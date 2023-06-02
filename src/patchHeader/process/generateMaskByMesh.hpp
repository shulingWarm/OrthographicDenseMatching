#pragma once

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/cameras/Cameras_Common_command_line_helper.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_report.hpp"
#include "openMVG/system/timer.hpp"

using namespace openMVG;
using namespace openMVG::sfm;

#include<array>

using Face=std::array<Eigen::Vector3f,3>;

Eigen::Vector3f ComputeNormal(const std::array<Eigen::Vector3f, 3>& points) {
  Eigen::Vector3f v1 = points[1] - points[0];
  Eigen::Vector3f v2 = points[2] - points[0];
  return v1.cross(v2).normalized();
}

void readMesh(const std::string& meshPath,
              std::vector<Face>& faceList)
{
    //打开文件流
    std::fstream fileHandle;
    fileHandle.open(meshPath,std::ios::in|std::ios::binary);
    //读取face的个数
    unsigned faceNum;
    fileHandle.read((char*)&faceNum,sizeof(unsigned));
    //提前给face开辟空间
    faceList.resize(faceNum);
    //遍历每个要读取的face
    for(auto& eachFace : faceList)
    {
        //依次读取面片上的三个点
        for(unsigned idVertex=0;idVertex<3;++idVertex)
        {
            //读取三个数
            fileHandle.read((char*)(eachFace[idVertex].data()),sizeof(float)*3);
        }
    }
}

//把点投影在图片上
Eigen::Vector2d getViewProjection(const Eigen::Vector3d& point3d, const SfM_Data& sfm_data, const openMVG::IndexT id_view)
{
    // 获取指定id_view的相机参数
    const openMVG::cameras::IntrinsicBase* intrinsic = sfm_data.intrinsics.at(sfm_data.GetViews().at(id_view)->id_intrinsic).get();
    const openMVG::geometry::Pose3 pose = sfm_data.GetPoseOrDie(sfm_data.views.at(id_view).get());
    const Eigen::Matrix<double, 3, 4> projection_matrix = intrinsic->get_projective_equivalent(pose);

    // 将三维点投影到相机平面
    const Eigen::Vector4d point3d_homogeneous = point3d.homogeneous();
    const Eigen::Vector3d point2d_homogeneous = projection_matrix * point3d_homogeneous;
    const Eigen::Vector2d point2d = point2d_homogeneous.hnormalized();
    return point2d;
}

//判断点是否在图片的范围上
bool isPointInRange(const Eigen::Vector2d& point2d, const SfM_Data& sfm_data, const openMVG::IndexT id_view)
{
    // 获取图像大小
    const View* view = sfm_data.GetViews().at(id_view).get();
    const int width = view->ui_width;
    const int height = view->ui_height;

    // 判断投影点是否在图像范围内
    const double margin = 0.01;  // 增加一定的margin
    if (point2d.x() >= -margin && point2d.x() <= width - 1 + margin &&
        point2d.y() >= -margin && point2d.y() <= height - 1 + margin)
    {
        return true;
    }
    else
    {
        return false;
    }
}

//根据三个投影点在图片上画三角形
void DrawTriangle(const std::array<Eigen::Vector2d, 3>& projFace, cv::Mat& img) {
  std::vector<cv::Point> points(3);
  for (int i = 0; i < 3; i++) {
    points[i] = cv::Point(projFace[i].x(), projFace[i].y());
  }
  cv::fillConvexPoly(img, points, cv::Scalar(0));
}

//把一个face在图片上做标记
void labelFaceAtImg(SfM_Data& sfm,Face& dstFace)
{
    //范围
    const float range[]={-0.8,20,-116,-96};
    //遍历所有的图片
    for(auto& eachView : sfm.views)
    {
        if(!MID_ID_JUDGE(eachView.first)) continue;
        if(sfm.poses.count(eachView.second->id_pose)==0) continue;
        //遍历每个点确保都在图片范围内
        bool inFlag=true;
        //三个点的投影位置
        std::array<Eigen::Vector2d,3> projFace;
        for(unsigned idVertex=0;idVertex<3;++idVertex)
        {
            auto& eachPoint = dstFace[idVertex];
            //先把点转换成double形式
            Eigen::Vector3d tempDouble={eachPoint[0],eachPoint[1],eachPoint[2]};
            //判断点是否在范围内
            if(tempDouble[0]> range[1] || tempDouble[0]<range[0] ||
                    tempDouble[1]>range[3] || tempDouble[1]<range[2])
            {
                inFlag=false;
                break;
            }
            //把点投影到图片上
            projFace[idVertex]=getViewProjection(tempDouble,sfm,eachView.first);
            //判断投影点是否在范围内
            if(!isPointInRange(projFace[idVertex],sfm,eachView.first))
            {
                inFlag=false;
                break;
            }
        }
        //把三角形画到对应的原图上
        if(inFlag)
        {
            DrawTriangle(projFace,eachView.second->cvImg_);
            eachView.second->usedFlag=true;
        }
    }
}

//把mask加在原图上然后覆盖
void exportMask(SfM_Data& sfm,
                const std::string& maskFolder)
{
    //遍历所有的view
    std::vector<View*> viewList;
    std::vector<unsigned> idList;
    for(auto& eachView : sfm.views)
    {
        if(!MID_ID_JUDGE(eachView.first)) continue;
        if(sfm.poses.count(eachView.second->id_pose)==0) continue;
        viewList.push_back(eachView.second.get());
        idList.push_back(eachView.first);
    }
    //遍历所有需要处理的图片
#pragma omp parallel for
    for(unsigned id=0;id<idList.size();++id)
    {
        if(viewList[id]->usedFlag==false) continue;
        //读取它的原图
        cv::Mat realImg=cv::imread(sfm.s_root_path+"/"+
                                   viewList[id]->s_Img_path);
        auto& maskImg=viewList[id]->cvImg_;
        //把mask中白色的部分弄成黑色
        realImg.setTo(cv::Scalar(0,0,0),maskImg);
        //把结果保存到目标目录中
        cv::imwrite(maskFolder+"/"+std::to_string(idList[id])+".JPG",realImg);
    }
}

//输入一个sfm,把每个位置的原图都记录下来
void labelMeshMask(SfM_Data& sfm,
                   const std::string& meshPath, //mesh path是一种自定义的逻辑
                   const std::string& maskFolder
)
{
    //把所有的中相机的图片都读取进来，但是完全标记为黑色
    for(auto& eachView : sfm.views)
    {
        if(!MID_ID_JUDGE(eachView.first)) continue;
        if(sfm.poses.count(eachView.second->id_pose)==0) continue;

        eachView.second->cvImg_.create(eachView.second->ui_height,
                                       eachView.second->ui_width,CV_8UC1);
        eachView.second->cvImg_.setTo(cv::Scalar(255));
    }
    //期望的法向量
    Eigen::Vector3f expNorm={0,0,1};
    //期望的法向量点乘的阈值
    const float expDotProduct=0.8;
    //读取mesh
    std::vector<Face> faceList;
    readMesh(meshPath,faceList);
    //遍历每个face
#pragma omp parallel for
    for(unsigned idFace=0;idFace<faceList.size();++idFace)
    {
        auto& eachFace=faceList[idFace];
        //计算face的法向量
        auto tempNorm=ComputeNormal(eachFace);
        //判断两个法向量之间的点乘
        auto tempDot=std::abs(tempNorm.transpose()*expNorm);
        //如果点乘的结果满足期望就把点投影到原始的图片上做mask的标记
        if(tempDot>expDotProduct)
        {
            labelFaceAtImg(sfm,eachFace);
        }
    }
    //对已经生成的mask区域做反向选择
    exportMask(sfm,maskFolder);
}
