#pragma once
#include"openMVG/sfm/sfm_data.hpp"
#include"openMVG/sfm/sfm_landmark.hpp"

using namespace openMVG;
using namespace openMVG::sfm;

//2023-3-8 这个代码还没有写完，因为老师说他想用mesh的三角面片的法向量来决定一个哪些地方被标记为mask

//把所有的原图弄成全黑的
void makeAllViewImageBlank(SfM_Data& sfm)
{
    //遍历所有的view
    for(auto& eachView : sfm.views)
    {
        //看当前的view是否有pose
        if(sfm.poses.count(eachView.second->id_pose)==0) continue;
        //看当前的view的图片是否载入过
        if(eachView.second->cvImg_.empty()) continue;
        //把图片弄成黑色
        eachView.second->cvImg_.setTo(cv::Scalar(0,0,0));
    }
}

////输入一个点和id_view 判断这个点是否在这个view的投影范围内
//bool is_point_in_view(const Eigen::Vector3d& point3d, const SfM_Data& sfm_data, const openMVG::IndexT id_view)
//{
//    // 获取指定id_view的相机参数
//    const openMVG::cameras::IntrinsicBase* intrinsic = sfm_data.GetViews().at(id_view)->id_intrinsic.get();
//    const openMVG::geometry::Pose3 pose = sfm_data.GetPoseOrDie(id_view);
//    const Eigen::Matrix<double, 3, 4> projection_matrix = intrinsic->get_projective_equivalent(pose);

//    // 将三维点投影到相机平面
//    const Eigen::Vector3d point3d_homogeneous = point3d.homogeneous();
//    const Eigen::Vector3d point2d_homogeneous = projection_matrix * point3d_homogeneous;
//    const Eigen::Vector2d point2d = point2d_homogeneous.hnormalized();

//    // 获取图像大小
//    const View* view = sfm_data.GetViews().at(id_view).get();
//    const int width = view->ui_width;
//    const int height = view->ui_height;

//    // 判断投影点是否在图像范围内
//    const double margin = 0.01;  // 增加一定的margin
//    if (point2d.x() >= -margin && point2d.x() <= width - 1 + margin &&
//        point2d.y() >= -margin && point2d.y() <= height - 1 + margin)
//    {
//        return true;
//    }
//    else
//    {
//        return false;
//    }
//}


//把当前的点在所有的图片上投影
void projectAllPoints(SfM_Data& sfm,
                      Landmark& point)
{
    //遍历所有的view,看能不能投到这个点上
    for(auto& eachView : sfm.views)
    {
        //判断view的pose是否存在
        if(sfm.poses.count(eachView.second->id_view)==0) continue;
        //判断这个view的图片是否加载过
        if(eachView.second->cvImg_.empty()) continue;
        //判断能否投影到这个图片上

    }
}

//根据sfm的网格单元的结果，给所有的原图生成mask
void generateMask(SfM_Data& sfm)
{
    //把所有已经读取过的图片弄成全黑的
    makeAllViewImageBlank(sfm);
    //把网格单元里面的所有点弄成点云
    sfm.getZAsCloud(-2);
    //遍历所有的点
    for(unsigned idPoint=0;idPoint<sfm.structure.size();++idPoint)
    {
        //当前的点
        Landmark& currPoint=sfm.structure[idPoint];
        //把当前的点在所有的图片上投影
        projectAllPoints(sfm,currPoint);
    }
}
