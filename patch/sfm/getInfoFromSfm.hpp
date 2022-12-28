#pragma once
#include"openMVG/sfm/sfm_data.hpp"

//从sfm_data里面获取基本信息的扩展
//只是为了方便从里面获取信息
class SfmInterface
{
    using Mat3=Eigen::Matrix3d;
    using Vec3=Eigen::Vector3d;
    using SfmData=openMVG::sfm::SfM_Data;
public:
    //从内参信息里面获取内参矩阵的方法
    Mat3 getIntrMatByParams(openMVG::cameras::IntrinsicBase* params)
    {
        //获取内参的参数
        auto intrParam=params->getParams();
        //新建一个临时的内参矩阵
        Mat3 tempMat;
        tempMat<<intrParam[0],  0,                          intrParam[1],
                                0,                          intrParam[0],   intrParam[2],
                                0,                          0,                          1;
        //返回计算出来的内参矩阵
        return tempMat;
    }

    //从指定的view里面获取内参矩阵
    Mat3 getIntrMatByView(unsigned idView,SfmData& sfm)
    {
        return getIntrMatByParams(sfm.intrinsics.at(
                                      sfm.views.at(idView)->id_intrinsic).get());
    }

    //获取指定view的pose里的center
    Vec3 getCenterByView(unsigned idView,SfmData& sfm)
    {
        return sfm.poses.at(sfm.views.at(idView)->id_pose).center();
    }

    //根据view的标号获取相机光心
    Vec3 getPriorCenterByView(unsigned idView,SfmData& sfm)
    {
        //prior形式的指针
        auto priorPtr=dynamic_cast<openMVG::sfm::ViewPriors*>(sfm.views.at(idView).get());
        //不做安全检查，直接返回对应的center
        return priorPtr->pose_center_;
    }

    //从指定的view标号里面获取相机外参
    void getRtByView(unsigned idView,SfmData& sfm,
                     Mat3& rMat,Vec3& center)
    {
        //目标的pose
        const auto& pose=sfm.poses.at(sfm.views.at(idView)->id_pose);
        //给rt赋值
        rMat=pose.rotation();
        center=pose.center();
    }

    //判断一个投影坐标是不是在图片的范围内
    bool isInImageRange(double* xy,SfmData& sfm,unsigned idView)
    {
        //目标的view
        auto viewPtr=sfm.views.at(idView);
        return xy[0]>=0 && xy[0]<viewPtr->ui_width &&
                xy[1]>=0 && xy[1]<viewPtr->ui_height;
    }

    //计算点在某个相机上的投影位置 最后得到的是一个齐次点的坐标
    Vec3 getProjectionByView(unsigned idView,SfmData& sfm,const Vec3& worldPoint)
    {
        //获取内参和外参
        auto intr=getIntrMatByView(idView,sfm);
        //获取外参
        Vec3 center;
        Mat3 rotMat;
        getRtByView(idView,sfm,rotMat,center);
        //计算投影的位置
        Vec3 projection=intr*rotMat*(worldPoint-center);
        //做成非齐次点
        if(projection[2]!=0)
        {
            for(int i=0;i<3;++i) projection[i]/=projection[2];
        }
        return projection;
    }

    //获取sfm指定view的旋转矩阵，不做存在检查
    Mat3& getRotByView(unsigned idView,SfmData& sfm)
    {
        return sfm.poses.at(
                    sfm.views.at(idView)->id_pose).rotation();
    }

    //根据view获取相机的方向标号
    unsigned getDirByView(unsigned idView)
    {
        return idView/100000;
    }

    //根据view获取架次的标号
    unsigned getIdFlightByView(unsigned idView)
    {
        return getPoseIdByView(idView)/10000;
    }

    //根据view的标号获取它的pose标号，其实就是把相机标号去掉
    unsigned getPoseIdByView(unsigned idView)
    {
        return idView%100000;
    }

    //获取同架次的其它方向的相机
    unsigned getOtherDirByIdView(unsigned idView,unsigned idDir)
    {
        return getPoseIdByView(idView)+idDir*100000;
    }

    //判断某个view是否有pose
    bool isPoseExist(SfmData& sfm,unsigned idView)
    {
        if(sfm.views.count(idView)==0) return false;
        return sfm.poses.count(sfm.views.at(idView)->id_pose)>0;
    }

    //判断某个pose是否具备所有方向的pose
    bool isAllDirPoseExist(SfmData& sfm,unsigned idPose,unsigned dirNum=5)
    {
        //遍历每个方向
        for(unsigned idDir=1;idDir<=dirNum;++idDir)
        {
            //生成当前位置的idView
            unsigned idView=idPose+idDir*100000;
            //判断是否pose存在
            if(!isPoseExist(sfm,idView)) return false;
        }
        return true;
    }
};
