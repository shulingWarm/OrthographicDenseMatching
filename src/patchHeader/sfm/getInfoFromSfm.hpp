#pragma once
#include"openMVG/sfm/sfm_data.hpp"
#include<unordered_map>

//从sfm_data里面获取基本信息的扩展
//只是为了方便从里面获取信息
class SfmInterface
{
    using Mat3=Eigen::Matrix3d;
    using Vec3=Eigen::Vector3d;
    using SfmData=openMVG::sfm::SfM_Data;
    //每个view能看到的点数统计
    using ViewObsCount=std::unordered_map<unsigned,unsigned>;
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

    //根据架次pose标号和相机方向，生成对应的idView
    inline unsigned getIdViewByPoseAndDir(unsigned flightIdPose,unsigned idDir)
    {
        return flightIdPose+idDir*100000;
    }

    //获取同架次的其它方向的相机
    unsigned getOtherDirByIdView(unsigned idView,unsigned idDir)
    {
        return getPoseIdByView(idView)+idDir*100000;
    }
};
