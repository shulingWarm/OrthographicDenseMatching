#pragma once
#include"patch/grid/gridInterface.hpp"
#include<Eigen/Core>

//判断点到相机是否存在遮挡
class OcclusionJudger
{
public:

    //DOM网格单元，但只是根据对应位置获取它的Z值
    using DomZMat=GridInterface<double>;
    using Vec3=Eigen::Vector3cd;

    //检查时的分段数
    int stepSegNum=5;
    //差别时使用的比例
    double lineRate_=0.5;//比如投影线是y=kx 那么如果一个线的高度越过了0.5kx就认为是遮挡了

    //根据起始位置，步长和步数，判断遮挡
    bool isOccluded(DomZMat& domZ,
                   Vec3 xyzWorld,
                   Vec3 eachStepVec,
                   unsigned stepNum)
    {
        //遍历步长中的每个位置
        for(unsigned idStep=0;idStep<=stepNum;++idStep)
        {
            //当前的判别位置
            auto currLocation=xyzWorld+(idStep+1)*eachStepVec;
            //获取当前位置的DOM高度
            double zValueOfCurrLocation=domZ.getDataBySpatialCoord(currLocation.data());
            //2022-10-26 投影线在此处的高度必须比此处的DOM高度高出这么多
            double stepHeightHold=(currLocation[2]-xyzWorld[2])*lineRate_;
            //2022-10-26 为了更安全的去遮挡算法，更新处理的判断逻辑
            if(currLocation[2]-zValueOfCurrLocation<stepHeightHold)
            {
                return true;
            }
        }
        return false;
    }

    bool isOccluded(DomZMat& domZ,
                    double* xyzWorld,//三维的世界点
                    double* xyzCamera, //三维的相机光心
                    double zMax //z的最大值
    )
    {
        //计算两个点之间的向量
        Vec3 vectorToCam;
        for(int i=0;i<3;++i) vectorToCam[i]=xyzCamera[i]-xyzWorld[i];
        //场景的最大高度与当前位置的高度差
        double heightDiff=zMax-xyzWorld[2];
        //计算判别过程中每一步的向量
        Vec3 eachStepVec=vectorToCam*(heightDiff/(stepSegNum+1.0)/vectorToCam[2]);
        Vec3 worldPoint;
        for(int i=0;i<3;++i) worldPoint[i]=xyzWorld[i];
        //根据步长和起始位置判断是否遮挡
        return isOccluded(domZ,worldPoint,eachStepVec,stepSegNum);
    }
};
