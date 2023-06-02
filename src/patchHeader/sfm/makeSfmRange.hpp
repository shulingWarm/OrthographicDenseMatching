#pragma once
#include"patchHeader/sfm/getInfoFromSfm.hpp"
#include<memory>

//传入sfm，生成它的viewPrior的范围
class SfmRangeMaker
{
public:
    using Mat3=Eigen::Matrix3d;
    using Vec3=Eigen::Vector3d;
    using SfmData=openMVG::sfm::SfM_Data;

    std::shared_ptr<SfmInterface> sfmInterface_;

    SfmRangeMaker()
    {
        sfmInterface_.reset(new SfmInterface());
    }

    //range是4个数 x_min x_max y_min y_max
    virtual void makeRange(SfmData& sfm,double* dstRange)
    {
        //初始化标志
        bool initFlag=true;
        //遍历范围
        for(auto& eachView : sfm.views)
        {
            //获取光心坐标
            auto tempCenter=sfmInterface_->getPriorCenterByView(eachView.first,sfm);
            //是否需要初始化
            if(initFlag)
            {
                initFlag=false;
                //记录初始坐标
                dstRange[0]=tempCenter[0];
                dstRange[1]=tempCenter[0];
                dstRange[2]=tempCenter[1];
                dstRange[3]=tempCenter[1];
            }
            else
            {
                //更新最值
                if(tempCenter[0]>dstRange[1]) dstRange[1]=tempCenter[0];
                else if(tempCenter[0]<dstRange[0]) dstRange[0]=tempCenter[0];
                if(tempCenter[1]>dstRange[3]) dstRange[3]=tempCenter[1];
                else if(tempCenter[1]<dstRange[2]) dstRange[2]=tempCenter[1];
            }
        }
    }
};

//生成点的范围，可以再扩展
class PointRangeMaker : public SfmRangeMaker
{
public:
    using Mat3=Eigen::Matrix3d;
    using Vec3=Eigen::Vector3d;
    using SfmData=openMVG::sfm::SfM_Data;

    //range是4个数 x_min x_max y_min y_max
    virtual void makeRange(SfmData& sfm,double* dstRange) override
    {
        //初始化标志
        bool initFlag=true;
        //遍历范围
        for(auto& eachView : sfm.structure)
        {
            //获取光心坐标
            auto tempCenter=eachView.second.X;
            //是否需要初始化
            if(initFlag)
            {
                initFlag=false;
                //记录初始坐标
                dstRange[0]=tempCenter[0];
                dstRange[1]=tempCenter[0];
                dstRange[2]=tempCenter[1];
                dstRange[3]=tempCenter[1];
            }
            else
            {
                //更新最值
                if(tempCenter[0]>dstRange[1]) dstRange[1]=tempCenter[0];
                else if(tempCenter[0]<dstRange[0]) dstRange[0]=tempCenter[0];
                if(tempCenter[1]>dstRange[3]) dstRange[3]=tempCenter[1];
                else if(tempCenter[1]<dstRange[2]) dstRange[2]=tempCenter[1];
            }
        }
    }
};
