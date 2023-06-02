#pragma once
#include"patchHeader/sfm/getInfoFromSfm.hpp"
#include<memory>

//根据范围删点
class PointDeleterByViewRange
{
public:
    using SfmData=openMVG::sfm::SfM_Data;

    std::shared_ptr<SfmInterface> sfmInterface_;

    //获取sfm的viewPrior范围
    void makeViewRange(SfmData& sfm,double* range)
    {
        //初始化标志
        bool initFlag=true;
        //遍历每个范围
        for(auto& eachView : sfm.views)
        {
            //当前位置的点
            auto tempCenter=sfmInterface_->getPriorCenterByView(eachView.first,sfm);
            //判断是否初始化过
            if(initFlag)
            {
                initFlag=false;
                range[0]=tempCenter[0];
                range[1]=tempCenter[0];
                range[2]=tempCenter[1];
                range[3]=tempCenter[1];
            }
            else
            {
                if(tempCenter[0]<range[0]) range[0]=tempCenter[0];
                else if(tempCenter[0]>range[1]) range[1]=tempCenter[0];
                if(tempCenter[1]<range[2]) range[2]=tempCenter[1];
                else if(tempCenter[1]>range[3]) range[3]=tempCenter[1];
            }
        }
    }

    //计算viewPrior的平均高度
    double priorAvgHeight(SfmData& sfm)
    {
        //初始化高度求和
        double heightSum=0;
        //遍历viewPrior的高度
        for(auto& eachView : sfm.views)
        {
            //叠加高度
            auto tempCenter=sfmInterface_->getPriorCenterByView(eachView.first,sfm);
            heightSum+=tempCenter[2];
        }
        return heightSum/sfm.views.size();
    }

    //根据范围删除点
    void deletePointByRange(SfmData& sfm,double* range)
    {
        //遍历删除点
        for(auto iter=sfm.structure.begin();iter!=sfm.structure.end();)
        {
            //当前的点坐标
            auto point=iter->second.X;
            //判断是否在范围内
            if(point[0]>=range[0] && point[0]<range[1] &&
                    point[1]>=range[2] && point[1]<range[3] &&
                    point[2]>=range[4] && point[2]<range[5])
            {
                ++iter;
            }
            else
            {
                iter=sfm.structure.erase(iter);
            }
        }
    }

    //删除点，自己找范围然后自己删除点
    virtual void deletePoints(SfmData& sfm)
    {
        //查找范围
        double range[6];
        makeViewRange(sfm,range);
        //记录高度的范围
        range[5]=priorAvgHeight(sfm);
        range[4]=range[5]-500;
        //根据范围删除点
        deletePointByRange(sfm,range);
    }
};
