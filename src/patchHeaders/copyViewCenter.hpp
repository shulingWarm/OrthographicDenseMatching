#pragma once
#include"judgeViewValid.hpp"
#include"changeViewCenter.hpp"

template<class SfmData,
         class ViewPrior,//先验的位姿的指针
         class ChangeViewFunctor=
         ChangeViewCenter<SfmData,SfmData::views::value_type,ViewPrior>
>
class CopyViewCenter
{
public:

    ChangeViewCenter changeCenter;

    //把view的光心坐标从sfm2复制到sfm1
    void operator()(SfmData& sfm1,
                    SfmData& sfm2)
    {
        //遍历sfm2的每个view
        for(auto& eachView : sfm2.views)
        {
            //修改当前的光心坐标
            changeCenter(sfm1,eachView);
        }
    }
};
