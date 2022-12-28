#pragma once
#include"judgeViewValid.hpp"

//传入sfm_data信息和一个view信息，把sfm_data里面的view_id改成对应的光心坐标
template<class SfmData,
         class ViewPair,
         class ViewPrior,
         class JudgeViewValidFunctor=JudgeViewValid<SfmData>
>
class ChangeViewCenter
{
public:

    JudgeViewValidFunctor isViewValid;

    //把view复制到sfm里面
    void operator ()(SfmData& dstSfm,
                     const ViewPair& viewInfo
                     )
    {
        //判断sfm里面有没有这个view
        if(isViewValid(dstSfm,viewInfo.first))
        {
            //获取光心坐标
            auto& priorCenter=dynamic_cast<ViewPrior*>(viewInfo.second.get())->pose_center;
            //存储目标的光心坐标
            auto& dstCenter=dynamic_cast<ViewPrior*>(
                        dstSfm.views.at(viewInfo.first).get())->pose_center;
            //对光心坐标赋值
            dstCenter=priorCenter;
        }
    }
};
