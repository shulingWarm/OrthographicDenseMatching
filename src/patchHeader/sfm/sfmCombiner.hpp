#pragma once
#include"patchHeader/sfm/getInfoFromSfm.hpp"
#include"patchHeader/sfm/poseCombiner.hpp"

//合并两个sfm的接口
class SfmCombiner
{
public:
    using SfmData=openMVG::sfm::SfM_Data;

    //复制view的工具
    std::shared_ptr<PoseCombiner> combine_;

    SfmCombiner()
    {
        combine_.reset(new PoseCombiner());
    }

    //把sfm2合并到sfm1里面
    void combineSfm(SfmData& sfm1,SfmData& sfm2)
    {
        //合并两个sfm的view和pose
        for(auto& eachView : sfm2.views)
        {

        }
    }
};
