#pragma once

//通常来说是需要结合sfm_data做操作的，因为这里仍然使用sfmdata
template<class SfmData,
         class ViewFunctor //对每个view的pair做的操作,同时还会传入sfm_data
>
class IterateView
{
public:
    //对每个view的具体操作
    ViewFunctor viewFunc;

    //传入需要被操作的view列表
    void operator()(SfmData& dstSfm)
    {
        //遍历每个view
        for(auto& eachView : dstSfm)
        {
            //调用函数对当前的view做修改
            viewFunc(eachView,dstSfm);
        }
    }
};
