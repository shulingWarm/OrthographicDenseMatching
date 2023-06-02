#pragma once
#include"findMapNewKey.hpp"
#include"judgeViewValid.hpp"

//向一个sfmData里面添加一个点
template<class SfM_Data,
         class Landmark, //点的数据类型
         class Landmarks,//sfm里面的所有点的数据类型
         class FindNewKeyFunctor=FindMapNewKey<Landmarks>,
         class JudgeViewFunctor=JudgeViewValid<SfM_Data>
>
class AddPointToSfm
{

public:

    FindNewKeyFunctor findNewKey;
    //判断sfm里面是否存在某个view
    JudgeViewFunctor judgeViewValid;

    //向sfmData里面添加一个点
    void operator()(SfM_Data& targetSfm,
                    const Landmark& point)
    {
        //获取要添加的键值
        auto newKey=findNewKey(targetSfm.structure);
        //添加一个新的键值
        auto& newPoint=targetSfm.structure[newKey];
        //记录坐标点
        newPoint.X=point.X;
        //遍历记录点的obs
        for(auto& eachObs : point.obs)
        {
            //判断是否有这个view可以使用
            if(judgeViewValid(targetSfm,eachObs.first))
            {
                //在新添加的点里面记录这个obs
                newPoint.obs[eachObs.first]=eachObs.second;
            }
        }
    }
};
