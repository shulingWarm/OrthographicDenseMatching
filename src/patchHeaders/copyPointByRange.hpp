#pragma once
#include"cutPoint.hpp"
#include"copyPoint.hpp"

//2022-10-17 通过一个范围来对点做复制
template<class SfM_Data,
         class Landmakrs=decltype(SfM_Data::structure),
         class Landmark=Landmakrs::key_type,
         class CutPointFunctor=CutPoint<Landmark,Landmakrs>,
         class CopyPointFunctor=CopySfmDataPoint<SfM_Data,Landmark,Landmakrs>
>
class CopyPointByRange
{
public:

    CutPointFunctor cutPoints;
    CopyPointFunctor copyPoints;

    void operator()(SfM_Data& targetSfm,//点会被复制到它的structure里面
                    const Landmakrs& refPoints,//需要被复制的点
                    const double* range//需要复制的点的范围
    )
    {
        //复制一份点准备用来切割
        auto cutedPoints=refPoints;
        //对点做切割
        cutPoints(cutedPoints,range);
        //把切割过的点复制到sfm里面
        copyPoints(targetSfm,cutedPoints);
    }
};
