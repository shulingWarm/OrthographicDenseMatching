#pragma once
#include"cutPoint.hpp"
#include"cmpLess.hpp"
#include"judgeRangeDim1.hpp"
#include"getCertainDim.hpp"

//根据指定的Z值作切割的函数
template<class Landmakrs,class T=double>
using CutPointLessZ=CutPoint<Landmakrs::mapped_type,Landmakrs,T,
JudgePointDim1<T>
>;
