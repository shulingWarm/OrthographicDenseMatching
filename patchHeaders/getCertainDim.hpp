#pragma once

//获取某个特定的维度
//从0开始索引的
template<class T=double,int dim=2>
class GetCertainDim
{
public:
    //传入一个地址的开头，获取特定的维度
    const T* operator()(const T* datas) const
    {
        return &datas[dim];
    }
};

using GetDimZ=GetCertainDim<double,2>;
using GetDimX=GetCertainDim<double,0>;
using GetDimY=GetCertainDim<double,1>;
