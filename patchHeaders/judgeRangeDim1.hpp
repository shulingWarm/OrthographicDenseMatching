#pragma once
#include"getCertainDim.hpp"
#include"cmpLess.hpp"

//仅通过一个维度的多或少来对范围作判断
template<class T=double,
         class CompareFunctor=CompareLess<T>, //基本的二值比较
         class GetDimFunctor=GetDimZ //根据一个单个的维度来获取相应的范围
>
class JudgePointDim1
{
public:

    //每个函数都需要有一个对象
    CompareFunctor cmpFunctor;
    GetDimFunctor getDim;

    //传入的范围数据只是一个Z值的地址，传入的点是一个三维的范围
    bool operator()(const T* zPtr,const T* point) const
    {
        //获取某个点对应的维度值，然后再作是否可保留的判断
        return cmpFunctor(getDim(point),zPtr);
    }
};
