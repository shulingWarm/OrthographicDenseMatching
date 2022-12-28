#pragma once

//获取网格范围的基类
template<class T>
class RangeInterface
{
public:
    //获取对应维度的范围数值
    const T& getRange(unsigned idDim,bool minFlag)=0;
};
