#pragma once
#include"patchHeader/grid/gridInterface.hpp"

//只有一个数据的网格，其实就是用来欺骗相关接口的
template<class Cell>
class OnlyValueGrid : public GridInterface<Cell>
{
public:
    //被操作的固定值
    Cell data_;

    //获取特定数据的接口
    //依次是xy
    virtual Cell& getDataAt(int* xy) override
    {
        return data_;
    }

    //获取指定维度的单元格分辨率，一般就是0或1
    virtual double getSpatialResolution(int dim) override
    {
        return 1;
    }

    //获取指定维度的最小值
    virtual double getDimMin(int dim) override
    {
        return 0;
    }

    //获取每个维度的离散长度
    virtual int getDimSize(int dim) override
    {
        return 1;
    }

    //用离散的方式获取数据
    virtual Cell& getDataBySpatialCoord(double* xy) override
    {
        return data_;
    }

    //判断一个网格单元是否在范围内
    virtual bool isGridDimInRange(int dim,int x) override
    {
        return true;
    }

    //判断某个空间维度是否在范围内
    virtual bool isSpatialDimInRange(int dim,double x) override
    {
        return true;
    }

    //判断一个空间点是否在范围内
    virtual bool isSpatialPointInRange(double* xy) override
    {
        return true;
    }
};
