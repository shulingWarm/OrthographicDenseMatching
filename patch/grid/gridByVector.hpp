#pragma once
#include"patch/grid/gridInterface.hpp"
#include<vector>
#include<array>

//用vector来记录数据的网格单元，这里面会实质性地做数据存储
//而不是像父类那样只做接口
template<class Cell>
class GridByVector : public GridInterface<Cell>
{
public:
    //网格单元里面的数据
    using CellRow=std::vector<Cell>;
    using CellData=std::vector<CellRow>;
    CellData data_;
    //行数和列数的存储
    std::array<unsigned,2> xySize_;
    //xy方向的最小值
    std::array<double,2> xyMin_;
    //行列的分辨率的存储
    double resolution_;

    //重新做的初始化
    void reset(double resolution,int xSize,int ySize,double xMin,double yMin)
    {
        resolution_=resolution;
        xySize_[0]=xSize;
        xySize_[1]=ySize;
        xyMin_[0]=xMin;
        xyMin_[1]=yMin;
    }

    //根据行列数和分辨率初始化网格单元
    GridByVector(double resolution,int xSize,int ySize,double xMin,double yMin)
    {
        reset(resolution,xSize,ySize,xMin,yMin);
    }

    //获取特定数据的接口
    //依次是xy
    virtual Cell& getDataAt(int* xy) override
    {
        return data_[xy[1]][xy[0]];
    }

    //获取指定维度的单元格分辨率，一般就是0或1
    virtual double getSpatialResolution(int dim) override
    {
        return resolution_;
    }

    //获取指定维度的最小值
    virtual double getDimMin(int dim) override
    {
        return xyMin_[dim];
    }

    //获取每个维度的离散长度
    virtual int getDimSize(int dim) override
    {
        return xySize_[dim];
    }
};
