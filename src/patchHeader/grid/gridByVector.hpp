#pragma once
#include"patchHeader/grid/gridInterface.hpp"
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
        //初始化data数据的大小
        data_.resize(ySize);
        for(auto& eachRow : data_)
        {
            eachRow.resize(xSize,0);
        }

    }

    //根据行列数和分辨率初始化网格单元
    GridByVector(double resolution=1,int xSize=0,int ySize=0,double xMin=0,double yMin=0)
    {
        reset(resolution,xSize,ySize,xMin,yMin);
    }

    //获取特定数据的接口
    //依次是xy
    virtual Cell& getDataAt(int* xy) override
    {
        //如果超过范围了就按最接近的来
        int x=xy[0];
        if(x>=getDimSize(0)) x=getDimSize(0)-1;
        else if(x<0) x=0;
        int y=xy[1];
        if(y>=getDimSize(1)) y=getDimSize(1)-1;
        else if(y<0) y=0;
        return data_[y][x];
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
