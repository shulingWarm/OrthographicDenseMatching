#pragma once
#include"patch/grid/gridInterface.hpp"
#include"openMVG/sfm/sfm_data.hpp"

//把一个DOM网格包装成网格的抽象接口
//这里专门用于获取Z值
template<class Cell>
class InterfaceByDom : public GridInterface<Cell>
{
public:
    using DomGrid=openMVG::sfm::DomInfo;

protected:
    //提供数据的指针
    DomGrid* domPtr_;

public:

    //传入一个DOM数据，后面获取数据都从这里面取
    InterfaceByDom(DomGrid& dom)
    {
        domPtr_=&dom;
    }

    //获取指定维度的单元格分辨率，一般就是0或1
    virtual double getSpatialResolution(int dim) override
    {
        return domPtr_->pixelLength_;
    }

    //获取指定维度的最小值
    virtual double getDimMin(int dim) override
    {
        return domPtr_->cloudRange_[dim];
    }

    //获取每个维度的离散长度
    virtual int getDimSize(int dim) override
    {
        if(dim==0) return domPtr_->domWidth_;
        return domPtr_->domHeight_;
    }
};
