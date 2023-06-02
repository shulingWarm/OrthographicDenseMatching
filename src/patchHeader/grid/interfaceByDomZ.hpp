#pragma once
#include"patchHeader/grid/gridInterface.hpp"
#include"openMVG/sfm/sfm_data.hpp"

//2023-2-22 过于紧急，代码复用还是他娘的复制粘贴更快
//把一个DOM网格包装成网格的抽象接口
//这里专门用于获取Z值
class InterfaceByDomZ : public GridInterface<double>
{
public:
    using DomGrid=openMVG::sfm::DomInfo;

protected:
    //提供数据的指针
    DomGrid* domPtr_;

    using Cell=double;

public:

    //传入一个DOM数据，后面获取数据都从这里面取
    InterfaceByDomZ(DomGrid& dom)
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

    //获取数据的z值
    //获取特定数据的接口
    //依次是xy
    virtual Cell& getDataAt(int* xy) override
    {
        if(xy[0]<0) xy[0]=0;
        if(xy[0]>= getDimSize(0)) xy[0]=getDimSize(0)-1;
        if(xy[1]<0) xy[1]=0;
        if(xy[1]>=getDimSize(1)) xy[1]=getDimSize(1)-1;
        //获取对应位置的网格单元
        auto& unit=domPtr_->getUnit(xy[0],xy[1]);
        //返回z值
        return unit.z_;
    }
};
