#pragma once
#include"patch/domGrid/getGridRange.hpp"
#include"patch/domGrid/gridWidthHeight.hpp"

template<class GridCell>
class DomInterfaceBase : public WidthHeightInterface,
        public RangeInterface<double>
{
public:
    //获取每个网格单元间距的接口
    virtual double getCellGap()=0;

    //获取指定位置的网格单元
    virtual GridCell& getCellAt(unsigned idRow,unsigned idCol)=0;

    //把坐标从连续域转换到离散域
    template<class WorldPoint,class GridPoint>
    void transWorld2Grid(const WorldPoint& srcPoint,GridPoint& dstPoint)
    {
        dstPoint[0]=(srcPoint[0]-this->getRange(0,0))/this->getCellGap();
        dstPoint[1]=(srcPoint[1]-this->getRange(1,0))/this->getCellGap();
    }

    //把坐标从离散域转到连续域
    template<class GridPoint,class WorldPoint>
    void transGrid2World(const GridPoint& srcPoint,WorldPoint& dstPoint)
    {
        dstPoint[0]=srcPoint[0]*this->getCellGap()+this->getRange(0,0);
        dstPoint[1]=srcPoint[1]*this->getCellGap()+this->getRange(1,0);
    }
};
