#pragma once

//DOM网格的外部调用接口
template<class DomGrid,
         class DomCell>
class DomExtInterface
{
public:
    //获取网格的宽度或高度
    //dim为0的时候返回宽度，否则返回高度
    virtual unsigned getSize(DomGrid& grid,unsigned dim)=0;

    //获取指定位置的网格单元
    virtual DomCell& getCellAt(DomGrid& grid,unsigned x,unsigned y)=0;

    //获取网格单元的间隔
    virtual double getCellGap(DomGrid& grid)=0;
};
