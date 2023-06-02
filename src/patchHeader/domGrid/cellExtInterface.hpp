#pragma once

//DOM的外置调用抽象接口，这样会更灵活一些
template<class DOMCell,class T=double>
class CellExtInterface
{
public:
    //获取z值
    virtual T& getZValue(DOMCell& cell)=0;

    //获取DOM的分数
    virtual T& getCellScore(DOMCell& cell)=0;
};
