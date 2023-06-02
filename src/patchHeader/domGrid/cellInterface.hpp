#pragma once

//DOM网格单元的抽象接口，可以获取Z值和目前的分数
template<class T>//一般不是double就是float
class CellInterface
{
public:
    //获取网格单元目前的高度
    virtual T& getCellHeight()=0;

    //获取网格单元目前的分数
    virtual T& getCellScore()=0;
};
