#pragma once

//DOM网格的基类
class WidthHeightInterface
{
public:
    //获取网格的高度
    virtual unsigned getHeight() const=0;
    //获取网格的宽度
    virtual unsigned getWidth() const =0;
};
