#pragma once

//把一个DOM转换成点云的形式
template<class DomGrid,class PointCloud>
class TransDom2Cloud
{
public:
    //把一个DOM转换成点云的形式
    virtual void transDom2Cloud(DomGrid& srcGrid,PointCloud& dstCloud)=0;
};
