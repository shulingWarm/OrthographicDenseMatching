#pragma once

//用一个点云来初始化一个DOM
template<class DomGrid,class PointCloud>
class InitDomByCloud
{
public:
    //用一个点云来初始化一个DOM
    virtual void initDomByCloud(DomGrid& grid,PointCloud& refCloud)=0;
};
