#pragma once
#include"initDomByExist.hpp"
#include"initDomByCloud.hpp"
#include"patchHeader/pointcloud/transDom2Pointcloud.hpp"

//把一个点云转换成点云，然后再初始化一个新的DOM
template<class DomGrid,class PointCloud>
class MakeDomByTransferCloud : virtual public InitDomByExist<DomGrid>,
        virtual public InitDomByCloud<DomGrid,PointCloud>,
        virtual public TransDom2Cloud<DomGrid,PointCloud>
{
public:
    //实现用已有DOM初始化一个新的DOM
    void initDom(DomGrid& srcGrid,DomGrid& dstGrid) override
    {
        //新建一个点云
        PointCloud cloud;
        //把DOM转换成cloud
        this->transDom2Cloud(srcGrid,cloud);
        //用得到的点云初始化目标网格单元
        this->initDomByCloud(dstGrid,cloud);
    }
};
