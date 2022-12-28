#pragma once
#include"patch/initDOM/transferDom2DomByCloud.hpp"
#include"patch/initDOM/mvgInitDomByCloud.hpp"
#include"patch/pointcloud/transDom2CloudWithScore.hpp"
#include"patch/initDOM/initDomWithCloudWithScoreDown.hpp"


//用mvg的方法从一个DOM初始化另一个DOM
class MvgInitDomByExistByCloud : virtual public MakeDomByTransferCloud<DomGrid,PointCloud>,
        virtual public DomInitByCloudWithScoreDown,//以前用的这个: AddCloudToInitedDom
        virtual public TransDom2CloudWithScore
{

};
