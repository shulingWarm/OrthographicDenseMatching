#pragma once
#include"patchHeader/initDOM/transferDom2DomByCloud.hpp"
#include"patchHeader/initDOM/mvgInitDomByCloud.hpp"
#include"patchHeader/pointcloud/transDom2CloudWithScore.hpp"
#include"patchHeader/initDOM/initDomWithCloudWithScoreDown.hpp"


//用mvg的方法从一个DOM初始化另一个DOM
class MvgInitDomByExistByCloud : virtual public MakeDomByTransferCloud<DomGrid,PointCloud>,
        virtual public DomInitByCloudWithScoreDown,//以前用的这个: AddCloudToInitedDom
        virtual public TransDom2CloudWithScore
{

};
