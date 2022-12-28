#include<iostream>
#include"patch/domGrid/domIO.hpp"
#include<random>

//sfm的数据
using SfmData=openMVG::sfm::SfM_Data;
using DomGrid=openMVG::sfm::DomInfo;
using PointCloud=openMVG::sfm::Landmarks;
using Landmark=typename PointCloud::mapped_type;

int main()
{
    //初始化一个DOM网格
    DomGrid grid;
    //随机做一些点云用来初始化网格
    PointCloud tempCloud;
    std::uniform_real_distribution<double> maker(0,100);
    std::default_random_engine engine;
    for(unsigned id=0;id<10000;++id)
    {
        //随机做一个点
        auto tempPoint=Landmark();
        //记录点坐标
        tempPoint.X[0]=maker(engine);
        tempPoint.X[1]=maker(engine);
        tempPoint.X[2]=0;
        //记录到哈希表中
        tempCloud[id]=tempPoint;
    }
    //用这个点云来初始化dom
    grid.pixelLength_=0.1;
    grid.initDom(tempCloud);
    //保存网格单元的信息
    DomIO ioTool;
    ioTool.saveDom(grid,"/media/cvlab/data/temp/dom.bin");
	return 0;
}
