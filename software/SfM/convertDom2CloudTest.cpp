#include<iostream>
#include"patch/domGrid/domIO.hpp"
#include"patch/pointcloud/transDom2CloudWithScore.hpp"

using DomGrid=openMVG::sfm::DomInfo;

int main()
{
    //读取一个网格数据
    DomGrid grid;
    DomIO ioTool;
    std::string filePath="/media/cvlab/data/workSpace/mainProject/topViewConstruct/workspace/hekouzhenDOMBlock/result012/grid.bin";
    ioTool.loadDom(grid,filePath);
    //把DOM转换成点云
    TransDom2CloudWithScore transTool;
    PointCloud cloud;
    transTool.transDom2Cloud(grid,cloud);
	return 0;
}
