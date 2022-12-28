#include<iostream>
#include"patch/initDOM/mvgInitDomByExistDomByCloud.hpp"
#include"openMVG/sfm/sfm_data_io.hpp"
#include"patch/domGrid/domIO.hpp"

int main()
{
    //正常的稀疏点云
    std::string sfmPath="/media/cvlab/data/workSpace/mainProject/topViewConstruct/workspace/hekouzhenDOMBlock/feat16_mask_1/sfm_data.bin";
    //保存的dom网格数据
    std::string gridPath="/media/cvlab/data/workSpace/mainProject/topViewConstruct/workspace/hekouzhenDOMBlock/result012/grid.bin";
    //读取sfm数据
    using SfmData=openMVG::sfm::SfM_Data;
    SfmData sfmData;
    openMVG::sfm::Load(sfmData,sfmPath,openMVG::sfm::ESfM_Data::ALL);
    //设置分辨率
    sfmData.domInfo_.pixelLength_=0.04;
    //初始化点云
    sfmData.domInfo_.initDom(sfmData.structure);
    //读取已经生成过的网格单元
    DomGrid grid;
    {
        DomIO ioTool;
        ioTool.loadDom(grid,gridPath);
    }
    //已有dom的叠加器
    MvgInitDomByExistByCloud mvgInitTool;
    mvgInitTool.initDom(grid,sfmData.domInfo_);
	return 0;
}
