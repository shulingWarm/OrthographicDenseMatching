#pragma once
#include"openMVG/sfm/sfm_data.hpp"
#include"openMVG/sfm/switch.hpp"

using namespace openMVG;
using namespace openMVG::sfm;

//输入一个点云，对里面的每个网格单元计算一个初始化的分数
//这个时候需要确保每个高程位置都已经有了一个高度值
void initDOMScore(SfM_Data& sfm)
{
    std::cout<<"多分辨率网格初始化分数"<<std::endl;
    //要求所有的网格单元都用现有的高度参与计算
    GridConfig::useCurrentHeightFlag()=1;
    //记录它需要考虑遮挡的问题
    sfm.isUsingOccludeDetection_=true;
    //调用sfm直接对全体数据作一次窗口操作
    sfm.optimizeDomPixel(0);
    //把固定z高度的标志弄成0
    GridConfig::useCurrentHeightFlag()=0;
}
