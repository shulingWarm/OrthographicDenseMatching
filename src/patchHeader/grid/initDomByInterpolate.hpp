#pragma once
#include"openMVG/sfm/sfm_data.hpp"
#include"patchHeader/grid/gridHeight/makeHeightGrid.hpp"
#include"patchHeader/grid/interfaceByDomZ.hpp"

using namespace openMVG;
using namespace openMVG::sfm;

//用插值的方法初始化DOM
void initDomByInterpolate(SfM_Data& sfm,bool reloadFlag=true)
{
    //把网格单元保存在点云里
    if(reloadFlag)
    {
        sfm.getZAsCloud(0.1);
    }
    //把dom包装成网格的接口
    InterfaceByDomZ zInterface(sfm.domInfo_);
    //从点云初始化网格单元的接口
    GridHeightMaker gridMaker;
    //双线性插值的方法改用直接访问的形式，其实是假的双线性插值
    //gridMaker.indicateHeightMake_->bilinearInterpolate_.reset(new FakeInterpolate<double>());
    gridMaker.makeHeight(zInterface,sfm.structure);
}
