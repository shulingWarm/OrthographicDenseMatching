#pragma once
#include"patch/grid/gridInterface.hpp"
#include"patch/imgSelect/getBestViewByKnownMat.hpp"
#include"patch/grid/gridByVector.hpp"

//获取每个网格单元的最佳view的总流程
class TotalBestViewMaker
{
public:
    //网格单元里面每个位置的最好的投影图片
    //键值存储的是一个图片标号
    using ViewMat=GridByVector<unsigned>;
    //DOM网格单元，但只是根据对应位置获取它的Z值
    using DomZMat=GridInterface<double>;
    using SfmData=openMVG::sfm::SfM_Data;
    using Vec3=Eigen::Vector3d;
    //DOM的网格数据
    using DomGrid=openMVG::sfm::DomInfo;

    //最后的结果会存储在dstViewMat里面
    void makeBestView(DomGrid& dom,ViewMat& dstViewMat)
    {

    }
};
