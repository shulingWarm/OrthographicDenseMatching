#pragma once
#include"patchHeader/grid/gridHeight/makeHeightGrid.hpp"
#include"patchHeader/grid/gridByVector.hpp"
#include"patchHeader/sfm/getInfoFromSfm.hpp"
#include"patchHeader/sfm/makeSfmRange.hpp"
#include<memory>

//通过sfm生成heightGrid
class GridHeightMakerBySfm
{
public:
    //带有高度信息的网格单元
    using HeightGrid=GridByVector<double>;
    //点云的数据类型
    using PointCloud=openMVG::sfm::Landmarks;
    //sfm
    using SfmData=openMVG::sfm::SfM_Data;

    std::shared_ptr<SfmInterface> sfmInterface_;
    std::shared_ptr<GridHeightMaker> makeByPoint_;
    std::shared_ptr<SfmRangeMaker> makeRange_;

    GridHeightMakerBySfm()
    {
        sfmInterface_.reset(new SfmInterface());
        makeByPoint_.reset(new GridHeightMaker());
        makeRange_.reset(new SfmRangeMaker());
    }

    //从sfm的数据生成高度网格信息
    void makeHeightGrid(SfmData& sfm,HeightGrid& grid,double spatialResolution)
    {
        //生成sfm的范围
        double range[4];
        makeRange_->makeRange(sfm,range);
        //初始化网格
        grid.reset(spatialResolution,(range[1]-range[0])/spatialResolution,
                (range[3]-range[2])/spatialResolution,
                range[0],range[2]);
        //用点云生成最终的网格数据
        makeByPoint_->makeHeight(grid,sfm.structure);
    }
};
