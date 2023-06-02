#pragma once
#include"patchHeader/sfm/deletePointByRange.hpp"
#include"patchHeader/grid/gridHeight/makeHeightBySfm.hpp"

//通过高度网格删除点
class PointDeleterByHeightGrid : public PointDeleterByViewRange
{
public:
    using SfmData=openMVG::sfm::SfM_Data;
    //带有高度信息的网格单元
    using HeightGrid=GridByVector<double>;

    //生成网格用的工具
    std::shared_ptr<GridHeightMakerBySfm> makeGrid_;
    //双线性插值的方法获取的每个位置的高度
    std::shared_ptr<BilinearInterpolate<double>> bilinearInterpolate_;

    PointDeleterByHeightGrid()
    {
        makeGrid_.reset(new GridHeightMakerBySfm());
        //直接从class里面复制就可以
        bilinearInterpolate_=makeGrid_->makeByPoint_->indicateHeightMake_->bilinearInterpolate_;
    }

    //删除点时使用的分辨率，每个网格的大小
    double resolution_=5;
    //每个点到平均高度的阈值
    double disThreshold_=1;

    //用网格删除点云
    void deleteByGrid(SfmData& sfm,HeightGrid& grid)
    {
        //遍历sfm里面的每个点
        for(auto iter=sfm.structure.begin();iter!=sfm.structure.end();)
        {
            //当前的点坐标
            double* coord=iter->second.X.data();
            //判断点是否在范围内
            if(!grid.isSpatialPointInRange(coord))
            {
                iter=sfm.structure.erase(iter);
                continue;
            }
            //获取当前位置的平均高度
            double avgZ=bilinearInterpolate_->getDataAt(grid,coord);
            //判断z和平均值的差值是否达标
            if(std::abs(avgZ-coord[2])>disThreshold_)
            {
                iter=sfm.structure.erase(iter);
            }
            else
            {
                ++iter;
            }
        }
    }

    //删除点，自己找范围然后自己删除点
    virtual void deletePoints(SfmData& sfm) override
    {
        //根据sfm做高度范围的网格
        HeightGrid grid;
        makeGrid_->makeHeightGrid(sfm,grid,resolution_);
        //用分辨率删除点云
        deleteByGrid(sfm,grid);
    }
};
