#pragma once
#include"patchHeader/grid/gridInterface.hpp"
#include"openMVG/sfm/sfm_data.hpp"

//对于sfm来说常用的网格操作
//其实是个函数大杂烩
class CommonSfmGrid
{
public:
    //每个网格点云数量的计数器
    using PointCounter=GridInterface<unsigned>;
    using SfmData=openMVG::sfm::SfM_Data;
    //点云位置的操作
    using PointCloud=openMVG::sfm::Landmarks;

    //统计每个网格里面的点云数量
    //传入的时候需要自己保证所有的点计数都归零
    void makePointCound(PointCloud& points,
                        PointCounter& counter)
    {
        //遍历所有的点
        for(auto& eachPoint : points)
        {
            //令对应的网格位置的数值+1
            counter.getDataBySpatialCoord(eachPoint.second.X.data())++;
        }
    }
};
