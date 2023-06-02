#pragma once
#include"patchHeader/grid/gridInterface.hpp"
#include"openMVG/sfm/sfm_data.hpp"

//把高度网格转换成用点云来表示
class CloudConverterByGrid
{
public:
    using HeightGrid=GridInterface<double>;
    using SfmData=openMVG::sfm::SfM_Data;
    using Vec3=Eigen::Vector3d;

    //整理点云的点标号
    void sortPointId(SfmData& sfm)
    {
        //遍历点
        openMVG::sfm::Landmarks points;
        points.swap(sfm.structure);
        //遍历记录点
        for(auto& eachPoint : points)
        {
            //要添加的标号
            unsigned nextId=sfm.structure.size();
            sfm.structure[nextId]=eachPoint.second;
        }
    }

    //把一个xy坐标从网格里面拿出来然后送到sfm里面做点云
    void makePointToSfm(SfmData& sfm,HeightGrid& grid,
                        int* xy)
    {
        //生成一个三维坐标
        Vec3 point;
        point[2]=grid.getDataAt(xy);
        point[0]=grid.toSpatialDim(0,xy[0]);
        point[1]=grid.toSpatialDim(1,xy[1]);
        //获取要添加的点坐标
        unsigned idPoint=sfm.structure.size();
        //判断点是否存在
        if(sfm.structure.count(idPoint))
        {
            //整理一下点云里面的点标号
            sortPointId(sfm);
        }
        //新建一个点坐标
        auto& landmark=sfm.structure[idPoint];
        //只需要记录点坐标，不用管obs
        landmark.X=point;
    }

    //用高度类型的网格生成一个点云用于展示
    //在传入的时候需要保证sfm是整理过的，不然会出问题
    void makeCloud(SfmData& sfm,HeightGrid& grid)
    {
        //遍历每个位置的数据
        int xy[2];
        for(xy[0]=0;xy[0]<grid.getDimSize(0);++xy[0])
        {
            for(xy[1]=0;xy[1]<grid.getDimSize(1);++xy[1])
            {
                //把xy的当前位置弄成点云加入到grid里面
                makePointToSfm(sfm,grid,xy);
            }
        }
    }
};
