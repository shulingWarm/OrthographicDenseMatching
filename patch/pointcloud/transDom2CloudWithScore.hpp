#pragma once
#include"openMVG/sfm/sfm_data.hpp"
#include"patch/pointcloud/transDom2Pointcloud.hpp"
#include"patch/pointcloud/3DPointWithScore.hpp"

//dom网格的数据类型
using DomGrid=openMVG::sfm::DomInfo;
//点云的数据类型
using Point3D=PointWithScore;
using PointCloud=std::vector<Point3D>;
//网格单元的数据类型
using DomCell=openMVG::sfm::DomUnit;

//关键是要使用openmvg形式的点，并且要保留可扩展性
class TransDom2CloudWithScore : virtual public TransDom2Cloud<DomGrid,PointCloud>
{
protected:



    //把一个网格单元的信息记录到一个点里面
    Point3D transCell2Point(DomCell& cell,DomGrid& grid,unsigned x,unsigned y)
    {
        //新建一个三维点
        Point3D tempPoint;
        //点坐标
        auto pointPtr=tempPoint.getCoord();
        //把xy坐标转换成连续的坐标
        auto convertedPoint=grid.convertToContiousPoint<std::array<double,2>>(y,x);
        pointPtr[0]=convertedPoint[0];
        pointPtr[1]=convertedPoint[1];
        //从cell里面记录z值
        pointPtr[2]=cell.z_;
        tempPoint.getScore()=cell.nccScore_;
        //返回记录到的点
        return tempPoint;
    }
public:
    //把DOM转换成带分数的点云
    void transDom2Cloud(DomGrid& grid,PointCloud& cloud) override
    {
        //初始化点云数据
        cloud.reserve(grid.domHeight_*grid.domWidth_);
        //遍历点云里面的每个单元
        for(unsigned x=0;x<grid.domWidth_;++x)
        {
            for(unsigned y=0;y<grid.domHeight_;++y)
            {
                //当前的cell值
                auto& currCell=grid.getUnit(x,y);
                //把DOM信息保存到点云里面
                cloud.push_back(transCell2Point(currCell,grid,x,y));
            }
        }
    }
};
