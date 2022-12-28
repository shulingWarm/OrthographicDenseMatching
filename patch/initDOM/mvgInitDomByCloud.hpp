#pragma once
#include"patch/initDOM/initDomByCloud.hpp"
#include"openMVG/sfm/sfm_data.hpp"
#include"patch/pointcloud/3DPointWithScore.hpp"

//dom的网格
using DomGrid=openMVG::sfm::DomInfo;
//三维的点
using Point3D=PointWithScore;
//三维点云
using PointCloud=std::vector<Point3D>;

//openmvg 在已经初始化过的DOM点云上叠加来自其它DOM的点云
class AddCloudToInitedDom : virtual public InitDomByCloud<DomGrid,PointCloud>
{
protected:
    //dom的网格单元
    using DomCell=openMVG::sfm::DomUnit;


protected:
    //判断一个分数是否远大于另一个
    bool judgeLargeBigger(double score1,double score2)
    {
        return score1-score2>0.15;
    }

    //用一个点的数据来代替cell的数据
    virtual void replaceCellData(DomCell& cell,Point3D& point)
    {
        cell.z_=point.getCoord()[2];
        cell.nccScore_=point.getScore();
    }

    //在一个网格单元里面记录z值
    //分数差不多的时候取z值较大的，否则按照分数来比较
    void addPointToCell(DomCell& cell,Point3D& point)
    {
        //看一下已有的cell的分数明显大于这个，那就算了
        if(!judgeLargeBigger(cell.nccScore_,point.getScore()))
        {
            //如果点的分数明显大于cell的分数，那也不看高度了
            if(judgeLargeBigger(point.getScore(),cell.nccScore_))
            {
                this->replaceCellData(cell,point);
            }
            else
            {
                //两个还需要比较一下高度
                if(point.getCoord()[2]>cell.z_)
                {
                    this->replaceCellData(cell,point);
                }
            }
        }
    }

    //在dom点云里面添加一个点
    void addPointToGrid(DomGrid& grid,Point3D& point)
    {
        //点的三维坐标
        auto coord=point.getCoord();
        //判断点是否在范围内
        if(!grid.judgeContiPointOutRange(coord[0],coord[1]))
        {
            //获取网格单元
            auto& currCell=grid.getUnitByContiPoint(coord);
            //在网格单元里面记录信息
            addPointToCell(currCell,point);
        }
    }
public:
    //用一个点云来初始化DOM
    void initDomByCloud(DomGrid& grid,PointCloud& refCloud) override
    {
        //遍历点云里面的每个点
        for(auto& eachPoint : refCloud)
        {
            //在dom里面添加一个点
            addPointToGrid(grid,eachPoint);
        }
    }
};
