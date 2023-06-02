#pragma once
#include"patchHeader/grid/gridInterface.hpp"
#include"openMVG/sfm/sfm_data.hpp"
#include"patchHeader/grid/gridByVector.hpp"
#include"patchHeader/grid/bilinearInterpolate.hpp"
#include<memory>

//根据已有的网格指导生成新的网格高度
class HeightMakerByIndicate
{
public:
    //带有高度信息的网格单元
    using HeightGrid=GridInterface<double>;
    //点云的数据类型
    using PointCloud=openMVG::sfm::Landmarks;
    //网格单元的数据统计
    using CountGrid=GridByVector<unsigned>;

    //双线性插值的访问
    std::shared_ptr<BilinearInterpolate<double>> bilinearInterpolate_;

    HeightMakerByIndicate()
    {
        bilinearInterpolate_.reset(new BilinearInterpolate<double>());
    }

    //输入单个点，做计数上的维护
    //coord是点云的三维点坐标
    void addOnePointForAvg(HeightGrid& dstGrid,
                           CountGrid& count,double* coord
    )
    {
        //判断是否在范围内，不在范围的话就算了
        if(!dstGrid.isSpatialPointInRange(coord))
        {
            return;
        }
        //当前位置的计数
        auto& currCount=count.getDataBySpatialCoord(coord);
        //当前位置的计数数据
        auto& dstSum=dstGrid.getDataBySpatialCoord(coord);
        //如果是第1次计数就直接记录，不需要做加法
        if(currCount==0)
        {
            dstSum=coord[2];
        }
        else
        {
            //正常情况下是把传入的z值加到目标的sum上面
            dstSum+=coord[2];
        }
        //计数增加
        ++currCount;
    }

    //准备每个网格单元位置的平均高度
    void makeFinalAvgHeight(HeightGrid& dstGrid,HeightGrid& refGrid,
                            CountGrid& count)
    {
        //遍历每个网格单元
        int xy[2];
        //std::cout<<"******grid size*******"<<dstGrid.getDimSize(0)<<"******"<<dstGrid.getDimSize(1)<<std::endl;
        for(xy[1]=0;xy[1]<dstGrid.getDimSize(1);++xy[1])
        {
            for(xy[0]=0;xy[0]<dstGrid.getDimSize(0);++xy[0])
            {
                //获取当前位置的数据
                auto& currData=dstGrid.getDataAt(xy);
                //获取当前位置的计数
                auto& currCount=count.getDataAt(xy);
                //如果没有计数信息就从参考网格里面获取信息
                if(currCount==0)
                {
                    //把xy转换成空间坐标
                    double xyWorld[2];
                    xyWorld[0]=dstGrid.toSpatialDim(0,xy[0]);
                    xyWorld[1]=dstGrid.toSpatialDim(1,xy[1]);
                    //currData=refGrid.getDataBySpatialCoord(xyWorld);
                    currData=bilinearInterpolate_->getDataAt(refGrid,xyWorld);
                }
                else
                {
                    //在这个位置保存成平均的高度
                    currData/=(double)currCount;
                }
                //std::cout<<"x y: "<<xy[0]<<" "<<xy[1]<<" data: "<<currData<<" count: "<<currCount<<std::endl;
            }
        }
    }

    //根据高度获取每个位置的网格单元
    void makeHeightGrid(HeightGrid& dstGrid,HeightGrid& refGrid,
                        PointCloud& points)
    {
        //初始化每个位置的网格单元计数
        CountGrid count(dstGrid.getSpatialResolution(0),
                        dstGrid.getDimSize(0),dstGrid.getDimSize(1),dstGrid.getDimMin(0),
                        dstGrid.getDimMin(1));
        //遍历点云做高度计数
        for(auto& eachPoint : points)
        {
            //把点加进去用于计算平均高度
            addOnePointForAvg(dstGrid,count,eachPoint.second.X.data());
        }
        //根据计数结果生成最后的平均网格
        makeFinalAvgHeight(dstGrid,refGrid,count);
    }
};
