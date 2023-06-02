#pragma once
#include"patchHeader/grid/gridHeight/getGridHeightByIndicate.hpp"
#include"patchHeader/grid/onlyValueGrid.hpp"
#include<algorithm>
#include<memory>
#include<cmath>

class GridHeightMaker
{
public:
    //带有高度信息的网格单元
    using HeightGrid=GridInterface<double>;
    //点云的数据类型
    using PointCloud=openMVG::sfm::Landmarks;
    //每一层的缩放率
    double scaleEach_=1.5;
    //利用已有数据生成新的高度网格的操作
    std::shared_ptr<HeightMakerByIndicate> indicateHeightMake_;

    GridHeightMaker()
    {
        indicateHeightMake_.reset(new HeightMakerByIndicate());
    }

    //临时计算点的平均高度
    double pointAvgHeight(PointCloud& points)
    {
        double sum=0;
        for(auto& eachPont : points)
        {
            sum+=eachPont.second.X[2];
        }
        return sum/(double)points.size();
    }

    //计算分层数
    unsigned getLayerNum(HeightGrid& dstGrid)
    {
        //网格的最小宽度
        double minDim=std::min<int>(dstGrid.getDimSize(0),dstGrid.getDimSize(1));
        //1.5的对数
        double logValue=std::log(minDim)/std::log(scaleEach_);
        //如果小于1就按0算
        if(logValue<1) return 0;
        //正常情况下取整数-1
        return (unsigned)(logValue-1);
    }

    //生成高度网格信息
   void makeHeight(HeightGrid& dstGrid,PointCloud& points)
   {
        //用点的平均高度做最顶层
       OnlyValueGrid<double> avgLayer;
       avgLayer.getDataAt(nullptr)=pointAvgHeight(points);
       //获取分层数
       unsigned layerNum=getLayerNum(dstGrid);
       //用来保存中间结果的两层
       std::shared_ptr<GridByVector<double>> lastLayer(new GridByVector<double>());
       std::shared_ptr<GridByVector<double>> nextLayer(new GridByVector<double>());
       //临时计算范围的最大值，中间会用到
       double xMax=dstGrid.getDimSize(0)*dstGrid.getSpatialResolution(0)+dstGrid.getDimMin(0);
       double yMax=dstGrid.getDimSize(1)*dstGrid.getSpatialResolution(1)+dstGrid.getDimMin(1);
       //依次分层获取每一层的平均高度
       for(int idLayer=layerNum;idLayer>=0;--idLayer)
       {
            //交换上次计算的两层
            std::swap(lastLayer,nextLayer);
            //上一层的指针
            HeightGrid* lastPtr=lastLayer.get();
            //如果是第1层则使用平均层
            if(idLayer==layerNum)
            {
                lastPtr=&avgLayer;
            }
            //下一层的指针
            HeightGrid* nextPtr=nextLayer.get();
            //如果是最后一层则直接使用目标结果
            if(idLayer==0)
            {
                nextPtr=&dstGrid;
            }
            else
            {
                //当前层的分辨率
                double currResolution=dstGrid.getSpatialResolution(0)*std::pow(scaleEach_,idLayer);
                //不使用最后一层的情况下需要先初始化当前层
                nextLayer->reset(currResolution,
                                 std::ceil((xMax-dstGrid.getDimMin(0))/currResolution),
                                 std::ceil((yMax-dstGrid.getDimMin(1))/currResolution),
                                 dstGrid.getDimMin(0),dstGrid.getDimMin(1));
            }
            //用上层初始化下层
            indicateHeightMake_->makeHeightGrid(*nextPtr,*lastPtr,points);
       }
   }
};
