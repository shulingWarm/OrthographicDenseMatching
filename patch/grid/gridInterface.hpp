#pragma once

//抽象的网格单元
template<class Cell>
class GridInterface
{
public:
    //获取特定数据的接口
    //依次是xy
    virtual Cell& getDataAt(int* xy)=0;

    //获取指定维度的单元格分辨率，一般就是0或1
    virtual double getSpatialResolution(int dim)=0;

    //获取指定维度的最小值
    virtual double getDimMin(int dim)=0;

    //获取每个维度的离散长度
    virtual int getDimSize(int dim)=0;

    //用离散的方式获取数据
    virtual Cell& getDataBySpatialCoord(double* xy)
    {
        //转换成离散的形式
        int xyGrid[2];
        xyGrid[0]=toGridDim(0,xy[0]);
        xyGrid[1]=toGridDim(1,xy[1]);
        //用离散形式获取数据
        return getDataAt(xyGrid);
    }

    //指定某个维度，把连续值变成离散值
    int toGridDim(int dim,double x)
    {
        return (x-getDimMin(dim))/getSpatialResolution(dim);
    }

    //从离散域转成连续域
    double toSpatialDim(int dim,int x)
    {
        return (0.5+x)*getSpatialResolution(dim)+getDimMin(dim);
    }

    //判断一个点是否在范围内
    virtual double isSpatialPointInRange(double* xy)
    {
        for(int i=0;i<2;++i)
        {
            if(toGridDim(i,xy[i])>=getDimSize(dim) ||
                    toGridDim(i,xy[i])<0)
            {
                return false;
            }
        }
        return true;
    }

    //判断一个网格单元是否在范围内
    bool isGridDimInRange(int dim,int x)
    {
        return x>=0 && x<getDimSize(dim);
    }

    //判断某个空间维度是否在范围内
    bool isSpatialDimInRange(int dim,double x)
    {
        return isGridDimInRange(dim,toGridDim(dim,x));
    }

    //判断一个空间点是否在范围内
    bool isSpatialPointInRange(double* xy)
    {
        return isSpatialDimInRange(0,xy[0]) &&
                isSpatialDimInRange(1,xy[1]);
    }
};
