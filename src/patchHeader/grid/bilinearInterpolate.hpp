#pragma once
#include"patchHeader/grid/gridInterface.hpp"
#include<cmath>

//对网格的双线性插值访问
//要求使用这种数据操作的cell是可以做四则运算的
template<class Cell>
class BilinearInterpolate
{
public:
    //网格单元的数据类型
    using GridType=GridInterface<Cell>;

    //获取xy的4个角点
    //cornerList分别对应x0 y0 x1 y1 x2 y2 x3 y3
    //顺序依次是x- y- x+ y- x- y+ x+ y+
    //xGrid和yGrid是离散化的单元格坐标,但还没有取整
    void getIntegerCorner(double xGrid,double yGrid,int* cornerList)
    {
        //依次记录8个离散值
        for(int idCorner=0;idCorner<4;++idCorner)
        {
            int* currCorner=cornerList+idCorner*2;
            //x的值
            if(idCorner%2)
            {
                currCorner[0]=std::ceil(xGrid);
            }
            else
            {
                currCorner[0]=std::floor(xGrid);
            }
            //y的值
            if(idCorner/2)
            {
                currCorner[1]=std::ceil(yGrid);
            }
            else
            {
                currCorner[1]=std::floor(yGrid);
            }
        }
    }

    //从4个角点坐标里面取4个坐标值
    void getDataByArray(int* corners,Cell* dstData,GridType& grid,unsigned dataNum=4)
    {
        //遍历取点
        for(unsigned idData=0;idData<dataNum;++idData)
        {
            //取数据
            dstData[idData]=grid.getDataAt(corners+idData*2);
        }
    }

    //两点之间指定比例的插值
    Cell interpolate(Cell minData,Cell maxData,double rate)
    {
        return minData+(maxData-minData)*rate;
    }

    //对4个数据做双线性插值
    //这里面传入的xRate,yRate默认是0~1的区间，其实就是一个比例
    Cell bilinearInterpolate(Cell* data4,double xRate,double yRate)
    {
        //x范围的两次插值
        Cell xMin=interpolate(data4[0],data4[1],xRate);
        Cell xMax=interpolate(data4[2],data4[3],xRate);
        //对y维度的插值
        return interpolate(xMin,xMax,yRate);
    }

    //对网格单元作双线性插值访问,不对传入的位置做范围内检查
    virtual Cell getDataAt(GridType& grid,double* xy)
    {
        //获取离散化后的浮点坐标
        double xGrid=(xy[0]-grid.getDimMin(0))/grid.getSpatialResolution(0);
        double yGrid=(xy[1]-grid.getDimMin(1))/grid.getSpatialResolution(1);
        //对于一组xy数据，获取它的4个角点
        int xyCorner[8];
        getIntegerCorner(xGrid,yGrid,xyCorner);
        //取对应位置的4个坐标值
        Cell cornerData[4];
        getDataByArray(xyCorner,cornerData,grid,4);
        //对4个数组做双线性插值
        return bilinearInterpolate(cornerData,xGrid-xyCorner[0],yGrid-xyCorner[1]);
    }
};

//对网格的双线性插值访问
//要求使用这种数据操作的cell是可以做四则运算的
template<class Cell>
class FakeInterpolate : public BilinearInterpolate<Cell>
{
public:
    //网格单元的数据类型
    using GridType=GridInterface<Cell>;

    //对网格单元作双线性插值访问,不对传入的位置做范围内检查
    virtual Cell getDataAt(GridType& grid,double* xy) override
    {
        //直接使用普通的数据访问，不考虑双线性插值
        return grid.getDataBySpatialCoord(xy);
    }
};
