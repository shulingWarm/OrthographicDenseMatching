#pragma once
#include<fstream>
#include"patchHeader/grid/gridInterface.hpp"

//网格单元数据的IO操作接口
template<class Cell>
class GridIOTool
{
public:
    //被操作的目标网格数据类型
    using GridType=GridInterface<Cell>;

    //对网格单元的操作函数
    typedef void(*FileOp)(Cell&,std::fstream&);

    //向数据流里面保存网格单元的操作
    static void saveOperate(Cell& cell,std::fstream& fileHandle)
    {
        //写入数据
        fileHandle.write((char*)&cell,sizeof(Cell));
    }

    //从数据流里面载入数据的操作
    static void loadOperate(Cell& cell,std::fstream& fileHandle)
    {
        //写入数据
        fileHandle.read((char*)&cell,sizeof(Cell));
    }

    //遍历网格单元并对数据做保存或载入
    //opFunction是一个函数指针
    void operateData(GridType& grid,std::fstream& fileHandle,FileOp opFunction)
    {
        //遍历网格里面的每个数据
        int xy[2];
        for(xy[0]=0;xy[0]<grid.getDimSize(0);++xy[0])
        {
            for(xy[1]=0;xy[1]<grid.getDimSize(1);++xy[1])
            {
                //操作当前位置的数据
                opFunction(grid.getDataAt(xy),fileHandle);
            }
        }
    }

    //保存网格
    //fileHandle是一个已经打开了的二进制文件
    void saveGrid(GridType& grid,std::fstream& fileHandle)
    {
        //保存数据
        operateData(grid,fileHandle,saveOperate);
    }

    //载入网格的操作
    void loadGrid(GridType& grid,std::fstream& fileHandle)
    {
        //保存数据
        operateData(grid,fileHandle,loadOperate);
    }
};
