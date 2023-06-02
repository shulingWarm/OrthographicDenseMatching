#pragma once
#include"patchHeader/grid/gridByVector.hpp"
#include"patchHeader/grid/gridIO.hpp"
#include<string>
#include<memory>

//对gridByVecotr的IO操作
template<class Cell>
class GridByVectorIO
{
public:
    //承载cell的网格
    using GridType=GridByVector<Cell>;
    //字节流的操作 对一个字节流保存或读取
    typedef void(*ByteFunctor)(std::fstream&,char*,unsigned);
    //保存网格数据的操作
    std::shared_ptr<GridIOTool<Cell>> saveData_;

    GridByVectorIO()
    {
        saveData_.reset(new GridIOTool<Cell>());
    }

    //保存字节流
    static void saveByte(std::fstream& fileHandle,char* data,unsigned byteLength)
    {
        fileHandle.write(data,byteLength);
    }

    //载入字节流
    static void loadByte(std::fstream& fileHandle,char* data,unsigned byteLength)
    {
        fileHandle.read(data,byteLength);
    }

    //保存或读取网格的基本属性，长度和分辨率之类的
    void operateGridProperty(GridType& grid,std::fstream& fileHandle,ByteFunctor byteFunc)
    {
        //保存或载入分辨率
        byteFunc(fileHandle,(char*)&grid.resolution_,sizeof(double));
        //保存xy的最小值
        byteFunc(fileHandle,(char*)grid.xyMin_.data(),sizeof(double)*2);
        //保存长度和宽度
        byteFunc(fileHandle,(char*)grid.xySize_.data(),sizeof(unsigned)*2);
    }

    //把网格保存在特定的路径中
    void saveGrid(GridType& grid,const std::string& filePath)
    {
        //打开输入流
        std::fstream fileHandle;
        fileHandle.open(filePath,std::ios::out|std::ios::binary);
        if(!fileHandle.is_open())
        {
            std::cerr<<"cannot save "<<filePath<<std::endl;
            return;
        }
        //在二进制流中保存网格的相关属性
        operateGridProperty(grid,fileHandle,saveByte);
        //保存二进制的数据
        saveData_->saveGrid(grid,fileHandle);
    }

    //从指定的路径里面保存网格
    void loadGrid(GridType& grid,const std::string& filePath)
    {
        //打开输入流
        std::fstream fileHandle;
        fileHandle.open(filePath,std::ios::in|std::ios::binary);
        if(!fileHandle.is_open())
        {
            std::cerr<<"cannot load "<<filePath<<std::endl;
            return;
        }
        //在二进制流中保存网格的相关属性
        operateGridProperty(grid,fileHandle,loadByte);
        //根据读到的属性初始化每个网格单元
        grid.reset(grid.resolution_,grid.xySize_[0],grid.xySize_[1],
                        grid.xyMin_[0],grid.xyMin_[1]);
        //保存二进制的数据
        saveData_->loadGrid(grid,fileHandle);
    }
};
