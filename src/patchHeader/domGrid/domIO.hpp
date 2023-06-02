#pragma once
#include"openMVG/sfm/sfm_data.hpp"
#include<string>
#include<fstream>


//二进制数据的存取
template<bool saveFlag=true>
class BinaryIO
{
protected:
    using FileHandle=std::fstream;
public:
    //对二进制数据的保存或读取
    template<class T>
    void operateData(FileHandle& stream,T& opData)
    {
        if(saveFlag)
        {
            stream.write((char*)&opData,sizeof(T));
        }
        else
        {
            stream.read((char*)&opData,sizeof(T));
        }
    }

    //对多个数据的操作
    template<class T>
    void operateManyData(FileHandle& stream,T* opData,unsigned length)
    {
        if(saveFlag)
        {
            stream.write((char*)opData,length*sizeof(T));
        }
        else
        {
            stream.read((char*)opData,length*sizeof(T));
        }
    }
};

//DOM网格单元的存储接口
class DomIO
{
protected:
    //DOM网格
    using DomGrid=openMVG::sfm::DomInfo;
    //DOM的网格单元
    using DomCell=openMVG::sfm::DomUnit;

    using FileHandle=std::fstream;

    //对基本信息的记录
    template<bool saveFlag>
    void operateBaseInfo(FileHandle& stream,DomGrid& grid)
    {
        //二进制数据的操作
        auto ioOperator=BinaryIO<saveFlag>();
        //记录列数
        ioOperator.operateData(stream,grid.domWidth_);
        //记录行数
        ioOperator.operateData(stream,grid.domHeight_);
        //记录间距
        ioOperator.operateData(stream,grid.pixelLength_);
        //记录点云范围
        ioOperator.operateManyData(stream,grid.cloudRange_.data(),2);
    }

    //单独记录一个网格单元的信息
    template<bool saveFlag>
    void operateEachCellInfo(FileHandle& stream,DomCell& cell)
    {
        //二进制数据的操作
        auto ioOperator=BinaryIO<saveFlag>();
        //记录z值
        ioOperator.operateData(stream,cell.z_);
        //记录分数
        ioOperator.operateData(stream,cell.nccScore_);
    }

    //处理每个网格单元的信息
    //对基本信息的记录
    template<bool saveFlag>
    void operateCellInfo(FileHandle& stream,DomGrid& grid)
    {
        //记录网格的每个cell
        for(unsigned x=0;x<grid.domWidth_;++x)
        {
            for(unsigned y=0;y<grid.domHeight_;++y)
            {
                //获取当前的cell
                auto& currCell=grid.getUnit(x,y);
                //记录每个cell的信息
                this->operateEachCellInfo<saveFlag>(stream,currCell);
            }
        }
    }

    //对dom数据的操作
    template<bool saveFlag>
    void operateData(FileHandle& stream,DomGrid& grid)
    {
          //记录基本信息的内容
        this->operateBaseInfo<saveFlag>(stream,grid);
        //对于读取的情况，需要先初始化一下DOM
        if(!saveFlag)
        {
            grid.initDom();
        }
        //记录每个网格单元的信息
        this->operateCellInfo<saveFlag>(stream,grid);
    }
public:
    //把一个DOM保存成指定的文件
    void saveDom(DomGrid& grid,const std::string& filePath)
    {
        //读取输出流
        FileHandle stream;
        stream.open(filePath,std::ios::out|std::ios::binary);
        operateData<true>(stream,grid);
    }

    //从指定文件中读取DOM网格数据
    void loadDom(DomGrid& grid,const std::string& filePath)
    {
        //准备输入流
        FileHandle stream;
        stream.open(filePath,std::ios::in|std::ios::binary);
        operateData<false>(stream,grid);
    }
};
