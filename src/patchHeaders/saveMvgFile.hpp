#pragma once
#include"mvgHeader.hpp"
#include<string>
#include<iostream>
#include"fileException.hpp"

//保存openmvg的文件
template<class ThrowFunctor=FilePathError<>>//抛出异常所用的函数
class SaveMvgFile
{
public:
    //用于抛出异常的函数
    ThrowFunctor errorHandle;

    //需要被保存的路径
    std::string filePath_;

    //保存mvg的点云数据
    void operator()(SfM_Data& sfmData)
    {
        if(!openMVG::sfm::Save(sfmData,filePath_,openMVG::sfm::ESfM_Data::ALL))
        {
            std::cerr<<"cannot save "<<filePath_<<std::endl;
            throw -1;
        }
    }
};
