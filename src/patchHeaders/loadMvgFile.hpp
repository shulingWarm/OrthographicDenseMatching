#pragma once
#include<string>
#include"mvgHeader.hpp"


//用于载入openmvg的文件
class LoadMvgFile
{
public:
    //文件的读取路径,需要在调用函数之前修改
    std::string filePath_;

    //需要把结果存储在dstSfm里面
    void operator()(SfM_Data& dstSfm)
    {
        openMVG::sfm::Load(dstSfm,filePath_,openMVG::sfm::ESfM_Data::ALL);
    }
};
