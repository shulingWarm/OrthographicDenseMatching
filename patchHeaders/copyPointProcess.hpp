#pragma once
#include"loadMvgFile.hpp"
#include"saveMvgFile.hpp"
#include"copyPoint.hpp"
#include<string>

//完整的复制点云的流程
template<class SfmData,//sfm的文件
        class ReadFunctor,//载入点云函数
         class SaveFunctor,//保存点云的函数
         class CopyPointFunctor //基本的复制点云的函数
>
class CopyProcess
{
public:
    ReadFunctor read;
    SaveFunctor save;
    CopyPointFunctor copyHandle;

    //把其中一个点云合并到另一个点云
    void operator()(const std::string& dstFile,//点保存到这里
                    const std::string& refFile,//点从这里来
                    const std::string& savePath //最终结果保存目录
    )
    {
        //读取两个点云
        SfmData dstSfm,fromSfm;
        read.filePath_=dstFile;
        read(dstSfm);
        read.filePath_=refFile;
        read(fromSfm);
        //合并点云
        copyHandle(dstSfm,fromSfm.structure);
        //保存点云
        save.filePath_=savePath;
        save(dstSfm);
    }
};

//完整复制流程的具体化
template<class SfmData,class Landmarks=decltype(SfmData::structure)>
using MvgCopyProcess=CopyProcess<SfmData,LoadMvgFile,SaveMvgFile<>,
    CopySfmDataPoint<SfmData,Landmarks::mapped_type,Landmarks>>;
