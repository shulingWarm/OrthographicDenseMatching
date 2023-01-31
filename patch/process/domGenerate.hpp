#pragma once
#include"openMVG/sfm/sfm_data.hpp"



//从一个dom管理器里面做相关的dom生成管理
class DomGenerateByManager
{
protected:
    //所谓dom的管理器
    using DomManager=openMVG::sfm::SfM_Data;

    //读取sfm_data的数据
    void readDomManager(const std::string& sfmPath,
                        DomManager& dstManager)
    {
        if (!Load(dstManager, sfmPath, openMVG::sfm::ESfM_Data::ALL)) {
          std::cerr << std::endl
            << "The input SfM_Data file \""<< sfmPath << "\" cannot be read." << std::endl;
        }
    }

    //记录配置信息
    void setConfigInfo(DomManager& manager,
                       const std::string& rangeStr, //用来表示范围的字符串
                       const std::string& interOutputFolder, //中间过程的保存位置
                       double pixelLength,//空间分辨率
                       const std::string& imageRootPath //图片的根目录

    )
    {
        //记录空间分辨率
        manager.domInfo_.pixelLength_=pixelLength;
        //记录表示范围的字符串
        manager.cloudRangeInfo_=rangeStr;
        //记录中间过程保存的位置
        manager.midProceeFolder_=interOutputFolder;
        if(imageRootPath.size())
        {
            manager.s_root_path=imageRootPath;
        }
    }

    //生成dom的核心流程
    virtual void coreGenerateProcess(DomManager& manager)
    {
        manager.loadViewImgs();
        try {
            manager.denseDomLikeMvs();
        } catch (int errorFlag) {
            //这里只负责显示拿到的错误信息
            std::cerr<<errorFlag<<std::endl;
            throw errorFlag;
        }
    }

    //保存dom相关的结果
    void saveDom(DomManager& manager,
                 const std::string& outFolder)
    {
        //保存最后的DOM结果
        manager.saveDomResult(stlplus::create_filespec(outFolder, "domResult", ".bmp"));
        //释放图片内存
        manager.releaseImgs();
        //根据dom图里面的坐标重构点云
        manager.getZAsCloud();
        //保存点云
        openMVG::sfm::Save(manager,
          stlplus::create_filespec(outFolder, "cloud_and_poses", ".ply"),
          openMVG::sfm::ESfM_Data::STRUCTURE);
    }
public:
    //生成dom的流程
    void domGenerate(const std::string& sfmPath,
                     const std::string& outputPath,
                     const std::string& rangeStr,
                     double pixelLength, //空间分辨率
                     const std::string& imageRootPath //图片的根目录
    )
    {
        //读取dom管理器
        DomManager manager;
        readDomManager(sfmPath,manager);
        //记录空间分辨率
        setConfigInfo(manager,rangeStr,outputPath,pixelLength,imageRootPath);
        //生成dom的核心流程
        coreGenerateProcess(manager);
        //保存dom的结果
        saveDom(manager,outputPath);
    }
};
