#pragma once
#include"patchHeader/process/domGenerate.hpp"
#include"patchHeader/process/orthographicDenseMatching.hpp"
#include"patchHeader/domGrid/domIO.hpp"
#include"patchHeader/initDOM/mvgInitDomByExistDomByCloud.hpp"
#include<string>
#include<memory>
#include"patchHeader/initDOM/addInitBackground.hpp"
#include"patchHeader/initDOM/addMaskToDom.hpp"

//dom生成过程，但是生成的过程中还会引入一个生成过的dom作为指导
class DomGenerateWithIndicate : public DomGenerateByManager
{
protected:
    using DomGrid=openMVG::sfm::DomInfo;

    //顶视稠密匹配的核心函数，独立成模块方便解耦合
    std::shared_ptr<OrthographicDenseMatching> matchingPtr_;
    //向一个已经初始化过的程序中加背景图
    std::shared_ptr<AddInitBackground> backgroundAddPtr_;
    //给DOM添加mask的操作
    std::shared_ptr<AddMaskToDOM> maskAdder_;

    //带有网格形式的点云初始化
    void initDomWithGrid(const std::string& gridPath,
                         DomManager& manager)
    {
        //先做DOM从已知点云的初始化
        matchingPtr_->beforeIterateProcess(manager);
        //判断是否有可用的网格数据
        if(gridPath.size()>0)
        {
            DomIO ioTool;
            DomGrid grid;
            ioTool.loadDom(grid,gridPath);
            MvgInitDomByExistByCloud initTool;
            initTool.initDom(grid,manager.domInfo_);
        }
    }

    //核心的生成过程，因为已经有了专门的初始化过程，因此这里不需要再初始化了
    void coreGenerateProcess(DomManager& manager) override
    {
        //直接开始迭代过程就可以了
        try
        {
            matchingPtr_->iterateProcess(manager);
        }
        catch(int errorFlag)
        {
            std::cerr<<errorFlag<<std::endl;
            throw errorFlag;
        }
    }

public:
    //构造函数，初始化稠密匹配器
    DomGenerateWithIndicate()
    {
        matchingPtr_.reset(new OrthographicDenseMatching());
        //初始化一个指针，用于向DOM里面加一个初始化的背景图
        backgroundAddPtr_.reset(new AddInitBackground());
        maskAdder_.reset(new AddMaskToDOM());
    }

    //生成dom的过程，但是允许在点云初始化之后载入一个辅助性的坐标
    void generateWithIndicate(
            const std::string& sfmPath,
             const std::string& outputPath,
             const std::string& rangeStr,
             double pixelLength, //空间分辨率
             const std::string& imageRootPath, //图片的根目录
             const std::string& gridPath, //已经生成过的DOM网格单元
              const std::string& backgroundPath="", //预添加的背景图路径
            const std::string& maskPath="" //对已经生成的场景添加的mask
            )
    {
        //读取dom管理器
        DomManager manager;
        readDomManager(sfmPath,manager);
        //记录空间分辨率
        setConfigInfo(manager,rangeStr,outputPath,pixelLength,imageRootPath);
        //借助网格单元辅助生成DOM
        initDomWithGrid(gridPath,manager);
        //初始化dom的背景图
        backgroundAddPtr_->addBackgroundToDom(manager.domInfo_,backgroundPath);
        //给dom添加mask
        maskAdder_->addMask(manager.domInfo_,maskPath);
        //生成dom的核心流程
        this->coreGenerateProcess(manager);
        //保存dom的结果
        saveDom(manager,outputPath);
    }
};
