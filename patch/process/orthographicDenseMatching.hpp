#pragma once
#include"openMVG/sfm/sfm_data.hpp"

//正射的密集匹配过程
class OrthographicDenseMatching
{
public:
    using DomManager=openMVG::sfm::SfM_Data;

    //循环前的准备
    void beforeIterateProcess(DomManager& manager)
    {
        //查找关注位置的点云
          manager.focusCloudPt();
          //判断是否需要删除中间位置
#ifdef USE_ONLY_MID_CAM
          manager.midCameraFilter();
#endif
          //先把已经有的稀疏点云的信息部署在dom图上
          manager.drawSparsPtAtDom();
          //删除没有载入过的相机光心
          openMVG::sfm::DeleteNonLoadedViews nonLoadDeleter;
          nonLoadDeleter(manager.views);
          //清空点云
          manager.structure.clear();
    }

    //循环迭代的过程
    void iterateProcess(DomManager& manager)
    {
        //记录开始的时间
        time_t bgTime=time(nullptr);
        //是否存在新的可优化的点
        bool haveUnRefine=true;
        //遍历掩膜步长的每一步
#ifdef FIX_MASK_TIMES//使用固定的掩膜大小循环若干次
        for(int iterTimes=0;iterTimes<FIX_MASK_TIMES;++iterTimes)
#else
        for(int maskSize=MASK_SIZE;maskSize>=3;maskSize-=2)
#endif
        {
            //判断是否需要每次都保存运行的dom结果
#ifdef SAVE_EACH_TIME
            //保存当前循环次数的快速dom结果
            if(NEED_SAVE(iterTimes))
            {
                if(haveUnRefine)//必须是有新的优化成果再保存
                {
                    //计算时间差
                    unsigned int seconds=time(nullptr)-bgTime;
                    //转换为需要记录的时间
                    std::string fileName=std::to_string(seconds/60)+"-"+std::to_string(seconds%60);
                   manager.saveDomResult(manager.midProceeFolder_+"/"+fileName+".bmp");
                    //临时保存点云，正常运行的时候不需要动这个地方
                    //saveDomCloud(midProceeFolder_+"/"+fileName+".ply");
                }
            }
#endif
            //当前循环的面片大小
            int maskSize=MASK_SIZE(iterTimes);
            //更新当前周期是否需要考虑遮挡问题
            manager.isUsingOccludeDetection_=isIterateDetectOcclude(iterTimes);
            std::cout<<iterTimes<<std::endl;
            //初始化是否有可优化的点
            haveUnRefine=false;
            for(int currStep=0;currStep<maskSize;currStep+=MASK_STEP)
            {
                 //根据目前的步长优化dom图的每个像素
#ifdef FIX_MASK_TIMES//同大小窗口固定迭代若干次
                //临时的结果
                bool resultFlag=manager.optimizeDomPixel(currStep,false,maskSize,1,THRE_ID(iterTimes),ALLOW_MIN_Z(iterTimes));
                //判断是否有未优化的点
                if(resultFlag) haveUnRefine=true;
#endif
            }
            //判断是否需要更新dom像素的可用次数
            if(NEED_INIT_REF(iterTimes))
            {
                //manager.domInfo_.initRefTime();
                manager.domInfo_.initRefTimeForBlackCell();
            }
            //如果没有待优化的点了结束循环
            if(iterTimes>CONSIDER_BREAK_TIME && haveUnRefine==false)
                break;
            //判断是否强制结束
            //if(forceEndJudge()) break;
        }
    }

    //正射匹配的过程
    void denseMatching(DomManager& manager)
    {
          //迭代前的准备
            beforeIterateProcess(manager);
            //正式的迭代过程
            iterateProcess(manager);
    }
};
