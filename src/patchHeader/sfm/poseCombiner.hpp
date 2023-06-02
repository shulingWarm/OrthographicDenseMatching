#pragma once
#include"openMVG/sfm/sfm_data.hpp"
#include"patchHeader/sfm/getInfoFromSfm.hpp"
#include<unordered_map>

//合并两个sfm的指定pose用的工具
class PoseCombiner
{
public:
    using SfmData=openMVG::sfm::SfM_Data;
    //每个view能看到的点数统计
    using ViewObsCount=std::unordered_map<unsigned,unsigned>;
    //上次处理的sfm的指针
    SfmData* sfm2Ptr_=nullptr;
    //上次处理的总的sfm的指针
    SfmData* sfm1Ptr_=nullptr;
    //用来存储每个sfm的每个view能看到的信息记录
    ViewObsCount count1;
    ViewObsCount count2;

    std::shared_ptr<SfmInterface> sfmInterface_;

    PoseCombiner()
    {
        sfmInterface_.reset(new SfmInterface());
    }

    //统计一个sfm里面每个view能看到的点
    //其实应该单独成立一个模块，但毕竟不是什么大函数
    void makeViewCount(SfmData& sfm,ViewObsCount& obsCount)
    {
        //遍历所有的点
        for(auto& eachPoint : sfm.structure)
        {
            //遍历所有的obs
            for(auto& eachObs : eachPoint.second.obs)
            {
                //判断是否有这个view
                if(obsCount.count(eachObs.first)==0)
                {
                    obsCount[eachObs.first]=1;
                }
                else
                {
                    //常规计数
                    ++obsCount[eachObs.first];
                }
            }
        }
    }

    //更新两个sfm的指针
    //如果操作的sfm换了需要及时更新
    void updateTargetSfm(SfmData& sfm1,SfmData& sfm2)
    {
        //检查sfm有没有被更新
        if(sfm1Ptr_!=&sfm1)
        {
            sfm1Ptr_=&sfm1;
            count1.clear();
            makeViewCount(sfm1,count1);
        }
        //检查sfm2有没有更新
        if(sfm2Ptr_!=&sfm2)
        {
            sfm2Ptr_=&sfm2;
            count2.clear();
            //重新计算sfm2的每个view能看到的点数
            makeViewCount(sfm2,count2);
        }
    }

    //从sfm2里面复制view
    void copyView(SfmData& sfm1,SfmData& sfm2,unsigned idView)
    {
        //拿对应的priorView
        openMVG::sfm::ViewPriors* priorPtr=dynamic_cast<openMVG::sfm::ViewPriors*>(sfm2.views.at(idView).get());
        //拷贝构造新的智能指针
        std::shared_ptr<openMVG::sfm::ViewPriors> newPtr;
        newPtr.reset(new openMVG::sfm::ViewPriors(*priorPtr));
        //记录到sfm1里面
        sfm1.views[idView]=newPtr;
    }

    //合并两个sfm的指定pose
    //这种设计要求的是同一对sfm需要被集中调用
    void copyPose(SfmData& sfm1,SfmData& sfm2,unsigned idView)
    {
        //更新目标的sfm
        updateTargetSfm(sfm1,sfm2);
        //pose的标号
        unsigned idPose=sfm2.views.at(idView)->id_pose;
        //pose不存在就直接跳过
        if(!sfmInterface_->isPoseExist(sfm2,idView))
        {
            return;
        }
        if(idPose!=idView)
        {
            std::cerr<<"warning: idPose: "<<idPose<<" != idView: "<<idView<<std::endl;
        }
        //需要确保sfm1里面有这个view
        if(sfm1.views.count(idView)==0)
        {
            copyView(sfm1,sfm2,idView);
        }
        //判断sfm1里面有没有这个view的pose记录
        if(count1.count(idView)==0)
        {
            count1[idView]=0;
        }
        //如果新的pose能看到更多的点就记录下来
        if(count2[idView]>count1[idView])
        {
            count1[idView]=count2[idView];
            sfm1.poses[idPose]=sfm2.poses[idPose];
            sfm1.views[idView]->id_pose=idPose;
        }
    }
};
