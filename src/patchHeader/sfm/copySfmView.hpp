#pragma once
#include"openMVG/sfm/sfm_data.hpp"

//复制sfm的view的操作，只有复制view的时候才需要这样
//因为存在指针深拷贝的问题
class CopySfmView
{
protected:
    using SfmData=openMVG::sfm::SfM_Data;
    //带有先验信息的view的数据类型
    using PriorView=openMVG::sfm::ViewPriors;

    //判断这个view是否需要被复制
    virtual bool judgeNeedCopy(unsigned idView,SfmData& sfm)
    {
        return true;
    }

    //把一个指定的view从一个sfm复制到另一个里面，复制的时候需要确保数据和粘连
    void copyOneView(unsigned idView,SfmData& srcSfm,
                     SfmData& dstSfm)
    {
        //目标的指针
        auto targetPtr=srcSfm.views.at(idView);
        //指针需要被复制到的位置
        auto& dstPtr=dstSfm.views[idView];
        //目前默认每个view都是有gps信息的，暂不考虑没有gps信息的情况
        dstPtr=std::shared_ptr<PriorView>(new PriorView());
        //把原本的指针也解析成prior指针
        auto srcPrior=dynamic_cast<PriorView*>(targetPtr.get());
        auto dstPrior=dynamic_cast<PriorView*>(dstPtr.get());
        //直接复制指向的内容 但最后是两个对象，也就是深拷贝
        *dstPrior=*srcPrior;
        //判断原本的view是否存在pose
        if(srcSfm.poses.count(srcPrior->id_pose)>0)
        {
            dstSfm.poses[srcPrior->id_pose]=srcSfm.poses.at(srcPrior->id_pose);
        }
    }
public:
    //复制sfm的view信息
    void copySfmView(SfmData& srcSfm,
                     SfmData& dstSfm)
    {
        //遍历每个view
        for(auto& eachView : srcSfm.views)
        {
            //判断这个view是否需要被删除
            if(this->judgeNeedCopy(eachView.first,srcSfm))
            {
                copyOneView(eachView.first,srcSfm,dstSfm);
            }
        }
    }
};
