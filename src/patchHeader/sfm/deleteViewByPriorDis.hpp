#pragma once
#include"patchHeader/sfm/getInfoFromSfm.hpp"
#include<memory>

//通过pose到prior的距离删除对应的pose
class PoseDeleter
{
public:
    using Mat3=Eigen::Matrix3d;
    using Vec3=Eigen::Vector3d;
    using SfmData=openMVG::sfm::SfM_Data;

    std::shared_ptr<SfmInterface> sfmInterface_;
    //删除点的阈值
    double deleteThreshold_=2;

    PoseDeleter()
    {
        sfmInterface_.reset(new SfmInterface());
    }

    //判断这个view的pose是否需要删除
    bool needDelete(SfmData& sfm,unsigned idView)
    {
        //判断是否存在pose,没有的话就不必了
        if(!sfmInterface_->isPoseExist(sfm,idView))
        {
            return false;
        }
        //获取光心
        Vec3 poseCenter=sfmInterface_->getPoseCenterByView(idView,sfm);
        //获取先验信息的光心
        Vec3 priorCenter=sfmInterface_->getPriorCenterByView(idView,sfm);
        //光心的距离
        double dis=(poseCenter-priorCenter).norm();
        //大于阈值的情况需要删除
        return dis>deleteThreshold_;
    }

    //删除这个view的pose
    void deleteViewPose(SfmData& sfm,unsigned idView)
    {
        //目标的view
        auto& view=sfm.views.at(idView);
        //删除pose
        sfm.poses.erase(view->id_pose);
        //把view的id_pose设置成初始值
        view->id_pose=openMVG::UndefinedIndexT;
    }

    void deletePose(SfmData& sfm)
    {
        //遍历每个view,删除pose
        for(auto& eachView : sfm.views)
        {
            //判断这个view的pose是否需要删除
            if(needDelete(sfm,eachView.first))
            {
                deleteViewPose(sfm,eachView.first);
            }
        }
    }
};
