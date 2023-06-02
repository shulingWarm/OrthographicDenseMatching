#pragma once
#include"patchHeader/sfm/getInfoFromSfm.hpp"

class InvalidObsDeleter
{
public:
    //过滤view丢掉的obs
    void removeInvalidObs(openMVG::sfm::SfM_Data& sfm)
    {
        //删除非法的pose
        if(sfm.poses.count(openMVG::UndefinedIndexT))
        {
            sfm.poses.erase(openMVG::UndefinedIndexT);
        }
        //遍历每个点
        for(auto pointIter=sfm.structure.begin();pointIter!=sfm.structure.end();)
        {
            //当前的obs
            auto& obs=pointIter->second.obs;
            //遍历每个观察点
            for(auto iter=obs.begin();iter!=obs.end();)
            {
                //判断view是否存在
                if(sfm.views.count(iter->first)>0 && sfm.poses.count(
                            sfm.views.at(iter->first)->id_pose)>0)
                {
                    ++iter;
                }
                else
                {
                    iter=obs.erase(iter);
                }
            }
            //判断点观察者是否足够
            if(obs.size()<3)
            {
                pointIter=sfm.structure.erase(pointIter);
            }
            else
            {
                ++pointIter;
            }
        }
    }
};
