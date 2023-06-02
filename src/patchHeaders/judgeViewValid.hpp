#pragma once

//判断一个sfm_data里面是否存在某个view
template<class SfM_Data,class ViewIndex=unsigned>
class JudgeViewValid
{
public:
    //判断某个view是否可用
    bool operator()(const SfM_Data& sfmData,ViewIndex idView)
    {
        //判断view是否存在
        if(sfmData.views.count(idView)==0)
        {
            return false;
        }
        //判断相应的pose是否存在
        return sfmData.poses.count(sfmData.views[idView]->id_pose)>0;
    }
};
