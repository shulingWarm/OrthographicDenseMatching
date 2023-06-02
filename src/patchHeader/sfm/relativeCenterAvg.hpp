#pragma once
#include"patchHeader/sfm/getInfoFromSfm.hpp"
#include<unordered_map>

//估计一个sfm的相对光心
class RelativeCenterMaker
{
public:
    using Mat3=Eigen::Matrix3d;
    using Vec3=Eigen::Vector3d;
    using SfmData=openMVG::sfm::SfM_Data;

    std::shared_ptr<SfmInterface> sfmInterface_;

    RelativeCenterMaker()
    {
        sfmInterface_.reset(new SfmInterface());
    }

    //相对光心的计数器
    class RelativeCenterCount
    {
    public:
        Vec3 centerSum_;
        unsigned count_=0;

        //向里面添加新的相对光心
        void add(Vec3 newCenter)
        {
            //如果还是空的，就先初始化
            if(count_==0)
            {
                centerSum_<<0,0,0;
            }
            //叠加光心
            for(int i=0;i<3;++i) centerSum_(i)+=newCenter(i);
            ++count_;
        }

        //获取平均光心
        Vec3 getAvg()
        {
            return centerSum_/count_;
        }
    };

    //边相机到中相机的相对光心
    using RelativeCenterMap=std::unordered_map<unsigned,Vec3>;
    //计数器的列表
    using CountMap=std::unordered_map<unsigned,RelativeCenterCount>;

    //计算view2的光心在view1的坐标系下的光心位置
    Vec3 getCenter2inCenter1Axis(unsigned idView1,unsigned idView2,
                                 SfmData& sfm)
    {
        Vec3 center1=sfmInterface_->getPriorCenterByView(idView1,sfm);
        Vec3 center2=sfmInterface_->getPriorCenterByView(idView2,sfm);
        //1相机的光心旋转，不做检查，必须有
        Mat3 rot1=sfmInterface_->getRotByView(idView1,sfm);
        //返回2在1坐标系下的位置
        return rot1*(center2-center1);
    }

    //用一个输入的中相机的view来维护它的每个边相机的相对光心
    void countRelaiveCenter(SfmData& sfm,unsigned idView,CountMap& dstMap)
    {
        //传入的相机是中相机
        if(sfmInterface_->getDirByView(idView)!=2)
        {
            std::cerr<<"invalid direction "<<std::endl;
            return;
        }
        //遍历其它的方向
        for(unsigned idDir=1;idDir<=5;++idDir)
        {
            //中相机不处理
            if(idDir==2) continue;
            //判断有没有这个view
            unsigned tempView=sfmInterface_->getOtherDirByIdView(idView,idDir);
            if(sfm.views.count(tempView))
            {
                //计数器
                RelativeCenterCount& counter=dstMap[idDir];
                //计算这两个view的相对光心
                counter.add(getCenter2inCenter1Axis(idView,tempView,sfm));
            }
        }
    }

    //输入sfm,获取每个边相机到中相机的相对光心
    void getRelativeCenters(SfmData& sfm,RelativeCenterMap& dstCenters)
    {
        //新建计数器的列表
        CountMap tempMap;
        //遍历每个有pose的view
        for(auto& eachView : sfm.views)
        {
            //只处理中相机
            if(sfmInterface_->getDirByView(eachView.first)!=2) continue;
            //判断中相机的pose是否存在
            if(sfmInterface_->isPoseExist(sfm,eachView.first))
            {
                //用这个中相机维护计数器
                countRelaiveCenter(sfm,eachView.first,tempMap);
            }
        }
        //遍历每个count的列表，从里面获取平均的相对光心
        for(auto& eachCount : tempMap)
        {
            dstCenters[eachCount.first]=eachCount.second.getAvg();
        }
    }
};
