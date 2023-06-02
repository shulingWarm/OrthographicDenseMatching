#pragma once
#include"patch/grid/gridInterface.hpp"
#include"openMVG/sfm/sfm_data.hpp"
#include"patch/sfm/getInfoFromSfm.hpp"
#include<memory>
#include"patch/imgSelect/occlusionJudger.hpp"

//通过sfm信息找一个指定的合理的view
//对于一个指定点， 投影到所有的图片上，获取一个积分最高且不遮挡的图片
class BestViewMakerBySfm
{
public:
    //DOM网格单元，但只是根据对应位置获取它的Z值
    using DomZMat=GridInterface<double>;
    using SfmData=openMVG::sfm::SfM_Data;
    using Vec3=Eigen::Vector3d;

    std::shared_ptr<SfmInterface> sfmInterface_;
    //用来判断是否存在遮挡的工具
    std::shared_ptr<OcclusionJudger> occlusionJudge_;

    BestViewMakerBySfm()
    {
        sfmInterface_.reset(new SfmInterface());
        occlusionJudge_.reset(new OcclusionJudger());
    }

    //获取点在某个相机上的投影分数
    double getProjectScoreByView(SfmData& sfm,
                                 Vec3 worldPoint,unsigned idView,double* zRange)
    {
        //把点放到最低值
        worldPoint[2]=zRange[0];
        //投影到图片上
        auto projLow=sfmInterface_->getProjectionByView(idView,sfm,worldPoint);
        //高位置的投影
        worldPoint[2]=zRange[1];
        auto projHigh=sfmInterface_->getProjectionByView(idView,sfm,worldPoint);
        //简单定义的投影分数，两个投影点的距离的相反数
        return 10000-(projLow-projHigh).norm();
    }

    //把点投影到每个图片上来获取最合适的view
    //成本很高，尽量少调用
    unsigned getBestView(DomZMat& domZ,SfmData& sfm,
                         Vec3 worldPoint,double* zRange //场景中z值的范围
    )
    {
        //目前已经获得的最佳分数
        double currMaxScore=-1e9;
        //目前选择的最佳图片
        unsigned bestView=0;
        //遍历每个view
        for(auto& eachView : sfm.views)
        {
            //当前点的投影位置
            auto projection=sfmInterface_->getProjectionByView(eachView.first,sfm,worldPoint);
            //判断投影位置是否在范围内
            if(sfmInterface_->isInImageRange(projection.data(),sfm,eachView.first))
            {
                //获取点在此处的投影分数
                double tempScore=getProjectScoreByView(sfm,worldPoint,eachView.first,zRange);
                //如果分数没超过就不用往下走了
                if(tempScore<=currMaxScore) continue;
                //获取相机光心
                auto viewCenter=sfmInterface_->getCenterByView(eachView.first,sfm);
                //如果超过了旧分数再检查一下是否有遮挡
                if(!occlusionJudge_->isOccluded(domZ,worldPoint.data(),viewCenter.data(),zRange[1]))
                {
                    //如果没有遮挡说明最佳的图片就可以更新了
                    currMaxScore=tempScore;
                    bestView=eachView.first;
                }
            }
        }
        //返回找到的最佳图片
        return bestView;
    }
};
