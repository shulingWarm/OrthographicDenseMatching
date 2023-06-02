#pragma once
#include"patch/grid/gridInterface.hpp"
#include"openMVG/sfm/sfm_data.hpp"
#include"patch/sfm/getInfoFromSfm.hpp"
#include<memory>
#include"patch/imgSelect/occlusionJudger.hpp"
#include"patch/imgSelect/getBestViewBySfm.hpp"

class BestViewMakerByKnownMat
{
public:
    //网格单元里面每个位置的最好的投影图片
    //键值存储的是一个图片标号
    using ViewMat=GridInterface<unsigned>;
    //DOM网格单元，但只是根据对应位置获取它的Z值
    using DomZMat=GridInterface<double>;
    using SfmData=openMVG::sfm::SfM_Data;
    using Vec3=Eigen::Vector3d;

    std::shared_ptr<SfmInterface> sfmInterface_;
    //用来判断是否存在遮挡的工具
    std::shared_ptr<OcclusionJudger> occlusionJudge_;
    //通过sfm信息获取到的最佳投影图片
    std::shared_ptr<BestViewMakerBySfm> getBestViewBySfm_;

    BestViewMakerByKnownMat()
    {
        sfmInterface_.reset(new SfmInterface());
        occlusionJudge_.reset(new OcclusionJudger());
        getBestViewBySfm_.reset(new BestViewMakerBySfm());
    }

    //参考的最佳图片是否可用
    bool isRefBestViewUsable(ViewMat& refViewMat,
                             DomZMat& domZ,//可获取每个位置Z值的方法
                             SfmData& sfm,//sfm的数据
                             Vec3 worldPoint,
                             double maxZ)
    {
        if(refViewMat.isSpatialPointInRange(xy))
        {
            //获取这个位置的最佳参考图片
            unsigned refBest=refViewMat.getDataBySpatialCoord(xy);
            //获取此处的投影点位置
            auto projection=sfmInterface_->getProjectionByView(refBest,sfm,worldPoint);
            //判断一下点是否在投影范围内
            if(sfmInterface_->isInImageRange(projection.data(),sfm,refBest))
            {
                //相机光心的三维坐标
                auto cameraCenter=sfmInterface_->getCenterByView(refBest,sfm);
                //判断从当前点到目标相机之间是否存在遮挡
                if(!occlusionJudge_->isOccluded(domZ,worldPoint.data(),
                                                cameraCenter.data(),maxZ))
                {
                    return true;
                }
            }
        }
        return false;
    }

    //根据指定的位置，结合可参考的数据，获取最佳的观测图片
    unsigned getBestViewByKnownMat(ViewMat& refViewMat,
                                   DomZMat& domZ,//可获取每个位置Z值的方法
                                   SfmData& sfm,//sfm的数据
                                   double* xy,
                                   double* zRange
    )
    {
        //计算点的三维坐标
        Vec3 worldPoint;
        worldPoint<<xy[0],xy[1],domZ.getDataBySpatialCoord(xy);
        //判断参考的view是否可用
        if(isRefBestViewUsable(refViewMat,domZ,sfm,worldPoint,zRange[1]))
        {
            return refViewMat.getDataBySpatialCoord(xy);
        }
        //重新从sfm里面所有的view里找出来一个最合适的view
        return getBestViewBySfm_->getBestView(domZ,sfm,worldPoint,zRange);
    }

    //根据操作生成新的view
    void makeNewBestView(ViewMat& dstView,
                         ViewMat& refView,DomZMat& domZ,
                         SfmData& sfm,
                         double* zRange //两个数字，依次是z_min 和z_max
    )
    {
        //遍历目标矩阵里面的每个位置
        for(int y=0;y<dstView.getDimSize(1);++y)
        {
            for(int x=0;x<dstView.getDimSize(0);++x)
            {
                //当前位置的坐标
                double xy[2]={dstView.toSpatialDim(0,x),
                             dstView.toSpatialDim(1,y)};
                int gridXy[]={x,y};
                //记录算出来的最佳结果
                dstView.getDataAt(gridXy)=getBestViewByKnownMat(refView,domZ,sfm,xy,zRange);
            }
        }
    }
};
