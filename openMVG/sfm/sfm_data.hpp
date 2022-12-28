// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_DATA_HPP
#define OPENMVG_SFM_SFM_DATA_HPP

#include <string>

#include"delaunator.hpp"
#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/sfm/sfm_landmark.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/sfm/sfm_view_priors.hpp"
#include "openMVG/types.hpp"
#include"pointcloudMatch.hpp"
#include<opencv2/opencv.hpp>
#include"triInterface.h"
#include<ctime>
#include"openMVG/sfm/sfm_data_io.hpp"


namespace openMVG {
namespace sfm {

/// Define a collection of Pose (indexed by View::id_pose)
using Poses = Hash_Map<IndexT, geometry::Pose3>;

/// Define a collection of IntrinsicParameter (indexed by View::id_intrinsic)
using Intrinsics = Hash_Map<IndexT, std::shared_ptr<cameras::IntrinsicBase>>;

#include"deleteViews.hpp"
#include"isCvViewEmpty.hpp"
#include"saveViewImgId.hpp"
#include"iterateAndSave.hpp"
#include"iterateView.hpp"
//2022-10-21 删除没有载入过相机的views的函数
using DeleteNonLoadedViews=DeleteViews<Views,
    isViewImageLoad<Views::value_type,true>>;

//判断一个数是否在两个数之间[left,right),左闭右开区间
template <typename T>
bool judgeInside(T midValue, T leftValue, T rightValue)
{
    return midValue<rightValue && midValue>=leftValue;
}

//权值使用的数据类型
typedef double WeightT;
//做德劳内三角化的时候用到的点坐标列表,里面存储的东西分别是x1,y1,x2,y2...
typedef std::vector<double> DelaunayPtList;
typedef std::map<IndexT,int> StaticInfo;//统计信息

//做快速DOM时候的一个信息单元
//必须保留一个无参数的构造函数，std::vector里面有用
class DomUnit
{
public:
    //总共的权值
    WeightT weight_=0;

    //dom在当前位置的z值
    double z_=0;

    //这个点是否来自稀疏点云
    bool isPriorPoint_=false;

    //该信息可以被参考的次数
    RefType referencedTime_=REF_TIME_INIT;

    //该像素位置的ncc分数
    double nccScore_=-2;

    //是否为一个正常范围之外的点，当没有三角面片涉及到这个位置的点的时候
    //它会被设置为true
    bool isOutRangePt_=true;

    //像素的法向量
    PixelDir pixDir_;

    //复制信息
    virtual void copyTo(DomUnit& otherUnit) const
    {
        otherUnit.weight_=weight_;
        otherUnit.z_=z_;
    }

    //设置这个网格单元是否可用
    void setValidFlag(bool flag)
    {
        //如果要设置成可用，就设置成正常的原始分数
        //不方便再设置新的变量了，这里面每加一个新的变量就会大量消耗内存
        if(flag==false)
        {
            nccScore_=100000;
            referencedTime_=0;//表示不可被参考
        }
    }

    //判断一个网格单元是否可用
    bool isValid()
    {
        return nccScore_<1000;
    }

    //往里面添加z值
    //如果添加成功返回true
    bool addZ(double zValue)
    {
        //判断当前位置的点是否为场景之外的点
        //能被调用这个函数的点肯定不是场景之外的点
        //这里判断一下是看一下这里是不是已经有了一个z值了
        if(isOutRangePt_==false)
        {
            //判断z值,仅保留比较大的z值，因为这里做的是俯视视角的dom
            if(z_>zValue) return false;
        }
        //记录z值
        z_=zValue;
        //标记场景内的点
        isOutRangePt_=false;
        //正常返回
        return true;
    }
};

//完全体的每个dom像素的信息，里面包括计算需要用到的所有信息
class DomUnitAll : public DomUnit
{
public:
    TempPtBase scenePt_;//对应于场景中的xy坐标
    //对应于dom图片中的像素坐标
    int pixelX_;
    int pixelY_;
    //当前位置对应的颜色指针
    cv::Vec3b* corresColor_;
    //存储在dom图里面的真正的数据
    //需要使用先验信息的情况下可能会用到这个功能
    //这里算是一个很不好看的补丁
    DomUnit* srcUnitPtr_=nullptr;

    //保留外来的颜色
    void saveColor(const UserColor& otherColor)
    {
        corresColor_->operator[](0)=otherColor.getB();
        corresColor_->operator[](1)=otherColor.getG();
        corresColor_->operator[](2)=otherColor.getR();
        //std::cout<<*corresColor_<<std::endl;
    }

    //从一个其它的对象里面复制信息
    //这里是用来从dom的class里面保存信息的，重要的是要保存一个指针
    void copyFrom(DomUnit& srcUnit)
    {
        //记录指针
        srcUnitPtr_=&srcUnit;
        //记录权值和z坐标
        z_=srcUnit.z_;
        weight_=srcUnit.weight_;
    }
};

//dom里面每个像素信息的二维列表，和DOM输出图一一对应
//虽然是二维的列表，但实际上是按照一维来存储的
typedef std::vector<DomUnit> UnitMap;
//用于表示点云范围的数据类型,左右上下低高
typedef std::array<double,6> RangeT;

//做快速DOM的时候用的信息
class DomInfo
{
protected:
    //每个像素对应位置的信息
    UnitMap unitMap_;

public:
    //dom的处理结果
    cv::Mat domResult_;
    //分辨率与点云数量的比值
    double resolutionRate_=100.0;
    //输出图片的最大分辨率
    long int maxResolution_=24000000;
    //指定的每个像素对应的长度，如果是负数的话那就参考分辨率的比值
    double pixelLength_=-1.f;
    //dom图片的宽和高
    std::size_t domWidth_;
    std::size_t domHeight_;
    //点云的坐标开始位置，左上角,分别表示左和上
    std::array<double,2> cloudRange_;
    //点云的z坐标的范围，分别表示低和高
    std::array<double,2> zRange_;
    //记录区分后的每个高度阶层的高度值
    double zLevelHeight_=0.f;
    //每个高度的分布，也就是记录一个每个高度分别有多少个点
    std::map<int,uint> zMap_;
    //是否使用权值,改成false的时候，权值仅仅会被比较，不会被叠加,权值高的被保留
    bool useWeight_=true;

    //初始化每个位置的先验点可用次数
    void initRefTime()
    {
        //遍历unitmap里面的每个位置
#ifdef USE_OMP
#pragma omp parallel for
#endif
        for(unsigned int unitCount=0;unitCount<unitMap_.size();++unitCount)
        {
            //更新当前位置的值
            unitMap_[unitCount].referencedTime_=REF_TIME_INIT;
        }
    }

    //恢复某个范围内的可参考次数
    void restoreRefTimeInRange(int xMin,int xMax,int yMin,int yMax)
    {
        for(int x=xMin;x<=xMax;++x)
        {
            for(int y=yMin;y<=yMax;++y)
            {
                auto& currUnit=getUnit(x,y);
                currUnit.referencedTime_=REF_TIME_INIT;
            }
        }
    }

    //把这个网格单元周围的单元恢复
    void restoreSourndCellRefTime(unsigned idRow,unsigned idCol)
    {
        //范围矫正
        int xMin=idCol-1;if(xMin<0) xMin=0;
        int xMax=idCol+1;if(xMax>=domWidth_) xMax=domWidth_-1;
        int yMin=idRow-1;if(yMin<0) yMin=0;
        int yMax=idRow+1;if(yMax>=domHeight_) yMax=domHeight_-1;
        restoreRefTimeInRange(xMin,xMax,yMin,yMax);
    }

    //只允许未优化过的网格单元被优化的策略
    void initRefTimeForBlackCell()
    {
        //遍历每一行和每一列
        for(unsigned idRow=0;idRow<domHeight_;++idRow)
        {
            for(unsigned idCol=0;idCol<domWidth_;++idCol)
            {
                //目前的单元格
                auto& currUnit=getUnit(idCol,idRow);
                //如果还没有被优化过
                if(currUnit.nccScore_<0)
                {
                    restoreSourndCellRefTime(idRow,idCol);
                }
            }
        }
    }

    //dom图片实际的分辨率
    double getResolution() const
    {
        return domWidth_*domHeight_;
    }

    //获取像素的长度
    double getPixelLength() const
    {
        return pixelLength_;
    }

    //根据点云范围，确定最终的分辨率,需要确定下来图像像素的长宽，每个像素单元对应的实际长度
    void makeResolution(double leftRange,double rightRange,
                        double upRange,double downRange,int pointNum)
    {
        //点云范围的宽和高
        double cloudWidth=rightRange-leftRange;
        double cloudHeight=downRange-upRange;
        //如果没有初始的像素长度
        if(pixelLength_<0)
        {
            //计算一个临时的分辨率
            long int tempResolution=pointNum*resolutionRate_;
            //判断是否超过了最大允许的分辨率
            tempResolution=tempResolution>maxResolution_?maxResolution_:tempResolution;
            //根据最终的分辨率确定每个像素单元的面积
            double pixelArea=cloudWidth*cloudHeight/tempResolution;
            //由此计算每个像素单元的边长
            pixelLength_=std::sqrt(pixelArea);
        }
        //计算最后的DOM图片的宽和高
        domWidth_=cloudWidth/pixelLength_;
        domHeight_=cloudHeight/pixelLength_;
        //记录点云范围的最左边和最上边
        cloudRange_[0]=leftRange;
        cloudRange_[1]=upRange;
    }

    //获取点云的范围
    static RangeT getCloudPointRange(Landmarks& cloudPoints)
    {
        //初始化用于返回的点云范围
        RangeT domRange;
        //范围已经被初始化过的标志
        bool initFlag=false;
        //遍历所有的点云
        for(auto& eachPoint : cloudPoints)
        {
            //取出当前位置的点坐标
            const double* const thisPoint=eachPoint.second.X.data();
            //判断是否已经初始化
            if(initFlag)
            {
                domRange[0]=thisPoint[0]<domRange[0]?thisPoint[0]:domRange[0];//左
                domRange[1]=thisPoint[0]>domRange[1]?thisPoint[0]:domRange[1];//右
                domRange[2]=thisPoint[1]<domRange[2]?thisPoint[1]:domRange[2];//上
                domRange[3]=thisPoint[1]>domRange[3]?thisPoint[1]:domRange[3];//下
                domRange[4]=thisPoint[2]<domRange[4]?thisPoint[2]:domRange[4];//低
                domRange[5]=thisPoint[2]>domRange[5]?thisPoint[2]:domRange[5];//高
            }
            else //直接把点赋值给范围
            {
                domRange[0]=thisPoint[0];
                domRange[1]=thisPoint[0];
                domRange[2]=thisPoint[1];
                domRange[3]=thisPoint[1];
                domRange[4]=thisPoint[2];
                domRange[5]=thisPoint[2];
                //设定已经初始化过的标志
                initFlag=true;
            }
        }
        //返回最后的范围
        return domRange;
    }

    //计算传入的z值所属的阶层
    int getZLevel(double zValue) const
    {
        //判断z的阶层高度是否初始化过
        if(zLevelHeight_==0) return ERROR_NO_INIT_ZRANGE;
        return (zValue-zRange_[0])/zLevelHeight_;
    }

    //获取指定level的底部z值
    //例如输入1的时候返回的是第1层的下限，也就是第2层的上限
    double getLevelHeight(int level)
    {
        //判断z的阶层高度是否初始化过
        if(zLevelHeight_==0) return ERROR_NO_INIT_ZRANGE;
        return zRange_[0]+zLevelHeight_*level;
    }

    //传入一堆点云，这个点云必须是快速DOM的时候用的点云。
    void initDom(Landmarks &cloudPoints,bool initFlag=true //是否需要初始化分辨率
    )
    {
        //初始化点云的范围,左右上下
        RangeT domRange=getCloudPointRange(cloudPoints);
        //记录最高和最低
        zRange_[0]=domRange[4];
        zRange_[1]=domRange[5];
        //输出点云的范围
        std::cout<<"left:"<<domRange[0]<<std::endl;
        std::cout<<"right:"<<domRange[1]<<std::endl;
        std::cout<<"up:"<<domRange[2]<<std::endl;
        std::cout<<"down:"<<domRange[3]<<std::endl;
        std::cout<<"low:"<<domRange[4]<<std::endl;
        std::cout<<"high:"<<domRange[5]<<std::endl;
        //记录z值每个高度的阶层的高度
        zLevelHeight_=(zRange_[1]-zRange_[0])*STEP_L1;
        //记录点云分辨率的相关数据
        if(initFlag)
        {
            makeResolution(domRange[0],domRange[1],domRange[2],domRange[3],
                    cloudPoints.size());
            initDom();
        }
    }

    //根据目前的网格大小，重新初始化网格单元和图片信息
    void initDom()
    {
        //初始化最后的输出图片
        domResult_.create(domHeight_,domWidth_,CV_8UC3);
        //初始化每个位置的信息表
        unitMap_.resize(domHeight_*domWidth_);
        //默认都是黑色的
        domResult_.setTo(cv::Scalar(0,0,0));
    }

    //初始化每个高度的分布
    void initZMap()
    {
        //遍历所有的dom信息
        for(UnitMap::iterator iter=unitMap_.begin();
            iter!=unitMap_.end();++iter)
        {
            //判断当前位置是否为先验点
            if(iter->isPriorPoint_)
            {
                //计算当前点所在的阶层
                int tempStepCount=getZLevel(iter->z_);
                //判断当前位置的阶层是否存在
                if(zMap_.count(tempStepCount))
                {
                    //直接计数增加
                    zMap_[tempStepCount]++;
                }
                else
                {
                    //把对应的数值初始化成1
                    zMap_[tempStepCount]=1;
                }
            }
        }
    }

    //判断传入的坐标是否在一个合理的范围内
    bool judgeInvalid(int pixelX,int pixelY) const
    {
        return !judgeInside<int>(pixelX,0,domWidth_) ||
                !judgeInside<int>(pixelY,0,domHeight_);
    }

    //把离散点转换成连续点
    //dim为0的时候表示转换x的维度
    int convertConti2Pixel(double x,unsigned dim=0)
    {
        return (x-cloudRange_[dim])/pixelLength_;
    }

    //判断传入的连续坐标是否在范围内
    bool judgeContiPointOutRange(double x,double y)
    {
        //把坐标转换成离散形式
        int pixelX=convertConti2Pixel(x);
        int pixelY=convertConti2Pixel(y,1);
        return judgeInvalid(pixelX,pixelY);
    }

    //把一个二维的坐标转换成unit的顺序存储单元的索引下标
    //不做任何错误校验
    IndexT convertLocalToIdLabal(int pixelX,int pixelY)
    {
        return pixelY*domWidth_+pixelX;
    }

    //获取一个确定位置的像素单元,自己负责访问范围的合理性
    DomUnit& getUnit(int pixelX,int pixelY)
    {
        //获取unit里面的索引位置
        IndexT idxLocal=pixelY*domWidth_+pixelX;
        //判断是否超过了索引长度
        if(idxLocal>=unitMap_.size()) throw ERROR_DOMINFO_OUTRANGE;
        return unitMap_[idxLocal];
    }

    //获取完全体的像素单元
    bool getUnit(int pixelX,int pixelY,DomUnitAll& dstUnit)
    {
        //如果不在一个合理的范围内就不做处理
        if(judgeInvalid(pixelX,pixelY)) return false;
        //获取目标像素单元
        DomUnit& tempUnit=getUnit(pixelX,pixelY);
        //在最后的unit里面写上像素坐标
        dstUnit.pixelX_=pixelX;
        dstUnit.pixelY_=pixelY;
        //在最后的结果里面记录属于父类的部分
        //重要的是这里面记录了一个指针
        dstUnit.copyFrom(tempUnit);
        //记录点云坐标
        dstUnit.scenePt_<<pixelX*pixelLength_+cloudRange_[0],
                pixelY*pixelLength_+cloudRange_[1];
        //获取颜色数据
        dstUnit.corresColor_=&(domResult_.at<cv::Vec3b>(domHeight_-pixelY-1,pixelX));
        //正常情况下返回true
        return true;
    }

    //获取dom图上的颜色指针
    CvColorPtr getCvColorPtr(int pixelX,int pixelY)
    {
        return domResult_.ptr<uchar>(domHeight_-pixelY-1)+3*pixelX;
    }

    //传入连续点获取三维的点坐标
    DomUnit& getUnitByContiPoint(const double* point)
    {
        //把点云坐标变换成像素点的坐标
        int pixelX=(point[0]-cloudRange_[0])/pixelLength_;
        int pixelY=(point[1]-cloudRange_[1])/pixelLength_;
        //获取目标位置的点云
        return getUnit(pixelX,pixelY);
    }

    //直接获取dom平面上目标位置的信息
    DomUnit& getUnitByCloud(const Point3DBase& cloudPt)
    {
        return getUnitByContiPoint(cloudPt.data());
    }

    //把离散点转换成连续坐标
    template<typename PointT>
    PointT convertToContiousPoint(int idRow,int idCol)
    {
        //新建目标点
        PointT pointForReturn;
        pointForReturn[0]=idCol*getPixelLength()+cloudRange_[0];
        pointForReturn[1]=idRow*getPixelLength()+cloudRange_[1];
        //返回结果
        return pointForReturn;
    }

    //根据点云的浮点型坐标获取完成的dom像素信息
    bool getUnitByCloud(const Vec3& cloudPt,DomUnitAll& dstUnit)
    {
        //把点云坐标变换成像素点的坐标
        int pixelX=convertConti2Pixel(cloudPt[0]);
        int pixelY=convertConti2Pixel(cloudPt[1],1);
        //获取目标的像素单元
        return getUnit(pixelX,pixelY,dstUnit);
    }

    //添加像素的颜色信息
    void addPixel(int pixelX,int pixelY,cv::Vec3i colorValue,WeightT weight=1)
    {
        //判断当前位置的像素是否合法
        if(judgeInvalid(pixelX,pixelY)) return;
        //获取当前位置的像素颜色
        cv::Vec3b &thisColor=domResult_.at<cv::Vec3b>(domHeight_-pixelY-1,pixelX);
        //获取当前位置的像素信息
        DomUnit &thisUnit=getUnit(pixelX,pixelY);
        //判断是否为一个空的颜色
        if(thisUnit.weight_==0)
        {
            thisUnit.weight_=weight;
            thisColor=colorValue;
            return;
        }
        //如果不使用权值，直接进行权值比较，然后保留数据即可
        if(!useWeight_)
        {
            //判断是否为更大的权值
            if(weight>thisUnit.weight_)
            {
                //保留权值和颜色
                thisUnit.weight_=weight;
                thisColor=colorValue;
            }
            return;
        }
        //构造一个带权值的颜色
        UserWColor mainColor(weight);
        //复制颜色
        mainColor.Vec::operator=(colorValue);
        //把目标存储的颜色保存成userColor的形式
        UserWColor dstColor(weight,thisColor);
        //两个颜色的比较，保留亮度和权值比较大的
        dstColor.cmpColor(mainColor);
        //记录最终的颜色和权值
        thisUnit.weight_=dstColor.weight_;
        thisColor=dstColor;
    }

    //在DOM图的对应位置上添加颜色信息,这里传入的是点云的精确坐标
    void addCloudPixel(double xLoc,double yLoc,cv::Vec3i colorValue,WeightT weight=1)
    {
        //计算点云所在的像素坐标
        int pixelX=std::round((xLoc-cloudRange_[0])/pixelLength_);
        int pixelY=std::round((yLoc-cloudRange_[1])/pixelLength_);
        //计算当前位置的像素中心
        std::pair<double,double> pixelCenter(pixelLength_*(0.5f+pixelX),
                              pixelLength_*(0.5f+pixelY));
        //计算像素坐标与像素中心的距离,为了防止除0操作，附加0.5个像素单元距离
        WeightT pixelDistance=0.5f*pixelLength_+std::sqrt(
                    std::pow(xLoc-pixelCenter.first,2)+std::pow(yLoc-pixelCenter.second,2));
        //添加当前位置的像素信息
        addPixel(pixelX,pixelY,colorValue,weight/pixelDistance);
    }

    //传入一个颜色序列，每个颜色对应一个坐标，把它们添加到dom图里面
    void addColorPtList(const ColorPtList &colorPtList)
    {
        //遍历每个颜色数据
        for(ColorPtList::const_iterator iter=colorPtList.begin();
            iter!=colorPtList.end();++iter)
        {
            if(iter->weight_<SIN_THRE) continue;
            //记录当前位置的颜色数据
            addCloudPixel(iter->xCoord(),iter->yCoord(),iter->ptColor_,iter->weight_);
        }
    }

    //批量添加三维点信息
    void addManyZValue(const Point3DList& ptList)
    {
        //遍历待添加的每个点，这里可以考虑使用多线程
        for(unsigned int ptId=0;ptId<ptList.size();++ptId)
        {
            //当前位置的点
            const Point3D& thisPt=ptList[ptId];
            //根据xy坐标获取当前位置的dom单元
            DomUnit& thisUnit=getUnitByCloud(thisPt);
            //判断当前dom单元是否已经包含了先验信息
            if(thisUnit.isPriorPoint_) continue;
            //记录里面的z值
            thisUnit.addZ(thisPt[2]);
        }

    }

    //用一系列三角面片初始化每个位置的z值
    void initGlobalZValue(const DomPatchList& patchList)
    {
        //遍历每个位置的三角面片
        //这里不宜使用多线程，不同面片之间的访问可能有冲突
        for(DomPatchList::const_iterator iter=patchList.begin();
            iter!=patchList.end();++iter)
        {
            Point3DList ptList;
            try {
                //获取当前位置的点列表
                iter->getPtInTriangle(ptList,pixelLength_/2);
                //计算三维点里面的z坐标
                iter->computeManyZ(ptList);
                //让z值适面片的范围，不然房沿的位置点会特别离谱
                iter->adjustZInRange(ptList);
            } catch (int errorFlag) {
                //的可能三个点太密集，导致它图片的范围变成了0
                if(errorFlag==ERROR_ZERO_IMG_SIZE)
                    iter->getPlanePts(ptList);//这种情况下只获取点云里面的3个点
            }
            //把三维点列表添加到dom信息集合中
            addManyZValue(ptList);
        }
    }

    //利用cv里面的点列表初始化优化器数据
    //当最后一个参数存在的时候，只取边缘点和已经优化过的点
    void makeOptimizerVector(const CvPointVector& srcPtList,
                             OptimizerVector& dstVector,const ImgRange* rangeRef=nullptr,double nccThre=NCC_THRE)
    {
        //根据点的个数，初始化优化器列表的个数
        dstVector.reserve(srcPtList.size());
        //遍历点的列表
        for(CvPointVector::const_iterator iter=srcPtList.begin();
            iter!=srcPtList.end();++iter)
        {
            //根据点坐标获取当前位置的unit
            DomUnit& tempUnit=getUnit(iter->x,iter->y);
            //如果是一个无效的单元，就跳过
#if defined (GLOBAL_HELP) || (!defined (USE_Z_PLANE))//使用z平面的算法的时候不需要考虑外点的问题
            if(tempUnit.isOutRangePt_) continue;
#endif
            //判断是否需要进一步筛选边缘点
            if(rangeRef)
            {
                //判断是否为优化过的点
                if(tempUnit.nccScore_<nccThre && (!rangeRef->judgeMidSide(iter->x,iter->y)))
                    continue;
            }
            //往优化器的列表里面记录这个数据，计算实际的点云位置
            dstVector.push_back(DomZOptimizer(iter->x*pixelLength_+cloudRange_[0],
                                iter->y*pixelLength_+cloudRange_[1],&tempUnit.z_));
            //获取刚刚添加的数据
            DomZOptimizer &newPt=dstVector[dstVector.size()-1];
            //记录它是否为先验点
            newPt.isPrior_=tempUnit.isPriorPoint_;
            //记录该点的ncc分数
            newPt.nccPtr_=&tempUnit.nccScore_;
            //记录该点对应的颜色数据
            newPt.colorPtr_=getCvColorPtr(iter->x,iter->y);
            //在它的对象信息里面记录z的范围
            newPt.zRange_=&zRange_;
            //记录可使用的次数
            newPt.refTimePtr_=&tempUnit.referencedTime_;
            //记录像素方向的指针
            newPt.pixDirPtr_=&tempUnit.pixDir_;
            //记录像素的顺序标号 20211126
            newPt.domUnitIdLabal_=convertLocalToIdLabal(iter->x,iter->y);
        }
    }

    //根据二维世界坐标获取当前DOM像素的位置
    //如果超出范围则返回最低位置
    double getZValueByCoordinate(double xValue,double yValue)
    {
        //把浮点型坐标转换成离散坐标
        int idCol=(xValue-cloudRange_[0])/pixelLength_;
        int idRow=(yValue-cloudRange_[1])/pixelLength_;
        //判断范围
        if(judgeInvalid(idCol,idRow))
        {
            return zRange_[0];
        }
        //获取目标位置的Z值
        DomUnit& dstUnit=getUnit(idCol,idRow);
        //如果没有被优化过，直接返回
        if(dstUnit.nccScore_<0.3) return zRange_[0];
        return dstUnit.z_;
    }

    //根据传入的图片的范围，获取z优化器的向量
    void getRect(const ImgRange& imgRange,OptimizerVector& dstVector,double nccThre=NCC_THRE)
    {
        //获取范围内的点列表
        CvPointVector cvPtList;
        imgRange.getPtsInRange(cvPtList);
#ifdef FIX_SAMPLE
        //传入范围，表示仅采样边缘点和优化过的点，采样优化过的点是为了辅助后面的优化
        makeOptimizerVector(cvPtList,dstVector,&imgRange,nccThre);
#else
        //根据点列表获取优化器列表
        makeOptimizerVector(cvPtList,dstVector,nullptr,nccThre);
#endif

    }

    //把传入的整数坐标转换成点云坐标
    //不负责z值
    void convertPixelLocal(int pixelX,int pixelY,Vec3& dstPt)
    {
        //分别赋值x和y
        dstPt[0]=pixelX*pixelLength_+cloudRange_[0];
        dstPt[1]=pixelY*pixelLength_+cloudRange_[1];
    }

    //获取指定位置的颜色数据，不是获取指针
    void getDomResultColor(int pixelX,int pixelY,cv::Vec3i& dstColor)
    {
        //获取指定位置的像素
        dstColor=domResult_.at<cv::Vec3b>(domHeight_-pixelY-1,pixelX);
    }

    //根据range目前的位置，获取它的下一个位置
    //直接在传入的位置上更新
    //如果已经遍历到头了，就返回true
    //这个东西也是可以考虑复用的，
    //如果将来有人愿意重构这个代码，请不要太恨我。。。
    bool getNextRange(ImgRange& currRange,int stepLen=MASK_STEP)
    {
        //将范围整体向右平移一步
        currRange.translate(stepLen,0);
        //判断右边界是否超了
        if(currRange.right()>=domWidth_)
        {
            //将范围向下平移一步
            currRange.moveTo(0,currRange.up()+stepLen);
            //判断下边界是否超了,这个时候说明已经遍历到头了
            if(currRange.down()>domHeight_)
                return true;
        }
        //记录目前的范围
#ifdef DEBUG_PRINT
        currRange.iterInfo_=std::to_string(currRange.left())+"/"+
                std::to_string(domWidth_)+"\t"+std::to_string(currRange.up())+"/"+
                std::to_string(domHeight_);
#endif
        //没有超过范围，正常返回
        return false;
    }
};//class DomInfo

//做德劳内三角化的工具，继承别人的做德劳内三角化的接口，然后封装一些别的功能
class MyDelaunator : public delaunator::Delaunator
{
    //点标号的列表，父类里面存储的三角面片对应的是点的顺序，但算法真正需要的并不是按照点的传入顺序的点标号
    const IndexList* idxListPtr_;
public:
    //构造函数，仿造父类 Delaunator(std::vector<double> const& in_coords)
    MyDelaunator(DelaunayPtList& ptList,const IndexList* idxList=nullptr) : delaunator::Delaunator(ptList)
    {
        //记录点传入顺序的标号
        idxListPtr_=idxList;
    }

    //获取三角面片的列表,传入的patchType必须满足传入三个点标号从而实现的初始化操作
    template<typename PatchType>
    void getPatchList(std::vector<PatchType>& dstPatchList)
    {
        //取出三角面片的点列表
        const std::vector<std::size_t> &patchCornerList=triangles;
        //提前开辟空间
        dstPatchList.reserve(patchCornerList.size()/3);
        //遍历已经找到的每个三角面片
        for(unsigned int ptCount=0;ptCount<patchCornerList.size();ptCount+=3)
        {
            //判断有没有可参考的标号列表
            if(idxListPtr_==nullptr)
            {
                dstPatchList.push_back(PatchType(patchCornerList[ptCount],
                        patchCornerList[ptCount+1],
                        patchCornerList[ptCount+2]));
            }
            else
            {
                const IndexList& idxList=*idxListPtr_;
                dstPatchList.push_back(PatchType(idxList[patchCornerList[ptCount]],
                        idxList[patchCornerList[ptCount+1]],
                        idxList[patchCornerList[ptCount+2]]));
            }
        }
    }
};

/// Generic SfM data container
/// Store structure and camera properties:
struct SfM_Data
{
  /// Considered views
  Views views;
  //图片之间需要用到的F矩阵是否已经被计算过了
  bool fPrepared_=false;
  //每个图片能看到的特征点列表是否已经初始化过了
  bool imgLookReady_=false;
  /// Considered poses (indexed by view.id_pose)
  Poses poses;
  /// Considered camera intrinsics (indexed by view.id_intrinsic)
  Intrinsics intrinsics;
  /// Structure (3D points with their 2D observations)
  Landmarks structure;
  //三角面片的列表
  DelaunayPatchList patchList_;
  //点云信息的指针列表,上面的structure哈希表会按照顺序存在这里面
    LandmarkPtrList strucPtrList_;

    //做快速DOM需要用到的相关信息
    DomInfo domInfo_;

  /// Controls points (stored as Landmarks (id_feat has no meaning here))
  Landmarks control_points;

  /// Root Views path
  std::string s_root_path;

  //可用图片的个数
  unsigned int imgNum_=0;

  //是否使用遮挡判断。
  bool isUsingOccludeDetection_=false;

  //是否锁定所有的Z值
  bool holdZValue_=false;

  //处理DOM时候的点云范围
  std::string cloudRangeInfo_="";

  //中间过程的输出目录 2022-8-18 添加的，为了适应其它机器上的运行
  std::string midProceeFolder_="";

  //--
  // Accessors
  //--
  const Views & GetViews() const {return views;}
  const Poses & GetPoses() const {return poses;}
  const Intrinsics & GetIntrinsics() const {return intrinsics;}
  const Landmarks & GetLandmarks() const {return structure;}
  const Landmarks & GetControl_Points() const {return control_points;}

  /// Check if the View have defined intrinsic and pose
  bool IsPoseAndIntrinsicDefined(const View * view) const
  {
    if (view == nullptr ) return false;
    return (
      view->id_intrinsic != UndefinedIndexT &&
      view->id_pose != UndefinedIndexT &&
      intrinsics.find(view->id_intrinsic) != intrinsics.end() &&
      poses.find(view->id_pose) != poses.end());
  }

  /// Get the pose associated to a view
  const geometry::Pose3 GetPoseOrDie(const View * view) const
  {
    return poses.at(view->id_pose);
  }

  //载入view的图片
  void loadViewImgs()
  {
      //如果使用动态载入图片的方案
      //那么在这里就暂时不载入图片
#ifdef DYNAMIC_LOAD_IMAGE
      return;
#endif
      //初始化可用图片的个数
      imgNum_=0;
      //遍历view
      for(auto& eachView : views)
      {
#ifdef USE_ONLY_MID_CAM
          if(!MID_ID_JUDGE(eachView.first)) continue;
#endif
          eachView.second->loadImg(s_root_path);
          //可用图片的个数增加
          imgNum_++;
      }
  }

  //调试的时候使用的东西
  void readRange(int& range1,int &range2)
  {
      //读取文本
      std::ifstream txtHandle;
      txtHandle.open("/media/cvlab/data/workSpace/mainProject/topViewConstruct/workspace/reconstruct/range.txt");
      //读取两行
      char buffer[200];
      txtHandle.getline(buffer,200);
      range1=atoi(buffer);
      txtHandle.getline(buffer,200);
      range2=atoi(buffer);
  }

  //依次输出点云里面每个点的坐标
  void printPointList()
  {
        //遍历sfm_data的structure属性
        for(auto &structureIter : structure)
        {
            //输出当前位置的标号
            std::cout<<structureIter.first<<" ";
            //获取当前位置的坐标
            Vec3 &thisPt=structureIter.second.X;
            std::cout<<thisPt[0]<<" "<<thisPt[1]<<" "<<thisPt[2];
            std::cout<<std::endl;
        }
  }

  //把structure里面的哈希表信息的指针存储在指针列表里面，方便按照顺序随机访问
  //这个东西并不随时更新，也不是经常使用，仅仅是需要的时候就用一下
  void makeStrucPtr()
  {
      //初始化指针列表
      strucPtrList_.clear();
      //提前开辟空间
      strucPtrList_.reserve(structure.size());
      //遍历点云里面的点
      for(auto &eachPt : structure)
      {
          //记录当前位置的指针
          strucPtrList_.push_back(&eachPt.second);
      }
  }

  //把点云信息保存成pointcloudMatch的形式，是自己写的一种文件保存规则
  //为了把点云中的点和图片上像素的对应关系保存下来，后面做快速DOM的时候要用
  void makePointcloudMatch(PointcloudMatch& tempMatch)
  {
        //图片的个数
        int imgAccount=views.size();
        //图片的大小
        int imgCols=views.begin()->second->ui_width;
        int imgRows=views.begin()->second->ui_height;
        //新建点云匹配器
        tempMatch.rebuild(imgAccount,imgCols,imgRows);
        //点的个数
        int pointAccount=structure.size();
        //往点云匹配器里面添加点的个数
        tempMatch.resizePoint(pointAccount);
        //初始化需要访问的点的id,它和哈希表里面的索引标号是两回事
        int ptIdx=0;
        //遍历每个点的信息
        for(auto &eachPair : structure)
        {
            //取出当前位置的pointmatch点容器
            PointcloudMatch::PointContainer &thisPtCter=
                    tempMatch.getPtContainer(ptIdx);
            //初始化图片的索引
            int imgIdx=0;
            //遍历能观察到当前点的所有图片
            for(const std::pair<IndexT, Observation> &eachObv : eachPair.second.obs)
            {
                //获取当前图片下的投影位置
                PointcloudMatch::PointCoord &thisUv=thisPtCter[eachObv.first];
                //当前位置点的投影坐标
                const double* const sfmUv=eachObv.second.x.data();
                //记录当前图片下的两个投影坐标
                thisUv.first=std::round(sfmUv[0]);
                thisUv.second=std::round(sfmUv[1]);
                //在图片相应的像素投影位置上记录下它的点云索引
                tempMatch.imgAt(imgIdx,thisUv.second,thisUv.first)=ptIdx;
                //图片的索引标号增加
                imgIdx++;
            }
            //标号增加
            ptIdx++;
        }
  }

  //对点云里面的点做德劳内三角化,保存在三角面片的列表里
  //最后的三角面片结果会被存储在patchList_属性里面
  void delaunayTriangulation2D()
  {
      //生成一个临时的点列表，用来适配第三方做德劳内三角化的接口
        std::vector<double> simplePtList;
        //提前开辟空间
        simplePtList.reserve(structure.size()*2);
        //对应地，记录每个位置的点的下标
        std::vector<IndexT> indexList;
        //给下标序列提前开辟空间
        indexList.reserve(structure.size());
        //遍历点云里面的点
        for(auto &structureIter : structure)
        {
            //获取当前位置的坐标
            const double* const thisPt=structureIter.second.X.data();
            //添加当前位置的坐标
            simplePtList.push_back(thisPt[0]);
            simplePtList.push_back(thisPt[1]);
            //记录当前点的下标
            indexList.push_back(structureIter.first);
        }
        //传入的德劳内三角化对象中，准备做三角化
        delaunator::Delaunator delaTool(simplePtList);
        //提前给三角面片的容器开辟空间
        patchList_.reserve(delaTool.triangles.size()/3);
        //遍历德劳内三角面片的结果
        for(std::size_t patchId=0;patchId<delaTool.triangles.size();patchId+=3)
        {
            //记录三角面片所属点的三个下标
            patchList_.push_back(TriangularPatch2D(
                                      indexList[delaTool.triangles[patchId]],
                                      indexList[delaTool.triangles[patchId+1]],
                                      indexList[delaTool.triangles[patchId+2]]));
        }
  }


  //制作易用接口,方便做三维的德劳内三角化
  //idList里面存储的是索引标号，方便从数据里面取回结果
  void make3DSimplePtList(std::vector<double> &dstPtList,std::vector<IndexT> &idList)
  {
      //提前开辟空间
      dstPtList.reserve(structure.size()*3);
      idList.reserve(structure.size());
      //遍历点云中所有的点
      for(auto &eachPt : structure)
      {
          //当前位置的点坐标
          Vec3 &thisPt=eachPt.second.X;
          //记录点坐标
          for(int i=0;i<3;++i) dstPtList.push_back(thisPt[i]);
          //记录村号
          idList.push_back(eachPt.first);
      }
  }

  //把patchList保存成自己需要的变量形式
  //第1个变量表示的是数据第次，第二个变量表示的是哈希索引
  void savePatchList(IndexList& cornerList, IndexList& idList)
  {
      //提前给最终的三角面片开辟空间
      patchList_.reserve(cornerList.size()/3);
      //遍历角点
      for(unsigned int cornerId=0;cornerId+2<cornerList.size();cornerId+=3)
      {
          //根据标号生成三角面片
          patchList_.push_back(TriangularPatch2D(idList[cornerList[cornerId]],
                               idList[cornerList[cornerId+1]],
                               idList[cornerList[cornerId+2]]
                                   ));
      }
  }

  //对所有的点云做三维德劳内三角化
  void delaunayTriangulation3D()
  {
        //制作简单的三维德劳内三角化接口
        std::vector<double> simplePtList;
        IndexList idList;
        make3DSimplePtList(simplePtList,idList);
        //获取三角面片的数据,cornerList里面存储的是从0开始的一组数，它表示struct里面的第几个点
        //这个"第几个点"指的不是first索引
        IndexList cornerList;
        delaunayTriangulation(simplePtList,cornerList);
        //保存面片的数据
        savePatchList(cornerList,idList);
  }

  //判断某个图片能不能看到某个点
  bool judgeImgLook(IndexT imgId,IndexT cloudId)
  {
      return structure.at(cloudId).obs.count(imgId);
  }


  //判断一个统计信息里面是否有2,仅仅是临时随便写的函数
  bool judgeStatHaveValue(const StaticInfo& statInfo, const int value)
  {
      //遍历统计信息
      for(StaticInfo::const_iterator iter=statInfo.begin();
          iter!=statInfo.end();++iter)
      {
          if(iter->second==value) return true;
      }
      return false;
  }


  //给一个点云中的点添加投影信息，具体的添加了指定点云的id和图片的id，然后尝试添加
  //添加的时候是通过重投影的方式添加的
  //如果添加成功就返回true
  bool addProjectPoint(IndexT cloudId,IndexT imgId)
  {
      //找到点云的引用
      Landmark& cloudPt=structure.at(cloudId);
      //判断这个图是否本来就能看到这个点
      if(cloudPt.obs.count(imgId)>0)
      {
          return false;
      }
      //把这个点投影到这个图上
      Vec2 projPt=projectCloudPoint(cloudPt,imgId);
      //判断投影点是否在图片范围内
      if(!views.at(imgId)->inImgRange(projPt[0],projPt[1])) return false;
      //获取这个点的平均颜色
      UserColor avgColor;
      //判断这个点是否已经有了平均颜色
      if(cloudPt.avgInited)
      {
          avgColor.Vec::operator=(cloudPt.avgColor);
      }
      else
      {
          pointAvgColor(cloudPt,avgColor);
          cloudPt.avgColor=avgColor;
          cloudPt.avgInited=true;
      }
      //获取点在当前位置的颜色
      UserColor thisColor(views.at(imgId)->colorAt(projPt[0],projPt[1]));
      //计算当前颜色与平均颜色的距离
      double colorDis=avgColor.colorDis(thisColor);
      //判断两个颜色是否足够接近
      if(colorDis<COL_THRE) return false;
      //在obs里面添加这个数据
      cloudPt.obs.insert(std::make_pair(imgId,Observation(projPt,imgId)));
      //正常情况下返回true
      return true;
  }

  //对三角面片的进一步保护,这种情况下，三角面片上的每个点都只有一个数据能看到它
  bool protectForward(DelaunayPatchList::iterator dyingPatchIter,StaticInfo& statInfo)
  {
      //遍历所有的统计信息
      for(StaticInfo::iterator iter=statInfo.begin();iter!=statInfo.end();++iter)
      {
          //添加新的观测图，并统计添加成功的次数
          for(int i=0;i<3;++i)
          {
              if(addProjectPoint(dyingPatchIter->corners_[i],iter->first))
              {
                  //添加成功了就把它计数
                  iter->second++;
              }
          }
          //判断是否达到了三个
          if(iter->second>=3)
          {
              //记录新的观测者
              dyingPatchIter->obvs_.push_back(iter->first);
              return true;
          }
      }
      //没救成的情况下返回false
      return false;
  }


  //拯救即将消失的三角面片,第二个参数是当前三角面片的统计信息
  bool protectPatch(DelaunayPatchList::iterator dyingPatchIter,StaticInfo& statInfo)
  {
      //遍历所有的统计信息，看看有没有哪个点能通过重投影获取一个新的观测点
        for(StaticInfo::iterator statIter=statInfo.begin();statIter!=statInfo.end();++statIter)
        {
            //判断是否已经看到了两个点
            if(statIter->second!=2) continue;
            //这个图没看到的是哪个点
            int lostId=dyingPatchIter->whoLostThisImg(structure,statIter->first);
            //在观测信息里面添加这个图片，如果可以添加成功的话，那就拯救成功
            if(addProjectPoint(dyingPatchIter->corners_[lostId],statIter->first))
            {
                //给当前的统计信息+1
                statIter->second++;
                //记录新的观测者
                dyingPatchIter->obvs_.push_back(statIter->first);
                //表示救成了
                return true;
            }
        }
        //效果还是不行，只好采取进一步的救援
      return protectForward(dyingPatchIter,statInfo);
  }


  //计算F矩阵,在调用这个之前需要先把三角面片给算出来
  //这里面同时还有标记三角面片的观测点的功能，所以就算不需要F，也是需要调用一下的
  void prepareFundamentalMatrixs()
  {
      //非法的面片统计
      unsigned int invPatchCount=0;
      //面片的标号
      unsigned int patchId=0;
      //遍历每个三角面片
      for(DelaunayPatchList::iterator iter=patchList_.begin();
          iter!=patchList_.end();++iter)
      {
          //取出当前位置正在迭代的三角面片
          TriangularPatch2D &thisPatch=*iter;
          //有哪些图能同时看到三角面片上的这三个点,做统计
          std::map<IndexT,int> obvStati;
          //遍历三角面片的三个角点
          for(int cornerId=0;cornerId<3;++cornerId)
          {
              //当前迭代位置的角点
              Landmark &thisPoint=structure[thisPatch.corners_[cornerId]];
              //遍历能看到当前点的所有图
              for(auto &eachObv : thisPoint.obs)
              {
                  //当前位置的图片标号计数
                    obvStati[eachObv.first]++;
                    //判断计数的标号是否达到了3个
                    if(obvStati[eachObv.first]==3)
                    {
                        //在当前三角面片里面记录一下，都有谁能完整地看到这个三角面片
                        thisPatch.obvs_.push_back(eachObv.first);
                    }
              }
          }
          //判断到底有没有图片可以完全看到这个三角面片
          if(thisPatch.obvs_.size()==0)
          {
              //尝试拯救一下这个面片
              if(!protectPatch(iter,obvStati))
              {
                  //标记为非法面片
                  thisPatch.isInvalid_=true;
                  //非法面片计数
                  invPatchCount++;
                  continue;
              }
          }
          patchId++;
#ifdef NEED_F
          //计算所有图两两之间的F矩阵
          for(IndexList::iterator firstIter=thisPatch.obvs_.begin();
              firstIter!=thisPatch.obvs_.end();++firstIter)
          {
              //取出第1个图
              View &firstView=*views[*firstIter];
              //第二层的遍历
              for(IndexList::iterator secondIter=thisPatch.obvs_.begin();
                  secondIter!=thisPatch.obvs_.end();++secondIter)
              {
                  //判断第1层和第2层是否一样
                    if(firstIter==secondIter) continue;
                    //判断第1个图到第2个图的F是否已经算过
                    if(firstView.fundmentals_.count(*secondIter)==1)
                        continue;
                    //初始化第1个图到第2个图的F矩阵
                    firstView.fundmentals_[*secondIter]=FundamentalMatrix();
                    //取出第2个图的view
                    View &secondView=*views[*secondIter];
                    //判断第2个图到第1个图的F是否已经计算过
                    if(secondView.fundmentals_.count(*firstIter)==1)
                    {
                        //第1个图到第2个图的F矩阵就是第2个图到第1个图F矩阵的转置
                        firstView.fundmentals_[*secondIter].fundData_=
                                secondView.fundmentals_[*firstIter].fundData_.transpose();
                        continue;
                    }
                    //取出两个图的旋转和相机光心
                    Pose3 &firstPose=poses[firstView.id_pose];
                    Pose3 &secondPose=poses[secondView.id_pose];
                    //计算第2个图到第1个图的相对旋转和相对的光心位置,相对旋转矩阵不太确定
                    Mat3 relateR=secondPose.rotation()*firstPose.rotation().transpose();
                    Vec3 relateC=secondPose.center()-firstPose.center();
                    //两个图的内参矩阵
                    Mat3 firstIntr=intrinsics[firstView.id_intrinsic]->getIntrinsicMatrix();
                    Mat3 secondIntr=intrinsics[secondView.id_intrinsic]->getIntrinsicMatrix();
                    //计算第1个图到第2个图的F矩阵
                    firstView.fundmentals_[*secondIter].computeF(relateR,relateC,firstIntr,secondIntr);
              }
          }
#endif
      }
        //标记一下，F矩阵已经算过了
      fPrepared_=true;
      //输出非法面片的统计
#ifdef DEBUG_PRINT
      std::cout<<"invalid patch count: "<<invPatchCount<<std::endl;
#endif
  }

  //利用XY平面的约束算法，计算出它在DOM平面上的约束直线
  LineT getDomLine(IndexT imgId,const TempPtBase& imgPt) const
  {
      //当前图片的旋转矩阵
      const Mat3& thisRMat=poses.at(views.at(imgId)->id_pose).rotation();
      //当前图片的光心坐标
      const Vec3& camCenter=poses.at(views.at(imgId)->id_pose).center();
      //获取像素点在相机坐标系下的投影
      Vec3 xcPoint=intrinsics.at(views.at(imgId)->id_intrinsic)->getXcPoint(imgPt[0],imgPt[1]);
      //r2和xc的点积
      double r2Xc=(thisRMat.col(1).transpose()*xcPoint)[0];
      //r1和xc的点积
      double r1Xc=(thisRMat.col(0).transpose()*xcPoint)[0];
      //生成用于返回的直线[(x-Cx)r2-(y-Cy)r1]*Xc=0
      LineT retLine;
    retLine<<r2Xc,-r1Xc,-camCenter[0]*r2Xc+camCenter[1]*r1Xc;
    return retLine;
  }

  //把这3个点画到dom图上
  void drawCornerPoints(const Pt3Group &ptIdList)
  {
        //遍历每个角点
      for(unsigned int ptId=0;ptId<ptIdList.size();++ptId)
      {
          //当前点的平均颜色值
          cv::Vec3i avgColor(0,0,0);
          //取出当前位置的点
          Landmark &thisCloudPoint=structure[ptIdList[ptId]];
          //遍历当前点的每个视角
          for(auto &eachId : thisCloudPoint.obs)
          {
              //根据eachId的first，在对应的view下取出它的图片
              cv::Mat &thisImg=views[eachId.first]->cvImg_;
              //判断图片是否读取过了
              if(thisImg.empty())
              {
                  //图片的路径
                  std::string imgPath=s_root_path+"/"+views[eachId.first]->s_Img_path;
                  //读取它的图片
                  thisImg=cv::imread(imgPath);
              }
              //获取当前位置的像素坐标
              const double* const pixelCoord=eachId.second.x.data();
              //当前位置的颜色
              cv::Vec3b tempColor=thisImg.at<cv::Vec3b>(cv::Point2i(pixelCoord[0],pixelCoord[1]));
              //累加点的平均颜色值
                avgColor+=tempColor;
          }
          //当前点共有多少个观察者
          int obvAccount=thisCloudPoint.obs.size();
          //计算点的颜色平均值
          avgColor=avgColor/obvAccount;
          //当前位置的点坐标
          const double* const pointCoord=thisCloudPoint.X.data();
          //把点的颜色平均值添加到DOM图的对应位置上
            domInfo_.addCloudPixel(pointCoord[0],pointCoord[1],avgColor,obvAccount);
      }
  }

  //填充三角面片的一个边,目前采用方案2,需要确保三角面片上的三个点处于一个真实存在的平面上
  //---------------------方案1------------------------
  //1.以最长的边为基准，最长边上的每个点与dom平面上直线的交点为新的稠密点
  //2.用F矩阵求对极线，对极线和其它图片上对应直线的交点被认为是上面那个稠密点的同名点
  //3.用所有的同名点的加权平均作为该稠密点的颜色
  //----------------------方案2-------------------
  //依次遍历每个图同名边的每个像素，像素点与dom平面直线的交点被作为新的稠密点
  //找到的若干稠密点融合起来作为最终的dom图
  //融合的时候像素值采用加权平均，权值与两直线的夹角有关，夹角越大权值越大
  //----------------------方案3-----------------------
  //1.以最长边为基准，用F矩阵计算对极线，在其它图片的同名直线上寻找同名点
  //2.把同名点集综合起来计算稠密点
  void fillALine(PatchIterInfo& iterInfo)
  {
        //遍历iterInfo里面保存的每个图片,其实是每个图片上的直线
      for(ViewImgLines::iterator imgIter=iterInfo.imgLines_.begin();
          imgIter!=iterInfo.imgLines_.end();++imgIter)
      {
          //对于长度小于1个像素的直线，不做处理,这里的像素就单纯是图片像素
          if(imgIter->getLength()<2) continue;
          //获取当前图片上的像素点集,点集里面每个点的父类表示的是在dom图上的坐标，
          //里面的imtPt属性表示的是在图片上的坐标
          imgIter->prepareColorPtList(true);
          ColorPtList &pixelPtList=imgIter->colorPtList_;
          //遍历图片上的每个点
#ifdef USE_OMP
#pragma omp parallel for
#endif
          for(unsigned int ptId=0;ptId<pixelPtList.size();++ptId)
          {
              //获取当前位置的图片点坐标
              const TempPtBase& thisPt=pixelPtList[ptId].imgPt_;
              //当前位置的目标颜色坐标
              ColorDomPt& dstColorPt=pixelPtList[ptId];
              //根据点坐标获取在DOM平面上的投影直线
              LineMid domLine=getDomLine(imgIter->imgId_,thisPt);
              //投影直线和dom平面上的直线的夹角正弦,后面会作为权值的参考
                double sinValue=iterInfo.domLine_.sinAngle(domLine);
                //计算两个dom直线的交点
                dstColorPt=iterInfo.domLine_.lineCross(domLine);
                //计算当前交点的颜色,另外需要判断点是否在图片的合理范围内
                if(!views.at(imgIter->imgId_)->inImgRange(thisPt[0],thisPt[1])) continue;
                dstColorPt.initColor(views.at(imgIter->imgId_)->cvImg_.at<cv::Vec3b>(
                                         cv::Point2i(thisPt[0],thisPt[1])));
                //当前位置的颜色权值,另外一个权值是图片上的直线和DOM直线的长度比值
                //if(sinValue>0.2) sinValue=1.f;
                dstColorPt.weight_=sinValue;
                //判断点是否在dom线段的范围内
                if(!iterInfo.domLine_.judgeBetween(dstColorPt))
                    dstColorPt.weight_=0;
          }
          //把记录好的颜色点列表添加到dom图里面
          domInfo_.addColorPtList(pixelPtList);
      }
  }

  //这个时候iterInfo的最长边已经被填充了，需要把另一个角与边上的每个点连接然后生成新的稠密点
  void filleCornerLine(PatchMainIterInfo& mainIterInfo)
  {
      //遍历里面的每个图片序列
      for(ViewImgLines::iterator iter=mainIterInfo.imgLines_.begin();
          iter!=mainIterInfo.imgLines_.end();++iter)
      {
          //第3个点在该图上的投影坐标
          const TempPtBase& otherImgPt=
                  structure.at(mainIterInfo.otherCorner_.idx_).obs.at(iter->imgId_).x;
          //遍历当前图上的每个点
          for(ColorPtList::iterator ptIter=iter->colorPtList_.begin();
              ptIter!=iter->colorPtList_.end();++ptIter)
          {
              //权值不够的数据不需要再处理
              if(ptIter->weight_<SIN_THRE) continue;
              //生成临时的迭代信息
              PatchIterInfo tempIter;
              //两个点在dom图上形成的直线
              tempIter.domLine_.initLine(mainIterInfo.otherCorner_,*ptIter);
              //在迭代信息里面加入一个空直线
              tempIter.imgLines_.push_back(DomImgLine());
              DomImgLine &thisImgLine=tempIter.imgLines_[0];
              //修改图片标号
              thisImgLine.imgId_=iter->imgId_;
              //调用它父类的接口，初始化直线
              thisImgLine.DomLine::initLine(ptIter->imgPt_,otherImgPt);
              //为了防止多线程的时候访问冲突，提前准备直线方程
              thisImgLine.getLineEqu();
              //填充这个新生成的直线
              fillALine(tempIter);
          }
      }
  }

  //把点云稠密化，我憋了好久才憋出来的
  //这里做的是按照三角面片的顺序，然后依次处理三角面片所属的图片
  //最后做出来的效果不好，准备换一下顺序，依次遍历每个图片，然后遍历这个图片能看到的所有三角面片
  void densifyPointcloud()
  {
        //判断F矩阵是否已经计算过了
      if(!fPrepared_)
      {
          //计算F矩阵
          prepareFundamentalMatrixs();
      }
      //根据目前的稀疏点云，确定最终的图片大小和分辨率，如果图片太大的话，需要卡一个阈值
      domInfo_.initDom(structure);
      //不使用权值了
      domInfo_.useWeight_=false;
      //读取范围,调试使用
      int range1,range2;
      readRange(range1,range2);
      //遍历每一个三角面片
      for(unsigned int patchIdx=0;patchIdx<patchList_.size();++patchIdx)
      {
          if(patchIdx>range2||patchIdx<range1) continue;
          if(patchIdx%1000==0)
            std::cout<<patchIdx<<"/"<<patchList_.size()<<std::endl;
          if(patchIdx%10000==0)
              cv::imwrite("/media/cvlab/data/workSpace/mainProject/topViewConstruct/workspace/reconstruct/domResult.bmp",domInfo_.domResult_);
            //当前位置正在遍历的面片
          TriangularPatch2D &thisPatch=patchList_[patchIdx];
          //把最基本的三个角点信息画到最后的输出图上,同时需要预加载一下每个view的图片信息
          drawCornerPoints(thisPatch.corners_);
          //判断当前面片是否合法
          if(thisPatch.isInvalid_) continue;
          //获取首次迭代需要做的迭代器
            PatchMainIterInfo mainIterInfo;
            thisPatch.getFirstPatchIterInfo(mainIterInfo,structure);
           //填充当前三角面片的主边
            fillALine(mainIterInfo);
            //填充三角面片的另一个角与其它边上的连线
            filleCornerLine(mainIterInfo);
      }
  }

  //输出最后的DOM图,建议保存BMP的图，质量会好一些
  void saveDomResult(std::string savePath)
  {
      cv::imwrite(savePath,domInfo_.domResult_);
  }

  //记录每个图片能看到的点的列表，理论上该函数只需要调用一次
  void prepareImgLookingPt()
  {
      //判断每个图片能看到的点是否已经被初始化过了
      if(imgLookReady_)
      {
          return;
      }
      else
      {
          imgLookReady_=true;
      }
      //遍历每个点
      for(const auto& eachPt : structure)
      {
          //能看到当前点的每个图片
          const Observations &obvs=eachPt.second.obs;
          //遍历能看到当前点的每个图片
          for(const auto& eachObv : obvs)
          {
              //在图片上记录当前点的标号
              views.at(eachObv.first)->ptIdList_.push_back(eachPt.first);
          }
      }
  }

  //制作适用于德劳内三角化的接口的点坐标列表
  void makeDelaunayInterfacePtList(const IndexList& ptIdList,//点的标号列表，标号指的是struct里面的点标号
                                   const IndexT imgId,//需要做三角化的图片ID
                                   DelaunayPtList& dstPtList //生成的点坐标列表会被存储的这里
                                   )
  {
      //提前给最后的点结果开辟空间
      dstPtList.reserve(ptIdList.size()*2);
        //遍历点标号列表里面的每个点
      for(IndexList::const_iterator iter=ptIdList.begin();iter!=ptIdList.end();++iter)
      {
          //取出当前位置的点坐标
          const Vec2& thisLocal=structure.at(*iter).obs.at(imgId).x;
          //记录两个点坐标
          dstPtList.push_back(thisLocal[0]);
          dstPtList.push_back(thisLocal[1]);
      }
  }

  //针对全局点去制作简单的点列表
  //和另外一个函数make3DSimplePtList的制作过程是类似的
  void make2DSimplePtList(IndexList& idxList,DelaunayPtList& dstPtList)
  {
      //提前给两个点集开辟空间
      idxList.reserve(structure.size());
      dstPtList.reserve(structure.size()*2);
      //遍历点云里面的每个点
      for(auto& eachPt : structure)
      {
          //记录当前位置的点坐标
          dstPtList.push_back(eachPt.second.X[0]);
          dstPtList.push_back(eachPt.second.X[1]);
          //记录点的标号
          idxList.push_back(eachPt.first);
      }
  }

  //模仿mvs的时候使用的三角面片信息
  //这个三角面片的特殊之处在于它还继承的平面和三角形的相关接口
  //主要是用来初始化DOM像素上每个点的z坐标的
  void makeAdvancedTriPatch(DomPatchList& dstList)
  {
      //生成适用于接口的点列表和点标号列表
      IndexList idxList;
      DelaunayPtList ptList;
      make2DSimplePtList(idxList,ptList);
      //使用点的标号列表和点列表初始化德劳内信息
      MyDelaunator delaunator(ptList,&idxList);
      //从德劳内信息里面保存最后的结果
      delaunator.getPatchList<DomTrianglePatch>(dstList);
      //帮助每个三角面片记录点云里面的坐标信息，从而初始化它的平面信息
      DomTrianglePatch::updatePatchVector(dstList,structure);
  }


  //对指定的图片进行德劳内三角化
  void triangularImgPoints(const IndexT imgId,//需要做三角化的图片ID
                           ImgPatchList& dstPatchList //获取到的三角面片会最终存储在这里
                           )
  {
        //按照德劳内三角化的接口，制作点坐标列表
        DelaunayPtList ptList;
        makeDelaunayInterfacePtList(views.at(imgId)->ptIdList_,imgId,ptList);
        //把点列表传入德劳内三角化的对象中，准备做德劳内三角化
        MyDelaunator trigulator(ptList,&(views.at(imgId)->ptIdList_));
        //记录找到的三角面片
        trigulator.getPatchList<ImageTriPatch>(dstPatchList);
  }

  //把德劳内三角化的情况画在图片上
  //传入这个图片三角面片的列表以及对应的图片
  void drawDelaunayImg(ImgPatchList& patches,IndexT imgId)
  {
      //获取图片
      cv::Mat img=views.at(imgId)->loadImg(s_root_path).clone();
      //遍历每个三角面片
      for(ImgPatchList::iterator patchIter=patches .begin();
          patchIter!=patches.end();++patchIter)
      {
          //把当前三角面片画在图上
          patchIter->drawPatch(img,structure);
      }
      //把画好的图存储下来
      cv::imwrite(s_root_path+"/"+std::to_string(imgId)+"_tri.jpg",img);
  }


  //上面的稠密化点云使用的是对整个DOM图做德劳内三角化，然后遍历每个三角面片
  //这里要做的算法是，遍历每个图片，然后对图片上能看到的点进行三角化
  void denseDomImageByImage()
  {
      //初始化最终的DOM图
      domInfo_.initDom(structure);
        //准备每个图片能看到的点列表
      prepareImgLookingPt();
      //false的时候不使用权值，直接保留数据
      domInfo_.useWeight_=true;
      //遍历每个能看到的图片
      for(auto& eachView : views)
      {
          //如果这个图片看到的点太少，也跳过
          if(eachView.second->ptIdList_.size()<10) continue;
          std::cout<<eachView.first<<std::endl;
          if(eachView.first>299999||eachView.first<200000) continue;
          //对当前图片做德劳内三角化，获取三角面片的列表
          ImgPatchList imgPatches;
          triangularImgPoints(eachView.first,imgPatches);
          //载入图片，以免多线程处理的时候出现访问冲突
          eachView.second->loadImg(s_root_path);
          int patchId=0;
          std::ifstream txtHandle;
          //遍历每个三角面片
          for(ImgPatchList::iterator patchIter=imgPatches.begin();
              patchIter!=imgPatches.end();++patchIter)
          {
              patchId++;
              //std::cout<<patchId<<" ";
              //记录三角面片的图片标号
              patchIter->imgId_=eachView.first;
              //获取当前三角面片的主边迭代器,也就是最长边的迭代信息
                PatchMainIterInfo mainIter;
                patchIter->getFirstPatchIterInfo(mainIter,structure);
                //填充它的主边
                fillALine(mainIter);
                //填充每个边和对角线上的像素
                filleCornerLine(mainIter);

          }
          cv::imwrite("/media/cvlab/data/workSpace/mainProject/topViewConstruct/workspace/reconstruct/domResult.bmp",domInfo_.domResult_);
          //画三角化后的图
          //drawDelaunayImg(imgPatches,eachView.first);
      }
  }

  //根据xy坐标，获取每个位置上的投影直线
  void getProjLines(std::vector<ProjectLine>& projLineList,double domX,double domY)
  {
      //提前开辟view空间
      projLineList.reserve(views.size());
      //遍历每个view
      for(auto &eachView : views)
      {
#ifdef USE_ONLY_MID_CAM
          if(!MID_ID_JUDGE(eachView.first)) continue;
#endif
          //判断图片是否载入过
          if(eachView.second->cvImg_.empty()) continue;
          //判断这个view是否有外参，并不是每个view都有外参
          if(poses.count(eachView.second->id_pose)==0) continue;
          //新建直线信息
          projLineList.push_back(ProjectLine(poses.at(eachView.second->id_pose).rotation(),//旋转矩阵
                                     poses.at(eachView.second->id_pose).center(),//光心
                                             intrinsics.at(eachView.second->id_intrinsic)->getParams(),//内参
                                             domX,domY,eachView.first));//dom图上的坐标位置
      }
  }

  //根据z值，把每个线上的颜色都记录下来
  void getZColors(std::vector<ProjectLine>& projLineList,//直线的列表
                  double z,//z坐标
                  std::vector<VectorInterface<UserColorType>*>& dstColList)//最后的颜色序列会被存储到这里
  {
        //提前给颜色序列开辟空间，反正都是指针，不费什么
        dstColList.reserve(projLineList.size());
        //遍历每个直线
        for(std::vector<ProjectLine>::const_iterator iter=projLineList.begin();
            iter!=projLineList.end();++iter)
        {
            //获取该直线上的投影坐标,注意，这里一步到位，直接得到的就是图片上的像素坐标
            double pixelX,pixelY;
            iter->getXy(z,pixelX,pixelY);
            //判断是否在图片上
            if(!views.at(iter->imgId_)->inImgRange(pixelX,pixelY)) continue;
            //记录颜色
            dstColList.push_back(new UserColor(views.at(iter->imgId_)->colorAt(pixelX,pixelY)));
        }
  }

  //根据传入的直线集获取合适的颜色
  void getFitColorZ(std::vector<ProjectLine>& projLineList,//投影直线的列表
                    DomUnitAll& dstColorInfo,//最后计算出来的合适的颜色会被保存在这里
                    std::array<double,2>& zRange,//z的变化范围,下限，上限
                    double step) //步长
  {
      //遍历每个z
      for(double z=zRange[1];z>zRange[0];z-=step)
      {
            //获取z在当前位置上每个图片上的颜色投影
            std::vector<VectorInterface<UserColorType>*> colorList;
            getZColors(projLineList,z,colorList);
            //如果没有找到平均颜色就直接跳过
            if(colorList.size()==0) continue;
            //计算所有的向量的平均模长
            double avgNorm=UserColorVecBase::getAvgNorm(colorList);
            //对找到的所有颜色做归一化
            UserColorVecBase::uniformVectors(colorList);
            //计算平均颜色
            UserColor avgColor;
            UserColorVecBase::getAverageValue(colorList,avgColor);
            //计算所有颜色的角度标准差,当传入的所有向量都是归一化过的，可以把最后一项写成true
            double cosVariance=UserColorVecBase::computeCosVairance(colorList,avgColor,true);
            //到后面这里的颜色信息已经没用了，可以把颜色列表释放掉了
            UserColorVecBase::freeEachPtr(colorList);
            //判断角度方差是否足够小
            if(colorList.size()<MIN_CORR||cosVariance<CORR_THRE)
            {
                continue;
            }
            //如果超过了阈值，那就直接保留这个数据了,保留之前需要把颜色的模长恢复到平均的模长
            avgColor.setNormAs(avgNorm);
            dstColorInfo.saveColor(avgColor);
            //记录z值和权值
            dstColorInfo.z_=z;
            dstColorInfo.weight_=cosVariance;
            //正常通过之后直接结束
            break;
      }
  }

  //把稀疏点云画在DOM图上
  //返回可用的图片的个数
  void drawSparsPtAtDom()
  {
      std::cout<<"draw sparse"<<std::endl;
      //遍历每个structure
      for(auto& eachPt : structure)
      {
          //获取当前点的平均颜色
          UserColor avgColor;
#ifdef USE_IMG_Z_CONSTRAIN
          uint ptCount=pointAvgColor(eachPt.second,avgColor,true);
#else
          uint ptCount=pointAvgColor(eachPt.second,avgColor,false);
#endif
          //获取当前位置的dom信息,有可能对于边界位置是访问越界的，那就直接不要了
          DomUnitAll tempUnit;
          if(!domInfo_.getUnitByCloud(eachPt.second.X,tempUnit)) continue;
          //判断这里面的z值是否已经有了一个数据,z的默认数据是0
          //有可能某个数据的数值真的是0,等出现了这种情况再说吧，反正目前是不会出现这种情况的
          //这里面需要使用的是原始数据，也就是类属性里面的那个指针
          if(!tempUnit.srcUnitPtr_->addZ(eachPt.second.X[2])) continue;
          //记录平均颜色
          if(ptCount>0)
          {
            //把nnc分数标记成最大值
            tempUnit.srcUnitPtr_->nccScore_=PRIOR_NCC;
            tempUnit.saveColor(avgColor);
            //将当前位置的点标记为先验的点
            tempUnit.srcUnitPtr_->isPriorPoint_=true;
            //标记一下，它不是范围之外的点
            tempUnit.srcUnitPtr_->isOutRangePt_=false;
          }
          else
          {
              //把ncc的分数记录成0
              tempUnit.srcUnitPtr_->nccScore_=0;
              //将当前位置的点标记为先验的点
              tempUnit.srcUnitPtr_->isPriorPoint_=false;
              //标记一下，它不是范围之外的点
              tempUnit.srcUnitPtr_->isOutRangePt_=true;
          }
      }
  }

  //判断某个点在某个图片上的投影是否会存在倾斜遮挡
  bool isOcclude(Point3D& domPoint,IndexT idView,
                 IndexT stepSegNum=OCCLUDE_DETECT_SEG_NUM)
  {
        //获取目标位置的光心坐标
        Point3DBase camCenter=poses.at(views.at(idView)->id_pose).center();
        //计算两个点之间的向量
        Point3DBase vectorToCam=camCenter-domPoint;
        //场景的最大高度与当前位置的高度差
        double heightDiff=domInfo_.zRange_[1]-domPoint[2];
        //2022-10-26 临界位置还是需要多判断一段距离的
        //heightDiff*=1.5;
        //计算判别过程中每一步的向量
        Point3DBase eachStepVec=vectorToCam*(heightDiff/(stepSegNum+1.0)/vectorToCam[2]);
        //根据步长与DOM像素的比值判断是否需要考虑遮挡问题
        double eachStepLength=std::sqrt(std::pow(eachStepVec[0],2)+std::pow(eachStepVec[1],2));
        if(!needJudgeOcclude(eachStepLength/domInfo_.pixelLength_)) return false;

        //遍历步长中的每个位置
        for(IndexT idStep=0;idStep<=stepSegNum;++idStep)
        {
            //当前的判别位置
            Point3DBase currLocation=domPoint+(idStep+1)*eachStepVec;
            //获取当前位置的DOM高度
            double zValueOfCurrLocation=domInfo_.getZValueByCoordinate(currLocation[0],
                    currLocation[1]);
            //2022-10-26 每一步要求的阈值都不一样 总是要求拉高两倍以上
            double stepHeightHold=(currLocation[2]-domPoint[2])*0.5;
            //判断两个高度的差是否超过了1个阶层
#ifdef OCCLUDE_DETECT_ANY_HIGHER
            if(isDOMHigher(domInfo_.getZLevel(domPoint[2]),
                           domInfo_.getZLevel(zValueOfCurrLocation)))
#else
//            if(isDOMHigher(domInfo_.getZLevel(currLocation[2]),
//                           domInfo_.getZLevel(zValueOfCurrLocation)))
            //2022-10-26 为了更安全的去遮挡算法，更新处理的判断逻辑
            if(currLocation[2]-zValueOfCurrLocation<stepHeightHold)
#endif
            {
                return true;
            }
        }
        return false;
  }

  //根据domPatch的每个位置的z值，获取它们目前的颜色信息
  void getColorPatch(DomRectPatch& srcPatch,ColorPatchGroup& dstPatch,
                     FileHandler* filePtr=nullptr //log管理文件的指针
          )
  {
        //获取图片里面一共有多少个投影直线
      unsigned int imgNum=srcPatch.optiVec_[0].viewLines_.size();
      //提前给最后的颜色面片开辟空间
      dstPatch.patchList_.reserve(imgNum);
      //遍历每个投影直线,因为添加的过程中可能涉及到增加和删除，不能使用多线程
        for(unsigned int imgId=0;imgId<imgNum;++imgId)
        {
            //获取当前位置的图片集的引用
            dstPatch.patchList_.push_back(ColorPatch());
            ColorPatch& thisPatch=dstPatch.patchList_[dstPatch.patchList_.size()-1];
            //根据投影点的数量，提前开辟空间
            thisPatch.patchColors_.reserve(srcPatch.optiVec_.size());
            //如果需要限制z值的话，还需要对每个点的投影位置做一下预处理
#ifdef USE_IMG_Z_CONSTRAIN
            thisPatch.projList_.reserve(srcPatch.optiVec_.size());
#endif
            //记录一个临时的投影点列表，用于观察点的位置,这个东西仅仅用于debug
#ifdef SAVE_PROJ_PT
            CvPointVector tempCvPtList;
            tempCvPtList.reserve(srcPatch.optiVec_.size());
#endif
            //记录该面片的投影直线的标号
            thisPatch.viewLineId_=imgId;
            //遍历每个点在这个图片上的投影，依次添加进去
            for(OptimizerVector::iterator iter=srcPatch.optiVec_.begin();
                iter!=srcPatch.optiVec_.end();++iter)
            {
                //判断是否为第1次迭代，如果是第1次迭代，那就记录一下图片的标号
                if(iter==srcPatch.optiVec_.begin())
                    thisPatch.imgId_=iter->viewLines_[imgId].imgId_;
                //获取当前位置在图片上的投影点
                TempPtBase projPt;
                iter->getProjLocal(imgId,projPt[0],projPt[1]);
                //获取点在目标图片上的像素
                UserColor currColor;
#ifdef USE_IMG_Z_CONSTRAIN
                //获取当前位置的z高度
                int currZLevel=domInfo_.getZLevel(*iter->zPtr_);
#ifdef BAN_CENTER_FOR_CONSTRAIN
                //获取颜色
                if(!views.at(thisPatch.imgId_)->colorAt(projPt[0],projPt[1],currColor)) break;
                //获取目标位置的Z值
                int dstZLevel=views.at(thisPatch.imgId_)->getZLevel(projPt[0],projPt[1]);
                //判断取得的数据是否更大
                if(Z_LEVEL_VISIT_TRANS(currZLevel)<dstZLevel)
                {
                    //将当前的颜色面片标记为不可作为中心
                    thisPatch.banCenter_=true;
                }
#else
                //如果发现矩形面片里面有一个访问越界，那这个图片上的像素信息就全都不要了
                //如果访问到了过高的已经锁定的z值，会返回一个纯黑色的颜色，这个时候也会循环退出
                if(!views.at(thisPatch.imgId_)->colorAt(projPt[0],projPt[1],currColor,currZLevel)) break;
#endif
#else
                if(!views.at(thisPatch.imgId_)->colorAt(projPt[0],projPt[1],currColor)) break;
#endif
                //记录它在这个图片上的投影位置
                thisPatch.projList_.push_back(projPt);
                //判断一下找到的是不是一个纯黑的颜色，那很可能是去畸变后留下的黑边
                //对于每个图片的z高度注册结果，它也有可能因为高度值不如它原有的高度值而被弃用了
                if(currColor.judgeAllValue(0)) break;
                //判断是否需要参考主亮度
                if(srcPatch.mainBright_>0)
                {
                    //判断亮度的差别是否超过阈值
                    if(std::abs(currColor.getBright()-srcPatch.mainBright_)>BRIGHT_MAX_DIS) break;
                }
#ifdef SET_COLOR_NORM
                currColor.setNormAs(SET_COLOR_NORM);
#endif
                //把找到的颜色记录下来
                thisPatch.patchColors_.push_back(currColor);
                //记录找到的投影点位置 debug的时候使用
#ifdef SAVE_PROJ_PT
                tempCvPtList.push_back(CvPoint(projPt[0],projPt[1]));
#endif
            }
            //判断找到的颜色是否足够
            if(thisPatch.patchColors_.size()<srcPatch.optiVec_.size())
            {
                //如果不够，就直接不要它在这个图片上的投影了
                dstPatch.patchList_.pop_back();
            }
            else
            {
                //判断是否需要记录log信息
                if(filePtr)
                {
                    //获取所有的颜色
                    ColorTVec tempColorList;
                    thisPatch.getCvColorList(tempColorList);
                    //获取当前面片的颜色信息
                    filePtr->putColorInfo(thisPatch.imgId_,
                                          thisPatch.projList_,tempColorList);
                }
                //判断是否需要计算每个面片的权值
#ifdef USE_NCC_WEIGHT
                //获取相机的光心坐标
                const Vec3& camCenter=poses.at(views.at(thisPatch.imgId_)->id_pose).center();
                //新建一个向量用于表示中心位置的Z优化器的三维坐标
                Vec3 patchCenter=srcPatch.getCenterZOpt().get3DPoint();
                //计算两个向量的差值
                Vec3 ptVector=camCenter-patchCenter;
                //给权值赋值
                thisPatch.weight_=ptVector[2]/ptVector.norm();
                //改成正弦
                //thisPatch.weight_=std::sqrt(1-std::pow(thisPatch.weight_,2));
#endif
                if(isUsingOccludeDetection_)
                {
                    //判断当前面片是否为srcPatch指定的最佳面片 这是进入此函数之前就已经指定好的了
                    if(imgId==srcPatch.bestProjectView_)
                    {
                        //指定colorGroup里面的最佳颜色
                        dstPatch.bestImgId_=dstPatch.patchList_.size()-1;
                    }
                }
                else
                {
#ifdef USE_PRI_AS_AVG
                 //使用遮挡判断的时候，每个图片的投影范围在前面已经计算过了
                if(isUsingOccludeDetection_==false && MID_ID_JUDGE(thisPatch.imgId_))
                {
                    //获取中心迭代器
                    DomZOptimizer& centerOpt=srcPatch.getCenterZOpt();
                    //计算在当前图片上的上下界的投影距离
                    thisPatch.rangeProjDis_=centerOpt.getProjRangeDis(imgId);
                }
#endif
                }
                //纯属是为了方便debug
#ifdef SAVE_PROJ_PT
                std::ifstream tempIn;
                tempIn.open("/home/cvlab/empty0");
                //在view里面把这些点画到图片上并存储下来
                if(tempIn.is_open())
                    views.at(thisPatch.imgId_)->drawPoints(tempCvPtList);
#endif
            }
        }
  }

  //将一条垂直线记录到z值图上
  void recordRectZLine(double zValue,const ProjectLine& projLine,IndexT viewId)
  {
      //判断当前位置所属的z阶层
      int topZLevel=domInfo_.getZLevel(zValue);
      //当前位置的投影点
      CvPoint topProj;
      double tempX,tempY;
      projLine.getXy(zValue,tempX,tempY);
      topProj.x=std::round(tempX);
      topProj.y=std::round(tempY);
      //需要被操作的view
      View& mainView=*views.at(viewId);
      //变换到Z值图上
      topProj=mainView.transToZCoord(topProj.x,topProj.y);
      //将z值向下迭代依次标记遇到的像素点
      while(true)
      {
            //获取当前阶层的底部z值
          double nextZValue=domInfo_.getLevelHeight(topZLevel);
          //下一个位置的投影点
          projLine.getXy(nextZValue,tempX,tempY);
          CvPoint nextPt;
          nextPt.x=std::round(tempX);
          nextPt.y=std::round(tempY);
          //变换到Z值图的坐标空间
          nextPt=mainView.transToZCoord(nextPt.x,nextPt.y);
          //获取当前点和下一个点之间的中间点
          CvPointVector internalPt;
          ImgRange::getInternalPt(topProj,nextPt,internalPt);
          if(mainView.registerZMap(internalPt,topZLevel)) break;
          //更新z阶层
          --topZLevel;
          //更新z值
          zValue=nextZValue;
          //更新下一个位置的投影点
          topProj=nextPt;
      }
  }

  //记录面片上每个位置的z值
  void recordPatchZ(DomRectPatch& srcPatch,ColorPatchGroup &colorPatches)
  {
      //遍历patchGroup里面的每个面片
      for(ColorPatchList::iterator iter=colorPatches.patchList_.begin();
          iter!=colorPatches.patchList_.end();++iter)
      {
          //判断当前位置里面的点个数是否和面片里面的z优化器个数相同
          //按说一定是相同的
         if(iter->projList_.size()!=srcPatch.optiVec_.size()) throw ERROR_PT_ZOPT_DIFF_SIZE;
         //遍历每个位置的z优化器
         for(uint ptCount=0;ptCount<iter->projList_.size();++ptCount)
         {
             //判断是否使用完全的平面约束，如果是的话，则使用另一套逻辑
#ifdef CONSTRAIN_Z_LINE
             //记录当前位置对应的z直线
             recordRectZLine(*srcPatch.optiVec_[ptCount].zPtr_,
                             srcPatch.optiVec_[ptCount].viewLines_[iter->viewLineId_],
                     iter->imgId_);
             continue;
#endif
             //当前位置的投影点
             const TempPtBase& projPt=iter->projList_[ptCount];
             //当前位置的z阶层
             int zLevel=domInfo_.getZLevel(*srcPatch.optiVec_[ptCount].zPtr_);
             //在当前位置的投影点中记录z值
             views.at(iter->imgId_)->registerZMap(projPt[0],projPt[1],zLevel);
         }
      }
  }

  //为了debug添加的一个小功能，强制结束并保存点云结果
  bool forceEndJudge()
  {
      //目标文件目录
      std::string fileName="/home/cvlab/end";
      //打开文件
      std::ifstream fileHandle;
      fileHandle.open(fileName);
      //判断是否打开
      if(fileHandle.is_open()) return true;
      return false;
  }

  //通过优化法向量，对dom的面片做优化
  //到这个位置的时候，平面里面的信息还没有被初始化过
  //到这个位置的时候，每个优化器已经载入完毕了，每个优化器的对应到图片上的投影直线也就绪了
  //如果优化成功就返回true
  //如果有需要的话，最后的ncc结果可以被存储在nccScorePtr里面
  //当最后一个参数为true的时候，只要遇到了更大的数值，就记录一次颜色，如果是false,遇到了达到阈值的ncc
  //的时候才记录一次颜色
  //filePtr是log文件管理器的指针，不用的情况下传一个空的指针就行了
  bool optimizePatchNorm(DomRectPatch& srcPatch,double* nccScorePtr=nullptr,bool findMax=false,double nccThre=NCC_THRE,
                         FileHandler* filePtr=nullptr)
  {
        //初始化最大的ncc分数
        double maxNcc=0;
        //判断传入的ncc分数是否已经有值，如果已经有值，就用它作为临时的最大值
        if(nccScorePtr!=nullptr) maxNcc=*nccScorePtr;
        //通知三角面片缓存当前位置的数据
        srcPatch.saveData();
        //遍历必要的次数，寻找一个合适的ncc分数
        for(unsigned int iterTime=0;iterTime<NORM_OPT_TIMES;++iterTime)
        {
            //获取面片在每个位置上的投影颜色
            ColorPatchGroup currColor;
            getColorPatch(srcPatch,currColor,filePtr);
            //查看预先选定的图片是否存在
            if(isUsingOccludeDetection_)
            {
                //记录预先选定的ID和真实的ID
                if(filePtr)
                {
                    filePtr->putSelectedImageIndex(srcPatch.getBestProjectViewSrcImageId(),
                                                   currColor.getBestImageViewIndex());
                }
                if(currColor.getBestImageViewIndex()!=srcPatch.getBestProjectViewSrcImageId())
                {
                    throw ERROR_NO_MID_FOR_PRI;
                }
            }
            //如果能看到这个矩阵面片的图片少于某个数,那就没必要做约束
            if(currColor.patchList_.size()>=NCC_MIN_IMG)
            {
                //根据目前的颜色情况计算ncc分数
#ifdef USE_KMEANS
                double tempScore=currColor.getKmeansMax();
#else
#ifdef NCC_WRITE_AVG_COLOR
#ifdef USE_EUCLID_DIS
                double tempScore=currColor.computeEuclidScore(tempPriPt,filePtr);
                //对欧氏距离做变换
                tempScore=10.f/tempScore;
#else
                double tempScore=currColor.computeNccByAvg(filePtr);
#endif
                //判断是否需要记录中心颜色
                if(filePtr)
                {
                    //获取中心颜色
                    ColorTVec avgColor;
                    currColor.getCenterPatch().getCvColorList(avgColor);
                    filePtr->putAvgColor(avgColor);
                    //记录ncc的匹配结果
                    filePtr->putNccScore(tempScore);
                }
#else
                double tempScore=currColor.computeNcc();
#endif
#endif
                //判断是否需要对分数做根号处理
#ifdef NCC_EVO
                tempScore=std::pow(tempScore,1.f/currColor.patchList_.size());
#endif
                //记录最大的分数
                if(tempScore>maxNcc)
                {
                    maxNcc=tempScore;
                    //通知三角面片记录当前的数据情况
                    srcPatch.saveData(maxNcc);
                    //判断是否有必要记录颜色和法向量
                    if(findMax)
                    {
                        Point3DBase tempNorm;
#ifdef USE_NORM_PROG
                        srcPatch.getNormal(tempNorm);
#endif

                        //如果是固定Z值的情况，直接把颜色记录下来就可以走了
                        if(holdZValue_)
                        {
                            srcPatch.drawColorPatch(currColor.getBestColorList(),
                                                                                    tempNorm,1);
                            return true;
                        }
                        else
                            srcPatch.drawColorPatch(currColor.getBestColorList(),
                                                                                tempNorm,maxNcc);
                    }
                }
                else
                {
                    //恢复到修改之前的数据
                    srcPatch.loadData();
                }
                //判断最大的数据是否已经达标了
                if(maxNcc>nccThre)
                {
                    //判断是否需要在这里记录一下颜色
                    if(!findMax)
                    {
                        Point3DBase tempNorm;
#ifdef USE_NORM_PROG
                        srcPatch.getNormal(tempNorm);
#endif
                        srcPatch.drawColorPatch(currColor.getBestColorList(),
                                                                                tempNorm,maxNcc);
                    }
#ifdef USE_IMG_Z_CONSTRAIN
                        //记录每个位置的z值
                        recordPatchZ(srcPatch,currColor);
#endif
                    break;
                }
            }
            else //如果能看到这个面片的图片太少，也需要恢复到以前的状态
            {
                //如果第1次迭代就看不见，那就算了
                if(iterTime==0) break;
                srcPatch.loadData();
            }
            //判断这是不是最后一次迭代
            if(iterTime+1<NORM_OPT_TIMES)
            {
                //如果不需要更新Z值，那就直接返回
                if(holdZValue_) break;
                //将平面的法向量随机初始化
                srcPatch.updateRandNorm();
                //根据平面方程，更新每个位置的点
                srcPatch.updateOptimizer();
            }
        }
        //判断是否需要记录最大的分数
        if(nccScorePtr!=nullptr) *nccScorePtr=maxNcc;
        //返回最大分数的比较结果，如果超过最大分数，则视为优化成功了
        return maxNcc>nccThre;
  }

  //对包含z信息的矩形dom信息块做优化
  //第2个参数表示是否对z进行优化
  //如果优化成功了返回true，不成功则返回false
  bool optimizeRectPatch(DomRectPatch& srcPatch,bool optiZ=false)
  {
        //遍历面片的每个位置，初始化它们的投影直线
        for(OptimizerVector::iterator iter=srcPatch.optiVec_.begin();
            iter!=srcPatch.optiVec_.end();++iter)
        {
            //获取当前位置的投影直线
            getProjLines(iter->viewLines_,iter->xCoord(),iter->yCoord());
        }
        //更新矩阵面片的平面信息
        srcPatch.updatePlane();
        //新建一个double变量用于获取最大的分数数值
        double maxScore=0.f;
        //通过对法向量的优化，找到一个合适的三角面片的z值,如果优化成功了就直接返回了
        if(optimizePatchNorm(srcPatch,&maxScore,true)) return true;
        //如果优化不成功，判断是否需要对z做优化
        if(optiZ)
        {
            //调用优化z的函数
            //iteratePatchZ(srcPatch,domInfo_.zRange_,maxScore,STEP_L1,true);
            optimizePatchZ(srcPatch,maxScore);
            //这里一定会想办法让它强行优化成功的
            return true;
        }
        else
            return false;
  }

  //在排除倾斜遮挡的情况下，选取一个最适合的中相机作为当前面片的主要参考相机
  //如果发现所有的相机都是遮挡的，那就直接不弄了
  void occlusionSelect(DomRectPatch& srcPatch,FileHandler* filePtr=nullptr)
  {
      //获取当前面片的中心点
      Point3D centerPoint;
      srcPatch.getCenterPt(centerPoint);
      //中心位置的Z优化器
      DomZOptimizer& centerZOpt=srcPatch.getCenterZOpt();
      //初始化最佳的投影直线范围
      double bestProjectRange=100000;
      //最佳的投影直线标号
      IndexT bestProjectLineId=0;
      //遍历所有的投影直线
      for(IndexT idLine=0;idLine<centerZOpt.viewLines_.size();++idLine)
      {
          //当前位置的投影直线
          ProjectLine& currProjLine=centerZOpt.viewLines_[idLine];
#ifdef USE_COMPLEX_OCCLUDE_DETECTION
          //判断是否存在遮挡
          if(isOcclude(centerPoint,currProjLine.imgId_))
          {
              //把当前直线标记成不可用
              currProjLine.isUsable_=false;
              continue;
          }
          else
          {
              currProjLine.isUsable_=true;
          }
#endif
          //只处理中相机
          if(!MID_ID_JUDGE(centerZOpt.viewLines_[idLine].imgId_)) continue;
          //获取中心像素的投影位置
          TempPtBase projectPixel;
          currProjLine.getXy(centerPoint[2],projectPixel[0],projectPixel[1]);
          //判断中心像素是否在图片的投影范围内
          if(!views.at(currProjLine.imgId_)->inImgRange(projectPixel[0],projectPixel[1]))
          {
              continue;
          }
          //记录当前位置的投影直线范围
          double tempProjectRange=centerZOpt.getProjRangeDis(idLine);
          //判断是否需要记录投影直线范围信息
          if(filePtr)
          {
              //记录图片名称和投影直线的范围
              filePtr->recordProjectRange(currProjLine.imgId_,tempProjectRange);
          }
          //判断是否为更小的投影直线范围
          if(tempProjectRange<bestProjectRange)
          {
              //判断是否存在遮挡
#ifndef USE_COMPLEX_OCCLUDE_DETECTION
              if(!isOcclude(centerPoint,centerZOpt.viewLines_[idLine].imgId_))
#endif
              {
                  //记录投影直线的标号
                  bestProjectLineId=idLine;
                  bestProjectRange=tempProjectRange;
              }
          }
      }
      //如果没有找到合适的图片，报错
      if(bestProjectRange==100000)
      {
          throw ERROR_NO_MID_FOR_PRI;
      }
      else
      {
          //记录最佳的标号
          srcPatch.bestProjectView_=bestProjectLineId;
      }
  }

  //z平面算法状态下优化每个面片的位置
  bool optimizePlaneZPatch(DomRectPatch& srcPatch,double nccThre=NCC_THRE,
                           FileHandler* filePtr=nullptr //log文件的指针
    )
  {
#ifdef SAVE_EXTERNAL_PROJECT_LINE
      //如果使用外部存储直线的方式，那么这里通过外部存储的方式载入直线
      srcPatch.initProjectLinesByExternalData();
#endif
      //遍历面片的每个位置，初始化它们的投影直线
      for(OptimizerVector::iterator iter=srcPatch.optiVec_.begin();
          iter!=srcPatch.optiVec_.end();++iter)
      {
#ifndef SAVE_EXTERNAL_PROJECT_LINE
          //获取当前位置的投影直线
          getProjLines(iter->viewLines_,iter->xCoord(),iter->yCoord());
#endif
          //判断是否需要写入log
          if(filePtr)
          {
              {
                  //记录先验点信息
                  filePtr->putPriorInfo(iter->get3DPoint(),*(iter->nccPtr_),
                                        iter-srcPatch.optiVec_.begin());
              }
          }
          //如果每个像素有最大的可优化次数，将每个优化器的优化次数下降
#ifdef USE_MAX_OPT_TIME
          iter->minusOptTime();
#endif
      }
      double maxScore=0.f;
      //记录每个位置的z值
      Point3DList savedPtList;
      srcPatch.getPlanePts(savedPtList);
      //依次遍历面片的每个可能的高度阶层信息
      for(uint zId=0;zId<srcPatch.getZStepNum();++zId)
      {
            //把平面设置为指定的值
            //如果当前位置不可用，则跳过该位置
            if(!srcPatch.updateAsZStep(zId,holdZValue_)) continue;
            try {
                if(isUsingOccludeDetection_)
                {
                    //根据遮挡情况和投影情况，选择最合适的投影直线
                    occlusionSelect(srcPatch,filePtr);
                }
                //开始优化
                bool optFlag=optimizePatchNorm(srcPatch,&maxScore,true,nccThre,filePtr);
                //如果锁定Z值的话就只优化一次
                if(holdZValue_) return optFlag;
                //给每个点一次恢复到以前状态的机会
                srcPatch.loadPointList(savedPtList,nullptr,maxScore);
                //如果优化成功就停止
                if(optFlag) return true;
                else//如果没有优化成功，判断这是不是这个点的最后一次优化机会
                {
                    //判断中心位置的Z优化器是否已经没机会了
#ifdef USE_REFINE_TIME
                    if(*(srcPatch.getCenterZOpt().refTimePtr_)==0)
                    {
                        //强行记录当前位置的Z信息
                        ColorPatchGroup currColor;
                        getColorPatch(srcPatch,currColor);
                        recordPatchZ(srcPatch,currColor);
                    }
#endif
                }
            } catch (int errorFlag) {
                if(errorFlag!=ERROR_NO_MID_FOR_PRI)
                {
                    throw errorFlag;
                }
            }
            //判断是否为最后一次优化
            if(zId+1<srcPatch.getZStepNum())
            {
                //重新记录当前位置的点
                savedPtList.clear();
                srcPatch.getPlanePts(savedPtList);
            }
      }
      return maxScore>nccThre;
  }

  //遍历z的可能的范围，对每个可能的位置依次做法向量的优化
  void iteratePatchZ(DomRectPatch& srcPatch,//矩形面片
                     const RangeType& zRange,//z坐标的迭代范围
                     double maxScore=0,//目前的最大分数，通常不是0
                     double stepRate=STEP_L1, //步长的变化率
                     bool optForward=false //是否需要进一步做优化
          )
  {
      //计算z的实际步长
      double zStep=(zRange[1]-zRange[0])*stepRate;
      //初始化最大值对应的z位置
      double maxScoreZ=srcPatch.getCenterZ();
      //记录目前的点坐标列表
      Point3DList savedPtList;
      srcPatch.getPlanePts(savedPtList);
      //记录目前位置的直线方程
      PlaneT savedPlane=srcPatch.getPlaneEqu();
      //遍历每个位置的z值
      for(double currZ=zRange[1];currZ>=zRange[0];currZ-=zStep)
      {
            //更新矩阵面片里面每个位置的z值
            srcPatch.setAllZ(currZ);
            //保留当前位置的分数
            double tempScore=maxScore;
            //通过随机更改法向量来优化出更好的NCC匹配值
            //传入true表示一旦找到更合适的数据就把颜色保存下来
            bool optFlag=optimizePatchNorm(srcPatch,&tempScore,true);
            //判断是否得到了更大的分数
            if(tempScore>maxScore)
            {
                //记录更大的分数
                maxScore=tempScore;
                //记录最大的分数对应的z的位置
                maxScoreZ=currZ;
                //清空目前位置的点列表
                savedPtList.clear();
                //记录目前状态下的点列表
                srcPatch.getPlanePts(savedPtList);
                //判断分数是否已经达到了阈值,已经达到阈值就不需要再优化了
                if(optFlag) return;
            }
            else
            {
                //恢复以前记录过的状态
                srcPatch.loadPointList(savedPtList,&savedPlane);
            }
      }
      //判断是否需要做进一步的优化
      if(optForward)
      {
          //新建下一次的优化范围
          RangeType nextRange;
          nextRange[0]=maxScoreZ-zStep;
          nextRange[1]=maxScoreZ+zStep;
          //执行下一次优化，但下一次优化就已经是最后一次优化了
          iteratePatchZ(srcPatch,nextRange,maxScore,STEP_L1,false);
      }
  }

  //对矩形面片里面每个位置的z值做优化，寻找不同位置的z值，看一下能不能找到更合适的z信息
  //当调用这个函数的时候，要保证patch里面的平面信息是已经初始化过的
  void optimizePatchZ(DomRectPatch& srcPatch,double maxScore=0)
  {
        //获取z的范围
        const RangeType& zRange=domInfo_.zRange_;
        //计算z的步长
        double zStep=(zRange[1]-zRange[0])*RAND_Z;
        //记录目前的点坐标列表
        Point3DList savedPtList;
        srcPatch.getPlanePts(savedPtList);
        //记录目前的直线方程
        PlaneT savedPlane=srcPatch.getPlaneEqu();
        //初始化opencv的随机数生成器
        cv::RNG randMaker;
        //循环若干次，用于更新z值
        for(unsigned int iterTime=0;iterTime<Z_OPT_TIMES;++iterTime)
        {
            //生成一个随机的z变化量
            double randValue=randMaker.uniform(-zStep,zStep);
            //将随机的z值叠加到旧的z值上,这个过程中需要重新更新平面中的每个点
            //由接口里面的点自动更新
            srcPatch.updatePriZ(randValue);
            //更新当前位置的法向量，并计算它临时的分数
            double tempScore=maxScore;
            //写入true表示不达到阈值也会把颜色记录下来，最终记录最大的颜色值
            bool optFlag=optimizePatchNorm(srcPatch,&tempScore,true);
            //每个z在这里都有一次机会恢复到自己以前的最佳状态
            //但在这里并不恢复平面方程
            srcPatch.loadPointList(savedPtList,nullptr,tempScore);
            //判断是否得到了更大的分数
            if(tempScore>maxScore)
            {
                //记录更大的分数
                maxScore=tempScore;
                //清空目前位置的点列表
                savedPtList.clear();
                //记录目前状态下的点列表
                srcPatch.getPlanePts(savedPtList);
                //判断是否约束成功
                if(optFlag) return;
            }
            else
            {
                //恢复到旧的平面方程，点在上面已经恢复过了
                srcPatch.setPlaneEqu(savedPlane);
            }
        }
  }

  //把dom图里面的每个位置的高度存储成点云
  void getZAsCloud()
  {
      std::cout<<"reload DOM in pointcloud"<<std::endl;
      //清空目前位置的点云
      structure.clear();
      //遍历dom图里面的每个位置
      for(unsigned int rowCount=0;rowCount<domInfo_.domHeight_;++rowCount)
      {
          for(unsigned int colCount=0;colCount<domInfo_.domWidth_;++colCount)
          {
              //获取dom信息
              DomUnit& getUnit=domInfo_.getUnit(colCount,rowCount);
              //没优化出来的点就不用了
              if(getUnit.nccScore_<0) continue;
#if defined (GLOBAL_HELP) || (!defined (USE_Z_PLANE))
              //判断是否为外点
              if(getUnit.isOutRangePt_) continue;
#endif
              //新建一个点对
              std::pair<IndexT,Landmark> tempPair;
              tempPair.first=structure.size();
              //记录点坐标
              domInfo_.convertPixelLocal(colCount,rowCount,tempPair.second.X);
              tempPair.second.X[2]=getUnit.z_;
              //记录点云的颜色
              domInfo_.getDomResultColor(colCount,rowCount,tempPair.second.avgColor);
              //标记点云的颜色已经初始化过了
              tempPair.second.avgInited=true;
              //把点添加到点云里面
              structure.insert(tempPair);
          }
      }
  }

  //把dom图里面每个位置的匹配分数保存成点云
  void getScoreAsCloud(double nccThre=NCC_THRE)
  {
      //清空目前位置的点云
      structure.clear();
      //遍历dom图里面的每个位置
      for(unsigned int rowCount=0;rowCount<domInfo_.domHeight_;++rowCount)
      {
          for(unsigned int colCount=0;colCount<domInfo_.domWidth_;++colCount)
          {
              //获取dom信息
              DomUnit& getUnit=domInfo_.getUnit(colCount,rowCount);
              //判断是否为优化完成的点
              if(getUnit.nccScore_<nccThre||getUnit.isOutRangePt_) continue;
              //新建一个点对
              std::pair<IndexT,Landmark> tempPair;
              tempPair.first=structure.size();
              //记录点坐标
              domInfo_.convertPixelLocal(colCount,rowCount,tempPair.second.X);
              tempPair.second.X[2]=getUnit.z_;
              //把点添加到点云里面
              structure.insert(tempPair);
          }
      }
  }

  void printProcess(std::string infoStr)
  {
      //输出信息
      std::cout<<infoStr<<std::endl;
  }

  //遍历dom图做优化
  //传入的是起始步长
  //第2个参数是true的时候，它会优化到底，否则每个面片的中心点的z不会被改变
  //sliderLen是窗口向下滑动的长度，为了填充稀疏点用的,使用德劳内的z值做辅助的时候，
  //sliderLen，它传入1的时候表示使用德劳内的z值做优化
  //minZStep是允许的最小的z阶梯，为了保证z值较高的点优先重建设置的一个参数
  bool optimizeDomPixel(int bgStep,bool optAnyway=false,int maskSize=MASK_SIZE(0),int sliderLen=MASK_SIZE(0),double nccThre=NCC_THRE,
                        int minZStep=0)
  {
      //dom图的宽度的迭代范围
      unsigned int domWidth=domInfo_.domWidth_-maskSize;
        //初始化是否使用德劳内的z值
      bool useGlobalHelp=false;
      //是否有新的待优化的点
      bool haveUnrefine=false;
#ifdef GLOBAL_HELP
      //根据滑动值是否为1来判断
      if(sliderLen==1)
      {
          //真实需要的滑动值其实还是窗口的宽度，这里仅仅是为了传递一个信息
          sliderLen=maskSize;
          //把使用德劳内z值的标志设为真
          useGlobalHelp=true;
      }
#endif
      //遍历dom图的每一行每一列
      for(unsigned int domRow=bgStep;domRow+maskSize<=domInfo_.domHeight_;domRow+=sliderLen)
      {
          //输出迭代位置
          //std::cout<<domRow<<"/"<<domInfo_.domHeight_<<std::endl;
#ifdef USE_OMP
#pragma omp parallel for
#endif
          for(unsigned int domCol=bgStep;domCol<=domWidth;domCol+=maskSize)
          {
              //生成当前位置的范围数据
              ImgRange iterRange;
              iterRange.initValue(domCol,domCol+maskSize,domRow,domRow+maskSize);
              //新建dom图里面的矩阵块
              DomRectPatch tempPatch;
              //从dom图需要获取矩阵块的每个位置的z优化器
              domInfo_.getRect(iterRange,tempPatch.optiVec_,nccThre);
              //判断面片内有效点的个数是否足够
#if defined (GLOBAL_HELP) || (!defined (USE_Z_PLANE))//使用z平台算法的时候，点的数量一定是够的
              if(tempPatch.optiVec_.size()>=MIN_MASK_PT)
#endif
              {
                  //判断是否存在未优化的点
                  if(!tempPatch.haveUnrefine(nccThre)) continue;
                  //更新确定还有待优化的点
                  haveUnrefine=true;
#ifdef USE_Z_PLANE
                  //判断是否有可使用的点，如果没有可使用的点那就先不管
                  if(tempPatch.initZSteps(domInfo_.zRange_,STEP_L1,useGlobalHelp,nccThre-REF_THRE_DIS,minZStep)==0) continue;
                  //判断当前位置是否需要输出信息
                  if(debugStop(domCol,domInfo_.domHeight_-domRow-1))
                  {
                      //新建一个输出流管理器
                      FileHandler logFile(logFileName(domCol,domInfo_.domHeight_-domRow-1));
                      optimizePlaneZPatch(tempPatch,nccThre,&logFile);
                  }
                  else
                  {
                      //使用z平面方案的专属优化方式
                      optimizePlaneZPatch(tempPatch,nccThre);
                  }
#else
                  //有了z优化器信息后，更新这个优化器的代表点
                  tempPatch.updateRepPoint(nccThre);
                  //完成对矩形面片里面内容的优化
                  //没有被优化成功的z，它上面是不会有标记的
                  optimizeRectPatch(tempPatch,optAnyway);
#endif
              }
          }
      }
      //返回是否还有待优化的点
      return  haveUnrefine;
  }

  //释放view里面的数据
  void releaseImgs()
  {
        //遍历每个view
      for(auto& eachView : views)
      {
          //判断是否只使用中间的图片
#ifdef USE_ONLY_MID_CAM
          if(!MID_ID_JUDGE(eachView.first)) continue;
#endif
          //释放图片
          eachView.second->releaseImg();
          //保存当前位置的z值图
#ifdef USE_IMG_Z_CONSTRAIN
#ifdef SAVE_Z_MAP
          cv::imwrite(SAVE_IMG_PATH(std::to_string(eachView.first)),
                                    eachView.second->imgZMap_);
#endif
#endif
      }
  }

  //获取从DOM像素到相机像素在某个高度值上的单应性变换
  Mat3 getHomographyDomToImage(View& srcView,float zValue)
  {
        //当前的view的内参矩阵
        Mat3 intrMat=intrinsics.at(srcView.id_intrinsic)->getIntrinsicMatrix();
        //旋转矩阵
        Mat3 rotMat=poses.at(srcView.id_pose).rotation();
        //相机光心
        Vec3 camCenter=poses.at(srcView.id_pose).center();
        //后置的变换矩阵，形成单应性矩阵的关键
        Mat3 postFix=Eigen::Matrix3d::Identity(3,3);
        postFix(0,2)=-camCenter[0];
        postFix(1,2)=-camCenter[1];
        postFix(2,2)=zValue-camCenter[2];
        //返回三个矩阵相乘的结果
        return intrMat*rotMat*postFix;
  }

  //判断一个图片能否看到DOM的场景
  bool isViewUsable(View& srcView)
  {
      //判断pose是否存在，如果不存在则返回false
      if(poses.count(srcView.id_pose)==0) return  false;
      //获取当前的view从dom像素到图片像素在某个高度上的单应性变换
      Mat3 homoDom2ImageLow=getHomographyDomToImage(srcView,domInfo_.zRange_[0]);
      //获取最高点平面的单应性变换
      Mat3 homoDom2ImageHigh=getHomographyDomToImage(srcView,domInfo_.zRange_[1]);
      //对涉及到的矩阵取逆
      Mat3 homoImage2DomLow=homoDom2ImageLow.inverse();
      Mat3 homoImage2DomHigh=homoDom2ImageHigh.inverse();
      //遍历原图的四个角判断是否能投影在DOM像素上
      for(IndexT idCorner=0;idCorner<4;++idCorner)
      {
          //获取当前位置的角
          Vec3 cornerPoint;
          cornerPoint[0]=idCorner/2?0:srcView.ui_width;
          cornerPoint[1]=idCorner%2?0:srcView.ui_height;
          cornerPoint[2]=1;
          //对当前的位置做投影
          Vec3 projectPointHigh=MatFunc::homoProject(cornerPoint,homoImage2DomHigh);
          //判断投影点是否在图片上
          if(projectPointHigh[2]!=0 &&
                  srcView.inImgRange(projectPointHigh[0],projectPointHigh[1]))
          {
              //对最低点做投影
              Vec3 projectPointLow=MatFunc::homoProject(cornerPoint,homoImage2DomLow);
              //判断投影点是否在图片上
              if(projectPointLow[2]!=0 &&
                      srcView.inImgRange(projectPointLow[0],projectPointLow[1]))
              {
                  return true;
              }
          }
      }
      return false;
  }

  //判断投影直线在最高点和最低点的投影是否有可能在目标图片上
  bool isProjectLineUsable(ProjectLine& srcProjectLine)
  {
        //目标图片
        View& currView=*views.at(srcProjectLine.imgId_);
        //投影点
        TempPtBase projectPixel;
        //获取最高点的投影位置
        srcProjectLine.getXy(domInfo_.zRange_[1],projectPixel[0],projectPixel[1]);
        //判断投影是否在范围内
        if(!currView.inImgRange(projectPixel[0],projectPixel[1])) return false;
        //获取最低点的投影位置
        srcProjectLine.getXy(domInfo_.zRange_[0],projectPixel[0],projectPixel[1]);
        return currView.inImgRange(projectPixel[0],projectPixel[1]);
  }

  //准备每个DOM像素的投影直线信息,把计算出来的DOM像素存储在外部存储里面
  //如果需要的话，顺便载入图片
  void prepareExternalProjectLine()
  {
      //遍历所有的DOM单元
      for(IndexT idRow=0;idRow<domInfo_.domHeight_;++idRow)
      {
          //当前的标号偏移量
          IndexT idOffset=idRow*domInfo_.domWidth_;
          for(IndexT idCol=0;idCol<domInfo_.domWidth_;++idCol)
          {
              //当前位置的坐标
              Vec3 currLocation;
              domInfo_.convertPixelLocal(idCol,idRow,currLocation);
              //当前DOM像素对应的二进制文件
              std::fstream unitBinaryFile;
                unitBinaryFile.open(projectBinFileName(idOffset+idCol),std::ios::out|std::ios::binary);
                if(!unitBinaryFile.is_open()) throw ERROR_FAIL_OPEN_PROJECT_BIN;
              //遍历所有的图片
              for(auto& eachView : views)
              {
                  //判断这个view是否有外参，并不是每个view都有外参
                  if(poses.count(eachView.second->id_pose)==0) continue;
                  //新建直线信息
                  ProjectLine tempLine(poses.at(eachView.second->id_pose).rotation(),//旋转矩阵
                            poses.at(eachView.second->id_pose).center(),//光心
                            intrinsics.at(eachView.second->id_intrinsic)->getParams(),//内参
                            currLocation[0],currLocation[1],eachView.first);//dom图上的坐标位置
                  //判断这个投影直线对当前图片是否有价值
                  if(isProjectLineUsable(tempLine))
                  {
                      //载入当前的图片
                      eachView.second->loadImg(s_root_path);
                      //把当前的投影直线记录下来
                     tempLine.saveBinaryFile(unitBinaryFile);
                  }
              }
              //关闭二进制文件
              unitBinaryFile.close();
          }
      }
  }

  //动态载入图片
  void prepareDynamicLoadImage()
  {
      imgNum_=0;
      //遍历所有的view
      for(auto& eachView : views)
      {
          //判断当前的view是否可用
          if(isViewUsable(*eachView.second))
          {
              //std::cout<<"load image:"<<eachView.first<<std::endl;
              imgNum_++;
              //载入图片
              eachView.second->loadImg(s_root_path);
          }
      }
  }

  //根据属性里面保存的点云信息，解析点云的范围
  //点云范围的格式样例："-20 120 0 50"
    void parseRangeInfo(RangeType& xRange,RangeType& yRange)
    {
        std::stringstream tempStream;
        tempStream<<cloudRangeInfo_;
        tempStream>>xRange[0]>>xRange[1]>>yRange[0]>>yRange[1];
    }

  //根据关注的窗口筛选点云
  //仅仅用于debug，做局部区域的处理
  //广东数据的原始范围x:-340,230 y:-510,-60,z:-271,-242
  //广东数据的x:-150~ 50 y: -350~ -150 是房子所在的区域
  //广东数据的x:-150~-50 y:-450~-350是用来比较速度的庄稼区域
  void focusCloudPt()
  {
      //判断是否存在有效的点云信息
      RangeType xRange,yRange;
      if(cloudRangeInfo_.size())
      {
            parseRangeInfo(xRange,yRange);
            domInfo_.makeResolution(xRange[0],xRange[1],yRange[0],yRange[1],structure.size());
            domInfo_.initDom();
            //遍历每个点云 把范围外的点删除
            for(Landmarks::iterator iter=structure.begin();iter!=structure.end();)
            {
                //获取当前位置的点
                Vec3& thisPt=iter->second.X;
                //判断范围
                if(thisPt[0]<xRange[0]) iter=structure.erase(iter);
                else if(thisPt[0]>xRange[1]) iter=structure.erase(iter);
                else if(thisPt[1]<yRange[0]) iter=structure.erase(iter);
                else if(thisPt[1]>yRange[1]) iter=structure.erase(iter);
                else ++iter;
            }
            //同样需要用点云来初始化，但不需要初始化网格了
            domInfo_.initDom(structure,false);
      }
      else
      {
          //否则用点云来做初始化
          domInfo_.initDom(structure);
      }
  }

  //根据要求的图片和dom分辨率的比例更改图片的分辨率
  //处理后，每个图片和dom图的分辨率的比值是dstRate
  void changeImgResolution(double dstRate)
  {
      //初始化的标志
      bool initFlag=true;
      //图片的缩放比例
      double scaleRate=-1;
      //遍历每个图片
      for(auto& eachView : views)
      {
#ifdef USE_ONLY_MID_CAM
          if(!MID_ID_JUDGE(eachView.first)) continue;
#endif
          //如果没有载入图片就放弃
          if(eachView.second->cvImg_.empty()) continue;
          //判断是否为第1次遍历
          if(initFlag)
          {
              //计算图片需要的缩放比例
              scaleRate=std::sqrt(dstRate*
                                  domInfo_.getResolution()/eachView.second->imgResolution());
              std::cout<<"scale rate:"<<scaleRate<<std::endl;
              //取消初始化的标志
              initFlag=false;
          }
          //更改图片的分辨率
          eachView.second->resizeScale(scaleRate);
      }
  }

  //中间位置相机的过滤
  //当算法中选择只使用中相机的时候，在这里先把中相机拍摄不到的点去除
  void midCameraFilter()
  {
      //遍历所有的稀疏点云
      for(Landmarks::iterator iter=structure.begin();
          iter!=structure.end();)
      {
          //观察点中是否存在中相机的标志
          bool haveMidCam=false;
          //遍历能看到当前点的每个相机
          for(auto& eachObv : iter->second.obs)
          {
              //判断当前位置的相机标号是否符合要求
              if(MID_ID_JUDGE(eachObv.first))
              {
                  //标记它具有中相机的观察者
                  haveMidCam=true;
                  break;
              }
          }
          //如果有中相机的观察者，就保留当前的稀疏点
          //如果没有的话就删除当前位置的点，删除后返回被删除位置的下一个位置
          if(haveMidCam) ++iter;
          else iter=structure.erase(iter);
      }
  }


  //存储sfm_data的点云
  void saveDomCloud(const std::string& filePath//存储的路径
                  )
  {
      //把DOM信息保存成cloud
      getZAsCloud();
      //保存sfm_data的点云信息
      openMVG::sfm::Save(*this,filePath,openMVG::sfm::ESfM_Data::STRUCTURE);
  }


  //20210719 与下面denseDomByDomPixel的算法类似
  //在下面那个算法的基础上结合了mvs论文里面的方法
  //如果选择onlyTestRange,则仅仅显示一下运行的范围，什么都不做
  void denseDomLikeMvs(bool onlyTestRange=false)
  {
      //查找关注位置的点云
        focusCloudPt();
        //判断是否需要删除中间位置
#ifdef USE_ONLY_MID_CAM
        midCameraFilter();
#endif
        if(onlyTestRange) return;
#ifdef SAVE_EXTERNAL_PROJECT_LINE
        prepareExternalProjectLine();//准备每个DOM像素的投影直线
#endif
        //先把已经有的稀疏点云的信息部署在dom图上
        drawSparsPtAtDom();
        //删除没有载入过的相机光心
        DeleteNonLoadedViews nonLoadDeleter;
        nonLoadDeleter(views);

        //如果需要缩放图片的比例的话，就在这里缩放,万万不可在这之前缩放
#ifdef OVERLAP_RATE
        changeImgResolution(OVERLAP_RATE(imgNum_));
#endif
        //根据目前的dom重新载入点云
        structure.clear();
#if defined (GLOBAL_HELP) || (!defined (USE_Z_PLANE))
        //做德劳内三角化，用三角化的结果初始化每个位置的z值
        {//带个大括号是为了让它早点释放掉内存
            DomPatchList patchList;
            makeAdvancedTriPatch(patchList);
            //用三角化的结果初始化dom里面的z值
            domInfo_.initGlobalZValue(patchList);
        }
#endif
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
                if(haveUnRefine)//必须是有新的优化成果再做优化
                {
                    //计算时间差
                    unsigned int seconds=time(nullptr)-bgTime;
                    //转换为需要记录的时间
                    std::string fileName=std::to_string(seconds/60)+"-"+std::to_string(seconds%60);
                    saveDomResult(midProceeFolder_+"/"+fileName+".bmp");
                    //临时保存点云，正常运行的时候不需要动这个地方
                    //saveDomCloud(midProceeFolder_+"/"+fileName+".ply");
                }
            }
#endif
            //当前循环的面片大小
            int maskSize=MASK_SIZE(iterTimes);
            //更新当前周期是否需要考虑遮挡问题
            isUsingOccludeDetection_=isIterateDetectOcclude(iterTimes);
            std::cout<<iterTimes<<std::endl;
            //初始化是否有可优化的点
            haveUnRefine=false;
            for(int currStep=0;currStep<maskSize;currStep+=MASK_STEP)
            {
                 //根据目前的步长优化dom图的每个像素
#ifdef FIX_MASK_TIMES//同大小窗口固定迭代若干次
                //临时的结果
                bool resultFlag=optimizeDomPixel(currStep,false,maskSize,1,THRE_ID(iterTimes),ALLOW_MIN_Z(iterTimes));
                //判断是否有未优化的点
                if(resultFlag) haveUnRefine=true;
#endif
            }
            //判断是否需要更新dom像素的可用次数
            if(NEED_INIT_REF(iterTimes))
            {
                domInfo_.initRefTime();
            }
            //如果没有待优化的点了结束循环
            if(iterTimes>CONSIDER_BREAK_TIME && haveUnRefine==false)
                break;
            //判断是否强制结束
            //if(forceEndJudge()) break;
        }
#ifdef END_WITH_OCCLUDE_DETECTION
        //更新可优化次数
        domInfo_.initRefTime();
        //锁定Z值
        holdZValue_=true;
        //使用遮挡判断
        isUsingOccludeDetection_=true;
        optimizeDomPixel(0,false,3,1,1,0);
#endif
        //做最后一次优化，寻找还是没有被优化成功的点
        //optimizeDomPixel(0,true);
  }

  //依次遍历dom平面上的每个位置
  //依次填充dom图上每个像素的颜色值
  void denseDomByDomPixel()
  {
        //初始化dom图
      domInfo_.initDom(structure);
      //z范围的长度
      double zLength=domInfo_.zRange_[1]-domInfo_.zRange_[0];
      //遍历z的时候合适的范围
      std::array<double,2> zRange={domInfo_.zRange_[0]-zLength*SEARCH_UP_RATE
                                  ,domInfo_.zRange_[1]+zLength*SEARCH_UP_RATE};
      //步长
      double zStep=zLength*STEP_L1;
        //遍历dom图上的每个位置
      for(unsigned int rowCount=0;rowCount<domInfo_.domHeight_;++rowCount)
      {
          std::cout<<rowCount<<"/"<<domInfo_.domHeight_<<std::endl;
#ifdef USE_OMP
#pragma omp parallel for
#endif
          for(unsigned int colCount=0;colCount<domInfo_.domWidth_;++colCount)
          {
              //获取当前位置的像素信息
              DomUnitAll allInfo;
              //正常情况下是不会continue的
              if(!domInfo_.getUnit(colCount,rowCount,allInfo)) continue;
              //获取当前位置下的所有直线
              std::vector<ProjectLine> projLineList;
              getProjLines(projLineList,allInfo.scenePt_[0],allInfo.scenePt_[1]);
              //在各个直线上寻找合适的颜色
              getFitColorZ(projLineList,allInfo,zRange,zStep);
          }
      }
  }

  //去除稀疏点云里面误差比较大的点
  void sparseFilter()
  {
      //获取点云的范围
      RangeT cloudRange=DomInfo::getCloudPointRange(structure);
      //计算平均每个点云占用的面积
    double avgArea=(cloudRange[1]-cloudRange[0])*(cloudRange[3]-cloudRange[2])/structure.size();
    //根据平均每个点云占用的面积，计算出单位点云占用的长度,这在下面是阈值的重要参考标准
    double avgLength=std::sqrt(avgArea);
      //遍历点云里面的点,由于里面涉及删除操作，因此需要在循环体里面自行移位
      for(Landmarks::iterator iter=structure.begin();iter!=structure.end();)
      {
            //获取当前的点云坐标
          const Vec3& cloudEigenPt=iter->second.X;
          CornerPt cloudPt(iter->first,cloudEigenPt[0],cloudEigenPt[1]);
          //找一下边缘点
          if(cloudPt[0]-cloudRange[0]<=avgLength)
          {
              std::cout<<cloudPt[0]<<" "<<cloudPt[1]<<std::endl;
          }
          //能看到这个点的所有图
          Observations& obvImgs=iter->second.obs;
          double errorSum=0;
          //累积计算它在每个图上的误差
          for(Observations::iterator imgIter=obvImgs.begin();imgIter!=obvImgs.end();++imgIter)
          {
              //当前位置的图片投影坐标
              const TempPtBase &imgPt=imgIter->second.x;
              //计算当前坐标在dom平面上的投影
              LineMid domLine=getDomLine(imgIter->first,imgPt);
              //累加点到投影直线的距离
            errorSum+=domLine.pointDistance(cloudPt);
          }
          //计算平均误差
          double avgError=errorSum/obvImgs.size();
          //std::cout<<"average error: "<<avgError<<"\t"<<"average len: "<<avgLength<<std::endl;
          //判断点到直线距离的平均值是否超过了阈值
          if(avgError>avgLength)
          {
              //删除当前位置的点云并指向下一个位置的点云
              iter=structure.erase(iter);
          }
          else
          {
              //单纯地移位
              ++iter;
          }
      }
  }//去除稀疏点云里面误差比较大的点

  //确保每个点云都有至少3个有效的观察者
  void filterLossObv()
  {
      //足够分离的角度范围暂且判断为10度
      static const double angleThrehold=(double)10/180/M_PI;
      //每个点最少要有3个有效观测图
      static const int leastObvCount=4;
      //遍历点云里面的每个点,迭代器的移位在循环体内部进行
      for(Landmarks::iterator iter=structure.begin();iter!=structure.end();)
      {
            //当前位置的点云坐标
          const Vec3 &cloudEigenPt=iter->second.X;
          //生成二维的点坐标
          TempPtMid cloudPt(cloudEigenPt[0],cloudEigenPt[1]);
          //能看到这个点的所有图
          Observations &obvs=iter->second.obs;
          //是否保留这个点云点
          bool saveFlag=true;
          //判断看到这个点的数量是否足够
          if(obvs.size()<leastObvCount)
          {
              saveFlag=false;
          }
          else //判断包含这个点的图片的dom直线的角度是否足够分离
          {
                //每个投影点生成的dom直线的角度情况
              std::map<double,bool> angleMap;
              //遍历每个观察图片，然后记录它们的角度分布
              for(Observations::iterator imgIter=obvs.begin();imgIter!=obvs.end();++imgIter)
              {
                  //当前位置的点坐标
                  const TempPtBase& imgPt=imgIter->second.x;
                  //当前投影点对应的dom直线
                  LineMid domLine=getDomLine(imgIter->first,imgPt);
                  //记录当前位置的角度数据
                    angleMap[domLine.absoluteAngle()]=true;
              }
              //有效的观测数据的统计,并不是一个很精妙的算法，但基本够用
              unsigned int obvCount=0;
              //上一次的统计角度
              double lastAngle=0;
              //遍历统计到的角度信息
              for(std::map<double,bool>::iterator mapIter=angleMap.begin();
                  mapIter!=angleMap.end();++mapIter)
              {
                  //判断是否存在记录过的角度
                  if(obvCount==0)
                  {
                      //记录角度然后下一次
                      obvCount++;
                      lastAngle=mapIter->first;
                      continue;
                  }
                  //和第1个角度太近也不行
                  if(angleMap.begin()->first-mapIter->first+M_PI<angleThrehold)
                  {
                      //说明已经到头儿了，不用再统计了
                      break;
                  }
                  //判断当前角度与上一次角度是否超过阈值
                  if(mapIter->first-lastAngle>angleThrehold)
                  {
                      //记录新的角度
                      obvCount++;
                      lastAngle=mapIter->first;
                  }
              }
              //判断有效的信息点是否足够多
               if(obvCount<leastObvCount)
                   saveFlag=false;
          }
          //判断是否保存当前的点
          if(saveFlag)
          {
              iter++;
          }
          else //不保存的情况下删除当前迭代器所在位置
          {
              iter=structure.erase(iter);
          }
      }
  }//void filterLossObv() 确保每个点云至少有3个观察者

  //计算某个特征点的平均颜色，投影到具体的图上然后看它们的颜色
  //最后计算出来的平均颜色会被存储到这个avgColor里面
  //当needRecordZ为true的时候，记录每次的投影结果，并且把z值记录到图片上
  //返回有效颜色的个数
  uint pointAvgColor(Landmark& cloudPt, UserColor& avgColor,bool needRecordZ=false)
  {
      //初始化颜色
      avgColor*=0;
      //有效颜色的计数
      uint colorCount=0;
      //遍历能看到这个点的每个图片
      for(auto &eachObv : cloudPt.obs)
      {
          //判断是否只使用中间的相机
#ifdef USE_ONLY_MID_CAM
          if(!MID_ID_JUDGE(eachObv.first)) continue;
#endif
          //判断图片是否载入过
          views.at(eachObv.first)->loadImg(s_root_path);
          //当前位置
          const Vec2 &thisLocal=eachObv.second.x;
          //读取颜色
          try {
              cv::Vec3b &thisColor=views.at(eachObv.first)->colorAt(thisLocal[0],thisLocal[1]);
              //用userColor来初始化这个颜色
              UserColor tempColor(thisColor);
              //把颜色叠加到平均颜色上,这个+=是自己单独实现的，不是opencv自带的东西
              avgColor+=tempColor;
              //判断是否需要记录z值
              if(needRecordZ)
              {
                  //计算临时的z所处的阶层
                  int tempZLevel=domInfo_.getZLevel(cloudPt.X[2]);
                  //把z所在的阶层添加到目标的像素位置上
                  views.at(eachObv.first)->registerZMap(thisLocal[0],thisLocal[1],tempZLevel);
              }
              //记录一个有效颜色
              colorCount++;
          } catch (int errorFlag) {
              //的可能出现恰好在边缘的情况，这种情况下直接就不要了
              if(errorFlag==ERROR_IMG_OUT_RANGE) continue;
              else
                  throw errorFlag;
          }
      }
      //除以总数获取平均值
      avgColor.VectorInterface::operator*=(1.f/colorCount);
      //avgColor.Vec::operator=(cv::Vec3b(255,255,255));
      //返回有效颜色的个数
      return colorCount;
  }

  //制作一体机需要使用的pose文件
  //这里制作的是直接的坐标数据，而下面制作的是经纬度坐标
  void makeAltizurePose(std::string poseFile)
  {
      //输出流
      std::fstream fileHandle;
      fileHandle.open(poseFile,std::ios::out);
      std::cout<<views.size()<<std::endl;
      //文件的开头 表示使用局部坐标系
      fileHandle<<"coordinatesystem local"<<std::endl;
      //遍历所有的view
      for(auto& eachView : views)
      {
          //判断图片是否载入过，图片载入条件会在调用前处理
          if(eachView.second->cvImg_.empty()) continue;
          //当前位置的pose
          Pose3& currPose=poses.at(eachView.second->id_pose);
          //光心坐标
          Vec3 camCenter=currPose.center();
          //输出信息
          fileHandle<<eachView.second->s_Img_path<<" "
                   <<camCenter[0]<<" "<<camCenter[1]<<" "<<camCenter[2]<<std::endl;
      }
      //关闭
      fileHandle.close();
  }

  //获取图片的pose
  //这个函数是用来针对一体机制作pose接口使用的
  void makePose(std::string poseFile)
  {
      //每个纬度的长度
      double latiDis=111000;
      //纬度中心
      long double latiCenter=34.58000000000;
      //经度中心
      long double longtiCenter=116.90000000;
      //高度中心
      double heightCenter=132.6;
      //一个经度的长度
      double longtiDis=latiDis*std::cos(latiCenter/180*M_PI);
      //打开pose文件
      std::ofstream poseHandle;
      poseHandle.open(poseFile,std::ios::out);
      poseHandle.precision(7);
      //遍历所有的view
      for(const auto& eachView : views)
      {
#ifdef USE_ONLY_MID_CAM
          if(!MID_ID_JUDGE(eachView.first)) continue;
#endif
          //获取当前位置的坐标
          const Vec3& thisPose=poses.at(eachView.second->id_pose).center();
          //记录文件名
          poseHandle<<eachView.first<<".JPG ";
          //记录坐标
          poseHandle<<latiCenter+(long double)thisPose[1]/latiDis<<" "
                    <<longtiCenter+(long double)thisPose[0]/longtiDis<<" "
                 <<thisPose[2]+heightCenter<<"\n";
      }
      //关闭文件
      poseHandle.close();
  }

  //一体机做三维重建还需要用到的group的信息
  void makeAltizureGroup(std::string groupFile)
  {
      std::fstream fileHandle;
      fileHandle.open(groupFile,std::ios::out);
      //遍历每个view
      for(auto& eachView : views)
      {
          //判断view是否载入过
          if(eachView.second->cvImg_.empty()) continue;
          //记录group内容
          fileHandle<<eachView.second->s_Img_path<<" "<<eachView.first/100000<<std::endl;
      }
      //关闭文件
      fileHandle.close();
  }

  //制作相机的内参信息
  //针对一体机制作输入数据接口用的
  void makeIntrinsicInfo(std::string camFile)
  {
      //打开输出流文件
      std::ofstream camHandle;
      camHandle.open(camFile,std::ios::out);
      //遍历所有的图片
      for(const auto& eachView : views)
      {
          //如果是没有载入过的图片就不要了
          if(eachView.second->cvImg_.empty()) continue;
#ifdef USE_ONLY_MID_CAM
          if(!MID_ID_JUDGE(eachView.first)) continue;
#endif
          //当前位置的内参
          std::vector<double> tempIntr=intrinsics.at(eachView.second->id_intrinsic)->getParams();
          //记录图片名称
          camHandle<<eachView.first<<".JPG ";
          //依次记录三个参数
          camHandle<<tempIntr[0]<<" "<<tempIntr[1]<<" "<<tempIntr[2]<<"\n";
      }
      //关闭文件
      camHandle.close();
  }

  //把某个点云中的点投影到图片上，记录它在图片上的投影位置
  Vec2 projectCloudPoint(const Vec3& cloudPoint,//被投影的点云
                         IndexT imgId //把点云投影到这个图片上
                         )
  {
        //把点转换到相机坐标系下,用pose里面的一个成员函数来实现
      Vec3 xcPoint=poses.at(views.at(imgId)->id_pose).toCameraCoord(cloudPoint);
      //把它投影到相机坐标系下,需要调用内参接口
        return intrinsics.at(views.at(imgId)->id_intrinsic)->projectXcPoint(xcPoint);
  }

  //把某个点云中的点投影到图片上，记录它在图片上的投影位置
  Vec2 projectCloudPoint(const Landmark& cloudPoint,//被投影的点云
                         IndexT imgId //把点云投影到这个图片上
                         )
  {
       return projectCloudPoint(cloudPoint.X,imgId);
  }

  //对场景做了稀疏的三维重建后，让不同的view之间共享特征点
  //毕竟寻找特征点的时候，有的特征点虽然被拍到了，但并没有被所有的图判断为特征点
  //让view之间可用的特征点共享后，方便后面做DOM稠密化的时候有更多的素材可以使用
    void shareViewFeature()
    {
#ifdef DEBUG_PRINT
        std::cout<<"point account: "<<structure.size()<<std::endl;
#endif
        //遍历所有的特征点
        for(auto& eachPt : structure)
        {
            //计算当前位置的平均颜色
            UserColor avgColor;
            pointAvgColor(eachPt.second,avgColor);
            //当前点的观测信息
            Observations &obvs=eachPt.second.obs;
#ifdef DEBUG_PRINT
            std::cout<<"obs: "<<obvs.size()<<" ";
#endif
            unsigned int giveUpCount=0;
            //遍历每个view
            for(auto& eachView : views)
            {
                //判断obv里面是否已经有了这个图片
                if(obvs.count(eachView.first)>0) continue;
                //判断这个相机是否存在位姿，并不是所有的图都存在位姿,有的位姿可能是约束不出来的
                if(poses.count(eachView.second->id_pose)==0) continue;
                //载入view的图片
                eachView.second->loadImg(s_root_path);
                //用于添加到obv里面的pair数据
                std::pair<IndexT,Observation> tempPair;
                tempPair.first=eachView.first;
                //复制投影点，用于后续操作
                Vec2& imgProj=tempPair.second.x;
                //当前点云在该图片上的投影点
                imgProj=projectCloudPoint(eachPt.second,eachView.first);
                //判断投影点是否在图片范围内
                if(!eachView.second->inImgRange(imgProj[0],imgProj[1])) continue;
                //获取当前图片在投影位置的颜色
                cv::Vec3b viewColor=eachView.second->colorAt(imgProj[0],imgProj[1]);
                //构造颜色数据
                UserColor userColor(viewColor);
                //判断颜色距离是否足够小,注意这里是余弦距离，先减去它们的一个先验的均值
                double cosDistance=avgColor.cosDis(userColor);
                if(cosDistance<COL_THRE)
                {
                    ++giveUpCount;
                    continue;
                }
                //检验通过，在哈希表里面添加这个投影数据
                obvs.insert(tempPair);
            }
#ifdef DEBUG_PRINT
            std::cout<<"give up point: "<<giveUpCount<<std::endl;
#endif
        }

    }
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_DATA_HPP
