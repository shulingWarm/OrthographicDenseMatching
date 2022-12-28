// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_LANDMARK_HPP
#define OPENMVG_SFM_SFM_LANDMARK_HPP

#include "openMVG/types.hpp"
#include<iostream>
#include<opencv2/opencv.hpp>
#include"switch.hpp"
#include"kmeans.hpp"
#include<Eigen/Eigen>
#include<Eigen/Dense>

namespace openMVG {
namespace sfm {

/// Define 3D-2D tracking data: 3D landmark with it's 2D observations
struct Observation
{
  Observation():id_feat(UndefinedIndexT) {  }
  Observation(const Vec2 & p, IndexT idFeat): x(p), id_feat(idFeat) {}

  Vec2 x;
  Vec3 xcPoint;
  IndexT id_feat;

  // Serialization
  template <class Archive>
  void save( Archive & ar) const;

  // Serialization
  template <class Archive>
  void load( Archive & ar);
};

/// Observations are indexed by their View_id
using Observations = Hash_Map<IndexT, Observation>;

/// Define a landmark (a 3D point, with it's 2d observations)
struct Landmark
{
  Vec3 X;
  Observations obs;

  //所属的三角面片序列
    std::vector<IndexT> patchIdList_;

    //该点的平均颜色
    cv::Vec3i avgColor;
    bool avgInited=false;

    //依次把每个观察点的坐标都打印一遍
    void printObvs()
    {
        for(auto& eachObv : obs)
        {
            std::cout<<eachObv.first<<" "<<eachObv.second.x<<std::endl;
        }
    }

  // Serialization
  template <class Archive>
  void save( Archive & ar) const;

  template <class Archive>
  void load( Archive & ar);
};

/// Define a collection of landmarks are indexed by their TrackId
using Landmarks = Hash_Map<IndexT, Landmark>;
//点云信息的顺序容器，哈希表遍历起来有点不方便
typedef std::vector<Landmark*> LandmarkPtrList;
//三个点的组合
typedef std::array<IndexT,3> Pt3Group;
typedef std::map<IndexT,bool> Pt3Map;
//标号的列表
typedef  std::vector<IndexT> IndexList;
//用来做统计计数的哈希表
typedef std::unordered_map<IndexT,IndexT> CountMap;
//直线方程的数据类型
typedef Eigen::Vector3d LineT;
//平面方程的数据类型
typedef Eigen::Vector4d PlaneT;
//齐次点的类型和直线保持一致
typedef LineT HomoPt;
//临时点的基类，typedef出来方便以后改,之前本来用的是std::array,我真傻，真的
typedef Eigen::Vector2d TempPtBase;
//opencv里面点的数据类型
typedef cv::Point2i CvPoint;
//opencv里面点数据类型的列表
typedef std::vector<CvPoint> CvPointVector;
//三维点的基类，同时还可以继承一个向量的接口
typedef Eigen::Vector3d Point3DBase;
//表示颜色的数据类型,里面的颜色固定是按照bgr来存储的
typedef cv::Vec3i ColorT;
//基本的颜色数据类型的向量
typedef std::vector<ColorT> ColorTVec;
//目前range使用的数据类型
typedef std::array<double,2> RangeType;
//用于临时表示主点坐标的数据类型
typedef std::array<double,2> PriPoint;
//适用于图片的数据类型，感觉和下面的range2d特别重复
//但当初没有使用泛型，现在只能凑活用了
//分别表示[左右)[上下),都是左闭右开区间
typedef std::array<int,4> ImgRangeBase;
//颜色的串行数据序列，每3个为一组颜色，这个粗犷的数据存储方式方便传输s
typedef std::vector<double> ColorDataSerial;
//先验点可参考次数的数据类型，为了节约内存，暂时使用char类型
typedef uchar RefType;
typedef Eigen::Vector3d EigenVec3;
typedef Eigen::Matrix3d EigenMat3;

//矩阵运算的相关接口
//里面只有静态函数
class MatFunc
{
public:
    //对一个向量做单应性变换并返回结果
    //如果计算出来的位置不是无穷远点，会把它的最后一项弄成1,否则什么都不做
    static EigenVec3 homoProject(EigenVec3& srcVector,EigenMat3& homoMat)
    {
        //计算相乘结果
        EigenVec3 finalAns=homoMat*srcVector;
        //判断是否为无穷远点
        if(finalAns[2]!=0)
        {
            finalAns/=finalAns[2];
        }
        //返回结果
        return  finalAns;
    }
};

//两两互相操作的接口
//本来应该单独放到一个文件里面的
//整个类就实现一个功能，传入一个向量，然后把向量里面的数据两两互相运算
//得到一个特定的结果
//这样的一个功能本来是可以通过做一个函数，然后通过函数传递的方式来做的
//但函数传递的方式写起来感觉还是太乱了，而且我还不会lamda表达式
template<typename SrcT,typename DstT>
class EachOperateInterface
{
protected:
    //分别对应第1层迭代和第2层迭代,方便子类访问必要的信息
    unsigned int vecCount=0;
    unsigned int vec2Count=0;

    //总共的数据个数
    int dataAccount_=0;
public:
    //根据两个位置的数据计算目标的结果
    virtual DstT compute(const SrcT& src1,const SrcT& src2)=0;


    //运行两两运算
    void run(const std::vector<SrcT>& srcVec,std::vector<DstT>& dstVec)
    {
        //判断传入的原始向量个数是否正常
        if(srcVec.size()<2) return;
        //记录数据的个数
        dataAccount_=srcVec.size();
        //提前给目标向量开辟空间
        dstVec.reserve(srcVec.size()*(srcVec.size()-1)/2);
        //遍历向量里面的数据，做两两运算
        for(vecCount=0;vecCount<srcVec.size();++vecCount)
        {
            //当前位置的数据
            const SrcT& thisData=srcVec[vecCount];
            //遍历后面的每个位置的数据
            for(vec2Count=vecCount+1;vec2Count<srcVec.size();++vec2Count)
            {
                //记录当前位置的两两计算结果
                dstVec.push_back(compute(thisData,srcVec[vec2Count]));
            }
        }
    }
};

//20211030为了研究纯色区域传播慢的问题
//特别设定的一个输出流的管理工具
class FileHandler
{
private:
    std::ofstream file_;
public:
    //构造函数，读取文件
    FileHandler(std::string filePath)
    {
        file_.open(filePath,std::ios::app);
    }

    //判断输出流是否可用
    inline bool isFileValid() const
    {
        return file_.is_open();
    }

    //获取关于先验点的信息
    void putPriorInfo(Point3DBase pt3D,double currScore,IndexT currId)
    {
        //判断输出流是否打开
        if(!isFileValid()) return;
        //记录坐标值
        file_<<currId<<"\t";
        for(IndexT dimCount=0;dimCount<3;++dimCount)
            file_<<pt3D[dimCount]<<"\t";
        //记录ncc分数
        file_<<"ncc:"<<currScore<<std::endl;
    }

    //记录颜色信息
    void putColorInfo(IndexT imgId,//图片标号
                      const std::vector<TempPtBase>& projList,//投影点的标号
                      const ColorTVec& colorList)//颜色的列表
    {
        //判断输出流是否可用
        if(!isFileValid()) return;
        //颜色的个数和投影点的个数必须一致
        if(projList.size()!=colorList.size())
        {
            std::cout<<"FileHandle putColorInfo: project list and color not equal"<<std::endl;
            throw;
        }
        //记录图片标号
        file_<<"----------------"<<imgId<<"----------------"<<std::endl;
        //遍历所有的投影点
        for(IndexT ptCount=0;ptCount<projList.size();++ptCount)
        {
            //记录投影点
            const TempPtBase& currPt=projList[ptCount];
            file_<<"("<<currPt[0]<<","<<currPt[1]<<") [";
            //记录颜色
            const ColorT& currColor=colorList[ptCount];
            file_<<currColor[0]<<","<<currColor[1]<<","<<currColor[2]<<"]\n";
        }
    }

    //记录一次ncc的匹配分数
    void putNccScore(double nccScore)
    {
        file_<<"ncc match score:"<<nccScore<<std::endl;
    }

    //记录预先选定的图片标号和实际找到的图片标号
    void putSelectedImageIndex(IndexT adviceId,IndexT realId)
    {
        file_<<"advice id:"<<adviceId<<std::endl;
        file_<<"real id:"<<realId<<std::endl;
    }

    //记录平均颜色
    void putAvgColor(const ColorTVec& colorList)
    {
        //判断是否可用
        if(!isFileValid()) return;
        //平均颜色的信息
        file_<<"-----------------average patch----------------"<<std::endl;
        //记录平均颜色
        for(IndexT idCount=0;idCount<colorList.size();++idCount)
        {
            //当前位置的颜色
            const ColorT& currColor=colorList[idCount];
            file_<<"["<<currColor[0]<<","<<currColor[1]<<","<<currColor[2]<<"]\n";
        }
    }

    //记录中心向量和每个图片向量的余弦距离
    void putCosDistance(IndexT imgId,double cosDis)
    {
        file_<<"cos distance with "<<imgId<<": "<<cosDis<<std::endl;
    }

    //传入一个容器，记录这个容器的向量信息
    //向量必须支持使用中括号做随机访问
    template<typename VecT>
    void recordVector(const VecT& vecInfo,IndexT vecLen,std::string headStr)
    {
        //写入标头信息
        file_<<"*********vector:"<<headStr<<std::endl;
        //遍历向量的每个值
        for(IndexT vecCount=0;vecCount<vecLen;++vecCount)
        {
            //写入当前位置的数值
            file_<<vecInfo[vecCount]<<" ";
        }
        //记录一个换行
        file_<<"\n";
    }

    //记录投影直线的范围
    void recordProjectRange(IndexT imageIndex,double projectRange)
    {
        file_<<"range: "<<imageIndex<<" "<<projectRange<<std::endl;
    }

    //析构函数
    virtual ~FileHandler()
    {
        //如果输出流是可用的，就释放掉它
        if(isFileValid()) file_.close();
    }
};

//一维的范围数据
class Range1D
{
private:
    //0是最小值 1是最大值
    RangeType rangeData_;//数据的范围
public:
    //初始化范围
    void initRange(double value)
    {
        //初始化范围
        rangeData_[0]=value;
        rangeData_[1]=value;
    }

    Range1D(double value=0)
    {
        initRange(value);
    }

    //往范围里面添加数据
    void addData(double value)
    {
        //判断是否小于最小值
        if(value<rangeData_[0]) rangeData_[0]=value;
        if(value>rangeData_[1]) rangeData_[1]=value;
    }

    //让传入的数据适应这个范围
    //如果传入的数据超过了这个范围，就让它等于一个最接近的值
    void adjustDataInRange(double& otherValue)
    {
        if(otherValue<rangeData_[0]) otherValue=rangeData_[0];
        if(otherValue>rangeData_[1]) otherValue=rangeData_[1];
    }
};

//表示坐标的范围
class Range2D
{
protected:
    //都是默认左小右大
    RangeType xRange_;
    RangeType yRange_;
public:
    //初始化坐标范围的比较
    void initRange(double x,double y)
    {
        xRange_.fill(x);
        yRange_.fill(y);
    }

    void addRange(RangeType& dstRange,double value)
    {
        if(value>dstRange[1]) dstRange[1]=value;
        if(value<dstRange[0]) dstRange[0]=value;
    }

    //添加坐标用以修正范围
    void addPoint(double x,double y)
    {
        //修正两个坐标范围
        addRange(xRange_,x);
        addRange(yRange_,y);
    }

    //判断一个数字是否在某个位置内
    bool judgeRange(const RangeType& range,double value) const
    {
        return value<=range[1] && value>=range[0];
    }

    //判断一个坐标是否在范围内
    bool judgeIn(double x, double y) const
    {
        return judgeRange(xRange_,x) && judgeRange(yRange_,y);
    }

    //无参数的构造函数，这里是为了保证通用性
    Range2D(){}

    //使用一系列点坐标来初始化这个范围,保留这个接口是为了防止已经写过的代码出问题
    //后面再有别的问题都去使用那个泛型类
    void initRange(const std::vector<TempPtBase>& ptList)
    {
        //判断传入的点是不是空的
        if(ptList.size()==0) throw ERROR_EMPTY_PT_LIST_MAKING_RANGE;
        //取出第1个点来做初始化
        const TempPtBase& firstPt=ptList[0];
        initRange(firstPt[0],firstPt[1]);
        //遍历剩下的点
        for(unsigned int ptId=1;ptId<ptList.size();++ptId)
            addPoint(ptList[ptId][0],ptList[ptId][1]);
    }

    //带模板的初始化坐标列表的方法
    //传入的点的模板类型必须满足可以使用[]来随机访问
    template<typename PointType>
    void initRangeT(const std::vector<PointType>& ptList)
    {
        //判断传入的点是不是空的
        if(ptList.size()==0) throw ERROR_EMPTY_PT_LIST_MAKING_RANGE;
        //取出第1个点来做初始化
        const PointType& firstPt=ptList[0];
        initRange(firstPt[0],firstPt[1]);
        //遍历剩下的点
        for(unsigned int ptId=1;ptId<ptList.size();++ptId)
            addPoint(ptList[ptId][0],ptList[ptId][1]);
    }

    //传入点列表的构造函数,根据点列表直接出范围
    Range2D(const std::vector<TempPtBase>& ptList)
    {
        //使用传入的点列表来初始化这个范围对象
        initRange(ptList);
    }
};

//带有分辨率的范围，仅仅是为了提供一个接口，方便从整数坐标和浮点型坐标的相互转换
class ResolutionRange2D : public Range2D
{
public:
    double pixelLen_=0;//每个像素分辨率的长度
public:
    //使用与父类类似的构造函数
    ResolutionRange2D() : Range2D(){}
    ResolutionRange2D(const std::vector<TempPtBase>& ptList,double pixelLen)
        : Range2D(ptList), pixelLen_(pixelLen){}

    //把浮点型的坐标转换成整型的坐标
    //这里并不负责转换的结果在合理的范围内
    void toInt(const double floatX,const double floatY,
                      int& intX,int& intY) const
    {
        //判断像素的长度是否为0
        if(pixelLen_==0) throw ERROR_ZERO_PIXEL_LEN;
        intX=(floatX-xRange_[0])/pixelLen_;
        intY=(floatY-yRange_[0])/pixelLen_;
    }

    //把点列表批量转换为opencv的点列表
    template<typename PointType>
    void toCvPt(const std::vector<PointType>& srcPtList,
                std::vector<cv::Point>& dstPtList) const
    {
        //给目标点列表开辟空间
        dstPtList.reserve(srcPtList.size());
        //遍历源点列表
        for(unsigned int ptId=0;ptId<srcPtList.size();++ptId)
        {
            //新建一个临时的点
            cv::Point tempPt;
            //当前位置的原始点
            const PointType& thisPt=srcPtList[ptId];
            //获取转换后的点坐标
            toInt(thisPt[0],thisPt[1],tempPt.x,tempPt.y);
            //把转换后的点坐标添加到列表里面
            dstPtList.push_back(tempPt);
        }
    }

    //把整型坐标转换成浮点型坐标
    void toFloat(const int intX,const int intY,double& floatX,double &floatY) const
    {
        //判断像素的长度是否为0
        if(pixelLen_==0) throw ERROR_ZERO_PIXEL_LEN;
        floatX=intX*pixelLen_+xRange_[0];
        floatY=intY*pixelLen_+yRange_[0];
    }

    //获取图片的大小
    cv::Size getImgRange() const
    {
        //新建opencv的范围数据
        cv::Size retSize;
        //计算长宽
        retSize.width=(xRange_[1]-xRange_[0])/pixelLen_;
        retSize.height=(yRange_[1]-yRange_[0])/pixelLen_;
        return retSize;
    }

    //使用范围数据初始化一个opencv的矩阵
    void initCvMat(cv::Mat &dstMat,int cvFlag) const
    {
        //获取图片的范围
        cv::Size imgSize=getImgRange();
        //判断是大小是否为0
        if(imgSize.width*imgSize.height==0) throw ERROR_ZERO_IMG_SIZE;
        dstMat.create(imgSize,cvFlag);
    }
};

//点坐标的中间层基类,方便以后扩展,先有的cornerPt才有的这个，cornerpt已经成熟了，暂时不动cornerpt
class TempPtMid : public TempPtBase
{
public:
    //初始化点坐标
    virtual void initPoint(double xLoc, double yLoc)
    {
        //记录点坐标
        TempPtBase::operator[](0)=xLoc;
        TempPtBase::operator[](1)=yLoc;
    }

    TempPtMid(double xLoc=0, double yLoc=0) : TempPtBase()
    {
        //初始化点坐标
        initPoint(xLoc,yLoc);
    }

    //重载等号操作
    virtual void operator=(const TempPtBase& otherPt)
    {
        //调用父类的等号操作
        TempPtBase::operator=(otherPt);
    }

    //获取x坐标
    double xCoord() const
    {
        return TempPtBase::operator[](0);
    }

    //获取y坐标
    double yCoord() const
    {
        return TempPtBase::operator[](1);
    }
};

//dom问题中，向量的数据类型
class DomDir : public TempPtMid
{
public:
    //构造函数
    DomDir(double xDir=0, double yDir=0) : TempPtMid(xDir,yDir)
    {

    }

    //重载等号
    void operator=(const TempPtBase& otherPt)
    {
        //调用父类的等号操作
        TempPtMid::operator=(otherPt);
    }

    //计算两个向量之间的余弦夹角
    virtual double cosAngle(const TempPtBase& otherDir) const
    {
        return (transpose()*otherDir)[0]/norm()/otherDir.norm();
    }

    //计算两个向量之间的正弦夹角
    virtual double sinAngle(const TempPtBase& otherDir) const
    {
        return std::sqrt(1-std::pow(cosAngle(otherDir),2));
    }

    //两个向量之间的夹角，夹角的范围是[0,pi)
    virtual double dirAngle(const DomDir& otherDir) const
    {
        //计算两个向量的余弦夹角
        double cosValue=cosAngle(otherDir);
        //返回夹角
        return std::acos(cosValue);
    }

    //计算向量的夹角的绝对值,夹角的范围是(-pi,pi]
    virtual double absoluteAngle() const
    {
        //生成主方向
       const DomDir mainDir(1,0);
        //判断y值
        if(yCoord()>0)
        {
            //直接计算
            return dirAngle(mainDir);
        }
        else if(yCoord()<0)
        {
            //计算后取相反数
            return -dirAngle(mainDir);
        }
        else if(xCoord()>0) //判断x的范围
        {
            return 0;
        }
        else if(xCoord()<0) //180度
        {
            return M_PI;
        }
         //异常情况
        std::cout<<"zero vector has no angle\n";
        return -M_PI;
    }

    //将向量标准化,y值不为0则令y值大于0,y值为0则令x大于0,该向量本身来会被改变
    DomDir vectorNorm() const
    {
        //用于返回的数据
        DomDir retDir;
        //判断y方向是否存在
        if(yCoord()>0)
        {
            retDir=*this;
        }
        else if(yCoord()<0)
        {
            //保存相反数的情况
            retDir=-(*this);
        }
        else if(xCoord()>0)
        {
            //记录x坐标
            retDir[0]=TempPtBase::operator[](0);
        }
        else if(xCoord()<0)
        {
            //记录x坐标的相反数
            retDir[0]=-TempPtBase::operator[](0);
        }
        //返回标准化后的数据
        return retDir;
    }
};

//表示图片范围的数据类型
class ImgRange : public ImgRangeBase
{
public:
    //分别获取左右上下
    int left() const{return operator[](0);}
    int right() const{return operator[](1);}
    int up() const{return operator[](2);}//对于图片范围来说，上小下大
    int down() const{return operator[](3);}

    //获取引用的形式
    int& left() {return operator[](0);}
    int& right() {return operator[](1);}
    int& up() {return operator[](2);}//对于图片范围来说，上小下大
    int& down() {return operator[](3);}

    //范围窗口的迭代信息,用于辅助观察
    std::string iterInfo_="";

    //计算宽度
    unsigned int getWidth() const
    {
        return right()-left();
    }

    //计算高度
    unsigned int getHeight() const
    {
        return down()-up();
    }

    //计算范围里面的面积
    unsigned int rangeArea() const
    {
        return getWidth()*getHeight();
    }

    //将范围整体平移
    void translate(int transX,int transY)
    {
        //分别把平移叠加上去
        if(transX)
        {
            left()+=transX;
            right()+=transX;
        }
        if(transY)
        {
            up()+=transY;
            down()+=transY;
        }
    }

    //将范围的左上角移动到某个位置
    void moveTo(int dstX,int dstY)
    {
        //记录范围的长度和高度
        int currWidth=getWidth();
        int currHeight=getHeight();
        //把坐标移动过去
        left()=dstX;
        up()=dstY;
        //根据当时的长度更新右下的信息
        right()=dstX+currWidth;
        down()=dstY+currHeight;
    }

    //获取范围内的所有点坐标
    void getPtsInRange(CvPointVector& dstPtList) const
    {
        //根据范围的面积，提前给点的列表开辟空间
        dstPtList.reserve(rangeArea());
        //遍历范围内的点
        for(int colCount=left();colCount<right();++colCount)
        {
            for(int rowCount=up();rowCount<down();++rowCount)
                dstPtList.push_back(CvPoint(colCount,rowCount));
        }
    }

    //获取采样点
    void getSamplePt(CvPointVector &dstPtList) const
    {
        //提前开辟空间
        dstPtList.reserve(9);
        //计算中间值
        int midX=left()+getWidth()/2;
        int midY=up()+getHeight()/2;
        dstPtList.push_back(CvPoint(left(),up()));//左上
        dstPtList.push_back(CvPoint(midX,up()));//中上
        dstPtList.push_back(CvPoint(right(),up()));//右上
        dstPtList.push_back(CvPoint(left(),midY));//左中
        dstPtList.push_back(CvPoint(midX,midY));//中中
        dstPtList.push_back(CvPoint(right(),midY));//右中
        dstPtList.push_back(CvPoint(left(),down()));//左下
        dstPtList.push_back(CvPoint(midX,down()));//中下
        dstPtList.push_back(CvPoint(right(),down()));//右下
    }

    //将像素的移动方向规约成适合移动的像素方向
    static void adjustImgDir(double &xDir,double &yDir)
    {
        //如果两个都是0,直接返回
        if(xDir==0 && yDir==0) return;
        //计算两个方向的绝对值
        double xDirAbs=std::abs(xDir);
        double yDirAbs=std::abs(yDir);
        //判断如果x等于0
        if(xDir==0)
        {
            //更新方向
            yDir=yDir/yDirAbs;
            return;
        }
        //判断如果y等于0
        if(yDir==0)
        {
            //更新方向
            xDir=xDir/xDirAbs;
            return;
        }
        //初始化两个要除以的数字
        double divValue=xDirAbs;
        //判断是否y更大
        if(yDirAbs>xDirAbs) divValue=yDirAbs;
        //两个数字做除法
        xDir/=divValue;
        yDir/=divValue;
    }

    //获取两个点之间的线段上的点
    //这个代码以前写过，但是接口不太匹配
    //这个代码的接口更通用一些
    static void getInternalPt(const CvPoint& pt1,
                              const CvPoint& pt2,
                              CvPointVector& dstPtList)
    {
        //计算x的增量
        double xStep=pt2.x-pt1.x;
        //计算y的增量
        double yStep=pt2.y-pt1.y;
        //两个点之间的距离
        int ptDistance=std::max(std::abs(pt1.x-pt2.x),std::abs(pt1.y-pt2.y));
        //根据距离判断一下是否有必要加入
        if(ptDistance<=1)
        {
            //只添加第1个点就够了
            dstPtList.push_back(pt1);
            return;
        }
        //提前开辟空间
        dstPtList.reserve(ptDistance);
        //调整方向
        adjustImgDir(xStep,yStep);
        //初始化xy目前的位置
        CvPoint currPt=pt1;
        //移动的步长
        int forwardStep=0;
        //向后迭代寻找
        while(true)
        {
            //记录当前位置的点
            dstPtList.push_back(currPt);
            //移动的步长增加
            ++forwardStep;
            //移动固定个数的步长
            if(forwardStep>=ptDistance) break;
            //更新下一个点的位置
            currPt.x=pt1.x+xStep*forwardStep;
            currPt.y=pt1.y+yStep*forwardStep;
        }
    }

    //判断当前位置是否位于范围的边角中位置
    bool judgeMidSide(int srcX,int srcY) const
    {
        //x的三个值
        std::map<int,bool> xList;
        xList[left()]=true;
        xList[right()]=true;
        xList[left()+getWidth()/2]=true;
        std::map<int,bool> yList;
        yList[up()]=true;
        yList[down()]=true;
        yList[up()+getHeight()/2]=true;
        //判断是否存在
        return (xList.count(srcX)>0) && (yList.count(srcY)>0);
    }

    //初始化范围的数据
    void initValue(int leftValue,int rightValue,int upValue,int downValue)
    {
        left()=leftValue;
        right()=rightValue;
        up()=upValue;
        down()=downValue;
    }
};
//表示图片范围的列表,做不成功的面片的buffer的时候会用到这个东西
typedef std::vector<ImgRange> ImgRangeList;

//每个三角面片的角点，sfm封装的这个vec2感觉用着不舒服
class   CornerPt : public TempPtBase
{
public:
    //在点云中该点的标号，如果是图片点，指的是它在点云中的对应点
    //如果是点云点，那指的就是点云点本身
    IndexT idx_=0;

    //重载一下等号，需要用
    void operator=(const CornerPt& fromPt)
    {
        //记录点的标号
        idx_=fromPt.idx_;
        //调用父类的等号函数,并把传入的东西强制转换为父类
        TempPtBase::operator=(static_cast<TempPtBase>(fromPt));
    }

    //传入普通点的等号操作
    void operator=(const TempPtBase& fromPt)
    {
        //调用父类的等号函数
        TempPtBase::operator=(fromPt);
    }

    //初始化当前的点,为了复用，单独拿出来弄成个函数
    void initPoint(IndexT cloudIdx,double xLoc,double yLoc)
    {
        //记录点云点的坐标
        idx_=cloudIdx;
        //写入坐标
        TempPtBase::operator[](0)=xLoc;
        TempPtBase::operator[](1)=yLoc;
    }

    //传入点云的构造方式
    void initPoint(const Landmarks& cloudPoints,IndexT cloudIdx)
    {
        //获取当前位置的点坐标
        const double* const ptCArray=cloudPoints.at(cloudIdx).X.data();
        //调用初始化函数
        initPoint(cloudIdx,ptCArray[0],ptCArray[1]);
    }

    //根据特定图片进行构造
    void initPoint(const Landmarks& cloudPoints, IndexT cloudIdx, IndexT imgIdx)
    {
        //获取当前位置的图片坐标
        const double* const ptCArray=cloudPoints.at(cloudIdx).obs.at(imgIdx).x.data();
        //调用初始化的构造函数
        initPoint(cloudIdx,ptCArray[0],ptCArray[1]);
    }

    //构造函数，构造的时候可以传入点坐标，点标号必须传
    CornerPt(IndexT cloudIdx=0,double xLoc=0,double yLoc=0) : TempPtBase()
    {
        //调用初始化函数
        initPoint(cloudIdx,xLoc,yLoc);
    }

    //直接传入一个点云列表，然后从点云列表里面根据点标号获取坐标
    CornerPt(const Landmarks& cloudPoints,IndexT cloudIdx)
    {
        initPoint(cloudPoints,cloudIdx);
    }

    //根据指定点在指定图片的投影点进行构造
    CornerPt(const Landmarks& cloudPoints, IndexT cloudIdx, IndexT imgIdx)
    {
        initPoint(cloudPoints,cloudIdx,imgIdx);
    }

    //获取点的齐次形式,这个功能本应该被放在基类里面,由于种种历史原因。。。
    void gethomoModel(HomoPt& dstPt)
    {
        dstPt<<operator[](0),operator[](1),1.f;
    }
};

//任意维度的向量接口,勉强实现一些简单的数学运算
//子类最好要能够支持随机访问，不能随机访问的话，时间消耗可能会比较大
template<typename T>
class VectorInterface
{
public:
    //用来实现对索引的访问
    virtual T& operator[](unsigned int idx) =0;
    virtual const T& operator[](unsigned int idx) const=0;

    //功能类似于一个索引,c++的enum感觉不好用
    struct IdList
    {
        const static int VecBase=0;//基类
        const static int WeightColor=1;//带权值的颜色
    };

    //为了方便判断指针或引用对象属于哪个基类
    virtual int selfId() const
    {
        return IdList::VecBase;
    }

    //向量的大小
    virtual std::size_t size() const=0;

    //转换成eigen向量
    virtual Eigen::Matrix<T,Eigen::Dynamic,1> toEigen() const
    {
        //新建用于返回的向量
        Eigen::Matrix<T,Eigen::Dynamic,1> retVec;
        //指定大小
        retVec.resize(size(),1);
        //遍历向量，赋值
        for(int dCount=0;dCount<size();++dCount)
        {
            retVec[dCount]=operator[](dCount);
        }
        //返回结果
        return retVec;
    }

    //转换成特定大小的eigen向量，只能与上层模板相配合的时候才能使用
    template<unsigned int chn>
    Eigen::Matrix<double,chn,1> getEigen() const
    {
        //判断指定的通道数和自己原来的通道数是否一致
        if(size()!=chn) throw ERROR_EIGEN_VEC_DIFF_DIM;
        //新建一个eigen数据用于返回
        Eigen::Matrix<double,chn,1> retEigen;
        //遍历向量，依次给每个位置赋值
        for(unsigned int dimCount=0;dimCount<chn;++dimCount)
        {
            retEigen[dimCount]=operator[](dimCount);
        }
        //返回最后的结果
        return retEigen;
    }

    //判断是否所有元素都是某个定值
    bool judgeAllValue(T value=0) const
    {
        //遍历整个向量
        for(unsigned int vecCount=0;vecCount<size();++vecCount)
            if(operator[](vecCount)!=value) return false;
        return true;
    }

    //计算向量的总和
    virtual T sum() const
    {
        //初始化加和
        T sumValue=0;
        //遍历数据
        for(int dimCount=0;dimCount<size();++dimCount)
            sumValue=sumValue+operator[](dimCount);
        //返回结果
        return sumValue;
    }

    //计算指数次方和
    virtual double indexSum(const double index=2) const
    {
        double sumValue=0;
        //遍历数据
        for(unsigned int dimCount=0;dimCount<size();++dimCount)
            sumValue=sumValue+std::pow(operator[](dimCount),index);
        //返回数据
        return sumValue;
    }

    //向量的二范数
    virtual double norm() const
    {
        //平方和再开根
        return std::sqrt(indexSum(2));
    }

    //计算两个向量的欧氏距离
    virtual double getEuclidDis(const VectorInterface &otherVec) const
    {
        //判断向量大小是否一致
        if(size()!=otherVec.size())
        {
            std::cout<<"getEuclidDis require equal dimesion vector\n";
            return -1;
        }
        //初始化差值平方和
        double disSum=0;
        //遍历每个元素
        for(IndexT dimCount=0;dimCount<size();++dimCount)
        {
            //计算差值
            double tempDis=this->operator[](dimCount)-otherVec[dimCount];
            //距离累计
            disSum+=tempDis*tempDis;
        }
        //开方返回距离
        return std::sqrt(disSum);
    }

    //计算两个向量的数量积
    virtual double dotProduct(const VectorInterface &otherVec) const
    {
        //判断向量大小是否一致
        if(size()!=otherVec.size())
        {
            std::cout<<"dotProduct require equal dimesion vector\n";
            return -1;
        }
        //初始化总和
        double sumValue=0;
        //检查是否有非0数字
        bool zeroFlag=true;
        //遍历向量
        for(unsigned int dimCount=0;dimCount<size();++dimCount)
        {
            if(otherVec[dimCount]!=0 || operator[](dimCount)!=0)
            {
                zeroFlag=false;
            }
            //记录相乘数据
            sumValue+=operator[](dimCount)*otherVec[dimCount];
        }
        //如果是全0的，返回1
        if(zeroFlag) return 1;
        //返回相加结果
        return sumValue;
    }

    //计算向量的余弦距离[-1,1],传入另一个向量,两个向量越接近，返回值越接近1
    virtual double cosDis(const VectorInterface& otherVec) const
    {
        //判断向量是否大小一致
        if(size()!=otherVec.size())
        {
            std::cout<<"cos distance cannot be computed by different dimesion vector\n";
            return -1;
        }
        //计算两个向量的数量积然后再除以它们各自的二范数
        return dotProduct(otherVec)/norm()/otherVec.norm();
    }

    //向量整体加上一个数字
    virtual void operator+=(const T& value)
    {
        //遍历向量
        for(unsigned int dimCount=0;dimCount<size();++dimCount)
        {
            //加法
            operator[](dimCount)+=value;
        }
    }

    //把向量整体乘上一个数字
    virtual void operator*=(const T& value)
    {
        //遍历向量
        for(unsigned int dimCount=0;dimCount<size();++dimCount)
        {
            //加法
            operator[](dimCount)*=value;
        }
    }

    //向量整体减去某个数字
    virtual void operator-=(const T &value)
    {
        //调用相反数的加号
        operator+=(-value);
    }

    //向量整体减去某个向量
    virtual void operator-=(VectorInterface<T>& otherVec)
    {
        //判断维度是否相等
        if(otherVec.size()!=size()) throw ERROR_VEC_SIZE_DIFF;
        //遍历向量里面的每个值
        for(unsigned int dimCount=0;dimCount<size();++dimCount)
        {
            //对应位置相减
            operator[](dimCount)-=otherVec[dimCount];
        }
    }

    //把向量的二范数弄成一个指定的数值，但保持比例不变
    virtual void setNormAs(double dstNorm)
    {
        //计算向量的范数
        double currNorm=norm();
        //如果目前的范数是0,那它不可能被修改，那就直接返回吧
        if(currNorm==0) return;
        //如果两个范数相同，就不用做了
        if(currNorm==dstNorm) return;
        //把向量整体乘它们的相对值
        operator*=(dstNorm/currNorm);
    }

    //把向量复制给另外一个向量
    //需要保证它们有相同的大小，这个函数不负责保证它能复制成功
    //因为子类可能是一个固定大小的向量
    //如果复制失败就返回false了
    bool copyTo(VectorInterface<T> &dstVector) const
    {
        //判断是否具有相同的大小
        if(dstVector.size()!=size())
        {
            std::cout<<"copy method of VectorInterface require a same size vector\n";
            return false;
        }
        //遍历向量的大小，把每个位置的数值都赋值给它
        for(unsigned int dimCount=0;dimCount<size();++dimCount)
            dstVector[dimCount]=operator[](dimCount);
        return true;
    }

    //复制一份向量，然后把向量的二范数弄成指定的数值
    //把最后的向量存储到dstVector里面
    //如果最后没有弄成就返回false了
    virtual bool getNormAs(double dstNorm, VectorInterface<T> &dstVector) const
    {
        //把向量复制给这个传入的向量
        bool copyStage=copyTo(dstVector);
        if(!copyStage) return false;
        //把复制后的向量转换成对应的范数
        dstVector.setNormAs(dstNorm);
        return true;
    }

    //把另外一个向量线性相加，加到当前向量上,第二个参数是相加的系数
    virtual bool linearAdd(const VectorInterface<T>& otherVec,double coef=1)
    {
        //确定两个向量长度是否相同
        if(otherVec.size()!=size())
        {
            std::cout<<"only vector with same size can be added";
            return false;
        }
        //遍历向量，按照比例相加
        for(unsigned int dimCount=0;dimCount<size();++dimCount)
        {
            //临时的数据
            T& thisData=operator[](dimCount);
            //在临时数据上根据系数做叠加
            thisData+=coef*otherVec[dimCount];
        }
        //正常情况下返回true
        return true;
    }

    //两个向量的直接相加
    virtual void operator+=(const VectorInterface<T>& otherVec)
    {
        //直接相加即可
        linearAdd(otherVec);
    }

    //判断一组向量的长度是否一致
    //如果一致的话返回大小,否则返回0
    static unsigned int judgeSizeEqual(const std::vector<VectorInterface<T>*>& vecList)
    {
        //判断是否为空向量
        if(vecList.size()==0) return 0;
        //第1组数据的长度
        std::size_t allSize=vecList[0]->size();
        //遍历所有的数据
        for(unsigned int vecCount=1;vecCount<vecList.size();++vecCount)
        {
            //存在不相同的数据的时候返回false
            if(vecList[vecCount]->size()!=allSize) return 0;
        }
        return allSize;
    }

    //判断一组向量的长度是否一致
    //如果一致的话返回大小,否则返回0
    static unsigned int judgeSizeEqual(const std::vector<const VectorInterface<T>*>& vecList)
    {
        //判断是否为空向量
        if(vecList.size()==0) return 0;
        //第1组数据的长度
        std::size_t allSize=vecList[0]->size();
        //遍历所有的数据
        for(unsigned int vecCount=1;vecCount<vecList.size();++vecCount)
        {
            //存在不相同的数据的时候返回false
            if(vecList[vecCount]->size()!=allSize) return 0;
        }
        return allSize;
    }

    //计算一组向量的互相关性
    //直接仿照两个向量做的互相关性的操作最后发现效果是不正常的
    static double crossCorrelation(const std::vector<VectorInterface<T>*>& vecList)
    {
        //判断是否为空向量
        if(vecList.size()==0) return 0;
        //向量的统一大小
        unsigned int allSize=judgeSizeEqual(vecList);
        //判断是否所有的数据维度都相同
        if(!allSize)
        {
            std::cout<<"only equal vector can be computed crossCorrelation\n";
            return 0;
        }
        //初始化最后的返回结果
        double retValue=0;
        //遍历每个维度
        for(unsigned int dimCount=0;dimCount<allSize;++dimCount)
        {
            //当前维度的乘积
            double dimProduct=1;
            //遍历每个向量在这个维度上的数值
            for(unsigned int vecCount=0;vecCount<vecList.size();++vecCount)
            {
                dimProduct*=vecList[vecCount]->operator[](dimCount);
            }
            //把相乘结果加在最后的数据上
            retValue+=dimProduct;
        }
        //依次除以每个向量的模
        for(unsigned int vecCount=0;vecCount<vecList.size();++vecCount)
            retValue/=vecList[vecCount]->norm();
        //返回最后的结果
        return retValue;
    }

    //把所有的向量正则化,把所有的向量的norm都弄成1
    static void uniformVectors(std::vector<VectorInterface<T>*> &vecList)
    {
        //遍历每个向量
        for(unsigned int vecCount=0;vecCount<vecList.size();++vecCount)
            vecList[vecCount]->setNormAs(1);
    }

    //计算传入的所有向量的平均模长
    static double getAvgNorm(const std::vector<VectorInterface<T>*> &vecList)
    {
        //初始化模长
        double normSum=0;
        //计算传入的所有向量的平均模长
        for(unsigned int vecCount=0;vecCount<vecList.size();++vecCount)
        {
            normSum+=vecList[vecCount]->norm();
        }
        //返回平均的模长
        return normSum/vecList.size();
    }

    //释放每个vector指针
    static void freeEachPtr(std::vector<VectorInterface<T>*>& vecList)
    {
        //遍历每个向量
        for(unsigned int vecCount=0;vecCount<vecList.size();++vecCount)
            delete vecList[vecCount];
    }

    //计算所有向量的角度标准差，每个向量计算与平均向量的角度，最后再取平均，从而评判一下聚合程度
    static double computeCosVairance(const std::vector<VectorInterface<T>*> &vecList,//所有向量的指针
                                const VectorInterface<T> &avgVec,//平均向量
                                   bool uniformFlag=false,//如果这里是true，就默认所有向量都是做过归一化的，这样计算角度距离的时候就直接点积了
                                     unsigned int* nearestVec=nullptr) //当这里传入一个具体的指针的时候，把最小的数据记录在里面
    {
        //初始化每个向量和平均向量的cos距离之和
        double cosDisSum=0;
        //初始化最近的距离
        double nearestDis=0;
        //遍历每个向量，计算与平均向量的距离
        for(unsigned int vecCount=0;vecCount<vecList.size();vecCount++)
        {
            //临时的余弦距离
            double tempCosDis=0;
            //判断是否做过归一化
            if(uniformFlag)
                tempCosDis=avgVec.dotProduct(*vecList[vecCount]);
            else //没做过归一化的时候，就老老实实使用每个位置的余弦距离依次算吧
                tempCosDis=avgVec.cosDis(*vecList[vecCount]);
            //判断是否找到了更接近的距离
            if(tempCosDis>nearestDis)
            {
                //记录最近的距离
                nearestDis=tempCosDis;
                //判断是否需要记录最近的标号
                if(nearestVec) *nearestVec=vecCount;
            }
            //在总数里面记录距离之和
            cosDisSum+=tempCosDis;
        }
        //返回平均值
        return cosDisSum/vecList.size();
    }

    //计算所有向量的平均向量
    //需要从外部保证所有的向量的维度都是相同的
    static void getAverageValue(const std::vector<VectorInterface<T>*>& vecList,
                                VectorInterface<T>& dstVec)
    {
        //判断是否所有的向量的维度都相同
        unsigned int allSize=judgeSizeEqual(vecList);
        //判断是否具有相同的维度
        if(dstVec.size()!=allSize)
        {
            std::cout<<"require equal vector size in getAverageValue\n";
            return;
        }
        //遍历向量的每个维度
        for(unsigned int dimCount=0;dimCount<allSize;++dimCount)
        {
            //初始化目标维度上的当前值
            dstVec[dimCount]=0;
            //遍历每个向量，把数值加上去
            for(unsigned int vecCount=0;vecCount<vecList.size();++vecCount)
                dstVec[dimCount]+=vecList[vecCount]->operator[](dimCount);
            //取平均值
            dstVec[dimCount]/=vecList.size();
        }
    }

    //计算所有向量的平均向量
    //需要从外部保证所有的向量的维度都是相同的
    static void getAverageValue(std::vector<const VectorInterface<T>*>& vecList,
                                VectorInterface<T>& dstVec)
    {
        //判断是否所有的向量的维度都相同
        unsigned int allSize=judgeSizeEqual(vecList);
        //判断是否具有相同的维度
        if(dstVec.size()!=allSize)
        {
            std::cout<<"require equal vector size in getAverageValue\n";
            return;
        }
        //遍历向量的每个维度
        for(unsigned int dimCount=0;dimCount<allSize;++dimCount)
        {
            //初始化目标维度上的当前值
            dstVec[dimCount]=0;
            //遍历每个向量，把数值加上去
            for(unsigned int vecCount=0;vecCount<vecList.size();++vecCount)
                dstVec[dimCount]+=vecList[vecCount]->operator[](dimCount);
            //取平均值
            dstVec[dimCount]/=vecList.size();
        }
    }

    //把vectorInterface的某个子类封装成vectorInterface的指针
    template<typename SonClass>
    static void makeInterfacePtr(std::vector<SonClass>& sonList,//子类的列表
            std::vector<VectorInterface<T>*>& dstInterfaceList) //目标的接口列表
    {
        //给接口的列表提前开辟空间
        dstInterfaceList.resize(sonList.size());
        //遍历子类的列表
        for(unsigned int idCount=0;idCount<sonList.size();++idCount)
            dstInterfaceList[idCount]=&(sonList[idCount]);
    }

    //把vectorInterface的某个子类封闭成const指针的列表
    template<typename SonClass>
    static void makeInterfaceConstPtr(const std::vector<SonClass>& sonList,
                                      std::vector<const VectorInterface<T>*>& dstInterfaceList)
    {
        //给接口的列表提前开辟空间
        dstInterfaceList.resize(sonList.size());
        //遍历子类的列表
        for(unsigned int idCount=0;idCount<sonList.size();++idCount)
            dstInterfaceList[idCount]=&(sonList[idCount]);
    }

    //把若干向量转换成一个eigen矩阵,n行,chn列，每个向量的具体维度是chn
    template<unsigned int chn>
    static void makeEigenMat(const std::vector<VectorInterface<T>*>& srcVecList,
                             Eigen::Matrix<double,Eigen::Dynamic,chn>& dstMat)
    {
        //提前给最终的矩阵开辟空间
        dstMat.resize(srcVecList.size(),chn);
        //遍历每个向量
        //这里是可以使用多线程的，但这是一个很细枝末节的地方
        for(unsigned int vecCount=0;vecCount<srcVecList.size();++vecCount)
        {
            //当前位置向量的引用
            const VectorInterface<T>& thisVec=*(srcVecList[vecCount]);
            //遍历每个维度，依次给矩阵的当前行赋值
            for(unsigned int dimCount=0;dimCount<chn;++dimCount)
            {
                dstMat(vecCount,dimCount)=thisVec[dimCount];
            }
        }
    }

    //根据一个n行,chn列的矩阵，计算这个矩阵的最小特征值对应的特征向量
    //这里说的最小特征值是绝对值最小的特征值
    template<unsigned int chn>
    static void getMinEigen(const Eigen::Matrix<double,Eigen::Dynamic,chn>& srcMat,//用于计算特征值和特征向量的矩阵
                            Eigen::Matrix<double,chn,1>& dstVec)//最小特征值对应的特征向量结果
    {
        //转置相乘，传统套路
        Eigen::Matrix<double,chn,chn> mainMat=srcMat.transpose()*srcMat;
        //计算特征值和特征向量
        Eigen::EigenSolver<Eigen::Matrix<double,chn,chn>> eigSolver(mainMat);
        //获取所有特征值的绝对值
        Eigen::Matrix<double,chn,1> eigList=eigSolver.eigenvalues().array().abs();
        //获取最小特征值的位置
        Eigen::MatrixXd::Index minIdx;
        eigList.minCoeff(&minIdx);
        //根据最小特征值的位置取出最小特征值
        dstVec=eigSolver.eigenvectors().col(minIdx).transpose().real();
    }

    //传入一组向量，将这些向量组合成一组eigen数据，然后计算这组eigen数据的最小特征值对应的特征向量
    template<unsigned int chn>
    static void getMinEigen(const std::vector<VectorInterface<T>*> &srcVec,
                            VectorInterface<double>& dstEigen)//最小特征值对应的特征向量会被存储在这里
    {
        //判断目标向量的维度和指定的维度是否一致
        if(dstEigen.size()!=chn) throw ERROR_EIGEN_VEC_DIFF_DIM;
        //用向量的列表生成一组矩阵
        Eigen::Matrix<double,Eigen::Dynamic,chn> cmbVec;
        makeEigenMat<chn>(srcVec,cmbVec);
        //根据生成的矩阵获取最小特征值
        Eigen::Matrix<double,chn,1> minEigen;
        getMinEigen<chn>(cmbVec,minEigen);
        //把eigen向量转换到目标向量中
        for(unsigned int dimCount=0;dimCount<chn;++dimCount)
            dstEigen[dimCount]=minEigen[dimCount];
    }
};//VectorInterface的结尾
//定义一下double的vector数据类型，方便使用
typedef VectorInterface<double> VecInterDouble;

//三维点的数据类型
class Point3D : public Point3DBase, public VectorInterface<double>
{
public:
    //用来实现对索引的访问
    virtual double& operator[](unsigned int idx) override
    {
        return Point3DBase::operator[](idx);
    }

    //实现对索引的访问
    virtual const double& operator[](unsigned int idx) const override
    {
        return Point3DBase::operator[](idx);
    }

    //空的构造函数，仅仅是一个可选的项目
    Point3D(){}

    //用一个eigen的数据类型来初始化这个三维点数据
    //直接调用基类在这方面的构造函数
    Point3D(const Point3DBase& otherPt) : Point3DBase(otherPt)
    {

    }

    //传入x,y,z实现的初始化方式
    Point3D(double xValue,double yValue,double zValue)
    {
        Point3DBase::operator[](0)=xValue;
        Point3DBase::operator[](1)=yValue;
        Point3DBase::operator[](2)=zValue;
    }

    //获取z坐标
    double getZ() const{return Point3DBase::operator[](2);}

    //实现父类的size函数
    virtual std::size_t size() const override {return 3;}
};
//三维数据点的列表
typedef std::vector<Point3D> Point3DList;

//三维向量的数据类型
//虽说这个东西很大程度上和以前的DomDir重复了，但这里的很多东西都可以复用vectorInterface
class Dir3D : public Point3D
{
public:
    //通过传入两个点来构造一个三维向量
    //从第1个点到第2个点
    Dir3D(const Point3D& fromPt,const Point3D& toPt)
    {
        //对两个点做减法
        this->Point3DBase::operator=(toPt-fromPt);
    }

    //允许使用无参数的构造函数
    Dir3D(){}
};
//定义点到向量的互相减数据类型
class PtDirEachMinus : public EachOperateInterface<Point3D,Dir3D>
{
public:
    //实现接口里面要求的减法操作
    Dir3D compute(const Point3D &src1, const Point3D &src2) override
    {
        //返回两个点的相减结果
        //向量数据类型的构造函数可以自动满足相减的操作
        return Dir3D(src1,src2);
    }
};

//自己用来搞事情的颜色数据接口，因为有的数据已经实现了颜色的功能，所以只好做一个上层接口了
//使用的T最好是一种常规的数据类型，不然可能会出问题
//仅限三通道的颜色，其它颜色模型暂时用不到就不做了，用到了再说
template<typename T>
class ColorInterface : public VectorInterface<T>
{
protected:
    //是否已经把每个颜色通道弄成了均值为0的形式
    bool avgZero_=false;
public:
    //获取每种颜色数据,并且是可以操作的
    virtual T& getR() =0;
    virtual T& getG() =0;
    virtual T& getB() =0;
    virtual const T& getR() const=0;
    virtual const T& getG() const=0;
    virtual const T& getB() const=0;

    //向量的大小
    virtual std::size_t size() const override
    {
        return 3;
    }

    //实现父类的纯虚函数
    virtual T& operator[](unsigned int idx) override
    {
        //按照bgr的顺序来返回
        if(idx==0) return getB();
        else if(idx==1) return getG();
        else if(idx==2) return getR();

        std::cout<<"color have no more than 3 value\n";
        return getR();
    }

    //实现父类的纯虚函数
    virtual const T& operator[](unsigned int idx) const override
    {
        //按照bgr的顺序来返回
        if(idx==0) return getB();
        else if(idx==1) return getG();
        else if(idx==2) return getR();

        std::cout<<"color have no more than 3 value\n";
        return getR();
    }

    //获取颜色的数据，按照opencv的接口来返回
    //按照BGR的格式来存储
    virtual ColorT getCvColor() const
    {
        return ColorT(getB(),getG(),getR());
    }

    //计算两个颜色的色度距离，为了消除亮度的影响，这里使用的是余弦距离
    virtual double colorDis(const ColorInterface<T>& otherColor)
    {
        //直接计算两个颜色的余弦距离，本质上是一样的
        return VectorInterface<T>::cosDis(otherColor);
    }

    //每个颜色通道的最大值
    virtual T maxValue() const=0;

    //复制内容
    bool copyTo(ColorInterface<T> &otherColor) const
    {
        //调用父类的复制
        bool copyStage=VectorInterface<T>::copyTo(otherColor);
        if(!copyStage) return false;
        //复制是否归0的信息
        otherColor.avgZero_=avgZero_;
        //默认返回true
        return true;
    }

    //把颜色的均值设置为0
    virtual void setAvgZero()
    {
        //如果已经设置过了，就不用设置了
        if(avgZero_) return;
        //每个通道都减去通道最大值的一半
        VectorInterface<T>::operator-=(maxValue()/2);
        //记录颜色模型
        avgZero_=true;
    }

    //把颜色的均值恢复为正常
    virtual void restoreAvgZero()
    {
        //如果本来就是正常的，就不用恢复了
        if(!avgZero_) return;
        //每个通道都增加最大值的一半
        VectorInterface<T>::operator+=(maxValue()/2);
        //标记颜色正常
        avgZero_=false;
    }

    //获取当前颜色的均值为零的版本
    virtual void getAvgZero(ColorInterface<T> &otherColor) const
    {
        //把当前的颜色复制给它,这里的copyTo不会复制子类的属性
        this->copyTo(otherColor);
        //把它的颜色信息设置成均值为0
        otherColor.setAvgZero();
    }

    //获取颜色的亮度，并且是那种模板的类型，传入的是一个指针，这应该属于通用性最强的接口了
    //颜色的默认排列顺序是BGR
    static double colorBrightness(const T* colorData,uint colorSize=3)
    {
        //初始化颜色通道的最大值
        double maxChn=0;
        //遍历每个颜色通道
        for(uint chnCount=0;chnCount<colorSize;++chnCount)
        {
            //判断是否为更大的颜色
            if(colorData[chnCount]>maxChn)
                maxChn=colorData[chnCount];
        }
        //返回统计得到的最大值
        return maxChn;
    }

    //把每个通道的颜色包装成指针
    double getBright() const
    {
        T colorData[3]={getB(),getG(),getR()};
        //传入到静态函数中用于计算亮度
        return colorBrightness(colorData);
    }
};

//为了方便计算颜色之间的距离，定义的一个类，同时继承opencv的颜色数据和自定义的虚接口
//加上一个T是为了区别后面的UserColor
template <typename T>
class UserColorT : public cv::Vec<T,3>, public ColorInterface<T>
{
public:
    //实现获取每种颜色数据
    virtual T & getB() override{return cv::Vec<T,3>::operator[](0);}
    virtual T & getG() override{return cv::Vec<T,3>::operator[](1);}
    virtual T & getR() override{return cv::Vec<T,3>::operator[](2);}
    virtual const T & getB() const override{return cv::Vec<T,3>::operator[](0);}
    virtual const T & getG() const override{return cv::Vec<T,3>::operator[](1);}
    virtual const T & getR() const override{return cv::Vec<T,3>::operator[](2);}

    //空的构造函数
    UserColorT(){}

    //传入一个vec3b的构造函数
    UserColorT(const cv::Vec3b& otherColor)
        : cv::Vec<T,3>(otherColor[0],otherColor[1],otherColor[2]){}

    //直接传入三个值从而实现的构造函数
    UserColorT(T data1,T data2,T data3) : cv::Vec<T,3>(data1,data2,data3){}

    //每个颜色通道的最大值
    virtual T maxValue() const override{return 255;}
};

//带权值的颜色
template<typename T>
class WeightColor : public UserColorT<T>
{
public:
    double weight_=0;//权值
public:
    //空的构造函数
    WeightColor(double weight=0) : UserColorT<T>()
    {
        //记录权值
        weight_=weight;
    }

    //传入颜色的构造函数
    WeightColor(double weight,const cv::Vec3b &otherColor) : UserColorT<T>(otherColor)
    {
        //记录权值
        weight_=weight;
    }

    //返回自己的身份，用于标识子类类型
    virtual int selfId() const override
    {
        return VectorInterface<T>::IdList::WeightColor;
    }

    //重载vector那一层的copyTo
    virtual bool copyTo(WeightColor<T>& otherColor) const
    {
        //调用vector那一层的复制
        bool copyStage=ColorInterface<T>::copyTo(otherColor);
        //如果没有复制成功就算了
        if(copyStage==false) return false;
        //复制权值
        otherColor.weight_=this->weight_;
        return true;
    }

    //带权值的颜色融合
    //仅仅融合色度而不融合亮度
    //有时候在外部已经比较过了两个norm，就不用再比较了
    virtual void fuseColor(const WeightColor<T>& otherColor, bool cmpNorm=true)
    {
        //判断是否需要比较两个范数
        if(cmpNorm)
        {
            //分别计算两个范数
            double norm1=this->norm();
            double norm2=otherColor.norm();
            //判断是否第2个norm较大
            if(norm2>norm1)
            {
                //把范数重新处理一下就行
                this->operator*=(norm2/norm1);
                fuseColor(otherColor,false);
                return;
            }
            else
            {
                //新建一个临时的权值颜色
                WeightColor tempColor;
                //把另外一个颜色复制成当前的颜色，并更改范数
                otherColor.copyTo(tempColor);
                tempColor*=(norm1/norm2);
                //用临时生成的另外一个颜色来做融合
                fuseColor(tempColor,false);
                return;
            }
        }
        //把当前的颜色整合成权值的形式
        this->operator*=(weight_);
        //权值相加
        weight_+=otherColor.weight_;
        //两个颜色线性相加
        this->linearAdd(otherColor,otherColor.weight_);
        //恢复到平均颜色的形式
        this->operator*=(1/weight_);
    }

    //比较两个颜色是否相同，如果两个颜色相同，就保存成范数比较大的那个，并且保存权值较大的那个
    //如果两个颜色不同，按照权值大的保存
    //如果当前属性发生了改变，返回false，否则返回true
    virtual bool cmpColor(WeightColor<T> otherColor)
    {
        //计算颜色的距离
        double colDistance=this->colorDis(otherColor);
        //判断两个颜色是否相同
        if(colDistance<COL_THRE)//颜色不同的情况
        {
            //比较权值
            if(otherColor.weight_>this->weight_)
            {
                //把它的属性复制到这里
                otherColor.copyTo(*this);
                return true;
            }
            //如果没有改变就返回false
            return false;
        }
        else//颜色相同的情况
        {
            //两个颜色的亮度
            double norm1=this->norm();
            double norm2=otherColor.norm();
            bool changeFlag=false;
            //比较亮度
            if(norm2>norm1)
            {
                //直接复制颜色
                otherColor.ColorInterface<T>::copyTo(*this);
                changeFlag=true;
            }
            //比较权值
            if(otherColor.weight_>this->weight_)
            {
                this->weight_=otherColor.weight_;
                changeFlag=true;
            }
            return changeFlag;
        }
        //正常情况下返回false
        return false;
    }
};

//目前主要使用的userColor的数据类型
typedef double UserColorType;
typedef UserColorT<UserColorType> UserColor;
//用户常用颜色的列表
typedef std::vector<UserColor> UserColorList;
typedef WeightColor<UserColorType> UserWColor;
typedef VectorInterface<UserColorType> UserColorVecBase;

//带有颜色信息的点,里面存储的坐标表示的是dom平面上的坐标，而不是它在图片上的坐标
class ColorDomPt : public TempPtMid
{
public:
    ColorT ptColor_;
    //权值
    double weight_=0;
    //对应点在图片上的坐标
    TempPtMid imgPt_;
public:
    //传入点坐标的构造函数
    ColorDomPt(ColorT ptColor,double xLoc=0,double yLoc=0) : TempPtMid(xLoc,yLoc)
    {
        //记录当前点的颜色
        ptColor_=ptColor;
    }

    //对颜色做初始化
    virtual void initColor()
    {
        for(int dimId=0;dimId<3;++dimId) ptColor_[dimId]=0;
    }

    //初始化当前对象的颜色
    virtual void initColor(const ColorT &otherColor)
    {
        //记录颜色
        ptColor_=otherColor;
    }

    //不传入颜色的构造函数
    ColorDomPt(double xLoc=0,double yLoc=0) : TempPtMid(xLoc,yLoc)
    {
        //对颜色做一个初始化
        initColor();
    }

    virtual void operator=(const TempPtBase& otherPt)
    {
        TempPtMid::operator=(otherPt);
    }

    //传入基类的点的构造函数,调用这个函数的时候通常指的是对应点在图片上的坐标
    ColorDomPt(const TempPtBase& otherPt)
    {
        //调用等号操作
        operator=(otherPt);
        //对应点在图片上的坐标也做同样的记录
        imgPt_=otherPt;
        //初始化颜色
        initColor();
    }
};

//颜色点信息的列表
typedef std::vector<ColorDomPt> ColorPtList;
//直线的两个端点
typedef std::array<CornerPt,2> EndPoints;
//两个直线的端点在点云点集中对应的点标号，有时候还是有用的
typedef std::array<IndexT,2> EndIndexes;

//实现直线相关运算需要使用的虚接口
//怕以后又要复用这个东西，还是把它独立出来吧
class LineInterface
{
public:
    //纯虚函数,获取直线的方程
    virtual const LineT& getLineEqu()=0;

    //传入y值获取x值
    virtual double getX(double fromY)
    {
        //获取直线方程
        const LineT& mainLine=getLineEqu();
        //判断这个方程能不能用来计算x
        if(mainLine(0)==0)
        {
            std::cout<<"a horizontal line cannot be used to compute x coordinate\n";
            return -1;
        }
        //根据直线方程的规则计算x的值
        return (-mainLine(2)-mainLine(1)*fromY)/mainLine(0);
    }

    //传入x值获取y值
    virtual double getY(double fromX)
    {
        //获取直线方程
        const LineT& mainLine=getLineEqu();
        //判断这个方程能不能用来计算y
        if(mainLine(1)==0)
        {
            std::cout<<"a vertical line cannot be used to compute y coordinate\n";
            return -1;
        }
        //根据直线方程计算y值
        return (-mainLine(2)-mainLine(0)*fromX)/mainLine(1);
    }

    //计算传入的直线与当前直线的交点
    virtual TempPtBase lineCross(const LineT &otherLine)
    {
        //获取当前直线的方程
        const LineT &thisLine=getLineEqu();
        //计算直线交点的齐次形式
        HomoPt homoPt=thisLine.cross(otherLine);
        //返回正常的坐标形式
        TempPtBase retPt;
        //判断是否为无穷远点
        if(homoPt[2]==0)
        {
            std::cout<<"the cross point is in infinte\n";
            retPt<<0,0;
        }
        else
        {
            //转换为正常的欧氏点
            retPt<<homoPt[0]/homoPt[2],homoPt[1]/homoPt[2];
        }
        return retPt;
    }

    //获取直线的方向向量，其实可以直接通过直线方程算出来的，暂时用不着这样，懒得写了
    virtual void getDir(DomDir& dstDir)
    {
        //获取直线方程
        const LineT& thisLine=getLineEqu();
        //判断是否为合法直线
        if(thisLine[0]==0&&thisLine[1]==0)
        {
            std::cout<<"this line has an error direction\n";
            return;
        }
        //判断x或y的系数是否存在0
        if(thisLine[0]==0)
        {
            dstDir.initPoint(1,0);
            return;
        }
        else if(thisLine[1]==0)
        {
            dstDir.initPoint(0,1);
            return;
        }
        //用于计算结果的两个点
        TempPtMid pt1,pt2;
        //判断是否过原点
        if(thisLine[2]==0)
        {
            pt1.initPoint(0,0);
            pt2.initPoint(1,getY(1));
        }
        else
        {
            //两个截距
            pt1.initPoint(0,getY(0));
            pt2.initPoint(getX(0),0);
        }
        //计算最后的结果
        dstDir=pt1-pt2;
    }

    //计算两个向量的余弦夹角，感觉没什么用
    virtual double cosAngle(LineInterface& anotherLine)
    {
        //获取两个向量的方向向量
        DomDir thisDir,otherDir;
        getDir(thisDir);
        anotherLine.getDir(otherDir);
        //返回两个向量的夹角余弦
        return thisDir.cosAngle(otherDir);
    }

    //与传入的直线计算夹角,传入的是一个实现了直线接口的东西，直接传入直线还得自己
    //实现一下通过直线计算方向向量
    //最后得到的正弦值是非负数
    virtual double sinAngle(LineInterface& anotherLine)
    {
        //获取两个向量的方向向量
        DomDir thisDir,otherDir;
        getDir(thisDir);
        anotherLine.getDir(otherDir);
        //计算两个向量之间的夹角正弦并返回
        return thisDir.sinAngle(otherDir);
    }

    //计算点到直线的距离
    virtual double pointDistance(const TempPtBase& otherPt)
    {
        //获取直线方程
        const LineT &thisLine=getLineEqu();
        //直线的前两个系数的模
        double lineNorm=std::sqrt(std::pow(thisLine[0],2)+std::pow(thisLine[1],2));
        //传入点的齐次形式
        HomoPt homoPt;
        homoPt<<otherPt[0],otherPt[1],1.f;
        //按照公式计算距离
        return std::abs((thisLine.transpose()*homoPt)[0]/lineNorm);
    }

    //直线的绝对角度
    virtual double absoluteAngle()
    {
        //获取直线的向量
        DomDir thisDir;
        getDir(thisDir);
        //返回标准化后的向量的绝对角度,有点eigen内味儿了
        return thisDir.vectorNorm().absoluteAngle();
    }
};

//一个中间层的继承，为了用上直线的运算接口，感觉是一个很不优雅的补丁
//双继承，继承直线，继承直线的运算接口
class LineMid : public LineT, public LineInterface
{
public:
    //重载等号
    void operator=(const LineT &otherLine)
    {
        LineT::operator=(otherLine);
    }

    //构造函数
    LineMid(const LineT& otherLine)
    {
        //调用等号
        LineMid::operator=(otherLine);
    }

    //获取直线方程
    virtual const LineT & getLineEqu() override
    {
        //直接返回自己
        return *this;
    }
};

//dom平面上的一个直线
//这里是存在可优化的空间的，三角面片是存在共用同一个边的情况的，但目前的算法没有考虑这个情况
//到后面看具体的运行速度再决定
class DomLine : public LineInterface
{
protected:
    //dom直线的长度,当被计算过一次之后，就不再计算了
    double lineLen_=-1;
    //当前直线对应的二维直线方程
    LineT lineEqu_;
    //直线方程初始化过的标志
    bool lineInitFlag=false;
public:
    //直线的两个端点
    EndPoints endPts_;

    //初始化这条直线,方便给构造函数调用,指定由点云中的哪两个点生成当前位置的直线
    virtual void initLine(const Landmarks& cloudPoints,IndexT idx1,IndexT idx2)
    {
        //分别构造第1个点和第2个点
        endPts_[0].initPoint(cloudPoints,idx1);
        endPts_[1].initPoint(cloudPoints,idx2);
    }

    //通过传入两个点来初始化直线
    virtual void initLine(const TempPtBase& pt1,const TempPtBase& pt2)
    {
        //初始化两个直线
        endPts_[0]=pt1;
        endPts_[1]=pt2;
    }

    //DOM直线的构造函数，需要传入两个点
    DomLine(CornerPt &pt1, CornerPt &pt2)
    {
        //分别给两个商战赋值
        endPts_[0]=pt1;
        endPts_[1]=pt2;
    }

    //也可以不初始化，那就直接是一个空的直线
    DomLine()
    {

    }

    //重载一下赋值操作
    virtual void operator=(const DomLine& anotherLine)
    {
        endPts_=anotherLine.endPts_;
        //这样访问都不报错的吗?
        lineLen_=anotherLine.lineLen_;
        //记录直线方程的相关内容
        lineEqu_=anotherLine.lineEqu_;
        lineInitFlag=anotherLine.lineInitFlag;
    }

    //获取线段的长度
    double getLength()
    {
        if(lineLen_>0)
            return lineLen_;
        lineLen_=std::sqrt(std::pow(endPts_[0][0]-endPts_[1][0],2)+
                std::pow(endPts_[0][1]-endPts_[1][1],2));
        return lineLen_;
    }

    //获取直线方程
    virtual const LineT& getLineEqu() override
    {
        //如果没计算过就先计算一下
        if(!lineInitFlag)
        {
            //获取两个端点的齐次形式
            HomoPt homoPt1,homoPt2;
            endPts_[0].gethomoModel(homoPt1);
            endPts_[1].gethomoModel(homoPt2);
            //计算直线方程的过程
            lineEqu_=homoPt1.cross(homoPt2);
            //标记已经计算过
            lineInitFlag=true;
        }
        //返回结果
        return lineEqu_;
    }

    //获取两个点之间的向量
    virtual void getDir(DomDir& dstDir) override
    {
        dstDir=static_cast<TempPtBase>(endPts_[1]-endPts_[0]);
    }

    //线段在xy坐标上的最大投影长度
    virtual double maxProjectXy()
    {
        //计算两个端点之间的差
        DomDir pointDiff;
        //这里禁止调用子类的getDir
        DomLine::getDir(pointDiff);
        //返回端点差的最大值
        return std::max(std::abs(pointDiff[0]),std::abs(pointDiff[1]));
    }

    //获取线段的范围,其实这个功能可以单独做一个range的类，懒得弄了
    virtual void getRange(Range2D& dstRange) const
    {
        //初始化坐标
        dstRange.initRange(endPts_[0][0],endPts_[0][1]);
        //添加坐标
        dstRange.addPoint(endPts_[1][0],endPts_[1][1]);
    }

    //判断传入的点是否在两个端点之间
    //不负责判断点是否在该直线上，仅仅从区间上判断在两个点之间的可能性
    virtual bool judgeBetween(TempPtBase& otherPt) const
    {
        //获取两个点的范围
        Range2D tempRange;
        getRange(tempRange);
        //判断是否在这个范围内
        return tempRange.judgeIn(otherPt[0],otherPt[1]);
    }
};

//处理dom问题时在图片上的一个直线
class DomImgLine : public DomLine
{
public:
    IndexT imgId_;//直线所属的图片标号
    //图片上的颜色点序列
    ColorPtList colorPtList_;
public:
    //空的构造函数
    DomImgLine() : DomLine()
    {
        //图片的标号初始化成0，虽然初始化成0它也不能用
        imgId_=0;
    }

    //最后一个参数传入的是图片标号
    DomImgLine(CornerPt &pt1,CornerPt& pt2,IndexT imgId) : DomLine(pt1,pt2)
    {
        //记录图片标号
        imgId_=imgId;
    }

    //初始化图片的构造函数,用点云中的两个点在特定图片下的投影来生成一条直线
    void initImgLine(const Landmarks& cloudPoints, IndexT pt1, IndexT pt2,IndexT imgIdx)
    {
        //记录图片的标号
        imgId_=imgIdx;
        //根据图片初始化线段的两个点
        endPts_[0].initPoint(cloudPoints,pt1,imgIdx);
        endPts_[1].initPoint(cloudPoints,pt2,imgIdx);
    }

    virtual void operator=(const DomImgLine& anotherLine)
    {
        //调用父类的部分
        DomLine::operator=(anotherLine);
        //复制图片的标号
        imgId_=anotherLine.imgId_;
    }

    //复制另一个直线的信息
    DomImgLine(const DomImgLine& anotherLine) : DomLine()
    {
        //调用赋值操作
        operator=(anotherLine);
    }

    //需要保证向量的x分量为1
    virtual void getDir(DomDir &dstDir) override
    {
        DomLine::getDir(dstDir);
        //需要考虑x=0的情况，不然会崩
        if(dstDir[0]==0)
        {
            //验证y值能否归一化
            if(dstDir[1]==0) return;
            else
            {
                dstDir/=std::abs(dstDir[1]);
                return;
            }
        }
        //确保第1个分量的绝对值是1
        dstDir/=std::abs(dstDir[0]);
    }

    //获取需要遍历的坐标列表,把坐标列表存储到传入的目标容器中
    //传入的数据类型PointT需要能够接受TempPtBase作为构造函数的传入变量
    //新找到的点会被添加到dstPtList上面，旧的点不会被删除
    template<typename PointT>
    void getIterPtList(std::vector<PointT>& dstPtList)
    {
        //获取从第1个点到第2个点的方向向量
        DomDir lineDir;
        getDir(lineDir);
        //判断直线方向是否合法
        if(lineDir[0]==0&&lineDir[1]==0)
        {
            std::cout<<"getIterPtList error, a zero dom vector\n";
            return;
        }
        //遍历y值的时候，y值的增加方向
        int yDir;
        //判断y方向是否有变化
        if(lineDir[1]==0) yDir=0;
        else
            yDir=lineDir[1]/std::abs(lineDir[1]);
        //获取第1个点
        TempPtBase firstPt=endPts_[0];
        //提前开辟空间
        dstPtList.reserve((int)maxProjectXy());
        //记录第1个点
        dstPtList.push_back(PointT(firstPt));
        while(true)
        {
            //获取下一个位置的点
            TempPtBase nextPt=firstPt+lineDir;
            //是否迭代完这次就结束
            bool outFlag=false;
            //判断是否越界,乘上lineDir[0]是为了保证方向一致
            if(nextPt[0]*lineDir[0]>=endPts_[1][0]*lineDir[0])
            {
                //令下一个点等于终点
                nextPt=endPts_[1];
                //迭代完这一次就结束
                outFlag=true;
            }
            //记录下一个点
            dstPtList.push_back(PointT(nextPt));
            //判断y方向是否有变化
            if(yDir!=0)
            {
                //y迭代范围,这个地方的写法犹豫了好久
                int yEnd=std::floor(nextPt[1]);
                int yBegin=std::floor(firstPt[1])+yDir;
                //遍历起点和终点之间的y值
                for(int yIter=yBegin;yIter*yDir<=yEnd*yDir;yIter+=yDir)
                {
                    //计算当前位置的x坐标
                    double xCoord=getX(yIter);
                    //点坐标形式
                    TempPtBase tempPt;
                    tempPt<<xCoord,yIter;
                    //记录坐标
                    dstPtList.push_back(PointT(tempPt));
                }
            }
            //是否结束
            if(outFlag) break;
            //更新点的位置
            firstPt=nextPt;
        }
    }

    //准备填充颜色点列表,传入true的时候无论如何都会被初始化
    void prepareColorPtList(bool prepareAnyway=false)
    {
        //判断是否准备过
        if(prepareAnyway||colorPtList_.size()==0)
        {
            colorPtList_.clear();
            getIterPtList<ColorDomPt>(colorPtList_);
        }

    }
};

//参数方程形式的直线
//目前仅限单参数的直线，多参数的直线好像没什么意义
//它的上层其它可以有一个父类的直线，也可以有一个父类的参数曲线，暂时用不到
//用到了再说
class ParameterLine
{
public:
    //传入参数，然后获取对应的xy坐标
    virtual bool getXy(double par,double &x,double &y) const=0;
};

//线性运算方程
//可以继承直线接口，暂时不继承，等用到了再继承也不妨
class LinearFunc
{
public:
    //斜率和截距
    double k_;
    double b_;
public:
    LinearFunc(double k=1,double b=0) : k_(k),b_(b)
    {

    }

    //获取y坐标
    double getY(double xValue) const
    {
        return k_*xValue+b_;
    }
};

//这是一个带参数的直线，这样的一个直线是由参数方向构成的
//它当然也可以继承line接口，但暂时没必要
//关于z的参数方程，都是常数除以一个与z相关的分母然后再加一个常数的操作
class ProjectLine : public ParameterLine
{
protected:
    //分母对应的线性函数
    LinearFunc denFunc_;
    std::array<double,2> xPar_;//x的运算参数，第1个参数是分子，第二个参数是常数加
    std::array<double,2> yPar_;

    //初始化分母
    void initDen(double r31,double r32,double r33,//旋转矩阵第3行
                 double cx,double cy,double cz,//光心坐标
                 double x,double y)//dom平面上的坐标
    {
        denFunc_.k_=r33;
        denFunc_.b_=r31*(x-cx)+r32*(y-cy)-r33*cz;
    }

    //初始化计算参数
    void initXyPar(double rx1,double rx2,double rx3,//旋转矩阵行向量x是第1行，y是第2行
                   double cx,double cy,double cz,//相机光心
                   double x,double y,//xy坐标和主点坐标
                   double pri,double f,//主点坐标,焦距
                   std::array<double,2> &dstPar)//参数的目标结果
    {
        //常数参数
        dstPar[1]=rx3/denFunc_.k_;
        //分子参数
        dstPar[0]=rx1*(x-cx)+rx2*(y-cy)-rx3*cz-dstPar[1]*denFunc_.b_;
        //把计算结果乘上f
        dstPar[0]*=f;
        //常数参数再另外加上主点
        dstPar[1]=f*dstPar[1]+pri;
    }

    //计算具体的x值或y值
    double computeByPar(const std::array<double,2>& par,double denValue) const
    {
        //分别是分子除以它再加个常数
        return par[0]/denValue+par[1];
    }
public:
    IndexT imgId_;//图片标号
    bool isUsable_=true;

    //传入参数，获取最后的xy结果
    //这里的参数就是z值
    virtual bool getXy(double par, double &x, double &y) const override
    {
        //计算分母数值
        double denValue=denFunc_.getY(par);
        //判断分母为0
        if(denValue==0)
        {
            return false;
        }
        //分别计算xy
        x=computeByPar(xPar_,denValue);
        y=computeByPar(yPar_,denValue);
        return true;
    }

    //空的构造函数，方便通过其它形式载入
    ProjectLine(){}

    //构造函数
    ProjectLine(const Eigen::Matrix3d& rMat,//旋转矩阵
                const Eigen::Vector3d &camCen,//光心
                const std::vector<double> &intrPar,//内参,焦距,主点x,主点y
                double domX,double domY,IndexT imgId) : imgId_(imgId)
    {
        //初始化分母
        initDen(rMat(2,0),rMat(2,1),rMat(2,2),camCen(0),camCen(1),camCen(2),domX,domY);
        //初始化分子的常数参数
        initXyPar(rMat(0,0),rMat(0,1),rMat(0,2),camCen(0),camCen(1),camCen(2),
                  domX,domY,intrPar[1],intrPar[0],xPar_);
        initXyPar(rMat(1,0),rMat(1,1),rMat(1,2),camCen(0),camCen(1),camCen(2),
                  domX,domY,intrPar[2],intrPar[0],yPar_);
    }

    //把当前直线的信息保存进二进制文件中
    void saveBinaryFile(std::fstream& fileHandle)
    {
        //记录图片标号
        fileHandle.write((char*)&imgId_,sizeof(imgId_));
        //记录分母线性函数的斜率和截距
        fileHandle.write((char*)&denFunc_.k_,sizeof(double));
        fileHandle.write((char*)&denFunc_.b_,sizeof(double));
        //记录决定x的参数
        fileHandle.write((char*)xPar_.data(),sizeof(double)*2);
        fileHandle.write((char*)yPar_.data(),sizeof(double)*2);
    }

    //从输入流中载入二进制文件信息
    void loadByBinaryStream(std::fstream& fileHandle)
    {
        //读取图片标号
        fileHandle.read((char*)&imgId_,sizeof(imgId_));
        //读取分母线性函数的斜率和截距
        fileHandle.read((char*)&denFunc_.k_,sizeof(double));
        fileHandle.read((char*)&denFunc_.b_,sizeof(double));
        //读取xy的参数
        fileHandle.read((char*)xPar_.data(),sizeof(double)*2);
        fileHandle.read((char*)yPar_.data(),sizeof(double)*2);
    }

    //从二进制文件中载入所有的直线数据
    static void loadProjectLines(std::vector<ProjectLine>& dstProjectLineGroup,IndexT domUnitLabel)
    {
        //获取二进制文件的路径
        std::string binFilePath=projectBinFileName(domUnitLabel);
        //打开二进制文件
        std::fstream binFileHandle;
        binFileHandle.open(binFilePath,std::ios::binary);
        //读取数据
        while(!binFileHandle.eof())
        {
            //添加一个新的直线
            dstProjectLineGroup.push_back(ProjectLine());
            //载入二进制文件信息
            dstProjectLineGroup.back().loadByBinaryStream(binFileHandle);
        }
        //关闭文件
        binFileHandle.close();
    }
};
//投影直线的列表
typedef std::vector<ProjectLine> ProjLineList;
typedef std::vector<ProjLineList> ProjectLineListGroup;
//适用于快速DOM的一个数据类型，每个图片标号对应一个直线
typedef std::vector<DomImgLine> ViewImgLines;

//对三角面片做稠密化的时候，需要往里面添加要遍历的点，要遍历的边在每个图以及DOM图中的直线方程
class PatchIterInfo
{
public:
    //需要被迭代的直线，每个图上都有一个直线,第1个图就是主图
    //找稠密点的时候都会从主图出发去寻找
    //图片的存储顺序没有意义，乱序存储的
    ViewImgLines imgLines_;
    //dom平面上的直线,每个图片的直线都和这个东西对应
    DomLine domLine_;

    //感觉c++默认的enum不好用，所以用了struct
    struct ClassEnum
    {
        //最原始的基类
        static const unsigned int Base=0;
        //每个patch的主迭代器
        static const unsigned int MainIter=1;
    };

    //用于标记当前的对象来自哪个class
    virtual unsigned int classId()
    {
        return ClassEnum::Base;
    }
};

//每个patch的全局迭代器，有一些专属的信息
//父类里面是最长的边以及它的对应边列表
class PatchMainIterInfo : public PatchIterInfo
{
public:
    //来自三角面片，除了一条直线，应该还有剩下的一个点
    CornerPt otherCorner_;

    //标记当前的对象来自哪一层的继承
    virtual unsigned int classId() override
    {
        return ClassEnum::MainIter;
    }
};

//二维的三角面片，里面存储的仅仅是点的索引，具体的点坐标需要去sfm的点里面找
class TriangularPatch2D
{
public:
    //点的索引列表
    Pt3Group corners_;
    //能观察到这三个点的图片ID的列表
    IndexList obvs_;
    //如果是一个已经被放弃治疗的面片，设置成true
    bool isInvalid_=false;
public:
    //构造函数，传入三角面片三个点在点云里面的下标
    TriangularPatch2D(IndexT pt1=0,IndexT pt2=0,IndexT pt3=0)
    {
        corners_[0]=pt1;
        corners_[1]=pt2;
        corners_[2]=pt3;
    }

    //maxImgLine表示的是最长的一条直线，里面包括它的两个端点以及它来自哪个图片
    virtual void getFirstPatchIterInfo(PatchMainIterInfo& dstInfo,const Landmarks& cloudPoints,
                               const DomImgLine& maxImgLine, IndexT maxLineIdx)
    {
        //确定迭代直线的对角点,初始化它的点
        dstInfo.otherCorner_.initPoint(cloudPoints,corners_[(maxLineIdx+2)%3]);
        //直线的另一端的点标号
        IndexT nextLinePtIdx=(maxLineIdx+1)%3;
        //记录它点云点中的主直线
        dstInfo.domLine_.initLine(cloudPoints,corners_[maxLineIdx],corners_[nextLinePtIdx]);
        //给每个线提前开辟空间
#ifdef USE_ONLY_IMG //仅使用最长的直线
        dstInfo.imgLines_.resize(1);
#else
        dstInfo.imgLines_.resize(obvs_.size());
#endif
        //记录第1个直线
        dstInfo.imgLines_[0]=maxImgLine;
#ifndef USE_ONLY_IMG
        //初始化当前的迭代位置
        unsigned int iterLoc=1;
        //遍历能观察到这三个点的其它图片，给每个图片注册一条直线
        for(IndexList::iterator iter=obvs_.begin();iter!=obvs_.end();++iter)
        {
            //已经添加过的第1个直线就不用添加了
            if(*iter==maxImgLine.imgId_) continue;
            //初始化当前位置的图片直线,每次都需要反复访问一个固定的index，存在优化空间
            dstInfo.imgLines_[iterLoc].initImgLine(cloudPoints,corners_[maxLineIdx],
                                                   corners_[nextLinePtIdx],*iter);
            //提前计算一下直线方程，防止后面多线程访问的时候冲突
            dstInfo.imgLines_[iterLoc].getLineEqu();
            //迭代位置计数增加
            ++iterLoc;
        }
#endif
    }

    //获取接下来要迭代的信息,传入的是点云信息
    virtual void getFirstPatchIterInfo(PatchMainIterInfo& dstInfo,const Landmarks& cloudPoints)
    {
        //初始化过的标志，如果还没有被初始化过，就直接记录下第1个点
        bool initFlag=false;
        //最长的一条线
        DomImgLine maxImgLine;
        //最长的线对应的线段标号,其实存储的是角点编号，当前角点和下一个角点的连线就是所谓的线
        IndexT maxLineIdx=0;
        //遍历三角面片的每个边
        for(int cornerId=0;cornerId<3;++cornerId)
        {
            //当前位置的点
            const Landmark &thisCloudPt=cloudPoints.at(corners_[cornerId]);
            //下一个位置的点标号
            IndexT nextPtIdx=corners_[(cornerId+1)%3];
            //下一个位置的点云点
            const Landmark &nextCloudPt=cloudPoints.at(nextPtIdx);
            //遍历能观察到这3个点的每个图片id
            for(IndexList::iterator iter=obvs_.begin();iter!=obvs_.end();++iter)
            {
                //当前位置的两个投影点
                const double* const thisObvCArray=thisCloudPt.obs.at(*iter).x.data();
                CornerPt thisObvPt(corners_[cornerId],thisObvCArray[0],thisObvCArray[1]);
                const double* const nextObvCArray=nextCloudPt.obs.at(*iter).x.data();
                CornerPt nextObvPt(nextPtIdx,nextObvCArray[0],nextObvCArray[1]);
                //用两个投影点构造一个直线
                DomImgLine tempLine(thisObvPt,nextObvPt,*iter);
                //获取线段的长度
                double lineLen=tempLine.getLength();
                //记录这个东西的标志
                bool saveFlag=false;
                //判断有没有初始化过
                if(!initFlag)
                {
                    //初始化的标记更改
                    initFlag=true;
                    //记录当前的信息
                    saveFlag=true;
                }
                else if(lineLen>maxImgLine.getLength())//判断这是不是最长的一个线
                {
                    saveFlag=true;
                }
                //是否记录
                if(saveFlag)
                {
                    //重载过等号运算符,直接记录直线的相关信息
                    maxImgLine=tempLine;
                    //最长的线对应的线段标号
                    maxLineIdx=cornerId;
                }
            }
        }
        //如果最长的直线长度太短，那就算了
        if(maxImgLine.getLength()<=2)
        {
            isInvalid_=true;
            return;
        }
        //找到最长的线和最长的线对应的直线标号之后，正式开始构造第一次的迭代信息
        getFirstPatchIterInfo(dstInfo,cloudPoints,maxImgLine,maxLineIdx);
    }

    //三个点，这个图只看到了两个点，这个图没有看到哪个点
    //返回0,1,2中的一个
    int whoLostThisImg(Landmarks& cloudPoints, IndexT lostImg)
    {
        //遍历3个点
        for(int ptId=0;ptId<3;++ptId)
        {
            //如果看不到，就返回这个标号
            if(cloudPoints.at(corners_[ptId]).obs.count(lostImg)==0)
                return ptId;
        }
        //不正常的情况
        std::cout<<"this img can see all of the patch\n";
        return -1;
    }

    //把当前面片能看到的三个点存储为map的形式，仅存储里面的标号
    void getCornerMap(Pt3Map& dstMap) const
    {
        //遍历三个角点
        for(int ptCount=0;ptCount<3;++ptCount)
            dstMap[corners_[ptCount]]=true;
    }
};

//20210926
//DOM像素的法向量方向
//DOM向量仅仅传播Z值似乎是不够的，还需要传播法向量
//这个法向量与其它的法向量不再发生交叉继承关系
//它这里表示的其实是一个半球上的各种数值
class PixelDir
{
private:
    //向量的x分量和y分量，最后其实是用来表示一个-1~1的数字
    char x_=0;
    char y_=0;

    //数据的最大值
    static const int maxValue_=128;
public:
    //获取向量
    void makeDir(Point3DBase& dstDir) const
    {
        //分别计算x值和y值
        dstDir[0]=(double)x_/maxValue_;
        dstDir[1]=(double)y_/maxValue_;
        //根据前两个数字计算第3个值
        //z值一定是正的
        dstDir[2]=std::sqrt(1.f-dstDir[0]*dstDir[0]-dstDir[1]*dstDir[1]);
    }

    //从一个向量中获取半球形式的向量
    void fromDir(const Point3DBase& otherDir)
    {
        //复制一份数据
        Point3DBase tempDir=otherDir;
        //正则化
        tempDir.normalize();
        //如果z值小于0,就取反
        if(tempDir[2]<0) tempDir=-tempDir;
        //对x,y数据赋值
        x_=tempDir[0]*maxValue_;
        y_=tempDir[1]*maxValue_;
    }

    //从其它的向量里面复制
    PixelDir& operator=(const PixelDir& otherDir)
    {
        //记录坐标值
        x_=otherDir.x_;
        y_=otherDir.y_;
        //返回自身
        return *this;
    }
};

//数学意义上的平面的接口
class PlaneInterface
{
protected:
    //平面方程
    PlaneT planeEqu_;

    //平面方程是否已经初始化过了
    //本来是想在构造函数里面调用update()的，但在构造函数里面调用虚函数是很危险的行为
    //所以设立了这样的一个flag，需要由子类来调用初始化函数
    bool initFlag_=false;
public:
    //需要从子类里面获取若干个点来初始化这个平面
    virtual void getPlanePts(std::vector<Point3D>& dstPtList) const=0;

    //获取点集的中心点
    //可以被子类重载然后用其它的方式获取
    virtual void getCenterPt(Point3D& dstPt)
    {
        //获取所有的点
        std::vector<Point3D> ptList;
        std::vector<VecInterDouble*> ptrList;
        VecInterDouble::makeInterfacePtr<Point3D>(ptList,ptrList);
        //获取中心点
        VecInterDouble::getAverageValue(ptrList,dstPt);
    }

    //用一个法向量和一个点来制作一个平面
    static void makePlaneByNormPt(PlaneT& dstPlane,const Dir3D& normal,const Point3D& pt)
    {
        //判断传入的向量是不是一个纯零的向量
        if(normal.judgeAllValue(0)) throw ERROR_PLANE_COLLINEATION;
        //将法向量的前3个数值写入到平面里面
        for(int vecCount=0;vecCount<3;++vecCount) dstPlane[vecCount]=normal[vecCount];
        //计算它的最后一个数值
        dstPlane[3]=-normal.dotProduct(pt);
    }

    //用svd的方法计算法向量，主要是为了配合下面的平面方程的计算
    static void computeNormalSvd(Dir3D& dstDir,const std::vector<Point3D>& ptList)
    {
        //计算点集们两两相减形成向量的过程
        std::vector<Dir3D> eachDir;
        PtDirEachMinus minusObj;
        minusObj.run(ptList,eachDir);
        //把向量列表转换成父类的指针列表
        std::vector<VecInterDouble*> ptrList;
        VecInterDouble::makeInterfacePtr<Dir3D>(eachDir,ptrList);
        //用向量的集合计算法向量
        VecInterDouble::getMinEigen<3>(ptrList,dstDir);
    }

    //使用SVD方法计算的平面方程
    static void computePlaneEquationSvd(PlaneT& dstPlane,const std::vector<Point3D>& ptList)
    {
        //判断是否超过了3个点
        if(ptList.size()<3) throw ERROR_PLANE_PT_LESS_3;
        //用svd的方法计算法向量
        Dir3D normDir;
        computeNormalSvd(normDir,ptList);
        //计算所有点的中心
        Point3D centerPt;
        std::vector<const VecInterDouble*> ptrList;
        VecInterDouble::makeInterfaceConstPtr<Point3D>(ptList,ptrList);
        VecInterDouble::getAverageValue(ptrList,centerPt);
        //用法向量和所有点的中心计算最终的平面方程
        makePlaneByNormPt(dstPlane,normDir,centerPt);
    }

    //根据传入的点列表计算平面方程
    //目前仅仅使用前三个点做一个简单的计算，如果真的需要特别多的点来拟合平面的时候
    //需要另外找一个最小二乘法来完成
    static void computePlaneEquation(PlaneT& dstPlane,const std::vector<Point3D>& ptList)
    {
        //判断是否超过了3个点
        if(ptList.size()<3) throw ERROR_PLANE_PT_LESS_3;
        //如果超过3个点，使用svd的方法计算
        if(ptList.size()>3)
        {
            computePlaneEquationSvd(dstPlane,ptList);
            return;
        }
        //根据传入的点新建两个三维向量
        Dir3D dir1(ptList[0],ptList[1]);
        Dir3D dir2(ptList[0],ptList[2]);
        //计算法向量
        Dir3D normDir;
        normDir.Point3DBase::operator=(dir1.cross(dir2));
        //用一个法向量和一个点来构造一个平面
        //随便找一个点就行
        makePlaneByNormPt(dstPlane,normDir,ptList[0]);
    }

    //将当前位置的法向量随机更新一个数值
    void updateRandNorm(double normStep=NORM_STEP)
    {
        //获取当前的法向量
        Dir3D normVec;
        getNormal(normVec);
        //随机生成一个三维的法向量
        Dir3D randVec;
        while(true)
        {
            //随便生成一个eigen向量
            Point3DBase tempRand=Point3DBase::Random(3,1);
            //记录到向量中
            randVec.Point3DBase::operator=(tempRand);
            //如果是全零向量，就继续随机
            if(!randVec.judgeAllValue(0)) break;
        }
        //将得到的随机向量归一化
        randVec.normalize();
        //把临时得到的向量叠加到目前的法向量上
        normVec.Point3DBase::operator+=(randVec*normStep);
        //按照eigen库里面的东西来操作吧
        Point3DBase& normBase=normVec;
        //把它归一化
        normBase.normalize();
        //判断x,y是否超界
        if(normBase[2]<0) normBase=-normBase;
        if(normBase[2]<MIN_Z_WEIGHT)
        {
            normBase<<0,0,1;
        }
        //获取平面的中心点
        Point3D centerPt;
        getCenterPt(centerPt);
        //使用新的法向量和中心点初始化平面
        updatePlane(normVec,centerPt);
    }


    //更新这个平面的数据
    void updatePlane()
    {
        //从子类里面获取若干平面上的点数据
        std::vector<Point3D> ptList;
        getPlanePts(ptList);
        //根据若干点计算的平面方程
        computePlaneEquation(planeEqu_,ptList);
    }

    //指定使用特定点和特定直线来更新平面方程
    void updatePlane(const Dir3D& normal,const Point3D& centerPt)
    {
        //判断是否为全零的法向量
        if(normal.Point3DBase::norm()==0)
        {
            throw ERROR_ZERO_NORMAL;
        }
        //调用计算接口
        makePlaneByNormPt(planeEqu_,normal,centerPt);
    }

    //给定另外两个维度的坐标，然后计算它在第3个维度上的值
    double computeValue(double xValue,double yValue,IndexT dim=2) const
    {
        //判断传入的维度是否正常
        if(dim>2) throw ERROR_PLANE_DIM_OUT;
        //判断目标维度的数据是否可解
        if(planeEqu_[dim]==0) throw ERROR_MANY_VALUE;
        //拿出来平面上的两个数据
        std::array<double,2> planeValue;
        //如果是第1个维度，需要把顺序反过来
        if(dim==1)
        {
            planeValue[0]=planeEqu_[0];
            planeValue[1]=planeEqu_[1];
        }
        else //按照取余的方式刚好可以做
        {
            planeValue[0]=planeEqu_[(dim+1)%3];
            planeValue[1]=planeEqu_[(dim+2)%3];
        }
        //计算最后的结果
        return (-planeEqu_[3]-planeValue[0]*xValue-planeValue[1]*yValue)/planeEqu_[dim];
    }

    //计算给定位置的z坐标
    double computeZ(double xValue,double yValue) const
    {
        //计算第3个维度上的结果
        return computeValue(xValue,yValue,2);
    }

    //批量计算z值，本来是可以弄一个批量计算某个维度的值的
    //但目前没必要，懒得弄了，追求代码复用性真的是个大坑
    void computeManyZ(Point3DList& dstPt) const
    {
        //遍历点列表
        //可以考虑使用多线程
        for(IndexT ptId=0;ptId<dstPt.size();++ptId)
        {
            //当前位置的点
            Point3D& thisPt=dstPt[ptId];
            //计算当前位置的z值
            thisPt[2]=computeZ(thisPt[0],thisPt[1]);
        }
    }

    //获取平面的法向量
    void getNormal(Point3DBase& dstDir) const
    {
        //记录平面方程里面的前3项
        for(int dimCount=0;dimCount<3;++dimCount)
            dstDir[dimCount]=planeEqu_[dimCount];
        //计算向量的模长
        double normValue=dstDir.norm();
        //如果模长是0,那就有问题了
        if(normValue==0) throw ERROR_GET_ZERO_NORM;
        //将向量归一化
        dstDir.normalize();
        //保证z值大于0
        if(dstDir[2]<0) dstDir=-dstDir;
    }

    //获取平面方程
    const PlaneT& getPlaneEqu(){return planeEqu_;}

    //直接设置平面的方程
    void setPlaneEqu(const PlaneT& otherPlane)
    {
        planeEqu_=otherPlane;
    }
};

//二维三角形的接口
template<typename PointType>//pointType必须满足使用operator[]来对其进行操作
class Triangle2DInterface
{
public:
    //获取子类三角形的三维二维点，用来服务相关的计算工作
    //这里面传回来的点是子类的点的拷贝
    virtual void getTrianglePts(std::vector<PointType>& triPts) const=0;

    //获取三角形的二维坐标范围
    void getRange(Range2D& dstRange,const std::vector<PointType>* ptListPtr=nullptr) const
    {
        //判断是否为空指针
        if(ptListPtr==nullptr)
        {
            //获取点的坐标列表
            std::vector<PointType> points;
            getTrianglePts(points);
            //递归调用
            getRange(dstRange,&points);
            return;
        }
        //利用得到的点集初始化点的范围
        dstRange.initRangeT<PointType>(*ptListPtr);
    }

    //将三角形画在图片上，opencv里面有现成的接口
    //画完之后默认是要填充整个区域的
    void drawTriangle(cv::Mat &dstImg,//图片存储的目标位置
                      double pixelLength, //每个像素对应的长度
                      const std::vector<PointType>* ptListPtr=nullptr, //点列表的指针
                      const ResolutionRange2D* rangePtr=nullptr //范围数据的指针
    ) const
    {
        //判断点列表是否为空指针
        if(ptListPtr==nullptr)
        {
            //获取点列表
            std::vector<PointType> points;
            getTrianglePts(points);
            //递归调用
            drawTriangle(dstImg,pixelLength,&points,rangePtr);
            return;
        }
        //判断点列表的个数是否正常
        if(ptListPtr->size()!=3) throw ERROR_TRI_PT_NOT_3;
        //判断是否传入的指针数据
        if(rangePtr==nullptr)
        {
            //根据像素的分布制作一个范围
            ResolutionRange2D triRange;
            //记录每个像素的长度
            triRange.pixelLen_=pixelLength;
            //获取像素的坐标范围
            getRange(triRange,ptListPtr);
            //递归调用
            drawTriangle(dstImg,pixelLength,ptListPtr,&triRange);
            return;
        }
        const ResolutionRange2D& triRange=*rangePtr;
        //用像素的坐标范围初始化图片
        triRange.initCvMat(dstImg,CV_8UC1);
        //把点列表转换为适用于opencv的点列表
        std::vector<std::vector<cv::Point>> corners(1);
        triRange.toCvPt<PointType>(*ptListPtr,corners[0]);
        //在最终的图片上画线
        cv::polylines(dstImg,corners,true,cv::Scalar(255));
        //在最终的图片上填充三角形
        cv::fillPoly(dstImg,corners,cv::Scalar(255));
    }

    //获取三角形区域内的所有点坐标
    //传入的像素长度是点密度的一个参考
    void getPtInTriangle(std::vector<PointType>& dstPtList,
                         double pixelLength) const
    {
        //获取点坐标的列表
        std::vector<PointType> points;
        getTrianglePts(points);
        //根据点坐标列表和对应的分辨率获取三角形的坐标范围
        ResolutionRange2D ptRange;
        ptRange.pixelLen_=pixelLength;
        ptRange.initRangeT<PointType>(points);
        //把三角形区域画在图片上，调用opencv的接口
        cv::Mat triImg;
        drawTriangle(triImg,pixelLength,&points,&ptRange);
        //遍历图片上的每个像素
        for(int rowCount=0;rowCount<triImg.rows;++rowCount)
        {
            //获取当前行的指针
            uchar* rowPtr=triImg.ptr<uchar>(rowCount);
            //遍历每一列
            for(int colCount=0;colCount<triImg.cols;++colCount)
            {
                //判断当前位置是否被写值
                if(*(rowPtr+colCount)==255)
                {
                    //计算当前位置对应的浮点型坐标
                    double srcX,srcY;
                    ptRange.toFloat(colCount,rowCount,srcX,srcY);
                    //生成目标形式的点
                    PointType tempPt;
                    tempPt[0]=srcX;tempPt[1]=srcY;
                    //把临时生成的点添加到目标列表里面
                    dstPtList.push_back(tempPt);
                }
            }
        }
    }


};//class Triangle2DInterface

//目前做快速DOM的时候需要使用到的数据类型
typedef Triangle2DInterface<Point3D> DomTriInterface;

//dom的颜色块信息，每个块上都有若干年颜色，它们都来自同一个图片
//它们到底来片哪个图片，似乎还挺重要的
class ColorPatch : public VectorInterface<double>
{
protected:
    //已经被归一化过的颜色列表
    std::vector<double> normedColor_;

    //每个颜色提供的向量中数值个数
    static const IndexT normBase_=3;

    //计算一个颜色的平均亮度
    static double avgBirght(const cv::Vec3d& srcColor)
    {
        //亮度值累计
        double sumData=0;
        for(IndexT dimCount=0;dimCount<3;++dimCount)
        {
            sumData+=srcColor[dimCount];
        }
        return sumData/3;
    }

    //获取一个颜色的最大绝对值
    static double getMaxAbs(const cv::Vec3d& srcColor)
    {
        //初始化数据
        double maxValue=0;
        for(IndexT dimCount=0;dimCount<3;++dimCount)
        {
            //判断是否为更大的值
            if(std::abs(srcColor[dimCount])>maxValue)
            {
                maxValue=std::abs(srcColor[dimCount]);
            }
        }
        return maxValue;
    }

    //向归一化过的颜色列表里面添加新的颜色
    void addNormColor(cv::Vec3d& addedColor)
    {
        //依次添加三个颜色
        for(int i=0;i<3;++i) normedColor_.push_back(addedColor[i]);
    }

    //原始颜色与均值颜色的运算规则
    //传入的颜色满足用中括号做随机访问就行
    cv::Vec3d transformColorByAvg(const cv::Vec3d& avgColor,
                                  const cv::Vec3d& mainColor)
    {
        //新建待返回的值
        cv::Vec3d retColor;
        //计算主颜色的平均亮度
        double avgValue=avgBirght(mainColor);
        //遍历每个通道
        for(IndexT chnCount=0;chnCount<3;++chnCount)
        {
            //将每个颜色通道与平均亮度相减
            retColor[chnCount]=(mainColor[chnCount]-avgValue);
        }
        //返回结果
        return retColor;
    }
public:
    //颜色数据的平均值被设置为0的标志
    //都弄成public是省得给自己找麻烦
    bool avgFlag_=false;

    //将当前的面片标记为不可作为中心
    bool banCenter_=false;

    //颜色面片的权值，当计算平均的NCC分数的时候，会用到这个权值
    double weight_=1;

    IndexT imgId_=0;//图片像素面片来自的图片标号,这个标号是哈希表里面的索引标号
    IndexT viewLineId_=0;//图片投影直线的标号，这个标号指的是viewLines的顺序标号

    //投影直线里面上界到下界的距离
    float rangeProjDis_=0;

    //点的颜色列表，每个位置对应于dom图上的一个颜色
    std::vector<UserColor> patchColors_;

    //它表示的是在这个图片上的若干投影位置
    //与上面的patchColors_一一对应
    std::vector<TempPtBase> projList_;

    //获取当前面片的颜色列表，需要的是正常的颜色列表
    void getCvColorList(ColorTVec& dstColorVec) const
    {
        dstColorVec.reserve(patchColors_.size());
        //记录所有的颜色
        for(IndexT colorCount=0;colorCount<patchColors_.size();++colorCount)
        {
            dstColorVec.push_back((ColorT)patchColors_[colorCount]);
        }
    }

    //将一个外部数据添加到这个颜色面片上
    void addPatch(const ColorPatch& otherPatch)
    {
        //判断两个颜色面片的大小是否相等
        if(patchColors_.size()!=otherPatch.patchColors_.size()) throw ERROR_ADD_PATCH_DIFF_SIZE;
        //依次遍历每个数据，然后依次添加
        for(uint colorCount=0;colorCount<patchColors_.size();++colorCount)
        {
            //实现两个颜色的相加
            patchColors_[colorCount].VectorInterface::operator+=(
                        otherPatch.patchColors_[colorCount]);
        }
    }

    //将每个颜色的幅值都乘以一个数字，这里是对原始的颜色做操作
    //不是对normedcolor做操作
    void multScale(double scale)
    {
        //遍历每个颜色面片
        for(uint colorCount=0;colorCount<patchColors_.size();++colorCount)
        {
            //做乘法操作
            patchColors_[colorCount].VectorInterface::operator*=(scale);
        }
    }

    //从其它颜色面片里面复制颜色信息
    //只复制每个patch的颜色
    void copyColorPatch(const ColorPatch& otherPatch)
    {
        //初始化颜色的数量
        patchColors_.clear();
        patchColors_.reserve(otherPatch.patchColors_.size());
        //记录权值
        weight_=otherPatch.weight_;
        //遍历每个颜色面片，把它们的颜色依次记录下来
        for(uint colorCount=0;colorCount<otherPatch.patchColors_.size();++colorCount)
        {
            //往面片里面添加颜色数据
            const cv::Vec<UserColorType,3> &currColor=otherPatch.patchColors_[colorCount];
            patchColors_.push_back(UserColor(currColor[0],currColor[1],currColor[2]));
        }
    }

    //用一个数据序列来初始化颜色数据,比如传入27个数据，它每3个是一组颜色，共9个颜色
    void initWithColorSerial(const ColorDataSerial& dataSerial)
    {
        //判断颜色的数据序列是否为3的整数倍，不然是有问题的
        if(dataSerial.size()%normBase_!=0) throw ERROR_NOT_3_TIMES;
        //根据颜色数据的大小提前开辟空间
        patchColors_.reserve(dataSerial.size()/normBase_);
        //遍历颜色序列，3个3个地往里面依次添加数据
        for(uint dataCount=0;dataCount<dataSerial.size();dataCount+=normBase_)
        {
            patchColors_.push_back(UserColor(
                dataSerial[dataCount],dataSerial[dataCount+1],dataSerial[dataCount+2]));
        }
    }

    //实现接口类里面要求的size的大小
    std::size_t size() const override
    {
        return normedColor_.size();
    }

    //实现接口的访问机制
    //用来实现对索引的访问
    virtual double& operator[](unsigned int idx) override
    {
            return normedColor_[idx];
    }

    //const类型的访问机制
    virtual const double& operator[](unsigned int idx) const override
    {
            return normedColor_[idx];
    }

    //必须保留一个无参数的构造函数，因为需要对它使用到vector里面的resize
    ColorPatch(){}

    //获取相机投影点到主点的距离
    double getPriDistance(PriPoint& imgCenter)
    {
        //判断投影点里面是否有值
        if(projList_.size()==0)
        {
            throw ERROR_PROJ_EMPTY;
        }
        //获取中间位置的投影点
        TempPtBase& midPt=projList_[projList_.size()/2];
        //计算与主点的距离
        return std::abs(midPt[0]-imgCenter[0])+std::abs(midPt[1]-imgCenter[1]);
    }

    //获取颜色列表的平均颜色
    //这里不复用VevtorInterface,因为有的加减法没有实现
    void getAvgColor(UserColor& dstColor)
    {
        //初始化平均值
        dstColor.operator*=(0);
        //遍历每个颜色
        for(UserColorList::iterator iter=patchColors_.begin();
            iter!=patchColors_.end();++iter)
            dstColor.operator+=(*iter);
        //除以颜色的总数
        dstColor.operator*=(1.f/patchColors_.size());
#ifdef NCC_AVG_AVG
        //计算颜色平均值
        int rgbAvg=0;
        for(int chnCount=0;chnCount<3;++chnCount) rgbAvg+=dstColor.Vec::operator[](chnCount);
        rgbAvg/=3;
        //赋值颜色平均值
        for(int chnCount=0;chnCount<3;++chnCount) dstColor.Vec::operator[](chnCount)=rgbAvg;
#endif
    }

    //颜色面片中总共的颜色的个数
    uint colorAccount() const
    {
        return patchColors_.size();
    }

    //获取的是归一化后的数据，如果还没有做过均值化，那会报错
    const double* getData() const
    {
        //判断当前的颜色是否已经归一化过
        if(!avgFlag_) throw ERROR_NOT_AVG_ZERO;
        //返回vector里面的数据
        return normedColor_.data();
    }

    //将颜色面片里面的东西填充到向量里面
    void fillNormVec()
    {
        //清空向量数据
        normedColor_.clear();
        //开辟空间
        normedColor_.reserve(normBase_*patchColors_.size());
        //遍历每个颜色
        for(IndexT colorCount=0;colorCount<patchColors_.size();++colorCount)
        {
            //获取当前位置的颜色
            cv::Vec3d& currColor=patchColors_[colorCount];
            //记录颜色
            normedColor_.push_back(currColor[0]);
            normedColor_.push_back(currColor[1]);
            normedColor_.push_back(currColor[2]);
        }
    }

    //把颜色的平均值设置为0
    //同时还把norm也设置成1了，本来是没有这一项的
    //弄上这个是为了减少后面的计算量
    void setAvgZero(const UserColor* avgPtr=nullptr,
                    FileHandler* filePtr=nullptr,//log文件管理器的指针
                    std::string headStr=""//log文件管理器需要的一个头部信息
            )
    {
        //如果已经是平均值为0了，那就不用做了
        if(avgFlag_) return;
        {
            //计算颜色的平均值
            UserColor avgColor;
    #ifdef NCC_RM_AVG
            //判断是否有可参考的平均颜色
            if(avgPtr)
            {
                avgColor.Vec::operator=(*avgPtr);
            }
            else
                getAvgColor(avgColor);
    #else
            avgColor.getR()=0;
            avgColor.getG()=0;
            avgColor.getB()=0;
    #endif
            //给归一化过的颜色列表开辟空间
            normedColor_.reserve(normBase_*patchColors_.size());
            //遍历每个颜色，让它们分别减去平均值
            for(std::vector<UserColor>::iterator iter=patchColors_.begin();
                iter!=patchColors_.end();++iter)
            {
                //当前位置的颜色
                cv::Vec3d& currBgr=*iter;
                //依次添加每个通道的差值
                normedColor_.push_back(currBgr[0]-avgColor.Vec::operator[](0));
                normedColor_.push_back(currBgr[1]-avgColor.Vec::operator[](1));
                normedColor_.push_back(currBgr[2]-avgColor.Vec::operator[](2));
                //再另外添加一个平均亮度的差值
                //normedColor_.push_back((tempAvg-avgAvgValue)*3);
            }
        }
        //判断是否需要记录向量
        if(filePtr)
        {
            filePtr->recordVector<VecInterDouble>(*this,this->size(),headStr);
        }
        //把方差设置为1
#ifdef NCC_NORM_1
        setNormAs(1);
#endif
        //处理完成后把标志更改回来
        avgFlag_=true;
    }
};
//颜色面片的列表
typedef std::vector<ColorPatch> ColorPatchList;

//数据的选取器，添加一堆输入流数据，里面的数据到底怎么处理由它里面自己弄
template <typename T>
class DataSelectorBase
{
public:
    //往数据窗口里面添加数据，具体怎么添加由子类负责
    virtual void addData(const T& otherData)=0;

    //获取最终值得选取的数据
    virtual T getData() const=0;

    //获取所有数据的总和
    virtual T getSumData() const =0;
};

//选取中位数的方案
//数据T必须自己实现大于号和小于号
//这里面的逻辑是准备采用插入排序
typedef double MidDataT;//使用泛型编译不过，无奈
class MediumSelector : public DataSelectorBase<MidDataT>
{
protected:
    //从大到小排列
    std::list<MidDataT> datas_;//数据的容器

    std::list<MidDataT>::iterator midIter_;

    //数据的总和
    MidDataT sumValue_=0;
public:
    //实现父类里面添加数据的逻辑规则
    void addData(const MidDataT& otherData) override
    {
        //给数据的总和计数
        sumValue_+=otherData;
        //判断窗口里面是否已经有数据
        if(datas_.size()==0)
        {
            datas_.push_back(otherData);
            midIter_=datas_.begin();
            return;
        }
        //记录数据的标志
        bool saveFlag=false;
        //遍历整个列表，寻找适合添加的位置
        for(std::list<MidDataT>::iterator iter=datas_.begin();
            iter!=datas_.end();++iter)
        {
            //判断传入的数据是否大于当前位置的数据
            if(otherData>*iter)
            {
                //标记找到了适合记录的位置
                saveFlag=true;
                datas_.insert(iter,otherData);
                //结束循环,后面不需要再找了
                break;
            }
        }
        //判断刚才是否找到了可以添加的位置,如果没找到把它添加到最后
        if(!saveFlag) datas_.push_back(otherData);
        //判断刚刚添加的数据是否添加在了中位数的前面
        if(otherData>*midIter_)
        {
            //如果添加完以后是奇数，就不需要移动
            if(datas_.size()%2==0)
                //中位数向前移位
                --midIter_;
        }
        else
        {
            //如果添加完以后是偶数，那就不用动
            if(datas_.size()%2!=0)
                //中位数向后移位
                ++midIter_;
        }
    }

    //获取最后的中位数
    double getData() const override
    {
        //判断数据是否为空
        if(datas_.size()==0) throw ERROR_VISIT_EMPTY_LIST;
        //如果是偶数，返回当前指针和下一个指针的平均值
        if(datas_.size()%2==0)
        {
            //下一个位置的迭代器
            std::list<MidDataT>::iterator nextIter=midIter_;
            ++nextIter;
            //返回两个迭代器数据的平均值
            return (*midIter_+*nextIter)/2;
        }
        return *midIter_;
    }

    //获取数据的总和
    MidDataT getSumData() const override
    {
        return sumValue_;
    }
};

//平均值选取器
template<typename T>
class MeanSelector : public DataSelectorBase<T>
{
private:
    //数据的总和的数据
    T sumData_=0;

    //数据的个数
    uint dataNum_=0;
public:
    //重载父类的接口
    virtual void addData(const T& otherData) override
    {
        //记录数据的累积
        sumData_+=otherData;
        //累积数据
        dataNum_++;
    }

    //重载父类的接口，获取最终计算结果的方式
    virtual T getData() const override
    {
        //返回平均值的计算结果
        return sumData_/dataNum_;
    }

    //重载父类获取数据总和的接口
    virtual T getSumData() const override
    {
        //返回数据的总和
        return sumData_;
    }
};

//用于计算一堆图片之间两两的NCC分数
class EachNcc : public EachOperateInterface<ColorPatch,double>
{
protected:

    //ncc数据的哈希表,用于计算每个图片和其它图片的所有ncc分数的平均值
    //第1个位置的索引表示的仅仅是它在vector里面的标号索引
    //第2个参数是数据选取器
#ifdef USE_NCC_MEDIUM
    std::map<IndexT,MediumSelector> nccMap_;
#else
    std::map<IndexT,MeanSelector<double>> nccMap_;
#endif
    //最大值，相比传入函数的指针，这大概就是做成接口方便的地方吧
    //这个最大值是每个图片和其它ncc图片的总和，但最后返回的是最大图片里面和另一个图片最大的分数
    double maxValue_=0;


    //维护ncc分数
    //把临时的ncc分数添加到所有的ncc分数的哈希表里面
    //并判断是否找到了更大的ncc分数
    void addNccScore(IndexT addId,double tempNcc)
    {
        //记录新传入的数据
        nccMap_[addId].addData(tempNcc);
        //获取目前位置的总和数据
        double currData=nccMap_[addId].getSumData();
        //判断目前位置的选取数据是否超过了最大数据
        if(currData>maxValue_)
        {
            //更新最大值
            maxValue_=currData;
            //更新最大值对应的索引
            bestImgId_=addId;
        }
    }
public:
    //最佳匹配颜色对应的图片索引标号，索引标号不是哈希表里面的键值，而是传入的vector的键值
    IndexT bestImgId_=0;

    //获取最大的分数
    double getMaxScore() const
    {
        //判断最大位置的分数是否存在
        if(nccMap_.count(bestImgId_)==0) throw ERROR_NO_MAP_IDX;
        //获取目标位置的按照规则选取的数据
        return nccMap_.at(bestImgId_).getData();
    }


    //从两个颜色面片中计算结果
    virtual double compute(const ColorPatch &src1, const ColorPatch &src2) override
    {
        //判断两个数据是否都已经设置成平均值为0了
        if(!(src1.avgFlag_&&src2.avgFlag_)) throw ERROR_NOT_AVG_ZERO;
        //计算ncc分数的结果
        double tempNcc=src1.dotProduct(src2);
        //为了防止亮度差异，取绝对值
#ifdef NCC_ABS
        tempNcc=std::abs(tempNcc);
#endif
        //把临时获取到的ncc分数添加到哈希表里面，用于寻找最大的ncc分数
        addNccScore(vecCount,tempNcc);
        addNccScore(vec2Count,tempNcc);
        //返回值直接使用它们两个的数量积
        return tempNcc;
    }
};

//其实就是有多个歧义颜色的情况下应该选择哪个的问题
//本来是应该继承dataSelector接口的，但真正用起来发现情况并不是很合适
class ColorSelector
{
private:
    //聚类算法中使用的点的类型，别人的接口
    typedef kmeans::Point KmPt;

    //待聚类的点的列表，勉强使用别人接口里面使用过的点的数据类型
    std::vector<KmPt> ptList_;

    //点目前的标号
    int pointId_=0;

    //聚类结果最多的一组数据的数量
    int maxDataAccount_=0;

    //聚类结果最多的类对应的颜色数据，每3个为一组，串行连接
    ColorDataSerial centerColor_;

    //数据的维度,0仅仅是默认值
    int dataDim_=0;

    //判断传入的数据维度是否合适
    bool judgeDim(int newDim)
    {
        //如果上一个维度是0,直接记录
        if(dataDim_==0)
        {
            dataDim_=newDim;
            return true;
        }
        //判断新传入的维度是否相同，不同的话要出大问题
        if(newDim!=dataDim_) return false;
        //相同的情况下返回true
        return true;
    }
public:
    //往颜色选取器里面添加待聚类的颜色
    void addData(const ColorPatch &otherData)
    {
        //判断传入的数据维度是否可以接受
        //如果传入的维度不一样，那是要出大问题的
        if(!judgeDim(otherData.patchColors_.size()*3))
            throw ERROR_POINT_DIFF_SIZE;
        //根据传入的颜色数据生成聚类接口的点数据
        ptList_.push_back(KmPt(pointId_,otherData.getData(),
                               otherData.colorAccount()*3));
        //点的标号增加
        ++pointId_;
    }

    //批量添加待聚类的颜色数据
    void addManyData(const ColorPatchList &colorList)
    {
        //遍历所有颜色数据，依次添加
        for(ColorPatchList::const_iterator iter=colorList.begin();
            iter!=colorList.end();++iter)
            addData(*iter);
    }

    //从此结束往点列表里面添加新的数据，并就此开始进行运算
    void endAdding()
    {
        //新建一个聚类算法的对象
        kmeans::KMeans kMeanComputer(CLASS_NUM,KMEAN_ITER_TIME);
        //运行聚类算法,传入已经接收到的点列表
        kMeanComputer.run(ptList_);
        //给颜色数据的容器开辟空间
        centerColor_.resize(dataDim_);
        //记录中心颜色和最大类的聚类数据个数,沙雕QT，能编译通过，但是显示了个报错
        //可能是沙雕openmp造成的
        maxDataAccount_=kMeanComputer.getMaxCluster<ColorDataSerial>(centerColor_);
    }

    //获取最终的颜色数据，其实就是一个颜色的序列
    int getData(ColorPatch& dstColorPatch) const
    {
        //用一个串行的颜色列表来初始化它全部的颜色
        dstColorPatch.initWithColorSerial(centerColor_);
        //返回最多聚类结果对应的数量
        return maxDataAccount_;
    }

    //仅仅获取最大数据的结果
    int getMaxAccount() const
    {
        return maxDataAccount_;
    }
};

//颜色面片的集合，里面包含来自每个面片的数据
class ColorPatchGroup
{
public:
    //中心位置的颜色，刚开始里面是个空的
    ColorPatch centerColor_;

    //将最佳位置的颜色面片设置成所有颜色的平均值
    void setBestAsAvg()
    {
        //最佳位置的颜色面片
        ColorPatch& bestPatch=patchList_[bestImgId_];
        //遍历每个位置的颜色面片
        for(uint imgId=0;imgId<patchList_.size();++imgId)
        {
            //如果是这个面片它自己，那就不用管
            if(imgId==bestImgId_) continue;
            //将当前面片的每个颜色都加到这个最佳的数据上
            bestPatch.addPatch(patchList_[imgId]);
        }
        //除以所有的颜色个数，得到颜色的平均值
        bestPatch.multScale(1.f/patchList_.size());
    }
public:
    //颜色面片的列表
    ColorPatchList patchList_;

    //记录匹配最佳的颜色标号
    //当选用距离主点最近的下视相机作为所谓的平均颜色的时候
    //这个地方存储的是距离主点最近的patch的标号
    int bestImgId_=-1;

    //获取中心颜色的面片
    const ColorPatch& getCenterPatch()
    {
        return centerColor_;
    }

    //使用中心面片颜色对其它面片做处理
    //注意，对中心面片赋值以后才可以做这个操作
    void initOtherPatchWithCenter(FileHandler* filePtr)
    {
        //将中心位置的颜色的方差设置成1方便后面使用
        centerColor_.setAvgZero(nullptr,filePtr,"center");
        //获取中心位置的颜色
        UserColor centerAvg;
        centerColor_.getAvgColor(centerAvg);
        //遍历面片列表里面的每个颜色，依次用现在的平均值来做处理
        for(ColorPatchList::iterator iter=patchList_.begin();
            iter!=patchList_.end();++iter)
        {
#ifdef NCC_CENTER_AVG
            //设置成平均值
            iter->setAvgZero(&centerAvg,filePtr,std::to_string(iter->imgId_));
#else
            //设置成平均值
            iter->setAvgZero(nullptr,filePtr,std::to_string(iter->imgId_));
#endif
        }
    }

    //获取最佳面片的原始图片标号
    IndexT getBestImageViewIndex()
    {
        return patchList_[bestImgId_].imgId_;
    }

    //选取投影最接近主点的图片作为平均颜色面片
    void selectPriNearest()
    {
        //检查是否已经存在了最佳的图片位置
        if(bestImgId_>=0)
        {
            centerColor_.copyColorPatch(patchList_[bestImgId_]);
            return;
        }
        //初始化上下界的距离
        double minDistance=99999;
        //目前最接近的点的索引
        IndexT minColorId=0;
        //判断是否找到过中相机的点
        bool existMid=false;
        //遍历所有的颜色面片
        for(IndexT colorCount=0;colorCount<patchList_.size();++colorCount)
        {
            //当前颜色的面片
            ColorPatch& currPatch=patchList_[colorCount];
            //如果当前的面片不可作为中心，则直接跳过
            if(currPatch.banCenter_) continue;
            //判断当前面片是否为中相机
            if(MID_ID_JUDGE(currPatch.imgId_))
            {
                //更新找到中相机的标志
                existMid=true;
                double tempDis=currPatch.rangeProjDis_;
                //判断是否为更小的距离
                if(tempDis<minDistance)
                {
                    //更新信息
                    minDistance=tempDis;
                    minColorId=colorCount;
                }
            }
        }
        //判断是否存在中相机
        if(existMid)
        {
            //将对应的权值改成余弦
            //patchList_[minColorId].weight_=std::sqrt(1-std::pow(patchList_[minColorId].weight_,2));
            //复制对应位置的面片
            centerColor_.copyColorPatch(patchList_[minColorId]);
            //记录最后被选用的中相机
            bestImgId_=minColorId;
        }
        else
        {
            throw ERROR_NO_MID_FOR_PRI;
        }
    }

    //将中心位置和颜色设置成所有面片的平均颜色
    void setCenterAsAvg()
    {
        //判断颜色的数量是否达标，不达标就不用再做了
        if(patchList_.size()<NCC_MIN_IMG) return;
        //用第1个面片的数据初始化中心位置的面片
        centerColor_.copyColorPatch(patchList_[0]);
        //遍历剩下的颜色面片，依次添加进去
        for(uint colorCount=1;colorCount<patchList_.size();++colorCount)
        {
            //往中心位置的数据里面添加颜色
            centerColor_.addPatch(patchList_[colorCount]);
        }
        //对整体的数据做一个除法，得到最后的平均数据
        centerColor_.multScale(1.f/patchList_.size());
    }

    //获取最适合的颜色值的序列
    const UserColorList& getBestColorList() const
    {
        //使用聚类算法的时候返回的是聚类结果的中心颜色
#ifdef USE_KMEANS
        return centerColor_.patchColors_;
#else //使用NCC算法的时候，返回的是
#ifdef NCC_WRITE_AVG_COLOR//使用平均颜色的时候，中心颜色就是它的平均颜色
        return centerColor_.patchColors_;
#else
        return patchList_[bestImgId_].patchColors_;
#endif
#endif
    }

    //获取聚类结果中最大的一组数据的个数
    double getKmeansMax()
    {
        //新建一个基于聚类算法的颜色选择器
        ColorSelector selector;
        //往里面指添加颜色
        selector.addManyData(patchList_);
        //结束添加颜色，开始进行聚类算法
        selector.endAdding();
        //获取它的中心颜色，并且返回最大的数据结果
        return (double)selector.getData(centerColor_)/patchList_.size();
    }

    //计算目前的颜色信息的ncc结果
    double computeNcc()
    {
        //图片之间两两相互计算NCC
        EachNcc nccComputer;
        std::vector<double> nccMap;//两两计算结果会被存储在这里
        nccComputer.run(patchList_,nccMap);
        //记录匹配最佳的颜色标号
        bestImgId_=nccComputer.bestImgId_;
        //返回计算结果里面的最大值
        return nccComputer.getMaxScore();
    }

    //使用中心位置的颜色来计算和每个颜色的ncc数据
    double computeNccByAvg(FileHandler* filePtr=nullptr)
    {
#ifdef USE_PRI_AS_AVG//如果选取最接近主点的投影点作为面片
        selectPriNearest();
#else
        //将中心颜色设置成所有面片的平均颜色
        setCenterAsAvg();
#endif
        //使用中心面片初始化其它面片
        initOtherPatchWithCenter(filePtr);
        //初始化平均ncc分数
        double avgScore=0;
        //初始化分数的最小值
        double minScore=1;
        //总共的权值叠加
        double weightSum=0;
        //遍历每个颜色,依次计算与平均值的ncc分数
        for(uint colorCount=0;colorCount<patchList_.size();++colorCount)
        {
            //判断是否为匹配最近的点标号
#ifdef USE_PRI_AS_AVG
            if(colorCount==bestImgId_) continue;
#endif
            //当前位置的颜色面片
            ColorPatch& currPatch=patchList_[colorCount];
            //计算当前位置的点与中心位置的平均颜色
            double tempScore=centerColor_.dotProduct(currPatch);
#ifdef NCC_ABS
            tempScore=std::abs(tempScore);
#endif
            //叠加当前位置的ncc分数
            //需要再乘上权值
#ifdef USE_EACH_NCC_WEIGHT
            avgScore+=tempScore*currPatch.weight_;
            //叠加权值
            weightSum+=currPatch.weight_;
#else
            avgScore+=tempScore;
            weightSum++;
#endif
            //判断是否为更小的分数
            if(tempScore<minScore) minScore=tempScore;
            //记录匹配分数
            if(filePtr)
            {
                filePtr->putCosDistance(currPatch.imgId_,tempScore);
            }
        }
        //返回ncc分数的平均值
#ifdef USE_PRI_AS_AVG
        return centerColor_.weight_*(avgScore)/weightSum;
#else
        return (avgScore)/weightSum;
#endif
    }

    //通过欧氏距离计算匹配分数
    double computeEuclidScore(FileHandler* filePtr=nullptr)
    {
#ifdef USE_PRI_AS_AVG//如果选取最接近主点的投影点作为面片
        selectPriNearest();
#else
        //将中心颜色设置成所有面片的平均颜色
        setCenterAsAvg();
#endif
        //填充中心位置的颜色
        centerColor_.fillNormVec();
        //初始化平均距离
        double sumDis=0;
        //遍历每个颜色，依次计算与中心颜色的欧氏距离
        for(IndexT patchCount=0;patchCount<patchList_.size();++patchCount)
        {
            //判断是否为匹配最近的点标号
#ifdef USE_PRI_AS_AVG
            if(patchCount==bestImgId_) continue;
#endif
            //当前位置的面片
            ColorPatch& currPatch=patchList_[patchCount];
            //填充当前位置的颜色
            currPatch.fillNormVec();
            //计算当前颜色与中心颜色的距离
            double tempDis=centerColor_.getEuclidDis(currPatch);
            //将分数除以维度数，最后计算的是每个维度的平均距离
            tempDis/=centerColor_.size();
            //判断是否需要记录距离
            if(filePtr)
            {
                filePtr->putCosDistance(currPatch.imgId_,tempDis);
            }
            //累加平均距离
            sumDis+=tempDis;
        }
        //返回最后的平均距离
#ifdef USE_PRI_AS_AVG
        return (sumDis)/(patchList_.size()-1);
#else
        return (sumDis)/(patchList_.size());
#endif
    }
};


//扩展了z功能的dom平面上的三角面片
//这里用来给整个dom图做初始化的
class DomTrianglePatch : public TriangularPatch2D, public PlaneInterface,public DomTriInterface
{
public:
    //点云数据
    Landmarks* cloudPtsPtr_;
public:

    //重写父类的函数，给接口类提供几个点，方便它计算平面方程
    virtual void getPlanePts(std::vector<Point3D> &dstPtList) const override
    {
        //提前给点开辟空间，这个东西仅仅会往里面加入3个点
        dstPtList.reserve(3);
        //遍历当前三角面片的三个角点
        for(int cornerId=0;cornerId<3;++cornerId)
        {
            //往vector里面添加这三个角点
            dstPtList.push_back(Point3D(cloudPtsPtr_->at(corners_[cornerId]).X));
        }
    }

    //继承三角形里面父类的函数
    virtual void getTrianglePts(std::vector<Point3D> &triPts) const override
    {
        //跟平面上的那个父类其实是同一个接口
        getPlanePts(triPts);
    }

    //构造函数，使用的时候需要调用父类的构造函数，传入三个点云坐标
    //传入点云的数据，做三维数据的时候这里有必要使用
    DomTrianglePatch(IndexT idx0,IndexT idx1,IndexT idx2,Landmarks& cloudPoints)
        : TriangularPatch2D(idx0,idx1,idx2), cloudPtsPtr_(&cloudPoints)
    {
        //需要从外部调用update函数更新平面信息
        //从构造函数里面调用虚函数是危险的
    }

    //仅提供父类信息的构造函数，做这样的一个构造函数是为了复用别的接口
    DomTrianglePatch(IndexT idx0,IndexT idx1,IndexT idx2) : TriangularPatch2D(idx0,idx1,idx2){}

    //给每个三角面片提供一个三维点云的原始数据，然后更新它们的平面信息
    static void updatePatchVector(std::vector<DomTrianglePatch>& patchList,
                                  Landmarks& cloudPts)
    {
        //遍历每个三角面片的信息，分别更新它们的数据
        //这里可以考虑使用多线程
        for(unsigned int patchId=0;patchId<patchList.size();++patchId)
        {
            //当前位置的三角面片
            DomTrianglePatch& thisPatch=patchList[patchId];
            //记录点云信息的指针
            thisPatch.cloudPtsPtr_=&cloudPts;
            //更新平面信息
            thisPatch.updatePlane();
        }
    }

    //把传入的若干点的数据规划在一个合理的范围内
    //这里面已经有若干个点了，但它的z值可能超过这个平面的z值的范围
    //如果有点的z值超过了平面的z值范围，把它规划在一个合理的范围内
    void adjustZInRange(std::vector<Point3D>& ptList) const
    {
        //获取平面内点的方程
        std::vector<Point3D> planePt;
        getPlanePts(planePt);
        //初始化一个范围
        Range1D zRange(planePt[0][2]);
        for(int ptCount=1;ptCount<3;++ptCount) zRange.addData(planePt[ptCount][2]);
        //遍历每个传入的点，让点的z值适应这个范围
        for(Point3DList::iterator iter=ptList.begin();
            iter!=ptList.end();++iter)
        {
            zRange.adjustDataInRange(iter->Point3DBase::operator[](2));
        }
    }
};
//dom问题下使用的三角面片的数据类型
typedef std::vector<DomTrianglePatch> DomPatchList;

//opencv里面使用的颜色数据，一个临时的定义
typedef uchar* CvColorPtr;
//opencv的颜色指针列表，用于在另一个地方给dom图的结果写入颜色
typedef std::vector<uchar*> CvColorPtrList;

//dom的每个单元的优化器
//当写下这个class之后，我也不知道我该干嘛了
//这一段儿算法好难
//父类里面的二维点表示的是当前位置的x,y信息,真实值，不是dom图的像素值
class DomZOptimizer : public TempPtMid
{
protected:
    double saveZ_;
public:
    //整个class里面这是最关键的属性，优化就是在优化这个东西
    //这个指针的数值可以修改，其它情况下建议不要修改这里的数值
    double* zPtr_;
    //是否为带有先验信息的z值
    bool isPrior_=false;
    //该点与dom图链接的ncc分数
    //这个东西是在类外被手动赋值的，domInfo那里
    double* nccPtr_=nullptr;
    //当前位置的(x,y)坐标在每个图片上的投影直线
    ProjLineList viewLines_;
    //待写入的颜色，与当前位置相对应,最好不要自己动它
    CvColorPtr colorPtr_;
    //表示z范围的全局变量
    const std::array<double,2>* zRange_=nullptr;
    //该点作为参考点时可使用的次数
    RefType* refTimePtr_=nullptr;
    //像素的方向数据
    PixelDir* pixDirPtr_=nullptr;
    //DOM像素在整体的DOM单元中的顺序标号
    //添加这个属性是为了方便做外部直线的存取
    Index domUnitIdLabal_=0;

    //减少当前DOM像素的可被优化次数
    void minusOptTime()
    {
        //如果已经是0了直接返回
        if(*refTimePtr_==0) return;
        --(*refTimePtr_);
    }

    //根据已知的信息获取当前Z优化器里面存储的一个三维点
    template<typename PointT>
    void getPoint3D(PointT& dstPoint)
    {
        //分别记录每个位置的值
        dstPoint[0]=xCoord();
        dstPoint[1]=yCoord();
        dstPoint[2]=*zPtr_;
    }

    //获取它在某个图片上的投影点
    //这里的viewId不是scene里面的viewId，仅仅是属性viewLines_的索引标号
    void getProjLocal(IndexT viewId,double& dstX,double& dstY)
    {
        //判断表示z范围的全局变量是否初始化过
        if(zRange_==nullptr) throw ERROR_EMPTY_ZRANGE;
        //判断目前的z坐标是否超出了范围 如果直接本身不可用也不行
        //直线不可用这个情况是20211203添加的，是一种新的处理遮挡问题的方案
        if(*zPtr_<zRange_->operator[](0)||*zPtr_>zRange_->operator[](1)||
                viewLines_[viewId].isUsable_==false)
        {
            //记录一个非法的坐标，让它变得不能用
            dstX=-1;
            dstY=-1;
            return;
        }
        //根据目前的z值获取x和y坐标
        viewLines_[viewId].getXy(*zPtr_,dstX,dstY);
    }

    //获取某个图片的投影界限距离
    //对于某一个图片，将上界投影到图片上获取一个点
    //再把Z的下界投影到这个图片上获取一个点
    //计算两个点的距离并返回
    //注意这里的viewid不是图片的索引标号，仅仅是直线里面的索引标号
    double getProjRangeDis(IndexT viewId)
    {
        //获取要处理的投影直线
        ProjectLine& currLine=viewLines_[viewId];
        //获取下界
        TempPtBase projDown;
        currLine.getXy(zRange_->operator[](0),projDown[0],projDown[1]);
        //获取上界
        TempPtBase projUp;
        currLine.getXy(zRange_->operator[](1),projUp[0],projUp[1]);
        //计算上下界的距离
        return (projDown-projUp).norm();
    }

    //获取颜色的亮度,调用颜色接口里面的操作
    double colorBrightness()
    {
        //调用颜色接口里面的静态函数
        return ColorInterface<uchar>::colorBrightness(colorPtr_);
    }

    //设置z值，同时还需要判断一下这里是不是一个先验的值
    //如果是先验的值，那就不理会传入的数据，也不报错
    //如果里面的z值不能被修改，返回false
    bool setZ(double zValue)
    {
        //判断这是否为一个先验的点
#ifndef PRIOR_REFINE
        if(isPrior_) return false;
#endif
#ifndef OPT_REFINE//使用z平面的状态下，找到了就不动了
        if(*nccPtr_>=NCC_THRE) return false;
#endif
#ifndef REFINE_REFINE_ZERO
        if(*refTimePtr_==0) return false;
#endif
        //正常情况下就直接赋值然后返回
        *zPtr_=zValue;
        return true;
    }

    //这个优化器刚刚参与了一次面片优化
    //这里传入的最大分数是刚刚那个面片的整体最佳分数，不一定适用于这个优化器
    //传入的z值是这个优化器历史上的最佳值，目前这个优化器里面的z值已经在优化面片的过程中被污染了
    bool selectZ(double zValue,double maxScore)
    {
        //判断传入的最佳分数是否大于自己的历史最佳分数
        //如果大于自己的历史最佳分数，那就不恢复自己以前的那个历史最佳z值了
        if(maxScore>=*nccPtr_) return false;
        else
        {
            //恢复到自己以前的历史最佳z值
            return setZ(zValue);
        }
    }

    //往数据里面写入颜色
    void writeColor(uchar bColor,uchar gColor,uchar rColor,double nccScore)
    {
        //检查黑色数值
        auto colorSum=colorPtr_[0]+colorPtr_[1]+colorPtr_[2];
        //判断传入的分数是否更大
        if(colorSum!=0&&
                ((*nccPtr_)!=PRIOR_NCC&&nccScore<*nccPtr_)) return;
        colorPtr_[0]=bColor;
        colorPtr_[1]=gColor;
        colorPtr_[2]=rColor;
    }

    //往数据里面记录法向量
    void writeNorm(const PixelDir& otherDir,double nccScore)
    {
        //判断传入的分数是否更大
        if(nccScore<*nccPtr_) return;
        //记录法向量
        pixDirPtr_->operator=(otherDir);
    }

    //对z坐标做叠加
    bool addZ(double zForward)
    {
        return setZ(*zPtr_+zForward);
    }

    //通过传入两个点坐标来实现的初始化操作
    DomZOptimizer(double xValue=0,double yValue=0,double* zPtr=nullptr)
        : TempPtMid(xValue,yValue)
    {
        //记录指针
        zPtr_=zPtr;
    }

    //转换成三维点数据
    void toPoint3D(Point3D& dstPoint) const
    {
        dstPoint[0]=xCoord();
        dstPoint[1]=yCoord();
        dstPoint[2]=*zPtr_;
    }

    //记录目前的z值,用于缓存以前的数据
    //如果传入的分数小于0,那就是一个强行记录的指令
    void saveData(double nccScore)
    {
        if(nccScore<0)
        {
            //强行记录z值
            saveZ_=*zPtr_;
            return;
        }
        //判断ncc分数是否大于这里已经保存的分数
        if(nccScore>*nccPtr_)
        {
            //记录z值
            saveZ_=*zPtr_;
            //更新ncc分数,这里很奇怪，我明明记得我单独写过一个位置让它去更新分数的
            //但我找不到我当时写过的那个更新分数的地方了
            //用QT的搜索功能又找了一下，确实没有了
            *nccPtr_=nccScore;
        }
        else//否则把目前的z值恢复到以前的z值
        {
            loadData();
        }
    }

    //载入z数值，用于恢复上次的计算结果
    void loadData()
    {
        *zPtr_=saveZ_;
    }

    //获取此点的3D信息
    Point3DBase get3DPoint() const
    {
        //新建一个三维点
        Point3D tempPt;
        tempPt<<xCoord(),yCoord(),*zPtr_;
        //返回三维点
        return tempPt;
    }
};
//z优化器的列表
typedef std::vector<DomZOptimizer> OptimizerVector;

//z值的直方图对儿
//表示z的value和它的代表标号
class ZHistogramPair
{
public:
    double value_;//z值
    IndexT repId_;//在原始的面片里面的代表标号
    double brightness_=-1;//当前位置的最高亮度

    //因为需要用到vector的resize，所以必须保留一个无参数的构造函数
    ZHistogramPair(){}

    //构造函数
    ZHistogramPair(double zValue,IndexT repId,double bright=-1)
    {
        value_=zValue;
        repId_=repId;
        //记录当前位置的最高亮度
        brightness_=bright;
    }

    //类似于构造函数的初始化各种数据
    void initValues(double zValue,IndexT repId,double bright=-1)
    {
        value_=zValue;
        repId_=repId;
        //记录当前位置的最高亮度
        brightness_=bright;
    }
};

//dom单元里面的一个矩形块
//继承平面的相关特性
class DomRectPatch : public PlaneInterface
{
protected:
    Dir3D saveNorm_;//缓存的平面法向量的信息

    //面片内可能出现的z的高度值
    std::vector<ZHistogramPair> zSteps_;
public:
    //不同的z优化器,每个优化器代表一个位置
    //这个class本身并不在乎一共有多少个位置
    //这里面的点是按照从左到右从上到下的顺序存储的
    //想要使用中心位置的点，那就用size()/2来找
    OptimizerVector optiVec_;

    //这里使用的是不遮挡并且投影范围最近的一个点
    //保存的是最合适的点在optiVec_里面的索引位置
    IndexT bestProjectView_=0;

    //当前面片的主亮度，面片的亮度基本一致是一个基本假设，不同的亮度不参与同一面片的计算
    //主亮度小于0的时候则主亮度不做参考
    double mainBright_=-1;

    //获取最佳射影直线的最佳投影位置的图片标号
    IndexT getBestProjectViewSrcImageId()
    {
        return getCenterZOpt().viewLines_[bestProjectView_].imgId_;
    }

    //判断是否存在没被优化过的点
    bool haveUnrefine(double nccThre=NCC_THRE) const
    {
        //遍历优化器
        for(OptimizerVector::const_iterator iter=optiVec_.begin();
            iter!=optiVec_.end();++iter)
        {
            //判断是否没优化过
            if(*(iter->nccPtr_)<=nccThre)
            {
#ifdef USE_MAX_OPT_TIME
                if(*(iter->refTimePtr_)>0)
#endif
                    return true;
            }
        }
        //没找到就返回false
        return false;
    }

    //区域里面的代表点的索引
    IndexT repOptIdx_=0;

    //寻找所有优化点里面的先验点，找到即返回
    IndexT findPriorPoint() const
    {
        //遍历所有的优化数据
        for(IndexT optId=0;optId<optiVec_.size();++optId)
        {
            //判断是否为先验点
            if(optiVec_[optId].isPrior_) return optId;
        }
        //如果都没有找到就返回整个向量的大小
        return optiVec_.size();
    }

    //初始化当前矩形面片里面可能出现的z阶层
    //并返回总共的z阶层的个数
    //最后一个参数是用德劳内三角化的z值做辅助判断的时候使用的
    //minZStep是所允许的最小的z值的阶梯，为了保证较高位置优先被重建所设置的参数
    uint initZSteps(const RangeType& zRange,double stepRate,bool useCenterZ=false,double nccThre=NCC_THRE
            ,int minZStep=0)
    {
        //分别记录都有哪个阶层的数据出现过
        //第1个参数表示z的哪个阶层，第2个参数表示这个阶层的代表标号
        std::map<int,IndexT> zRecorder;
        //每个阶层的长度
        double stepHeight=(zRange[1]-zRange[0])*stepRate;
        //遍历每个z优化器
        for(uint optCount=0;optCount<optiVec_.size();++optCount)
        {
            //当前的z优化器
            const DomZOptimizer& iter=optiVec_[optCount];
            //判断是否为优化过的点
            if(iter.isPrior_ || *(iter.nccPtr_)>=nccThre)
            {
                //判断是否还有可用的次数
#ifdef USE_REFINE_TIME //限制每个点的可参考次数的时候需要禁用不可用的点，否则一个点是一直可被参考的
                if(*(iter.refTimePtr_)==0) continue;
#endif
                //判断当前的z值所属的阶层
                int currStep=(*iter.zPtr_-zRange[0])/stepHeight;
                //判断是否超过了所允许的最小的高度阶层
                if(currStep<minZStep) continue;
                //判断这个阶层是否存在
                if(zRecorder.count(currStep)==0)
                {
                    //添加这个位置的数值
                    zRecorder[currStep]=optCount;
                }
            }
        }
        //初始化z阶层的记录
        zSteps_.resize(zRecorder.size());
        //记录目前写入的位置
        IndexT currId=zRecorder.size()-1;
        //初始化目前的位置，保证它是从大到小排列的
        //遍历z阶层的记录值
        for(auto eachStep : zRecorder)
        {
            //当前阶层的z数值
            double zValue=*(optiVec_[eachStep.second].zPtr_);
            //添加z阶层的数值
            zSteps_[currId].initValues(zValue,eachStep.second
#ifdef USE_BRIGHT_CONSTRAIN
                                       ,optiVec_[eachStep.second].colorBrightness());
#else
                                  );
#endif
            //更新接下来要写入的位置
            --currId;
        }
        //判断是否存在z的阶层数
        if(zSteps_.size()==0)
        {
            //判断是否需要借助中间位置的点
            if(useCenterZ)
            {
                //计算中间位置的索引
                IndexT centerId=optiVec_.size()/2;
                //把中间位置的z值添加到z阶层表里面
                zSteps_.push_back(ZHistogramPair(*optiVec_[centerId].zPtr_,centerId));
            }
        }
        //返回总共的z阶层个数
        return zSteps_.size();
    }

    //寻找所有优化点里面已经被优化过的点
    IndexT findRefinePt(double nccThre=NCC_THRE) const
    {
        //遍历所有数据
        for(IndexT optId=0;optId<optiVec_.size();++optId)
        {
            //判断是否为已经优化过的点
            if(*(optiVec_[optId].nccPtr_)>nccThre) return optId;
        }
        //如果妹找到，返回最后一位
        return optiVec_.size();
    }

    //更新优化器里面的代表点
    //这个函数是需要在外部调用的，最好更新了优化器信息后立刻就调用这个函数
    void updateRepPoint(double nccThre=NCC_THRE)
    {
        //判断点的个数是否超过了要求的最少的点
        if(optiVec_.size()<MIN_MASK_PT) throw ERROR_MASK_PT_TOO_LESS;
        //寻找点列表里面有没有先验的点
        repOptIdx_=findPriorPoint();
        if(repOptIdx_<optiVec_.size()) return;
        //寻找列表里面有没有已经优化过的点
        repOptIdx_=findRefinePt(nccThre);
        if(repOptIdx_<optiVec_.size()) return;
        //把中心位置的点作为代表点
        repOptIdx_=optiVec_.size()/2;
    }

    //获取平面高度的阶层的个数
    uint getZStepNum() const
    {
        return zSteps_.size();
    }

    //更新面片里面的z值，让它们等于特定阶层的高度
    //holdZValue为true的时候仅选择中心像素，不对其它像素做修改
    bool updateAsZStep(IndexT zStepId,bool holdZValue=false)
    {
        //判断阶层的高度是否超过了索引的范围
        if(zStepId>=getZStepNum()) throw ERROR_VEC_OUT_RANGE;
        //更新代表点的标号
        repOptIdx_=zSteps_[zStepId].repId_;
#ifdef USE_REFINE_TIME
        //判断这个位置的z是否可用
        if((*optiVec_[repOptIdx_].refTimePtr_)==0) return false;
        //根据代表点的标号，将它的可用次数递减
        --(*optiVec_[repOptIdx_].refTimePtr_);
#endif
        //计算该面片的法向量
        Dir3D tempDir;
#ifdef USE_NORM_PROG
        optiVec_[repOptIdx_].pixDirPtr_->makeDir(tempDir);
#else
        tempDir<<0,0,1;
#endif
        //获取中心点
        Point3D centerPt;
        getCenterPt(centerPt);
        //使用法向量和中心点生成新的平面
        PlaneInterface::makePlaneByNormPt(planeEqu_,tempDir,centerPt);
        if(holdZValue) return true;
        //使用已经更新过的点和法向量
        updateOptimizer();
#ifdef USE_BRIGHT_CONSTRAIN
        //设置当前位置的主亮度
        //最高位置的点是不需要做这样的参考的，因为不存在遮挡
        if(zStepId+1<getZStepNum())
            mainBright_=zSteps_[zStepId].brightness_;
#endif
        //正常情况下返回可用
        return true;
    }

    //实现父类里面获取点列表的接口
    void getPlanePts(std::vector<Point3D> &dstPtList) const override
    {
        //提前开辟空间
        dstPtList.reserve(optiVec_.size());
        //遍历每个z优化器
        for(OptimizerVector::const_iterator iter=optiVec_.begin();
            iter!=optiVec_.end();++iter)
        {
            //记录当前位置的点
            dstPtList.push_back(Point3D(iter->xCoord(),iter->yCoord(),*iter->zPtr_));
        }
    }

    //获取中心位置的Z优化器
    DomZOptimizer& getCenterZOpt()
    {
        return  optiVec_[repOptIdx_];
    }

    //获取中心点
    void getCenterPt(Point3D &dstPt) override
    {
        //获取中心位置的点
        optiVec_[repOptIdx_].toPoint3D(dstPt);
    }

    //获取中心位置点的z值
    double getCenterZ() const
    {
        return *(optiVec_[repOptIdx_].zPtr_);
    }

    //保存面片里面当前的数据信息，万一后面的优化结果不行了
    //需要在这里把优化数据重新取出来,做优化的时候需要使用这样的操作
    void saveData(double nccScore=-2.f)
    {
        //记录父类平面的法向量信息
        getNormal(saveNorm_);
        //遍历每个位置的z优化器，记录它们此时的z值
        for(OptimizerVector::iterator iter=optiVec_.begin();
            iter!=optiVec_.end();++iter)
                iter->saveData(nccScore);
    }

    //载入之前记录过的三角面片的数据,这是属于优化失败了采取的操作
    void loadData()
    {
        //获取法向量和中心点
        Point3D centerPt;
        getCenterPt(centerPt);
        //用上次迭代的法向量和中心点更新平面
        updatePlane(saveNorm_,centerPt);
        //遍历每个位置的z优化器，恢复它们上次的z值
        for(OptimizerVector::iterator iter=optiVec_.begin();
            iter!=optiVec_.end();++iter)
            iter->loadData();
    }

    //根据新的平面方程，更新z优化器
    void updateOptimizer()
    {
        //遍历每个z优化器
        for(OptimizerVector::iterator iter=optiVec_.begin();
            iter!=optiVec_.end();++iter)
        {
            //从平面信息里面计算z值并对优化器里面的z赋值
            iter->setZ(computeZ(iter->xCoord(),iter->yCoord()));
        }
    }

    //根据传入位置的点列表依次把每个位置的点都记录一下
    //需要自行和z优化器里面每个位置是一一对应的
    //这里面传入的是它以前的点坐标，如果有更好的分数，就不需要再拿这里的数据了
    //这里传入的分数是刚刚约束出来的分数，每个优化器里面还存储着自己历史上的最佳分数
    //现在给每个优化器一个机会，选择是否恢复到自己以前的历史最佳分数
    //如果发现传入的分数还不如自己的历史最佳分数，就恢复成自己的历史最佳分数
    //把maxScore的默认值设置成0,这样如果不写这个数值，那么所有的数据都会选择恢复到自己以前的最佳数据
    void loadPointList(const Point3DList& ptList, const PlaneT* planePtr=nullptr,double maxScore=0)
    {
        //判断传入的点的个数和优化器的个数是否相同
        if(ptList.size()!=optiVec_.size()) throw ERROR_POINT_DIFF_SIZE;
        //记录每个点的坐标
        //这里可以考虑使用多线程
        for(unsigned int ptId=0;ptId<ptList.size();++ptId)
            optiVec_[ptId].selectZ(ptList[ptId].getZ(),maxScore);
        //判断是否需要更新平面
        if(planePtr!=nullptr) planeEqu_=*planePtr;
    }


    //把最终选取出来的颜色序列记录到每个对应位置的z上
    //传入ncc分数是用于判断是否有必要保存
    //normVec表示的是法向量
    void drawColorPatch(const UserColorList& colorList,const Point3DBase& normVec,double nccScore)
    {
        //判断颜色列表和内部存储的颜色个数是否一致
        if(colorList.size()!=optiVec_.size()) throw ERROR_POINT_DIFF_SIZE;
        //从法向量里面初始化一个适用于存储格式的法向量
        PixelDir tempDir;
        tempDir.fromDir(normVec);
        //可以考虑使用多线程
        for(unsigned int ptCount=0;ptCount<optiVec_.size();++ptCount)
        {
            //当前位置的颜色
           const cv::Vec<UserColorType,3>& thisColor=colorList[ptCount];
            //在z优化器里面写入对应的颜色
            optiVec_[ptCount].writeColor(thisColor[0],thisColor[1],thisColor[2],nccScore);
            //在z优化器里面记录对应的法向量
#ifdef USE_NORM_PROG
            optiVec_[ptCount].writeNorm(tempDir,nccScore);
#endif
        }
    }

    //更新主点的z值,这里的更新值是增量而不是直接的变量
    void updatePriZ(double zForward)
    {
        //遍历迭代器，将每个位置都做一个平移
        for(OptimizerVector::iterator iter=optiVec_.begin();
            iter!=optiVec_.end();++iter)
            iter->addZ(zForward);
    }

    //设置平面里面每个位置的z值，把它们设置成一个水平的平面
    void setAllZ(double zValue)
    {
        //遍历迭代器，设置每个位置的z值
        for(OptimizerVector::iterator iter=optiVec_.begin();
            iter!=optiVec_.end();++iter)
            iter->setZ(zValue);
    }

    //利用已经读取下来的直线数据列表初始化当前patch里面的Z优化器里面的直线数据
    //一个图片必须能够看到所有的Z优化器才可以，所以统计信息里面必须和优化器的个数一致
    void initProjectLinesByListGroup(ProjectLineListGroup& linesForEachOptimizer,
                                     CountMap& idCountMap)
    {
        if(linesForEachOptimizer.size()!=optiVec_.size())
        {
            throw ERROR_INIT_PROJECT_LINES;
        }
        //遍历每个Z优化器
        for(IndexT idZOptimizer=0;idZOptimizer<optiVec_.size();++idZOptimizer)
        {
            //当前位置的Z优化器
            DomZOptimizer& currOptimizer=optiVec_[idZOptimizer];
            //当前位置的直线列表
            ProjLineList& currProjectLineVec=linesForEachOptimizer[idZOptimizer];
            //开辟空间
            currOptimizer.viewLines_.reserve(currProjectLineVec.size());
            //遍历直线列表里面的每个直线
            for(IndexT idLine=0;idLine<currProjectLineVec.size();++idLine)
            {
                //当前位置的直线
                ProjectLine& currLine=currProjectLineVec[idLine];
                //判断数量是否足够
                if(idCountMap[currLine.imgId_]==optiVec_.size())
                {
                    //把直线加入到列表中
                    currOptimizer.viewLines_.push_back(currLine);
                }
            }
        }
    }

    //利用外部数据初始化每个DOM像素的投影直线
    void initProjectLinesByExternalData()
    {
        //能看到当前面片的相机数据统计
        CountMap idCountMap;
        //新建临时的向量，用于临时存储每条线
        ProjectLineListGroup linesForEachOptimizer;
        //给每个DOM的Z优化器提前开辟相应的位置
        linesForEachOptimizer.resize(optiVec_.size());
        //当前面片的所有Z优化器
        for(IndexT idZOptimizer=0;idZOptimizer<optiVec_.size();++idZOptimizer)
        {
            //当前的Z优化器对应的投影直线的列表
            ProjLineList& currProjectLineList=linesForEachOptimizer[idZOptimizer];
            //当前的Z优化器
            DomZOptimizer& currOptimizer=optiVec_[idZOptimizer];
            //获取一个Z优化器里面全部的直线
            ProjectLine::loadProjectLines(currProjectLineList,currOptimizer.domUnitIdLabal_);
            //对所有的直线信息做统计
            for(IndexT idLine=0;idLine<currProjectLineList.size();++idLine)
            {
                //当前位置的直线
                ProjectLine& currLine=currProjectLineList[idLine];
                if(idCountMap.count(currLine.imgId_)==0)
                {
                    idCountMap[currLine.imgId_]=1;
                }
                else
                {
                    ++idCountMap[currLine.imgId_];
                }
            }
        }
        //用直线列表和统计信息初始化所有的Z优化器里面的直线数据
        initProjectLinesByListGroup(linesForEachOptimizer,idCountMap);
    }
};

//基于图片的三角面片,表示的是一个图片上三角面片的坐标
//继承了一些不需要的属性，毕竟那个基于DOM点云的三角面片的类是先被设计出来的
//没想到有一天居然还会用到别的形式
//注意imgId_需要特地被初始化一下
class ImageTriPatch : public TriangularPatch2D
{
public:
    //当前的图片三角面片所属的图片标号
    IndexT imgId_;//需要专门被初始化，由于应用场景，不得不这样

public:
    //构造函数，调用父类
    ImageTriPatch(IndexT pt1=0,IndexT pt2=0,IndexT pt3=0) : TriangularPatch2D(pt1,pt2,pt3)
    {
    }

    //获取三角面片的最长边,比如它返回了一个0,那么代表的就是0~1两个点形成的边
    IndexT getLongestSide(const Landmarks &cloudPoints)
    {
        //初始化最长的边长度
        double maxLen=0;
        //最长边对应的标号
        IndexT maxLineId=4;
        //遍历三角面片的每个边
        for(int cornerId=0;cornerId<3;++cornerId)
        {
            //当前位置的直线
            DomImgLine tempLine;
            tempLine.initImgLine(cloudPoints,corners_[cornerId],corners_[(cornerId+1)%3],imgId_);
            //判断直线的长度是否够大
            if(tempLine.getLength()>maxLen)
            {
                maxLen=tempLine.getLength();
                maxLineId=cornerId;
            }
        }
        //返回最长边对应的直线标号
        return maxLineId;
    }

    //获取最终的主迭代信息
    virtual void getFirstPatchIterInfo(PatchMainIterInfo &dstInfo, const Landmarks &cloudPoints,
            const DomImgLine &maxImgLine, IndexT maxLineIdx) override
    {
        //在主迭代信息的父类里面写上图片的直线
        dstInfo.imgLines_.push_back(maxImgLine);
        //确定迭代直线的对角点,初始化它的点
        dstInfo.otherCorner_.initPoint(cloudPoints,corners_[(maxLineIdx+2)%3]);
        //直线的另一端的点标号
        IndexT nextLinePtIdx=(maxLineIdx+1)%3;
        //记录它点云点中的主直线
        dstInfo.domLine_.initLine(cloudPoints,corners_[maxLineIdx],corners_[nextLinePtIdx]);
    }

    //获取主迭代信息
    virtual void getFirstPatchIterInfo(PatchMainIterInfo &dstInfo, const Landmarks &cloudPoints) override
    {
        //获取三角面片的最长边对应的直线标号
        IndexT longestLineId=getLongestSide(cloudPoints);
        //生成最长的直线
        DomImgLine maxLine;
        maxLine.initImgLine(cloudPoints,corners_[longestLineId],corners_[(longestLineId+1)%3],imgId_);
        //调用上层的生成主迭代信息的接口，获取最终的主迭代信息
        getFirstPatchIterInfo(dstInfo,cloudPoints,maxLine,longestLineId);
    }

    //获取cv形式的坐标列表
    void getCvPtList(const Landmarks &cloudPoints, std::vector<cv::Point2i> &dstPtList)
    {
        //开辟空间
        dstPtList.reserve(3);
        //取出每个点
        for(int ptCount=0;ptCount<3;++ptCount)
        {
            //当前位置的点坐标
            const TempPtBase& obvPt=cloudPoints.at(corners_[ptCount]).obs.at(imgId_).x;
            //记录坐标
            dstPtList.push_back(cv::Point2i(obvPt[0],obvPt[1]));
        }
    }

    //把当前的三角面片画在图上
    void drawPatch(cv::Mat &img,const Landmarks &cloudPoints)
    {
        //获取opencv形式的坐标列表
        std::vector<cv::Point2i> cvPtList;
        getCvPtList(cloudPoints,cvPtList);
        //在图上画线
        cv::polylines(img,cvPtList,true,cv::Scalar(0,0,255));
    }
};


//德劳内的三角面片序列
typedef std::vector<TriangularPatch2D> DelaunayPatchList;
//适用于图片的三角面片列表
typedef std::vector<ImageTriPatch> ImgPatchList;

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_LANDMARK_HPP
