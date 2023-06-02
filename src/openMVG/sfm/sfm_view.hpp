// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_VIEW_HPP
#define OPENMVG_SFM_SFM_VIEW_HPP

#include <string>

#include "openMVG/types.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include<opencv2/opencv.hpp>
#include"switch.hpp"

namespace openMVG {
namespace sfm {

//计算传入的向量的反对称矩阵,全局函数不带static会报错，openmvg不知道这是个什么编译环境
static Eigen::Matrix3d makeSkew(Eigen::Vector3d fromVec)
{
    Eigen::Matrix3d finalMat;
    finalMat<<0,-fromVec(2),fromVec(1),
            fromVec(2),0,-fromVec(0),
            -fromVec(1),fromVec(0),0;
    return finalMat;
}

//基础矩阵的数据类型
typedef Eigen::Matrix3d FundData;
//基础矩阵,之所以要做一个类，是为了方便以后对传入的点做投影
class FundamentalMatrix
{
public:
    //基础矩阵的内容,自己负责赋值,本class不管
    FundData fundData_;

    //根据传入的外参，计算单应性矩阵
    void computeF(const Eigen::Matrix3d &relaR,//相对旋转
                  const Eigen::Vector3d &relaC, //相对光心
                  const Eigen::Matrix3d &intr1, //第1个相机的内参
                  const Eigen::Matrix3d &intr2 //第2个相机的内参
                  )
    {
        //根据计算机视觉中的多视力几何163页
        fundData_=intr2.inverse().transpose()*relaR*intr1.transpose()*
                makeSkew(intr1*relaC);
    }
};

//存储基础矩阵的哈希表
typedef std::map<IndexT,FundamentalMatrix> Fundamentals;
//存储点标号的列表
typedef std::vector<IndexT> PtIdList;
//表示图片上4个角的范围,表示左右上下,上小下大，左小右大
typedef std::array<int,4> ImgCornerRange;
//图片颜色指针的列表
typedef std::array<const cv::Vec3b*,4> ColorPtrArray;
//z值图使用的数据类型
typedef uchar ZType;


/// A view define an image by a string and unique indexes for the view, the camera intrinsic & the pose
struct View
{

    //获取view总共的载入了的图片个数
    static unsigned& getLoadedImageCount()
    {
        //添加它对应的静态变量
        static unsigned imageCount;
        return imageCount;
    }

  // image path on disk
  std::string s_Img_path;

  //读取后的图片，感觉又要费好多内存
  cv::Mat cvImg_;

  //图片的高度图，与上面的图片一一对应，分别对应于它的高度
  cv::Mat imgZMap_;

  //2022-3-11 生成原图的mesh用的一个小变量，平常没用
  bool usedFlag=false;

  //一个全零的颜色，纯属就是为了给当初不好的设计擦屁股
  //这个写法实在是太low了，害，能用就行，真是难为后来人了
  cv::Vec3b zeroColor_;

  // Id of the view
  IndexT id_view;

  // Index of intrinsics and the pose
  IndexT id_intrinsic, id_pose;

  // image size
  IndexT ui_width, ui_height;

  //当前的图片和其它的图片的基础矩阵的列表
  Fundamentals fundmentals_;

  //获取图片的分辨率
  double imgResolution() const
  {
      return ui_width*ui_height;
  }

  //根据传入的缩放比例将图片缩放
  void resizeScale(double scale)
  {
      //对于z约束的情况，把z图更改一下大小
#ifdef USE_IMG_Z_CONSTRAIN
      cv::Mat tempSrcImg=imgZMap_.clone();
      cv::resize(tempSrcImg,imgZMap_,cv::Size(scale*ui_width,scale*ui_height));
#endif
  }

  //载入图片并返回图片的引用
  cv::Mat &loadImg(std::string rootPath="")
  {
      //判断是否需要载入
      if(cvImg_.empty())
      {
          cvImg_=cv::imread(rootPath+"/"+s_Img_path);
          //判断是否仍然没有图片
          if(cvImg_.empty())
          {
              throw ERROR_IMREAD_FAIL;
          }
          //如果需要使用图片高度注册的方式，那么这里需要提前注册一个每个像素的高度大小
#ifdef USE_IMG_Z_CONSTRAIN
          imgZMap_.create(cvImg_.rows,cvImg_.cols,CV_8UC1);
          imgZMap_.setTo(0);
#endif
          //记录成功添加过的图片个数
          getLoadedImageCount()++;
      }
      //返回图片
      return cvImg_;
  }

  //判断点坐标是否在图片的范围内
  bool inImgRange(double pixelX,double pixelY) const
  {
      int xLoc=std::round(pixelX);
      int yLoc=std::round(pixelY);
        return xLoc>=0 && yLoc>= 0 && xLoc<ui_width && yLoc<ui_height;
  }

  //判断传入的点坐标是否在Z值图的范围内
  bool inZMapRange(double pixelX,double pixelY) const
  {
      int xLoc=std::round(pixelX);
      int yLoc=std::round(pixelY);
        return xLoc>=0 && yLoc>= 0 && xLoc<imgZMap_.cols && yLoc<imgZMap_.rows;
  }

  //把传入的点转换到Z值图的坐标系下，它们可以差一个缩放比例
  cv::Point2i transToZCoord(int xValue,int yValue)
  {
      //初始化待返回的值
      cv::Point2i retValue;
      retValue.x=((double)xValue/cvImg_.cols)*imgZMap_.cols;
      retValue.y=((double)yValue/cvImg_.rows)*imgZMap_.rows;
      //返回结果
      return retValue;
  }

  //访问一个特定位置的颜色，返回它的引用
  //调用之前必须保证图片已经被初始化
  //最后一个参数是z的参考高度，如果传入小于0的数字表示不使用这个参数
  cv::Vec3b& colorAt(double pixelX,double pixelY,int zLevel=-1)
  {
      //对Z值做变换
      zLevel=Z_LEVEL_TRANS(zLevel);
      //读取图片
      cv::Mat &img=loadImg();
      //如果不在范围内，就返回这个颜色在
      if(!inImgRange(pixelX,pixelY))
      {
          throw ERROR_IMG_OUT_RANGE;
      }
      int xRd=std::round(pixelX);
      int yRd=std::round(pixelY);
#ifdef USE_IMG_Z_CONSTRAIN
      //如果是更小的z直接返回一个全零内容
      //如果传入的z是一个小于零的数字，说明不用比较
      if(zLevel>=0&&zLevel<imgZMap_.at<ZType>(transToZCoord(xRd,yRd)))
      {
            return zeroColor_;
      }
#endif
      //返回图片的颜色
      return img.at<cv::Vec3b>(yRd,xRd);
  }

  //在指定的投影位置注册一个高度
  void registerZMap(double pixelX,double pixelY,int zLevel)
  {
      //对Z值做变换
      zLevel=Z_LEVEL_TRANS(zLevel);
      //判断z高度是否为空
      if(imgZMap_.empty()) throw ERROR_IMREAD_FAIL;
      //判断图片是否在范围内
      if(!inImgRange(pixelX,pixelY)) throw  ERROR_IMG_OUT_RANGE;
      //获取指定位置的z值
      ZType& dstZValue=imgZMap_.at<ZType>(
                  transToZCoord(std::round(pixelX),std::round(pixelY)));
      //判断是否传入的是更大的z值
      if(zLevel>dstZValue)
      {
          dstZValue=zLevel;
      }
  }

  //在图上记录两个点之间的z值
  //如果遇到了结束条件，就return true
  bool registerZMap(const std::vector<cv::Point2i>& ptList,int zLevel)
  {
        //临时存储一份旧的Z值
      int srcLevel=zLevel;
      //对Z阶层做变换
      zLevel=Z_LEVEL_TRANS(zLevel);
        //遍历所有的点
        for(std::vector<cv::Point2i>::const_iterator iter=ptList.begin();
            iter!=ptList.end();++iter)
        {
            //判断z值是否在范围内
            if(!inZMapRange(iter->x,iter->y)) return true;
            //获取当前位置的z值
            uchar& dstZLevel=imgZMap_.at<uchar>(*iter);
            //判断它的下一个位置是否更大
            //这里是大于还是大于等于，区别会很大，如果是大于的话，填充应该会比较紧密，但是时间开销也大
            if(dstZLevel>zLevel)
                continue;
            //将当前位置记录下来
            dstZLevel=zLevel;
        }
        //如果z值已经是0了，也返回true
        return srcLevel==0;
  }

  //获取某个位置的z值
  ZType& getZLevel(double pixelX,double pixelY)
  {
      //判断图片是否在范围内
      if(!inImgRange(pixelX,pixelY)) throw  ERROR_IMG_OUT_RANGE;
        //获取指定位置的z值
      return imgZMap_.at<ZType>(std::round(pixelY),std::round(pixelX));
  }

  //判断传入的z值是否为可以访问的
  bool judgeZAvailable(double pixelX,double pixelY,int zLevel)
  {
      //返回比较结果
      return zLevel>=getZLevel(pixelX,pixelY);
  }

  //传入一个浮点数，然后计算它附近的整数范围
  static void floatSournd(double xValue,double yValue,ImgCornerRange& dstRange)
  {
        //左右
        dstRange[0]=std::floor(xValue);
        dstRange[1]=std::ceil(xValue);
        //上下
        dstRange[2]=std::floor(yValue);
        dstRange[3]=std::ceil(yValue);
  }

  //根据传入的范围获取4个角位置上的像素
  void getSourndColor(const ImgCornerRange& sourndRange,
                      ColorPtrArray& colorList)//从左上开始顺时针记录
  {
        colorList[0]=&(colorAt(sourndRange[0],sourndRange[2]));//左上
        colorList[1]=&(colorAt(sourndRange[0],sourndRange[2]));//右上
        colorList[2]=&(colorAt(sourndRange[0],sourndRange[2]));//右下
        colorList[3]=&(colorAt(sourndRange[0],sourndRange[2]));//左下
  }

  //不使用双线性内插是因为图片的像素已经够高了
  //使用双线性内插会增大很多计算量，但最后效果的增长量可能不会很多
  //这里的函数不保证坐标的像素范围内
  //如果坐标不在像素范围内，尽管throw
  bool colorAt(double pixelX,double pixelY,cv::Vec3d& dstColor,int zLevel=-1)
  {
        //判断坐标是否在图片的范围内
        if(!inImgRange(pixelX,pixelY)) return false;
        //获取目标位置的颜色
        const cv::Vec3b& tempColor=colorAt(pixelX,pixelY,zLevel);
        //遍历颜色，复制每个通道
        for(int i=0;i<3;++i) dstColor[i]=tempColor[i];
        //正常返回
        return true;
  }

  //释放图片内存
  void releaseImg()
  {
      //判断图片是否为空
      if(cvImg_.empty()) return;
      //释放内存
      cvImg_.release();
  }

  //把若干点画到图片上，然后输出图片的结果
  void drawPoints(const std::vector<cv::Point2i> &ptList)
  {
      //获取图片
      cv::Mat &mainImg=loadImg();
      //目标图片
      cv::Mat dstImg=mainImg.clone();
      //把线画在图片上，画线是为了方便观察
      cv::polylines(dstImg,ptList,false,cv::Scalar(255,0,0),3);
      //保存图片
      cv::imwrite(s_Img_path,dstImg);
  }



  //当前图片所能看到的特征点列表，仅仅记录每个点云中点的标号
  //想要获取这个点在图片上的投影，需要去点列表里面的obvs属性里面寻找.
    PtIdList ptIdList_;

  // Constructor (use unique index for the view_id)
  View(
    const std::string & sImgPath = "",
    IndexT view_id = UndefinedIndexT,
    IndexT intrinsic_id = UndefinedIndexT,
    IndexT pose_id = UndefinedIndexT,
    IndexT width = UndefinedIndexT, IndexT height = UndefinedIndexT)
    :s_Img_path(sImgPath), id_view(view_id), id_intrinsic(intrinsic_id),
    id_pose(pose_id), ui_width(width), ui_height(height)
    {
      //初始化全零的颜色
      zeroColor_=cv::Vec3b(0,0,0);
  }

  virtual ~View() = default;

  /**
  * Serialization out
  * @param ar Archive
  */
  template <class Archive>
  void save( Archive & ar ) const;

  /**
  * @brief Serialization in
  * @param ar Archive
  */
  template <class Archive>
  void load( Archive & ar );
};

/// Define a collection of View
using Views = Hash_Map<IndexT, std::shared_ptr<View> >;

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_VIEW_HPP
