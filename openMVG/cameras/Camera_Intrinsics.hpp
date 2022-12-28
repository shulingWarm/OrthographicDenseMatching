// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERAS_CAMERA_INTRINSICS_HPP
#define OPENMVG_CAMERAS_CAMERA_INTRINSICS_HPP

#include <vector>

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/numeric/numeric.h"
#include "openMVG/stl/hash.hpp"

#include "ceres/problem.h"
#include "ceres/solver.h"
#include"ceres/autodiff_cost_function.h"
//#include "openMVG/sfm/sfm_data_BA_ceres_camera_functor.hpp"
#include <ceres/rotation.h>
#include <ceres/types.h>


namespace openMVG
{
namespace cameras
{

/**
* @brief Struct used to force "clonability"
*/
template< typename T>
struct Clonable
{
  virtual T * clone() const = 0;
};

//针对ceres的用于约束xc的约束函数
struct XcConstrain
{
public:
    const double* imgPt_;//图像中的点坐标
    const std::vector<double>& intrinsicPar_;//内参

    //构造函数，初始化点坐标,记录内参
    XcConstrain(const double* const imgPoint,const std::vector<double> &intrList)
        : imgPt_(imgPoint) , intrinsicPar_(intrList)
    {

    }

    // 每个内参的索引
    enum : uint8_t {
      OFFSET_FOCAL_LENGTH = 0,
      OFFSET_PRINCIPAL_POINT_X = 1,
      OFFSET_PRINCIPAL_POINT_Y = 2,
      OFFSET_DISTO_K1 = 3,
      OFFSET_DISTO_K2 = 4,
      OFFSET_DISTO_K3 = 5,
      OFFSET_DISTO_T1 = 6,
      OFFSET_DISTO_T2 = 7,
    };

    //约束函数
    template <typename T>
    bool operator()(
      const T* const pos_proj,//相机坐标系下的投影点，也就是xcPoint
      T* out_residuals) const
    {
        const T& focal = T(intrinsicPar_[OFFSET_FOCAL_LENGTH]);
        const T& principal_point_x = T(intrinsicPar_[OFFSET_PRINCIPAL_POINT_X]);
        const T& principal_point_y = T(intrinsicPar_[OFFSET_PRINCIPAL_POINT_Y]);

        // Transform the point from homogeneous to euclidean (undistorted point)
        T x_u = pos_proj[0];
        T y_u = pos_proj[1];
        const T& k1 = T(intrinsicPar_[OFFSET_DISTO_K1]);
        const T& k2 = T(intrinsicPar_[OFFSET_DISTO_K2]);
        const T& k3 = T(intrinsicPar_[OFFSET_DISTO_K3]);
        const T& t1 = T(intrinsicPar_[OFFSET_DISTO_T1]);
        const T& t2 = T(intrinsicPar_[OFFSET_DISTO_T2]);

        // Apply distortion (xd,yd) = disto(x_u,y_u)
        const T r2 = x_u*x_u + y_u*y_u;
        const T r4 = r2 * r2;
        const T r6 = r4 * r2;
        const T r_coeff = (1.0 + k1*r2 + k2*r4 + k3*r6);
        const T t_x = t2 * (r2 + 2.0 * x_u*x_u) + 2.0 * t1 * x_u * y_u;
        const T t_y = t1 * (r2 + 2.0 * y_u*y_u) + 2.0 * t2 * x_u * y_u;
        x_u = x_u * r_coeff + t_x;
        y_u = y_u * r_coeff + t_y;

        // Apply focal length and principal point to get the final image coordinates
        const T projected_x = principal_point_x + focal * x_u;
        const T projected_y = principal_point_y + focal * y_u;

        out_residuals[0] = projected_x - imgPt_[0];
        out_residuals[1] = projected_y - imgPt_[1];

        //另外一个强制为1
        out_residuals[2]=pos_proj[2]-T(1);
        return true;
    }
};

/**
* @brief Base class used to store common intrinsics parameters
*/
struct IntrinsicBase : public Clonable<IntrinsicBase>
{
  /// Width of image
  unsigned int w_;
  /// Height of image
  unsigned int h_;

  /**
  * @brief Constructor
  * @param w Width of the image
  * @param h Height of the image
  */
  IntrinsicBase( unsigned int w = 0, unsigned int h = 0 )
    : w_( w ),
      h_( h )
  {

  }

  /**
  * @brief Destructor
  */
  virtual ~IntrinsicBase() = default;

  /**
  * @brief Get width of the image
  * @return width of the image
  */
  unsigned int w() const
  {
    return w_;
  }

  /**
  * @brief Get height of the image
  * @return height of the image
  */
  unsigned int h() const
  {
    return h_;
  }

  /**
  * @brief Compute projection of a 3D point into the image plane
  * (Apply pose, disto (if any) and Intrinsics)
  * @param pose Pose used to compute projection
  * @param pt3D 3D-point to project on image plane
  * @return Projected (2D) point on image plane
  */
  virtual Vec2 project(
    const geometry::Pose3 & pose,
    const Vec3 & pt3D ) const
  {
    const Vec3 X = pose( pt3D ); // apply pose
    if ( this->have_disto() ) // apply disto & intrinsics
    {
      return this->cam2ima( this->add_disto( X.hnormalized() ) );
    }
    else // apply intrinsics
    {
      return this->cam2ima( X.hnormalized() );
    }
  }

  /**
  * @brief Compute the residual between the 3D projected point and an image observation
  * @param pose Pose used to project point on camera plane
  * @param X 3d point to project on camera plane
  * @param x image observation
  * @brief Relative 2d distance between projected and observed points
  */
  Vec2 residual(
    const geometry::Pose3 & pose,
    const Vec3 & X,
    const Vec2 & x ) const
  {
    const Vec2 proj = this->project( pose, X );
    return x - proj;
  }

  // --
  // Virtual members
  // --

  /**
  * @brief Tell from which type the embed camera is
  * @return Corresponding intrinsic
  */
  virtual EINTRINSIC getType() const = 0;

  /**
  * @brief Data wrapper for non linear optimization (get data)
  * @return vector of parameter of this intrinsic
  */
  virtual std::vector<double> getParams() const = 0;

  //获取传入的点在相机坐标系下的坐标
  virtual Vec3 getXcPoint(double xPixel,double yPixel)
  {
      //获取当前的内参序列
      std::vector<double> thisParList=getParams();
      //初始化用于返回的坐标值
      Vec3 retValue;
      //根据内参数据计算
      retValue<<(xPixel-thisParList[1])/thisParList[0],
              (yPixel-thisParList[2])/thisParList[0],1.f;
      //实时去畸变太慢了
      return retValue;
      //新建ceres的待优化模型
    ceres::Problem problem;
    //为了构造约束函数，临时弄成数组的形式
    double imgPoint[2]={xPixel,yPixel};
    //新建损失函数
    ceres::CostFunction *costFunc=
            new ceres::AutoDiffCostFunction<XcConstrain,3,3>(new XcConstrain(imgPoint,thisParList));
    //添加待求解的参数
    problem.AddResidualBlock(costFunc,nullptr,retValue.data());
    //求解时的选项
    ceres::Solver::Options ceres_config_options;
    ceres_config_options.max_num_iterations = 100;
    ceres_config_options.logging_type = ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
    //求解的时候不能没有summary，不能给空指针
    ceres::Solver::Summary summary;
    //准备求解
    ceres::Solve(ceres_config_options,&problem,&summary);
      return retValue;
  }

  //获取内参矩阵,仅仅适用于普通的小孔成像模型
  virtual Mat3 getIntrinsicMatrix() const
  {
      //获取内参序列
      std::vector<double> inP=getParams();
      //新建用于返回的内参矩阵
      Mat3 matReturn;
      //记录焦距和主点坐标
      matReturn<<inP[0],0,inP[1],
                                  0,inP[0],inP[2],
                                  0,0,1;
      return matReturn;
  }

  //把传入的相机坐标系下的点投影到图片上
  virtual Vec2 projectXcPoint(const Vec3 &xcPoint)
  {
      //如果点不在相机的前面，返回-1 -1
      if(xcPoint[2]<=0)
      {
          Vec2 retValue;
          retValue<<-1,-1;
          return retValue;
      }
      //获取内参列表
      std::vector<double> intrVec=getParams();
      //初始化返回的结果
      Vec2 retValue;
      //计算具体的坐标
      retValue<<(xcPoint[0]*intrVec[0]/xcPoint[2]+intrVec[1]),
              (xcPoint[1]*intrVec[0]/xcPoint[2]+intrVec[2]);
      //返回结果
      return retValue;
  }

  /**
  * @brief Data wrapper for non linear optimization (update from data)
  * @param params List of params used to update this intrinsic
  * @retval true if update is correct
  * @retval false if there was an error during update
  */
  virtual bool updateFromParams( const std::vector<double> & params ) = 0;

  //图片的大小更改后，需要根据图片的大小更改相应地把内参也更改
  void resizeIntrByImg(double scaleRate)
  {
      //获取目前的内参
      std::vector<double> params=getParams();
      //将前3个数值缩放一下，分别是焦距和x,y的主点
      for(int parCount=0;parCount<3;++parCount) params[parCount]*=scaleRate;
      //根据这里的数值更新内参
      updateFromParams(params);
  }

  /**
  * @brief Return the list of parameter indexes that must be held constant
  * @param parametrization The given parametrization
  */
  virtual std::vector<int> subsetParameterization(
    const Intrinsic_Parameter_Type & parametrization) const = 0;

  /**
  * @brief Get bearing vector of a point given an image coordinate
  * @return bearing vector
  */
  virtual Vec3 operator () ( const Vec2& p ) const = 0;

  /**
  * @brief Transform a point from the camera plane to the image plane
  * @param p Camera plane point
  * @return Point on image plane
  */
  virtual Vec2 cam2ima( const Vec2& p ) const = 0;

  /**
  * @brief Transform a point from the image plane to the camera plane
  * @param p Image plane point
  * @return camera plane point
  */
  virtual Vec2 ima2cam( const Vec2& p ) const = 0;

  /**
  * @brief Does the camera model handle a distortion field?
  * @retval true if intrinsic holds distortion
  * @retval false if intrinsic does not hold distortion
  */
  virtual bool have_disto() const
  {
    return false;
  }

  /**
  * @brief Add the distortion field to a point (that is in normalized camera frame)
  * @param p Point before distortion computation (in normalized camera frame)
  * @return point with distortion
  */
  virtual Vec2 add_disto( const Vec2& p ) const = 0;

  /**
  * @brief Remove the distortion to a camera point (that is in normalized camera frame)
  * @param p Point with distortion
  * @return Point without distortion
  */
  virtual Vec2 remove_disto( const Vec2& p ) const  = 0;

  /**
  * @brief Return the un-distorted pixel (with removed distortion)
  * @param p Input distorted pixel
  * @return Point without distortion
  */
  virtual Vec2 get_ud_pixel( const Vec2& p ) const = 0;

  /**
  * @brief Return the distorted pixel (with added distortion)
  * @param p Input pixel
  * @return Distorted pixel
  */
  virtual Vec2 get_d_pixel( const Vec2& p ) const = 0;

  /**
  * @brief Normalize a given unit pixel error to the camera plane
  * @param value Error in image plane
  * @return error of passing from the image plane to the camera plane
  */
  virtual double imagePlane_toCameraPlaneError( double value ) const = 0;

  /**
  * @brief Return the projection matrix (interior & exterior) as a simplified projective projection
  * @param pose Extrinsic matrix
  * @return Concatenation of intrinsic matrix and extrinsic matrix
  */
  virtual Mat34 get_projective_equivalent( const geometry::Pose3 & pose ) const = 0;
  
  /**
  * @brief Serialization out
  * @param ar Archive
  */
  template <class Archive>
  void save( Archive & ar ) const;


  /**
  * @brief  Serialization in
  * @param ar Archive
  */
  template <class Archive>
  void load( Archive & ar );

  /**
  * @brief Generate a unique Hash from the camera parameters (used for grouping)
  * @return Hash value
  */
  virtual std::size_t hashValue() const
  {
    size_t seed = 0;
    stl::hash_combine( seed, static_cast<int>( this->getType() ) );
    stl::hash_combine( seed, w_ );
    stl::hash_combine( seed, h_ );
    const std::vector<double> params = this->getParams();
    for ( const auto & param : params )
      stl::hash_combine( seed , param );
    return seed;
  }
};


/**
* @brief Compute angle between two bearing rays
* Bearing rays are computed from position on image plane in each cameras
*
* @param pose1 Pose of the first camera
* @param intrinsic1 Intrinsic of the first camera
* @param pose2 Pose of the second camera
* @param intrinsic2 Intrinsic of the second camera
* @param x1 Image coordinate of a point in first camera
* @param x2 Image coordinate of a point in the second camera
*
* @return Angle (in degree) between the two rays
*/
inline double AngleBetweenRay(
  const geometry::Pose3 & pose1,
  const IntrinsicBase * intrinsic1,
  const geometry::Pose3 & pose2,
  const IntrinsicBase * intrinsic2,
  const Vec2 & x1, const Vec2 & x2 )
{
  // x = (u, v, 1.0)  // image coordinates
  // X = R.t() * K.inv() * x + C // Camera world point
  // getting the ray:
  // ray = X - C = R.t() * K.inv() * x
  const Vec3 ray1 = ( pose1.rotation().transpose() * intrinsic1->operator()( x1 ) ).normalized();
  const Vec3 ray2 = ( pose2.rotation().transpose() * intrinsic2->operator()( x2 ) ).normalized();
  const double mag = ray1.norm() * ray2.norm();
  const double dotAngle = ray1.dot( ray2 );
  return R2D( acos( clamp( dotAngle / mag, -1.0 + 1.e-8, 1.0 - 1.e-8 ) ) );
}

} // namespace cameras
} // namespace openMVG

#endif // #ifndef OPENMVG_CAMERAS_CAMERA_INTRINSICS_HPP
