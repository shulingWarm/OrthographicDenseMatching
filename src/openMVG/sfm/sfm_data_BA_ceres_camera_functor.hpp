// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_DATA_BA_CERES_CAMERA_FUNCTOR_HPP
#define OPENMVG_SFM_SFM_DATA_BA_CERES_CAMERA_FUNCTOR_HPP

#include <ceres/ceres.h>
#include <ceres/rotation.h>

#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/cameras/Camera_Pinhole_Radial.hpp"
#include "openMVG/cameras/Camera_Pinhole_Brown.hpp"
#include "openMVG/cameras/Camera_Pinhole_Fisheye.hpp"
#include "openMVG/cameras/Camera_Spherical.hpp"

//--
//- Define ceres Cost_functor for each OpenMVG camera model
//--

namespace openMVG {
namespace sfm {


/// Decorator used to Weight a given cost camera functor
/// i.e useful to weight GCP (Ground Control Points)
template <typename CostFunctor>
struct WeightedCostFunction
{
  WeightedCostFunction() :weight_(1.0) {}

  explicit WeightedCostFunction
  (
    CostFunctor * func,
    const double weight
  )
    :functor_(func), weight_(weight)
  {}

  template <typename T>
  bool operator()
  (
    const T* const cam_intrinsic,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals
  ) const
  {
    if (functor_->operator()(cam_intrinsic, cam_extrinsics, pos_3dpoint, out_residuals))
    {
      // Reweight the residual values
      for (int i = 0; i < CostFunctor::num_residuals(); ++i)
      {
        out_residuals[i] *= T(weight_);
      }
      return true;
    }
    return false;
  }

  template <typename T>
  bool operator()
  (
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals
  ) const
  {
    if (functor_->operator()(cam_extrinsics, pos_3dpoint, out_residuals))
    {
      // Reweight the residual values
      for (int i = 0; i < CostFunctor::num_residuals(); ++i)
      {
        out_residuals[i] *= T(weight_);
      }
      return true;
    }
    return false;
  }

  ceres::internal::scoped_ptr<CostFunctor> functor_;
  const double weight_;
};

/**
 * @brief Ceres functor to use a Pinhole_Intrinsic (pinhole camera model K[R[t]) and a 3D point.
 *
 *  Data parameter blocks are the following <2,3,6,3>
 *  - 2 => dimension of the residuals,
 *  - 3 => the intrinsic data block [focal, principal point x, principal point y],
 *  - 6 => the camera extrinsic data block (camera orientation and position) [R;t],
 *         - rotation(angle axis), and translation [rX,rY,rZ,tx,ty,tz].
 *  - 3 => a 3D point data block.
 *
 */
struct ResidualErrorFunctor_Pinhole_Intrinsic
{
  ResidualErrorFunctor_Pinhole_Intrinsic(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum : uint8_t {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2
  };

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y] )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--



    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];

    T pos_proj[3];
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_proj);

    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    // Transform the point from homogeneous to euclidean (undistorted point)
    const T x_u = pos_proj[0] / pos_proj[2];
    const T y_u = pos_proj[1] / pos_proj[2];

    //--
    // Apply intrinsic parameters
    //--

    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal * x_u;
    const T projected_y = principal_point_y + focal * y_u;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - m_pos_2dpoint[0];
    out_residuals[1] = projected_y - m_pos_2dpoint[1];

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic, 2, 3, 6, 3>(
            new ResidualErrorFunctor_Pinhole_Intrinsic(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic>, 2, 3, 6, 3>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic>
            (new ResidualErrorFunctor_Pinhole_Intrinsic(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};

struct ResidualErrorFunctor_Pinhole_Intrinsic_WaterXY
{
  ResidualErrorFunctor_Pinhole_Intrinsic_WaterXY(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum : uint8_t {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2
  };

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y] )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_extrinsics;
    const T * cam_C = &cam_extrinsics[3];

    //prepare l
    T R[9];
    ceres::AngleAxisToRotationMatrix(cam_R, R);

    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];

    T x1_c[3];
    x1_c[0] = (m_pos_2dpoint[0] - principal_point_x) /focal;
    x1_c[1] = (m_pos_2dpoint[1] - principal_point_y) /focal;
    x1_c[2] = (T) 1.0;

    const T Xrep_x = R[0]*x1_c[0] + R[1]*x1_c[1] + R[2] + cam_C[0];
    const T Xrep_y = R[3]*x1_c[0] + R[4]*x1_c[1] + R[5] + cam_C[1];
    //const T Xrep_z = R[6]*x1_c0 + R[7]*x1_c1 + R[8] + cam_C[2];


    const T L_XC_0 = cam_C[1] - pos_3dpoint[1];
    const T L_XC_1 = pos_3dpoint[0] - cam_C[0];
    const T L_XC_2 = cam_C[0]*pos_3dpoint[1] - cam_C[1]*pos_3dpoint[0];

    //const T d_Xrep_L = ceres::abs(Xrep_x*L_XC_0 + Xrep_y*L_XC_1 + L_XC_2);// / ceres::sqrt((L_XC_0*L_XC_0 + L_XC_1*L_XC_1));
    const T d_Xrep_L = ceres::abs(Xrep_x*L_XC_0 + Xrep_y*L_XC_1 + L_XC_2) / ceres::sqrt((L_XC_0*L_XC_0 + L_XC_1*L_XC_1));

    out_residuals[0] = d_Xrep_L;

    return true;
  }

  static int num_residuals() { return 1; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_WaterXY, 1, 3, 6, 3>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_WaterXY(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_WaterXY>, 1, 3, 6, 3>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_WaterXY>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_WaterXY(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};

//赵志豪 20210119 单独对Z做约束的时候使用的一种用来约束先验位姿的东西
struct ConstrainXYPrior
{
  Vec3 weight_;
  Vec3 pose_center_constraint_;
  Vec3 rotAxis_;//旋转角
  bool useRotAxis_=false;//是否对旋转角做约束
    bool needZ_;//是否需要对z做约束
  ConstrainXYPrior
  (
    const Vec3 & center,
    const Vec3 & weight,
    bool needConstrainZ=true //是否对z做约束
  ): weight_(weight), pose_center_constraint_(center),needZ_(needConstrainZ)
  {
      //标记不对旋转角做约束
      useRotAxis_=false;
  }

  //对旋转角做约束的情况
  ConstrainXYPrior
  (
    const Vec3 & center,
    const Vec3 & weight,
    const Vec3 &rotAxis,//旋转角
    bool needConstrainZ=true //是否对z做约束
  ): weight_(weight), pose_center_constraint_(center),needZ_(needConstrainZ)
  {
      //对旋转角做约束
      useRotAxis_=true;
      //记录旋转角
      rotAxis_=rotAxis;
  }

  template <typename T> bool
  operator()
  (
    const T* const extrainsicInfo, // 外参,3,4是光芯XY
    T* residuals
  )
  const
  {
    //计算光芯坐标和先验的光芯坐标的Z值的差
      residuals[0]=(extrainsicInfo[3]-pose_center_constraint_[0]);
      residuals[1]=extrainsicInfo[4]-pose_center_constraint_[1];
      //判断是否需要约束z
      if(needZ_)
        residuals[2]=(extrainsicInfo[5]-pose_center_constraint_[2]);
      else residuals[2]=T(0);
      //判断是否对旋转角做约束
      if(useRotAxis_)
      {
          //比较旋转角的差
          residuals[3]=extrainsicInfo[0]-rotAxis_[0];
          residuals[4]=extrainsicInfo[1]-rotAxis_[1];
          residuals[5]=extrainsicInfo[2]-rotAxis_[2];
      }

    return true;
  }
};


//用来对折射率做先验约束的优化结构体
struct ConstrinRefN
{
    double priorIdx_;//先验的折射率

    //构造函数
    ConstrinRefN(double refN):priorIdx_(refN){}

    template <typename T> bool
    operator()
    (
      const T* const refIdx, // 正在被优化的折射率
      T* residuals
    )
    const
    {

        //误差即为目前的量和先验量的区别
        residuals[0]=refIdx[0]-T(priorIdx_);
      return true;
    }
};


//赵志豪 20210119 单独对Z做约束的时候使用的一种用来约束先验位姿的东西
struct ConstrainZPrior
{
  Vec3 weight_;
  Vec3 pose_center_constraint_;

  ConstrainZPrior
  (
    const Vec3 & center,
    const Vec3 & weight
  ): weight_(weight), pose_center_constraint_(center)
  {
  }

  template <typename T> bool
  operator()
  (
    const T* const camZ, // 光芯坐标Z
    T* residuals
  )
  const
  {
    //计算光芯坐标和先验的光芯坐标的Z值的差
      residuals[0]=(camZ[0]-T(pose_center_constraint_[2]));

    return true;
  }
};



//赵志豪 20210116 对于已经约束了XY的坐标,进一步约束Z
struct WaterPlaneConstrainCostFunctionZ
{
    double oceanHeight_;//海面高度
    double posXy_[2];//已经经过约束的XY坐标
    double camXy_[2];//已经约束过的相机光芯的XY坐标
    double planeDis_;//世界坐标和光芯坐标的距离在水平面的投影
    double xcPoint_[3];//成像点在相机坐标系下的坐标
  WaterPlaneConstrainCostFunctionZ(const double* const pos_2dpoint,
                                    double oceanHeight,//海面高度
                                   Eigen::Matrix3d &rotaMat,//旋转矩阵
                                   Eigen::Vector3d &camPos,//光芯坐标
                                   Eigen::Vector3d &pointPos, //重建点坐标
                                   std::vector<double> &intrinsicInfo //内参信息
                                   )
  :m_pos_2dpoint(pos_2dpoint)
  {
      oceanHeight_=oceanHeight;//记录海面高度
      //成像点在相机坐标系下的坐标
      Eigen::Vector3d xcPoint;
      xcPoint<<pos_2dpoint[0]-intrinsicInfo[1],
              pos_2dpoint[1]-intrinsicInfo[2],intrinsicInfo[0];
      //记录成像点在相机坐标系下的坐标
      xcPoint_[0]=xcPoint[0];
      xcPoint_[1]=xcPoint[1];
      xcPoint_[2]=xcPoint[2];
      //记录已经约束过的XY坐标
      posXy_[0]=pointPos[0];
      posXy_[1]=pointPos[1];
      //记录相机光芯坐标
      camXy_[0]=camPos[0];
      camXy_[1]=camPos[1];
      //记录世界坐标和光芯坐标距离在水平面的投影
      planeDis_=(posXy_[0]-camXy_[0])*(posXy_[0]-camXy_[0])+
              (posXy_[1]-camXy_[1])*(posXy_[1]-camXy_[1]);
      planeDis_=ceres::sqrt(planeDis_);
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum : uint8_t {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2
  };

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y] )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const camZ,//相机光芯的Z坐标
    const T* const worldZ,//成像点Z坐标
    const T* const rotAxis,//旋转角
    T* out_residuals) const
  {
      //把旋转角转换为旋转矩阵
      //旋转矩阵的3个列向量分别是rotMat rotMat+3,rotMat+6
      T rotMat[9];
      ceres::AngleAxisToRotationMatrix(rotAxis,rotMat);
      //取成像点在相机坐标系下的坐标
      T xcPoint[3];
      xcPoint[0]=T(xcPoint_[0]);
      xcPoint[1]=T(xcPoint_[1]);
      xcPoint[2]=T(xcPoint_[2]);
      //旋转向量的3个列向量
      T* rotCol1=rotMat;
      T* rotCol2=rotMat+3;
      T* rotCol3=rotMat+6;
      //复制Z值
      T pointZ=worldZ[0];
      //判断当前Z的位置
      //if(pointZ>T(oceanHeight_))
      if(false)
      {
          //图片反向投影的方向向量
          T dirVec[3];
          dirVec[0]=ceres::DotProduct(rotCol1,xcPoint);
          dirVec[1]=ceres::DotProduct(rotCol2,xcPoint);
          dirVec[2]=ceres::DotProduct(rotCol3,xcPoint);
          //方向向量与法向量的夹角
          T cosValue=dirVec[2]/ceres::sqrt(ceres::DotProduct(dirVec,dirVec));
          T inputAngle=ceres::acos(ceres::abs(cosValue));
          //折射角
          T outputAngle=ceres::asin(ceres::sin(inputAngle)/1.33);
          //低于海面的高度
          T underHeight=T(oceanHeight_)-pointZ;
          //更改低于海面的高度
          underHeight=underHeight*ceres::tan(outputAngle)/ceres::tan(inputAngle);
          //更新Z值
          pointZ=T(oceanHeight_)-underHeight;
      }
      //r1和另外两个旋转向量的叉积
      T r12Cross[3],r13Cross[3];
      ceres::CrossProduct(rotCol1,rotCol2,r12Cross);
      ceres::CrossProduct(rotCol1,rotCol3,r13Cross);
      //按照比例加起来的向量
      T sumVec[3];
      for(int dimCount=0;dimCount<3;++dimCount)
      {
          sumVec[dimCount]=(posXy_[1]-camXy_[1])*r12Cross[dimCount]+
                  (pointZ-camZ[0])*r13Cross[dimCount];
      }
      //约束YZ
      out_residuals[0]=ceres::DotProduct(sumVec,xcPoint);
      //r2和另外两个旋转向量的叉乘
      T r21Cross[3],r23Cross[3];
      ceres::CrossProduct(rotCol2,rotCol1,r21Cross);
      ceres::CrossProduct(rotCol2,rotCol3,r23Cross);
      for(int dimCount=0;dimCount<3;++dimCount)
      {
          sumVec[dimCount]=(posXy_[0]-camXy_[0])*r21Cross[dimCount]+
                  (pointZ-camZ[0])*r23Cross[dimCount];
      }
      //约束XZ
      out_residuals[1]=ceres::DotProduct(sumVec,xcPoint);
      return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,double oceanHeight,
  Eigen::Matrix3d &rotaMat,//旋转矩阵
  Eigen::Vector3d &camPos,//光芯坐标
  Eigen::Vector3d &pointPos, //重建点坐标
  std::vector<double> &intrinsicInfo, //内参信息
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <WaterPlaneConstrainCostFunctionZ, 2,1,1,3>(
            new WaterPlaneConstrainCostFunctionZ(observation.data(),oceanHeight,
                rotaMat,camPos,pointPos,intrinsicInfo)));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<WaterPlaneConstrainCostFunctionZ>, 2,1,1,3>
          (new WeightedCostFunction<WaterPlaneConstrainCostFunctionZ>
            (new WaterPlaneConstrainCostFunctionZ(observation.data(),oceanHeight,
                rotaMat,camPos,pointPos,intrinsicInfo),weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};

//赵志豪 20210120
//用重投影的方案约束Z
struct ReProjectConstrainZ
{
    double oceanHeight_;//海面高度
  ReProjectConstrainZ(const double* const pos_2dpoint,
                                    double oceanHeight)
  :m_pos_2dpoint(pos_2dpoint)
  {
      oceanHeight_=oceanHeight;
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum : uint8_t {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2
  };

  //ceres形式的叉积,因为涉及到模板类的问题,cost函数里面不方便用eigen
  //第1个向量和第2个向量的叉积,计算结果保存到dstVec里面
  //向量默认是已经开辟好空间的
  template<typename T>
  static void ceresCross(const T* vec1,const T* vec2,T* dstVec)
  {
    //按照叉积的公式来做计算
      dstVec[0]=vec1[1]*vec2[2]-vec1[2]*vec2[1];
      dstVec[1]=vec1[2]*vec2[0]-vec1[0]*vec2[2];
      dstVec[2]=vec1[0]*vec2[1]-vec1[1]*vec2[0];
  }

  //在计算残差的时候使用的点积,做这个函数的原因和上面的叉积一样
  template<typename T>
  static T ceresDot(const T* vec1,const T* vec2)
  {
      T sum=T(0);
      for(int vecCount=0;vecCount<3;++vecCount)
          sum=sum+vec1[vecCount]*vec2[vecCount];
      return sum;
  }

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y] )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,//内参
    const T* const cam_extrinsics,//外参,外参中的后三个是平移量
    const T* const pos_3dpoint,//重建后的点的位置
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--
    //取出外参的数值
    const T * cam_R = cam_extrinsics;
    const T * cam_C = &cam_extrinsics[3];
    T R[9];//P矩阵中的旋转部分
    ceres::AngleAxisToRotationMatrix(cam_R, R);
    //取出旋转矩阵的每个列向量
    T rotCol1[3],rotCol2[3],rotCol3[3];
    rotCol1[0]=R[0];rotCol1[1]=R[1];rotCol1[2]=R[2];
    rotCol2[0]=R[3];rotCol2[1]=R[4];rotCol2[2]=R[5];
    rotCol3[0]=R[6];rotCol3[1]=R[7];rotCol3[2]=R[8];
    //取出内参向量中的数值
    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];
    //把成像点转换为摄像机矩阵下的坐标
    T xcPoint[3];
    xcPoint[0]=(m_pos_2dpoint[0] - principal_point_x) /focal;
    xcPoint[1]=(m_pos_2dpoint[1] - principal_point_y) /focal;
    xcPoint[2]=T(1.0f);
    //XY坐标在水平面的投影距离
    T planeDis=(cam_C[0]-pos_3dpoint[0])*(cam_C[0]-pos_3dpoint[0])+
            (cam_C[1]-pos_3dpoint[1])*(cam_C[1]-pos_3dpoint[1]);
    planeDis=ceres::sqrt(planeDis);
    //反向投影的方向向量
    T dirVec[3];
    dirVec[0]=ceresDot(rotCol1,xcPoint);
    dirVec[1]=ceresDot(rotCol2,xcPoint);
    dirVec[2]=ceresDot(rotCol3,xcPoint);
    //方向向量和海平面法向量的夹角余弦
    T cosValue=dirVec[2]/ceres::sqrt(ceresDot(dirVec,dirVec));
    //方向向量与海平面法向量的夹角
    T inputAngle=ceres::acos(ceres::abs(cosValue));
    //期望的Z坐标
    T exceptionZ=cam_C[2]-planeDis/ceres::tan(inputAngle);
    //根据期望的Z坐标做出来重建点坐标
    T tempPos[3];
    tempPos[0]=pos_3dpoint[0];
    tempPos[1]=pos_3dpoint[1];
    tempPos[2]=exceptionZ;
    //把旋转向量应用到临时坐标上
    T afterRotate[3];
    ceres::AngleAxisRotatePoint(cam_R,tempPos,afterRotate);
    //把相机的光芯坐标应用到旋转上
    T transPos[3];
    ceres::AngleAxisRotatePoint(cam_R,cam_C,transPos);
    //临时坐标的投影结果
    T projectPos[3];
    projectPos[0]=afterRotate[0]-transPos[0];
    projectPos[1]=afterRotate[1]-transPos[1];
    projectPos[2]=afterRotate[2]-transPos[2];
    //计算重投影误差
    out_residuals[1]=projectPos[0]/projectPos[2]-xcPoint[0];
    out_residuals[2]=projectPos[1]/projectPos[2]-xcPoint[1];
    //判断期望的Z是否低于海平面
    if(false)
    {
        //低于海平面的高度
        T underHeight=oceanHeight_-exceptionZ;
        //折射角
        T fractionAngle=ceres::asin(ceres::sin(inputAngle)/1.33);
        //实际的低于海平面的高度
        underHeight=underHeight*ceres::tan(inputAngle)/ceres::sin(fractionAngle);
        //实际的Z期望值
        exceptionZ=oceanHeight_-underHeight;
    }
    //误差为Z和期望Z的差
    out_residuals[0]=pos_3dpoint[2]-exceptionZ;
    return true;
  }

  static int num_residuals() { return 3; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,double oceanHeight,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ReProjectConstrainZ, 3, 3, 6, 3>(
            new ReProjectConstrainZ(observation.data(),oceanHeight)));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ReProjectConstrainZ>, 3, 3, 6, 3>
          (new WeightedCostFunction<ReProjectConstrainZ>
            (new ReProjectConstrainZ(observation.data(),oceanHeight),weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};


//用类似约束XY的方法来约束Z,因为使用反射投影的方式可能不太稳定
struct ConstrainZSimilarXY
{
    double oceanHeight_;//海面高度
    //反方向的Z,有时候海面的法向可能是0,0,-1,适用于这种情况
    const static bool oppositeZ=true;
  ConstrainZSimilarXY(const double* const pos_2dpoint,
                                    double oceanHeight)
  :m_pos_2dpoint(pos_2dpoint)
  {
      oceanHeight_=oceanHeight;
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum : uint8_t {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2
  };

  //ceres形式的叉积,因为涉及到模板类的问题,cost函数里面不方便用eigen
  //第1个向量和第2个向量的叉积,计算结果保存到dstVec里面
  //向量默认是已经开辟好空间的
  template<typename T>
  static void ceresCross(const T* vec1,const T* vec2,T* dstVec)
  {
    //按照叉积的公式来做计算
      dstVec[0]=vec1[1]*vec2[2]-vec1[2]*vec2[1];
      dstVec[1]=vec1[2]*vec2[0]-vec1[0]*vec2[2];
      dstVec[2]=vec1[0]*vec2[1]-vec1[1]*vec2[0];
  }

  //在计算残差的时候使用的点积,做这个函数的原因和上面的叉积一样
  template<typename T>
  static T ceresDot(const T* vec1,const T* vec2)
  {
      T sum=T(0);
      for(int vecCount=0;vecCount<3;++vecCount)
          sum=sum+vec1[vecCount]*vec2[vecCount];
      return sum;
  }

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y] )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,//内参
    const T* const cam_extrinsics,//外参,外参中的后三个是平移量
    const T* const pos_3dpoint,//重建后的点的位置
    const T* const refIndex,//折射率
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--
    //取出外参的数值
    const T * cam_R = cam_extrinsics;
    const T * cam_C = &cam_extrinsics[3];
    T R[9];//P矩阵中的旋转部分
    ceres::AngleAxisToRotationMatrix(cam_R, R);
    //取出旋转矩阵的每个列向量
    T rotCol1[3],rotCol2[3],rotCol3[3];
    rotCol1[0]=R[0];rotCol1[1]=R[1];rotCol1[2]=R[2];
    rotCol2[0]=R[3];rotCol2[1]=R[4];rotCol2[2]=R[5];
    rotCol3[0]=R[6];rotCol3[1]=R[7];rotCol3[2]=R[8];
    //取出内参向量中的数值
    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];
    //把成像点转换为摄像机矩阵下的坐标
    T xcPoint[3];
    xcPoint[0]=(m_pos_2dpoint[0] - principal_point_x) /focal;
    xcPoint[1]=(m_pos_2dpoint[1] - principal_point_y) /focal;
    xcPoint[2]=T(1.0f);
    //复制Z的数据
    T pointZ=pos_3dpoint[2];
    //判断Z是否小于海面高度
//    if((pointZ>T(oceanHeight_)&&oppositeZ)||
//            (pointZ<T(oceanHeight_)&&(!oppositeZ)))
    if(false)
    {
        //反向投影的方向向量
        T dirVec[3];
        dirVec[0]=ceresDot(rotCol1,xcPoint);
        dirVec[1]=ceresDot(rotCol2,xcPoint);
        dirVec[2]=ceresDot(rotCol3,xcPoint);
        //方向向量和海平面法向量的夹角余弦
        T cosValue=dirVec[2]/ceres::sqrt(ceresDot(dirVec,dirVec));
        //方向向量与海平面法向量的夹角
        T inputAngle=ceres::acos(ceres::abs(cosValue));
        //低于海面的高度
        T underHeight=T(oceanHeight_)-pointZ;
        //折射角
        T outputAngle=ceres::asin(ceres::sin(inputAngle)/refIndex[0]);
        //更新低于海面的高度
        underHeight=underHeight*ceres::tan(outputAngle)/ceres::tan(inputAngle);
        //更新Z值
        pointZ=T(oceanHeight_)-underHeight;
    }
    //用R1和另外两个旋转向量叉乘
    T r12Cross[3],r13Cross[3];
    ceresCross(rotCol1,rotCol2,r12Cross);
    ceresCross(rotCol1,rotCol3,r13Cross);
    //叉乘结果按照比较相加,根据约束公式
    T sumVec[3];
    for(int vecCount=0;vecCount<3;++vecCount)
    {
        sumVec[vecCount]=(pos_3dpoint[1]-cam_C[1])*r12Cross[vecCount]+
                (pointZ-cam_C[2])*r13Cross[vecCount];
    }
    //约束YZ
    out_residuals[0]=ceresDot(sumVec,xcPoint);
    //用r2和另外两个旋转向量相乘
    T r21Cross[3],r23Cross[3];
    ceresCross(rotCol2,rotCol1,r21Cross);
    ceresCross(rotCol2,rotCol3,r23Cross);
    //按照叉乘结果相加
    for(int vecCount=0;vecCount<3;++vecCount)
    {
        sumVec[vecCount]=(pos_3dpoint[0]-cam_C[0])*r21Cross[vecCount]+
                (pointZ-cam_C[2])*r23Cross[vecCount];
    }
    //约束XZ
    out_residuals[1]=ceresDot(sumVec,xcPoint);
    return true;
    //根据pointZ做出来的虚拟点
    T pointPos[3];
    pointPos[0]=pos_3dpoint[0];
    pointPos[1]=pos_3dpoint[1];
    pointPos[2]=pointZ;
    //把旋转应用到临时点上
    T afterRotate[3];
    ceres::AngleAxisRotatePoint(cam_R,pointPos,afterRotate);
    //把旋转应用到光芯坐标上
    T rotateCam[3];
    ceres::AngleAxisRotatePoint(cam_R,cam_C,rotateCam);
    //最终的重投影坐标
    T reprojectPos[3];
    reprojectPos[0]=afterRotate[0]-rotateCam[0];
    reprojectPos[1]=afterRotate[1]-rotateCam[1];
    reprojectPos[2]=afterRotate[2]-rotateCam[2];
    //计算重投影误差
    out_residuals[0]=reprojectPos[0]/reprojectPos[2]-xcPoint[0];
    out_residuals[1]=reprojectPos[1]/reprojectPos[2]-xcPoint[1];
    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,double oceanHeight,
    const double weight = 0.0
  )
  {
      return
        (new ceres::AutoDiffCostFunction
          <ConstrainZSimilarXY, 2, 3, 6, 3,1>(
            new ConstrainZSimilarXY(observation.data(),oceanHeight)));
  }

  const double * m_pos_2dpoint; // The 2D observation
};

//相机坐标系投影到图片上时的约束关系，为了去畸变做的权宜之计
struct XcPointProjectError
{
    //false是时候不使用去畸变
    const static bool useUndistortion=true;
  XcPointProjectError(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {

  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
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

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y] )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,//相机内参
    const T* const pos_proj,//相机坐标系下的投影点
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--
      const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
      const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
      const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];

      // Transform the point from homogeneous to euclidean (undistorted point)
      T x_u = pos_proj[0];
      T y_u = pos_proj[1];

    //如果不使用去畸变算法，就直接计算最后的投影结果
      if(useUndistortion)
      {
          const T& k1 = cam_intrinsics[OFFSET_DISTO_K1];
          const T& k2 = cam_intrinsics[OFFSET_DISTO_K2];
          const T& k3 = cam_intrinsics[OFFSET_DISTO_K3];
          const T& t1 = cam_intrinsics[OFFSET_DISTO_T1];
          const T& t2 = cam_intrinsics[OFFSET_DISTO_T2];

          // Apply distortion (xd,yd) = disto(x_u,y_u)
          const T r2 = x_u*x_u + y_u*y_u;
          const T r4 = r2 * r2;
          const T r6 = r4 * r2;
          const T r_coeff = (1.0 + k1*r2 + k2*r4 + k3*r6);
          const T t_x = t2 * (r2 + 2.0 * x_u*x_u) + 2.0 * t1 * x_u * y_u;
          const T t_y = t1 * (r2 + 2.0 * y_u*y_u) + 2.0 * t2 * x_u * y_u;
          x_u = x_u * r_coeff + t_x;
          y_u = y_u * r_coeff + t_y;


      }

      // Apply focal length and principal point to get the final image coordinates
      const T projected_x = principal_point_x + focal * x_u;
      const T projected_y = principal_point_y + focal * y_u;


      // Compute and return the error is the difference between the predicted
      //  and observed position
      out_residuals[0] = projected_x - m_pos_2dpoint[0];
      out_residuals[1] = projected_y - m_pos_2dpoint[1];

      //另外一个强制为1
      out_residuals[2]=pos_proj[2]-T(1);
      return true;
  }

  static int num_residuals() {
      return 2;
  }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <XcPointProjectError, 3,8,3>(
            new XcPointProjectError(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<XcPointProjectError>, 3,8,3>
          (new WeightedCostFunction<XcPointProjectError>
            (new XcPointProjectError(observation.data()),weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};



struct XyConstrainWithVaribleXc
{
  XyConstrainWithVaribleXc()
  {

  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum : uint8_t {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2
  };

  //ceres形式的叉积,因为涉及到模板类的问题,cost函数里面不方便用eigen
  //第1个向量和第2个向量的叉积,计算结果保存到dstVec里面
  //向量默认是已经开辟好空间的
  template<typename T>
  static void ceresCross(const T* vec1,const T* vec2,T* dstVec)
  {
    //按照叉积的公式来做计算
      dstVec[0]=vec1[1]*vec2[2]-vec1[2]*vec2[1];
      dstVec[1]=vec1[2]*vec2[0]-vec1[0]*vec2[2];
      dstVec[2]=vec1[0]*vec2[1]-vec1[1]*vec2[0];
  }

  //在计算残差的时候使用的点积,做这个函数的原因和上面的叉积一样
  template<typename T>
  static T ceresDot(const T* vec1,const T* vec2)
  {
      T sum=T(0);
      for(int vecCount=0;vecCount<3;++vecCount)
          sum=sum+vec1[vecCount]*vec2[vecCount];
      return sum;
  }

  template <typename T>
  bool operator()(
    const T* const xcPoint,//成像点在相机坐标系下的位置
    const T* const cam_extrinsics,//外参,外参中的后三个是平移量
    const T* const pos_3dpoint,//重建后的点的位置
    T* out_residuals) const
  {
    //取出外参的数值
    const T * cam_R = cam_extrinsics;
    const T * cam_C = &cam_extrinsics[3];
    T R[9];//P矩阵中的旋转部分
    ceres::AngleAxisToRotationMatrix(cam_R, R);
    //取出旋转矩阵的每个列向量
    T rotCol1[3],rotCol2[3],rotCol3[3];
    rotCol1[0]=R[0];rotCol1[1]=R[1];rotCol1[2]=R[2];
    rotCol2[0]=R[3];rotCol2[1]=R[4];rotCol2[2]=R[5];
    rotCol3[0]=R[6];rotCol3[1]=R[7];rotCol3[2]=R[8];
    //根据公式,计算旋转向量之间的叉积,因为对正负不确定,所以就老老实实地做叉积吧
    T r1Cross[3],r2Cross[3];
    ceresCross<T>(rotCol3,rotCol1,r1Cross);
    ceresCross<T>(rotCol3,rotCol2,r2Cross);
    //把xc归一化
    T xcScale=ceres::sqrt(ceresDot(xcPoint,xcPoint));
    //把3个叉积后的向量按照比例加起来
    T sumVec[3];
    for(int vecCount=0;vecCount<3;++vecCount)
    {
        sumVec[vecCount]=(pos_3dpoint[0]-cam_C[0])*r1Cross[vecCount]+
                (pos_3dpoint[1]-cam_C[1])*r2Cross[vecCount];
    }
    //计算向量的幅值，把sumVec也归一化
    T sumScale=ceres::sqrt(ceresDot(sumVec,sumVec));
    //判断幅值为0
    if(sumScale==T(0))
    {
        out_residuals[0]=T(0);
        return true;
    }
    //根据公式计算最终的残差
    out_residuals[0]=ceresDot(sumVec,xcPoint)/sumScale/xcScale;
    return true;
  }

  static int num_residuals() {
      return 1;
  }


  static ceres::CostFunction* Create
  (
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <XyConstrainWithVaribleXc, 1, 3, 6, 3>(
            new XyConstrainWithVaribleXc()));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<XyConstrainWithVaribleXc>, 1, 3, 6, 3>
          (new WeightedCostFunction<XyConstrainWithVaribleXc>
            (new XyConstrainWithVaribleXc(),weight)));
    }
  }
};





//赵志豪 20201226 从论文里面引进的一种新算法,这个算法是用来约束XY的
//,约束Z的方法会在另一个地方体现
struct WaterPlaneConstrainCostFunctionXy
{
    double oceanHeight_;//海面高度
    const static bool useReProject_=false;//是否使用重投影误差

    //相机坐标系下的图像点
    double xcValue[2];
    //pos_2dpoint目前传入的是三维坐标，表示的是它在相机坐标系下的位置
  WaterPlaneConstrainCostFunctionXy(const double* const pos_2dpoint,
                                    double oceanHeight)
  :m_pos_2dpoint(pos_2dpoint)
  {
      oceanHeight_=oceanHeight;
      //记录相机坐标系下的点坐标
      xcValue[0]=pos_2dpoint[0];
      xcValue[1]=pos_2dpoint[1];
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum : uint8_t {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2
  };

  //ceres形式的叉积,因为涉及到模板类的问题,cost函数里面不方便用eigen
  //第1个向量和第2个向量的叉积,计算结果保存到dstVec里面
  //向量默认是已经开辟好空间的
  template<typename T>
  static void ceresCross(const T* vec1,const T* vec2,T* dstVec)
  {
    //按照叉积的公式来做计算
      dstVec[0]=vec1[1]*vec2[2]-vec1[2]*vec2[1];
      dstVec[1]=vec1[2]*vec2[0]-vec1[0]*vec2[2];
      dstVec[2]=vec1[0]*vec2[1]-vec1[1]*vec2[0];
  }

  //在计算残差的时候使用的点积,做这个函数的原因和上面的叉积一样
  template<typename T>
  static T ceresDot(const T* vec1,const T* vec2)
  {
      T sum=T(0);
      for(int vecCount=0;vecCount<3;++vecCount)
          sum=sum+vec1[vecCount]*vec2[vecCount];
      return sum;
  }

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y] )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,//内参
    const T* const cam_extrinsics,//外参,外参中的后三个是平移量
    const T* const pos_3dpoint,//重建后的点的位置
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--
    //取出外参的数值
    const T * cam_R = cam_extrinsics;
    const T * cam_C = &cam_extrinsics[3];
    T R[9];//P矩阵中的旋转部分
    ceres::AngleAxisToRotationMatrix(cam_R, R);
    //取出旋转矩阵的每个列向量
    T rotCol1[3],rotCol2[3],rotCol3[3];
    rotCol1[0]=R[0];rotCol1[1]=R[1];rotCol1[2]=R[2];
    rotCol2[0]=R[3];rotCol2[1]=R[4];rotCol2[2]=R[5];
    rotCol3[0]=R[6];rotCol3[1]=R[7];rotCol3[2]=R[8];
    //根据公式,计算旋转向量之间的叉积,因为对正负不确定,所以就老老实实地做叉积吧
    T r1Cross[3],r2Cross[3];
    ceresCross<T>(rotCol3,rotCol1,r1Cross);
    ceresCross<T>(rotCol3,rotCol2,r2Cross);
    //ceresCross<T>(rotCol3,cam_C,tCross);
    //相机坐标系下的坐标在这里是可以直接使用的 20210624
    T xcPoint[3];
    xcPoint[0]=T(xcValue[0]);
    xcPoint[1]=T(xcValue[1]);
    xcPoint[2]=T(1);
    //计算两个旋转向量与xc的点积
    T xcDot1=ceresDot(xcPoint,rotCol1);
    T xcDot2=ceresDot(xcPoint,rotCol2);
    //计算向量的幅值，用来计算点到直线的距离
    T sumScale=ceres::sqrt(xcDot1*xcDot1+xcDot2*xcDot2);
    //判断幅值为0
    if(sumScale==T(0))
    {
        out_residuals[0]=T(0);
        return true;
    }
    //根据公式计算最终的残差
    out_residuals[0]=((pos_3dpoint[0]-cam_C[0])*xcDot2-
            (pos_3dpoint[1]-cam_C[1])*xcDot1)/sumScale;
    //如果不使用重投影误差,就此返回
    if(!useReProject_) return true;
  }

  static int num_residuals() {
      return useReProject_?3:1;
  }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec3 & observation,double oceanHeight,//这里被临时换成了Vec3，表示的是相机坐标系下的坐标
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
        if(useReProject_)
        {
            return
              (new ceres::AutoDiffCostFunction
                <WaterPlaneConstrainCostFunctionXy, 3, 8, 6, 3>(
                  new WaterPlaneConstrainCostFunctionXy(observation.data(),oceanHeight)));
        }
      return
        (new ceres::AutoDiffCostFunction
          <WaterPlaneConstrainCostFunctionXy, 1, 8, 6, 3>(
            new WaterPlaneConstrainCostFunctionXy(observation.data(),oceanHeight)));
    }
    else
    {
        if(useReProject_)
        {
            return
              (new ceres::AutoDiffCostFunction
                <WeightedCostFunction<WaterPlaneConstrainCostFunctionXy>, 3, 8, 6, 3>
                (new WeightedCostFunction<WaterPlaneConstrainCostFunctionXy>
                  (new WaterPlaneConstrainCostFunctionXy(observation.data(),oceanHeight),weight)));
        }
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<WaterPlaneConstrainCostFunctionXy>, 1, 8, 6, 3>
          (new WeightedCostFunction<WaterPlaneConstrainCostFunctionXy>
            (new WaterPlaneConstrainCostFunctionXy(observation.data(),oceanHeight),weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};


struct ResidualErrorFunctor_Pinhole_Intrinsic_WaterZ
{
  ResidualErrorFunctor_Pinhole_Intrinsic_WaterZ(const double* const pos_2dpoint, const double* const _water_plane)
  :m_pos_2dpoint(pos_2dpoint), water_plane(_water_plane)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum : uint8_t {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2
  };

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y] )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_extrinsics;
    const T * cam_C = &cam_extrinsics[3];

    //prepare l
    T R[9];
    ceres::AngleAxisToRotationMatrix(cam_R, R);

    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];

    T x1_c[3];
    x1_c[0] = (m_pos_2dpoint[0] - principal_point_x) /focal;
    x1_c[1] = (m_pos_2dpoint[1] - principal_point_y) /focal;
    x1_c[2] = (T) 1.0;

    const T Xrep_x = R[0]*x1_c[0] + R[1]*x1_c[1] + R[2] + cam_C[0];
    const T Xrep_y = R[3]*x1_c[0] + R[4]*x1_c[1] + R[5] + cam_C[1];
    //const T Xrep_z = R[6]*x1_c0 + R[7]*x1_c1 + R[8] + cam_C[2];


    const T L_XC_0 = cam_C[1] - pos_3dpoint[1];
    const T L_XC_1 = pos_3dpoint[0] - cam_C[0];
    const T L_XC_2 = cam_C[0]*pos_3dpoint[1] - cam_C[1]*pos_3dpoint[0];

    //const T d_Xrep_L = ceres::abs(Xrep_x*L_XC_0 + Xrep_y*L_XC_1 + L_XC_2);// / ceres::sqrt((L_XC_0*L_XC_0 + L_XC_1*L_XC_1));
    const T d_Xrep_L = ceres::abs(Xrep_x*L_XC_0 + Xrep_y*L_XC_1 + L_XC_2) / ceres::sqrt((L_XC_0*L_XC_0 + L_XC_1*L_XC_1));

    out_residuals[0] = d_Xrep_L;

    //
    //
    const T n = (T) 1.33;
    const T * Za = (T*) water_plane;

    T nl1[3];
    nl1[0] = R[0]*x1_c[0] + R[1]*x1_c[1] + R[2];
    nl1[1] = R[3]*x1_c[0] + R[4]*x1_c[1] + R[5];
    nl1[2] = R[6]*x1_c[0] + R[7]*x1_c[1] + R[8];

    T tan_a1 = ceres::tan(ceres::acos( nl1[2] / ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1] + nl1[2]*nl1[2])));
    T a1 = ceres::acos((nl1[2]) / ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1] + nl1[2]*nl1[2]));
    if (a1 > ceres::acos(-1)/2)
    {
        a1 = ceres::acos(-1) - a1;
        tan_a1 = ceres::tan(a1);
    }
    const T tan_a2 = ceres::tan(ceres::asin(ceres::sin(a1)/n));

    const T Zp = pos_3dpoint[2];
    const T r_in = tan_a1 * (Za[0]-cam_C[2]);

    const T r_xy = ceres::sqrt((pos_3dpoint[0]-cam_C[0])*(pos_3dpoint[0]-cam_C[0]) + (pos_3dpoint[1]-cam_C[1])*(pos_3dpoint[1]-cam_C[1]) );
    const T Xq[3] = {(Za[0]-cam_C[2])/nl1[2]*nl1[0]+cam_C[0], (Za[0]-cam_C[2])/nl1[2]*nl1[1]+cam_C[1], Za[0]};
    const T l_ref[3] = {nl1[0]/ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1]), nl1[1]/ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1]), (T)1.0/tan_a2*(ceres::sin(a1)/n)*nl1[2]/ceres::abs(nl1[2])};

//    if (Zp <= Za[0])
    if (r_in >= r_xy)
    {
        const T X_subtract_C[3] = {pos_3dpoint[0]-cam_C[0], pos_3dpoint[1]-cam_C[1], pos_3dpoint[2]-cam_C[2]};
        const T crossXC[3] = {X_subtract_C[1]*nl1[2]-X_subtract_C[2]*nl1[1], X_subtract_C[2]*nl1[0]-X_subtract_C[0]*nl1[2], X_subtract_C[0]*nl1[1]-X_subtract_C[1]*nl1[0]};
        const T D_in = ceres::sqrt(crossXC[0]*crossXC[0] + crossXC[1]*crossXC[1] + crossXC[2]*crossXC[2]) / ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1] + nl1[2]*nl1[2]);

        //out_residuals[1] = r_xy-(r_in-tan_a1*(Za[0]-Zp));
        out_residuals[1] = r_in - r_xy-tan_a1*(Za[0]-Zp);

    }else{
        const T n_ref[3] = {nl1[0]/ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1]), nl1[1]/ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1]), (T)1.0/tan_a2*(ceres::sin(a1)/n)*nl1[2]/ceres::abs(nl1[2])};
        const T X_subtract_C[3] = {pos_3dpoint[0]-Xq[0], pos_3dpoint[1]-Xq[1], pos_3dpoint[2]-Xq[2]};
        const T crossXC[3] = {X_subtract_C[1]*n_ref[2]-X_subtract_C[2]*n_ref[1], X_subtract_C[2]*n_ref[0]-X_subtract_C[0]*n_ref[2], X_subtract_C[0]*n_ref[1]-X_subtract_C[1]*n_ref[0]};
        const T D_ref = ceres::sqrt(crossXC[0]*crossXC[0] + crossXC[1]*crossXC[1] + crossXC[2]*crossXC[2]) / ceres::sqrt(n_ref[0]*n_ref[0] + n_ref[1]*n_ref[1] + n_ref[2]*n_ref[2]);
        const T r_ref = tan_a2 * (Zp-Za[0]);

        out_residuals[1] = r_xy-r_in-r_ref;

    }

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double & water_plane,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_WaterZ, 2, 3, 6, 3>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_WaterZ(observation.data(), &water_plane)));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_WaterZ>, 2, 3, 6, 3>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_WaterZ>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_WaterZ(observation.data(), &water_plane), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
  const double * water_plane;
};


struct ResidualErrorFunctor_Pinhole_Intrinsic_WaterZN
{
  ResidualErrorFunctor_Pinhole_Intrinsic_WaterZN(const double* const pos_2dpoint, const double* const _water_plane)
  :m_pos_2dpoint(pos_2dpoint), water_plane(_water_plane)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum : uint8_t {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2
  };

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y] )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    //const T* const ref_N,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_extrinsics;
    const T * cam_C = &cam_extrinsics[3];

    //prepare l
    T R[9];
    ceres::AngleAxisToRotationMatrix(cam_R, R);

    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];

    T x1_c[3];
    x1_c[0] = (m_pos_2dpoint[0] - principal_point_x) /focal;
    x1_c[1] = (m_pos_2dpoint[1] - principal_point_y) /focal;
    x1_c[2] = (T) 1.0;

    const T Xrep_x = R[0]*x1_c[0] + R[1]*x1_c[1] + R[2] + cam_C[0];
    const T Xrep_y = R[3]*x1_c[0] + R[4]*x1_c[1] + R[5] + cam_C[1];
    //const T Xrep_z = R[6]*x1_c0 + R[7]*x1_c1 + R[8] + cam_C[2];


    const T L_XC_0 = cam_C[1] - pos_3dpoint[1];
    const T L_XC_1 = pos_3dpoint[0] - cam_C[0];
    const T L_XC_2 = cam_C[0]*pos_3dpoint[1] - cam_C[1]*pos_3dpoint[0];

    //const T d_Xrep_L = ceres::abs(Xrep_x*L_XC_0 + Xrep_y*L_XC_1 + L_XC_2);// / ceres::sqrt((L_XC_0*L_XC_0 + L_XC_1*L_XC_1));
    const T d_Xrep_L = ceres::abs(Xrep_x*L_XC_0 + Xrep_y*L_XC_1 + L_XC_2) / ceres::sqrt((L_XC_0*L_XC_0 + L_XC_1*L_XC_1));

    out_residuals[0] = d_Xrep_L;

    //
    //
//    const T n = ref_N[0];
    const T& n = cam_intrinsics[3];
    const T * Za = (T*) water_plane;

    T nl1[3];
    nl1[0] = R[0]*x1_c[0] + R[1]*x1_c[1] + R[2];
    nl1[1] = R[3]*x1_c[0] + R[4]*x1_c[1] + R[5];
    nl1[2] = R[6]*x1_c[0] + R[7]*x1_c[1] + R[8];

    T tan_a1 = ceres::tan(ceres::acos( nl1[2] / ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1] + nl1[2]*nl1[2])));
    T a1 = ceres::acos((nl1[2]) / ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1] + nl1[2]*nl1[2]));
    if (a1 > ceres::acos(-1)/2)
    {
        a1 = ceres::acos(-1) - a1;
        tan_a1 = ceres::tan(a1);
    }
    const T tan_a2 = ceres::tan(ceres::asin(ceres::sin(a1)/n));

    const T Zp = pos_3dpoint[2];
    const T r_in = tan_a1 * (Za[0]-cam_C[2]);

    const T r_xy = ceres::sqrt((pos_3dpoint[0]-cam_C[0])*(pos_3dpoint[0]-cam_C[0]) + (pos_3dpoint[1]-cam_C[1])*(pos_3dpoint[1]-cam_C[1]) );
    const T Xq[3] = {(Za[0]-cam_C[2])/nl1[2]*nl1[0]+cam_C[0], (Za[0]-cam_C[2])/nl1[2]*nl1[1]+cam_C[1], Za[0]};
    const T l_ref[3] = {nl1[0]/ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1]), nl1[1]/ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1]), (T)1.0/tan_a2*(ceres::sin(a1)/n)*nl1[2]/ceres::abs(nl1[2])};

//    if (Zp <= Za[0])
    if (r_in >= r_xy)
    {
        const T X_subtract_C[3] = {pos_3dpoint[0]-cam_C[0], pos_3dpoint[1]-cam_C[1], pos_3dpoint[2]-cam_C[2]};
        const T crossXC[3] = {X_subtract_C[1]*nl1[2]-X_subtract_C[2]*nl1[1], X_subtract_C[2]*nl1[0]-X_subtract_C[0]*nl1[2], X_subtract_C[0]*nl1[1]-X_subtract_C[1]*nl1[0]};
        const T D_in = ceres::sqrt(crossXC[0]*crossXC[0] + crossXC[1]*crossXC[1] + crossXC[2]*crossXC[2]) / ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1] + nl1[2]*nl1[2]);

        //out_residuals[1] = r_xy-(r_in-tan_a1*(Za[0]-Zp));
        out_residuals[1] = r_in - r_xy-tan_a1*(Za[0]-Zp);

    }else{
        const T n_ref[3] = {nl1[0]/ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1]), nl1[1]/ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1]), (T)1.0/tan_a2*(ceres::sin(a1)/n)*nl1[2]/ceres::abs(nl1[2])};
        const T X_subtract_C[3] = {pos_3dpoint[0]-Xq[0], pos_3dpoint[1]-Xq[1], pos_3dpoint[2]-Xq[2]};
        const T crossXC[3] = {X_subtract_C[1]*n_ref[2]-X_subtract_C[2]*n_ref[1], X_subtract_C[2]*n_ref[0]-X_subtract_C[0]*n_ref[2], X_subtract_C[0]*n_ref[1]-X_subtract_C[1]*n_ref[0]};
        const T D_ref = ceres::sqrt(crossXC[0]*crossXC[0] + crossXC[1]*crossXC[1] + crossXC[2]*crossXC[2]) / ceres::sqrt(n_ref[0]*n_ref[0] + n_ref[1]*n_ref[1] + n_ref[2]*n_ref[2]);
        const T r_ref = tan_a2 * (Zp-Za[0]);

        out_residuals[1] = r_xy-r_in-r_ref;

    }

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double & water_plane,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_WaterZN, 2, 4, 6, 3>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_WaterZN(observation.data(), &water_plane)));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_WaterZN>, 2, 4, 6, 3>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_WaterZN>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_WaterZN(observation.data(), &water_plane), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
  const double * water_plane;
};

struct ResidualErrorFunctor_Pinhole_Intrinsic_Water_minz
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Water_minz(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum : uint8_t {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2
  };

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y] )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    //const T* const ref_N,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_extrinsics;
    const T * cam_C = &cam_extrinsics[3];

    //prepare l
    T R[9];
    ceres::AngleAxisToRotationMatrix(cam_R, R);

    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];

    T x1_c[3];
    x1_c[0] = (m_pos_2dpoint[0] - principal_point_x) /focal;
    x1_c[1] = (m_pos_2dpoint[1] - principal_point_y) /focal;
    x1_c[2] = (T) 1.0;

    const T Xrep_x = R[0]*x1_c[0] + R[1]*x1_c[1] + R[2] + cam_C[0];
    const T Xrep_y = R[3]*x1_c[0] + R[4]*x1_c[1] + R[5] + cam_C[1];
    //const T Xrep_z = R[6]*x1_c0 + R[7]*x1_c1 + R[8] + cam_C[2];


    const T L_XC_0 = cam_C[1] - pos_3dpoint[1];
    const T L_XC_1 = pos_3dpoint[0] - cam_C[0];
    const T L_XC_2 = cam_C[0]*pos_3dpoint[1] - cam_C[1]*pos_3dpoint[0];

    //const T d_Xrep_L = ceres::abs(Xrep_x*L_XC_0 + Xrep_y*L_XC_1 + L_XC_2);// / ceres::sqrt((L_XC_0*L_XC_0 + L_XC_1*L_XC_1));
    const T d_Xrep_L = ceres::abs(Xrep_x*L_XC_0 + Xrep_y*L_XC_1 + L_XC_2) / ceres::sqrt((L_XC_0*L_XC_0 + L_XC_1*L_XC_1));

    out_residuals[0] = d_Xrep_L;

    //
    //
//    const T n = ref_N[0];
    const T & water_plane = cam_intrinsics[3];
    const T n = (T)1.33;
    const T *Za =  &water_plane;

    T nl1[3];
    nl1[0] = R[0]*x1_c[0] + R[1]*x1_c[1] + R[2];
    nl1[1] = R[3]*x1_c[0] + R[4]*x1_c[1] + R[5];
    nl1[2] = R[6]*x1_c[0] + R[7]*x1_c[1] + R[8];

    T tan_a1 = ceres::tan(ceres::acos( nl1[2] / ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1] + nl1[2]*nl1[2])));
    T a1 = ceres::acos((nl1[2]) / ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1] + nl1[2]*nl1[2]));
    if (a1 > ceres::acos(-1)/2)
    {
        a1 = ceres::acos(-1) - a1;
        tan_a1 = ceres::tan(a1);
    }
    const T tan_a2 = ceres::tan(ceres::asin(ceres::sin(a1)/n));

    const T Zp = pos_3dpoint[2];
    const T r_in = tan_a1 * (Za[0]-cam_C[2]);

    const T r_xy = ceres::sqrt((pos_3dpoint[0]-cam_C[0])*(pos_3dpoint[0]-cam_C[0]) + (pos_3dpoint[1]-cam_C[1])*(pos_3dpoint[1]-cam_C[1]) );
    const T Xq[3] = {(Za[0]-cam_C[2])/nl1[2]*nl1[0]+cam_C[0], (Za[0]-cam_C[2])/nl1[2]*nl1[1]+cam_C[1], Za[0]};
    const T l_ref[3] = {nl1[0]/ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1]), nl1[1]/ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1]), (T)1.0/tan_a2*(ceres::sin(a1)/n)*nl1[2]/ceres::abs(nl1[2])};

//    if (Zp <= Za[0])
    if (r_in >= r_xy)
    {
        const T X_subtract_C[3] = {pos_3dpoint[0]-cam_C[0], pos_3dpoint[1]-cam_C[1], pos_3dpoint[2]-cam_C[2]};
        const T crossXC[3] = {X_subtract_C[1]*nl1[2]-X_subtract_C[2]*nl1[1], X_subtract_C[2]*nl1[0]-X_subtract_C[0]*nl1[2], X_subtract_C[0]*nl1[1]-X_subtract_C[1]*nl1[0]};
        const T D_in = ceres::sqrt(crossXC[0]*crossXC[0] + crossXC[1]*crossXC[1] + crossXC[2]*crossXC[2]) / ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1] + nl1[2]*nl1[2]);

        //out_residuals[1] = r_xy-(r_in-tan_a1*(Za[0]-Zp));
        out_residuals[1] = r_in - r_xy-tan_a1*(Za[0]-Zp);

    }else{
        const T n_ref[3] = {nl1[0]/ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1]), nl1[1]/ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1]), (T)1.0/tan_a2*(ceres::sin(a1)/n)*nl1[2]/ceres::abs(nl1[2])};
        const T X_subtract_C[3] = {pos_3dpoint[0]-Xq[0], pos_3dpoint[1]-Xq[1], pos_3dpoint[2]-Xq[2]};
        const T crossXC[3] = {X_subtract_C[1]*n_ref[2]-X_subtract_C[2]*n_ref[1], X_subtract_C[2]*n_ref[0]-X_subtract_C[0]*n_ref[2], X_subtract_C[0]*n_ref[1]-X_subtract_C[1]*n_ref[0]};
        const T D_ref = ceres::sqrt(crossXC[0]*crossXC[0] + crossXC[1]*crossXC[1] + crossXC[2]*crossXC[2]) / ceres::sqrt(n_ref[0]*n_ref[0] + n_ref[1]*n_ref[1] + n_ref[2]*n_ref[2]);
        const T r_ref = tan_a2 * (Zp-Za[0]);

        out_residuals[1] = r_xy-r_in-r_ref;

    }

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_Water_minz, 2, 4, 6, 3>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Water_minz(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Water_minz>, 2, 4, 6, 3>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Water_minz>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Water_minz(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};





//struct ResidualErrorFunctor_Pinhole_Intrinsic_WaterXY
//{
//  ResidualErrorFunctor_Pinhole_Intrinsic_WaterXY(const double* const pos_2dpoint)
//  :m_pos_2dpoint(pos_2dpoint)
//  {
//  }

//  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
//  enum : uint8_t {
//    OFFSET_FOCAL_LENGTH = 0,
//    OFFSET_PRINCIPAL_POINT_X = 1,
//    OFFSET_PRINCIPAL_POINT_Y = 2
//  };

//  /**
//   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y] )
//   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
//   *   - 3 for rotation(angle axis), 3 for translation
//   * @param[in] pos_3dpoint
//   * @param[out] out_residuals
//   */
//  template <typename T>
//  bool operator()(
//    const T* const cam_intrinsics,
//    const T* const cam_extrinsics,
//    const T* const pos_3dpoint,
//    T* out_residuals) const
//  {
//    //--
//    // Apply external parameters (Pose)
//    //--

//    const T * cam_R = cam_extrinsics;
//    const T * cam_C = &cam_extrinsics[3];

//    //prepare l
//    T R[9];
//    ceres::AngleAxisToRotationMatrix(cam_R, R);

//    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
//    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
//    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];

//    T x1_c[3];
//    x1_c[0] = (m_pos_2dpoint[0] - principal_point_x) /focal;
//    x1_c[1] = (m_pos_2dpoint[1] - principal_point_y) /focal;
//    x1_c[2] = (T) 1.0;

//    const T Xrep_x = R[0]*x1_c[0] + R[1]*x1_c[1] + R[2] + cam_C[0];
//    const T Xrep_y = R[3]*x1_c[0] + R[4]*x1_c[1] + R[5] + cam_C[1];
//    //const T Xrep_z = R[6]*x1_c0 + R[7]*x1_c1 + R[8] + cam_C[2];


//    const T L_XC_0 = cam_C[1] - pos_3dpoint[1];
//    const T L_XC_1 = pos_3dpoint[0] - cam_C[0];
//    const T L_XC_2 = cam_C[0]*pos_3dpoint[1] - cam_C[1]*pos_3dpoint[0];

//    //const T d_Xrep_L = ceres::abs(Xrep_x*L_XC_0 + Xrep_y*L_XC_1 + L_XC_2);// / ceres::sqrt((L_XC_0*L_XC_0 + L_XC_1*L_XC_1));
//    const T d_Xrep_L = ceres::abs(Xrep_x*L_XC_0 + Xrep_y*L_XC_1 + L_XC_2) / ceres::sqrt((L_XC_0*L_XC_0 + L_XC_1*L_XC_1));

//    out_residuals[0] = d_Xrep_L;

//    return true;
//  }

//  static int num_residuals() { return 1; }

//  // Factory to hide the construction of the CostFunction object from
//  // the client code.
//  static ceres::CostFunction* Create
//  (
//    const Vec2 & observation,
//    const double weight = 0.0
//  )
//  {
//    if (weight == 0.0)
//    {
//      return
//        (new ceres::AutoDiffCostFunction
//          <ResidualErrorFunctor_Pinhole_Intrinsic_WaterXY, 1, 3, 6, 2>(
//            new ResidualErrorFunctor_Pinhole_Intrinsic_WaterXY(observation.data())));
//    }
//    else
//    {
//      return
//        (new ceres::AutoDiffCostFunction
//          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_WaterXY>, 1, 3, 6, 2>
//          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_WaterXY>
//            (new ResidualErrorFunctor_Pinhole_Intrinsic_WaterXY(observation.data()), weight)));
//    }
//  }

//  const double * m_pos_2dpoint; // The 2D observation
//};



struct ResidualErrorFunctor_Pinhole_Intrinsic_water_prepare
{
  ResidualErrorFunctor_Pinhole_Intrinsic_water_prepare(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum : uint8_t {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2
  };

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y] )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];

    T pos_proj[3];
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_proj);

    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    // Transform the point from homogeneous to euclidean (undistorted point)
    const T x_u = pos_proj[0] / pos_proj[2];
    const T y_u = pos_proj[1] / pos_proj[2];

    //--
    // Apply intrinsic parameters
    //--

    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal * x_u;
    const T projected_y = principal_point_y + focal * y_u;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - m_pos_2dpoint[0];
    out_residuals[1] = projected_y - m_pos_2dpoint[1];

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::NumericDiffCostFunction//ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_water_prepare, ceres::RIDDERS,2, 3, 6, 3>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_water_prepare(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic>, 2, 3, 6, 3>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic>
            (new ResidualErrorFunctor_Pinhole_Intrinsic(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};


/**
 * @brief Ceres functor to use a Pinhole_Intrinsic_Radial_K1
 *
 *  Data parameter blocks are the following <2,4,6,3>
 *  - 2 => dimension of the residuals,
 *  - 4 => the intrinsic data block [focal, principal point x, principal point y, K1],
 *  - 6 => the camera extrinsic data block (camera orientation and position) [R;t],
 *         - rotation(angle axis), and translation [rX,rY,rZ,tx,ty,tz].
 *  - 3 => a 3D point data block.
 *
 */
struct ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum : uint8_t {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2,
    OFFSET_DISTO_K1 = 3
  };

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], K1 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];

    T pos_proj[3];
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_proj);

    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    // Transform the point from homogeneous to euclidean (undistorted point)
    const T x_u = pos_proj[0] / pos_proj[2];
    const T y_u = pos_proj[1] / pos_proj[2];

    //--
    // Apply intrinsic parameters
    //--

    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];
    const T& k1 = cam_intrinsics[OFFSET_DISTO_K1];

    // Apply distortion (xd,yd) = disto(x_u,y_u)
    const T r2 = x_u*x_u + y_u*y_u;
    const T r_coeff = 1.0 + k1*r2;
    const T x_d = x_u * r_coeff;
    const T y_d = y_u * r_coeff;

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal * x_d;
    const T projected_y = principal_point_y + focal * y_d;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - m_pos_2dpoint[0];
    out_residuals[1] = projected_y - m_pos_2dpoint[1];

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1, 2, 4, 6, 3>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1>, 2, 4, 6, 3>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};

/**
 * @brief Ceres functor to use a Pinhole_Intrinsic_Radial_K3
 *
 *  Data parameter blocks are the following <2,6,6,3>
 *  - 2 => dimension of the residuals,
 *  - 6 => the intrinsic data block [focal, principal point x, principal point y, K1, K2, K3],
 *  - 6 => the camera extrinsic data block (camera orientation and position) [R;t],
 *         - rotation(angle axis), and translation [rX,rY,rZ,tx,ty,tz].
 *  - 3 => a 3D point data block.
 *
 */
struct ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum : uint8_t {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2,
    OFFSET_DISTO_K1 = 3,
    OFFSET_DISTO_K2 = 4,
    OFFSET_DISTO_K3 = 5,
  };

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--
    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];

    T pos_proj[3];
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_proj);

    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    // Transform the point from homogeneous to euclidean (undistorted point)
    const T x_u = pos_proj[0] / pos_proj[2];
    const T y_u = pos_proj[1] / pos_proj[2];

    //--
    // Apply intrinsic parameters
    //--

    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];
    const T& k1 = cam_intrinsics[OFFSET_DISTO_K1];
    const T& k2 = cam_intrinsics[OFFSET_DISTO_K2];
    const T& k3 = cam_intrinsics[OFFSET_DISTO_K3];

    // Apply distortion (xd,yd) = disto(x_u,y_u)
    const T r2 = x_u*x_u + y_u*y_u;
    const T r4 = r2 * r2;
    const T r6 = r4 * r2;
    const T r_coeff = (1.0 + k1*r2 + k2*r4 + k3*r6);
    const T x_d = x_u * r_coeff;
    const T y_d = y_u * r_coeff;

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal * x_d;
    const T projected_y = principal_point_y + focal * y_d;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - m_pos_2dpoint[0];
    out_residuals[1] = projected_y - m_pos_2dpoint[1];

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3, 2, 6, 6, 3>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3>, 2, 6, 6, 3>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};

/**
 * @brief Ceres functor with constrained 3D points to use a Pinhole_Intrinsic_Brown_T2
 *
 *  Data parameter blocks are the following <2,8,6,3>
 *  - 2 => dimension of the residuals,
 *  - 8 => the intrinsic data block [focal, principal point x, principal point y, K1, K2, K3, T1, T2],
 *  - 6 => the camera extrinsic data block (camera orientation and position) [R;t],
 *         - rotation(angle axis), and translation [rX,rY,rZ,tx,ty,tz].
 *  - 3 => a 3D point data block.
 *
 */
struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint),useFraction_(true)
  {
  }
    const bool useFraction_;//是否考虑海面折射
  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
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

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];
    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];
    T pos_proj[3];
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_proj);

    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    // Transform the point from homogeneous to euclidean (undistorted point)
    const T x_u = pos_proj[0] / pos_proj[2];
    const T y_u = pos_proj[1] / pos_proj[2];

    //--
    // Apply intrinsic parameters
    //--


    const T& k1 = cam_intrinsics[OFFSET_DISTO_K1];
    const T& k2 = cam_intrinsics[OFFSET_DISTO_K2];
    const T& k3 = cam_intrinsics[OFFSET_DISTO_K3];
    const T& t1 = cam_intrinsics[OFFSET_DISTO_T1];
    const T& t2 = cam_intrinsics[OFFSET_DISTO_T2];

    // Apply distortion (xd,yd) = disto(x_u,y_u)
    const T r2 = x_u*x_u + y_u*y_u;
    const T r4 = r2 * r2;
    const T r6 = r4 * r2;
    const T r_coeff = (1.0 + k1*r2 + k2*r4 + k3*r6);
    const T t_x = t2 * (r2 + 2.0 * x_u*x_u) + 2.0 * t1 * x_u * y_u;
    const T t_y = t1 * (r2 + 2.0 * y_u*y_u) + 2.0 * t2 * x_u * y_u;
    const T x_d = x_u * r_coeff + t_x;
    const T y_d = y_u * r_coeff + t_y;

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal * x_d;
    const T projected_y = principal_point_y + focal * y_d;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - m_pos_2dpoint[0];
    out_residuals[1] = projected_y - m_pos_2dpoint[1];
    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2, 2, 8, 6, 3>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2>, 2, 8, 6, 3>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};

struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_prepare
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_prepare(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
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

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];

    T pos_proj[3];
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_proj);

    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    // Transform the point from homogeneous to euclidean (undistorted point)
    const T x_u = pos_proj[0] / pos_proj[2];
    const T y_u = pos_proj[1] / pos_proj[2];

    //--
    // Apply intrinsic parameters
    //--

    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];
    const T& k1 = cam_intrinsics[OFFSET_DISTO_K1];
    const T& k2 = cam_intrinsics[OFFSET_DISTO_K2];
    const T& k3 = cam_intrinsics[OFFSET_DISTO_K3];
    const T& t1 = cam_intrinsics[OFFSET_DISTO_T1];
    const T& t2 = cam_intrinsics[OFFSET_DISTO_T2];

    // Apply distortion (xd,yd) = disto(x_u,y_u)
    const T r2 = x_u*x_u + y_u*y_u;
    const T r4 = r2 * r2;
    const T r6 = r4 * r2;
    const T r_coeff = (1.0 + k1*r2 + k2*r4 + k3*r6);
    const T t_x = t2 * (r2 + 2.0 * x_u*x_u) + 2.0 * t1 * x_u * y_u;
    const T t_y = t1 * (r2 + 2.0 * y_u*y_u) + 2.0 * t2 * x_u * y_u;
    const T x_d = x_u * r_coeff + t_x;
    const T y_d = y_u * r_coeff + t_y;

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal * x_d;
    const T projected_y = principal_point_y + focal * y_d;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - m_pos_2dpoint[0];
    out_residuals[1] = projected_y - m_pos_2dpoint[1];

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::NumericDiffCostFunction//ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_prepare, ceres::RIDDERS, 2, 8, 6, 3>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_prepare(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2>, 2, 8, 6, 3>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};



///
struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_ThreeCamera_IMU_GPS_trans_change
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_ThreeCamera_IMU_GPS_trans_change(const double* const pos_2dpoint)
  {
    m_pos_2dpoint[0] = pos_2dpoint[0];
    m_pos_2dpoint[1] = pos_2dpoint[1];
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2,
    OFFSET_DISTO_K1 = 3,
    OFFSET_DISTO_K2 = 4,
    OFFSET_DISTO_K3 = 5,
    OFFSET_DISTO_T1 = 6,
    OFFSET_DISTO_T2 = 7,
  };

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 ) 8 parameters for three camera
   * @param[in] cam_R : R_x 3 paramerters
   * @param[in] transformation_R: Rotation of coordinate system ,use to correct world to camera rotation (imu)
   *                              3 paramerters for rotation(angle axis)
   * @param[in] transformation_x: Translation of coordinate system ,use to correct world to camera rotation (gps)
   * @param[in] gps_x: data from gps, keep content
   * @param[in] reference : R_x t_x 6 paramerters
   * @param[in] pos_3dpoint for each of three camera
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_R,
    const T* const gps_X,
    const T* const transformation_R,
    const T* const transformation_x,
    const T* const cam_transformation,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {

    //--
    // Apply external parameters (Pose)
    //--
    T point_gpsX_camR[3];
    ceres::AngleAxisRotatePoint(cam_R, gps_X, point_gpsX_camR);
    T point_gpsX_imuR[3];
    ceres::AngleAxisRotatePoint(transformation_R, point_gpsX_camR, point_gpsX_imuR);

    T pos_camR[3];
    ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_camR);
    T pos_imuR[3];
    ceres::AngleAxisRotatePoint(transformation_R, pos_camR, pos_imuR);

    T pos_proj_temp[3];
    pos_proj_temp[0] = pos_imuR[0] + transformation_x[0] - point_gpsX_imuR[0];
    pos_proj_temp[1] = pos_imuR[1] + transformation_x[1] - point_gpsX_imuR[1];
    pos_proj_temp[2] = pos_imuR[2] + transformation_x[2] - point_gpsX_imuR[2];

    const T *trans_R = cam_transformation;
    const T *trans_t = &cam_transformation[3];

    T pos_proj[3];
    ceres::AngleAxisRotatePoint(trans_R, pos_proj_temp, pos_proj);
    pos_proj[0] += trans_t[0];
    pos_proj[1] += trans_t[1];
    pos_proj[2] += trans_t[2];




    // Transform the point from homogeneous to euclidean (undistorted point)
    const T x_u = pos_proj[0] / pos_proj[2];
    const T y_u = pos_proj[1] / pos_proj[2];

    //--
    // Apply intrinsic parameters
    //--

    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];
    const T& k1 = cam_intrinsics[OFFSET_DISTO_K1];
    const T& k2 = cam_intrinsics[OFFSET_DISTO_K2];
    const T& k3 = cam_intrinsics[OFFSET_DISTO_K3];
    const T& t1 = cam_intrinsics[OFFSET_DISTO_T1];
    const T& t2 = cam_intrinsics[OFFSET_DISTO_T2];

    // Apply distortion (xd,yd) = disto(x_u,y_u)
    const T r2 = x_u*x_u + y_u*y_u;
    const T r4 = r2 * r2;
    const T r6 = r4 * r2;
    const T r_coeff = (T(1) + k1*r2 + k2*r4 + k3*r6);
    const T t_x = t2 * (r2 + T(2) * x_u*x_u) + T(2) * t1 * x_u * y_u;
    const T t_y = t1 * (r2 + T(2) * y_u*y_u) + T(2) * t2 * x_u * y_u;
    const T x_d = x_u * r_coeff + t_x;
    const T y_d = y_u * r_coeff + t_y;

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal * x_d;
    const T projected_y = principal_point_y + focal * y_d;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - T(m_pos_2dpoint[0]);
    out_residuals[1] = projected_y - T(m_pos_2dpoint[1]);

    return true;
  }

  static const int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
      if(weight == 0.0)
      {
          return
              (new ceres::AutoDiffCostFunction
                <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_ThreeCamera_IMU_GPS_trans_change, 2, 8, 3, 3, 3, 3, 6, 3>(
                  new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_ThreeCamera_IMU_GPS_trans_change(observation.data())));
      }else{
          return
              (new ceres::AutoDiffCostFunction
                <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_ThreeCamera_IMU_GPS_trans_change, 2, 8, 3, 3, 3, 3, 6, 3>(
                  new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_ThreeCamera_IMU_GPS_trans_change(observation.data()), weight));
      }

  }

  double m_pos_2dpoint[2]; // The 2D observation
};

///c1 gps imu
struct ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE
{
  ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(const double* const pos_2dpoint)
  {
    m_pos_2dpoint[0] = pos_2dpoint[0];
    m_pos_2dpoint[1] = pos_2dpoint[1];
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2
  };



  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y]) 3 parameters for three camera
   * @param[in] cam_R : R_x 3 paramerters
   * @param[in] gps_x: data from gps, keep content
   * @param[in] transformation_x: Translation of coordinate system ,use to correct world to camera rotation (gps)
   * @param[in] transformation_R: Rotation of coordinate system ,use to correct world to camera rotation (imu)
   *                              3 paramerters for rotation(angle axis)
   * @param[in] reference : R_x t_x 6 paramerters
   * @param[in] pos_3dpoint for each of three camera
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_R,
    const T* const gps_X,
    const T* const transformation_R,
    const T* const transformation_x,
    const T* const cam_transformation,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {

      //--
      // Apply external parameters (Pose)
      //--
      T point_gpsX_camR[3];
      ceres::AngleAxisRotatePoint(cam_R, gps_X, point_gpsX_camR);
      T point_gpsX_imuR[3];
      ceres::AngleAxisRotatePoint(transformation_R, point_gpsX_camR, point_gpsX_imuR);

      T pos_camR[3];
      ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_camR);
      T pos_imuR[3];
      ceres::AngleAxisRotatePoint(transformation_R, pos_camR, pos_imuR);

      T pos_proj_temp[3];
      pos_proj_temp[0] = pos_imuR[0] + transformation_x[0] - point_gpsX_imuR[0];
      pos_proj_temp[1] = pos_imuR[1] + transformation_x[1] - point_gpsX_imuR[1];
      pos_proj_temp[2] = pos_imuR[2] + transformation_x[2] - point_gpsX_imuR[2];

      const T *trans_R = cam_transformation;
      const T *trans_t = &cam_transformation[3];

      T pos_proj[3];
      ceres::AngleAxisRotatePoint(trans_R, pos_proj_temp, pos_proj);
      pos_proj[0] += trans_t[0];
      pos_proj[1] += trans_t[1];
      pos_proj[2] += trans_t[2];

    // Transform the point from homogeneous to euclidean (undistorted point)
    const T x_u = pos_proj[0] / pos_proj[2];
    const T y_u = pos_proj[1] / pos_proj[2];

    //--
    // Apply intrinsic parameters
    //--

    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal * x_u;
    const T projected_y = principal_point_y + focal * y_u;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - T(m_pos_2dpoint[0]);
    out_residuals[1] = projected_y - T(m_pos_2dpoint[1]);

    return true;
  }

  static const int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
//    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(observation.data())));
    }

  }

  double m_pos_2dpoint[2]; // The 2D observation
};

///c4 gps imu
struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_ThreeCamera_IMU_GPS
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_ThreeCamera_IMU_GPS(const double* const pos_2dpoint)
  {
    m_pos_2dpoint[0] = pos_2dpoint[0];
    m_pos_2dpoint[1] = pos_2dpoint[1];
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2,
    OFFSET_DISTO_K1 = 3,
    OFFSET_DISTO_K2 = 4,
    OFFSET_DISTO_K3 = 5,
    OFFSET_DISTO_T1 = 6,
    OFFSET_DISTO_T2 = 7,
  };

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 ) 8 parameters for three camera
   * @param[in] cam_R : R_x 3 paramerters
   * @param[in] transformation_R: Rotation of coordinate system ,use to correct world to camera rotation (imu)
   *                              3 paramerters for rotation(angle axis)
   * @param[in] transformation_x: Translation of coordinate system ,use to correct world to camera rotation (gps)
   * @param[in] gps_x: data from gps, keep content
   * @param[in] reference : R_x t_x 6 paramerters
   * @param[in] pos_3dpoint for each of three camera
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_R,
    const T* const gps_X,
    const T* const transformation_R,
    const T* const transformation_x,
    const T* const cam_transformation,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {

    //--
    // Apply external parameters (Pose)
    //--
    T point_gpsX_camR[3];
    ceres::AngleAxisRotatePoint(cam_R, gps_X, point_gpsX_camR);
    T point_gpsX_imuR[3];
    ceres::AngleAxisRotatePoint(transformation_R, point_gpsX_camR, point_gpsX_imuR);

    T pos_camR[3];
    ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_camR);
    T pos_imuR[3];
    ceres::AngleAxisRotatePoint(transformation_R, pos_camR, pos_imuR);

    T pos_proj_temp[3];
    pos_proj_temp[0] = pos_imuR[0] + transformation_x[0] - point_gpsX_imuR[0];
    pos_proj_temp[1] = pos_imuR[1] + transformation_x[1] - point_gpsX_imuR[1];
    pos_proj_temp[2] = pos_imuR[2] + transformation_x[2] - point_gpsX_imuR[2];

    const T *trans_R = cam_transformation;
    const T *trans_t = &cam_transformation[3];

    T pos_proj[3];
    ceres::AngleAxisRotatePoint(trans_R, pos_proj_temp, pos_proj);
    pos_proj[0] += trans_t[0];
    pos_proj[1] += trans_t[1];
    pos_proj[2] += trans_t[2];




    // Transform the point from homogeneous to euclidean (undistorted point)
    const T x_u = pos_proj[0] / pos_proj[2];
    const T y_u = pos_proj[1] / pos_proj[2];

    //--
    // Apply intrinsic parameters
    //--

    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];
    const T& k1 = cam_intrinsics[OFFSET_DISTO_K1];
    const T& k2 = cam_intrinsics[OFFSET_DISTO_K2];
    const T& k3 = cam_intrinsics[OFFSET_DISTO_K3];
    const T& t1 = cam_intrinsics[OFFSET_DISTO_T1];
    const T& t2 = cam_intrinsics[OFFSET_DISTO_T2];

    // Apply distortion (xd,yd) = disto(x_u,y_u)
    const T r2 = x_u*x_u + y_u*y_u;
    const T r4 = r2 * r2;
    const T r6 = r4 * r2;
    const T r_coeff = (T(1) + k1*r2 + k2*r4 + k3*r6);
    const T t_x = t2 * (r2 + T(2) * x_u*x_u) + T(2) * t1 * x_u * y_u;
    const T t_y = t1 * (r2 + T(2) * y_u*y_u) + T(2) * t2 * x_u * y_u;
    const T x_d = x_u * r_coeff + t_x;
    const T y_d = y_u * r_coeff + t_y;

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal * x_d;
    const T projected_y = principal_point_y + focal * y_d;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - T(m_pos_2dpoint[0]);
    out_residuals[1] = projected_y - T(m_pos_2dpoint[1]);

    return true;
  }

  static const int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
      if(weight == 0.0)
      {
          return
              (new ceres::AutoDiffCostFunction
                <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_ThreeCamera_IMU_GPS, 2, 8, 3, 3, 3, 3, 6, 3>(
                  new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_ThreeCamera_IMU_GPS(observation.data())));
      }else{
          return
              (new ceres::AutoDiffCostFunction
                <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_ThreeCamera_IMU_GPS, 2, 8, 3, 3, 3, 3, 6, 3>(
                  new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_ThreeCamera_IMU_GPS(observation.data()), weight));
      }

  }

  double m_pos_2dpoint[2]; // The 2D observation
};

///c1 water
struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
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

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--
//#pragma omp critical
//      {
    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];


    T pos_proj[3];
    T X[3];
    X[0] = pos_3dpoint[0];
    X[1] = pos_3dpoint[1];
    X[2] = T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);
    // Rotate the point according the camera rotation
//    ceres::AngleAxisRotatePoint(cam_R, X, pos_proj);

//    //prepare R
    const T theta = sqrt(cam_R[0]*cam_R[0] + cam_R[1]*cam_R[1] + cam_R[2]*cam_R[2]);
    const T cosTheta = cos(theta);
    const T sinTheta = sin(theta);
    const T r[3] = {cam_R[0]/theta, cam_R[1]/theta, cam_R[2]/theta};
    const T one_cosTheta = T(1.0) - cosTheta;

    T R[9];
    R[0] = cosTheta + one_cosTheta*r[0]*r[0];
    R[1] = one_cosTheta*r[0]*r[1] + sinTheta*r[2];
    R[2] = one_cosTheta*r[0]*r[2] - sinTheta*r[1];
    R[3] = one_cosTheta*r[0]*r[1] - sinTheta*r[2];
    R[4] = cosTheta + one_cosTheta*r[1]*r[1];
    R[5] = one_cosTheta*r[1]*r[2] + sinTheta*r[0];
    R[6] = one_cosTheta*r[0]*r[2] + sinTheta*r[1];
    R[7] = one_cosTheta*r[1]*r[2] - sinTheta*r[0];
    R[8] = cosTheta + one_cosTheta*r[2]*r[2];

    const T& f = cam_intrinsics[OFFSET_FOCAL_LENGTH]; //focal
    const T& u = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
    const T& v = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y
    //
    const T P1 = f*R[0] + u*R[2];
    const T P2 = f*R[3] + u*R[5];
    const T P3 = f*R[6] + u*R[8];
    const T P4 = f*cam_t[0] + u*cam_t[2];
    const T P5 = f*R[1] + v*R[2];
    const T P6 = f*R[4] + v*R[5];
    const T P7 = f*R[7] + v*R[8];
    const T P8 = f*cam_t[1] + v*cam_t[2];
    const T P9 = R[2];
    const T P10 = R[5];
    const T P11 = R[8];
    const T P12 = cam_t[2];
    //
    const T x_i = (P1*X[0] + P2*X[1] + P3*X[2] + P4) / (P9*X[0] + P10*X[1] + P11*X[2] + P12);
    const T y_i = (P5*X[0] + P6*X[1] + P7*X[2] + P8) / (P9*X[0] + P10*X[1] + P11*X[2] + P12);
    //
    const T xu = T(m_pos_2dpoint[0]);
    const T yu = T(m_pos_2dpoint[1]);
    //
    const T nc_x = P3/P11;
    const T nc_y = P7/P11;
    //
//    const T l_ncxi_x = x_i - nc_x;
//    const T l_ncxi_y = y_i - nc_y;
//    const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
//    const T normal_ncxi = (l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
    //
    const T l_ncxu_x = xu - nc_x;
    const T l_ncxu_y = yu - nc_y;
    const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);
//    const T normal_ncxu = (l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);

    const T A = yu - nc_y;
    const T B = nc_x - xu;
    const T C = xu*nc_y - yu*nc_x;

    const T D = abs(A*x_i + B*y_i + C) / sqrt(A*A+B*B);
    ///
    {
        out_residuals[0] = D * normal_ncxu / T(2.0);
    }

    {
//        const T l_ncxi_x = x_i - nc_x;
//        const T l_ncxi_y = y_i - nc_y;
//        const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);

//        out_residuals[0] = abs(l_ncxi_x/normal_ncxi * l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi * l_ncxu_x/normal_ncxu) / T(2.0);

    }
    return true;
  }

  static int num_residuals() { return 1; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
          (new ceres::NumericDiffCostFunction//ceres::AutoDiffCostFunction//
            <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1, ceres::RIDDERS, 1, 3, 6, 2>(
//        (new ceres::AutoDiffCostFunction
//          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1, 1, 3, 6, 2>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1(observation.data())));
    }
    else
    {
      return
        (new ceres::NumericDiffCostFunction//ceres::AutoDiffCostFunction
//          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3, ceres::RIDDERS, 1, 3, 6, 2>(
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1>, ceres::RIDDERS, 1, 3, 6, 2>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};

///c1 water
struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_2_old
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_2_old(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
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

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--
//#pragma omp critical
//      {
//          std::cout << "error 1" << std::endl;
    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];


    T pos_proj[3];
    T X[3];
    X[0] = pos_3dpoint[0];
    X[1] = pos_3dpoint[1];
    X[2] = T(0.0);
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, X, pos_proj);

    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    //prepare l
    T R[9];
    ceres::AngleAxisToRotationMatrix(cam_R, R);

    {
//        const T theta = sqrt(cam_R[0]*cam_R[0] + cam_R[1]*cam_R[1] + cam_R[2]*cam_R[2]);

//        const T cosTheta = cos(theta);
//        const T sinTheta = sin(theta);
//        const T r0 = cam_R[0] / theta;
//        const T r1 = cam_R[1] / theta;
//        const T r2 = cam_R[2] / theta;

//        const T one_cosTheta = T(1.0) - cosTheta;
//        const T r0_2 = r0 * r0;
//        const T r1_2 = r1 * r1;
//        const T r2_2 = r2 * r2;

//        R[0] = cosTheta + one_cosTheta * r0_2;
//        R[1] = one_cosTheta * r0 * r1 + sinTheta * r2;
//        R[2] = one_cosTheta * r0 * r2 - sinTheta + r1;
//        R[3] = one_cosTheta * r0 * r1 - sinTheta * r2;
//        R[4] = cosTheta + one_cosTheta * r1_2;
//        R[5] = one_cosTheta * r1 * r2 + sinTheta * r0;
//        R[6] = one_cosTheta * r0 * r2 + sinTheta * r1;
//        R[7] = one_cosTheta * r1 * r2 - sinTheta * r0;
//        R[8] = cosTheta + one_cosTheta * r2_2;

    }


    const T& f = cam_intrinsics[OFFSET_FOCAL_LENGTH]; //focal
    const T& u = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
    const T& v = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y

    // Transform the point from homogeneous to euclidean (undistorted point)
//    const T x_i = u + f * (pos_proj[0] / pos_proj[2]);
//    const T y_i = v + f * (pos_proj[1] / pos_proj[2]);

    const T P1 = f*R[0] + u*R[2];
    const T P2 = f*R[3] + u*R[5];
    const T P3 = f*R[6] + u*R[8];
    const T P4 = f*cam_t[0] + u*cam_t[2];
    const T P5 = f*R[1] + v*R[2];
    const T P6 = f*R[4] + v*R[5];
    const T P7 = f*R[7] + v*R[8];
    const T P8 = f*cam_t[1] + v*cam_t[2];
    const T P9 = R[2];
    const T P10 = R[5];
    const T P11 = R[8];
    const T P12 = cam_t[2];

    const T x_i = (P1*pos_3dpoint[0] + P2*pos_3dpoint[1] + P4) / (P9*pos_3dpoint[0] + P10*pos_3dpoint[1] + P12);
    const T y_i = (P5*pos_3dpoint[0] + P6*pos_3dpoint[1] + P8) / (P9*pos_3dpoint[0] + P10*pos_3dpoint[1] + P12);



    const T xu = T(m_pos_2dpoint[0]);
    const T yu = T(m_pos_2dpoint[1]);

//    const T l0 = y_i*P11 - P7;//P7 - P11*yu;
//    const T l1 = P3 - P11*x_i;//P11*xu - P3;
//    const T l2 = P7*x_i - P3*y_i;//P3*yu - P7*xu;
//    out_residuals[0] = abs(xu*l0 + yu*l1 + l2)/sqrt(l0*l0 + l1*l1);
////    out_residuals[0] = xu*l0 + yu*l1 + l2;

//    const T l0_ = P7 - P11*yu;
//    const T l1_ = P11*xu - P3;
//    const T l2_ = P3*yu - P7*xu;
//    out_residuals[1] = abs(x_i*l0_ + y_i*l1_ + l2_)/sqrt(l0_*l0_ + l1_*l1_);


//    std::cout << "error 2" << std::endl;

    /// angle
//    out_residuals[0] = acos((l_ncxi_x*l_ncxu_x + l_ncxi_y*l_ncxu_y) / (sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y)*sqrt(l_ncxu_x*l_ncxu_x+l_ncxu_y*l_ncxu_y)));


    const T nc_x = P3/P11;
    const T nc_y = P7/P11;
    //
    const T l_ncxi_x = x_i - nc_x;
    const T l_ncxi_y = y_i - nc_y;
    const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
    //
    const T l_ncxu_x = xu - nc_x;
    const T l_ncxu_y = yu - nc_y;
    const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);
    ///
    out_residuals[0] = l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi;
    out_residuals[1] = l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi;


//    std::cout << "error 6" << std::endl;
//    {
//        std::cout << "cam_R : " << cam_R[0] << " " << cam_R[1] << " " << cam_R[2] << std::endl;
//        std::cout << "R : " << R[0] << " " << R[1] << " " << R[2]
//                            << R[3] << " " << R[4] << " " << R[5]
//                            << R[6] << " " << R[7] << " " << R[8] << std::endl;
//        std::cout << "cam_t : " << cam_t[0] << " " << cam_t[1] << " " << cam_t[2] << std::endl;
//        std::cout << "f : " << f << " u : " << u << " v : " << v << std::endl;
//        std::cout << "xu : " << xu << " yu : " << yu << std::endl;
//        std::cout << "pos_3dpoint[0]: " << pos_3dpoint[0] << " pos_3dpoint[1]: " << pos_3dpoint[1] << " pos_3dpoint[2]: " << pos_3dpoint[2] << std::endl;
//        std::cout << "nc_x : " << nc_x << " nc_y : " << nc_y << std::endl;
//        std::cout << "l_ncxi_x : " << l_ncxi_x << " l_ncxi_y : " << l_ncxi_y << std::endl;
//        std::cout << "l_ncxu_x : " << l_ncxu_x << " l_ncxu_y : " << l_ncxu_y << std::endl;
//        std::cout << "out_residuals[0]: "  << out_residuals[0] << std::endl;
//        std::cout << "****-----------------****" << std::endl;
//    }

//    std::cout << "error 7" << std::endl;
//      }

//    out_residuals[0] = sqrt((x_i-P3/P11)*(x_i-P3/P11) + (y_i-P7/P11)*(y_i-P7/P11));

//    out_residuals[1] = x_i*l0_ + y_i*l1_ + l2_;
//    if(l0 != T(0.0))
//    {
//        out_residuals[0] = (- y_i*l1 - l2)/l0 - x_i;

//    }else{
//        out_residuals[0] = T(0.0);
//    }
//    if(l1 != T(0.0))
//    {
//        out_residuals[1] = ( -x_i*l0 - l2)/l1 - y_i;
//    }else{
//        out_residuals[1] = T(0.0);

//    }
//    out_residuals[0] = ( x_i*l0 + y_i*l1 + l2);
//    out_residuals[0] = abs(x_i*l0 + y_i*l1 + l2)/sqrt(l0*l0 + l1*l1);
//    out_residuals[0] = abs(x_i*l0 + y_i*l1 + l2)/sqrt(l0*l0 + l1*l1);
//    out_residuals[1] = T(0.0);
//    out_residuals[0] = abs(x_i*l0 + y_i*l1 + l2)/sqrt(l0*l0 + l1*l1);

//    if(out_residuals[0] > T(0.1))
//    {
//        std::cout << "cam_R : " << cam_R[0] << " " << cam_R[1] << " " << cam_R[2] << std::endl;
//        std::cout << "R : " << R[0] << " " << R[1] << " " << R[2]
//                            << R[3] << " " << R[4] << " " << R[5]
//                            << R[6] << " " << R[7] << " " << R[8] << std::endl;
//        std::cout << "cam_t : " << cam_t[0] << " " << cam_t[1] << " " << cam_t[2] << std::endl;
//        std::cout << "f : " << f << " u : " << u << " v : " << v << std::endl;
//        std::cout << "xu : " << xu << " yu : " << yu << std::endl;
//        std::cout << "pos_3dpoint[0]: " << pos_3dpoint[0] << " pos_3dpoint[1]: " << pos_3dpoint[1] << " pos_3dpoint[2]: " << pos_3dpoint[2] << std::endl;
//        std::cout << "out_residuals[0]: "  << out_residuals[0] << std::endl;
//        std::cout << "****-----------------****" << std::endl;
//    }

//      }

//    T pos_proj[3];
//    // Rotate the point according the camera rotation
//    ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_proj);

//    // Apply the camera translation
//    pos_proj[0] += cam_t[0];
//    pos_proj[1] += cam_t[1];
//    pos_proj[2] += cam_t[2];

//    // Transform the point from homogeneous to euclidean (undistorted point)
//    const T x_u = pos_proj[0] / pos_proj[2];
//    const T y_u = pos_proj[1] / pos_proj[2];

//    //--
//    // Apply intrinsic parameters
//    //--


//    const T& k1 = cam_intrinsics[OFFSET_DISTO_K1];
//    const T& k2 = cam_intrinsics[OFFSET_DISTO_K2];
//    const T& k3 = cam_intrinsics[OFFSET_DISTO_K3];
//    const T& t1 = cam_intrinsics[OFFSET_DISTO_T1];
//    const T& t2 = cam_intrinsics[OFFSET_DISTO_T2];

//    // Apply distortion (xd,yd) = disto(x_u,y_u)
//    const T r2 = x_u*x_u + y_u*y_u;
//    const T r4 = r2 * r2;
//    const T r6 = r4 * r2;
//    const T r_coeff = (1.0 + k1*r2 + k2*r4 + k3*r6);
//    const T t_x = t2 * (r2 + 2.0 * x_u*x_u) + 2.0 * t1 * x_u * y_u;
//    const T t_y = t1 * (r2 + 2.0 * y_u*y_u) + 2.0 * t2 * x_u * y_u;
//    const T x_d = x_u * r_coeff + t_x;
//    const T y_d = y_u * r_coeff + t_y;

//    // Apply focal length and principal point to get the final image coordinates
//    const T projected_x = principal_point_x + focal * x_d;
//    const T projected_y = principal_point_y + focal * y_d;

//    // Compute and return the error is the difference between the predicted
//    //  and observed position
//    out_residuals[0] = projected_x - m_pos_2dpoint[0];
//    out_residuals[1] = projected_y - m_pos_2dpoint[1];

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_2_old, 2, 3, 6, 2>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_2_old(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_2_old>, 2, 3, 6, 2>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_2_old>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_2_old(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};

struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_1
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_1(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
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

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--
//#pragma omp critical
//      {
    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];


    T pos_proj[3];
    T X[3];
    X[0] = pos_3dpoint[0];
    X[1] = pos_3dpoint[1];
    X[2] = pos_3dpoint[2];//T(0.0);
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, X, pos_proj);

    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    //prepare l
    T R[9];
//    ceres::MatrixAdapter<T, 3, 3>(R);
    ceres::AngleAxisToRotationMatrix(cam_R, R);

    const T& f = cam_intrinsics[OFFSET_FOCAL_LENGTH]; //focal
    const T& u = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
    const T& v = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y

    // Transform the point from homogeneous to euclidean (undistorted point)
    const T x_i = u + f * (pos_proj[0] / pos_proj[2]);
    const T y_i = v + f * (pos_proj[1] / pos_proj[2]);

    const T P1 = f*R[0] + u*R[2];
    const T P2 = f*R[3] + u*R[5];
    const T P3 = f*R[6] + u*R[8];
    const T P4 = f*cam_t[0] + u*cam_t[2];
    const T P5 = f*R[1] + v*R[2];
    const T P6 = f*R[4] + v*R[5];
    const T P7 = f*R[7] + v*R[8];
    const T P8 = f*cam_t[1] + v*cam_t[2];
    const T P9 = R[2];
    const T P10 = R[5];
    const T P11 = R[8];
    const T P12 = cam_t[2];

    const T xu = T(m_pos_2dpoint[0]);
    const T yu = T(m_pos_2dpoint[1]);

    const T l0 = P7 - P11*yu;
    const T l1 = P11*xu - P3;
    const T l2 = P3*yu - P7*xu;

//    if(l0 != T(0.0))
//    {
//        out_residuals[0] = (- y_i*l1 - l2)/l0 - x_i;

//    }else{
//        out_residuals[0] = T(0.0);
//    }
//    if(l1 != T(0.0))
//    {
//        out_residuals[1] = ( -x_i*l0 - l2)/l1 - y_i;
//    }else{
//        out_residuals[1] = T(0.0);

//    }
//    out_residuals[0] = ( x_i*l0 + y_i*l1 + l2);
    out_residuals[0] = (x_i*l0 + y_i*l1 + l2)/sqrt(l0*l0 + l1*l1);
    out_residuals[1] = (x_i*l0 + y_i*l1 + l2)/sqrt(l0*l0 + l1*l1);//sqrt((xu - x_i)*(xu - x_i)+(yu - y_i)*(yu - y_i));
//    out_residuals[0] = abs(x_i*l0 + y_i*l1 + l2)/sqrt(l0*l0 + l1*l1);
//    out_residuals[1] = T(0.0);
//    out_residuals[0] = abs(x_i*l0 + y_i*l1 + l2)/sqrt(l0*l0 + l1*l1);



//    if(out_residuals[0] > T(10.0) )// && f > T(8000))
//    {
//        std::cout << "cam_R : " << cam_R[0] << " " << cam_R[1] << " " << cam_R[2] << std::endl;
//        std::cout << "R : " << R[0] << " " << R[1] << " " << R[2]
//                            << R[3] << " " << R[4] << " " << R[5]
//                            << R[6] << " " << R[7] << " " << R[8] << std::endl;
//        std::cout << "cam_t : " << cam_t[0] << " " << cam_t[1] << " " << cam_t[2] << std::endl;
//        std::cout << "f : " << f << " u : " << u << " v : " << v << std::endl;
//        std::cout << "xu : " << xu << " yu : " << yu << std::endl;
//        std::cout << "pos_3dpoint[0]: " << pos_3dpoint[0] << " pos_3dpoint[1]: " << pos_3dpoint[1] << " pos_3dpoint[2]: " << pos_3dpoint[2] << std::endl;
//        std::cout << "P3 : " << P3 << std::endl;
//        std::cout << "P7 : " << P7 << std::endl;
//        std::cout << "P11 : " << P11 << std::endl;
//        std::cout << "l0 : " << l0 << std::endl;
//        std::cout << "l1 : " << l1 << std::endl;
//        std::cout << "l2 : " << l2 << std::endl;
//        std::cout << "x_i : " << x_i << std::endl;
//        std::cout << "y_i : " << y_i << std::endl;
//        std::cout << "pos_proj[0] : " << pos_proj[0] << std::endl;
//        std::cout << "pos_proj[1] : " << pos_proj[1] << std::endl;
//        std::cout << "pos_proj[2] : " << pos_proj[2] << std::endl;
//        std::cout << "out_residuals[0]: "  << out_residuals[0] << std::endl;
//        std::cout << "****-----------------****" << std::endl;
//    }
//      }
//    if(out_residuals[0] > T(0.1))
//    {
//        std::cout << "cam_R : " << cam_R[0] << " " << cam_R[1] << " " << cam_R[2] << std::endl;
//        std::cout << "R : " << R[0] << " " << R[1] << " " << R[2]
//                            << R[3] << " " << R[4] << " " << R[5]
//                            << R[6] << " " << R[7] << " " << R[8] << std::endl;
//        std::cout << "cam_t : " << cam_t[0] << " " << cam_t[1] << " " << cam_t[2] << std::endl;
//        std::cout << "f : " << f << " u : " << u << " v : " << v << std::endl;
//        std::cout << "xu : " << xu << " yu : " << yu << std::endl;
//        std::cout << "pos_3dpoint[0]: " << pos_3dpoint[0] << " pos_3dpoint[1]: " << pos_3dpoint[1] << " pos_3dpoint[2]: " << pos_3dpoint[2] << std::endl;
//        std::cout << "out_residuals[0]: "  << out_residuals[0] << std::endl;
//        std::cout << "****-----------------****" << std::endl;
//    }

//      }

//    T pos_proj[3];
//    // Rotate the point according the camera rotation
//    ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_proj);

//    // Apply the camera translation
//    pos_proj[0] += cam_t[0];
//    pos_proj[1] += cam_t[1];
//    pos_proj[2] += cam_t[2];

//    // Transform the point from homogeneous to euclidean (undistorted point)
//    const T x_u = pos_proj[0] / pos_proj[2];
//    const T y_u = pos_proj[1] / pos_proj[2];

//    //--
//    // Apply intrinsic parameters
//    //--


//    const T& k1 = cam_intrinsics[OFFSET_DISTO_K1];
//    const T& k2 = cam_intrinsics[OFFSET_DISTO_K2];
//    const T& k3 = cam_intrinsics[OFFSET_DISTO_K3];
//    const T& t1 = cam_intrinsics[OFFSET_DISTO_T1];
//    const T& t2 = cam_intrinsics[OFFSET_DISTO_T2];

//    // Apply distortion (xd,yd) = disto(x_u,y_u)
//    const T r2 = x_u*x_u + y_u*y_u;
//    const T r4 = r2 * r2;
//    const T r6 = r4 * r2;
//    const T r_coeff = (1.0 + k1*r2 + k2*r4 + k3*r6);
//    const T t_x = t2 * (r2 + 2.0 * x_u*x_u) + 2.0 * t1 * x_u * y_u;
//    const T t_y = t1 * (r2 + 2.0 * y_u*y_u) + 2.0 * t2 * x_u * y_u;
//    const T x_d = x_u * r_coeff + t_x;
//    const T y_d = y_u * r_coeff + t_y;

//    // Apply focal length and principal point to get the final image coordinates
//    const T projected_x = principal_point_x + focal * x_d;
//    const T projected_y = principal_point_y + focal * y_d;

//    // Compute and return the error is the difference between the predicted
//    //  and observed position
//    out_residuals[0] = projected_x - m_pos_2dpoint[0];
//    out_residuals[1] = projected_y - m_pos_2dpoint[1];

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_1, 2, 3, 6, 3>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_1(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_1>, 2, 3, 6, 3>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_1>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_1(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};

struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_2
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_2(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
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

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--
//#pragma omp critical
//      {
    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];


    T pos_proj[3];
    T X[3];
    X[0] = pos_3dpoint[0];
    X[1] = pos_3dpoint[1];
    X[2] = T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);

    //prepare R
    const T theta = sqrt(cam_R[0]*cam_R[0] + cam_R[1]*cam_R[1] + cam_R[2]*cam_R[2]);
    const T cosTheta = cos(theta);
    const T sinTheta = sin(theta);
    const T r[3] = {cam_R[0]/theta, cam_R[1]/theta, cam_R[2]/theta};
    const T one_cosTheta = T(1.0) - cosTheta;

    T R[9];
    R[0] = cosTheta + one_cosTheta*r[0]*r[0];
    R[1] = one_cosTheta*r[0]*r[1] + sinTheta*r[2];
    R[2] = one_cosTheta*r[0]*r[2] - sinTheta*r[1];
    R[3] = one_cosTheta*r[0]*r[1] - sinTheta*r[2];
    R[4] = cosTheta + one_cosTheta*r[1]*r[1];
    R[5] = one_cosTheta*r[1]*r[2] + sinTheta*r[0];
    R[6] = one_cosTheta*r[0]*r[2] + sinTheta*r[1];
    R[7] = one_cosTheta*r[1]*r[2] - sinTheta*r[0];
    R[8] = cosTheta + one_cosTheta*r[2]*r[2];

    const T& f = cam_intrinsics[OFFSET_FOCAL_LENGTH]; //focal
    const T& u = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
    const T& v = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y

    // Transform the point from homogeneous to euclidean (undistorted point)
//    const T x_i = u + f * (pos_proj[0] / pos_proj[2]);
//    const T y_i = v + f * (pos_proj[1] / pos_proj[2]);

    const T P1 = f*R[0] + u*R[2];
    const T P2 = f*R[3] + u*R[5];
    const T P3 = f*R[6] + u*R[8];
    const T P4 = f*cam_t[0] + u*cam_t[2];
    const T P5 = f*R[1] + v*R[2];
    const T P6 = f*R[4] + v*R[5];
    const T P7 = f*R[7] + v*R[8];
    const T P8 = f*cam_t[1] + v*cam_t[2];
    const T P9 = R[2];
    const T P10 = R[5];
    const T P11 = R[8];
    const T P12 = cam_t[2];

    const T depth_X = (P9*X[0] + P10*X[1] + P12);

    const T x_i = (P1*X[0] + P2*X[1] + P4) / (P9*X[0] + P10*X[1] + P12);
    const T y_i = (P5*X[0] + P6*X[1] + P8) / (P9*X[0] + P10*X[1] + P12);

    const T xu = T(m_pos_2dpoint[0]);
    const T yu = T(m_pos_2dpoint[1]);

    const T nc_x = P3/P11;
    const T nc_y = P7/P11;
    //
    const T l_ncxi_x = x_i - nc_x;
    const T l_ncxi_y = y_i - nc_y;
//    const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
    const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
    //
    const T l_ncxu_x = xu - nc_x;
    const T l_ncxu_y = yu - nc_y;
//    const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);
    const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);
    const T A = yu - nc_y;
    const T B = nc_x - xu;
    const T C = xu*nc_y - yu*nc_x;

    const T normalXi = (x_i-nc_x) / sqrt((x_i-nc_x)*(x_i-nc_x) + (y_i-nc_y)*(y_i-nc_y));
    const T normalYi = (y_i-nc_y) / sqrt((x_i-nc_x)*(x_i-nc_x) + (y_i-nc_y)*(y_i-nc_y));
    const T D = abs(A*x_i + B*y_i + C) / sqrt(A*A+B*B);
    ///
    out_residuals[1] = D * normal_ncxu / T(2.0);
//    out_residuals[1] = acos((l_ncxi_x*l_ncxu_x + l_ncxi_y*l_ncxu_y) / (normal_ncxi*normal_ncxu));//* ((l_ncxi_x*l_ncxu_x + l_ncxi_y*l_ncxu_y) / (normal_ncxi*normal_ncxu));
//    out_residuals[2] = depth_X;

    {
        const T xu_cx = (xu - u)*depth_X/f - cam_t[0];
        const T xu_cy = (yu - v)*depth_X/f - cam_t[1];
        const T xu_cz = depth_X - cam_t[2];
        //
        const T xu_lx = R[0]*xu_cx + R[1]*xu_cy + R[2]*xu_cz;
        const T xu_ly = R[3]*xu_cx + R[4]*xu_cy + R[5]*xu_cz;
        const T xu_lz = R[6]*xu_cx + R[7]*xu_cy + R[8]*xu_cz;
//        const T xi_lx = (T)invR(0,0)*xi_cx + (T)invR(0,1)*xi_cy + (T)invR(0,2)*xi_cz;
//        const T xi_ly = (T)invR(1,0)*xi_cx + (T)invR(1,1)*xi_cy + (T)invR(1,2)*xi_cz;
//        const T xi_lz = (T)invR(2,0)*xi_cx + (T)invR(2,1)*xi_cy + (T)invR(2,2)*xi_cz;

        const T nc_cx = (nc_x - u)*depth_X/f - cam_t[0];
        const T nc_cy = (nc_y - v)*depth_X/f - cam_t[1];
        const T nc_cz = depth_X - cam_t[2];
        //
        const T nc_lx = R[0]*nc_cx + R[1]*nc_cy + R[2]*nc_cz;
        const T nc_ly = R[3]*nc_cx + R[4]*nc_cy + R[5]*nc_cz;
        const T nc_lz = R[6]*nc_cx + R[7]*nc_cy + R[8]*nc_cz;
//        const T nc_lx = (T)invR(0,0)*nc_cx + (T)invR(0,1)*nc_cy + (T)invR(0,2)*nc_cz;
//        const T nc_ly = (T)invR(1,0)*nc_cx + (T)invR(1,1)*nc_cy + (T)invR(1,2)*nc_cz;
//        const T nc_lz = (T)invR(2,0)*nc_cx + (T)invR(2,1)*nc_cy + (T)invR(2,2)*nc_cz;


        const T L_ncxu_x = nc_lx - xu_lx;
        const T L_ncxu_y = nc_ly - xu_ly;
        const T L_ncxu_z = nc_lz - xu_lz;

        const T L_Xxu_x = X[0] - xu_lx;
        const T L_Xxu_y = X[1] - xu_ly;
        const T L_Xxu_z = X[2] - xu_lz;

        const T cross_x = L_ncxu_y*L_Xxu_z - L_ncxu_z*L_Xxu_y;
        const T cross_y = L_ncxu_z*L_Xxu_x - L_ncxu_x*L_Xxu_z;
        const T cross_z = L_ncxu_x*L_Xxu_y - L_ncxu_y*L_Xxu_x;
        const T S = T(1.0/2.0) * sqrt((cross_x)*(cross_x) +(cross_y)*(cross_y) +(cross_z)*(cross_z));
//        std::cout << "S : " << S << std::endl;
        const T V = T(1.0/3.0) * S * depth_X;
//        std::cout << "V : " << V << std::endl;
//        std::cout << "d : " << depth_X << std::endl;
        out_residuals[0] = V;
        if(S == T(0.0))
        {
            out_residuals[0] = depth_X;
        }
//        out_residuals[1] = D * normal_ncxu / T(2.0);
        out_residuals[0] = S;

//        out_residuals[0] = sqrt((L_Xxu_x)*(L_Xxu_x) + (L_Xxu_y)*(L_Xxu_y) + (L_Xxu_z)*(L_Xxu_z));





//        out_residuals[0] = (l_ncxu_x*l_ncxu_x/normal_ncxu - l_ncxi_x*l_ncxi_x/normal_ncxi);
//        out_residuals[1] = (l_ncxu_y*l_ncxu_y/normal_ncxu - l_ncxi_y*l_ncxi_y/normal_ncxi);

//        out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi) * (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi);
//        out_residuals[1] = (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi) * (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);

//        out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi) * (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi)
//                         + (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi) * (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);

//        out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi) * (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi);
//        out_residuals[1] = (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi) * (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);

//        out_residuals[0] = (l_ncxu_x/normal_ncxu*l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi*l_ncxi_x/normal_ncxi)
//                         + (l_ncxu_y/normal_ncxu*l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi*l_ncxi_y/normal_ncxi);
    }

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::NumericDiffCostFunction//ceres::AutoDiffCostFunction//
           <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_2, ceres::CENTRAL, 2, 3, 6, 2>(
//        (new ceres::AutoDiffCostFunction
//          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_2, 2, 3, 6, 2>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_2(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_2>, 2, 3, 6, 2>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_2>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_2(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};

struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_S
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_S(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
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

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--
//#pragma omp critical
//      {
    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];


    T pos_proj[3];
    T X[3];
    X[0] = pos_3dpoint[0];
    X[1] = pos_3dpoint[1];
    X[2] = T(-0.613);//T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);

    //prepare R
    const T theta = sqrt(cam_R[0]*cam_R[0] + cam_R[1]*cam_R[1] + cam_R[2]*cam_R[2]);
    const T cosTheta = cos(theta);
    const T sinTheta = sin(theta);
    const T r[3] = {cam_R[0]/theta, cam_R[1]/theta, cam_R[2]/theta};
    const T one_cosTheta = T(1.0) - cosTheta;

    T R[9];
    R[0] = cosTheta + one_cosTheta*r[0]*r[0];
    R[1] = one_cosTheta*r[0]*r[1] + sinTheta*r[2];
    R[2] = one_cosTheta*r[0]*r[2] - sinTheta*r[1];
    R[3] = one_cosTheta*r[0]*r[1] - sinTheta*r[2];
    R[4] = cosTheta + one_cosTheta*r[1]*r[1];
    R[5] = one_cosTheta*r[1]*r[2] + sinTheta*r[0];
    R[6] = one_cosTheta*r[0]*r[2] + sinTheta*r[1];
    R[7] = one_cosTheta*r[1]*r[2] - sinTheta*r[0];
    R[8] = cosTheta + one_cosTheta*r[2]*r[2];

    const T& f = cam_intrinsics[OFFSET_FOCAL_LENGTH]; //focal
    const T& u = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
    const T& v = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y

    // Transform the point from homogeneous to euclidean (undistorted point)
//    const T x_i = u + f * (pos_proj[0] / pos_proj[2]);
//    const T y_i = v + f * (pos_proj[1] / pos_proj[2]);

    const T P1 = f*R[0] + u*R[2];
    const T P2 = f*R[3] + u*R[5];
    const T P3 = f*R[6] + u*R[8];
    const T P4 = f*cam_t[0] + u*cam_t[2];
    const T P5 = f*R[1] + v*R[2];
    const T P6 = f*R[4] + v*R[5];
    const T P7 = f*R[7] + v*R[8];
    const T P8 = f*cam_t[1] + v*cam_t[2];
    const T P9 = R[2];
    const T P10 = R[5];
    const T P11 = R[8];
    const T P12 = cam_t[2];

    const T depth_X = (P9*X[0] + P10*X[1] + P12);

    const T x_i = (P1*X[0] + P2*X[1] + P4) / (P9*X[0] + P10*X[1] + P12);
    const T y_i = (P5*X[0] + P6*X[1] + P8) / (P9*X[0] + P10*X[1] + P12);

    const T xu = T(m_pos_2dpoint[0]);
    const T yu = T(m_pos_2dpoint[1]);

    const T nc_x = P3/P11;
    const T nc_y = P7/P11;
    //
    const T l_ncxi_x = x_i - nc_x;
    const T l_ncxi_y = y_i - nc_y;
//    const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
    const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
    //
    const T l_ncxu_x = xu - nc_x;
    const T l_ncxu_y = yu - nc_y;
//    const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);
    const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);
    const T A = yu - nc_y;
    const T B = nc_x - xu;
    const T C = xu*nc_y - yu*nc_x;

    const T normalXi = (x_i-nc_x) / sqrt((x_i-nc_x)*(x_i-nc_x) + (y_i-nc_y)*(y_i-nc_y));
    const T normalYi = (y_i-nc_y) / sqrt((x_i-nc_x)*(x_i-nc_x) + (y_i-nc_y)*(y_i-nc_y));
    const T D = abs(A*x_i + B*y_i + C) / sqrt(A*A+B*B);
    ///

//    out_residuals[1] = acos((l_ncxi_x*l_ncxu_x + l_ncxi_y*l_ncxu_y) / (normal_ncxi*normal_ncxu));//* ((l_ncxi_x*l_ncxu_x + l_ncxi_y*l_ncxu_y) / (normal_ncxi*normal_ncxu));
//    out_residuals[2] = depth_X;

    {
        const T xu_cx = (xu - u)*depth_X/f - cam_t[0];
        const T xu_cy = (yu - v)*depth_X/f - cam_t[1];
        const T xu_cz = depth_X - cam_t[2];
        //
        const T xu_lx = R[0]*xu_cx + R[1]*xu_cy + R[2]*xu_cz;
        const T xu_ly = R[3]*xu_cx + R[4]*xu_cy + R[5]*xu_cz;
        const T xu_lz = R[6]*xu_cx + R[7]*xu_cy + R[8]*xu_cz;
//        const T xi_lx = (T)invR(0,0)*xi_cx + (T)invR(0,1)*xi_cy + (T)invR(0,2)*xi_cz;
//        const T xi_ly = (T)invR(1,0)*xi_cx + (T)invR(1,1)*xi_cy + (T)invR(1,2)*xi_cz;
//        const T xi_lz = (T)invR(2,0)*xi_cx + (T)invR(2,1)*xi_cy + (T)invR(2,2)*xi_cz;

        const T nc_cx = (nc_x - u)*depth_X/f - cam_t[0];
        const T nc_cy = (nc_y - v)*depth_X/f - cam_t[1];
        const T nc_cz = depth_X - cam_t[2];
        //
        const T nc_lx = R[0]*nc_cx + R[1]*nc_cy + R[2]*nc_cz;
        const T nc_ly = R[3]*nc_cx + R[4]*nc_cy + R[5]*nc_cz;
        const T nc_lz = R[6]*nc_cx + R[7]*nc_cy + R[8]*nc_cz;
//        const T nc_lx = (T)invR(0,0)*nc_cx + (T)invR(0,1)*nc_cy + (T)invR(0,2)*nc_cz;
//        const T nc_ly = (T)invR(1,0)*nc_cx + (T)invR(1,1)*nc_cy + (T)invR(1,2)*nc_cz;
//        const T nc_lz = (T)invR(2,0)*nc_cx + (T)invR(2,1)*nc_cy + (T)invR(2,2)*nc_cz;


        const T L_ncxu_x = nc_lx - xu_lx;
        const T L_ncxu_y = nc_ly - xu_ly;
        const T L_ncxu_z = nc_lz - xu_lz;

        const T L_Xxu_x = X[0] - xu_lx;
        const T L_Xxu_y = X[1] - xu_ly;
        const T L_Xxu_z = X[2] - xu_lz;

        const T cross_x = L_ncxu_y*L_Xxu_z - L_ncxu_z*L_Xxu_y;
        const T cross_y = L_ncxu_z*L_Xxu_x - L_ncxu_x*L_Xxu_z;
        const T cross_z = L_ncxu_x*L_Xxu_y - L_ncxu_y*L_Xxu_x;
        const T S = T(1.0/2.0) * sqrt((cross_x)*(cross_x) +(cross_y)*(cross_y) +(cross_z)*(cross_z));
//        std::cout << "S : " << S << std::endl;
        const T V = T(1.0/3.0) * S * depth_X;
//        std::cout << "V : " << V << std::endl;
//        std::cout << "d : " << depth_X << std::endl;
//        out_residuals[0] = V;
//        if(S == T(0.0))
//        {
//            out_residuals[0] = depth_X;
//        }
//        out_residuals[1] = D * normal_ncxu / T(2.0);
//        out_residuals[0] = S;
//        out_residuals[1] = D * normal_ncxu / T(2.0);

        const T l_xi_nc_x = nc_x - x_i;
        const T l_xi_nc_y = nc_y - y_i;
        const T l_xi_xu_x = xu - x_i;
        const T l_xi_xu_y = yu - y_i;
        const T cross_l = (l_xi_nc_x*l_xi_xu_y - l_xi_nc_y*l_xi_xu_x);
        out_residuals[0] = ((cross_l)/T(2.0));

//        out_residuals[0] = cross_x ;
//        out_residuals[1] = cross_y ;
//        out_residuals[2] = cross_z ;

//        out_residuals[4] = (cross_x)*(cross_x) +(cross_y)*(cross_y) +(cross_z)*(cross_z);

    }

    return true;
  }

  static int num_residuals() { return 1; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
//        (new ceres::NumericDiffCostFunction//ceres::AutoDiffCostFunction//
//           <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_S, ceres::RIDDERS, 1, 3, 6, 2>(
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_S, 1, 3, 6, 2>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_S(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_S>, 1, 3, 6, 2>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_S>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_S(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};

struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_pi
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_pi(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
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

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics1,
    const T* const cam_extrinsics1,
    const T* const cam_extrinsics2,
    const T* const cam_extrinsics3,
    T* out_residuals) const
  {
//      const T* const cam_intrinsics1,
//      const T* const cam_extrinsics1,
//      const T* const cam_intrinsics2,
//      const T* const cam_extrinsics2,
//      const T* const cam_intrinsics3,
//      const T* const cam_extrinsics3,
    //--
    // Apply external parameters (Pose)
    //--
//#pragma omp critical
      {
    const T * cam_R1 = cam_extrinsics1;
    const T * cam_t1 = &cam_extrinsics1[3];
    const T * cam_R2 = cam_extrinsics2;
    const T * cam_t2 = &cam_extrinsics2[3];
    const T * cam_R3 = cam_extrinsics3;
    const T * cam_t3 = &cam_extrinsics3[3];


//    T pos_proj[3];
//    T X[3];
//    X[0] = pos_3dpoint[0];
//    X[1] = pos_3dpoint[1];
//    X[2] = T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);


    //prepare R3
    const T theta1 = sqrt(cam_R1[0]*cam_R1[0] + cam_R1[1]*cam_R1[1] + cam_R1[2]*cam_R1[2]);
    const T cosTheta1 = cos(theta1);
    const T sinTheta1 = sin(theta1);
    const T r1[3] = {cam_R1[0]/theta1, cam_R1[1]/theta1, cam_R1[2]/theta1};
    const T one_cosTheta1 = T(1.0) - cosTheta1;

    T R1[9];
    R1[0] = cosTheta1 + one_cosTheta1*r1[0]*r1[0];
    R1[1] = one_cosTheta1*r1[0]*r1[1] + sinTheta1*r1[2];
    R1[2] = one_cosTheta1*r1[0]*r1[2] - sinTheta1*r1[1];
    R1[3] = one_cosTheta1*r1[0]*r1[1] - sinTheta1*r1[2];
    R1[4] = cosTheta1 + one_cosTheta1*r1[1]*r1[1];
    R1[5] = one_cosTheta1*r1[1]*r1[2] + sinTheta1*r1[0];
    R1[6] = one_cosTheta1*r1[0]*r1[2] + sinTheta1*r1[1];
    R1[7] = one_cosTheta1*r1[1]*r1[2] - sinTheta1*r1[0];
    R1[8] = cosTheta1 + one_cosTheta1*r1[2]*r1[2];

    const T& f1 = cam_intrinsics1[OFFSET_FOCAL_LENGTH]; //focal
    const T& u1 = cam_intrinsics1[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
    const T& v1 = cam_intrinsics1[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y

    //prepare R3
    const T theta2 = sqrt(cam_R2[0]*cam_R2[0] + cam_R2[1]*cam_R2[1] + cam_R2[2]*cam_R2[2]);
    const T cosTheta2 = cos(theta2);
    const T sinTheta2 = sin(theta2);
    const T r2[3] = {cam_R2[0]/theta2, cam_R2[1]/theta2, cam_R2[2]/theta2};
    const T one_cosTheta2 = T(1.0) - cosTheta2;

    T R2[9];
    R2[0] = cosTheta2 + one_cosTheta2*r2[0]*r2[0];
    R2[1] = one_cosTheta2*r2[0]*r2[1] + sinTheta2*r2[2];
    R2[2] = one_cosTheta2*r2[0]*r2[2] - sinTheta2*r2[1];
    R2[3] = one_cosTheta2*r2[0]*r2[1] - sinTheta2*r2[2];
    R2[4] = cosTheta2 + one_cosTheta2*r2[1]*r2[1];
    R2[5] = one_cosTheta2*r2[1]*r2[2] + sinTheta2*r2[0];
    R2[6] = one_cosTheta2*r2[0]*r2[2] + sinTheta2*r2[1];
    R2[7] = one_cosTheta2*r2[1]*r2[2] - sinTheta2*r2[0];
    R2[8] = cosTheta2 + one_cosTheta2*r2[2]*r2[2];

//    const T& f2 = cam_intrinsics2[OFFSET_FOCAL_LENGTH]; //focal
//    const T& u2 = cam_intrinsics2[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
//    const T& v2 = cam_intrinsics2[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y
    const T& f2 = cam_intrinsics1[OFFSET_FOCAL_LENGTH]; //focal
    const T& u2 = cam_intrinsics1[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
    const T& v2 = cam_intrinsics1[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y

    //prepare R3
    const T theta3 = sqrt(cam_R3[0]*cam_R3[0] + cam_R3[1]*cam_R3[1] + cam_R3[2]*cam_R3[2]);
    const T cosTheta3 = cos(theta3);
    const T sinTheta3 = sin(theta3);
    const T r3[3] = {cam_R3[0]/theta3, cam_R3[1]/theta3, cam_R3[2]/theta3};
    const T one_cosTheta3 = T(1.0) - cosTheta3;

    T R3[9];
    R3[0] = cosTheta3 + one_cosTheta3*r3[0]*r3[0];
    R3[1] = one_cosTheta3*r3[0]*r3[1] + sinTheta3*r3[2];
    R3[2] = one_cosTheta3*r3[0]*r3[2] - sinTheta3*r3[1];
    R3[3] = one_cosTheta3*r3[0]*r3[1] - sinTheta3*r3[2];
    R3[4] = cosTheta3 + one_cosTheta3*r3[1]*r3[1];
    R3[5] = one_cosTheta3*r3[1]*r3[2] + sinTheta3*r3[0];
    R3[6] = one_cosTheta3*r3[0]*r3[2] + sinTheta3*r3[1];
    R3[7] = one_cosTheta3*r3[1]*r3[2] - sinTheta3*r3[0];
    R3[8] = cosTheta3 + one_cosTheta3*r3[2]*r3[2];

//    const T& f3 = cam_intrinsics3[OFFSET_FOCAL_LENGTH]; //focal
//    const T& u3 = cam_intrinsics3[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
//    const T& v3 = cam_intrinsics3[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y
    const T& f3 = cam_intrinsics1[OFFSET_FOCAL_LENGTH]; //focal
    const T& u3 = cam_intrinsics1[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
    const T& v3 = cam_intrinsics1[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y

    // Transform the point from homogeneous to euclidean (undistorted point)
//    const T x_i = u + f * (pos_proj[0] / pos_proj[2]);
//    const T y_i = v + f * (pos_proj[1] / pos_proj[2]);

    const T P1_1 = f1*R1[0] + u1*R1[2];
    const T P2_1 = f1*R1[3] + u1*R1[5];
    const T P3_1 = f1*R1[6] + u1*R1[8];
    const T P4_1 = f1*cam_t1[0] + u1*cam_t1[2];
    const T P5_1 = f1*R1[1] + v1*R1[2];
    const T P6_1 = f1*R1[4] + v1*R1[5];
    const T P7_1 = f1*R1[7] + v1*R1[8];
    const T P8_1 = f1*cam_t1[1] + v1*cam_t1[2];
    const T P9_1 = R1[2];
    const T P10_1 = R1[5];
    const T P11_1 = R1[8];
    const T P12_1 = cam_t1[2];

    const T P1_2 = f2*R2[0] + u2*R2[2];
    const T P2_2 = f2*R2[3] + u2*R2[5];
    const T P3_2 = f2*R2[6] + u2*R2[8];
    const T P4_2 = f2*cam_t2[0] + u2*cam_t2[2];
    const T P5_2 = f2*R2[1] + v2*R2[2];
    const T P6_2 = f2*R2[4] + v2*R2[5];
    const T P7_2 = f2*R2[7] + v2*R2[8];
    const T P8_2 = f2*cam_t2[1] + v2*cam_t2[2];
    const T P9_2 = R2[2];
    const T P10_2 = R2[5];
    const T P11_2 = R2[8];
    const T P12_2 = cam_t2[2];

    const T P1_3 = f3*R3[0] + u3*R3[2];
    const T P2_3 = f3*R3[3] + u3*R3[5];
    const T P3_3 = f3*R3[6] + u3*R3[8];
    const T P4_3 = f3*cam_t3[0] + u3*cam_t3[2];
    const T P5_3 = f3*R3[1] + v3*R3[2];
    const T P6_3 = f3*R3[4] + v3*R3[5];
    const T P7_3 = f3*R3[7] + v3*R3[8];
    const T P8_3 = f3*cam_t3[1] + v3*cam_t3[2];
    const T P9_3 = R3[2];
    const T P10_3 = R3[5];
    const T P11_3 = R3[8];
    const T P12_3 = cam_t3[2];

    //
    const T nc_x_1 = P3_1/P11_1;
    const T nc_y_1 = P7_1/P11_1;

    const T nc_x_2 = P3_2/P11_2;
    const T nc_y_2 = P7_2/P11_2;

    const T nc_x_3 = P3_3/P11_3;
    const T nc_y_3 = P7_3/P11_3;

    //
    const T xu_1 = T(m_pos_2dpoint[0]);
    const T yu_1 = T(m_pos_2dpoint[1]);

    const T xu_2 = T(m_pos_2dpoint[2]);
    const T yu_2 = T(m_pos_2dpoint[3]);

    const T xu_3 = T(m_pos_2dpoint[4]);
    const T yu_3 = T(m_pos_2dpoint[5]);

    //
    const T corss_nc_xu_x_1 = nc_y_1 - yu_1;
    const T corss_nc_xu_y_1 = xu_1 - nc_x_1;
    const T corss_nc_xu_z_1 = nc_x_1*yu_1 - nc_y_1*xu_1;

    const T corss_nc_xu_x_2 = nc_y_2 - yu_2;
    const T corss_nc_xu_y_2 = xu_2 - nc_x_2;
    const T corss_nc_xu_z_2 = nc_x_2*yu_2 - nc_y_2*xu_2;

    const T corss_nc_xu_x_3 = nc_y_3 - yu_3;
    const T corss_nc_xu_y_3 = xu_3 - nc_x_3;
    const T corss_nc_xu_z_3 = nc_x_3*yu_3 - nc_y_3*xu_3;

    //
    T P_pi_a_1 = P1_1 * corss_nc_xu_x_1 + P5_1 * corss_nc_xu_y_1 + P9_1 * corss_nc_xu_z_1;
    T P_pi_b_1 = P2_1 * corss_nc_xu_x_1 + P6_1 * corss_nc_xu_y_1 + P10_1 * corss_nc_xu_z_1;
    T P_pi_c_1 = P3_1 * corss_nc_xu_x_1 + P7_1 * corss_nc_xu_y_1 + P11_1 * corss_nc_xu_z_1;
    T P_pi_d_1 = P4_1 * corss_nc_xu_x_1 + P8_1 * corss_nc_xu_y_1 + P12_1 * corss_nc_xu_z_1;
    T p1 = sqrt(P_pi_a_1*P_pi_a_1 + P_pi_b_1*P_pi_b_1 + P_pi_c_1*P_pi_c_1 + P_pi_d_1*P_pi_d_1);
    P_pi_a_1 /= p1;
    P_pi_b_1 /= p1;
    P_pi_c_1 /= p1;
    P_pi_d_1 /= p1;

    T P_pi_a_2 = P1_2 * corss_nc_xu_x_2 + P5_2 * corss_nc_xu_y_2 + P9_2 * corss_nc_xu_z_2;
    T P_pi_b_2 = P2_2 * corss_nc_xu_x_2 + P6_2 * corss_nc_xu_y_2 + P10_2 * corss_nc_xu_z_2;
    T P_pi_c_2 = P3_2 * corss_nc_xu_x_2 + P7_2 * corss_nc_xu_y_2 + P11_2 * corss_nc_xu_z_2;
    T P_pi_d_2 = P4_2 * corss_nc_xu_x_2 + P8_2 * corss_nc_xu_y_2 + P12_2 * corss_nc_xu_z_2;
    T p2 = sqrt(P_pi_a_2*P_pi_a_2 + P_pi_b_2*P_pi_b_2 + P_pi_c_2*P_pi_c_2 + P_pi_d_2*P_pi_d_2);
    P_pi_a_2 /= p2;
    P_pi_b_2 /= p2;
    P_pi_c_2 /= p2;
    P_pi_d_2 /= p2;

    T P_pi_a_3 = P1_3 * corss_nc_xu_x_3 + P5_3 * corss_nc_xu_y_3 + P9_3 * corss_nc_xu_z_3;
    T P_pi_b_3 = P2_3 * corss_nc_xu_x_3 + P6_3 * corss_nc_xu_y_3 + P10_3 * corss_nc_xu_z_3;
    T P_pi_c_3 = P3_3 * corss_nc_xu_x_3 + P7_3 * corss_nc_xu_y_3 + P11_3 * corss_nc_xu_z_3;
    T P_pi_d_3 = P4_3 * corss_nc_xu_x_3 + P8_3 * corss_nc_xu_y_3 + P12_3 * corss_nc_xu_z_3;
    T p3 = sqrt(P_pi_a_3*P_pi_a_3 + P_pi_b_3*P_pi_b_3 + P_pi_c_3*P_pi_c_3 + P_pi_d_3*P_pi_d_3);
    P_pi_a_3 /= p3;
    P_pi_b_3 /= p3;
    P_pi_c_3 /= p3;
    P_pi_d_3 /= p3;

    const T r_1 = (P_pi_a_1 * P_pi_b_2 * P_pi_c_3 + P_pi_a_3 * P_pi_b_1 * P_pi_c_2 + P_pi_a_2 * P_pi_b_3 * P_pi_c_1 - P_pi_a_3 * P_pi_b_2 * P_pi_c_1 - P_pi_a_1 * P_pi_b_3 * P_pi_c_2 - P_pi_a_2 * P_pi_b_1 * P_pi_c_3);
    const T r_2 = (P_pi_a_1 * P_pi_b_2 * P_pi_d_3 + P_pi_a_3 * P_pi_b_1 * P_pi_d_2 + P_pi_a_2 * P_pi_b_3 * P_pi_d_1 - P_pi_a_3 * P_pi_b_2 * P_pi_d_1 - P_pi_a_1 * P_pi_b_3 * P_pi_d_2 - P_pi_a_2 * P_pi_b_1 * P_pi_d_3);
    const T r_3 = (P_pi_a_1 * P_pi_c_2 * P_pi_d_3 + P_pi_a_3 * P_pi_c_1 * P_pi_d_2 + P_pi_a_2 * P_pi_c_3 * P_pi_d_1 - P_pi_a_3 * P_pi_c_2 * P_pi_d_1 - P_pi_a_1 * P_pi_c_3 * P_pi_d_2 - P_pi_a_2 * P_pi_c_1 * P_pi_d_3);
    const T r_4 = (P_pi_b_1 * P_pi_c_2 * P_pi_d_3 + P_pi_b_3 * P_pi_c_1 * P_pi_d_2 + P_pi_b_2 * P_pi_c_3 * P_pi_d_1 - P_pi_b_3 * P_pi_c_2 * P_pi_d_1 - P_pi_b_1 * P_pi_c_3 * P_pi_d_2 - P_pi_b_2 * P_pi_c_1 * P_pi_d_3);


    out_residuals[0] = r_1;
    out_residuals[1] = r_2;
    out_residuals[2] = r_3;
    out_residuals[3] = r_4;
    //
//    out_residuals[0] = (P_pi_a_1 * P_pi_b_2 * P_pi_c_3 + P_pi_a_3 * P_pi_b_1 * P_pi_c_2 + P_pi_a_2 * P_pi_b_3 * P_pi_c_1 - P_pi_a_3 * P_pi_b_2 * P_pi_c_1 - P_pi_a_1 * P_pi_b_3 * P_pi_c_2 - P_pi_a_2 * P_pi_b_1 * P_pi_c_3)
//                     + (P_pi_a_1 * P_pi_b_2 * P_pi_d_3 + P_pi_a_3 * P_pi_b_1 * P_pi_d_2 + P_pi_a_2 * P_pi_b_3 * P_pi_d_1 - P_pi_a_3 * P_pi_b_2 * P_pi_d_1 - P_pi_a_1 * P_pi_b_3 * P_pi_d_2 - P_pi_a_2 * P_pi_b_1 * P_pi_d_3)
//                     + (P_pi_a_1 * P_pi_c_2 * P_pi_d_3 + P_pi_a_3 * P_pi_c_1 * P_pi_d_2 + P_pi_a_2 * P_pi_c_3 * P_pi_d_1 - P_pi_a_3 * P_pi_c_2 * P_pi_d_1 - P_pi_a_1 * P_pi_c_3 * P_pi_d_2 - P_pi_a_2 * P_pi_c_1 * P_pi_d_3)
//                     + (P_pi_b_1 * P_pi_c_2 * P_pi_d_3 + P_pi_b_3 * P_pi_c_1 * P_pi_d_2 + P_pi_b_2 * P_pi_c_3 * P_pi_d_1 - P_pi_b_3 * P_pi_c_2 * P_pi_d_1 - P_pi_b_1 * P_pi_c_3 * P_pi_d_2 - P_pi_b_2 * P_pi_c_1 * P_pi_d_3);


//    std::cout << "r1 " << r_1 << std::endl;
//    std::cout << "r2 " << r_2 << std::endl;
//    std::cout << "r3 " << r_3 << std::endl;
//    std::cout << "r4 " << r_4 << std::endl;
//    std::cout << "----------------" << std::endl;
      }
//        out_residuals[0] = all;

    return true;
  }

  static int num_residuals() { return 4; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec6 & observation,
    const double weight = 0.0
  )
  {
//    if (weight == 0.0)
    {
      return
//        (new ceres::NumericDiffCostFunction//ceres::AutoDiffCostFunction//
//           <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_pi, ceres::RIDDERS, 1, 3, 6, 2>(
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_pi, 4, 3, 6, 6, 6>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_pi(observation.data())));
    }
//    else
//    {
//      return
//        (new ceres::AutoDiffCostFunction
//          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_pi>, 1, 3, 6, 3, 6, 3, 6>
//          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_pi>
//            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_pi(observation.data()), weight)));
//    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};

struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_D
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_D(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
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

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_extrinsics1,
    const T* const cam_extrinsics2,
    const T* const cam_extrinsics3,
    T* out_residuals) const
  {
//      const T* const cam_intrinsics1,
//      const T* const cam_extrinsics1,
//      const T* const cam_intrinsics2,
//      const T* const cam_extrinsics2,
//      const T* const cam_intrinsics3,
//      const T* const cam_extrinsics3,
    //--
    // Apply external parameters (Pose)
    //--
//#pragma omp critical
      {
    const T * cam_R1 = cam_extrinsics1;
    const T * cam_t1 = &cam_extrinsics1[3];
    const T * cam_R2 = cam_extrinsics2;
    const T * cam_t2 = &cam_extrinsics2[3];
    const T * cam_R3 = cam_extrinsics3;
    const T * cam_t3 = &cam_extrinsics3[3];

    const T * x1 = (const T *)&m_pos_2dpoint[0];
    const T * x2 = (const T *)&m_pos_2dpoint[3];
    const T * x3 = (const T *)&m_pos_2dpoint[6];


//    T pos_proj[3];
//    T X[3];
//    X[0] = pos_3dpoint[0];
//    X[1] = pos_3dpoint[1];
//    X[2] = T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);


    //prepare R3
    const T theta1 = sqrt(cam_R1[0]*cam_R1[0] + cam_R1[1]*cam_R1[1] + cam_R1[2]*cam_R1[2]);
    const T cosTheta1 = cos(theta1);
    const T sinTheta1 = sin(theta1);
    const T r1[3] = {cam_R1[0]/theta1, cam_R1[1]/theta1, cam_R1[2]/theta1};
    const T one_cosTheta1 = T(1.0) - cosTheta1;

    T R1[9];
    R1[0] = cosTheta1 + one_cosTheta1*r1[0]*r1[0];
    R1[1] = one_cosTheta1*r1[0]*r1[1] + sinTheta1*r1[2];
    R1[2] = one_cosTheta1*r1[0]*r1[2] - sinTheta1*r1[1];
    R1[3] = one_cosTheta1*r1[0]*r1[1] - sinTheta1*r1[2];
    R1[4] = cosTheta1 + one_cosTheta1*r1[1]*r1[1];
    R1[5] = one_cosTheta1*r1[1]*r1[2] + sinTheta1*r1[0];
    R1[6] = one_cosTheta1*r1[0]*r1[2] + sinTheta1*r1[1];
    R1[7] = one_cosTheta1*r1[1]*r1[2] - sinTheta1*r1[0];
    R1[8] = cosTheta1 + one_cosTheta1*r1[2]*r1[2];

//    const T& f1 = cam_intrinsics1[OFFSET_FOCAL_LENGTH]; //focal
//    const T& u1 = cam_intrinsics1[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
//    const T& v1 = cam_intrinsics1[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y

    //prepare R3
    const T theta2 = sqrt(cam_R2[0]*cam_R2[0] + cam_R2[1]*cam_R2[1] + cam_R2[2]*cam_R2[2]);
    const T cosTheta2 = cos(theta2);
    const T sinTheta2 = sin(theta2);
    const T r2[3] = {cam_R2[0]/theta2, cam_R2[1]/theta2, cam_R2[2]/theta2};
    const T one_cosTheta2 = T(1.0) - cosTheta2;

    T R2[9];
    R2[0] = cosTheta2 + one_cosTheta2*r2[0]*r2[0];
    R2[1] = one_cosTheta2*r2[0]*r2[1] + sinTheta2*r2[2];
    R2[2] = one_cosTheta2*r2[0]*r2[2] - sinTheta2*r2[1];
    R2[3] = one_cosTheta2*r2[0]*r2[1] - sinTheta2*r2[2];
    R2[4] = cosTheta2 + one_cosTheta2*r2[1]*r2[1];
    R2[5] = one_cosTheta2*r2[1]*r2[2] + sinTheta2*r2[0];
    R2[6] = one_cosTheta2*r2[0]*r2[2] + sinTheta2*r2[1];
    R2[7] = one_cosTheta2*r2[1]*r2[2] - sinTheta2*r2[0];
    R2[8] = cosTheta2 + one_cosTheta2*r2[2]*r2[2];

//    const T& f2 = cam_intrinsics2[OFFSET_FOCAL_LENGTH]; //focal
//    const T& u2 = cam_intrinsics2[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
//    const T& v2 = cam_intrinsics2[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y
//    const T& f2 = cam_intrinsics1[OFFSET_FOCAL_LENGTH]; //focal
//    const T& u2 = cam_intrinsics1[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
//    const T& v2 = cam_intrinsics1[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y

    //prepare R3
    const T theta3 = sqrt(cam_R3[0]*cam_R3[0] + cam_R3[1]*cam_R3[1] + cam_R3[2]*cam_R3[2]);
    const T cosTheta3 = cos(theta3);
    const T sinTheta3 = sin(theta3);
    const T r3[3] = {cam_R3[0]/theta3, cam_R3[1]/theta3, cam_R3[2]/theta3};
    const T one_cosTheta3 = T(1.0) - cosTheta3;

    T R3[9];
    R3[0] = cosTheta3 + one_cosTheta3*r3[0]*r3[0];
    R3[1] = one_cosTheta3*r3[0]*r3[1] + sinTheta3*r3[2];
    R3[2] = one_cosTheta3*r3[0]*r3[2] - sinTheta3*r3[1];
    R3[3] = one_cosTheta3*r3[0]*r3[1] - sinTheta3*r3[2];
    R3[4] = cosTheta3 + one_cosTheta3*r3[1]*r3[1];
    R3[5] = one_cosTheta3*r3[1]*r3[2] + sinTheta3*r3[0];
    R3[6] = one_cosTheta3*r3[0]*r3[2] + sinTheta3*r3[1];
    R3[7] = one_cosTheta3*r3[1]*r3[2] - sinTheta3*r3[0];
    R3[8] = cosTheta3 + one_cosTheta3*r3[2]*r3[2];

//    const T& f3 = cam_intrinsics3[OFFSET_FOCAL_LENGTH]; //focal
//    const T& u3 = cam_intrinsics3[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
//    const T& v3 = cam_intrinsics3[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y
//    const T& f3 = cam_intrinsics1[OFFSET_FOCAL_LENGTH]; //focal
//    const T& u3 = cam_intrinsics1[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
//    const T& v3 = cam_intrinsics1[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y

    // Transform the point from homogeneous to euclidean (undistorted point)
//    const T x_i = u + f * (pos_proj[0] / pos_proj[2]);
//    const T y_i = v + f * (pos_proj[1] / pos_proj[2]);

    const T H1_1 = R1[0];
    const T H2_1 = R1[1];
    const T H3_1 = R1[2];
    const T H4_1 = R1[3];
    const T H5_1 = R1[4];
    const T H6_1 = R1[5];
    const T H7_1 = cam_t1[0];
    const T H8_1 = cam_t1[1];
    const T H9_1 = cam_t1[2];

    const T H1_2 = R2[0];
    const T H2_2 = R2[1];
    const T H3_2 = R2[2];
    const T H4_2 = R2[3];
    const T H5_2 = R2[4];
    const T H6_2 = R2[5];
    const T H7_2 = cam_t2[0];
    const T H8_2 = cam_t2[1];
    const T H9_2 = cam_t2[2];

    const T H1_3 = R3[0];
    const T H2_3 = R3[1];
    const T H3_3 = R3[2];
    const T H4_3 = R3[3];
    const T H5_3 = R3[4];
    const T H6_3 = R3[5];
    const T H7_3 = cam_t3[0];
    const T H8_3 = cam_t3[1];
    const T H9_3 = cam_t3[2];

    //
    const T cross_x_r3_x_1 = x1[1]*R1[8] - x1[2]*R1[7];
    const T cross_x_r3_y_1 = x1[2]*R1[6] - x1[0]*R1[8];
    const T cross_x_r3_z_1 = x1[0]*R1[7] - x1[1]*R1[6];

    const T cross_x_r3_x_2 = x2[1]*R2[8] - x2[2]*R2[7];
    const T cross_x_r3_y_2 = x2[2]*R2[6] - x2[0]*R2[8];
    const T cross_x_r3_z_2 = x2[0]*R2[7] - x2[1]*R2[6];

    const T cross_x_r3_x_3 = x3[1]*R3[8] - x3[2]*R3[7];
    const T cross_x_r3_y_3 = x3[2]*R3[6] - x3[0]*R3[8];
    const T cross_x_r3_z_3 = x3[0]*R3[7] - x3[1]*R3[6];

    //
    const T lz_x_1 = H1_1 *  cross_x_r3_x_1 + H2_1 *  cross_x_r3_y_1 + H3_1 *  cross_x_r3_z_1;
    const T lz_y_1 = H4_1 *  cross_x_r3_x_1 + H5_1 *  cross_x_r3_y_1 + H6_1 *  cross_x_r3_z_1;
    const T lz_z_1 = H7_1 *  cross_x_r3_x_1 + H8_1 *  cross_x_r3_y_1 + H9_1 *  cross_x_r3_z_1;

    const T lz_x_2 = H1_2 *  cross_x_r3_x_2 + H2_2 *  cross_x_r3_y_2 + H3_2 *  cross_x_r3_z_2;
    const T lz_y_2 = H4_2 *  cross_x_r3_x_2 + H5_2 *  cross_x_r3_y_2 + H6_2 *  cross_x_r3_z_2;
    const T lz_z_2 = H7_2 *  cross_x_r3_x_2 + H8_2 *  cross_x_r3_y_2 + H9_2 *  cross_x_r3_z_2;

    const T lz_x_3 = H1_3 *  cross_x_r3_x_3 + H2_3 *  cross_x_r3_y_3 + H3_3 *  cross_x_r3_z_3;
    const T lz_y_3 = H4_3 *  cross_x_r3_x_3 + H5_3 *  cross_x_r3_y_3 + H6_3 *  cross_x_r3_z_3;
    const T lz_z_3 = H7_3 *  cross_x_r3_x_3 + H8_3 *  cross_x_r3_y_3 + H9_3 *  cross_x_r3_z_3;

    //
    const T p_12_x = lz_y_1 * lz_z_2 - lz_z_1 * lz_y_2;
    const T p_12_y = lz_z_1 * lz_x_2 - lz_x_1 * lz_z_2;
    const T p_12_z = lz_x_1 * lz_y_2 - lz_z_1 * lz_x_2;

    const T p_13_x = lz_y_1 * lz_z_3 - lz_z_1 * lz_y_3;
    const T p_13_y = lz_z_1 * lz_x_3 - lz_x_1 * lz_z_3;
    const T p_13_z = lz_x_1 * lz_y_3 - lz_z_1 * lz_x_3;

    const T p_23_x = lz_y_2 * lz_z_3 - lz_z_2 * lz_y_3;
    const T p_23_y = lz_z_2 * lz_x_3 - lz_x_2 * lz_z_3;
    const T p_23_z = lz_x_2 * lz_y_3 - lz_z_2 * lz_x_3;

    const T d1 = (lz_x_1*p_23_x + lz_y_1*p_23_y + lz_z_1*p_23_z)*(lz_x_1*p_23_x + lz_y_1*p_23_y + lz_z_1*p_23_z) / (lz_x_1*lz_x_1 + lz_y_1*lz_y_1)
       + (lz_x_2*p_13_x + lz_y_2*p_13_y + lz_z_2*p_13_z)*(lz_x_2*p_13_x + lz_y_2*p_13_y + lz_z_2*p_13_z) / (lz_x_2*lz_x_2 + lz_y_2*lz_y_2)
       + (lz_x_3*p_12_x + lz_y_3*p_12_y + lz_z_3*p_12_z)*(lz_x_3*p_12_x + lz_y_3*p_12_y + lz_z_3*p_12_z) / (lz_x_3*lz_x_3 + lz_y_3*lz_y_3);

    out_residuals[0] = d1;

    //
//    out_residuals[0] = (P_pi_a_1 * P_pi_b_2 * P_pi_c_3 + P_pi_a_3 * P_pi_b_1 * P_pi_c_2 + P_pi_a_2 * P_pi_b_3 * P_pi_c_1 - P_pi_a_3 * P_pi_b_2 * P_pi_c_1 - P_pi_a_1 * P_pi_b_3 * P_pi_c_2 - P_pi_a_2 * P_pi_b_1 * P_pi_c_3)
//                     + (P_pi_a_1 * P_pi_b_2 * P_pi_d_3 + P_pi_a_3 * P_pi_b_1 * P_pi_d_2 + P_pi_a_2 * P_pi_b_3 * P_pi_d_1 - P_pi_a_3 * P_pi_b_2 * P_pi_d_1 - P_pi_a_1 * P_pi_b_3 * P_pi_d_2 - P_pi_a_2 * P_pi_b_1 * P_pi_d_3)
//                     + (P_pi_a_1 * P_pi_c_2 * P_pi_d_3 + P_pi_a_3 * P_pi_c_1 * P_pi_d_2 + P_pi_a_2 * P_pi_c_3 * P_pi_d_1 - P_pi_a_3 * P_pi_c_2 * P_pi_d_1 - P_pi_a_1 * P_pi_c_3 * P_pi_d_2 - P_pi_a_2 * P_pi_c_1 * P_pi_d_3)
//                     + (P_pi_b_1 * P_pi_c_2 * P_pi_d_3 + P_pi_b_3 * P_pi_c_1 * P_pi_d_2 + P_pi_b_2 * P_pi_c_3 * P_pi_d_1 - P_pi_b_3 * P_pi_c_2 * P_pi_d_1 - P_pi_b_1 * P_pi_c_3 * P_pi_d_2 - P_pi_b_2 * P_pi_c_1 * P_pi_d_3);


//    std::cout << "r1 " << r_1 << std::endl;
//    std::cout << "r2 " << r_2 << std::endl;
//    std::cout << "r3 " << r_3 << std::endl;
//    std::cout << "r4 " << r_4 << std::endl;
//    std::cout << "----------------" << std::endl;
      }
//        out_residuals[0] = all;

    return true;
  }

  static int num_residuals() { return 1; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec9 & observation,
    const double weight = 0.0
  )
  {
//    if (weight == 0.0)
    {
      return
        (new ceres::NumericDiffCostFunction//ceres::AutoDiffCostFunction//
           <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_D, ceres::RIDDERS, 1, 6, 6, 6>(
//        (new ceres::AutoDiffCostFunction
//          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_D, 1, 6, 6, 6>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_D(observation.data())));
    }
//    else
//    {
//      return
//        (new ceres::AutoDiffCostFunction
//          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_pi>, 1, 3, 6, 3, 6, 3, 6>
//          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_pi>
//            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_pi(observation.data()), weight)));
//    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};


struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C4_S
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C4_S(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
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

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--
//#pragma omp critical
//      {
    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];


    T pos_proj[3];
    T X[3];
    X[0] = pos_3dpoint[0];
    X[1] = pos_3dpoint[1];
    X[2] = T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);

    //prepare R
    const T theta = sqrt(cam_R[0]*cam_R[0] + cam_R[1]*cam_R[1] + cam_R[2]*cam_R[2]);
    const T cosTheta = cos(theta);
    const T sinTheta = sin(theta);
    const T r[3] = {cam_R[0]/theta, cam_R[1]/theta, cam_R[2]/theta};
    const T one_cosTheta = T(1.0) - cosTheta;

    T R[9];
    R[0] = cosTheta + one_cosTheta*r[0]*r[0];
    R[1] = one_cosTheta*r[0]*r[1] + sinTheta*r[2];
    R[2] = one_cosTheta*r[0]*r[2] - sinTheta*r[1];
    R[3] = one_cosTheta*r[0]*r[1] - sinTheta*r[2];
    R[4] = cosTheta + one_cosTheta*r[1]*r[1];
    R[5] = one_cosTheta*r[1]*r[2] + sinTheta*r[0];
    R[6] = one_cosTheta*r[0]*r[2] + sinTheta*r[1];
    R[7] = one_cosTheta*r[1]*r[2] - sinTheta*r[0];
    R[8] = cosTheta + one_cosTheta*r[2]*r[2];

    const T& f = cam_intrinsics[OFFSET_FOCAL_LENGTH]; //focal
    const T& u = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
    const T& v = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y
    const T& k1 = cam_intrinsics[OFFSET_DISTO_K1];
    const T& k2 = cam_intrinsics[OFFSET_DISTO_K2];
    const T& k3 = cam_intrinsics[OFFSET_DISTO_K3];
    const T& t1 = cam_intrinsics[OFFSET_DISTO_T1];
    const T& t2 = cam_intrinsics[OFFSET_DISTO_T2];

    // Transform the point from homogeneous to euclidean (undistorted point)
//    const T x_i = u + f * (pos_proj[0] / pos_proj[2]);
//    const T y_i = v + f * (pos_proj[1] / pos_proj[2]);

    const T P1 = f*R[0] + u*R[2];
    const T P2 = f*R[3] + u*R[5];
    const T P3 = f*R[6] + u*R[8];
    const T P4 = f*cam_t[0] + u*cam_t[2];
    const T P5 = f*R[1] + v*R[2];
    const T P6 = f*R[4] + v*R[5];
    const T P7 = f*R[7] + v*R[8];
    const T P8 = f*cam_t[1] + v*cam_t[2];
    const T P9 = R[2];
    const T P10 = R[5];
    const T P11 = R[8];
    const T P12 = cam_t[2];

    const T depth_X = (P9*X[0] + P10*X[1] + P12);

    const T x_i = (P1*X[0] + P2*X[1] + P4) / (P9*X[0] + P10*X[1] + P12);
    const T y_i = (P5*X[0] + P6*X[1] + P8) / (P9*X[0] + P10*X[1] + P12);

    const T xd = T(m_pos_2dpoint[0]);//sift
    const T yd = T(m_pos_2dpoint[1]);

    // Apply distortion (xd,yd) = disto(x_u,y_u)
    const T r2 = xd*xd + yd*yd;
    const T r4 = r2 * r2;
    const T r6 = r4 * r2;
    const T r_coeff = (1.0 + k1*r2 + k2*r4 + k3*r6);
    const T t_x = t2 * (r2 + 2.0 * xd*xd) + 2.0 * t1 * xd * yd;
    const T t_y = t1 * (r2 + 2.0 * yd*yd) + 2.0 * t2 * xd * yd;
    const T xu = xd * r_coeff + t_x;
    const T yu = yd * r_coeff + t_y;



    const T nc_x = P3/P11;
    const T nc_y = P7/P11;
    //
    const T l_ncxi_x = x_i - nc_x;
    const T l_ncxi_y = y_i - nc_y;
//    const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
    const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
    //
    const T l_ncxu_x = xu - nc_x;
    const T l_ncxu_y = yu - nc_y;
//    const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);
    const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);
    const T A = yu - nc_y;
    const T B = nc_x - xu;
    const T C = xu*nc_y - yu*nc_x;

    const T normalXi = (x_i-nc_x) / sqrt((x_i-nc_x)*(x_i-nc_x) + (y_i-nc_y)*(y_i-nc_y));
    const T normalYi = (y_i-nc_y) / sqrt((x_i-nc_x)*(x_i-nc_x) + (y_i-nc_y)*(y_i-nc_y));
    const T D = abs(A*x_i + B*y_i + C) / sqrt(A*A+B*B);
    ///
//    out_residuals[1] = acos((l_ncxi_x*l_ncxu_x + l_ncxi_y*l_ncxu_y) / (normal_ncxi*normal_ncxu));//* ((l_ncxi_x*l_ncxu_x + l_ncxi_y*l_ncxu_y) / (normal_ncxi*normal_ncxu));
//    out_residuals[2] = depth_X;

    {
        const T xu_cx = (xu - u)*depth_X/f - cam_t[0];
        const T xu_cy = (yu - v)*depth_X/f - cam_t[1];
        const T xu_cz = depth_X - cam_t[2];
        //
        const T xu_lx = R[0]*xu_cx + R[1]*xu_cy + R[2]*xu_cz;
        const T xu_ly = R[3]*xu_cx + R[4]*xu_cy + R[5]*xu_cz;
        const T xu_lz = R[6]*xu_cx + R[7]*xu_cy + R[8]*xu_cz;
//        const T xi_lx = (T)invR(0,0)*xi_cx + (T)invR(0,1)*xi_cy + (T)invR(0,2)*xi_cz;
//        const T xi_ly = (T)invR(1,0)*xi_cx + (T)invR(1,1)*xi_cy + (T)invR(1,2)*xi_cz;
//        const T xi_lz = (T)invR(2,0)*xi_cx + (T)invR(2,1)*xi_cy + (T)invR(2,2)*xi_cz;

        const T nc_cx = (nc_x - u)*depth_X/f - cam_t[0];
        const T nc_cy = (nc_y - v)*depth_X/f - cam_t[1];
        const T nc_cz = depth_X - cam_t[2];
        //
        const T nc_lx = R[0]*nc_cx + R[1]*nc_cy + R[2]*nc_cz;
        const T nc_ly = R[3]*nc_cx + R[4]*nc_cy + R[5]*nc_cz;
        const T nc_lz = R[6]*nc_cx + R[7]*nc_cy + R[8]*nc_cz;
//        const T nc_lx = (T)invR(0,0)*nc_cx + (T)invR(0,1)*nc_cy + (T)invR(0,2)*nc_cz;
//        const T nc_ly = (T)invR(1,0)*nc_cx + (T)invR(1,1)*nc_cy + (T)invR(1,2)*nc_cz;
//        const T nc_lz = (T)invR(2,0)*nc_cx + (T)invR(2,1)*nc_cy + (T)invR(2,2)*nc_cz;


        const T L_ncxu_x = nc_lx - xu_lx;
        const T L_ncxu_y = nc_ly - xu_ly;
        const T L_ncxu_z = nc_lz - xu_lz;

        const T L_Xxu_x = X[0] - xu_lx;
        const T L_Xxu_y = X[1] - xu_ly;
        const T L_Xxu_z = X[2] - xu_lz;

        const T cross_x = L_ncxu_y*L_Xxu_z - L_ncxu_z*L_Xxu_y;
        const T cross_y = L_ncxu_z*L_Xxu_x - L_ncxu_x*L_Xxu_z;
        const T cross_z = L_ncxu_x*L_Xxu_y - L_ncxu_y*L_Xxu_x;
        const T S = T(1.0/2.0) * sqrt((cross_x)*(cross_x) +(cross_y)*(cross_y) +(cross_z)*(cross_z));
//        std::cout << "S : " << S << std::endl;
        const T V = T(1.0/3.0) * S * depth_X;

//        out_residuals[0] = V;
//        if(S == T(0.0))
//        {
//            out_residuals[0] = depth_X;
//        }
//        out_residuals[1] = D * normal_ncxu / T(2.0);
//        out_residuals[0] = S;

//        out_residuals[1] = D * normal_ncxu / T(2.0);
//        out_residuals[0] = D * normal_ncxu / T(2.0);

        const T l_xi_nc_x = nc_x - x_i;
        const T l_xi_nc_y = nc_y - y_i;
        const T l_xi_xu_x = xu - x_i;
        const T l_xi_xu_y = yu - y_i;
        const T cross_l = (l_xi_nc_x*l_xi_xu_y - l_xi_nc_y*l_xi_xu_x);

        out_residuals[0] = cross_l/T(2.0);// / T(2.0);//*cross_l;// s
//        out_residuals[1] = D*D;// (cross_x)*(cross_x) +(cross_y)*(cross_y) +(cross_z)*(cross_z);
//        out_residuals[1] = S;

//        out_residuals[0] = cross_l*cross_l;//s
//        out_residuals[1] = ((cross_x)*(cross_x) +(cross_y)*(cross_y) +(cross_z)*(cross_z)); //S

//        out_residuals[0] = (cross_x)*(cross_x) +(cross_y)*(cross_y) +(cross_z)*(cross_z);//S;//abs(cross_l);//*cross_l;//s


    }

    return true;
  }

  static int num_residuals() { return 1; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
//        (new ceres::AutoDiffCostFunction//
//            <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C4_S, 1, 8, 6, 2>(
        (new ceres::NumericDiffCostFunction//ceres::AutoDiffCostFunction//
           <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C4_S, ceres::RIDDERS, 1, 8, 6, 2>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C4_S(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C4_S>, 1, 8, 6, 2>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C4_S>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C4_S(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};


struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3_old
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3_old(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
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

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--
//#pragma omp critical
//      {
    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];


    T pos_proj[3];
    T X[3];
    X[0] = pos_3dpoint[0];
    X[1] = pos_3dpoint[1];
    X[2] = T(0.0);// pos_3dpoint[2];//T(0.0);
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, X, pos_proj);



    //prepare l
    T R[9];
//    ceres::MatrixAdapter<T, 3, 3>(R);
    ceres::AngleAxisToRotationMatrix(cam_R, R);

//    pos_proj[0] = R[0]*pos_3dpoint[0] + R[3]*pos_3dpoint[1];
//    pos_proj[1] = R[1]*pos_3dpoint[0] + R[4]*pos_3dpoint[1];
//    pos_proj[2] = R[2]*pos_3dpoint[0] + R[5]*pos_3dpoint[1];

    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    T pos_proj2[3];
    pos_proj2[0] = R[0]*pos_3dpoint[0] + R[3]*pos_3dpoint[1];
    pos_proj2[1] = R[1]*pos_3dpoint[0] + R[4]*pos_3dpoint[1];
    pos_proj2[2] = R[2]*pos_3dpoint[0] + R[5]*pos_3dpoint[1];
    pos_proj2[0] += cam_t[0];
    pos_proj2[1] += cam_t[1];
    pos_proj2[2] += cam_t[2];

//    if(pos_proj[0] != pos_proj2[0] || pos_proj[1] != pos_proj2[1] || pos_proj[2] != pos_proj2[2])
//    {
//        std::cout << "pos_proj : " << pos_proj[0] << std::endl << pos_proj[1] << std::endl << pos_proj[2] << std::endl;
//        std::cout << "pos_proj2 : " << pos_proj2[0] << std::endl << pos_proj2[1] << std::endl << pos_proj2[2] << std::endl;
//        std::cout << "**------**" << std::endl;
//    }



    const T& f = cam_intrinsics[OFFSET_FOCAL_LENGTH]; //focal
    const T& u = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
    const T& v = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y

    // Transform the point from homogeneous to euclidean (undistorted point)
    const T x_i = u + f * (pos_proj[0] / pos_proj[2]);
    const T y_i = v + f * (pos_proj[1] / pos_proj[2]);

    const T P1 = f*R[0] + u*R[2];
    const T P2 = f*R[3] + u*R[5];
    const T P3 = f*R[6] + u*R[8];
    const T P4 = f*cam_t[0] + u*cam_t[2];
    const T P5 = f*R[1] + v*R[2];
    const T P6 = f*R[4] + v*R[5];
    const T P7 = f*R[7] + v*R[8];
    const T P8 = f*cam_t[1] + v*cam_t[2];
    const T P9 = R[2];
    const T P10 = R[5];
    const T P11 = R[8];
    const T P12 = cam_t[2];

    const T xu = T(m_pos_2dpoint[0]);
    const T yu = T(m_pos_2dpoint[1]);

    const T l0 = P7 - P11*yu;
    const T l1 = P11*xu - P3;
    const T l2 = P3*yu - P7*xu;
//    const T l0 = y_i*P11 - P7;//P7 - P11*yu;
//    const T l1 = P3 - P11*x_i;//P11*xu - P3;
//    const T l2 = P7*x_i - P3*y_i;//P3*yu - P7*xu;

    const T A = l0*P1 + l1*P5 + l2*P9 ;
    const T B = l0*P2 + l1*P6 + l2*P10 ;
    const T C = l0*P3 + l1*P7 + l2*P11 ;
    const T D = l0*P4 + l1*P8 + l2*P12 ;

    out_residuals[0] = abs(pos_3dpoint[0]*A + pos_3dpoint[1]*B + pos_3dpoint[2]*C + D)/sqrt(A*A + B*B + C*C);

//    out_residuals[0] = abs(xu*l0 + yu*l1 + l2)/sqrt(l0*l0 + l1*l1);
//    out_residuals[1] = abs(xu*l0 + yu*l1 + l2)/sqrt(l0*l0 + l1*l1);

//    if(l0 != T(0.0))
//    {
//        out_residuals[0] = (- y_i*l1 - l2)/l0 - x_i;

//    }else{
//        out_residuals[0] = T(0.0);
//    }
//    if(l1 != T(0.0))
//    {
//        out_residuals[1] = ( -x_i*l0 - l2)/l1 - y_i;
//    }else{
//        out_residuals[1] = T(0.0);

//    }
//    out_residuals[0] = ( x_i*l0 + y_i*l1 + l2);
//    out_residuals[0] = abs(x_i*l0 + y_i*l1 + l2)/sqrt(l0*l0 + l1*l1);
//    out_residuals[1] = T(0.0);

//    out_residuals[0] = abs(x_i*l0 + y_i*l1 + l2)/sqrt(l0*l0 + l1*l1);

//    if(out_residuals[0] > T(10.0) && f > T(8000))
//    {
//        std::cout << "cam_R : " << cam_R[0] << " " << cam_R[1] << " " << cam_R[2] << std::endl;
//        std::cout << "R : " << R[0] << " " << R[1] << " " << R[2]
//                            << R[3] << " " << R[4] << " " << R[5]
//                            << R[6] << " " << R[7] << " " << R[8] << std::endl;
//        std::cout << "cam_t : " << cam_t[0] << " " << cam_t[1] << " " << cam_t[2] << std::endl;
//        std::cout << "f : " << f << " u : " << u << " v : " << v << std::endl;
//        std::cout << "xu : " << xu << " yu : " << yu << std::endl;
//        std::cout << "pos_3dpoint[0]: " << pos_3dpoint[0] << " pos_3dpoint[1]: " << pos_3dpoint[1] << " pos_3dpoint[2]: " << pos_3dpoint[2] << std::endl;
//        std::cout << "P3 : " << P3 << std::endl;
//        std::cout << "P7 : " << P7 << std::endl;
//        std::cout << "P11 : " << P11 << std::endl;
//        std::cout << "l0 : " << l0 << std::endl;
//        std::cout << "l1 : " << l1 << std::endl;
//        std::cout << "l2 : " << l2 << std::endl;
//        std::cout << "x_i : " << x_i << std::endl;
//        std::cout << "y_i : " << y_i << std::endl;
//        std::cout << "pos_proj[0] : " << pos_proj[0] << std::endl;
//        std::cout << "pos_proj[1] : " << pos_proj[1] << std::endl;
//        std::cout << "pos_proj[2] : " << pos_proj[2] << std::endl;
//        std::cout << "out_residuals[0]: "  << out_residuals[0] << std::endl;
//        std::cout << "****-----------------****" << std::endl;
//    }

//      }


    return true;
  }

  static int num_residuals() { return 1; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3_old, 1, 3, 6, 2>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3_old(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3_old>, 1, 3, 6, 2>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3_old>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3_old(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};

//struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3
//{
//  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3(const double* const pos_2dpoint)
//  :m_pos_2dpoint(pos_2dpoint)
//  {
//  }

//  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
//  enum : uint8_t {
//    OFFSET_FOCAL_LENGTH = 0,
//    OFFSET_PRINCIPAL_POINT_X = 1,
//    OFFSET_PRINCIPAL_POINT_Y = 2,
//    OFFSET_DISTO_K1 = 3,
//    OFFSET_DISTO_K2 = 4,
//    OFFSET_DISTO_K3 = 5,
//    OFFSET_DISTO_T1 = 6,
//    OFFSET_DISTO_T2 = 7,
//  };

//  /**
//   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
//   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
//   *   - 3 for rotation(angle axis), 3 for translation
//   * @param[in] pos_3dpoint
//   * @param[out] out_residuals
//   */
//  template <typename T>
//  bool operator()(
//    const T* const cam_intrinsics,
//    const T* const cam_extrinsics,
//    const T* const pos_3dpoint,
//    T* out_residuals) const
//  {
//    //--
//    // Apply external parameters (Pose)
//    //--
////#pragma omp critical
////      {
//    const T * cam_R = cam_extrinsics;
//    const T * cam_t = &cam_extrinsics[3];


//    T pos_proj[3];
//    T X[3];
//    X[0] = pos_3dpoint[0];
//    X[1] = pos_3dpoint[1];
//    X[2] = T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);
//    // Rotate the point according the camera rotation
////    ceres::AngleAxisRotatePoint(cam_R, X, pos_proj);

////    // Apply the camera translation
////    pos_proj[0] += cam_t[0];
////    pos_proj[1] += cam_t[1];
////    pos_proj[2] += cam_t[2];

////    //prepare R
//    const T theta = sqrt(cam_R[0]*cam_R[0] + cam_R[1]*cam_R[1] + cam_R[2]*cam_R[2]);
//    const T cosTheta = cos(theta);
//    const T sinTheta = sin(theta);
//    const T r[3] = {cam_R[0]/theta, cam_R[1]/theta, cam_R[2]/theta};
//    const T one_cosTheta = T(1.0) - cosTheta;

//    T R[9];
//    R[0] = cosTheta + one_cosTheta*r[0]*r[0];
//    R[1] = one_cosTheta*r[0]*r[1] + sinTheta*r[2];
//    R[2] = one_cosTheta*r[0]*r[2] - sinTheta*r[1];
//    R[3] = one_cosTheta*r[0]*r[1] - sinTheta*r[2];
//    R[4] = cosTheta + one_cosTheta*r[1]*r[1];
//    R[5] = one_cosTheta*r[1]*r[2] + sinTheta*r[0];
//    R[6] = one_cosTheta*r[0]*r[2] + sinTheta*r[1];
//    R[7] = one_cosTheta*r[1]*r[2] - sinTheta*r[0];
//    R[8] = cosTheta + one_cosTheta*r[2]*r[2];

////    ceres::MatrixAdapter<T, 3, 3>(R);
////    ceres::AngleAxisToRotationMatrix(cam_R, R);

//    const T& f = cam_intrinsics[OFFSET_FOCAL_LENGTH]; //focal
//    const T& u = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
//    const T& v = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y

//    // Transform the point from homogeneous to euclidean (undistorted point)
////    const T x_i = u + f * (pos_proj[0] / pos_proj[2]);
////    const T y_i = v + f * (pos_proj[1] / pos_proj[2]);

//    const T P1 = f*R[0] + u*R[2];
//    const T P2 = f*R[3] + u*R[5];
//    const T P3 = f*R[6] + u*R[8];
//    const T P4 = f*cam_t[0] + u*cam_t[2];
//    const T P5 = f*R[1] + v*R[2];
//    const T P6 = f*R[4] + v*R[5];
//    const T P7 = f*R[7] + v*R[8];
//    const T P8 = f*cam_t[1] + v*cam_t[2];
//    const T P9 = R[2];
//    const T P10 = R[5];
//    const T P11 = R[8];
//    const T P12 = cam_t[2];

//    const T x_i = (P1*X[0] + P2*X[1] + P3*X[2] + P4) / (P9*X[0] + P10*X[1] + P11*X[2] + P12);
//    const T y_i = (P5*X[0] + P6*X[1] + P7*X[2] + P8) / (P9*X[0] + P10*X[1] + P11*X[2] + P12);

//    const T xu = T(m_pos_2dpoint[0]);
//    const T yu = T(m_pos_2dpoint[1]);

////    const T l0 = P7 - P11*yu;
////    const T l1 = P11*xu - P3;
////    const T l2 = P3*yu - P7*xu;

////    const T l0 = y_i*P11 - P7;//P7 - P11*yu;
////    const T l1 = P3 - P11*x_i;//P11*xu - P3;
////    const T l2 = P7*x_i - P3*y_i;//P3*yu - P7*xu;
////    const T l0_ = P7 - P11*yu;
////    const T l1_ = P11*xu - P3;
////    const T l2_ = P3*yu - P7*xu;

////    const T P3_ = P3 / P11;
////    const T P7_ = P7 / P11;
////    const T l0 = y_i - P7_;//P7 - P11*yu;
////    const T l1 = P3_ - x_i;//P11*xu - P3;
////    const T l2 = P7_*x_i - P3_*y_i;//P3*yu - P7*xu;
////    out_residuals[0] = abs(xu*l0 + yu*l1 + l2)/sqrt(l0*l0 + l1*l1);
////    const T l0_ = yu - P7_;
////    const T l1_ = P3_ - xu;
////    const T l2_ = P7_*xu - P3_*yu;
////    out_residuals[1] = abs(x_i*l0_ + y_i*l1_ + l2_)/sqrt(l0_*l0_ + l1_*l1_);

//    const T nc_x = P3/P11;
//    const T nc_y = P7/P11;
//    //
//    const T l_ncxi_x = x_i - nc_x;
//    const T l_ncxi_y = y_i - nc_y;
//    const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
////    const T normal_ncxi = (l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
//    //
//    const T l_ncxu_x = xu - nc_x;
//    const T l_ncxu_y = yu - nc_y;
//    const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);
////    const T normal_ncxu = (l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);

//    const T A = yu - nc_y;
//    const T B = nc_x - xu;
//    const T C = xu*nc_y - yu*nc_x;

//    const T normalXi = (x_i-nc_x) / sqrt((x_i-nc_x)*(x_i-nc_x) + (y_i-nc_y)*(y_i-nc_y));
//    const T normalYi = (y_i-nc_y) / sqrt((x_i-nc_x)*(x_i-nc_x) + (y_i-nc_y)*(y_i-nc_y));
//    const T D = abs(A*x_i + B*y_i + C) / sqrt(A*A+B*B);
//    ///
////    const T x2 = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi)*(l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi);
////    const T y2 = (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi)*(l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);
//    {
////        out_residuals[0] = (l_ncxu_x*l_ncxi_y - l_ncxu_y*l_ncxi_x) / T(2.0)  / normal_ncxi / normal_ncxu;// / D;
//        out_residuals[0] = D * normal_ncxu / T(2.0);



////        T X[3];
////        X[0] = pos_3dpoint[0];
////        X[1] = pos_3dpoint[1];
////        X[2] = pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);


////        const T x_i = (P1*X[0] + P2*X[1] + P3*X[2] + P4) / (P9*X[0] + P10*X[1] + P11*X[2] + P12);
////        const T y_i = (P5*X[0] + P6*X[1] + P7*X[2] + P8) / (P9*X[0] + P10*X[1] + P11*X[2] + P12);


////        const T nc_x = P3/P11;
////        const T nc_y = P7/P11;
////        //
////        const T l_ncxi_x = x_i - nc_x;
////        const T l_ncxi_y = y_i - nc_y;
////        const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
////    //    const T normal_ncxi = (l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
////        //
////        const T l_ncxu_x = xu - nc_x;
////        const T l_ncxu_y = yu - nc_y;
////        const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);
////    //    const T normal_ncxu = (l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);

////        out_residuals[1] = (x_i - xu);// sqrt(*(l_ncxi_x/normal_ncxi-l_ncxu_x/normal_ncxu) + *(l_ncxi_y/normal_ncxi-l_ncxu_y/normal_ncxu));
////        out_residuals[2] = (y_i - yu);

////        out_residuals[1] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi)* (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi);
////        out_residuals[2] = (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi)* (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);

////          out_residuals[0] = (l_ncxu_x*l_ncxi_y - l_ncxu_y*l_ncxi_x) / T(2.0);//


////        out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi);
////        out_residuals[1] = (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);

////        out_residuals[0] = (l_ncxu_x*l_ncxu_x/normal_ncxu - l_ncxi_x*l_ncxi_x/normal_ncxi);
////        out_residuals[1] = (l_ncxu_y*l_ncxu_y/normal_ncxu - l_ncxi_y*l_ncxi_y/normal_ncxi);

////        out_residuals[2] = sqrt((x_i - xu)*(x_i - xu) + (y_i - yu)*(y_i - yu));


////        out_residuals[0] = (l_ncxu_x*l_ncxu_x/normal_ncxu - l_ncxi_x*l_ncxi_x/normal_ncxi);
////        out_residuals[1] = (l_ncxu_y*l_ncxu_y/normal_ncxu - l_ncxi_y*l_ncxi_y/normal_ncxi);

////        out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi) * (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi);
////        out_residuals[1] = (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi) * (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);

////        out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi) * (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi)
////                         + (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi) * (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);

////        out_residuals[0] = (l_ncxu_x/normal_ncxu*l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi*l_ncxi_x/normal_ncxi)
////                         + (l_ncxu_y/normal_ncxu*l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi*l_ncxi_y/normal_ncxi);
//    }

//    return true;
//  }

//  static int num_residuals() { return 1; }

//  // Factory to hide the construction of the CostFunction object from
//  // the client code.
//  static ceres::CostFunction* Create
//  (
//    const Vec2 & observation,
//    const double weight = 0.0
//  )
//  {
//    if (weight == 0.0)
//    {
//      return
//          (new ceres::NumericDiffCostFunction//ceres::AutoDiffCostFunction//
//            <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3, ceres::RIDDERS, 1, 3, 6, 2>(
////        (new ceres::AutoDiffCostFunction
////          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3, 1, 3, 6, 3>(
//            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3(observation.data())));
//    }
//    else
//    {
//      return
//        (new ceres::NumericDiffCostFunction//ceres::AutoDiffCostFunction
////          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3, ceres::RIDDERS, 1, 3, 6, 2>(
//          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3>, ceres::RIDDERS, 1, 3, 6, 2>
//          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3>
//            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3(observation.data()), weight)));
//    }
//  }

//  const double * m_pos_2dpoint; // The 2D observation
//};

struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
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

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--
//#pragma omp critical
//      {
    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];


    T pos_proj[3];
    T X[3];
    X[0] = pos_3dpoint[0];
    X[1] = pos_3dpoint[1];
    X[2] = T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);

    //prepare R
    const T theta = sqrt(cam_R[0]*cam_R[0] + cam_R[1]*cam_R[1] + cam_R[2]*cam_R[2]);
    const T cosTheta = cos(theta);
    const T sinTheta = sin(theta);
    const T r[3] = {cam_R[0]/theta, cam_R[1]/theta, cam_R[2]/theta};
    const T one_cosTheta = T(1.0) - cosTheta;

    T R[9];
    R[0] = cosTheta + one_cosTheta*r[0]*r[0];
    R[1] = one_cosTheta*r[0]*r[1] + sinTheta*r[2];
    R[2] = one_cosTheta*r[0]*r[2] - sinTheta*r[1];
    R[3] = one_cosTheta*r[0]*r[1] - sinTheta*r[2];
    R[4] = cosTheta + one_cosTheta*r[1]*r[1];
    R[5] = one_cosTheta*r[1]*r[2] + sinTheta*r[0];
    R[6] = one_cosTheta*r[0]*r[2] + sinTheta*r[1];
    R[7] = one_cosTheta*r[1]*r[2] - sinTheta*r[0];
    R[8] = cosTheta + one_cosTheta*r[2]*r[2];

    const T& f = cam_intrinsics[OFFSET_FOCAL_LENGTH]; //focal
    const T& u = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
    const T& v = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y

    // Transform the point from homogeneous to euclidean (undistorted point)
//    const T x_i = u + f * (pos_proj[0] / pos_proj[2]);
//    const T y_i = v + f * (pos_proj[1] / pos_proj[2]);

    const T P1 = f*R[0] + u*R[2];
    const T P2 = f*R[3] + u*R[5];
    const T P3 = f*R[6] + u*R[8];
    const T P4 = f*cam_t[0] + u*cam_t[2];
    const T P5 = f*R[1] + v*R[2];
    const T P6 = f*R[4] + v*R[5];
    const T P7 = f*R[7] + v*R[8];
    const T P8 = f*cam_t[1] + v*cam_t[2];
    const T P9 = R[2];
    const T P10 = R[5];
    const T P11 = R[8];
    const T P12 = cam_t[2];

    const T depth_X = (P9*X[0] + P10*X[1] + P12);

    const T x_i = (P1*X[0] + P2*X[1] + P4) / (P9*X[0] + P10*X[1] + P12);
    const T y_i = (P5*X[0] + P6*X[1] + P8) / (P9*X[0] + P10*X[1] + P12);

    const T xu = T(m_pos_2dpoint[0]);
    const T yu = T(m_pos_2dpoint[1]);

    const T nc_x = P3/P11;
    const T nc_y = P7/P11;
    //
    const T l_ncxi_x = x_i - nc_x;
    const T l_ncxi_y = y_i - nc_y;
//    const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
    const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
    //
    const T l_ncxu_x = xu - nc_x;
    const T l_ncxu_y = yu - nc_y;
//    const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);
    const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);
    const T A = yu - nc_y;
    const T B = nc_x - xu;
    const T C = xu*nc_y - yu*nc_x;

    const T normalXi = (x_i-nc_x) / sqrt((x_i-nc_x)*(x_i-nc_x) + (y_i-nc_y)*(y_i-nc_y));
    const T normalYi = (y_i-nc_y) / sqrt((x_i-nc_x)*(x_i-nc_x) + (y_i-nc_y)*(y_i-nc_y));
    const T D = abs(A*x_i + B*y_i + C) / sqrt(A*A+B*B);
    ///
//    out_residuals[0] = D * normal_ncxu / T(2.0);
//    out_residuals[1] = ((l_ncxi_x*l_ncxu_x + l_ncxi_y*l_ncxu_y) / (normal_ncxi*normal_ncxu))* ((l_ncxi_x*l_ncxu_x + l_ncxi_y*l_ncxu_y) / (normal_ncxi*normal_ncxu));
//    out_residuals[1] = depth_X;

    {
        const T xu_cx = (xu - u)*depth_X/f - cam_t[0];
        const T xu_cy = (yu - v)*depth_X/f - cam_t[1];
        const T xu_cz = depth_X - cam_t[2];
        //
        const T xu_lx = R[0]*xu_cx + R[1]*xu_cy + R[2]*xu_cz;
        const T xu_ly = R[3]*xu_cx + R[4]*xu_cy + R[5]*xu_cz;
        const T xu_lz = R[6]*xu_cx + R[7]*xu_cy + R[8]*xu_cz;
//        const T xi_lx = (T)invR(0,0)*xi_cx + (T)invR(0,1)*xi_cy + (T)invR(0,2)*xi_cz;
//        const T xi_ly = (T)invR(1,0)*xi_cx + (T)invR(1,1)*xi_cy + (T)invR(1,2)*xi_cz;
//        const T xi_lz = (T)invR(2,0)*xi_cx + (T)invR(2,1)*xi_cy + (T)invR(2,2)*xi_cz;

        const T nc_cx = (nc_x - u)*depth_X/f - cam_t[0];
        const T nc_cy = (nc_y - v)*depth_X/f - cam_t[1];
        const T nc_cz = depth_X - cam_t[2];
        //
        const T nc_lx = R[0]*nc_cx + R[1]*nc_cy + R[2]*nc_cz;
        const T nc_ly = R[3]*nc_cx + R[4]*nc_cy + R[5]*nc_cz;
        const T nc_lz = R[6]*nc_cx + R[7]*nc_cy + R[8]*nc_cz;
//        const T nc_lx = (T)invR(0,0)*nc_cx + (T)invR(0,1)*nc_cy + (T)invR(0,2)*nc_cz;
//        const T nc_ly = (T)invR(1,0)*nc_cx + (T)invR(1,1)*nc_cy + (T)invR(1,2)*nc_cz;
//        const T nc_lz = (T)invR(2,0)*nc_cx + (T)invR(2,1)*nc_cy + (T)invR(2,2)*nc_cz;


        const T L_ncxu_x = nc_lx - xu_lx;
        const T L_ncxu_y = nc_ly - xu_ly;
        const T L_ncxu_z = nc_lz - xu_lz;

        const T L_Xxu_x = X[0] - xu_lx;
        const T L_Xxu_y = X[1] - xu_ly;
        const T L_Xxu_z = X[2] - xu_lz;

        const T cross_x = L_ncxu_y*L_Xxu_z - L_ncxu_z*L_Xxu_y;
        const T cross_y = L_ncxu_z*L_Xxu_x - L_ncxu_x*L_Xxu_z;
        const T cross_z = L_ncxu_x*L_Xxu_y - L_ncxu_y*L_Xxu_x;
        const T S = T(1.0/2.0) * sqrt((cross_x)*(cross_x) +(cross_y)*(cross_y) +(cross_z)*(cross_z));
//        std::cout << "S : " << S << std::endl;
        const T V = T(1.0/3.0) * S * depth_X;
//        std::cout << "V : " << V << std::endl;
//        std::cout << "d : " << depth_X << std::endl;
        out_residuals[0] = V;

        out_residuals[1] = D;
//        out_residuals[1] = V;
//        if(S == T(0.0))
//        {
//            out_residuals[1] = depth_X;
//        }



//        out_residuals[0] = (l_ncxu_x*l_ncxu_x/normal_ncxu - l_ncxi_x*l_ncxi_x/normal_ncxi);
//        out_residuals[1] = (l_ncxu_y*l_ncxu_y/normal_ncxu - l_ncxi_y*l_ncxi_y/normal_ncxi);

//        out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi) * (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi);
//        out_residuals[1] = (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi) * (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);

//        out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi) * (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi)
//                         + (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi) * (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);

//        out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi) * (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi);
//        out_residuals[1] = (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi) * (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);

//        out_residuals[0] = (l_ncxu_x/normal_ncxu*l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi*l_ncxi_x/normal_ncxi)
//                         + (l_ncxu_y/normal_ncxu*l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi*l_ncxi_y/normal_ncxi);
    }

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
//        (new ceres::NumericDiffCostFunction//ceres::AutoDiffCostFunction//
//           <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3, ceres::RIDDERS, 2, 3, 6, 2>(
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3, 2, 3, 6, 2>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3>, 1, 3, 6, 2>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};


//struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3
//{
//  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3(const double* const pos_2dpoint)
//  :m_pos_2dpoint(pos_2dpoint)
//  {
//  }

//  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
//  enum : uint8_t {
//    OFFSET_FOCAL_LENGTH = 0,
//    OFFSET_PRINCIPAL_POINT_X = 1,
//    OFFSET_PRINCIPAL_POINT_Y = 2,
//    OFFSET_DISTO_K1 = 3,
//    OFFSET_DISTO_K2 = 4,
//    OFFSET_DISTO_K3 = 5,
//    OFFSET_DISTO_T1 = 6,
//    OFFSET_DISTO_T2 = 7,
//  };

//  /**
//   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
//   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
//   *   - 3 for rotation(angle axis), 3 for translation
//   * @param[in] pos_3dpoint
//   * @param[out] out_residuals
//   */
//  template <typename T>
//  bool operator()(
//    const T* const cam_intrinsics,
//    const T* const cam_extrinsics,
//    const T* const pos_3dpoint,
//    T* out_residuals) const
//  {
//    //--
//    // Apply external parameters (Pose)
//    //--
////#pragma omp critical
////      {
//    const T * cam_R = cam_extrinsics;
//    const T * cam_t = &cam_extrinsics[3];


//    T pos_proj[3];
//    T X[3];
//    X[0] = pos_3dpoint[0];
//    X[1] = pos_3dpoint[1];
//    X[2] = T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);
//    // Rotate the point according the camera rotation
////    ceres::AngleAxisRotatePoint(cam_R, X, pos_proj);

////    //prepare R
//    const T theta = sqrt(cam_R[0]*cam_R[0] + cam_R[1]*cam_R[1] + cam_R[2]*cam_R[2]);
//    const T cosTheta = cos(theta);
//    const T sinTheta = sin(theta);
//    const T r[3] = {cam_R[0]/theta, cam_R[1]/theta, cam_R[2]/theta};
//    const T one_cosTheta = T(1.0) - cosTheta;

//    T R[9];
//    R[0] = cosTheta + one_cosTheta*r[0]*r[0];
//    R[1] = one_cosTheta*r[0]*r[1] + sinTheta*r[2];
//    R[2] = one_cosTheta*r[0]*r[2] - sinTheta*r[1];
//    R[3] = one_cosTheta*r[0]*r[1] - sinTheta*r[2];
//    R[4] = cosTheta + one_cosTheta*r[1]*r[1];
//    R[5] = one_cosTheta*r[1]*r[2] + sinTheta*r[0];
//    R[6] = one_cosTheta*r[0]*r[2] + sinTheta*r[1];
//    R[7] = one_cosTheta*r[1]*r[2] - sinTheta*r[0];
//    R[8] = cosTheta + one_cosTheta*r[2]*r[2];

//    const T& f = cam_intrinsics[OFFSET_FOCAL_LENGTH]; //focal
//    const T& u = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
//    const T& v = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y
//    //
//    const T P1 = f*R[0] + u*R[2];
//    const T P2 = f*R[3] + u*R[5];
//    const T P3 = f*R[6] + u*R[8];
//    const T P4 = f*cam_t[0] + u*cam_t[2];
//    const T P5 = f*R[1] + v*R[2];
//    const T P6 = f*R[4] + v*R[5];
//    const T P7 = f*R[7] + v*R[8];
//    const T P8 = f*cam_t[1] + v*cam_t[2];
//    const T P9 = R[2];
//    const T P10 = R[5];
//    const T P11 = R[8];
//    const T P12 = cam_t[2];
//    //
//    const T x_i = (P1*X[0] + P2*X[1] + P3*X[2] + P4) / (P9*X[0] + P10*X[1] + P11*X[2] + P12);
//    const T y_i = (P5*X[0] + P6*X[1] + P7*X[2] + P8) / (P9*X[0] + P10*X[1] + P11*X[2] + P12);
//    //
//    const T xu = T(m_pos_2dpoint[0]);
//    const T yu = T(m_pos_2dpoint[1]);
//    //
//    const T nc_x = P3/P11;
//    const T nc_y = P7/P11;
//    //
////    const T l_ncxi_x = x_i - nc_x;
////    const T l_ncxi_y = y_i - nc_y;
////    const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
////    const T normal_ncxi = (l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
//    //
//    const T l_ncxu_x = xu - nc_x;
//    const T l_ncxu_y = yu - nc_y;
//    const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);
////    const T normal_ncxu = (l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);

//    const T A = yu - nc_y;
//    const T B = nc_x - xu;
//    const T C = xu*nc_y - yu*nc_x;

//    const T D = abs(A*x_i + B*y_i + C) / sqrt(A*A+B*B);
//    ///
//    {
//        out_residuals[0] = D * normal_ncxu / T(2.0);
//    }

//    {
////        const T l_ncxi_x = x_i - nc_x;
////        const T l_ncxi_y = y_i - nc_y;
////        const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);

////        out_residuals[0] = abs(l_ncxi_x/normal_ncxi * l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi * l_ncxu_x/normal_ncxu) / T(2.0);

//    }
//    return true;
//  }

//  static int num_residuals() { return 1; }

//  // Factory to hide the construction of the CostFunction object from
//  // the client code.
//  static ceres::CostFunction* Create
//  (
//    const Vec2 & observation,
//    const double weight = 0.0
//  )
//  {
//    if (weight == 0.0)
//    {
//      return
////          (new ceres::NumericDiffCostFunction//ceres::AutoDiffCostFunction//
////            <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1, ceres::RIDDERS, 1, 3, 6, 2>(
//        (new ceres::AutoDiffCostFunction
//          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3, 1, 3, 6, 3>(
//            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3(observation.data())));
//    }
//    else
//    {
//      return
//        (new ceres::NumericDiffCostFunction//ceres::AutoDiffCostFunction
////          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3, ceres::RIDDERS, 1, 3, 6, 2>(
//          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3>, ceres::RIDDERS, 1, 3, 6, 3>
//          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3>
//            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3(observation.data()), weight)));
//    }
//  }

//  const double * m_pos_2dpoint; // The 2D observation
//};


class Rat43Analytic : public ceres::SizedCostFunction<1,4> {
   public:
     Rat43Analytic(const double x, const double y) : x_(x), y_(y) {}
     virtual ~Rat43Analytic() {}
     virtual bool Evaluate(double const* const* parameters,
                           double* residuals,
                           double** jacobians) const {
       const double b1 = parameters[0][0];
       const double b2 = parameters[0][1];
       const double b3 = parameters[0][2];
       const double b4 = parameters[0][3];

       residuals[0] = b1 *  pow(1 + exp(b2 -  b3 * x_), -1.0 / b4) - y_;

       if (!jacobians) return true;
       double* jacobian = jacobians[0];
       if (!jacobian) return true;

       jacobian[0] = pow(1 + exp(b2 - b3 * x_), -1.0 / b4);
       jacobian[1] = -b1 * exp(b2 - b3 * x_) *
                     pow(1 + exp(b2 - b3 * x_), -1.0 / b4 - 1) / b4;
       jacobian[2] = x_ * b1 * exp(b2 - b3 * x_) *
                     pow(1 + exp(b2 - b3 * x_), -1.0 / b4 - 1) / b4;
       jacobian[3] = b1 * log(1 + exp(b2 - b3 * x_)) *
                     pow(1 + exp(b2 - b3 * x_), -1.0 / b4) / (b4 * b4);
       return true;
     }

    private:
     const double x_;
     const double y_;
 };

template <typename CostFunctor,
          int kNumResiduals,  // Number of residuals, or ceres::DYNAMIC.
          int N0,       // Number of parameters in block 0.
          int N1 = 0,   // Number of parameters in block 1.
          int N2 = 0,   // Number of parameters in block 2.
          int N3 = 0,   // Number of parameters in block 3.
          int N4 = 0,   // Number of parameters in block 4.
          int N5 = 0,   // Number of parameters in block 5.
          int N6 = 0,   // Number of parameters in block 6.
          int N7 = 0,   // Number of parameters in block 7.
          int N8 = 0,   // Number of parameters in block 8.
          int N9 = 0>   // Number of parameters in block 9.
class NewDiffCostFunction : public ceres::SizedCostFunction<kNumResiduals,
                                                      N0, N1, N2, N3, N4,
                                                      N5, N6, N7, N8, N9> {
 public:
  // Takes ownership of functor. Uses the template-provided value for the
  // number of residuals ("kNumResiduals").
  explicit NewDiffCostFunction(CostFunctor* functor)
      : functor_(functor) {
    CHECK_NE(kNumResiduals, ceres::DYNAMIC)
        << "Can't run the fixed-size constructor if the "
        << "number of residuals is set to ceres::DYNAMIC.";
  }

    NewDiffCostFunction(CostFunctor* functor, Vec2 _m_pos_2dpoint)
          : functor_(functor), m_pos_2dpoint(_m_pos_2dpoint) {
        CHECK_NE(kNumResiduals, ceres::DYNAMIC)
            << "Can't run the fixed-size constructor if the "
            << "number of residuals is set to ceres::DYNAMIC.";
      }

  // Takes ownership of functor. Ignores the template-provided
  // kNumResiduals in favor of the "num_residuals" argument provided.
  //
  // This allows for having autodiff cost functions which return varying
  // numbers of residuals at runtime.
  NewDiffCostFunction(CostFunctor* functor, int num_residuals)
      : functor_(functor) {
    CHECK_EQ(kNumResiduals, ceres::DYNAMIC)
        << "Can't run the dynamic-size constructor if the "
        << "number of residuals is not ceres::DYNAMIC.";
    ceres::SizedCostFunction<kNumResiduals,
                      N0, N1, N2, N3, N4,
                      N5, N6, N7, N8, N9>
        ::set_num_residuals(num_residuals);
  }

  virtual ~NewDiffCostFunction() {}

  // Implementation details follow; clients of the autodiff cost function should
  // not have to examine below here.
  //
  // To handle varardic cost functions, some template magic is needed. It's
  // mostly hidden inside autodiff.h.
  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) const {
//#pragma omp critical
//      {
//      if (!jacobians)
//          return true;
//      double *jacobian = jacobians[0];
//      if (!jacobian)
//          return true;
//      ceres::internal::AutoDiff<CostFunctor, double,
//             N0, N1, N2, N3, N4, N5, N6, N7, N8, N9>::Differentiate(
//                 *functor_,
//                 parameters,
//                 ceres::SizedCostFunction<kNumResiduals,
//                                   N0, N1, N2, N3, N4,
//                                   N5, N6, N7, N8, N9>::num_residuals(),
//                 residuals,
//                 jacobians);

      const double *cam_R = parameters[1];
      const double *cam_t = &(parameters[1][3]);

      double pos_proj[3];
      double X[3];
      X[0] = parameters[2][0];
      X[1] = parameters[2][1];
      X[2] = (0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);

//      std::cout << "error 4" << std::endl;
      double theta = sqrt(cam_R[0]*cam_R[0] + cam_R[1]*cam_R[1] + cam_R[2]*cam_R[2]);
      double r[3] = {cam_R[0]/theta, cam_R[1]/theta, cam_R[2]/theta};

      double R[9];
      R[0] = cos(theta) + r[0]*r[0];
      R[1] = r[0]*r[1] + sin(theta)*r[2];
      R[2] = r[0]*r[2] - sin(theta)*r[1];
      R[3] = r[0]*r[1] - sin(theta)*r[2];
      R[4] = cos(theta) + r[1]*r[1];
      R[5] = r[1]*r[2] + sin(theta)*r[0];
      R[6] = r[0]*r[2] + sin(theta)*r[1];
      R[7] = r[1]*r[2] - sin(theta)*r[0];
      R[8] = cos(theta) + r[2]*r[2];

      const double& f = parameters[0][0]; //focal
      const double& u = parameters[0][1]; //principal_point_x
      const double& v = parameters[0][2]; //principal_point_y

      // Apply the camera translation
      pos_proj[0] = R[0]*X[0] + R[3]*X[1] + cam_t[0];
      pos_proj[1] = R[1]*X[0] + R[4]*X[1] + cam_t[1];
      pos_proj[2] = R[2]*X[0] + R[5]*X[1] + cam_t[2];

//      std::cout << "error 6" << std::endl;

      // Transform the point from homogeneous to euclidean (undistorted point)
      const double x_i = u + f * (pos_proj[0] / pos_proj[2]);
      const double y_i = v + f * (pos_proj[1] / pos_proj[2]);

      const double P1 = f*R[0] + u*R[2];
      const double P2 = f*R[3] + u*R[5];
      const double P3 = f*R[6] + u*R[8];
      const double P4 = f*cam_t[0] + u*cam_t[2];
      const double P5 = f*R[1] + v*R[2];
      const double P6 = f*R[4] + v*R[5];
      const double P7 = f*R[7] + v*R[8];
      const double P8 = f*cam_t[1] + v*cam_t[2];
      const double P9 = R[2];
      const double P10 = R[5];
      const double P11 = R[8];
      const double P12 = cam_t[2];



      const double xu = m_pos_2dpoint[0];
      const double yu = m_pos_2dpoint[1];

//      m_pos_2dpoint
      const double nc_x = P3/P11;
      const double nc_y = P7/P11;
      //
      const double l_ncxi_x = x_i - nc_x;
      const double l_ncxi_y = y_i - nc_y;
      const double normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
  //    const T normal_ncxi = (l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
      //
      const double l_ncxu_x = xu - nc_x;
      const double l_ncxu_y = yu - nc_y;
      const double normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);

      residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi);
      residuals[1] = (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);


//      std::cout << "error 0" << std::endl;
//      if (!jacobians) {
//        return ceres::internal::VariadicEvaluate<
//            CostFunctor, double, N0, N1, N2, N3, N4, N5, N6, N7, N8, N9>
//            ::Call(*functor_, parameters, residuals);
//      }
//      return internal::AutoDiff<CostFunctor, double,
//             N0, N1, N2, N3, N4, N5, N6, N7, N8, N9>::Differentiate(
//                 *functor_,
//                 parameters,
//                 SizedCostFunction<kNumResiduals,
//                                   N0, N1, N2, N3, N4,
//                                   N5, N6, N7, N8, N9>::num_residuals(),
//                 residuals,
//                 jacobians);

//      if (!jacobians) {
//          std::cout << "error 1" << std::endl;
////        return ceres::internal::VariadicEvaluate<
////            CostFunctor, double, N0, N1, N2, N3, N4, N5, N6, N7, N8, N9>
////            ::Call(*functor_, parameters, residuals);
//      }

//      std::cout << "j0 " << ((R[6]/P11)/normal_ncxu + l_ncxu_x*(R[6]/P11*(l_ncxu_x/normal_ncxu) + R[7]/P11*(l_ncxu_x/normal_ncxu)))
//              - ((pos_proj[0]/pos_proj[2] - R[6]/P11)/normal_ncxi + l_ncxi_x*((pos_proj[0]/pos_proj[2] - R[6]/P11)*(l_ncxi_x/normal_ncxi) + (pos_proj[1]/pos_proj[2] - R[7]/P11)*(l_ncxi_x/normal_ncxi))) << std::endl;
//     bool flag = true;
      if (jacobians == NULL)
     {
//         std::cout << "error 1 1" << std::endl;
//         flag = false;
        return true;
     }
//    double *jacobian = jacobians[0];
//    if (jacobian == NULL)
//    {
////        std::cout << "error 1 2" << std::endl;
////        flag = false;
//        return true;
//    }
    // 0,f



      if(jacobians != NULL)// && jacobians[0] != NULL)
      {
          std::cout << "error 1 3" << std::endl;
          jacobians[0][0] = -1.0;
          std::cout << "error 1 4" << std::endl;


          std::cout << "error 1 5" << std::endl;
        std::cout << "error 4 2 " << jacobians[0][0] << std::endl;
        std::cout << "error 1 6" << std::endl;
      }
      //    std::cout << "error 1 3" << std::endl;
//    ceres::internal::AutoDiff<CostFunctor, double,
//           N0, N1, N2, N3, N4, N5, N6, N7, N8, N9>::Differentiate(
//               *functor_,
//               parameters,
//               ceres::SizedCostFunction<kNumResiduals,
//                                 N0, N1, N2, N3, N4,
//                                 N5, N6, N7, N8, N9>::num_residuals(),
//               residuals,
//               jacobians);


    //if(flag)
    {
////      std::cout << "error 3" << std::endl;
//      if(P11 == 0.0 || normal_ncxu == 0.0 || pos_proj[2] == 0.0 || normal_ncxi == 0.0)
//      {
//          std::cout << "error 2" << std::endl;
//      }
////      std::cout << "error 4 " << normal_ncxu << std::endl;
////      std::cout << "size " << jacobians.size() << std::endl;
////      std::cout << "size " << jacobians[0].size() << std::endl;
////      std::cout << "error 4 1 " << jacobians[1] << std::endl;
////      std::cout << "error 4 2 " << jacobian[3] << std::endl;
////      jacobian[3] = normal_ncxu;//((R[6]/P11)/normal_ncxu + l_ncxu_x*(R[6]/P11*(l_ncxu_x/normal_ncxu) + R[7]/P11*(l_ncxu_x/normal_ncxu)))
////                - ((pos_proj[0]/pos_proj[2] - R[6]/P11)/normal_ncxi + l_ncxi_x*((pos_proj[0]/pos_proj[2] - R[6]/P11)*(l_ncxi_x/normal_ncxi) + (pos_proj[1]/pos_proj[2] - R[7]/P11)*(l_ncxi_x/normal_ncxi)));
////      std::cout << "error 5 " << jacobian[3] << std::endl;

//    // 0,u
//    jacobian[0] = 0.0;
    }
//}
    return true;


  }

 private:
  ceres::internal::scoped_ptr<CostFunctor> functor_;
  Vec2 m_pos_2dpoint;
};


struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_4
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_4(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
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

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--
//#pragma omp critical
//      {
    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];


    T pos_proj[3];
    T X[3];
    X[0] = pos_3dpoint[0];
    X[1] = pos_3dpoint[1];
    X[2] = T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);

    //prepare R
    const T theta = sqrt(cam_R[0]*cam_R[0] + cam_R[1]*cam_R[1] + cam_R[2]*cam_R[2]);
    const T cosTheta = cos(theta);
    const T sinTheta = sin(theta);
    const T r[3] = {cam_R[0]/theta, cam_R[1]/theta, cam_R[2]/theta};
    const T one_cosTheta = T(1.0) - cosTheta;

    T R[9];
    R[0] = cosTheta + one_cosTheta*r[0]*r[0];
    R[1] = one_cosTheta*r[0]*r[1] + sinTheta*r[2];
    R[2] = one_cosTheta*r[0]*r[2] - sinTheta*r[1];
    R[3] = one_cosTheta*r[0]*r[1] - sinTheta*r[2];
    R[4] = cosTheta + one_cosTheta*r[1]*r[1];
    R[5] = one_cosTheta*r[1]*r[2] + sinTheta*r[0];
    R[6] = one_cosTheta*r[0]*r[2] + sinTheta*r[1];
    R[7] = one_cosTheta*r[1]*r[2] - sinTheta*r[0];
    R[8] = cosTheta + one_cosTheta*r[2]*r[2];

    const T& f = cam_intrinsics[OFFSET_FOCAL_LENGTH]; //focal
    const T& u = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
    const T& v = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y

    // Transform the point from homogeneous to euclidean (undistorted point)
//    const T x_i = u + f * (pos_proj[0] / pos_proj[2]);
//    const T y_i = v + f * (pos_proj[1] / pos_proj[2]);

    const T P1 = f*R[0] + u*R[2];
    const T P2 = f*R[3] + u*R[5];
    const T P3 = f*R[6] + u*R[8];
    const T P4 = f*cam_t[0] + u*cam_t[2];
    const T P5 = f*R[1] + v*R[2];
    const T P6 = f*R[4] + v*R[5];
    const T P7 = f*R[7] + v*R[8];
    const T P8 = f*cam_t[1] + v*cam_t[2];
    const T P9 = R[2];
    const T P10 = R[5];
    const T P11 = R[8];
    const T P12 = cam_t[2];

    const T depth_X = (P9*X[0] + P10*X[1] + P12);

    const T x_i = (P1*X[0] + P2*X[1] + P4) / (P9*X[0] + P10*X[1] + P12);
    const T y_i = (P5*X[0] + P6*X[1] + P8) / (P9*X[0] + P10*X[1] + P12);

    const T xu = T(m_pos_2dpoint[0]);
    const T yu = T(m_pos_2dpoint[1]);

    const T nc_x = P3/P11;
    const T nc_y = P7/P11;
    //
    const T l_ncxi_x = x_i - nc_x;
    const T l_ncxi_y = y_i - nc_y;
//    const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
    const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
    //
    const T l_ncxu_x = xu - nc_x;
    const T l_ncxu_y = yu - nc_y;
//    const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);
    const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);
    const T A = yu - nc_y;
    const T B = nc_x - xu;
    const T C = xu*nc_y - yu*nc_x;

    const T normalXi = (x_i-nc_x) / sqrt((x_i-nc_x)*(x_i-nc_x) + (y_i-nc_y)*(y_i-nc_y));
    const T normalYi = (y_i-nc_y) / sqrt((x_i-nc_x)*(x_i-nc_x) + (y_i-nc_y)*(y_i-nc_y));
    const T D = abs(A*x_i + B*y_i + C) / sqrt(A*A+B*B);
    ///
//    out_residuals[0] = D * normal_ncxu / T(2.0);
//    out_residuals[1] = ((l_ncxi_x*l_ncxu_x + l_ncxi_y*l_ncxu_y) / (normal_ncxi*normal_ncxu))* ((l_ncxi_x*l_ncxu_x + l_ncxi_y*l_ncxu_y) / (normal_ncxi*normal_ncxu));
//    out_residuals[1] = depth_X;

    {
//        Eigen::Matrix3Xd tempR, invR;
//        tempR << R[0], R[3], R[6],
//                 R[1], R[4], R[7],
//                 R[2], R[5], R[8];
//        invR = tempR.inverse();

        const T xu_cx = (xu - u)*depth_X/f - cam_t[0];
        const T xu_cy = (yu - v)*depth_X/f - cam_t[1];
        const T xu_cz = depth_X - cam_t[2];
        //
        const T xu_lx = R[0]*xu_cx + R[1]*xu_cy + R[2]*xu_cz;
        const T xu_ly = R[3]*xu_cx + R[4]*xu_cy + R[5]*xu_cz;
        const T xu_lz = R[6]*xu_cx + R[7]*xu_cy + R[8]*xu_cz;
//        const T xi_lx = (T)invR(0,0)*xi_cx + (T)invR(0,1)*xi_cy + (T)invR(0,2)*xi_cz;
//        const T xi_ly = (T)invR(1,0)*xi_cx + (T)invR(1,1)*xi_cy + (T)invR(1,2)*xi_cz;
//        const T xi_lz = (T)invR(2,0)*xi_cx + (T)invR(2,1)*xi_cy + (T)invR(2,2)*xi_cz;

        const T nc_cx = (nc_x - u)*depth_X/f - cam_t[0];
        const T nc_cy = (nc_y - v)*depth_X/f - cam_t[1];
        const T nc_cz = depth_X - cam_t[2];
        //
        const T nc_lx = R[0]*nc_cx + R[1]*nc_cy + R[2]*nc_cz;
        const T nc_ly = R[3]*nc_cx + R[4]*nc_cy + R[5]*nc_cz;
        const T nc_lz = R[6]*nc_cx + R[7]*nc_cy + R[8]*nc_cz;
//        const T nc_lx = (T)invR(0,0)*nc_cx + (T)invR(0,1)*nc_cy + (T)invR(0,2)*nc_cz;
//        const T nc_ly = (T)invR(1,0)*nc_cx + (T)invR(1,1)*nc_cy + (T)invR(1,2)*nc_cz;
//        const T nc_lz = (T)invR(2,0)*nc_cx + (T)invR(2,1)*nc_cy + (T)invR(2,2)*nc_cz;


        const T L_ncxu_x = nc_lx - xu_lx;
        const T L_ncxu_y = nc_ly - xu_ly;
        const T L_ncxu_z = nc_lz - xu_lz;

        const T L_Xxu_x = X[0] - xu_lx;
        const T L_Xxu_y = X[1] - xu_ly;
        const T L_Xxu_z = X[2] - xu_lz;

        const T cross_x = L_ncxu_y*L_Xxu_z - L_ncxu_z*L_Xxu_y;
        const T cross_y = L_ncxu_z*L_Xxu_x - L_ncxu_x*L_Xxu_z;
        const T cross_z = L_ncxu_x*L_Xxu_y - L_ncxu_y*L_Xxu_x;
        const T S = T(1.0/2.0) * sqrt((cross_x)*(cross_x) +(cross_y)*(cross_y) +(cross_z)*(cross_z));
//        std::cout << "S : " << S << std::endl;
        const T V = T(1.0/3.0) * S * depth_X;
//        std::cout << "V : " << V << std::endl;
//        std::cout << "d : " << depth_X << std::endl;
//        out_residuals[0] = V;
//        if(S == T(0.0))
//        {
//            out_residuals[0] = depth_X;
//        }
//        const T S2d = D * normal_ncxu / T(2.0);

//        out_residuals[1] = pow(sqrt(L_Xxu_x*L_Xxu_x + L_Xxu_y*L_Xxu_y),2);// + L_Xxu_z*L_Xxu_z), 3);

//        out_residuals[1] = S2d*S2d*S2d;



//        out_residuals[0] = V;
        out_residuals[0] = S;
        out_residuals[1] = (L_Xxu_x*L_Xxu_x);
        out_residuals[2] = (L_Xxu_y*L_Xxu_y);




//        out_residuals[0] = (l_ncxu_x*l_ncxu_x/normal_ncxu - l_ncxi_x*l_ncxi_x/normal_ncxi);
//        out_residuals[1] = (l_ncxu_y*l_ncxu_y/normal_ncxu - l_ncxi_y*l_ncxi_y/normal_ncxi);

//        out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi) * (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi);
//        out_residuals[1] = (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi) * (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);

//        out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi) * (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi)
//                         + (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi) * (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);

//        out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi) * (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi);
//        out_residuals[1] = (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi) * (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);

//        out_residuals[0] = (l_ncxu_x/normal_ncxu*l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi*l_ncxi_x/normal_ncxi)
//                         + (l_ncxu_y/normal_ncxu*l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi*l_ncxi_y/normal_ncxi);
    }

//}
    return true;
  }

  static int num_residuals() { return 3; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
//        (new ceres::NumericDiffCostFunction//ceres::AutoDiffCostFunction//
//           <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_4, ceres::RIDDERS, 1, 3, 6, 2>(
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_4, 3, 3, 6, 2>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_4(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_4>, 1, 3, 6, 2>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_4>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_4(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};

struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
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

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--
//#pragma omp critical
//      {
    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];


    T pos_proj[3];
    T X[3];
    X[0] = pos_3dpoint[0];
    X[1] = pos_3dpoint[1];
    X[2] = T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);
    // Rotate the point according the camera rotation
//    ceres::AngleAxisRotatePoint(cam_R, X, pos_proj);

//    //prepare R
    const T theta = sqrt(cam_R[0]*cam_R[0] + cam_R[1]*cam_R[1] + cam_R[2]*cam_R[2]);
    const T cosTheta = cos(theta);
    const T sinTheta = sin(theta);
    const T r[3] = {cam_R[0]/theta, cam_R[1]/theta, cam_R[2]/theta};
    const T one_cosTheta = T(1.0) - cosTheta;

    T R[9];
    R[0] = cosTheta + one_cosTheta*r[0]*r[0];
    R[1] = one_cosTheta*r[0]*r[1] + sinTheta*r[2];
    R[2] = one_cosTheta*r[0]*r[2] - sinTheta*r[1];
    R[3] = one_cosTheta*r[0]*r[1] - sinTheta*r[2];
    R[4] = cosTheta + one_cosTheta*r[1]*r[1];
    R[5] = one_cosTheta*r[1]*r[2] + sinTheta*r[0];
    R[6] = one_cosTheta*r[0]*r[2] + sinTheta*r[1];
    R[7] = one_cosTheta*r[1]*r[2] - sinTheta*r[0];
    R[8] = cosTheta + one_cosTheta*r[2]*r[2];

    const T& f = cam_intrinsics[OFFSET_FOCAL_LENGTH]; //focal
    const T& u = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
    const T& v = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y
    //
    const T P1 = f*R[0] + u*R[2];
    const T P2 = f*R[3] + u*R[5];
    const T P3 = f*R[6] + u*R[8];
    const T P4 = f*cam_t[0] + u*cam_t[2];
    const T P5 = f*R[1] + v*R[2];
    const T P6 = f*R[4] + v*R[5];
    const T P7 = f*R[7] + v*R[8];
    const T P8 = f*cam_t[1] + v*cam_t[2];
    const T P9 = R[2];
    const T P10 = R[5];
    const T P11 = R[8];
    const T P12 = cam_t[2];
    //
    const T depth_X = (P9*X[0] + P10*X[1] + P12);
    //
    const T x_i = (P1*X[0] + P2*X[1] + P3*X[2] + P4) / depth_X;//(P9*X[0] + P10*X[1] + P11*X[2] + P12);
    const T y_i = (P5*X[0] + P6*X[1] + P7*X[2] + P8) / depth_X;//(P9*X[0] + P10*X[1] + P11*X[2] + P12);
    //
    const T xu = T(m_pos_2dpoint[0]);
    const T yu = T(m_pos_2dpoint[1]);
    //
    const T nc_x = P3/P11;
    const T nc_y = P7/P11;
    //
//    const T l_ncxi_x = x_i - nc_x;
//    const T l_ncxi_y = y_i - nc_y;
//    const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
//    const T normal_ncxi = (l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
    //
    const T l_ncxu_x = xu - nc_x;
    const T l_ncxu_y = yu - nc_y;
    const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);
//    const T normal_ncxu = (l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);

    const T A = yu - nc_y;
    const T B = nc_x - xu;
    const T C = xu*nc_y - yu*nc_x;

    const T D = abs(A*x_i + B*y_i + C) / sqrt(A*A+B*B);
    ///
    {
//        out_residuals[0] = D * normal_ncxu / T(2.0);
//        out_residuals[1] = ((x_i-xu)*(x_i-xu) + (y_i-yu)*(y_i-yu));
    }

    {
        const T xu_cx = (xu - u)*depth_X/f - cam_t[0];
        const T xu_cy = (yu - v)*depth_X/f - cam_t[1];
        const T xu_cz = depth_X - cam_t[2];
        //
        const T xu_lx = R[0]*xu_cx + R[1]*xu_cy + R[2]*xu_cz;
        const T xu_ly = R[3]*xu_cx + R[4]*xu_cy + R[5]*xu_cz;
        const T xu_lz = R[6]*xu_cx + R[7]*xu_cy + R[8]*xu_cz;
//        const T xi_lx = (T)invR(0,0)*xi_cx + (T)invR(0,1)*xi_cy + (T)invR(0,2)*xi_cz;
//        const T xi_ly = (T)invR(1,0)*xi_cx + (T)invR(1,1)*xi_cy + (T)invR(1,2)*xi_cz;
//        const T xi_lz = (T)invR(2,0)*xi_cx + (T)invR(2,1)*xi_cy + (T)invR(2,2)*xi_cz;

        const T nc_cx = (nc_x - u)*depth_X/f - cam_t[0];
        const T nc_cy = (nc_y - v)*depth_X/f - cam_t[1];
        const T nc_cz = depth_X - cam_t[2];
        //
        const T nc_lx = R[0]*nc_cx + R[1]*nc_cy + R[2]*nc_cz;
        const T nc_ly = R[3]*nc_cx + R[4]*nc_cy + R[5]*nc_cz;
        const T nc_lz = R[6]*nc_cx + R[7]*nc_cy + R[8]*nc_cz;
//        const T nc_lx = (T)invR(0,0)*nc_cx + (T)invR(0,1)*nc_cy + (T)invR(0,2)*nc_cz;
//        const T nc_ly = (T)invR(1,0)*nc_cx + (T)invR(1,1)*nc_cy + (T)invR(1,2)*nc_cz;
//        const T nc_lz = (T)invR(2,0)*nc_cx + (T)invR(2,1)*nc_cy + (T)invR(2,2)*nc_cz;


        const T L_ncxu_x = nc_lx - xu_lx;
        const T L_ncxu_y = nc_ly - xu_ly;
        const T L_ncxu_z = nc_lz - xu_lz;

        const T L_Xxu_x = X[0] - xu_lx;
        const T L_Xxu_y = X[1] - xu_ly;
        const T L_Xxu_z = X[2] - xu_lz;

        const T cross_x = L_ncxu_y*L_Xxu_z - L_ncxu_z*L_Xxu_y;
        const T cross_y = L_ncxu_z*L_Xxu_x - L_ncxu_x*L_Xxu_z;
        const T cross_z = L_ncxu_x*L_Xxu_y - L_ncxu_y*L_Xxu_x;
        const T S = T(1.0/2.0) * sqrt((cross_x)*(cross_x) +(cross_y)*(cross_y) +(cross_z)*(cross_z));
//        std::cout << "S : " << S << std::endl;
        const T V = T(1.0/3.0) * S * depth_X;
//        std::cout << "V : " << V << std::endl;
//        std::cout << "d : " << depth_X << std::endl;
//        out_residuals[0] = V*V;
//        if(S == T(0.0))
        {
//            out_residuals[0] = depth_X;
        }
//        const T S2d = D * normal_ncxu / T(2.0);

//        out_residuals[1] = S2d*S2d*S2d;

//        out_residuals[0] = sqrt(S);//D * normal_ncxu / T(2.0);
//        out_residuals[1] = depth_X;

//        out_residuals[0] = V;//D * normal_ncxu / T(2.0);
//        out_residuals[1] = pow(sqrt(L_Xxu_x*L_Xxu_x + L_Xxu_y*L_Xxu_y),3);// + L_Xxu_z*L_Xxu_z), 3);

//        out_residuals[0] = V;
        out_residuals[0] = S;
//        out_residuals[1] = abs(L_Xxu_x*L_Xxu_x*L_Xxu_x);
//        out_residuals[2] = abs(L_Xxu_y*L_Xxu_y*L_Xxu_y);

//        const T l_ncxi_x = x_i - nc_x;
//        const T l_ncxi_y = y_i - nc_y;
//        const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);

//        out_residuals[0] = abs(l_ncxi_x/normal_ncxi * l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi * l_ncxu_x/normal_ncxu) / T(2.0);

    }
    return true;
  }

  static int num_residuals() { return 1; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
          (new ceres::NumericDiffCostFunction//ceres::AutoDiffCostFunction//
            <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6, ceres::CENTRAL, 1, 3, 6, 2>(
//        (new ceres::AutoDiffCostFunction
//          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6, 1, 3, 6, 2>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6(observation.data())));
    }
    else
    {
      return
        (new ceres::NumericDiffCostFunction//ceres::AutoDiffCostFunction
//          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3, ceres::RIDDERS, 1, 3, 6, 2>(
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6>, ceres::RIDDERS, 1, 3, 6, 2>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};


struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_5
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_5(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
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

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--
//#pragma omp critical
//      {
    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];


    T pos_proj[3];
    T X[3];
    X[0] = pos_3dpoint[0];
    X[1] = pos_3dpoint[1];
    X[2] = T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, X, pos_proj);

    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    //prepare l

    const T theta = sqrt(cam_R[0]*cam_R[0] + cam_R[1]*cam_R[1] + cam_R[2]*cam_R[2]);
    const T cosTheta = cos(theta);
    const T sinTheta = sin(theta);
    const T r[3] = {cam_R[0]/theta, cam_R[1]/theta, cam_R[2]/theta};
    const T one_cosTheta = T(1.0) - cosTheta;

    T R[9];
    R[0] = cosTheta + one_cosTheta*r[0]*r[0];
    R[1] = one_cosTheta*r[0]*r[1] + sinTheta*r[2];
    R[2] = one_cosTheta*r[0]*r[2] - sinTheta*r[1];
    R[3] = one_cosTheta*r[0]*r[1] - sinTheta*r[2];
    R[4] = cosTheta + one_cosTheta*r[1]*r[1];
    R[5] = one_cosTheta*r[1]*r[2] + sinTheta*r[0];
    R[6] = one_cosTheta*r[0]*r[2] + sinTheta*r[1];
    R[7] = one_cosTheta*r[1]*r[2] - sinTheta*r[0];
    R[8] = cosTheta + one_cosTheta*r[2]*r[2];

//    T R[9];
////    ceres::MatrixAdapter<T, 3, 3>(R);
//    ceres::AngleAxisToRotationMatrix(cam_R, R);

    const T& f = cam_intrinsics[OFFSET_FOCAL_LENGTH]; //focal
    const T& u = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
    const T& v = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y

    // Transform the point from homogeneous to euclidean (undistorted point)
//    const T depth_X = pos_proj[2];//
//    const T x_i = u + f * (pos_proj[0] / pos_proj[2]);
//    const T y_i = v + f * (pos_proj[1] / pos_proj[2]);

    const T P1 = f*R[0] + u*R[2];
    const T P2 = f*R[3] + u*R[5];
    const T P3 = f*R[6] + u*R[8];
    const T P4 = f*cam_t[0] + u*cam_t[2];
    const T P5 = f*R[1] + v*R[2];
    const T P6 = f*R[4] + v*R[5];
    const T P7 = f*R[7] + v*R[8];
    const T P8 = f*cam_t[1] + v*cam_t[2];
    const T P9 = R[2];
    const T P10 = R[5];
    const T P11 = R[8];
    const T P12 = cam_t[2];

    const T depth_X = (P9*X[0] + P10*X[1] + P12);
    const T x_i = (P1*X[0] + P2*X[1] + P3*X[2] + P4) / depth_X;
    const T y_i = (P5*X[0] + P6*X[1] + P7*X[2] + P8) / depth_X;

    const T xu = T(m_pos_2dpoint[0]);
    const T yu = T(m_pos_2dpoint[1]);
    //
    const T nc_x = P3/P11;
    const T nc_y = P7/P11;

    const T l_ncxu_x = xu - nc_x;
    const T l_ncxu_y = yu - nc_y;
    const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);
//    const T normal_ncxu = (l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);

    const T A = yu - nc_y;
    const T B = nc_x - xu;
    const T C = xu*nc_y - yu*nc_x;

    const T D = abs(A*x_i + B*y_i + C) / sqrt(A*A+B*B);
    ///
//    {
//        out_residuals[1] = D * normal_ncxu / T(2.0);

    {
        const T xu_cx = (xu - u)*depth_X/f - cam_t[0];
        const T xu_cy = (yu - v)*depth_X/f - cam_t[1];
        const T xu_cz = depth_X - cam_t[2];
        //
        const T xu_lx = R[0]*xu_cx + R[1]*xu_cy + R[2]*xu_cz;
        const T xu_ly = R[3]*xu_cx + R[4]*xu_cy + R[5]*xu_cz;
        const T xu_lz = R[6]*xu_cx + R[7]*xu_cy + R[8]*xu_cz;
    //        const T xi_lx = (T)invR(0,0)*xi_cx + (T)invR(0,1)*xi_cy + (T)invR(0,2)*xi_cz;
    //        const T xi_ly = (T)invR(1,0)*xi_cx + (T)invR(1,1)*xi_cy + (T)invR(1,2)*xi_cz;
    //        const T xi_lz = (T)invR(2,0)*xi_cx + (T)invR(2,1)*xi_cy + (T)invR(2,2)*xi_cz;

        const T nc_cx = (nc_x - u)*depth_X/f - cam_t[0];
        const T nc_cy = (nc_y - v)*depth_X/f - cam_t[1];
        const T nc_cz = depth_X - cam_t[2];
        //
        const T nc_lx = R[0]*nc_cx + R[1]*nc_cy + R[2]*nc_cz;
        const T nc_ly = R[3]*nc_cx + R[4]*nc_cy + R[5]*nc_cz;
        const T nc_lz = R[6]*nc_cx + R[7]*nc_cy + R[8]*nc_cz;
    //        const T nc_lx = (T)invR(0,0)*nc_cx + (T)invR(0,1)*nc_cy + (T)invR(0,2)*nc_cz;
    //        const T nc_ly = (T)invR(1,0)*nc_cx + (T)invR(1,1)*nc_cy + (T)invR(1,2)*nc_cz;
    //        const T nc_lz = (T)invR(2,0)*nc_cx + (T)invR(2,1)*nc_cy + (T)invR(2,2)*nc_cz;


        const T L_ncxu_x = nc_lx - xu_lx;
        const T L_ncxu_y = nc_ly - xu_ly;
        const T L_ncxu_z = nc_lz - xu_lz;

        const T L_Xxu_x = X[0] - xu_lx;
        const T L_Xxu_y = X[1] - xu_ly;
        const T L_Xxu_z = X[2] - xu_lz;

        const T cross_x = L_ncxu_y*L_Xxu_z - L_ncxu_z*L_Xxu_y;
        const T cross_y = L_ncxu_z*L_Xxu_x - L_ncxu_x*L_Xxu_z;
        const T cross_z = L_ncxu_x*L_Xxu_y - L_ncxu_y*L_Xxu_x;
        const T S = T(1.0/2.0) * sqrt((cross_x)*(cross_x) +(cross_y)*(cross_y) +(cross_z)*(cross_z));
    //        std::cout << "S : " << S << std::endl;
        const T V = T(1.0/3.0) * S * depth_X;
//        if(S == T(0.0))
        {
//            out_residuals[0] = depth_X;
        }
//        out_residuals[0] = pow(V,1.0/3.0);

//        out_residuals[0] = V;
//        out_residuals[1] = abs(L_Xxu_x*L_Xxu_x*L_Xxu_x);
//        out_residuals[2] = abs(L_Xxu_y*L_Xxu_y*L_Xxu_y);

//        out_residuals[0] = V;
        out_residuals[0] = S;
        out_residuals[1] = (L_Xxu_x*L_Xxu_x);
        out_residuals[2] = (L_Xxu_y*L_Xxu_y);

//        out_residuals[0] = abs(L_Xxu_x*L_Xxu_x*L_Xxu_x / V);
//        out_residuals[1] = abs(L_Xxu_y*L_Xxu_y*L_Xxu_y / V);
//        {
//            std::cout <<"*****************************" << std::endl;
//            std::cout << "X[0] " << X[0] << std::endl;
//            std::cout << "X[1] " << X[1] << std::endl;
//            std::cout << "X[2] " << X[2] << std::endl;
//            std::cout << "xu_lx  " << xu_lx << std::endl;
//            std::cout << "xu_ly  " << xu_ly << std::endl;
//            std::cout << "xu_lz  " << xu_lz << std::endl;
//            std::cout << "nc_lx " << nc_lx << std::endl;
//            std::cout << "nc_ly " << nc_ly << std::endl;
//            std::cout << "nc_lz " << nc_lz << std::endl;
//            std::cout << "L_Xxu_x " << L_Xxu_x << std::endl;
//            std::cout << "L_Xxu_y " << L_Xxu_y << std::endl;
//            std::cout << "L_Xxu_z " << L_Xxu_z << std::endl;
//            std::cout << "L_ncxu_x " << L_ncxu_x << std::endl;
//            std::cout << "L_ncxu_y " << L_ncxu_y << std::endl;
//            std::cout << "L_ncxu_z " << L_ncxu_z << std::endl;
//            std::cout << "depth_X " << depth_X << std::endl;
//            std::cout << "S " << S << std::endl;
//            std::cout << "V " << V << std::endl;
//            std::cout << "cross_x " << cross_x << std::endl;
//            std::cout << "cross_y " << cross_y << std::endl;
//            std::cout << "cross_z " << cross_z << std::endl;
//        }

//        out_residuals[1] = -L_Xxu_x*L_Xxu_x*L_Xxu_x;//pow(sqrt(L_Xxu_x*L_Xxu_x + L_Xxu_y*L_Xxu_y + L_Xxu_z*L_Xxu_z),3);// + L_Xxu_z*L_Xxu_z), 3);
//        out_residuals[2] = -L_Xxu_y*L_Xxu_y*L_Xxu_y;//pow(sqrt(L_Xxu_x*L_Xxu_x + L_Xxu_y*L_Xxu_y + L_Xxu_z*L_Xxu_z),3);// + L_Xxu_z*L_Xxu_z), 3);


//        out_residuals[0] = depth_X;//S;
//        out_residuals[1] = sqrt(S);
//        out_residuals[0] = L_Xxu_x;//S;

//        //plane
//        const T C_x = -(R[0]*cam_t[0] + R[1]*cam_t[1] + R[2]*cam_t[2]);
//        const T C_y = -(R[3]*cam_t[0] + R[4]*cam_t[1] + R[5]*cam_t[2]);
//        const T C_z = -(R[6]*cam_t[0] + R[7]*cam_t[1] + R[8]*cam_t[2]);
//        const T L_ncC_x = nc_lx - C_x;
//        const T L_ncC_y = nc_ly - C_y;
//        const T L_ncC_z = nc_lz - C_z;
//        const T L_xuC_x = xu_lx - C_x;
//        const T L_xuC_y = xu_ly - C_y;
//        const T L_xuC_z = xu_lz - C_z;
//        const T Plane_A = L_ncC_y*L_xuC_z - L_ncC_z*L_xuC_y;
//        const T Plane_B = L_ncC_z*L_xuC_x - L_ncC_x*L_xuC_z;
//        const T Plane_C = L_ncC_x*L_xuC_y - L_ncC_y*L_xuC_x;
//        const T Plane_D = -(cross_x*xu_lx + cross_y*xu_ly + cross_z*xu_lz);
//        const T D3d = abs(Plane_A*X[0] + Plane_B*X[1] + Plane_C*X[2] + Plane_D) / sqrt(Plane_A*Plane_A + Plane_B*Plane_B + Plane_C*Plane_C);
//        out_residuals[1] = D3d * D3d * D3d;

//        out_residuals[0] = D3d;// / depth_X;//* D3d * D3d / V;



//        out_residuals[2] = L_Xxu_z;
    }
//}
    return true;
  }

  static int num_residuals() { return 3; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
//        (new ceres::NumericDiffCostFunction//ceres::AutoDiffCostFunction//
//          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_5, ceres::RIDDERS, 1, 3, 6, 2>(
        (new ceres::AutoDiffCostFunction//
          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_5, 3, 3, 6, 2>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_5(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_5>, 2, 3, 6, 2>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_5>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_5(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};


struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_5_old
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_5_old(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
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

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--
//#pragma omp critical
//      {
    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];


    T pos_proj[3];
    T X[3];
    X[0] = pos_3dpoint[0];
    X[1] = pos_3dpoint[1];
    X[2] = T(0.0);//pos_3dpoint[2];//T(0.0);
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, X, pos_proj);

    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    //prepare l
    T R[9];
//    ceres::MatrixAdapter<T, 3, 3>(R);
    ceres::AngleAxisToRotationMatrix(cam_R, R);

    const T& f = cam_intrinsics[OFFSET_FOCAL_LENGTH]; //focal
    const T& u = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
    const T& v = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y

    // Transform the point from homogeneous to euclidean (undistorted point)
    const T x_i = u + f * (pos_proj[0] / pos_proj[2]);
    const T y_i = v + f * (pos_proj[1] / pos_proj[2]);

    const T P1 = f*R[0] + u*R[2];
    const T P2 = f*R[3] + u*R[5];
    const T P3 = f*R[6] + u*R[8];
    const T P4 = f*cam_t[0] + u*cam_t[2];
    const T P5 = f*R[1] + v*R[2];
    const T P6 = f*R[4] + v*R[5];
    const T P7 = f*R[7] + v*R[8];
    const T P8 = f*cam_t[1] + v*cam_t[2];
    const T P9 = R[2];
    const T P10 = R[5];
    const T P11 = R[8];
    const T P12 = cam_t[2];

    const T xu = T(m_pos_2dpoint[0]);
    const T yu = T(m_pos_2dpoint[1]);

//    const T l0 = P7 - P11*yu;
//    const T l1 = P11*xu - P3;
//    const T l2 = P3*yu - P7*xu;
    const T l0 = y_i*P11 - P7;//P7 - P11*yu;
    const T l1 = P3 - P11*x_i;//P11*xu - P3;
    const T l2 = P7*x_i - P3*y_i;//P3*yu - P7*xu;
    out_residuals[0] = abs(xu*l0 + yu*l1 + l2)/sqrt(l0*l0 + l1*l1);

    const T l0_ = P7 - P11*yu;
    const T l1_ = P11*xu - P3;
    const T l2_ = P3*yu - P7*xu;
    out_residuals[1] = abs(x_i*l0_ + y_i*l1_ + l2_)/sqrt(l0_*l0_ + l1_*l1_);

    out_residuals[2] = sqrt((x_i-P3/P11)*(x_i-P3/P11) + (y_i-P7/P11)*(y_i-P7/P11));//abs(x_i*l0_ + y_i*l1_ + l2_)/sqrt(l0_*l0_ + l1_*l1_);

//    if(l0 != T(0.0))
//    {
//        out_residuals[0] = (- y_i*l1 - l2)/l0 - x_i;

//    }else{
//        out_residuals[0] = T(0.0);
//    }
//    if(l1 != T(0.0))
//    {
//        out_residuals[1] = ( -x_i*l0 - l2)/l1 - y_i;
//    }else{
//        out_residuals[1] = T(0.0);

//    }
//    out_residuals[0] = ( x_i*l0 + y_i*l1 + l2);
//    out_residuals[0] = (x_i*l0 + y_i*l1 + l2)/sqrt(l0*l0 + l1*l1);

//    out_residuals[0] = (x_i*l0 + y_i*l1 + l2)/sqrt(l0*l0 + l1*l1);
//    out_residuals[1] = T(0.0);
//    out_residuals[0] = abs(x_i*l0 + y_i*l1 + l2)/sqrt(l0*l0 + l1*l1);

//    if(out_residuals[0] > T(0.1))
//    {
//        std::cout << "cam_R : " << cam_R[0] << " " << cam_R[1] << " " << cam_R[2] << std::endl;
//        std::cout << "R : " << R[0] << " " << R[1] << " " << R[2]
//                            << R[3] << " " << R[4] << " " << R[5]
//                            << R[6] << " " << R[7] << " " << R[8] << std::endl;
//        std::cout << "cam_t : " << cam_t[0] << " " << cam_t[1] << " " << cam_t[2] << std::endl;
//        std::cout << "f : " << f << " u : " << u << " v : " << v << std::endl;
//        std::cout << "xu : " << xu << " yu : " << yu << std::endl;
//        std::cout << "pos_3dpoint[0]: " << pos_3dpoint[0] << " pos_3dpoint[1]: " << pos_3dpoint[1] << " pos_3dpoint[2]: " << pos_3dpoint[2] << std::endl;
//        std::cout << "out_residuals[0]: "  << out_residuals[0] << std::endl;
//        std::cout << "****-----------------****" << std::endl;
//    }

//      }

//    T pos_proj[3];
//    // Rotate the point according the camera rotation
//    ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_proj);

//    // Apply the camera translation
//    pos_proj[0] += cam_t[0];
//    pos_proj[1] += cam_t[1];
//    pos_proj[2] += cam_t[2];

//    // Transform the point from homogeneous to euclidean (undistorted point)
//    const T x_u = pos_proj[0] / pos_proj[2];
//    const T y_u = pos_proj[1] / pos_proj[2];

//    //--
//    // Apply intrinsic parameters
//    //--


//    const T& k1 = cam_intrinsics[OFFSET_DISTO_K1];
//    const T& k2 = cam_intrinsics[OFFSET_DISTO_K2];
//    const T& k3 = cam_intrinsics[OFFSET_DISTO_K3];
//    const T& t1 = cam_intrinsics[OFFSET_DISTO_T1];
//    const T& t2 = cam_intrinsics[OFFSET_DISTO_T2];

//    // Apply distortion (xd,yd) = disto(x_u,y_u)
//    const T r2 = x_u*x_u + y_u*y_u;
//    const T r4 = r2 * r2;
//    const T r6 = r4 * r2;
//    const T r_coeff = (1.0 + k1*r2 + k2*r4 + k3*r6);
//    const T t_x = t2 * (r2 + 2.0 * x_u*x_u) + 2.0 * t1 * x_u * y_u;
//    const T t_y = t1 * (r2 + 2.0 * y_u*y_u) + 2.0 * t2 * x_u * y_u;
//    const T x_d = x_u * r_coeff + t_x;
//    const T y_d = y_u * r_coeff + t_y;

//    // Apply focal length and principal point to get the final image coordinates
//    const T projected_x = principal_point_x + focal * x_d;
//    const T projected_y = principal_point_y + focal * y_d;

//    // Compute and return the error is the difference between the predicted
//    //  and observed position
//    out_residuals[0] = projected_x - m_pos_2dpoint[0];
//    out_residuals[1] = projected_y - m_pos_2dpoint[1];

    return true;
  }

  static int num_residuals() { return 3; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_5_old, 3, 3, 6, 2>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_5_old(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_5_old>, 3, 3, 6, 2>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_5_old>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_5_old(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};

struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6_old
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6_old(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
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

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--
//#pragma omp critical
//      {
    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];


    T pos_proj[3];
    T X[3];
    X[0] = pos_3dpoint[0];
    X[1] = pos_3dpoint[1];
    X[2] = T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, X, pos_proj);

    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    //prepare l
    T R[9];
//    ceres::MatrixAdapter<T, 3, 3>(R);
    ceres::AngleAxisToRotationMatrix(cam_R, R);

    const T& f = cam_intrinsics[OFFSET_FOCAL_LENGTH]; //focal
    const T& u = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
    const T& v = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y

    // Transform the point from homogeneous to euclidean (undistorted point)
    const T x_i = u + f * (pos_proj[0] / pos_proj[2]);
    const T y_i = v + f * (pos_proj[1] / pos_proj[2]);

    const T P1 = f*R[0] + u*R[2];
    const T P2 = f*R[3] + u*R[5];
    const T P3 = f*R[6] + u*R[8];
    const T P4 = f*cam_t[0] + u*cam_t[2];
    const T P5 = f*R[1] + v*R[2];
    const T P6 = f*R[4] + v*R[5];
    const T P7 = f*R[7] + v*R[8];
    const T P8 = f*cam_t[1] + v*cam_t[2];
    const T P9 = R[2];
    const T P10 = R[5];
    const T P11 = R[8];
    const T P12 = cam_t[2];

    const T xu = T(m_pos_2dpoint[0]);
    const T yu = T(m_pos_2dpoint[1]);


    const T l0 = y_i*P11 - P7;//P7 - P11*yu;
    const T l1 = P3 - P11*x_i;//P11*xu - P3;
    const T l2 = P7*x_i - P3*y_i;//P3*yu - P7*xu;
    out_residuals[0] = abs(xu*l0 + yu*l1 + l2)/sqrt(l0*l0 + l1*l1);

    const T l0_ = P7 - P11*yu;
    const T l1_ = P11*xu - P3;
    const T l2_ = P3*yu - P7*xu;
    out_residuals[1] = abs(x_i*l0_ + y_i*l1_ + l2_)/sqrt(l0_*l0_ + l1_*l1_);

//    const T nc0 = P1*pos_3dpoint[0] + P2*pos_3dpoint[1];
//    const T nc1 = P5*pos_3dpoint[0] + P6*pos_3dpoint[1];
//    const T nc2 = P9*pos_3dpoint[0] + P10*pos_3dpoint[1];

//    const T l0 = y_i*nc2 - nc1;//P7 - P11*yu;
//    const T l1 = nc0 - nc2*x_i;//P11*xu - P3;
//    const T l2 = nc1*x_i - nc0*y_i;//P3*yu - P7*xu;
//    out_residuals[0] = abs(xu*l0 + yu*l1 + l2)/sqrt(l0*l0 + l1*l1);

//    const T l0_ = nc1 - nc2*yu;
//    const T l1_ = nc2*xu - nc0;
//    const T l2_ = nc0*yu - nc2*xu;
//    out_residuals[1] = abs(x_i*l0_ + y_i*l1_ + l2_)/sqrt(l0_*l0_ + l1_*l1_);

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6_old, 2, 3, 6, 2>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6_old(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6_old>, 2, 3, 6, 2>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6_old>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6_old(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};

//struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6
//{
//  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6(const double* const pos_2dpoint)
//  :m_pos_2dpoint(pos_2dpoint)
//  {
//  }

//  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
//  enum : uint8_t {
//    OFFSET_FOCAL_LENGTH = 0,
//    OFFSET_PRINCIPAL_POINT_X = 1,
//    OFFSET_PRINCIPAL_POINT_Y = 2,
//    OFFSET_DISTO_K1 = 3,
//    OFFSET_DISTO_K2 = 4,
//    OFFSET_DISTO_K3 = 5,
//    OFFSET_DISTO_T1 = 6,
//    OFFSET_DISTO_T2 = 7,
//  };

//  /**
//   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
//   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
//   *   - 3 for rotation(angle axis), 3 for translation
//   * @param[in] pos_3dpoint
//   * @param[out] out_residuals
//   */
//  template <typename T>
//  bool operator()(
//    const T* const cam_intrinsics,
//    const T* const cam_extrinsics,
//    const T* const pos_3dpoint,
//    T* out_residuals) const
//  {
//    //--
//    // Apply external parameters (Pose)
//    //--
////#pragma omp critical
////      {
//    const T * cam_R = cam_extrinsics;
//    const T * cam_t = &cam_extrinsics[3];


//    T pos_proj[3];
//    T X[3];
//    X[0] = pos_3dpoint[0];
//    X[1] = pos_3dpoint[1];
//    X[2] = T(-1200000.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);
//    // Rotate the point according the camera rotation
//    ceres::AngleAxisRotatePoint(cam_R, X, pos_proj);

////    // Apply the camera translation
//    pos_proj[0] += cam_t[0];
//    pos_proj[1] += cam_t[1];
//    pos_proj[2] += cam_t[2];

////    //prepare l
//    T R[9];
////    ceres::MatrixAdapter<T, 3, 3>(R);
//    ceres::AngleAxisToRotationMatrix(cam_R, R);

//    const T& f = cam_intrinsics[OFFSET_FOCAL_LENGTH]; //focal
//    const T& u = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
//    const T& v = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y

//    // Transform the point from homogeneous to euclidean (undistorted point)
//    const T x_i = u + f * (pos_proj[0] / pos_proj[2]);
//    const T y_i = v + f * (pos_proj[1] / pos_proj[2]);

//    const T P1 = f*R[0] + u*R[2];
//    const T P2 = f*R[3] + u*R[5];
//    const T P3 = f*R[6] + u*R[8];
//    const T P4 = f*cam_t[0] + u*cam_t[2];
//    const T P5 = f*R[1] + v*R[2];
//    const T P6 = f*R[4] + v*R[5];
//    const T P7 = f*R[7] + v*R[8];
//    const T P8 = f*cam_t[1] + v*cam_t[2];
//    const T P9 = R[2];
//    const T P10 = R[5];
//    const T P11 = R[8];
//    const T P12 = cam_t[2];

//    {//for nc_ P11
////        const T xu = T(m_pos_2dpoint[0]);
////        const T yu = T(m_pos_2dpoint[1]);

////        const T nc_x = P3;
////        const T nc_y = P7;
////        //P11
////    //    const T nc_x = P3;
////    //    const T nc_y = P7;
////        const T A = yu*P11*P11 - nc_y*P11;
////        const T B = nc_x*P11 - xu*P11*P11;
////        const T C = xu*P11*nc_y - yu*P11*nc_x;

////        const T D = abs(A*x_i + B*y_i + C) / sqrt(A*A+B*B);

////    //    out_residuals[0] = (xu - x_i) / D;
////    //    out_residuals[1] = (yu - y_i) / D;

////        //
////        const T l_ncxi_x = x_i*P11 - nc_x;
////        const T l_ncxi_y = y_i*P11 - nc_y;
//////        const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
////        const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
////    //    const T normal_ncxi = (l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
////        //
////        const T l_ncxu_x = xu*P11 - nc_x;
////        const T l_ncxu_y = yu*P11 - nc_y;
//////        const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);
////        const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);

//////        out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi);
//////        out_residuals[1] = (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);

//////        const T l_nc = sqrt(normal_ncxi*normal_ncxi - D*D);//need change
////        const T p_x = (B*(A*y_i-B*x_i)+A*C) / (-A*A-B*B);
////        const T p_y = (B*C-A*(A*y_i-B*x_i)) / (-A*A-B*B);
////        const T l_nc = sqrt((p_x-nc_x/P11)*(p_x-nc_x/P11) + (p_y-nc_y/P11)*(p_y-nc_y/P11));
////        out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/l_nc);
////        out_residuals[1] = (l_ncxu_y/normal_ncxu - l_ncxi_y/l_nc);

//    }

////    const T x_i = (P1*X[0] + P2*X[1] + P3*X[2] + P4) / (P9*X[0] + P10*X[1] + P11*X[2] + P12);
////    const T y_i = (P5*X[0] + P6*X[1] + P7*X[2] + P8) / (P9*X[0] + P10*X[1] + P11*X[2] + P12);


//    const T xu = T(m_pos_2dpoint[0]);
//    const T yu = T(m_pos_2dpoint[1]);

//    const T nc_x = P3/P11;
//    const T nc_y = P7/P11;
////    const T nc_x = P3;
////    const T nc_y = P7;
//    const T A = yu - nc_y;
//    const T B = nc_x - xu;
//    const T C = xu*nc_y - yu*nc_x;

//    const T D = abs(A*x_i + B*y_i + C) / sqrt(A*A+B*B);

////    out_residuals[0] = (xu - x_i) / D;
////    out_residuals[1] = (yu - y_i) / D;

//    //
//    const T l_ncxi_x = x_i - nc_x;
//    const T l_ncxi_y = y_i - nc_y;
//    const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
////    const T normal_ncxi = (l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
//    //
//    const T l_ncxu_x = xu - nc_x;
//    const T l_ncxu_y = yu - nc_y;
//    const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);

//    out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi);
//    out_residuals[1] = (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);

////    if(X[2] != T(0.0))
//    {
////        std::cout << "***----------------------***" << std::endl;
////        std::cout << "R : " << cam_R[0] << " " << cam_R[1] << " " << cam_R[2] << std::endl;
////        std::cout << "t : " << cam_t[0] << " " << cam_t[1] << " " << cam_t[2] << std::endl;
////        std::cout << "f : " << f << " u : " << u << " v : " << v << std::endl;
////        std::cout << "x_i : " << x_i << " y_i : " << y_i << std::endl;
////        std::cout << "xu : " << xu << " yu : " << yu << std::endl;
////        std::cout << "nc_x : " << nc_x << " nc_y : " << nc_y << std::endl;
////        std::cout << "l_ncxi_x : " << l_ncxi_x << " l_ncxi_y : " << l_ncxi_y << std::endl;
////        std::cout << "l_ncxu_x : " << l_ncxu_x << " l_ncxu_y : " << l_ncxu_y << std::endl;
////        std::cout << "normal_ncxi : " << normal_ncxi << " normal_ncxu : " << normal_ncxu << std::endl;
////        std::cout << "out_residuals[0] : " << out_residuals[0] << " out_residuals[1] : " << out_residuals[1] << std::endl;
////        std::cout << "X : " << X[0] << " " << X[1] << " " << X[2] << std::endl;
////        std::cout << "----------" << std::endl;

////        {
////        T pos_proj[3];
////        T X[3];
////        X[0] = pos_3dpoint[0];
////        X[1] = pos_3dpoint[1];
////        X[2] = T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);
////        // Rotate the point according the camera rotation
////        ceres::AngleAxisRotatePoint(cam_R, X, pos_proj);

////    //    // Apply the camera translation
////        pos_proj[0] += cam_t[0];
////        pos_proj[1] += cam_t[1];
////        pos_proj[2] += cam_t[2];

////    //    //prepare l
////        T R[9];
////    //    ceres::MatrixAdapter<T, 3, 3>(R);
////        ceres::AngleAxisToRotationMatrix(cam_R, R);

////        const T& f = cam_intrinsics[OFFSET_FOCAL_LENGTH]; //focal
////        const T& u = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
////        const T& v = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y

////        // Transform the point from homogeneous to euclidean (undistorted point)
////        const T x_i = u + f * (pos_proj[0] / pos_proj[2]);
////        const T y_i = v + f * (pos_proj[1] / pos_proj[2]);

////        const T P1 = f*R[0] + u*R[2];
////        const T P2 = f*R[3] + u*R[5];
////        const T P3 = f*R[6] + u*R[8];
////        const T P4 = f*cam_t[0] + u*cam_t[2];
////        const T P5 = f*R[1] + v*R[2];
////        const T P6 = f*R[4] + v*R[5];
////        const T P7 = f*R[7] + v*R[8];
////        const T P8 = f*cam_t[1] + v*cam_t[2];
////        const T P9 = R[2];
////        const T P10 = R[5];
////        const T P11 = R[8];
////        const T P12 = cam_t[2];

////        const T xu = T(m_pos_2dpoint[0]);
////        const T yu = T(m_pos_2dpoint[1]);

////        const T nc_x = P3/P11;
////        const T nc_y = P7/P11;
////    //    const T nc_x = P3;
////    //    const T nc_y = P7;
////        const T A = yu - nc_y;
////        const T B = nc_x - xu;
////        const T C = xu*nc_y - yu*nc_x;

////        const T D = abs(A*x_i + B*y_i + C) / sqrt(A*A+B*B);

////    //    out_residuals[0] = (xu - x_i) / D;
////    //    out_residuals[1] = (yu - y_i) / D;

////        //
////        const T l_ncxi_x = x_i - nc_x;
////        const T l_ncxi_y = y_i - nc_y;
////        const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
////    //    const T normal_ncxi = (l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
////        //
////        const T l_ncxu_x = xu - nc_x;
////        const T l_ncxu_y = yu - nc_y;
////        const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);

////        out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi);
////        out_residuals[1] = (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);

////        std::cout << "R : " << cam_R[0] << " " << cam_R[1] << " " << cam_R[2] << std::endl;
////        std::cout << "t : " << cam_t[0] << " " << cam_t[1] << " " << cam_t[2] << std::endl;
////        std::cout << "f : " << f << " u : " << u << " v : " << v << std::endl;
////        std::cout << "x_i : " << x_i << " y_i : " << y_i << std::endl;
////        std::cout << "xu : " << xu << " yu : " << yu << std::endl;
////        std::cout << "nc_x : " << nc_x << " nc_y : " << nc_y << std::endl;
////        std::cout << "l_ncxi_x : " << l_ncxi_x << " l_ncxi_y : " << l_ncxi_y << std::endl;
////        std::cout << "l_ncxu_x : " << l_ncxu_x << " l_ncxu_y : " << l_ncxu_y << std::endl;
////        std::cout << "normal_ncxi : " << normal_ncxi << " normal_ncxu : " << normal_ncxu << std::endl;
////        std::cout << "out_residuals[0] : " << out_residuals[0] << " out_residuals[1] : " << out_residuals[1] << std::endl;
////        std::cout << "X : " << X[0] << " " << X[1] << " " << X[2] << std::endl;
////        }
//    }
////    const T normal_ncxu = (l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);
//    ///
////    const T x2 = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi)*(l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi);
////    const T y2 = (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi)*(l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);
//    {
////        const T l_nc = sqrt(normal_ncxi*normal_ncxi - D*D);
////        out_residuals[0] = D*D / (l_nc*l_nc);

////        out_residuals[0] = D;///*D / (l_nc*l_nc);
////        out_residuals[1] = sqrt((x_i-xu)*(x_i-xu)+(y_i-yu)*(y_i-yu));///*D / (l_nc*l_nc);

////        const T l_nc = sqrt(normal_ncxi*normal_ncxi - D*D);//need change
////        const T p_x = (-l_ncxu_x - nc_y*(l_ncxu_x*y_i - l_ncxu_y*x_i)) / (l_ncxu_y*nc_y + l_ncxu_x*nc_x);
////        const T p_y = (nc_x*(l_ncxu_x*y_i - l_ncxu_y*x_i) - l_ncxu_y) / (l_ncxu_y*nc_y + l_ncxu_x*nc_x);
////        const T l_nc = sqrt((p_x-nc_x)*(p_x-nc_x) + (p_y-nc_y)*(p_y-nc_y));
////        out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/l_nc);
////        out_residuals[1] = (l_ncxu_y/normal_ncxu - l_ncxi_y/l_nc);

//        { // can be used
////        const T l1_x = x_i - xu;
////        const T l1_y = y_i - yu;
////        const T l2_x = xu - nc_x;
////        const T l2_y = yu - nc_y;
////        const T norml1 = sqrt(l1_x*l1_x + l1_y*l1_y);
////        const T norml2 = sqrt(l2_x*l2_x + l2_y*l2_y);
////        out_residuals[0] = -(l1_x/ norml1-l2_x/ norml2);
////        out_residuals[1] = -(l1_y/ norml1-l2_y/ norml2);
//        }

////        T l_ncxu_x2, l_ncxi_x2;
////        if(l_ncxu_x < T(0.0))
////        {
////            l_ncxu_x2 = -l_ncxu_x*l_ncxu_x;
////        }else{
////            l_ncxu_x2 = l_ncxu_x*l_ncxu_x;
////        }
////        if(l_ncxi_x < T(0.0))
////        {
////            l_ncxi_x2 = -l_ncxi_x*l_ncxi_x;
////        }else{
////            l_ncxi_x2 = l_ncxi_x*l_ncxi_x;
////        }
////        //
////        T l_ncxu_y2, l_ncxi_y2;
////        if(l_ncxu_y < T(0.0))
////        {
////            l_ncxu_y2 = -l_ncxu_y*l_ncxu_y;
////        }else{
////            l_ncxu_y2 = l_ncxu_y*l_ncxu_y;
////        }
////        if(l_ncxi_y < T(0.0))
////        {
////            l_ncxi_y2 = -l_ncxi_y*l_ncxi_y;
////        }else{
////            l_ncxi_y2 = l_ncxi_y*l_ncxi_y;
////        }
////        //
////        out_residuals[0] = (l_ncxu_x2/normal_ncxu - l_ncxi_x2/normal_ncxi);
////        out_residuals[1] = (l_ncxu_y2/normal_ncxu - l_ncxi_y2/normal_ncxi);

////        out_residuals[0] = (l_ncxu_x/(l_ncxi_x*l_ncxu_y - l_ncxi_y*l_ncxu_x) - l_ncxi_x/(l_ncxi_x*l_ncxu_y - l_ncxi_y*l_ncxu_x));
////        out_residuals[1] = (l_ncxu_y/(l_ncxi_x*l_ncxu_y - l_ncxi_y*l_ncxu_x) - l_ncxi_y/(l_ncxi_x*l_ncxu_y - l_ncxi_y*l_ncxu_x));

////        const T dx = xu - x_i;//l_ncxu_x - l_ncxi_x;
////        const T dy = yu - y_i;//l_ncxu_y - l_ncxi_y;
////        out_residuals[0] = dx / sqrt(dx*dx + dy*dy);
////        out_residuals[1] = dy / sqrt(dx*dx + dy*dy);


////        const T dx = l_ncxu_x*l_ncxu_x*normal_ncxi - l_ncxi_x*l_ncxi_x*normal_ncxu;
////        const T dy = l_ncxu_y*l_ncxu_y*normal_ncxi - l_ncxi_y*l_ncxi_y*normal_ncxu;
////        out_residuals[0] = dx*dx / (dx*dx + dy*dy);
////        out_residuals[1] = dy*dy / (dx*dx + dy*dy);


////        const T dx = xu*xu*normal_ncxi - x_i*x_i*normal_ncxu;
////        const T dy = yu*yu*normal_ncxi - y_i*y_i*normal_ncxu;
////        out_residuals[0] = dx*dx / (dx*dx + dy*dy);
////        out_residuals[1] = dy*dy / (dx*dx + dy*dy);

////        out_residuals[0] = dx*dx / (dx*dx + dy*dy);
////        out_residuals[1] = dy*dy / (dx*dx + dy*dy);

////        out_residuals[0] = (l_ncxu_x - l_ncxi_x)/(l_ncxi_x*l_ncxu_x - l_ncxi_y*l_ncxu_y);
////        out_residuals[1] = (l_ncxu_y - l_ncxi_y)/(l_ncxi_x*l_ncxu_x - l_ncxi_y*l_ncxu_y);


////        out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi) * (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi);
////        out_residuals[1] = (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi) * (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);

////        out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi) * (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi)
////                         + (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi) * (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);

////        out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi) * (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi);
////        out_residuals[2] = (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi) * (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);

////        out_residuals[0] = (l_ncxu_x/normal_ncxu*l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi*l_ncxi_x/normal_ncxi)
////                         + (l_ncxu_y/normal_ncxu*l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi*l_ncxi_y/normal_ncxi);
//    }

////      }
//    return true;
//  }

//  static int num_residuals() { return 2; }

//  // Factory to hide the construction of the CostFunction object from
//  // the client code.
//  static ceres::CostFunction* Create
//  (
//    const Vec2 & observation,
//    const double weight = 0.0
//  )
//  {
//    if (weight == 0.0)
//    {
//      return
//          (new ceres::NumericDiffCostFunction//ceres::AutoDiffCostFunction//
//            <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6, ceres::FORWARD, 2, 3, 6, 2>(
////        (new ceres::AutoDiffCostFunction
////          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6, 2, 3, 6, 2>(
//            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6(observation.data())));
//    }
//    else
//    {
//      return
//        (new ceres::AutoDiffCostFunction
//          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6>, 2, 3, 6, 2>
//          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6>
//            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6(observation.data()), weight)));
//    }
//  }

//  const double * m_pos_2dpoint; // The 2D observation
//};


struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C4_1
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C4_1(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
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

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--
//#pragma omp critical
//      {
    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];


    T pos_proj[3];
    T X[3];
    X[0] = pos_3dpoint[0];
    X[1] = pos_3dpoint[1];
    X[2] = T(0.0);
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, X, pos_proj);

    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    //prepare l
    T R[9];
//    ceres::MatrixAdapter<T, 3, 3>(R);
    ceres::AngleAxisToRotationMatrix(cam_R, R);

    const T& f = cam_intrinsics[OFFSET_FOCAL_LENGTH]; //focal
    const T& u = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
    const T& v = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y
    const T& k1 = cam_intrinsics[OFFSET_DISTO_K1];
    const T& k2 = cam_intrinsics[OFFSET_DISTO_K2];
    const T& k3 = cam_intrinsics[OFFSET_DISTO_K3];
    const T& t1 = cam_intrinsics[OFFSET_DISTO_T1];
    const T& t2 = cam_intrinsics[OFFSET_DISTO_T2];

    const T x_u = u + f * (pos_proj[0] / pos_proj[2]);
    const T y_u = v + f * (pos_proj[1] / pos_proj[2]);

    // Apply distortion (xd,yd) = disto(x_u,y_u)
    const T r2 = x_u*x_u + y_u*y_u;
    const T r4 = r2 * r2;
    const T r6 = r4 * r2;
    const T r_coeff = (1.0 + k1*r2 + k2*r4 + k3*r6);
    const T t_x = t2 * (r2 + 2.0 * x_u*x_u) + 2.0 * t1 * x_u * y_u;
    const T t_y = t1 * (r2 + 2.0 * y_u*y_u) + 2.0 * t2 * x_u * y_u;
    const T x_d = x_u * r_coeff + t_x;
    const T y_d = y_u * r_coeff + t_y;

    // Apply focal length and principal point to get the final image coordinates
    const T x_i = u + f * x_d;
    const T y_i = v + f * y_d;


    const T P1 = f*R[0] + u*R[2];
    const T P2 = f*R[3] + u*R[5];
    const T P3 = f*R[6] + u*R[8];
    const T P4 = f*cam_t[0] + u*cam_t[2];
    const T P5 = f*R[1] + v*R[2];
    const T P6 = f*R[4] + v*R[5];
    const T P7 = f*R[7] + v*R[8];
    const T P8 = f*cam_t[1] + v*cam_t[2];
    const T P9 = R[2];
    const T P10 = R[5];
    const T P11 = R[8];
    const T P12 = cam_t[2];

//    const T x_i = (P1*pos_3dpoint[0] + P2*pos_3dpoint[1] + P4) / (P9*pos_3dpoint[0] + P10*pos_3dpoint[1] + P12);
//    const T y_i = (P5*pos_3dpoint[0] + P6*pos_3dpoint[1] + P8) / (P9*pos_3dpoint[0] + P10*pos_3dpoint[1] + P12);



    const T xu = T(m_pos_2dpoint[0]);
    const T yu = T(m_pos_2dpoint[1]);

    const T l0 = P7 - P11*yu;
    const T l1 = P11*xu - P3;
    const T l2 = P3*yu - P7*xu;

//    if(l0 != T(0.0))
//    {
//        out_residuals[0] = (- y_i*l1 - l2)/l0 - x_i;

//    }else{
//        out_residuals[0] = T(0.0);
//    }
//    if(l1 != T(0.0))
//    {
//        out_residuals[1] = ( -x_i*l0 - l2)/l1 - y_i;
//    }else{
//        out_residuals[1] = T(0.0);

//    }
//    out_residuals[0] = ( x_i*l0 + y_i*l1 + l2);
//    out_residuals[0] = abs(x_i*l0 + y_i*l1 + l2)/sqrt(l0*l0 + l1*l1);
    out_residuals[0] = (x_i*l0 + y_i*l1 + l2)/sqrt(l0*l0 + l1*l1);
//    out_residuals[1] = T(0.0);
//    out_residuals[0] = abs(x_i*l0 + y_i*l1 + l2)/sqrt(l0*l0 + l1*l1);

//    if(out_residuals[0] > T(0.1))
//    {
//        std::cout << "cam_R : " << cam_R[0] << " " << cam_R[1] << " " << cam_R[2] << std::endl;
//        std::cout << "R : " << R[0] << " " << R[1] << " " << R[2]
//                            << R[3] << " " << R[4] << " " << R[5]
//                            << R[6] << " " << R[7] << " " << R[8] << std::endl;
//        std::cout << "cam_t : " << cam_t[0] << " " << cam_t[1] << " " << cam_t[2] << std::endl;
//        std::cout << "f : " << f << " u : " << u << " v : " << v << std::endl;
//        std::cout << "xu : " << xu << " yu : " << yu << std::endl;
//        std::cout << "pos_3dpoint[0]: " << pos_3dpoint[0] << " pos_3dpoint[1]: " << pos_3dpoint[1] << " pos_3dpoint[2]: " << pos_3dpoint[2] << std::endl;
//        std::cout << "out_residuals[0]: "  << out_residuals[0] << std::endl;
//        std::cout << "****-----------------****" << std::endl;
//    }

//      }

//    T pos_proj[3];
//    // Rotate the point according the camera rotation
//    ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_proj);

//    // Apply the camera translation
//    pos_proj[0] += cam_t[0];
//    pos_proj[1] += cam_t[1];
//    pos_proj[2] += cam_t[2];

//    // Transform the point from homogeneous to euclidean (undistorted point)
//    const T x_u = pos_proj[0] / pos_proj[2];
//    const T y_u = pos_proj[1] / pos_proj[2];

//    //--
//    // Apply intrinsic parameters
//    //--


//    const T& k1 = cam_intrinsics[OFFSET_DISTO_K1];
//    const T& k2 = cam_intrinsics[OFFSET_DISTO_K2];
//    const T& k3 = cam_intrinsics[OFFSET_DISTO_K3];
//    const T& t1 = cam_intrinsics[OFFSET_DISTO_T1];
//    const T& t2 = cam_intrinsics[OFFSET_DISTO_T2];

//    // Apply distortion (xd,yd) = disto(x_u,y_u)
//    const T r2 = x_u*x_u + y_u*y_u;
//    const T r4 = r2 * r2;
//    const T r6 = r4 * r2;
//    const T r_coeff = (1.0 + k1*r2 + k2*r4 + k3*r6);
//    const T t_x = t2 * (r2 + 2.0 * x_u*x_u) + 2.0 * t1 * x_u * y_u;
//    const T t_y = t1 * (r2 + 2.0 * y_u*y_u) + 2.0 * t2 * x_u * y_u;
//    const T x_d = x_u * r_coeff + t_x;
//    const T y_d = y_u * r_coeff + t_y;

//    // Apply focal length and principal point to get the final image coordinates
//    const T projected_x = principal_point_x + focal * x_d;
//    const T projected_y = principal_point_y + focal * y_d;

//    // Compute and return the error is the difference between the predicted
//    //  and observed position
//    out_residuals[0] = projected_x - m_pos_2dpoint[0];
//    out_residuals[1] = projected_y - m_pos_2dpoint[1];

    return true;
  }

  static int num_residuals() { return 1; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C4_1, 1, 8, 6, 2>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C4_1(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C4_1>, 1, 8, 6, 2>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C4_1>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C4_1(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};



/**
 * @brief Ceres functor with constrained 3D points to use a Pinhole_Intrinsic_Fisheye
 *
 *  Data parameter blocks are the following <2,8,6,3>
 *  - 2 => dimension of the residuals,
 *  - 7 => the intrinsic data block [focal, principal point x, principal point y, K1, K2, K3, K4],
 *  - 6 => the camera extrinsic data block (camera orientation and position) [R;t],
 *         - rotation(angle axis), and translation [rX,rY,rZ,tx,ty,tz].
 *  - 3 => a 3D point data block.
 *
 */

struct ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum : uint8_t {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2,
    OFFSET_DISTO_K1 = 3,
    OFFSET_DISTO_K2 = 4,
    OFFSET_DISTO_K3 = 5,
    OFFSET_DISTO_K4 = 6,
  };

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, k4 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];

    T pos_proj[3];
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_proj);

    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    // Transform the point from homogeneous to euclidean (undistorted point)
    const T x_u = pos_proj[0] / pos_proj[2];
    const T y_u = pos_proj[1] / pos_proj[2];

    //--
    // Apply intrinsic parameters
    //--
    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];
    const T& k1 = cam_intrinsics[OFFSET_DISTO_K1];
    const T& k2 = cam_intrinsics[OFFSET_DISTO_K2];
    const T& k3 = cam_intrinsics[OFFSET_DISTO_K3];
    const T& k4 = cam_intrinsics[OFFSET_DISTO_K4];

    // Apply distortion (xd,yd) = disto(x_u,y_u)
    const T r2 = x_u*x_u + y_u*y_u;
    const T r = sqrt(r2);
    const T
      theta = atan(r),
      theta2 = theta*theta,
      theta3 = theta2*theta,
      theta4 = theta2*theta2,
      theta5 = theta4*theta,
      theta7 = theta3*theta3*theta, //thetha6*theta
      theta8 = theta4*theta4,
      theta9 = theta8*theta;
    const T theta_dist = theta + k1*theta3 + k2*theta5 + k3*theta7 + k4*theta9;
    const T inv_r = r > T(1e-8) ? T(1.0)/r : T(1.0);
    const T cdist = r > T(1e-8) ? theta_dist * inv_r : T(1.0);

    const T x_d = x_u * cdist;
    const T y_d = y_u * cdist;

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal * x_d;
    const T projected_y = principal_point_y + focal * y_d;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - m_pos_2dpoint[0];
    out_residuals[1] = projected_y - m_pos_2dpoint[1];

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye, 2, 7, 6, 3>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye>, 2, 7, 6, 3>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};

struct ResidualErrorFunctor_Intrinsic_Spherical
{
  ResidualErrorFunctor_Intrinsic_Spherical
  (
    const double* const pos_2dpoint,
    const size_t imageSize_w,
    const size_t imageSize_h
  )
  : m_pos_2dpoint(pos_2dpoint),
    m_imageSize{imageSize_w, imageSize_h}
  {
  }

  /**
  * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
  *   - 3 for rotation(angle axis), 3 for translation
  * @param[in] pos_3dpoint
  * @param[out] out_residuals
  */
  template <typename T>
  bool operator()
  (
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals
  )
  const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];

    T pos_proj[3];
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_proj);

    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    // Transform the coord in is Image space
    const T lon = ceres::atan2(pos_proj[0] , pos_proj[2]); // Horizontal normalization of the  X-Z component
    const T lat = ceres::atan2(-pos_proj[1], ceres::sqrt(pos_proj[0] * pos_proj[0]  + pos_proj[2] * pos_proj[2])); // Tilt angle
    const T coord[] = {lon / (2 * M_PI), lat / (2 * M_PI)}; // normalization

    const T size ( std::max(m_imageSize[0], m_imageSize[1]) );
    const T projected_x = coord[0] * size - 0.5 + m_imageSize[0] / 2.0;
    const T projected_y = coord[1] * size - 0.5 + m_imageSize[1] / 2.0;

    out_residuals[0] = projected_x - m_pos_2dpoint[0];
    out_residuals[1] = projected_y - m_pos_2dpoint[1];

    return true;
  }

  static const int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const cameras::IntrinsicBase * cameraInterface,
    const Vec2 & observation,
    const double weight = 0.0
  )
  {

    if (weight == 0.0)
    {
      return
          (new ceres::AutoDiffCostFunction
              <ResidualErrorFunctor_Intrinsic_Spherical, 2, 6, 3>(
              new ResidualErrorFunctor_Intrinsic_Spherical(
                observation.data(), cameraInterface->w(), cameraInterface->h())));
    }
    else
    {
      return
          (new ceres::AutoDiffCostFunction
              <WeightedCostFunction<ResidualErrorFunctor_Intrinsic_Spherical>, 2, 6, 3>
              (new WeightedCostFunction<ResidualErrorFunctor_Intrinsic_Spherical>
                   (new ResidualErrorFunctor_Intrinsic_Spherical(
                     observation.data(), cameraInterface->w(), cameraInterface->h()),
                     weight)));
    }
  }

  const double * m_pos_2dpoint;  // The 2D observation
  size_t         m_imageSize[2]; // The image width and height
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_DATA_BA_CERES_CAMERA_FUNCTOR_HPP
