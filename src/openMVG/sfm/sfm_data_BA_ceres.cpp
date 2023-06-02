// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data_BA_ceres.hpp"

#ifdef OPENMVG_USE_OPENMP
#include </usr/lib/gcc/x86_64-linux-gnu/9/include/omp.h>
#endif

#include "ceres/problem.h"
#include "ceres/solver.h"
#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/geometry/Similarity3.hpp"
#include "openMVG/geometry/Similarity3_Kernel.hpp"
//- Robust estimation - LMeds (since no threshold can be defined)
#include "openMVG/robust_estimation/robust_estimator_LMeds.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres_camera_functor.hpp"
#include "openMVG/sfm/sfm_data_transform.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/types.hpp"

#include <ceres/rotation.h>
#include <ceres/types.h>

#include <iostream>
#include <limits>
#include<fstream>

namespace openMVG {
namespace sfm {

using namespace openMVG::cameras;
using namespace openMVG::geometry;

// Ceres CostFunctor used for SfM pose center to GPS pose center minimization
struct PoseCenterConstraintCostFunction
{
  Vec3 weight_;
  Vec3 pose_center_constraint_;

  PoseCenterConstraintCostFunction
  (
    const Vec3 & center,
    const Vec3 & weight
  ): weight_(weight), pose_center_constraint_(center)
  {
  }

  template <typename T> bool
  operator()
  (
    const T* const cam_extrinsics, // R_t
    T* residuals
  )
  const
  {
    const T * cam_R = &cam_extrinsics[0];
    const T * cam_t = &cam_extrinsics[3];
    const T cam_R_transpose[3] = {-cam_R[0], -cam_R[1], -cam_R[2]};

    T pose_center[3];
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R_transpose, cam_t, pose_center);
    pose_center[0] *= T(-1);
    pose_center[1] *= T(-1);
    pose_center[2] *= T(-1);

    residuals[0] = T(weight_[0]) * (pose_center[0] - T(pose_center_constraint_[0]));
    residuals[1] = T(weight_[1]) * (pose_center[1] - T(pose_center_constraint_[1]));
    residuals[2] = T(weight_[2]) * (pose_center[2] - T(pose_center_constraint_[2]));

    return true;
  }
};

/// Create the appropriate cost functor according the provided input camera intrinsic model.
/// The residual can be weighetd if desired (default 0.0 means no weight).
ceres::CostFunction * IntrinsicsToCostFunction
(
  IntrinsicBase * intrinsic,
  const Vec2 & observation,
  const double weight
)
{
  switch (intrinsic->getType())
  {
    case PINHOLE_CAMERA:
      return ResidualErrorFunctor_Pinhole_Intrinsic::Create(observation, weight);
    break;
    case PINHOLE_CAMERA_RADIAL1:
      return ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1::Create(observation, weight);
    break;
    case PINHOLE_CAMERA_RADIAL3:
      return ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3::Create(observation, weight);
    break;
    case PINHOLE_CAMERA_BROWN:
      return ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2::Create(observation, weight);
    break;
    case PINHOLE_CAMERA_FISHEYE:
      return ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye::Create(observation, weight);
    case CAMERA_SPHERICAL:
      return ResidualErrorFunctor_Intrinsic_Spherical::Create(intrinsic, observation, weight);
    default:
      return nullptr;
  }
}

Bundle_Adjustment_Ceres::BA_Ceres_options::BA_Ceres_options
(
  const bool bVerbose,
  bool bmultithreaded
)
: bVerbose_(bVerbose),
  nb_threads_(1),
  parameter_tolerance_(1e-8), //~= numeric_limits<float>::epsilon()
  bUse_loss_function_(true)
{
  #ifdef OPENMVG_USE_OPENMP
    nb_threads_ = omp_get_max_threads();
  #endif // OPENMVG_USE_OPENMP
  if (!bmultithreaded)
    nb_threads_ = 1;

  bCeres_summary_ = false;

  // Default configuration use a DENSE representation
  linear_solver_type_ = ceres::DENSE_SCHUR;
  preconditioner_type_ = ceres::JACOBI;
  // If Sparse linear solver are available
  // Descending priority order by efficiency (SUITE_SPARSE > CX_SPARSE > EIGEN_SPARSE)
  if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE))
  {
    sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
    linear_solver_type_ = ceres::SPARSE_SCHUR;
  }
  else
  {
    if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE))
    {
      sparse_linear_algebra_library_type_ = ceres::CX_SPARSE;
      linear_solver_type_ = ceres::SPARSE_SCHUR;
    }
    else
    if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::EIGEN_SPARSE))
    {
      sparse_linear_algebra_library_type_ = ceres::EIGEN_SPARSE;
      linear_solver_type_ = ceres::SPARSE_SCHUR;
    }
  }
}


Bundle_Adjustment_Ceres::Bundle_Adjustment_Ceres
(
  Bundle_Adjustment_Ceres::BA_Ceres_options options
)
: ceres_options_(options)
{}

Bundle_Adjustment_Ceres::BA_Ceres_options &
Bundle_Adjustment_Ceres::ceres_options()
{
  return ceres_options_;
}

bool Bundle_Adjustment_Ceres::Adjust
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;
    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
        else
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
      }
    }
    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          IntrinsicsToCostFunction(
            sfm_data.intrinsics[view->id_intrinsic].get(),
            obs_it.second.x,
            options.control_point_opt.weight);

        if (cost_function)
        {
            if (!map_intrinsics.at(view->id_intrinsic).empty())
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_intrinsics.at(view->id_intrinsic)[0],
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
            }
            else
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
  }
        }

      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;

  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " BriefReport : " << summary.BriefReport() << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}


//根据已经约束过的XY,获取当前位置期望的Z坐标
double getExceptionZ(
        const openMVG::geometry::Pose3 &pose,//相机位姿
        const std::vector<double> &intrinInfo,//相机内参
        const Eigen::Vector2d &imgUvData,//图片上的成像坐标
        const Eigen::Vector3d &xyData //已经约束过的XY数据,第3个是Z,不动它
)
{
    //计算摄像机矩阵下的成像点坐标
    Eigen::Vector3d xcPoint;
    xcPoint<<imgUvData(0)-intrinInfo[1],
            imgUvData(1)-intrinInfo[2],intrinInfo[0];
    //取出位姿中的旋转矩阵
    Eigen::Matrix3d rotationMat=pose.rotation();
    //根据成像点坐标和旋转矩阵,计算反射投影射线的方向向量
    Eigen::Vector3d dirVec=rotationMat.transpose()*xcPoint;
    //海平面的法向量
    Eigen::Vector3d planeNorm;
    planeNorm<<0,0,1;
    //反投影方向向量延伸出来的平面的法向量
    Eigen::Vector3d crossPlaneNorm=
            dirVec.cross(dirVec.cross(planeNorm));
    //取出相机位姿中光芯的位置
    Eigen::Vector3d camPosition=pose.center();
    //根据光芯也在这个延伸的平面上,得到平面方程中的D
    double equationD=-camPosition.dot(crossPlaneNorm);
    //根据已经约束的XY,求出当前位置期望的z坐标
    return -(xyData(0)*crossPlaneNorm(0)
             +xyData(1)*crossPlaneNorm(1)+equationD)
            /crossPlaneNorm(2);
}



//赵志豪 20201230
//对三维重建结果做Z坐标的约束,需要先约束XY才能调用这个函数,否则会出问题
bool Bundle_Adjustment_Ceres::adjustWaterZ(
        SfM_Data & sfm_data,
        const Optimize_Options & options
)
{
    //----------
    // Add camera parameters
    // - intrinsics
    // - poses [R|t]

    // Create residuals for each observation in the bundle adjustment problem. The
    // parameters for cameras and points are added automatically.
    //----------


    double pose_center_robust_fitting_error = 0.0;
    openMVG::geometry::Similarity3 sim_to_center;
    bool b_usable_prior = false;
    if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
    {
      // - Compute a robust X-Y affine transformation & apply it
      // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
      {
        // Collect corresponding camera centers
        std::vector<Vec3> X_SfM, X_GPS;
        for (const auto & view_it : sfm_data.GetViews())
        {
          const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
          if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
          {
            X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
            X_GPS.push_back( prior->pose_center_ );
          }
        }
        openMVG::geometry::Similarity3 sim;

        // Compute the registration:
        if (X_GPS.size() > 3)
        {
          const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
          const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
          geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
          const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
          if (lmeds_median != std::numeric_limits<double>::max())
          {
            b_usable_prior = true; // PRIOR can be used safely

            // Compute the median residual error once the registration is applied
            for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
            {
              pos = sim(pos);
            }
            Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
            std::sort(residual.data(), residual.data() + residual.size());
            pose_center_robust_fitting_error = residual(residual.size()/2);

            // Apply the found transformation to the SfM Data Scene
            openMVG::sfm::ApplySimilarity(sim, sfm_data);

            // Move entire scene to center for better numerical stability
            Vec3 pose_centroid = Vec3::Zero();
            for (const auto & pose_it : sfm_data.poses)
            {
              pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
            }
            sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
            openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
          }
        }
      }
    }

    //true的时候表示仅仅对Z做优化,否则约束Z的时候会对XY做优化
    bool singleZ=false;
    const bool adjustN=true;//是否对折射率做优化
    ceres::Problem problem;
    //初始化折射率
    double refIndex=1.33;
    //把折射率添加到待优化的问题中
    problem.AddParameterBlock(&refIndex,1);
    //判断是否对折射率做优化
    if(!adjustN)//不做优化的情况设置为常量
        problem.SetParameterBlockConstant(&refIndex);

    // Data wrapper for refinement:
    Hash_Map<IndexT, std::vector<double> > map_poses, init_map_poses;

    // Setup Poses data & subparametrization
    for (const auto & pose_it : sfm_data.poses)
    {
      const IndexT indexPose = pose_it.first;
      const Pose3 & pose = pose_it.second;
      const Mat3 R = pose.rotation();
      const Vec3 C = pose.center();

      double angleAxis[3];
      ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
      // angleAxis + translation
      //这个地方添加外参本来是仅仅添加了c0和c1的,但后面损失函数还是会用到第3个维度的,所以又添加上了
      map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(0), C(1), C(2)};
      //map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(0), C(1)};
      //double  parameter_block[6] = {angleAxis[0], angleAxis[1], angleAxis[2], C(0), C(1), C(2)};
      double *blockPtr;
      //仅仅优化Z的方案
      if(singleZ)
      {
          blockPtr=&(map_poses[indexPose][5]);
          problem.AddParameterBlock(blockPtr, 1);
          problem.AddParameterBlock(&(map_poses[indexPose][0]), 3);
          //problem.SetParameterBlockConstant(&(map_poses[indexPose][0]));
      }
      else //和重投影误差相结合的方案
      {
          blockPtr=&(map_poses[indexPose][0]);
          problem.AddParameterBlock(blockPtr, 6);
      }
      //判断外参是否设为常量
      if(options.extrinsics_opt==Extrinsic_Parameter_Type::NONE)
          problem.SetParameterBlockConstant(blockPtr);

  //    std::vector<int> vec_constant_extrinsic;
  //    vec_constant_extrinsic.push_back(5);
  //    if (!vec_constant_extrinsic.empty())
  //    {
  //        ceres::SubsetParameterization *subset_parameterization =
  //          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
  //        problem.SetParameterization(parameter_block, subset_parameterization);
  //    }

    }
    init_map_poses = map_poses;


    Hash_Map<IndexT, std::vector<double> > map_intrinsics;
    // Setup Intrinsics data & subparametrization
    for (const auto & intrinsic_it : sfm_data.intrinsics)
    {
      const IndexT indexCam = intrinsic_it.first;
      if (isValid(intrinsic_it.second->getType()))
      {
          //获取内参数据
          std::vector<double> intrinsicVecData=intrinsic_it.second->getParams();
          //把数据保存到哈希表中
          map_intrinsics[indexCam]={intrinsicVecData[0],intrinsicVecData[1],intrinsicVecData[2]};
        //map_intrinsics[indexCam] = intrinsic_it.second->getParams();
        if (!map_intrinsics[indexCam].empty())
        {
            double *blockPtr=&(map_intrinsics[indexCam][0]);
          problem.AddParameterBlock(blockPtr, 3);
          if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
          {
            // set the whole parameter block as constant for best performance
            problem.SetParameterBlockConstant(blockPtr);
          }
          else
          {
            const std::vector<int> vec_constant_intrinsic =
              intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
            if (!vec_constant_intrinsic.empty())
            {
              ceres::SubsetParameterization *subset_parameterization =
                new ceres::SubsetParameterization(
                  map_intrinsics[indexCam].size(), vec_constant_intrinsic);
              problem.SetParameterization(blockPtr, subset_parameterization);
            }
          }
        }
      }
      else
      {
        std::cerr << "Unsupported camera type." << std::endl;
      }
    }

    // Set a LossFunction to be less penalized by false measurements
    //  - set it to nullptr if you don't want use a lossFunction.
    ceres::LossFunction * p_LossFunction =
      ceres_options_.bUse_loss_function_ ?
        new ceres::HuberLoss(Square(4.0))
        : nullptr;
    //约束光心坐标的时候使用的损失函数
    ceres::LossFunction* priorFunction=
            nullptr;

    // For all visibility add reprojections errors:
    Landmarks save_structure;//sfm_data.structure
    for (auto & structure_landmark_it : sfm_data.structure)
    {
      const Observations & obs = structure_landmark_it.second.obs;

  //    if (structure_landmark_it.second.X(2) > 1.039592 || structure_landmark_it.second.X(1) < -0.16)
  //        continue;
  //    if (structure_landmark_it.second.X(1) < 0.185018 && structure_landmark_it.second.X(2) > 2.5108049)
  //        continue;
  //    if (structure_landmark_it.second.X(2) > 1.062)
  //        continue;

      if (obs.size() < 3)
          continue;

      save_structure.insert(structure_landmark_it);
    }
    sfm_data.structure = save_structure;
    for (auto & structure_landmark_it : sfm_data.structure)
    {
      const Observations & obs = structure_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        ceres::CostFunction* cost_function;
        //仅仅约束Z的方案
        //约束Z的算法
        if(singleZ)
        {
            cost_function =
              WaterPlaneConstrainCostFunctionZ::Create(obs_it.second.x,4855,
                        sfm_data.poses[view->id_pose].rotation(),
                    sfm_data.poses[view->id_pose].center(),
                    structure_landmark_it.second.X,
                    map_intrinsics[view->id_intrinsic]);
            problem.AddResidualBlock(cost_function,
              p_LossFunction,
              &map_poses[view->id_pose][5],
              &structure_landmark_it.second.X.data()[2],
            &map_poses[view->id_pose][0]);
        }
        else //约束Z的同时优化XY
        {
            //做水盆重建的时候临时把高度写成4855,重建海面的时候使用外部输入的水面高度
            cost_function =
                    ConstrainZSimilarXY::Create(obs_it.second.x,4855);
            problem.AddResidualBlock(cost_function,
              p_LossFunction,
                &map_intrinsics[view->id_intrinsic][0],
              &map_poses[view->id_pose][0],
              structure_landmark_it.second.X.data(),&refIndex);
        }
      }
      if (options.structure_opt == Structure_Parameter_Type::NONE)
        problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
    }
    //如果对折射率做优化，对折射率做一个基本的约束，让它在1.33附近浮动
//    if(adjustN)
//    {
//        //新建误差函数
//        ceres::CostFunction *cost_function=new ceres::AutoDiffCostFunction<ConstrinRefN, 1, 1>(
//            new ConstrinRefN(1.33));
//        //把误差函数添加到problem里面
//        problem.AddResidualBlock(cost_function,priorFunction,
//                                 &refIndex);
//    }

    // Add Pose prior constraints if any
    if (b_usable_prior)
    //if(false)
    {
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
            ceres::CostFunction * cost_function;
            //仅仅约束Z的算法
            if(singleZ)
            {
                cost_function =
                  new ceres::AutoDiffCostFunction<ConstrainZPrior, 1, 1>(
                    new ConstrainZPrior(prior->pose_center_, prior->center_weight_));
                problem.AddResidualBlock(cost_function, priorFunction, &map_poses[prior->id_view][5]);
            }
            else //和重投影误差相结合的方法
            {
                cost_function =
                  new ceres::AutoDiffCostFunction<ConstrainXYPrior, 3, 6>(
                    new ConstrainXYPrior(prior->pose_center_, prior->center_weight_));
                problem.AddResidualBlock(cost_function, priorFunction, &map_poses[prior->id_view][0]);
            }



        }
      }
    }

    // Configure a BA engine and run it
    //  Make Ceres automatically detect the bundle structure.
  //  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
  //  linear_solver_type_ = ceres::SPARSE_SCHUR;
  //  linear_solver_type_ = ceres::DENSE_SCHUR;
  //  preconditioner_type_ = ceres::JACOBI;
    ceres::Solver::Options ceres_config_options;
    ceres_config_options.max_num_iterations = 500;
    ceres_config_options.preconditioner_type =
      static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
    ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
      static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
    ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
      static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
    ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
    ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
    ceres_config_options.num_threads = ceres_options_.nb_threads_;
    ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
    ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

  //  std::cout << "start ceres solver" << std::endl;
    // Solve BA
    ceres::Solver::Summary summary;
    ceres::Solve(ceres_config_options, &problem, &summary);
    if (ceres_options_.bCeres_summary_)
      std::cout << summary.FullReport() << std::endl;

    // If no error, get back refined parameters
    if (!summary.IsSolutionUsable())
    {
      if (ceres_options_.bVerbose_)
        std::cout << "Bundle Adjustment failed." << std::endl;
      return false;
    }
    else // Solution is usable
    {
      if (ceres_options_.bVerbose_)
      {
        // Display statistics about the minimization
        std::cout << std::endl
          << "Bundle Adjustment statistics (approximated RMSE):\n"
          << " #views: " << sfm_data.views.size() << "\n"
          << " #poses: " << sfm_data.poses.size() << "\n"
          << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
          << " #tracks: " << sfm_data.structure.size() << "\n"
          << " #residuals: " << summary.num_residuals << "\n"
          << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
          << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
          << " BriefReport : " << summary.BriefReport() << "\n"
          << " Time (s): " << summary.total_time_in_seconds << "\n"
          << std::endl;
        if (options.use_motion_priors_opt)
          std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
      }

      // Update camera poses with refined data
      //if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
      {
        for (auto & pose_it : sfm_data.poses)
        {
          const IndexT indexPose = pose_it.first;

          Mat3 R_refined;
          ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
          Pose3 & pose = pose_it.second;
          Vec3 C_init = pose.center();
          //Vec3 C_refined(map_poses[indexPose][3], map_poses[indexPose][4], C_init(2));
          Vec3 C_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
          // Update the pose
  //        Pose3 & pose = pose_it.second;
          pose = Pose3(R_refined, C_refined);
        }
      }

      // Update camera intrinsics with refined data
      //if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
      {
        for (auto & intrinsic_it : sfm_data.intrinsics)
        {
          const IndexT indexCam = intrinsic_it.first;

          const std::vector<double> & vec_params = map_intrinsics[indexCam];
          intrinsic_it.second->updateFromParams(vec_params);
        }
      }
        //输出优化后的折射率结果
      if(adjustN) std::cout<<"refractive index: "<<refIndex<<std::endl;
      // Structure is already updated directly if needed (no data wrapping)

      if (b_usable_prior)
      {
        // set back to the original scene centroid
        openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

        //--
        // - Compute some fitting statistics
        //--

        // Collect corresponding camera centers
        std::vector<Vec3> X_SfM, X_GPS;
        for (const auto & view_it : sfm_data.GetViews())
        {
          const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
          if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
          {
            X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
            X_GPS.push_back( prior->pose_center_ );
          }
        }
        // Compute the registration fitting error (once BA with Prior have been used):
        if (X_GPS.size() > 3)
        {
          // Compute the median residual error
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::cout
            << "Pose prior statistics (user units):\n"
            << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
            << " - Final fitting error:";
          minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
        }
      }
      return true;
    }
}




bool Bundle_Adjustment_Ceres::Adjust_water_xy
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options,
    bool getPriorDiff //为true的时候,不对点云做优化,仅仅是计算和先验光心位置的差
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }


  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_poses, init_map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;
    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 C = pose.center();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    //这个地方添加外参本来是仅仅添加了c0和c1的,但后面损失函数还是会用到第3个维度的,所以又添加上了
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(0), C(1), C(2)};
    //map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(0), C(1)};
    //double  parameter_block[6] = {angleAxis[0], angleAxis[1], angleAxis[2], C(0), C(1), C(2)};
    double *blockPtr=&(map_poses[indexPose][0]);
    problem.AddParameterBlock(blockPtr, 6);
    //判断外参是否设为常量
    if(options.extrinsics_opt==Extrinsic_Parameter_Type::NONE)
        problem.SetParameterBlockConstant(blockPtr);

//    std::vector<int> vec_constant_extrinsic;
//    vec_constant_extrinsic.push_back(5);
//    if (!vec_constant_extrinsic.empty())
//    {
//        ceres::SubsetParameterization *subset_parameterization =
//          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
//        problem.SetParameterization(parameter_block, subset_parameterization);
//    }

  }
  init_map_poses = map_poses;


  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
        //获取内参数据
        std::vector<double> intrinsicVecData=intrinsic_it.second->getParams();
        //把数据保存到哈希表中
        //map_intrinsics[indexCam]={intrinsicVecData[0],intrinsicVecData[1],intrinsicVecData[2]};
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
          double *blockPtr=&(map_intrinsics[indexCam][0]);
        problem.AddParameterBlock(blockPtr, 8);
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(blockPtr);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(blockPtr, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
          nullptr;

  //先验位姿的损失函数
  ceres::LossFunction* priorFunction=
          //new ceres::HuberLoss(4);
          nullptr;

  //确保每个点云至少有3个图片能看到它，这里的功能在函数外面被单独实现了
  // For all visibility add reprojections errors:
//  Landmarks save_structure;//sfm_data.structure
//  for (auto & structure_landmark_it : sfm_data.structure)
//  {
//    const Observations & obs = structure_landmark_it.second.obs;

//    if (obs.size() < 3)
//        continue;

//    save_structure.insert(structure_landmark_it);
//  }
//  //如果仅仅是看和先验光心的区别,不做优化
//  sfm_data.structure = save_structure;
  const bool useVaribleXc=false;//是否使用可变的xc
  if(!getPriorDiff)
  {
      for (auto & structure_landmark_it : sfm_data.structure)
      {
        Observations & obs = structure_landmark_it.second.obs;

        for (auto & obs_it : obs)
        {
          // Build the residual block corresponding to the track observation:
          const View * view = sfm_data.views.at(obs_it.first).get();

          //如果使用可变xc的算法
          if(useVaribleXc)
          {
              //获取内参
              std::vector<double> &tempIntrinsic=
                      map_intrinsics[view->id_intrinsic];
              //获取成像点坐标
              Vec2 &imgX=obs_it.second.x;
              //获取xc
              Vec3 &tempXc=obs_it.second.xcPoint;
              //初始化xc
              tempXc[0]=(imgX[0]-tempIntrinsic[1])/tempIntrinsic[0];
              tempXc[1]=(imgX[1]-tempIntrinsic[2])/tempIntrinsic[0];
              tempXc[2]=1;
              //建立校正xc用的误差函数
              ceres::CostFunction *cost_function=
                      XcPointProjectError::Create(imgX);
              problem.AddResidualBlock(cost_function,
                   nullptr,
                   &tempIntrinsic[0],
                      tempXc.data());
              //建立针对可变误差的XY坐标约束函数
              ceres::CostFunction *xyCostFunction=
                      XyConstrainWithVaribleXc::Create(obs.size());
              problem.AddResidualBlock(xyCostFunction,
                                       p_LossFunction,
                                       tempXc.data(),
                                       &map_poses[view->id_pose][0],
                      structure_landmark_it.second.X.data());
          }
          else
          {
              //获取当前位置相机坐标系下的坐标
              const Vec2 &imgPt=obs_it.second.x;
              Vec3 xcPoint=sfm_data.intrinsics.at(view->id_intrinsic)->getXcPoint(imgPt[0],imgPt[1]);
              ceres::CostFunction* cost_function =//权值10目前是最好的结果
                WaterPlaneConstrainCostFunctionXy::Create(xcPoint,-295);
              problem.AddResidualBlock(cost_function,
                p_LossFunction,
                &map_intrinsics[view->id_intrinsic][0],
                &map_poses[view->id_pose][0],
                structure_landmark_it.second.X.data());
          }
        }
        if (options.structure_opt == Structure_Parameter_Type::NONE)
          problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
      }
  }






  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
          //当前需要记录的外参信息
          std::vector<double> &tempExt=map_poses[prior->id_view];
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<ConstrainXYPrior, 3, 6>(
            new ConstrainXYPrior(prior->pose_center_,
                                 prior->center_weight_,
                                 false));
                                 //Vec3(tempExt[0],tempExt[1],tempExt[2]),false));
        problem.AddResidualBlock(cost_function, priorFunction, &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;

  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 300;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " BriefReport : " << summary.BriefReport() << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    //if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Pose3 & pose = pose_it.second;
        Vec3 C_init = pose.center();
        //Vec3 C_refined(map_poses[indexPose][3], map_poses[indexPose][4], C_init(2));
        Vec3 T_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
//        Pose3 & pose = pose_it.second;
        //如果仅仅是看看差别,不保留新的光心位置
        if(!getPriorDiff)
            pose = Pose3(R_refined, T_refined);

        std::cout << "diff pose : " << init_map_poses[indexPose][3] - map_poses[indexPose][3] << " " << init_map_poses[indexPose][4] - map_poses[indexPose][4]
                <<" "<<init_map_poses[indexPose][5]-map_poses[indexPose][5]<< std::endl;
      }
    }

    // Update camera intrinsics with refined data
    //if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}

bool Bundle_Adjustment_Ceres::Adjust_water_min_z
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  double & water_plane,
  //double & ref_N,
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses, init_map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 C = pose.center();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(0), C(1), C(2)};
    //map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    problem.SetParameterBlockConstant(parameter_block);

//    std::vector<int> vec_constant_extrinsic;
//    vec_constant_extrinsic.push_back(5);
//    if (!vec_constant_extrinsic.empty())
//    {
//        ceres::SubsetParameterization *subset_parameterization =
//          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
//        problem.SetParameterization(parameter_block, subset_parameterization);
//    }

  }
  init_map_poses = map_poses;

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
          std::cout << "error 1" << std::endl;
        map_intrinsics[indexCam].push_back(water_plane);
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, 4);//map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
            std::vector<int> vec_constant_extrinsic;
            vec_constant_extrinsic.push_back(0);
            vec_constant_extrinsic.push_back(1);
            vec_constant_extrinsic.push_back(2);
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(4, vec_constant_extrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);

          // set the whole parameter block as constant for best performance
          //problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  std::cout << "error 2" << std::endl;
  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  Landmarks save_structure;//sfm_data.structure
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

//    if (structure_landmark_it.second.X(0) < -5.4575 || structure_landmark_it.second.X(0) > 2.059387 || structure_landmark_it.second.X(1) < -2.205668 || structure_landmark_it.second.X(1) > 1.797887 )
//        continue;
    if (obs.size() < 3)
        continue;

    save_structure.insert(structure_landmark_it);
  }
  sfm_data.structure = save_structure;
  std::cout << "error 3" << std::endl;
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        ResidualErrorFunctor_Pinhole_Intrinsic_Water_minz::Create(obs_it.second.x);
      problem.AddResidualBlock(cost_function,
        p_LossFunction,
        &map_intrinsics[view->id_intrinsic][0],
        &map_poses[view->id_pose][0],
        //&ref_N,
        structure_landmark_it.second.X.data());

//      double * parameter_block = structure_landmark_it.second.X.data();
//      std::vector<int> vec_constant_extrinsic;
//      vec_constant_extrinsic.push_back(0);
//      vec_constant_extrinsic.push_back(1);
//      ceres::SubsetParameterization *subset_parameterization =
//        new ceres::SubsetParameterization(3, vec_constant_extrinsic);
//      problem.SetParameterization(parameter_block, subset_parameterization);

    }
//    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }


  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;

  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 1500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_/10000;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " BriefReport : " << summary.BriefReport() << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    //if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Pose3 & pose = pose_it.second;
        Vec3 C_init = pose.center();
        //Vec3 C_refined(C_init(0), C_init(1), map_poses[indexPose][3]);
        Vec3 C_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
//        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, C_refined);

        //std::cout << "diff pose : " << init_map_poses[indexPose][3] - map_poses[indexPose][3] << " " << init_map_poses[indexPose][4] - map_poses[indexPose][4] << std::endl;
      }
    }

    std::cout << "error 2 " << std::endl;
    std::cout << "error 2 " << map_intrinsics[0][3] << std::endl;

    // Update camera intrinsics with refined data
    //if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
          std::cout << "error 2.1 " << std::endl;
        const IndexT indexCam = intrinsic_it.first;
        //
        Hash_Map<IndexT, std::vector<double> > map_intrinsics_res;
        map_intrinsics_res[indexCam].push_back(map_intrinsics[indexCam][0]);
        map_intrinsics_res[indexCam].push_back(map_intrinsics[indexCam][1]);
        map_intrinsics_res[indexCam].push_back(map_intrinsics[indexCam][2]);
//        std::cout << "error 2.2 " << std::endl;
//        map_intrinsics_res[indexCam][1] = map_intrinsics[indexCam][1];
//        std::cout << "error 2.2 " << std::endl;
//        map_intrinsics_res[indexCam][2] = map_intrinsics[indexCam][2];
//        std::cout << "error 2.2 " << std::endl;
        //
        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }
    water_plane = map_intrinsics[0][3];

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}


bool Bundle_Adjustment_Ceres::Adjust_water_z
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  double & water_plane,
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses, init_map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 C = pose.center();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(0), C(1), C(2)};
    //map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);

//    std::vector<int> vec_constant_extrinsic;
//    vec_constant_extrinsic.push_back(5);
//    if (!vec_constant_extrinsic.empty())
//    {
//        ceres::SubsetParameterization *subset_parameterization =
//          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
//        problem.SetParameterization(parameter_block, subset_parameterization);
//    }

  }
  init_map_poses = map_poses;

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  Landmarks save_structure;//sfm_data.structure
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

//        if (structure_landmark_it.second.X(2) > 2.53322799) //72
//            continue;
//        if (structure_landmark_it.second.X(2) > 1.0416) //73
//            continue;

//        if (structure_landmark_it.second.X(2) > 1.063) //71_1
//            continue;

    if (obs.size() < 3)
        continue;

    save_structure.insert(structure_landmark_it);
  }
  sfm_data.structure = save_structure;
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        ResidualErrorFunctor_Pinhole_Intrinsic_WaterZ::Create(obs_it.second.x, water_plane);
      problem.AddResidualBlock(cost_function,
        p_LossFunction,
        &map_intrinsics[view->id_intrinsic][0],
        &map_poses[view->id_pose][0],
        structure_landmark_it.second.X.data());

//      double * parameter_block = structure_landmark_it.second.X.data();
//      std::vector<int> vec_constant_extrinsic;
//      vec_constant_extrinsic.push_back(0);
//      vec_constant_extrinsic.push_back(1);
//      ceres::SubsetParameterization *subset_parameterization =
//        new ceres::SubsetParameterization(3, vec_constant_extrinsic);
//      problem.SetParameterization(parameter_block, subset_parameterization);

    }
//    if (options.structure_opt == Structure_Parameter_Type::NONE)
//      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }


  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;

  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 1500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_/10000;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " BriefReport : " << summary.BriefReport() << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    //if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Pose3 & pose = pose_it.second;
        Vec3 C_init = pose.center();
        //Vec3 C_refined(C_init(0), C_init(1), map_poses[indexPose][3]);
        Vec3 C_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
//        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, C_refined);

        //std::cout << "diff pose : " << init_map_poses[indexPose][3] - map_poses[indexPose][3] << " " << init_map_poses[indexPose][4] - map_poses[indexPose][4] << std::endl;
      }
    }

    // Update camera intrinsics with refined data
    //if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}


bool Bundle_Adjustment_Ceres::Adjust_water_z_fixC1
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  double & water_plane,
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses, init_map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 C = pose.center();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(0), C(1), C(2)};
    //map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);

    if (indexPose == 0 || indexPose == 1)
        problem.SetParameterBlockConstant(parameter_block);
//    std::vector<int> vec_constant_extrinsic;
//    vec_constant_extrinsic.push_back(5);
//    if (!vec_constant_extrinsic.empty())
//    {
//        ceres::SubsetParameterization *subset_parameterization =
//          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
//        problem.SetParameterization(parameter_block, subset_parameterization);
//    }

  }
  init_map_poses = map_poses;

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  Landmarks save_structure;//sfm_data.structure
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

//        if (structure_landmark_it.second.X(2) > 2.53322799) //72
//            continue;
//        if (structure_landmark_it.second.X(2) > 1.0416) //73
//            continue;

//        if (structure_landmark_it.second.X(2) > 1.063) //71_1
//            continue;

    if (obs.size() < 3)
        continue;

    save_structure.insert(structure_landmark_it);
  }
  sfm_data.structure = save_structure;
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        ResidualErrorFunctor_Pinhole_Intrinsic_WaterZ::Create(obs_it.second.x, water_plane);
      problem.AddResidualBlock(cost_function,
        p_LossFunction,
        &(map_intrinsics[view->id_intrinsic][0]),
        &map_poses[view->id_pose][0],
        structure_landmark_it.second.X.data());

//      double * parameter_block = structure_landmark_it.second.X.data();
//      std::vector<int> vec_constant_extrinsic;
//      vec_constant_extrinsic.push_back(0);
//      vec_constant_extrinsic.push_back(1);
//      ceres::SubsetParameterization *subset_parameterization =
//        new ceres::SubsetParameterization(3, vec_constant_extrinsic);
//      problem.SetParameterization(parameter_block, subset_parameterization);

    }
//    if (options.structure_opt == Structure_Parameter_Type::NONE)
//      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }


  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;

  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 1500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_/10000;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " BriefReport : " << summary.BriefReport() << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    //if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Pose3 & pose = pose_it.second;
        Vec3 C_init = pose.center();
        //Vec3 C_refined(C_init(0), C_init(1), map_poses[indexPose][3]);
        Vec3 C_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
//        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, C_refined);

        //std::cout << "diff pose : " << init_map_poses[indexPose][3] - map_poses[indexPose][3] << " " << init_map_poses[indexPose][4] - map_poses[indexPose][4] << std::endl;
      }
    }

    // Update camera intrinsics with refined data
    //if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}

bool Bundle_Adjustment_Ceres::Adjust_water_z_n
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  double & water_plane,
  double & ref_N,
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses, init_map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 C = pose.center();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(0), C(1), C(2)};
    //map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);

//    std::vector<int> vec_constant_extrinsic;
//    vec_constant_extrinsic.push_back(5);
//    if (!vec_constant_extrinsic.empty())
//    {
//        ceres::SubsetParameterization *subset_parameterization =
//          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
//        problem.SetParameterization(parameter_block, subset_parameterization);
//    }

  }
  init_map_poses = map_poses;

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
          std::cout << "error 1" << std::endl;
        map_intrinsics[indexCam].push_back(ref_N);
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, 4);//map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
            std::vector<int> vec_constant_extrinsic;
            vec_constant_extrinsic.push_back(0);
            vec_constant_extrinsic.push_back(1);
            vec_constant_extrinsic.push_back(2);
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(4, vec_constant_extrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);

          // set the whole parameter block as constant for best performance
          //problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  std::cout << "error 2" << std::endl;
  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  Landmarks save_structure;//sfm_data.structure
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

//    if (structure_landmark_it.second.X(0) < -5.4575 || structure_landmark_it.second.X(0) > 2.059387 || structure_landmark_it.second.X(1) < -2.205668 || structure_landmark_it.second.X(1) > 1.797887 )
//        continue;
    if (obs.size() < 3)
        continue;

    save_structure.insert(structure_landmark_it);
  }
  sfm_data.structure = save_structure;
  std::cout << "error 3" << std::endl;
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        ResidualErrorFunctor_Pinhole_Intrinsic_WaterZN::Create(obs_it.second.x, water_plane);
      problem.AddResidualBlock(cost_function,
        p_LossFunction,
        &map_intrinsics[view->id_intrinsic][0],
        &map_poses[view->id_pose][0],
        //&ref_N,
        structure_landmark_it.second.X.data());

//      double * parameter_block = structure_landmark_it.second.X.data();
//      std::vector<int> vec_constant_extrinsic;
//      vec_constant_extrinsic.push_back(0);
//      vec_constant_extrinsic.push_back(1);
//      ceres::SubsetParameterization *subset_parameterization =
//        new ceres::SubsetParameterization(3, vec_constant_extrinsic);
//      problem.SetParameterization(parameter_block, subset_parameterization);

    }
//    if (options.structure_opt == Structure_Parameter_Type::NONE)
//      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }


  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;

  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 800;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_/10000;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " BriefReport : " << summary.BriefReport() << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    //if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Pose3 & pose = pose_it.second;
        Vec3 C_init = pose.center();
        //Vec3 C_refined(C_init(0), C_init(1), map_poses[indexPose][3]);
        Vec3 C_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
//        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, C_refined);

        //std::cout << "diff pose : " << init_map_poses[indexPose][3] - map_poses[indexPose][3] << " " << init_map_poses[indexPose][4] - map_poses[indexPose][4] << std::endl;
      }
    }

    std::cout << "error 2 " << std::endl;
    std::cout << "error 2 " << map_intrinsics[0][3] << std::endl;

    // Update camera intrinsics with refined data
    //if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
          std::cout << "error 2.1 " << std::endl;
        const IndexT indexCam = intrinsic_it.first;
        //
        Hash_Map<IndexT, std::vector<double> > map_intrinsics_res;
        map_intrinsics_res[indexCam].push_back(map_intrinsics[indexCam][0]);
        map_intrinsics_res[indexCam].push_back(map_intrinsics[indexCam][1]);
        map_intrinsics_res[indexCam].push_back(map_intrinsics[indexCam][2]);
//        std::cout << "error 2.2 " << std::endl;
//        map_intrinsics_res[indexCam][1] = map_intrinsics[indexCam][1];
//        std::cout << "error 2.2 " << std::endl;
//        map_intrinsics_res[indexCam][2] = map_intrinsics[indexCam][2];
//        std::cout << "error 2.2 " << std::endl;
        //
        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }
    ref_N = map_intrinsics[0][3];

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}


bool Bundle_Adjustment_Ceres::Adjust_water_c1_D
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;
    if(obs.size() <3)
        continue;

    std::vector<const View*> view_list;
    std::vector<Vec3> x_list;
    //
    for (const auto & obs_it : obs)
    {
        const View * view = sfm_data.views.at(obs_it.first).get();
        view_list.push_back(view);
        std::vector<double> cam_intrinsics = sfm_data.intrinsics[view->id_intrinsic]->getParams();
        const double focal = cam_intrinsics[0];
        const double principal_point_x = cam_intrinsics[1];
        const double principal_point_y = cam_intrinsics[2];
        Vec3 thisXc = Vec3( (obs_it.second.x(0)-principal_point_x) / focal , (obs_it.second.x(1)-principal_point_y) / focal, 1.0);

        x_list.push_back(thisXc);
    }
    //
    for(std::size_t idI = 0; idI < view_list.size()-2; ++idI)
    {
        for(std::size_t idJ = idI+1; idJ < view_list.size()-1; ++idJ)
        {
            for(std::size_t idK = idJ+1; idK < view_list.size(); ++idK)
            {
                const View * view1 = view_list[idI];
                const View * view2 = view_list[idJ];
                const View * view3 = view_list[idK];
                Vec9 x;/* = Vec6(x_list[idI](0), x_list[idI](1), x_list[idI](2),
                              x_list[idJ](0), x_list[idJ](1), x_list[idJ](2),
                              x_list[idK](0), x_list[idK](1), x_list[idK](2));*/
                x << x_list[idI](0), x_list[idI](1), x_list[idI](2),
                     x_list[idJ](0), x_list[idJ](1), x_list[idJ](2),
                     x_list[idK](0), x_list[idK](1), x_list[idK](2);
                ceres::CostFunction* cost_function =
                  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_D::Create(x);
                problem.AddResidualBlock(cost_function,
                  p_LossFunction,
                  &map_poses[view1->id_pose][0],
                  //&map_intrinsics[view2->id_intrinsic][0],
                  &map_poses[view2->id_pose][0],
                  //&map_intrinsics[view3->id_intrinsic][0],
                  &map_poses[view3->id_pose][0]);

            }

        }

    }

  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          IntrinsicsToCostFunction(
            sfm_data.intrinsics[view->id_intrinsic].get(),
            obs_it.second.x,
            options.control_point_opt.weight);

        if (cost_function)
        {
            if (!map_intrinsics.at(view->id_intrinsic).empty())
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_intrinsics.at(view->id_intrinsic)[0],
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
            }
            else
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
  }
        }

      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;

  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " BriefReport : " << summary.BriefReport() << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}


bool Bundle_Adjustment_Ceres::Adjust_water_c1_S
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & pose_it : sfm_data.GetPoses())
      {
        {
          X_SfM.push_back( pose_it.second.center() );
          X_GPS.push_back( pose_it.second.center());
        }
      }
//      for (const auto & view_it : sfm_data.GetViews())
//      {
//        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
//        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
//        {
//          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
//          X_GPS.push_back( prior->pose_center_ );
//        }
//      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, 3);//map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                3, vec_constant_intrinsic);
//                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_S::Create(obs_it.second.x);
//        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
        else
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
      }
    }
    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_S::Create(obs_it.second.x);
//          IntrinsicsToCostFunction(
//            sfm_data.intrinsics[view->id_intrinsic].get(),
//            obs_it.second.x,
//            options.control_point_opt.weight);

        if (cost_function)
        {
            if (!map_intrinsics.at(view->id_intrinsic).empty())
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_intrinsics.at(view->id_intrinsic)[0],
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
            }
            else
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
  }
        }

      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
////    for (const auto & pose_it : sfm_data.GetPoses())
//    for (const auto & view_it : sfm_data.GetViews())
//    {
//      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
//      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
//      {
////          std::cout << "prior" << std::endl;
////          ceres::CostFunction * cost_function =
////            new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
////              new PoseCenterConstraintCostFunction( pose_it.second.center(), Vec3(1000.0,1000.0,1000.0)));//prior->pose_center_, prior->center_weight_));

////        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
//        ceres::CostFunction * cost_function =
//          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
//            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

//        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);

////        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[pose_it.first][0]);
////        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(4.0)), &map_poses[pose_it.first][0]);
////        problem.AddResidualBlock(cost_function, nullptr, &map_poses[pose_it.first][0]);

//      }
//    }
        for (const auto & pose_it : sfm_data.GetPoses())
        {

          {
//              std::cout << "prior" << std::endl;
              ceres::CostFunction * cost_function =
                new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
                  new PoseCenterConstraintCostFunction( pose_it.second.center(), Vec3(1.0,1.0,1.0)));//prior->pose_center_, prior->center_weight_));


//            problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[pose_it.first][0]);
            problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(4.0)), &map_poses[pose_it.first][0]);
//            problem.AddResidualBlock(cost_function, nullptr, &map_poses[pose_it.first][0]);

          }
        }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;

  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_ / 10000000;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " BriefReport : " << summary.BriefReport() << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & pose_it : sfm_data.GetPoses())
      {
        {
          X_SfM.push_back( pose_it.second.center() );
          X_GPS.push_back( pose_it.second.center());
        }
      }
//      for (const auto & view_it : sfm_data.GetViews())
//      {
//        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
//        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
//        {
//          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
//          X_GPS.push_back( prior->pose_center_ );
//        }
//      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}

bool Bundle_Adjustment_Ceres::Adjust_water_c1_pi
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & pose_it : sfm_data.GetPoses())
      {
        {
          X_SfM.push_back( pose_it.second.center() );
          X_GPS.push_back( pose_it.second.center());
        }
      }
//      for (const auto & view_it : sfm_data.GetViews())
//      {
//        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
//        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
//        {
//          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
//          X_GPS.push_back( prior->pose_center_ );
//        }
//      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, 3);//map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                3, vec_constant_intrinsic);
//                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;
    if(obs.size() <3)
        continue;

    std::vector<const View*> view_list;
    std::vector<Vec2> x_list;
    //
    for (const auto & obs_it : obs)
    {
        const View * view = sfm_data.views.at(obs_it.first).get();
        view_list.push_back(view);
        x_list.push_back(obs_it.second.x);
    }
    //
    for(std::size_t idI = 0; idI < view_list.size()-2; ++idI)
    {
        for(std::size_t idJ = idI+1; idJ < view_list.size()-1; ++idJ)
        {
            for(std::size_t idK = idJ+1; idK < view_list.size(); ++idK)
            {
                const View * view1 = view_list[idI];
                const View * view2 = view_list[idJ];
                const View * view3 = view_list[idK];
                Vec6 x ;//= Vec6(x_list[idI](0), x_list[idI](1), x_list[idJ](0), x_list[idJ](1), x_list[idK](0), x_list[idK](1));
                x << x_list[idI](0), x_list[idI](1), x_list[idJ](0), x_list[idJ](1), x_list[idK](0), x_list[idK](1);
                ceres::CostFunction* cost_function =
                  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_pi::Create(x);
                problem.AddResidualBlock(cost_function,
                  p_LossFunction,
                  &map_intrinsics[view1->id_intrinsic][0],
                  &map_poses[view1->id_pose][0],
                  //&map_intrinsics[view2->id_intrinsic][0],
                  &map_poses[view2->id_pose][0],
                  //&map_intrinsics[view3->id_intrinsic][0],
                  &map_poses[view3->id_pose][0]);

            }

        }

    }

  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_S::Create(obs_it.second.x);
//          IntrinsicsToCostFunction(
//            sfm_data.intrinsics[view->id_intrinsic].get(),
//            obs_it.second.x,
//            options.control_point_opt.weight);

        if (cost_function)
        {
            if (!map_intrinsics.at(view->id_intrinsic).empty())
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_intrinsics.at(view->id_intrinsic)[0],
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
            }
            else
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
  }
        }

      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
////    for (const auto & pose_it : sfm_data.GetPoses())
//    for (const auto & view_it : sfm_data.GetViews())
//    {
//      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
//      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
//      {
////          std::cout << "prior" << std::endl;
////          ceres::CostFunction * cost_function =
////            new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
////              new PoseCenterConstraintCostFunction( pose_it.second.center(), Vec3(1000.0,1000.0,1000.0)));//prior->pose_center_, prior->center_weight_));

////        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
//        ceres::CostFunction * cost_function =
//          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
//            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

//        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);

////        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[pose_it.first][0]);
////        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(4.0)), &map_poses[pose_it.first][0]);
////        problem.AddResidualBlock(cost_function, nullptr, &map_poses[pose_it.first][0]);

//      }
//    }
        for (const auto & pose_it : sfm_data.GetPoses())
        {

          {
//              std::cout << "prior" << std::endl;
              ceres::CostFunction * cost_function =
                new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
                  new PoseCenterConstraintCostFunction( pose_it.second.center(), Vec3(1.0,1.0,1.0)));//prior->pose_center_, prior->center_weight_));


//            problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[pose_it.first][0]);
            problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(4.0)), &map_poses[pose_it.first][0]);
//            problem.AddResidualBlock(cost_function, nullptr, &map_poses[pose_it.first][0]);

          }
        }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;

  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_ ;// / 1000000000;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " BriefReport : " << summary.BriefReport() << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & pose_it : sfm_data.GetPoses())
      {
        {
          X_SfM.push_back( pose_it.second.center() );
          X_GPS.push_back( pose_it.second.center());
        }
      }
//      for (const auto & view_it : sfm_data.GetViews())
//      {
//        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
//        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
//        {
//          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
//          X_GPS.push_back( prior->pose_center_ );
//        }
//      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}


bool Bundle_Adjustment_Ceres::Adjust_water_c4_S
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
//      for (const auto & view_it : sfm_data.GetViews())
      for (const auto & pose_it : sfm_data.GetPoses())
      {
//        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
//        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( pose_it.second.center() );
          X_GPS.push_back( pose_it.second.center());
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C4_S::Create(obs_it.second.x);
//        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
        else
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
      }
    }
    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C4_S::Create(obs_it.second.x);
//          IntrinsicsToCostFunction(
//            sfm_data.intrinsics[view->id_intrinsic].get(),
//            obs_it.second.x,
//            options.control_point_opt.weight);

        if (cost_function)
        {
            if (!map_intrinsics.at(view->id_intrinsic).empty())
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_intrinsics.at(view->id_intrinsic)[0],
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
            }
            else
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
  }
        }

      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & pose_it : sfm_data.GetPoses())
//    for (const auto & view_it : sfm_data.GetViews())
    {
//      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
//      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction( pose_it.second.center(), Vec3(1.0,1.0,1.0)));//prior->pose_center_, prior->center_weight_));

        std::cout << "prior " << std::endl;
//        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[pose_it.first][0]);
        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(4.0)), &map_poses[pose_it.first][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;

  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 20000;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = true;//ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;// / 100000;//00000;//10000000

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " BriefReport : " << summary.BriefReport() << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & pose_it : sfm_data.GetPoses())
      {
//        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
//        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( pose_it.second.center() );
          X_GPS.push_back( pose_it.second.center());
        }
      }
//      for (const auto & view_it : sfm_data.GetViews())
//      {
//        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
//        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
//        {
//          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
//          X_GPS.push_back( prior->pose_center_ );
//        }
//      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}


bool Bundle_Adjustment_Ceres::Adjust_step23
(
  sfm::SfM_Data & sfm_data,
  const std::map<int, bool>& fixedImgIdList,
  const std::map<int, bool>& fixedXIdList,
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    bool needFixedExtri = (fixedImgIdList.find(indexPose) != fixedImgIdList.end());
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE || needFixedExtri)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
        else
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
      }
    }
    bool needFixedX = (fixedXIdList.find(structure_landmark_it.first) != fixedXIdList.end());
    if (options.structure_opt == Structure_Parameter_Type::NONE || needFixedX)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          IntrinsicsToCostFunction(
            sfm_data.intrinsics[view->id_intrinsic].get(),
            obs_it.second.x,
            options.control_point_opt.weight);

        if (cost_function)
        {
            if (!map_intrinsics.at(view->id_intrinsic).empty())
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_intrinsics.at(view->id_intrinsic)[0],
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
            }
            else
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
  }
        }

      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;

  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " BriefReport : " << summary.BriefReport() << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}


bool Bundle_Adjustment_Ceres::Adjust_BAStep
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  std::cout << "sfm_data.poses" << std::endl;
  // Setup Poses data & subparametrization
  #pragma omp parallel for
  for (std::size_t tPoseListId = 0; tPoseListId < sfm_data.poses.size(); ++tPoseListId)
  {
      Poses::const_iterator pose_it = sfm_data.poses.begin();
      std::advance(pose_it, tPoseListId);

    const IndexT indexPose = pose_it->first;

    const Pose3 & pose = pose_it->second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    #pragma omp critical
    {
        // angleAxis + translation
        map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};
    }
        double * parameter_block = &map_poses[indexPose][0];
    #pragma omp critical
    {
        problem.AddParameterBlock(parameter_block, 6);
    }
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
        #pragma omp critical
        {
            // set the whole parameter block as constant for best performance
            problem.SetParameterBlockConstant(parameter_block);
        }
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
          #pragma omp critical
          {
              ceres::SubsetParameterization *subset_parameterization =
                new ceres::SubsetParameterization(6, vec_constant_extrinsic);
              problem.SetParameterization(parameter_block, subset_parameterization);
          }
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  Hash_Map<IndexT, std::vector<double> > X;
  std::cout << "sfm_data.structure" << std::endl;
  // For all visibility add reprojections errors:
  #pragma omp parallel for
  for (std::size_t tSListId = 0; tSListId < sfm_data.structure.size(); ++tSListId)
  {
      Landmarks::const_iterator structure_landmark_it = sfm_data.structure.begin();
      std::advance(structure_landmark_it, tSListId);
    const Observations & obs = structure_landmark_it->second.obs;
    X[structure_landmark_it->first] = {structure_landmark_it->second.X(0), structure_landmark_it->second.X(1), structure_landmark_it->second.X(2)};

    #pragma omp parallel for
    for (std::size_t tObsListId = 0; tObsListId < obs.size(); ++tObsListId)
    {
        Observations::const_iterator obs_it = obs.begin();
        std::advance(obs_it, tObsListId);
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it->first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it->second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
            #pragma omp critical
            {
                problem.AddResidualBlock(cost_function,
                    p_LossFunction,
                    &map_intrinsics[view->id_intrinsic][0],
                    &map_poses[view->id_pose][0],
                    &X[structure_landmark_it->first][0]);
//                    structure_landmark_it->second.X.data());
            }
        }
        else
        {
            #pragma omp critical
            {
                problem.AddResidualBlock(cost_function,
                    p_LossFunction,
                    &map_poses[view->id_pose][0],
                    &X[structure_landmark_it->first][0]);
//                    (structure_landmark_it->second).X.data());
            }
        }
      }
    }
    if (options.structure_opt == Structure_Parameter_Type::NONE)
    {
        #pragma omp critical
        {
//            problem.SetParameterBlockConstant(structure_landmark_it->second.X.data());
            problem.SetParameterBlockConstant(&X[structure_landmark_it->first][0]);
        }
    }
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          IntrinsicsToCostFunction(
            sfm_data.intrinsics[view->id_intrinsic].get(),
            obs_it.second.x,
            options.control_point_opt.weight);

        if (cost_function)
        {
            if (!map_intrinsics.at(view->id_intrinsic).empty())
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_intrinsics.at(view->id_intrinsic)[0],
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
            }
            else
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
  }
        }

      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;

  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 1;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " BriefReport : " << summary.BriefReport() << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      #pragma omp parallel for
      for (std::size_t tPoseListId = 0; tPoseListId < sfm_data.poses.size(); ++tPoseListId)
      {
          Poses::iterator pose_it = sfm_data.poses.begin();
          std::advance(pose_it, tPoseListId);

        const IndexT indexPose = pose_it->first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = (pose_it->second);
        #pragma omp critical
        {
            pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
        }

      }
    }

    #pragma omp parallel for
    for (std::size_t tSListId = 0; tSListId < sfm_data.structure.size(); ++tSListId)
    {
        Landmarks::iterator structure_landmark_it = sfm_data.structure.begin();
        std::advance(structure_landmark_it, tSListId);
        //
        #pragma omp critical
        {
            structure_landmark_it->second.X = Vec3(X[structure_landmark_it->first][0],
                                                  X[structure_landmark_it->first][1],
                                                  X[structure_landmark_it->first][2]);
        }

    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}


bool Bundle_Adjustment_Ceres::Adjust_water_prepare
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function;
      if(map_intrinsics[view->id_intrinsic].size() == 3)
      {
            cost_function = ResidualErrorFunctor_Pinhole_Intrinsic_water_prepare::Create(obs_it.second.x);
      }else{
            cost_function = ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_prepare::Create(obs_it.second.x);
      }


        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
        else
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
      }
    }
    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          IntrinsicsToCostFunction(
            sfm_data.intrinsics[view->id_intrinsic].get(),
            obs_it.second.x,
            options.control_point_opt.weight);

        if (cost_function)
        {
            if (!map_intrinsics.at(view->id_intrinsic).empty())
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_intrinsics.at(view->id_intrinsic)[0],
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
            }
            else
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
  }
        }

      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;

  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_ / 1000.0;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " BriefReport : " << summary.BriefReport() << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}



bool Bundle_Adjustment_Ceres::Adjust_water_c1
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
//        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        problem.AddParameterBlock(parameter_block, 3);// map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
//                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
                3, vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

//      if(obs_it.first / 100000 == 2)
//      {
//          continue;
//      }

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
              ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1::Create(obs_it.second.x);
//        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
        else
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
      }
    }

    {
//        problem.SetParameterBlockConstant(parameter_block);

//        double * parameter_block = &map_poses[indexPose][0];
//        problem.AddParameterBlock(parameter_block, 6);

//        std::vector<int> vec_constant_extrinsic;
//        // If we adjust only the translation, we must set ROTATION as constant
//        if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
//        {
//          // Subset rotation parametrization
//          vec_constant_extrinsic.push_back(0);
//          vec_constant_extrinsic.push_back(1);
//          vec_constant_extrinsic.push_back(2);
//        }
//        // If we adjust only the rotation, we must set TRANSLATION as constant
//        if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
//        {
//          // Subset translation parametrization
//          vec_constant_extrinsic.push_back(3);
//          vec_constant_extrinsic.push_back(4);
//          vec_constant_extrinsic.push_back(5);
//        }
//        if (!vec_constant_extrinsic.empty())
//        {
//          ceres::SubsetParameterization *subset_parameterization =
//            new ceres::SubsetParameterization(6, vec_constant_extrinsic);
//          problem.SetParameterization(parameter_block, subset_parameterization);
//        }
    }

//    std::vector<int> vec_constant_extrinsic;
//    vec_constant_extrinsic.push_back(2);
//    double * parameter_block = structure_landmark_it.second.X.data();
//    problem.AddParameterBlock(parameter_block, 3);
//    ceres::SubsetParameterization *subset_parameterization =
//      new ceres::SubsetParameterization(3, vec_constant_extrinsic);
//    problem.SetParameterization(parameter_block, subset_parameterization);

    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
                ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1::Create(obs_it.second.x, options.control_point_opt.weight);
//          IntrinsicsToCostFunction(
//            sfm_data.intrinsics[view->id_intrinsic].get(),
//            obs_it.second.x,
//            options.control_point_opt.weight);

        if (cost_function)
        {
            if (!map_intrinsics.at(view->id_intrinsic).empty())
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_intrinsics.at(view->id_intrinsic)[0],
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
            }
            else
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
  }
        }

      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;


  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;



//  ceres::Solver::Options ceres_config_options;
//  ceres_config_options.max_num_iterations = 500;
//  ceres_config_options.preconditioner_type =
//    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
//  ceres_config_options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
////    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
//  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
//    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
//  ceres_config_options.minimizer_progress_to_stdout = true;//ceres_options_.bVerbose_;
//  ceres_config_options.logging_type = ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
//  ceres_config_options.num_threads = ceres_options_.nb_threads_;
//  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
//  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << " BriefReport : " << summary.BriefReport() << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}


bool Bundle_Adjustment_Ceres::Adjust_water2_c1
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  Hash_Map<IndexT, std::vector<double> > X;
  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;
    X[structure_landmark_it.first] = {structure_landmark_it.second.X(0), structure_landmark_it.second.X(1)};//, structure_landmark_it.second.X(2)};

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
              ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_2::Create(obs_it.second.x);
//        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            &X[structure_landmark_it.first][0]);
//            structure_landmark_it.second.X.data());
        }
        else
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_poses[view->id_pose][0],
            &X[structure_landmark_it.first][0]);
//            structure_landmark_it.second.X.data());
        }
      }
    }

    {
//        problem.SetParameterBlockConstant(parameter_block);

//        double * parameter_block = &map_poses[indexPose][0];
//        problem.AddParameterBlock(parameter_block, 6);

//        std::vector<int> vec_constant_extrinsic;
//        // If we adjust only the translation, we must set ROTATION as constant
//        if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
//        {
//          // Subset rotation parametrization
//          vec_constant_extrinsic.push_back(0);
//          vec_constant_extrinsic.push_back(1);
//          vec_constant_extrinsic.push_back(2);
//        }
//        // If we adjust only the rotation, we must set TRANSLATION as constant
//        if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
//        {
//          // Subset translation parametrization
//          vec_constant_extrinsic.push_back(3);
//          vec_constant_extrinsic.push_back(4);
//          vec_constant_extrinsic.push_back(5);
//        }
//        if (!vec_constant_extrinsic.empty())
//        {
//          ceres::SubsetParameterization *subset_parameterization =
//            new ceres::SubsetParameterization(6, vec_constant_extrinsic);
//          problem.SetParameterization(parameter_block, subset_parameterization);
//        }
    }

//    std::vector<int> vec_constant_extrinsic;
//    vec_constant_extrinsic.push_back(2);
//    double * parameter_block = structure_landmark_it.second.X.data();
//    problem.AddParameterBlock(parameter_block, 3);
//    ceres::SubsetParameterization *subset_parameterization =
//      new ceres::SubsetParameterization(3, vec_constant_extrinsic);
//    problem.SetParameterization(parameter_block, subset_parameterization);

    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
                ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_2::Create(obs_it.second.x, options.control_point_opt.weight);
//          IntrinsicsToCostFunction(
//            sfm_data.intrinsics[view->id_intrinsic].get(),
//            obs_it.second.x,
//            options.control_point_opt.weight);

        if (cost_function)
        {
            if (!map_intrinsics.at(view->id_intrinsic).empty())
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_intrinsics.at(view->id_intrinsic)[0],
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
            }
            else
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
  }
        }

      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;


  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 5000;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;


//  ceres::Solver::Options ceres_config_options;
//  ceres_config_options.max_num_iterations = 500;
//  ceres_config_options.preconditioner_type =
//    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
//  ceres_config_options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
////    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
//  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
//    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
//  ceres_config_options.minimizer_progress_to_stdout = true;//ceres_options_.bVerbose_;
//  ceres_config_options.logging_type = ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
//  ceres_config_options.num_threads = ceres_options_.nb_threads_;
//  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
//  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << " BriefReport : " << summary.BriefReport() << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    for (auto & structure_landmark_it : sfm_data.structure)
    {
        structure_landmark_it.second.X = Vec3(X[structure_landmark_it.first][0], X[structure_landmark_it.first][1], structure_landmark_it.second.X(2));//X[structure_landmark_it.first][2]);//structure_landmark_it.second.X(2));//X[structure_landmark_it.first][2]);
//      const Observations & obs = structure_landmark_it.second.obs;
//      X[structure_landmark_it.first] = {structure_landmark_it.second.X(0), structure_landmark_it.second.X(1), structure_landmark_it.second.X(2)};
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}

bool Bundle_Adjustment_Ceres::Adjust_water3_c1
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;



  Hash_Map<IndexT, std::vector<double> > X;

  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;
    X[structure_landmark_it.first] = {structure_landmark_it.second.X(0), structure_landmark_it.second.X(1)};//, structure_landmark_it.second.X(2)};

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
              ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3::Create(obs_it.second.x);
//        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            &X[structure_landmark_it.first][0]);
//            structure_landmark_it.second.X.data());
        }
        else
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_poses[view->id_pose][0],
            &X[structure_landmark_it.first][0]);
//            structure_landmark_it.second.X.data());
        }
      }
    }

    {
//        problem.SetParameterBlockConstant(parameter_block);

//        double * parameter_block = &map_poses[indexPose][0];
//        problem.AddParameterBlock(parameter_block, 6);

//        std::vector<int> vec_constant_extrinsic;
//        // If we adjust only the translation, we must set ROTATION as constant
//        if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
//        {
//          // Subset rotation parametrization
//          vec_constant_extrinsic.push_back(0);
//          vec_constant_extrinsic.push_back(1);
//          vec_constant_extrinsic.push_back(2);
//        }
//        // If we adjust only the rotation, we must set TRANSLATION as constant
//        if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
//        {
//          // Subset translation parametrization
//          vec_constant_extrinsic.push_back(3);
//          vec_constant_extrinsic.push_back(4);
//          vec_constant_extrinsic.push_back(5);
//        }
//        if (!vec_constant_extrinsic.empty())
//        {
//          ceres::SubsetParameterization *subset_parameterization =
//            new ceres::SubsetParameterization(6, vec_constant_extrinsic);
//          problem.SetParameterization(parameter_block, subset_parameterization);
//        }
    }

//    std::vector<int> vec_constant_extrinsic;
//    vec_constant_extrinsic.push_back(2);
//    double * parameter_block = structure_landmark_it.second.X.data();
//    problem.AddParameterBlock(parameter_block, 3);
//    ceres::SubsetParameterization *subset_parameterization =
//      new ceres::SubsetParameterization(3, vec_constant_extrinsic);
//    problem.SetParameterization(parameter_block, subset_parameterization);

    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
                ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_3::Create(obs_it.second.x, options.control_point_opt.weight);
//          IntrinsicsToCostFunction(
//            sfm_data.intrinsics[view->id_intrinsic].get(),
//            obs_it.second.x,
//            options.control_point_opt.weight);

        if (cost_function)
        {
            if (!map_intrinsics.at(view->id_intrinsic).empty())
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_intrinsics.at(view->id_intrinsic)[0],
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
            }
            else
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
  }
        }

      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;


  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 1000;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;
  std::cout << "parameter_tolerance_: " << ceres_options_.parameter_tolerance_ << std::endl;


//  ceres::Solver::Options ceres_config_options;
//  ceres_config_options.max_num_iterations = 500;
//  ceres_config_options.preconditioner_type =
//    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
//  ceres_config_options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
////    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
//  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
//    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
//  ceres_config_options.minimizer_progress_to_stdout = true;//ceres_options_.bVerbose_;
//  ceres_config_options.logging_type = ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
//  ceres_config_options.num_threads = ceres_options_.nb_threads_;
//  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
//  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << " BriefReport : " << summary.BriefReport() << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    for (auto & structure_landmark_it : sfm_data.structure)
    {
        structure_landmark_it.second.X = Vec3(X[structure_landmark_it.first][0], X[structure_landmark_it.first][1], structure_landmark_it.second.X(2));//X[structure_landmark_it.first][2]);
//      const Observations & obs = structure_landmark_it.second.obs;
//      X[structure_landmark_it.first] = {structure_landmark_it.second.X(0), structure_landmark_it.second.X(1), structure_landmark_it.second.X(2)};
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}

bool Bundle_Adjustment_Ceres::Adjust_water5_c1
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  Hash_Map<IndexT, std::vector<double> > X;
  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;
    X[structure_landmark_it.first] = {structure_landmark_it.second.X(0), structure_landmark_it.second.X(1)};//, structure_landmark_it.second.X(2)};
    double * parameter_block = &X[structure_landmark_it.first][0];
    problem.AddParameterBlock(parameter_block, 2);

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
              ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_5::Create(obs_it.second.x);
//        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            &X[structure_landmark_it.first][0]);
//            structure_landmark_it.second.X.data());
        }
        else
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_poses[view->id_pose][0],
            &X[structure_landmark_it.first][0]);
//            structure_landmark_it.second.X.data());
        }
      }
    }

    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
                ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_5::Create(obs_it.second.x, options.control_point_opt.weight);
//          IntrinsicsToCostFunction(
//            sfm_data.intrinsics[view->id_intrinsic].get(),
//            obs_it.second.x,
//            options.control_point_opt.weight);

        if (cost_function)
        {
            if (!map_intrinsics.at(view->id_intrinsic).empty())
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_intrinsics.at(view->id_intrinsic)[0],
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
            }
            else
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
  }
        }

      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;


  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 5000;
  ceres_config_options.preconditioner_type = //ceres::CLUSTER_JACOBI;
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =  //ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type = //ceres::SUITE_SPARSE;//ceres::NO_SPARSE;//Eceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  std::cout << "ceres_options_.sparse_linear_algebra_library_type_ " << ceres_options_.sparse_linear_algebra_library_type_ << std::endl;
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;// / 10000.0;

//  ceres_config_options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;

  ceres_config_options.use_approximate_eigenvalue_bfgs_scaling = 1;
  std::cout << "ceres_config_options.use_approximate_eigenvalue_bfgs_scaling " << ceres_config_options.use_approximate_eigenvalue_bfgs_scaling << std::endl;



//  ceres::Solver::Options ceres_config_options;
//  ceres_config_options.max_num_iterations = 500;
//  ceres_config_options.preconditioner_type =
//    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
//  ceres_config_options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
////    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
//  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
//    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
//  ceres_config_options.minimizer_progress_to_stdout = true;//ceres_options_.bVerbose_;
//  ceres_config_options.logging_type = ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
//  ceres_config_options.num_threads = ceres_options_.nb_threads_;
//  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
//  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << " BriefReport : " << summary.BriefReport() << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    for (auto & structure_landmark_it : sfm_data.structure)
    {
        structure_landmark_it.second.X = Vec3(X[structure_landmark_it.first][0], X[structure_landmark_it.first][1], structure_landmark_it.second.X(2));//X[structure_landmark_it.first][2]);//structure_landmark_it.second.X(2));//X[structure_landmark_it.first][2]);
//      const Observations & obs = structure_landmark_it.second.obs;
//      X[structure_landmark_it.first] = {structure_landmark_it.second.X(0), structure_landmark_it.second.X(1), structure_landmark_it.second.X(2)};
    }
    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}

bool Bundle_Adjustment_Ceres::Adjust_water6_c1
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  Hash_Map<IndexT, std::vector<double> > X;
  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;
    X[structure_landmark_it.first] = {structure_landmark_it.second.X(0), structure_landmark_it.second.X(1), structure_landmark_it.second.X(2)};

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
              ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6::Create(obs_it.second.x);
//        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            &X[structure_landmark_it.first][0]);
//            structure_landmark_it.second.X.data());
        }
        else
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_poses[view->id_pose][0],
            &X[structure_landmark_it.first][0]);
//            structure_landmark_it.second.X.data());
        }
      }
    }

    if (options.structure_opt == Structure_Parameter_Type::NONE)
        problem.SetParameterBlockConstant(&X[structure_landmark_it.first][0]);
//      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
                ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_6::Create(obs_it.second.x, options.control_point_opt.weight);
//          IntrinsicsToCostFunction(
//            sfm_data.intrinsics[view->id_intrinsic].get(),
//            obs_it.second.x,
//            options.control_point_opt.weight);

        if (cost_function)
        {
            if (!map_intrinsics.at(view->id_intrinsic).empty())
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_intrinsics.at(view->id_intrinsic)[0],
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
            }
            else
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
  }
        }

      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;


  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 1000;//500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =// ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_ / 100000.0;

  ceres_config_options.gradient_tolerance /= 100.0;
  std::cout << "gradient_tolerance " << ceres_config_options.gradient_tolerance << std::endl;
  std::cout << "parameter_tolerance " << ceres_config_options.parameter_tolerance << std::endl;
//  std::cout << "ceres_config_options.trust_region_strategy_type  " << ceres_config_options.trust_region_strategy_type << std::endl;
//  ceres_config_options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
//  std::cout << "ceres_config_options.trust_region_strategy_type  " << ceres_config_options.trust_region_strategy_type << std::endl;


//  ceres::Solver::Options ceres_config_options;
//  ceres_config_options.max_num_iterations = 500;
//  ceres_config_options.preconditioner_type =
//    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
//  ceres_config_options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
////    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
//  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
//    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
//  ceres_config_options.minimizer_progress_to_stdout = true;//ceres_options_.bVerbose_;
//  ceres_config_options.logging_type = ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
//  ceres_config_options.num_threads = ceres_options_.nb_threads_;
//  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
//  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << " BriefReport : " << summary.BriefReport() << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    for (auto & structure_landmark_it : sfm_data.structure)
    {
        structure_landmark_it.second.X = Vec3(X[structure_landmark_it.first][0], X[structure_landmark_it.first][1], X[structure_landmark_it.first][2]);
//      const Observations & obs = structure_landmark_it.second.obs;
//      X[structure_landmark_it.first] = {structure_landmark_it.second.X(0), structure_landmark_it.second.X(1), structure_landmark_it.second.X(2)};
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}


bool Bundle_Adjustment_Ceres::Adjust_water4_c1
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  Hash_Map<IndexT, std::vector<double> > X;
  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;
    X[structure_landmark_it.first] = {structure_landmark_it.second.X(0), structure_landmark_it.second.X(1)};//, structure_landmark_it.second.X(2)};
    double * parameter_block = &X[structure_landmark_it.first][0];
    problem.AddParameterBlock(parameter_block, 2);

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
              ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_4::Create(obs_it.second.x);
//        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            &X[structure_landmark_it.first][0]);
//            structure_landmark_it.second.X.data());
        }
        else
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_poses[view->id_pose][0],
            &X[structure_landmark_it.first][0]);
//            structure_landmark_it.second.X.data());
        }
      }
    }
    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
                ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_4::Create(obs_it.second.x, options.control_point_opt.weight);
//          IntrinsicsToCostFunction(
//            sfm_data.intrinsics[view->id_intrinsic].get(),
//            obs_it.second.x,
//            options.control_point_opt.weight);

        if (cost_function)
        {
            if (!map_intrinsics.at(view->id_intrinsic).empty())
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_intrinsics.at(view->id_intrinsic)[0],
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
            }
            else
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
  }
        }

      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;


  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 5000;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =  //ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;// / 1000.0;// / 1000.0;

  ceres_config_options.use_approximate_eigenvalue_bfgs_scaling = 1;

  std::cout << "use_approximate_eigenvalue_bfgs_scaling " << ceres_config_options.use_approximate_eigenvalue_bfgs_scaling << std::endl;
  std::cout << "gradient_tolerance " << ceres_config_options.gradient_tolerance << std::endl;
  std::cout << "parameter_tolerance " << ceres_config_options.parameter_tolerance << std::endl;



//  ceres::Solver::Options ceres_config_options;
//  ceres_config_options.max_num_iterations = 500;
//  ceres_config_options.preconditioner_type =
//    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
//  ceres_config_options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
////    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
//  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
//    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
//  ceres_config_options.minimizer_progress_to_stdout = true;//ceres_options_.bVerbose_;
//  ceres_config_options.logging_type = ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
//  ceres_config_options.num_threads = ceres_options_.nb_threads_;
//  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
//  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << " BriefReport : " << summary.BriefReport() << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    for (auto & structure_landmark_it : sfm_data.structure)
    {
        structure_landmark_it.second.X = Vec3(X[structure_landmark_it.first][0], X[structure_landmark_it.first][1], structure_landmark_it.second.X(2));//X[structure_landmark_it.first][2]);
//      const Observations & obs = structure_landmark_it.second.obs;
//      X[structure_landmark_it.first] = {structure_landmark_it.second.X(0), structure_landmark_it.second.X(1), structure_landmark_it.second.X(2)};
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}


bool Bundle_Adjustment_Ceres::Adjust_water4_c1_2
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      if(view->id_pose / 100000 != 2)
      {
          continue;

      }
      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
              ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_4::Create(obs_it.second.x);
//        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
        else
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
      }
    }

//    std::vector<int> vec_constant_extrinsic;
//    vec_constant_extrinsic.push_back(2);
//    double * parameter_block = structure_landmark_it.second.X.data();
//    problem.AddParameterBlock(parameter_block, 3);
//    ceres::SubsetParameterization *subset_parameterization =
//      new ceres::SubsetParameterization(3, vec_constant_extrinsic);
//    problem.SetParameterization(parameter_block, subset_parameterization);

    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
                ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_4::Create(obs_it.second.x, options.control_point_opt.weight);
//          IntrinsicsToCostFunction(
//            sfm_data.intrinsics[view->id_intrinsic].get(),
//            obs_it.second.x,
//            options.control_point_opt.weight);

        if (cost_function)
        {
            if (!map_intrinsics.at(view->id_intrinsic).empty())
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_intrinsics.at(view->id_intrinsic)[0],
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
            }
            else
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
  }
        }

      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;


  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;


//  ceres::Solver::Options ceres_config_options;
//  ceres_config_options.max_num_iterations = 500;
//  ceres_config_options.preconditioner_type =
//    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
//  ceres_config_options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
////    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
//  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
//    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
//  ceres_config_options.minimizer_progress_to_stdout = true;//ceres_options_.bVerbose_;
//  ceres_config_options.logging_type = ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
//  ceres_config_options.num_threads = ceres_options_.nb_threads_;
//  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
//  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << " BriefReport : " << summary.BriefReport() << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;
        if(indexPose / 100000 != 2)
        {
            continue;

        }

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}



bool Bundle_Adjustment_Ceres::Adjust_water4_c4
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
              ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C4_1::Create(obs_it.second.x);
//        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
        else
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
      }
    }

//    std::vector<int> vec_constant_extrinsic;
//    vec_constant_extrinsic.push_back(2);
//    double * parameter_block = structure_landmark_it.second.X.data();
//    problem.AddParameterBlock(parameter_block, 3);
//    ceres::SubsetParameterization *subset_parameterization =
//      new ceres::SubsetParameterization(3, vec_constant_extrinsic);
//    problem.SetParameterization(parameter_block, subset_parameterization);

    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
                ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C4_1::Create(obs_it.second.x, options.control_point_opt.weight);
//          IntrinsicsToCostFunction(
//            sfm_data.intrinsics[view->id_intrinsic].get(),
//            obs_it.second.x,
//            options.control_point_opt.weight);

        if (cost_function)
        {
            if (!map_intrinsics.at(view->id_intrinsic).empty())
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_intrinsics.at(view->id_intrinsic)[0],
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
            }
            else
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
  }
        }

      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;


  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;


//  ceres::Solver::Options ceres_config_options;
//  ceres_config_options.max_num_iterations = 500;
//  ceres_config_options.preconditioner_type =
//    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
//  ceres_config_options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
////    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
//  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
//    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
//  ceres_config_options.minimizer_progress_to_stdout = true;//ceres_options_.bVerbose_;
//  ceres_config_options.logging_type = ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
//  ceres_config_options.num_threads = ceres_options_.nb_threads_;
//  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
//  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << " BriefReport : " << summary.BriefReport() << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}

bool Bundle_Adjustment_Ceres::Adjust_6
(
  const std::map<int, bool> & controlImgMap,
  const std::map<int, bool> & controlStructureMap,
  const std::map<int, bool> & computeImageIdMap,
  std::map<int, bool> & computeStrctureIdMap,
  Landmarks & newStructure,
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
//      std::cout << "error 01" << std::endl;
    const IndexT indexPose = pose_it.first;

    bool inControlLiar = (controlImgMap.find(indexPose) != controlImgMap.end());
    if(computeImageIdMap.find(indexPose) == computeImageIdMap.end()
    && (!inControlLiar))
    {
        continue;
    }

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if(inControlLiar)
    {
        problem.SetParameterBlockConstant(parameter_block);
        continue;
    }
    //
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  //find fix intrinsics
  std::map<int, bool> fixIList;
  for(std::size_t idCI = 0; idCI < controlImgMap.size(); ++idCI)
  {
      std::map<int, bool>::const_iterator itCI = controlImgMap.begin();
      std::advance(itCI, idCI);
      //
      if(sfm_data.views.find(itCI->first) == sfm_data.views.end())
      {
          continue;
      }
      const View * view = sfm_data.views.at(itCI->first).get();
      if(fixIList.find(view->id_intrinsic) == fixIList.end())
      {
          fixIList.insert(std::pair<int, bool>(view->id_intrinsic, true));
      }

  }
  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
//      std::cout << "error 02" << std::endl;
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if(fixIList.find(indexCam) != fixIList.end())
        {
            problem.SetParameterBlockConstant(parameter_block);
        }
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

//  std::cout << "error 021" << std::endl;
  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {

//      std::cout << "error 03" << std::endl;
    const Observations & obs = structure_landmark_it.second.obs;
    const int XId = structure_landmark_it.first;

//    std::cout << "error 03 1" << std::endl;
    bool needSet = false;
    int obsImgCount = 0;
//    std::cout << "error 03 2" << std::endl;
    #pragma omp parallel for
    for (std::size_t idObs = 0; idObs < obs.size(); ++idObs)
    {
        Observations::const_iterator itObs = obs.begin();
        std::advance(itObs, idObs);
        if(sfm_data.views.find(itObs->first) == sfm_data.views.end())
        {
            continue;
        }
        const View * view = sfm_data.views.at(itObs->first).get();
        if(computeImageIdMap.find(view->id_pose) != computeImageIdMap.end()
        || controlImgMap.find(view->id_pose) != controlImgMap.end())
        {
            #pragma omp critical
            {
                ++obsImgCount;
            }
        }
    }
//    std::cout << "error 03 3" << std::endl;
    if(obsImgCount < 3 )
    {
        continue;
    }else{
        needSet = true;
    }
//    std::cout << "error 04" << std::endl;
    //
    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      if(sfm_data.views.find(obs_it.first) == sfm_data.views.end())
      {
          continue;
      }
      const View * view = sfm_data.views.at(obs_it.first).get();
      if(computeImageIdMap.find(view->id_pose) == computeImageIdMap.end()
      && controlImgMap.find(view->id_pose) == controlImgMap.end())
      {
          continue;
      }
      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
        else
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
      }
    }
    //
    if(controlStructureMap.find(XId) != controlStructureMap.end())
    {
        problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
    }else if(computeStrctureIdMap.find(XId) == computeStrctureIdMap.end())
    {
        computeStrctureIdMap.insert(std::pair<int, bool>(XId, true));
    }
    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          IntrinsicsToCostFunction(
            sfm_data.intrinsics[view->id_intrinsic].get(),
            obs_it.second.x,
            options.control_point_opt.weight);

        if (cost_function)
        {
            if (!map_intrinsics.at(view->id_intrinsic).empty())
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_intrinsics.at(view->id_intrinsic)[0],
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
            }
            else
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
            }
        }

      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;

  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;
        if(computeImageIdMap.find(indexPose) == computeImageIdMap.end()
        && controlImgMap.find(indexPose) == controlImgMap.end())
        {
            continue;
        }

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)
    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }

    {
        for (auto & structure_landmark_it : sfm_data.structure)
        {
          if(computeStrctureIdMap.find(structure_landmark_it.first) == computeStrctureIdMap.end()
          && controlStructureMap.find(structure_landmark_it.first) == controlStructureMap.end() )
          {
              continue;
          }
          newStructure[structure_landmark_it.first] = structure_landmark_it.second;
        }
    }

    return true;
  }
}


bool Bundle_Adjustment_Ceres::Adjust_6_water
(
  const std::map<int, bool> & controlImgMap,
  const std::map<int, bool> & controlStructureMap,
  const std::map<int, bool> & computeImageIdMap,
  std::map<int, bool> & computeStrctureIdMap,
  Landmarks & newStructure,
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
//      std::cout << "error 01" << std::endl;
    const IndexT indexPose = pose_it.first;

    bool inControlLiar = (controlImgMap.find(indexPose) != controlImgMap.end());
    if(computeImageIdMap.find(indexPose) == computeImageIdMap.end()
    && (!inControlLiar))
    {
        continue;
    }

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if(inControlLiar)
    {
        problem.SetParameterBlockConstant(parameter_block);
        continue;
    }
    //
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  //find fix intrinsics
  std::map<int, bool> fixIList;
  for(std::size_t idCI = 0; idCI < controlImgMap.size(); ++idCI)
  {
      std::map<int, bool>::const_iterator itCI = controlImgMap.begin();
      std::advance(itCI, idCI);
      //
      if(sfm_data.views.find(itCI->first) == sfm_data.views.end())
      {
          continue;
      }
      const View * view = sfm_data.views.at(itCI->first).get();
      if(fixIList.find(view->id_intrinsic) == fixIList.end())
      {
          fixIList.insert(std::pair<int, bool>(view->id_intrinsic, true));
      }

  }
  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
//      std::cout << "error 02" << std::endl;
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if(fixIList.find(indexCam) != fixIList.end())
        {
            problem.SetParameterBlockConstant(parameter_block);
        }
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

//  std::cout << "error 021" << std::endl;
  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {

//      std::cout << "error 03" << std::endl;
    const Observations & obs = structure_landmark_it.second.obs;
    const int XId = structure_landmark_it.first;

//    std::cout << "error 03 1" << std::endl;
    bool needSet = false;
    int obsImgCount = 0;
//    std::cout << "error 03 2" << std::endl;
    #pragma omp parallel for
    for (std::size_t idObs = 0; idObs < obs.size(); ++idObs)
    {
        Observations::const_iterator itObs = obs.begin();
        std::advance(itObs, idObs);
        if(sfm_data.views.find(itObs->first) == sfm_data.views.end())
        {
            continue;
        }
        const View * view = sfm_data.views.at(itObs->first).get();
        if(computeImageIdMap.find(view->id_pose) != computeImageIdMap.end()
        || controlImgMap.find(view->id_pose) != controlImgMap.end())
        {
            #pragma omp critical
            {
                ++obsImgCount;
            }
        }
    }
//    std::cout << "error 03 3" << std::endl;
    if(obsImgCount < 3 )
    {
        continue;
    }else{
        needSet = true;
    }
//    std::cout << "error 04" << std::endl;
    //
    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      if(sfm_data.views.find(obs_it.first) == sfm_data.views.end())
      {
          continue;
      }
      const View * view = sfm_data.views.at(obs_it.first).get();
      if(computeImageIdMap.find(view->id_pose) == computeImageIdMap.end()
      && controlImgMap.find(view->id_pose) == controlImgMap.end())
      {
          continue;
      }
      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
        else
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
      }
    }
    //
    if(controlStructureMap.find(XId) != controlStructureMap.end())
    {
        problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
    }else if(computeStrctureIdMap.find(XId) == computeStrctureIdMap.end())
    {
        computeStrctureIdMap.insert(std::pair<int, bool>(XId, true));
    }
    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          IntrinsicsToCostFunction(
            sfm_data.intrinsics[view->id_intrinsic].get(),
            obs_it.second.x,
            options.control_point_opt.weight);

        if (cost_function)
        {
            if (!map_intrinsics.at(view->id_intrinsic).empty())
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_intrinsics.at(view->id_intrinsic)[0],
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
            }
            else
            {
              problem.AddResidualBlock(cost_function,
                                       nullptr,
                                       &map_poses.at(view->id_pose)[0],
                                       gcp_landmark_it.second.X.data());
            }
        }

      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;

  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;
        if(computeImageIdMap.find(indexPose) == computeImageIdMap.end()
        && controlImgMap.find(indexPose) == controlImgMap.end())
        {
            continue;
        }

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)
    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }

    {
        for (auto & structure_landmark_it : sfm_data.structure)
        {
          if(computeStrctureIdMap.find(structure_landmark_it.first) == computeStrctureIdMap.end()
          && controlStructureMap.find(structure_landmark_it.first) == controlStructureMap.end() )
          {
              continue;
          }
          newStructure[structure_landmark_it.first] = structure_landmark_it.second;
        }
    }

    return true;
  }
}


bool Bundle_Adjustment_Ceres::Adjust_fix_some
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options,
  std::map<int, Pose3>& fixIdList
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    if(fixIdList.find(indexPose) == fixIdList.end())
    {
        const Pose3 & pose = pose_it.second;
        const Mat3 R = pose.rotation();
        const Vec3 t = pose.translation();

        double angleAxis[3];
        ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
        // angleAxis + translation
        map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

        double * parameter_block = &map_poses[indexPose][0];
        problem.AddParameterBlock(parameter_block, 6);
        if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else  // Subset parametrization
        {
          std::vector<int> vec_constant_extrinsic;
          // If we adjust only the translation, we must set ROTATION as constant
          if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
          {
            // Subset rotation parametrization
            vec_constant_extrinsic.push_back(0);
            vec_constant_extrinsic.push_back(1);
            vec_constant_extrinsic.push_back(2);
          }
          // If we adjust only the rotation, we must set TRANSLATION as constant
          if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
          {
            // Subset translation parametrization
            vec_constant_extrinsic.push_back(3);
            vec_constant_extrinsic.push_back(4);
            vec_constant_extrinsic.push_back(5);
          }
          if (!vec_constant_extrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(6, vec_constant_extrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }

    }else{
        const Pose3 & pose = fixIdList.find(indexPose)->second;
        const Mat3 R = pose.rotation();
        const Vec3 t = pose.translation();

        double angleAxis[3];
        ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
        // angleAxis + translation
        map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

        double * parameter_block = &map_poses[indexPose][0];
        problem.AddParameterBlock(parameter_block, 6);
//        if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
    }

  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
        else
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
      }
    }
    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          IntrinsicsToCostFunction(
            sfm_data.intrinsics[view->id_intrinsic].get(),
            obs_it.second.x,
            options.control_point_opt.weight);

        if (cost_function)
          problem.AddResidualBlock(cost_function,
            nullptr,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            gcp_landmark_it.second.X.data());
      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;

  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}

bool Bundle_Adjustment_Ceres::Adjust_down
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    if((indexPose/100000)!= 2)
    {
        continue;
    }

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      if(((view->id_pose)/100000)!= 2)
      {
          continue;
      }

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
        else
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
      }
    }
    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          IntrinsicsToCostFunction(
            sfm_data.intrinsics[view->id_intrinsic].get(),
            obs_it.second.x,
            options.control_point_opt.weight);

        if (cost_function)
          problem.AddResidualBlock(cost_function,
            nullptr,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            gcp_landmark_it.second.X.data());
      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;

  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        if((indexPose/100000)!= 2)
        {
            continue;
        }

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}

bool Bundle_Adjustment_Ceres::Adjust_front
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    if((indexPose/100000)!= 1)
    {
        continue;
    }

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      if(((view->id_pose)/100000)!= 1)
      {
          continue;
      }

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
        else
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
      }
    }
    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          IntrinsicsToCostFunction(
            sfm_data.intrinsics[view->id_intrinsic].get(),
            obs_it.second.x,
            options.control_point_opt.weight);

        if (cost_function)
          problem.AddResidualBlock(cost_function,
            nullptr,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            gcp_landmark_it.second.X.data());
      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;

  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        if((indexPose/100000)!= 2)
        {
            continue;
        }

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}


bool Bundle_Adjustment_Ceres::Adjust_back
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    if((indexPose/100000)!= 3)
    {
        continue;
    }

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      if(((view->id_pose)/100000)!= 3)
      {
          continue;
      }

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
        else
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
      }
    }
    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          IntrinsicsToCostFunction(
            sfm_data.intrinsics[view->id_intrinsic].get(),
            obs_it.second.x,
            options.control_point_opt.weight);

        if (cost_function)
          problem.AddResidualBlock(cost_function,
            nullptr,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            gcp_landmark_it.second.X.data());
      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;

  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        if((indexPose/100000)!= 2)
        {
            continue;
        }

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}

bool Bundle_Adjustment_Ceres::Adjust_down_front
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    if((indexPose/100000) == 3)
    {
        continue;
    }

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      if(((view->id_pose)/100000) == 3)
      {
          continue;
      }

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
        else
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
      }
    }
    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          IntrinsicsToCostFunction(
            sfm_data.intrinsics[view->id_intrinsic].get(),
            obs_it.second.x,
            options.control_point_opt.weight);

        if (cost_function)
          problem.AddResidualBlock(cost_function,
            nullptr,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            gcp_landmark_it.second.X.data());
      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;

  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        if((indexPose/100000) == 3)
        {
            continue;
        }

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}



bool Bundle_Adjustment_Ceres::Adjust_fix2
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options & options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------


  double pose_center_robust_fitting_error = 0.0;
  openMVG::geometry::Similarity3 sim_to_center;
  bool b_usable_prior = false;
  if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
  {
    // - Compute a robust X-Y affine transformation & apply it
    // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
    {
      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      openMVG::geometry::Similarity3 sim;

      // Compute the registration:
      if (X_GPS.size() > 3)
      {
        const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
        const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
        geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
        const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
        if (lmeds_median != std::numeric_limits<double>::max())
        {
          b_usable_prior = true; // PRIOR can be used safely

          // Compute the median residual error once the registration is applied
          for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
          {
            pos = sim(pos);
          }
          Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
          std::sort(residual.data(), residual.data() + residual.size());
          pose_center_robust_fitting_error = residual(residual.size()/2);

          // Apply the found transformation to the SfM Data Scene
          openMVG::sfm::ApplySimilarity(sim, sfm_data);

          // Move entire scene to center for better numerical stability
          Vec3 pose_centroid = Vec3::Zero();
          for (const auto & pose_it : sfm_data.poses)
          {
            pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
          }
          sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
          openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
        }
      }
    }
  }

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if(indexPose/100000 == 2)
    {
        problem.SetParameterBlockConstant(parameter_block);
    }
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
      {
        if (!map_intrinsics[view->id_intrinsic].empty())
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
        else
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_poses[view->id_pose][0],
            structure_landmark_it.second.X.data());
        }
      }
    }
    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          IntrinsicsToCostFunction(
            sfm_data.intrinsics[view->id_intrinsic].get(),
            obs_it.second.x,
            options.control_point_opt.weight);

        if (cost_function)
          problem.AddResidualBlock(cost_function,
            nullptr,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            gcp_landmark_it.second.X.data());
      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Add Pose prior constraints if any
  if (b_usable_prior)
  {
    for (const auto & view_it : sfm_data.GetViews())
    {
      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
      if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
      {
        // Add the cost functor (distance from Pose prior to the SfM_Data Pose center)
        ceres::CostFunction * cost_function =
          new ceres::AutoDiffCostFunction<PoseCenterConstraintCostFunction, 3, 6>(
            new PoseCenterConstraintCostFunction(prior->pose_center_, prior->center_weight_));

        problem.AddResidualBlock(cost_function, new ceres::HuberLoss(Square(pose_center_robust_fitting_error)), &map_poses[prior->id_view][0]);
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;

  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =// ceres::SUITE_SPARSE;
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
      if (options.use_motion_priors_opt)
        std::cout << "Usable motion priors: " << (int)b_usable_prior << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)

    if (b_usable_prior)
    {
      // set back to the original scene centroid
      openMVG::sfm::ApplySimilarity(sim_to_center.inverse(), sfm_data, true);

      //--
      // - Compute some fitting statistics
      //--

      // Collect corresponding camera centers
      std::vector<Vec3> X_SfM, X_GPS;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
        {
          X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
          X_GPS.push_back( prior->pose_center_ );
        }
      }
      // Compute the registration fitting error (once BA with Prior have been used):
      if (X_GPS.size() > 3)
      {
        // Compute the median residual error
        Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
        std::cout
          << "Pose prior statistics (user units):\n"
          << " - Starting median fitting error: " << pose_center_robust_fitting_error << "\n"
          << " - Final fitting error:";
        minMaxMeanMedian<Vec::Scalar>(residual.data(), residual.data() + residual.size());
      }
    }
    return true;
  }
}

bool Bundle_Adjustment_Ceres::Adjust_old
(
 SfM_Data & sfm_data,     // the SfM scene to refine
 const Optimize_Options & options
)
{

    //----------
    // Add camera parameters
    // - intrinsics
    // - poses [R|t]

    // Create residuals for each observation in the bundle adjustment problem. The
    // parameters for cameras and points are added automatically.
    //----------

    ceres::Problem problem;

    // Data wrapper for refinement:
    Hash_Map<IndexT, std::vector<double> > map_intrinsics;
    Hash_Map<IndexT, std::vector<double> > map_poses;

    // Setup Poses data & subparametrization
   for (Poses::const_iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose)
    {
      const IndexT indexPose = itPose->first;

      const Pose3 & pose = itPose->second;
      const Mat3 R = pose.rotation();
      const Vec3 t = pose.translation();

      double angleAxis[3];
      ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
      // angleAxis + translation
      map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

      double * parameter_block = &map_poses[indexPose][0];
      problem.AddParameterBlock(parameter_block, 6);
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
      {
        // set the whole parameter block as constant for best performance
        problem.SetParameterBlockConstant(parameter_block);
      }
      else  // Subset parametrization
      {
        std::vector<int> vec_constant_extrinsic;
        // If we adjust only the translation, we must set ROTATION as constant
        if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
        {
          // Subset rotation parametrization
          vec_constant_extrinsic.push_back(0);
          vec_constant_extrinsic.push_back(1);
          vec_constant_extrinsic.push_back(2);
        }
        // If we adjust only the rotation, we must set TRANSLATION as constant
        if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
        {
          // Subset translation parametrization
          vec_constant_extrinsic.push_back(3);
          vec_constant_extrinsic.push_back(4);
          vec_constant_extrinsic.push_back(5);
        }
        if (!vec_constant_extrinsic.empty())
        {
          ceres::SubsetParameterization *subset_parameterization =
            new ceres::SubsetParameterization(6, vec_constant_extrinsic);
          problem.SetParameterization(parameter_block, subset_parameterization);
        }
      }
    }

    // Setup Intrinsics data & subparametrization
    for (Intrinsics::const_iterator itIntrinsic = sfm_data.intrinsics.begin();
      itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
    {
      const IndexT indexCam = itIntrinsic->first;

      if (isValid(itIntrinsic->second->getType()))
      {
        map_intrinsics[indexCam] = itIntrinsic->second->getParams();

        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            itIntrinsic->second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
      else
      {
        std::cerr << "Unsupported camera type." << std::endl;
      }
    }

    // Set a LossFunction to be less penalized by false measurements
    //  - set it to NULL if you don't want use a lossFunction.
    ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(Square(4.0));
    // TODO: make the LOSS function and the parameter an option

    // For all visibility add reprojections errors:
    for (Landmarks::iterator iterTracks = sfm_data.structure.begin();
      iterTracks!= sfm_data.structure.end(); ++iterTracks)
    {
      const Observations & obs = iterTracks->second.obs;

      for (Observations::const_iterator itObs = obs.begin();
        itObs != obs.end(); ++itObs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(itObs->first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), itObs->second.x);

        if (cost_function)
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            iterTracks->second.X.data());
      }
      if (options.structure_opt == Structure_Parameter_Type::NONE)
        problem.SetParameterBlockConstant(iterTracks->second.X.data());
    }

    if (options.control_point_opt.bUse_control_points)
    {
      // Use Ground Control Point:
      // - fixed 3D points with weighted observations
      for (Landmarks::iterator iterGCPTracks = sfm_data.control_points.begin();
        iterGCPTracks!= sfm_data.control_points.end(); ++iterGCPTracks)
      {
        const Observations & obs = iterGCPTracks->second.obs;

        for (Observations::const_iterator itObs = obs.begin();
          itObs != obs.end(); ++itObs)
        {
          // Build the residual block corresponding to the track observation:
          const View * view = sfm_data.views.at(itObs->first).get();

          // Each Residual block takes a point and a camera as input and outputs a 2
          // dimensional residual. Internally, the cost function stores the observed
          // image location and compares the reprojection against the observation.
          ceres::CostFunction* cost_function =
            IntrinsicsToCostFunction(
              sfm_data.intrinsics[view->id_intrinsic].get(),
              itObs->second.x,
              options.control_point_opt.weight);

          if (cost_function)
            problem.AddResidualBlock(cost_function,
              nullptr,
              &map_intrinsics[view->id_intrinsic][0],
              &map_poses[view->id_pose][0],
              iterGCPTracks->second.X.data());
        }
        // Set the 3D point as FIXED (it's a GCP)
        problem.SetParameterBlockConstant(iterGCPTracks->second.X.data());
      }
    }


    // Configure a BA engine and run it
    //  Make Ceres automatically detect the bundle structure.
    ceres::Solver::Options ceres_config_options;
    ceres_config_options.max_num_iterations = 500;
    ceres_config_options.preconditioner_type =
      static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
    ceres_config_options.linear_solver_type =
      static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
    ceres_config_options.sparse_linear_algebra_library_type =
      static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
    ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
    ceres_config_options.logging_type = ceres::SILENT;
    ceres_config_options.num_threads = ceres_options_.nb_threads_;
    ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
    ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;


    // Solve BA
    ceres::Solver::Summary summary;
    ceres::Solve(ceres_config_options, &problem, &summary);
    if (ceres_options_.bCeres_summary_)
      std::cout << summary.FullReport() << std::endl;

    // If no error, get back refined parameters
    if (!summary.IsSolutionUsable())
    {
      if (ceres_options_.bVerbose_)
        std::cout << "Bundle Adjustment failed." << std::endl;
      return false;
    }
    else // Solution is usable
    {
      if (ceres_options_.bVerbose_)
      {
        // Display statistics about the minimization
        std::cout << std::endl
          << "Bundle Adjustment statistics (approximated RMSE):\n"
          << " #views: " << sfm_data.views.size() << "\n"
          << " #poses: " << sfm_data.poses.size() << "\n"
          << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
          << " #tracks: " << sfm_data.structure.size() << "\n"
          << " #residuals: " << summary.num_residuals << "\n"
          << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
          << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
          << " Time (s): " << summary.total_time_in_seconds << "\n"
          << std::endl;
      }

      // Update camera poses with refined data
      if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
      {
        for (Poses::iterator itPose = sfm_data.poses.begin();
          itPose != sfm_data.poses.end(); ++itPose)
        {
          const IndexT indexPose = itPose->first;

          Mat3 R_refined;
          ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
          Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
          // Update the pose
          Pose3 & pose = itPose->second;
          pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
        }
      }

      // Update camera intrinsics with refined data
      if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
      {
        for (Intrinsics::iterator itIntrinsic = sfm_data.intrinsics.begin();
          itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
        {
          const IndexT indexCam = itIntrinsic->first;

          const std::vector<double> & vec_params = map_intrinsics[indexCam];
          itIntrinsic->second.get()->updateFromParams(vec_params);
        }
      }
      // Structure is already updated directly if needed (no data wrapping)
      return true;
    }

}



bool
Bundle_Adjustment_Ceres::filiter_ununited_pose(Poses& poses, Views& views, std::map<int, int>& resMap )
{
    //initialize resMap first using views id
    for (Views::const_iterator it = views.begin(); it != views.end(); ++it)
    {
        std::pair<int, int> tempPair;
        tempPair.first = it->first;
        tempPair.second = 0;
        resMap.insert(tempPair);
    }

    //build hash table
    for (Poses::const_iterator itPose = poses.begin(); itPose != poses.end(); ++itPose)
    {
        std::map<int, int>::iterator itMap = resMap.find(itPose->first);
        itMap->second = 1;
    }

    for (std::map<int, int>::iterator it = resMap.begin();
         it != resMap.end(); ++it)
    {
        if(it->second == 0)
        {
            int poseId = it->first;

            if(poseId < views.size()/3)  //back
            {
                resMap.find(poseId + views.size()/3)->second = 0;
                resMap.find(poseId + views.size()/3*2)->second = 0;

            }else if(poseId >= views.size()/3*2)  //front
            {
                resMap.find(poseId - views.size()/3)->second = 0;
                resMap.find(poseId - views.size()/3*2)->second = 0;

            }else{  //down
                resMap.find(poseId + views.size()/3)->second = 0;
                resMap.find(poseId - views.size()/3)->second = 0;
            }
        }

    }

}

bool
Bundle_Adjustment_Ceres::filiter_ununited_pose_test(Poses& poses, Views& views, std::map<int, int>& resMap )
{
    //initialize resMap first using views id
    for (Views::const_iterator it = views.begin(); it != views.end(); ++it)
    {
        std::pair<int, int> tempPair;
        tempPair.first = it->first;
        tempPair.second = 0;
        resMap.insert(tempPair);
    }

    //build hash table
    for (Poses::const_iterator itPose = poses.begin(); itPose != poses.end(); ++itPose)
    {
        std::map<int, int>::iterator itMap = resMap.find(itPose->first);
        itMap->second = 1;
    }

    for (std::map<int, int>::iterator it = resMap.begin();
         it != resMap.end(); ++it)
    {
        {
            int poseId = it->first;

            if(poseId < views.size()/3)  //back
            {
                resMap.find(poseId + views.size()/3)->second = 0;
                resMap.find(poseId + views.size()/3*2)->second = 0;

            }else if(poseId >= views.size()/3*2)  //front
            {
                resMap.find(poseId - views.size()/3)->second = 0;
                resMap.find(poseId - views.size()/3*2)->second = 0;

            }else{  //down
                resMap.find(poseId + views.size()/3)->second = 0;
                resMap.find(poseId - views.size()/3)->second = 0;
            }
        }

    }

}


bool Bundle_Adjustment_Ceres::Adjust_single_fix_fuv
(SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options options

)
{
        //----------
        // Add camera parameters
        // - intrinsics
        // - poses [R|t]

        // Create residuals for each observation in the bundle adjustment problem. The
        // parameters for cameras and points are added automatically.
        //----------

        ceres::Problem problem;

        // Data wrapper for refinement:
        Hash_Map<IndexT, std::vector<double> > map_intrinsics;
        Hash_Map<IndexT, std::vector<double> > map_poses;


        //std::cout << "error 1" << std::endl;
        // Setup Poses data & subparametrization
       for (Poses::const_iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose)
        {
          const IndexT indexPose = itPose->first;

          const Pose3 & pose = itPose->second;
          const Mat3 R = pose.rotation();
          const Vec3 t = pose.translation();

          double angleAxis[3];
          ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
          // angleAxis + translation
          map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

          double * parameter_block = &map_poses[indexPose][0];
          problem.AddParameterBlock(parameter_block, map_poses[indexPose].size());
          if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
          {
            // set the whole parameter block as constant for best performance
            problem.SetParameterBlockConstant(parameter_block);
          }
          else  // Subset parametrization
          {
            std::vector<int> vec_constant_extrinsic;
            // If we adjust only the translation, we must set ROTATION as constant
            if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
            {
              // Subset rotation parametrization
              vec_constant_extrinsic.push_back(0);
              vec_constant_extrinsic.push_back(1);
              vec_constant_extrinsic.push_back(2);
            }
            // If we adjust only the rotation, we must set TRANSLATION as constant
            if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
            {
              // Subset translation parametrization
              vec_constant_extrinsic.push_back(3);
              vec_constant_extrinsic.push_back(4);
              vec_constant_extrinsic.push_back(5);
            }
            if (!vec_constant_extrinsic.empty())
            {
              ceres::SubsetParameterization *subset_parameterization =
                new ceres::SubsetParameterization(6, vec_constant_extrinsic);
              problem.SetParameterization(parameter_block, subset_parameterization);
            }
          }
        }

       //std::cout << "error 2" << std::endl;
       // Setup Intrinsics data & subparametrization
        for (Intrinsics::const_iterator itIntrinsic = sfm_data.intrinsics.begin();
          itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
        {
          const IndexT indexCam = itIntrinsic->first;

          if (isValid(itIntrinsic->second->getType()))
          {
            map_intrinsics[indexCam] = itIntrinsic->second->getParams();

            double * parameter_block = &map_intrinsics[indexCam][0];
            problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());


            if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
            {
              // set the whole parameter block as constant for best performance
              problem.SetParameterBlockConstant(parameter_block);
            }
            else
            {
              const std::vector<int> vec_constant_intrinsic =
                itIntrinsic->second->subsetParameterization(options.intrinsics_opt);
              if (!vec_constant_intrinsic.empty())
              {
                ceres::SubsetParameterization *subset_parameterization =
                  new ceres::SubsetParameterization(
                    map_intrinsics[indexCam].size(), vec_constant_intrinsic);
                problem.SetParameterization(parameter_block, subset_parameterization);
              }

//              //fix fuv
//              std::vector<int> vec_constant_fuv;
//              vec_constant_fuv.push_back(3);
//              vec_constant_fuv.push_back(4);
//              vec_constant_fuv.push_back(5);
//              vec_constant_fuv.push_back(6);
//              vec_constant_fuv.push_back(7);
//              ceres::SubsetParameterization *subset_parameterization =
//                new ceres::SubsetParameterization(map_intrinsics[indexCam].size(), vec_constant_fuv);
//              problem.SetParameterization(parameter_block, subset_parameterization);

            }
          }
          else
          {
            std::cerr << "Unsupported camera type." << std::endl;
          }
        }





        // Set a LossFunction to be less penalized by false measurements
        //  - set it to NULL if you don't want use a lossFunction.
        ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(Square(4.0));
        // TODO: make the LOSS function and the parameter an option

        //std::cout << "error sfm_data.structure size" << sfm_data.structure.size() << std::endl;
        // For all visibility add reprojections errors:
        for (Landmarks::iterator iterTracks = sfm_data.structure.begin();
          iterTracks!= sfm_data.structure.end(); ++iterTracks)
        {
          const Observations & obs = iterTracks->second.obs;

          for (Observations::const_iterator itObs = obs.begin();
            itObs != obs.end(); ++itObs)
          {
            // Build the residual block corresponding to the track observation:
            const View * view = sfm_data.views.at(itObs->first).get();

            // Each Residual block takes a point and a camera as input and outputs a 2
            // dimensional residual. Internally, the cost function stores the observed
            // image location and compares the reprojection against the observation.
            ceres::CostFunction* cost_function =
              IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), itObs->second.x);

            if (cost_function)
              problem.AddResidualBlock(cost_function,
                p_LossFunction,
                &map_intrinsics[view->id_intrinsic][0],
                &map_poses[view->id_pose][0],
                iterTracks->second.X.data());

          }
          if (options.structure_opt == Structure_Parameter_Type::NONE)
            problem.SetParameterBlockConstant(iterTracks->second.X.data());
        }

        if (options.control_point_opt.bUse_control_points)
        {
          // Use Ground Control Point:
          // - fixed 3D points with weighted observations
          for (Landmarks::iterator iterGCPTracks = sfm_data.control_points.begin();
            iterGCPTracks!= sfm_data.control_points.end(); ++iterGCPTracks)
          {
            const Observations & obs = iterGCPTracks->second.obs;

            for (Observations::const_iterator itObs = obs.begin();
              itObs != obs.end(); ++itObs)
            {
              // Build the residual block corresponding to the track observation:
              const View * view = sfm_data.views.at(itObs->first).get();

              // Each Residual block takes a point and a camera as input and outputs a 2
              // dimensional residual. Internally, the cost function stores the observed
              // image location and compares the reprojection against the observation.
              ceres::CostFunction* cost_function =
                IntrinsicsToCostFunction(
                  sfm_data.intrinsics[view->id_intrinsic].get(),
                  itObs->second.x,
                  options.control_point_opt.weight);

              if (cost_function)
                problem.AddResidualBlock(cost_function,
                  nullptr,
                  &map_intrinsics[view->id_intrinsic][0],
                  &map_poses[view->id_pose][0],
                  iterGCPTracks->second.X.data());
            }
            // Set the 3D point as FIXED (it's a GCP)
            problem.SetParameterBlockConstant(iterGCPTracks->second.X.data());
          }
        }

        // Configure a BA engine and run it
        //  Make Ceres automatically detect the bundle structure.
        ceres::Solver::Options ceres_config_options;
        ceres_config_options.max_num_iterations = 500;
        ceres_config_options.preconditioner_type =
          static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
        ceres_config_options.linear_solver_type =
          static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
        ceres_config_options.sparse_linear_algebra_library_type =
          static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
        ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
        ceres_config_options.logging_type = ceres::SILENT;
        ceres_config_options.num_threads = ceres_options_.nb_threads_;
        ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
        ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

        //std::cout << "error 4" << std::endl;
        // Solve BA
        ceres::Solver::Summary summary;
        ceres::Solve(ceres_config_options, &problem, &summary);
        if (ceres_options_.bCeres_summary_)
          std::cout << summary.FullReport() << std::endl;

        // If no error, get back refined parameters
        if (!summary.IsSolutionUsable())
        {
          if (ceres_options_.bVerbose_)
            std::cout << "Bundle Adjustment failed." << std::endl;
          return false;
        }
        else // Solution is usable
        {
          if (ceres_options_.bVerbose_)
          {
            // Display statistics about the minimization
            std::cout << std::endl
              << "Bundle Adjustment statistics (approximated RMSE):\n"
              << " #views: " << sfm_data.views.size() << "\n"
              << " #poses: " << sfm_data.poses.size() << "\n"
              << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
              << " #tracks: " << sfm_data.structure.size() << "\n"
              << " #residuals: " << summary.num_residuals << "\n"
              << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
              << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
              << " Time (s): " << summary.total_time_in_seconds << "\n"
              << std::endl;
          }

          // Update camera poses with refined data
          if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
          {
            for (Poses::iterator itPose = sfm_data.poses.begin();
              itPose != sfm_data.poses.end(); ++itPose)
            {
              const IndexT indexPose = itPose->first;

              Mat3 R_refined;
              ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
              Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
              // Update the pose
              Pose3 & pose = itPose->second;
              pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
            }
          }

          // Update camera intrinsics with refined data
          if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
          {
            for (Intrinsics::iterator itIntrinsic = sfm_data.intrinsics.begin();
              itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
            {
              const IndexT indexCam = itIntrinsic->first;

              const std::vector<double> & vec_params = map_intrinsics[indexCam];
              itIntrinsic->second.get()->updateFromParams(vec_params);
            }
          }
          // Structure is already updated directly if needed (no data wrapping)
          return true;
        }
}

bool Bundle_Adjustment_Ceres::Adjust_down_gps_cail_first_change
(SfM_Data & sfm_data,     // the SfM scene to refine
 const Optimize_Options options,
 std::vector<double>& out_transformation_R_imu,
 std::vector<double>& out_tramsformation_x_gps
)
{
    //----------
    // Add camera parameters
    // - intrinsics
    // - poses [R|t]

    // Create residuals for each observation in the bundle adjustment problem. The
    // parameters for cameras and points are added automatically.
    //----------

    ceres::Problem problem;

    //set parameters about gps imu
    Hash_Map<IndexT, std::vector<double> > R_imu;
    Hash_Map<IndexT, std::vector<double> > C_gps;
    std::vector<double> transformation_R_imu = {0.0,0.0,0.0};//{-0.02074, 0.0199567, 0.0271603}; //{-0.0177615, 0.0177593, 0.0196154};  //{0.0,0.0,0.0}; //{-0.000427019, 0.0125482, 0.0134619};//{-0.00594065, 0.0111603, 0.0170654}; //{0.0,0.0,0.0};//{-0.00986379, 0.00976848, 0.0253123};
    std::vector<double> tramsformation_x_gps = {0.0,0.0,0.0}; //{-0.82937, 0.462316, -2.59147}; //{0.0,0.0,0.0};//{-0.281802, 0.574862, 0.064961};

    double * parameter_transformation_R_imu = &transformation_R_imu[0];
    problem.AddParameterBlock(parameter_transformation_R_imu, 3);

    double * parameter_tramsformation_x_gps = &tramsformation_x_gps[0];
    problem.AddParameterBlock(parameter_tramsformation_x_gps, 3);


    //problem.SetParameterBlockConstant(parameter_transformation_R_imu);
    //problem.SetParameterBlockConstant(parameter_tramsformation_x_gps);


    //get gps and imu data
  {
        {
            //read file and initialise external
            std::vector<double> c(3);
            std::string tFilePath = sfm_data.s_root_path + "/latlon_c.res";  //!!!
            std::ifstream in;
            in.open(tFilePath);
            if(!in)
            {
                std::cout << "open tFilePath file error!" << std::endl;
                return false;
            }
            std::size_t line = 0;
            const int count = sfm_data.views.size();
            while((!in.eof()) && (line < count))
            {
                 in >> c[0] >> c[1] >> c[2];
                 C_gps[line] = {c[0], c[1], c[2]};
                 double * parameter_block = &C_gps[line][0];
                 problem.AddParameterBlock(parameter_block, C_gps[line].size());
                 problem.SetParameterBlockConstant(parameter_block);  ///needn't change
                ++line;
            }
            in.close();

        }

        //read imu
        {
                std::ifstream in;
                in.open(sfm_data.s_root_path + "/IMU.res"); //imu file

                int line = 0;
                const int count = sfm_data.views.size();

        //        std::ofstream out_show_camera;
        //        out_show_camera.open(sfm_data_.s_root_path + "/show_imu_R.txt");

                while((!in.eof()) && (line < count))
                {
                    double roll, pitch, yaw;
                    in >> roll >> pitch >> yaw;  //z x y

                    Mat3 Rz, Ry, Rx;
                    Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
                          sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
                          0.0, 0.0, 1.0;
                    Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
                          0.0, 1.0, 0.0,
                          -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
                    Rx << 1.0, 0.0, 0.0,
                          0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
                          0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);
                    Mat3 R_imu_in =  Rz*Rx*Ry;
                    Mat3 changeAxis_M;
                    changeAxis_M << 1.0, 0.0, 0.0,
                                    0.0, -1.0, 0.0,
                                    0.0, 0.0, -1.0;
                    R_imu_in = changeAxis_M * R_imu_in;


                    Vec3 temp_rd;
                    ceres::RotationMatrixToAngleAxis(R_imu_in.data(), temp_rd.data());
                    R_imu[line] = {temp_rd(0), temp_rd(1), temp_rd(2)};

                    double * parameter_block = &R_imu[line][0];
                    problem.AddParameterBlock(parameter_block, R_imu[line].size());
                    problem.SetParameterBlockConstant(parameter_block);

                    ++line;
                }
                in.close();


            }


  }


///-----------------------------------------
///
///
        //set reference
        double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};
        Vec3 t_dd;
        t_dd << 0.0, 0.0, 0.0;


        std::vector<double> transformation_dr(6);
        transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                             t_dd(0), t_dd(1), t_dd(2)};


    std::cout << std::endl << std::endl << "initial value: " << std::endl;

        std::cout << "transformation_dr : " << transformation_dr[0] << " "
                                            << transformation_dr[1] << " "
                                            << transformation_dr[2] << " "
                                            << transformation_dr[3] << " "
                                            << transformation_dr[4] << " "
                                            << transformation_dr[5] << " " << std::endl;

        std::cout << "transformation_R_imu : " << transformation_R_imu[0] << " "
                                               << transformation_R_imu[1] << " "
                                               << transformation_R_imu[2] << std::endl;

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;


        double * parameter_transformation_dr = &transformation_dr[0];
        problem.AddParameterBlock(parameter_transformation_dr, 6);

        problem.SetParameterBlockConstant(parameter_transformation_dr);  ///needn't change


#ifndef CALIBRATION_THREE_CAMERAS
        //problem.SetParameterBlockConstant(parameter_transformation_br);
        //problem.SetParameterBlockConstant(parameter_transformation_fr);
#endif

        if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_transformation_dr);
        }
        else  // Subset parametrization
        {
          std::vector<int> vec_constant_extrinsic_dr;
          // If we adjust only the translation, we must set ROTATION as constant
          if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
          {
            // Subset rotation parametrization
            vec_constant_extrinsic_dr.push_back(0);
            vec_constant_extrinsic_dr.push_back(1);
            vec_constant_extrinsic_dr.push_back(2);
          }
          // If we adjust only the rotation, we must set TRANSLATION as constant
          if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
          {
            // Subset translation parametrization
              vec_constant_extrinsic_dr.push_back(3);
              vec_constant_extrinsic_dr.push_back(4);
              vec_constant_extrinsic_dr.push_back(5);
          }
          if (!vec_constant_extrinsic_dr.empty())
          {
            ceres::SubsetParameterization *subset_parameterization_dr =
              new ceres::SubsetParameterization(6, vec_constant_extrinsic_dr);
            problem.SetParameterization(parameter_transformation_dr, subset_parameterization_dr);
          }
        }


    // Data wrapper for refinement:
    Hash_Map<IndexT, std::vector<double> > map_intrinsics;

   // Setup Intrinsics data & subparametrization
    for (Intrinsics::const_iterator itIntrinsic = sfm_data.intrinsics.begin();
      itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
    {
      const IndexT indexCam = itIntrinsic->first;

      if (isValid(itIntrinsic->second->getType()))
      {
        map_intrinsics[indexCam] = itIntrinsic->second->getParams();

        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            itIntrinsic->second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
      else
      {
        std::cerr << "Unsupported camera type." << std::endl;
      }
    }


    // Set a LossFunction to be less penalized by false measurements
    //  - set it to NULL if you don't want use a lossFunction.
    ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(Square(4.0));
    // TODO: make the LOSS function and the parameter an option

std::cout << "start." << std::endl;

    // For all visibility add reprojections errors:
    for (Landmarks::iterator iterTracks = sfm_data.structure.begin();
         iterTracks!= sfm_data.structure.end(); ++iterTracks)
    {
      const Observations & obs = iterTracks->second.obs;

      for (Observations::const_iterator itObs = obs.begin();
        itObs != obs.end(); ++itObs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(itObs->first).get();

        //decide camera pose through pose id
        double * parameter_transformation;
        IndexT poseId = view->id_pose;


        //single
//        if (hashMap.find(poseId)->second == 0)
//        {

//            // Each Residual block takes a point and a camera as input and outputs a 2
//            // dimensional residual. Internally, the cost function stores the observed
//            // image location and compares the reprojection against the observation.
//            ceres::CostFunction* cost_function =
//              IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), itObs->second.x);

//            if (cost_function)
//              problem.AddResidualBlock(cost_function,
//                p_LossFunction,
//                &map_intrinsics[view->id_intrinsic][0],
//                &single_poses[view->id_pose][0],
//                iterTracks->second.X.data());

//        }else
        {  //
                ceres::CostFunction* cost_function2 =
                        new ceres::AutoDiffCostFunction
                          <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
                            new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

                if (cost_function2)
                  problem.AddResidualBlock(cost_function2,
                    p_LossFunction,
                    &map_intrinsics[view->id_intrinsic][0],
                    &R_imu[poseId][0],
                    &C_gps[poseId][0],
                    &transformation_R_imu[0],
                    &tramsformation_x_gps[0],
                    &transformation_dr[0],
                    iterTracks->second.X.data());  //

        }

      }
            if (options.structure_opt == Structure_Parameter_Type::NONE)
              problem.SetParameterBlockConstant(iterTracks->second.X.data());
    }



    // Configure a BA engine and run it
    //  Make Ceres automatically detect the bundle structure.
    ceres::Solver::Options ceres_config_options;
    ceres_config_options.max_num_iterations = 500;
    ceres_config_options.preconditioner_type =
      static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
    ceres_config_options.linear_solver_type =
      static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
    ceres_config_options.sparse_linear_algebra_library_type =
      static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
    ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
    ceres_config_options.logging_type = ceres::SILENT;
    ceres_config_options.num_threads = ceres_options_.nb_threads_;
    ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
    ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;


    // Solve BA
    ceres::Solver::Summary summary;
    ceres::Solve(ceres_config_options, &problem, &summary);
    if (ceres_options_.bCeres_summary_)
      std::cout << summary.FullReport() << std::endl;


    // If no error, get back refined parameters
    if (!summary.IsSolutionUsable())
    {
      if (ceres_options_.bVerbose_)
        std::cout << "Bundle Adjustment failed." << std::endl;
      return false;
    }
    else // Solution is usable
    {
      if (ceres_options_.bVerbose_)
      {
        // Display statistics about the minimization
        std::cout << std::endl
          << "Bundle Adjustment statistics (approximated RMSE) Adjust_threeCameras:\n"
          << " #views: " << sfm_data.views.size() << "\n"
          << " #poses: " << sfm_data.poses.size() << "\n"
          << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
          << " #tracks: " << sfm_data.structure.size() << "\n"
          << " #residuals: " << summary.num_residuals << "\n"
          << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
          << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
          << " Time (s): " << summary.total_time_in_seconds << "\n"
          << std::endl;


        //show reslut
        std::cout << "show reslut : " << std::endl;
        std::cout << "transformation_dr : " << transformation_dr[0] << " "
                                            << transformation_dr[1] << " "
                                            << transformation_dr[2] << " "
                                            << transformation_dr[3] << " "
                                            << transformation_dr[4] << " "
                                            << transformation_dr[5] << std::endl;

        std::cout << "transformation_R_imu : " << transformation_R_imu[0] << " "
                                               << transformation_R_imu[1] << " "
                                               << transformation_R_imu[2] << std::endl;

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;

        out_transformation_R_imu[0] = transformation_R_imu[0];
        out_transformation_R_imu[1] = transformation_R_imu[1];
        out_transformation_R_imu[2] = transformation_R_imu[2];

        out_tramsformation_x_gps[0] = tramsformation_x_gps[0];
        out_tramsformation_x_gps[1] = tramsformation_x_gps[1];
        out_tramsformation_x_gps[2] = tramsformation_x_gps[2];
      }

      // Update camera poses with refined data
      if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
      {
        for (Poses::iterator itPose = sfm_data.poses.begin();
          itPose != sfm_data.poses.end(); ++itPose)
        {
          const IndexT indexPose = itPose->first;
          IndexT camId = indexPose;


          {
              Mat3 Rr_trans;
              Vec3 tr_trans;
              {
                  ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rr_trans.data());

                  tr_trans(0) = transformation_dr[3];
                  tr_trans(1) = transformation_dr[4];
                  tr_trans(2) = transformation_dr[5];

              }

              Vec3 gps_trans;
              gps_trans(0) = tramsformation_x_gps[0];
              gps_trans(1) = tramsformation_x_gps[1];
              gps_trans(2) = tramsformation_x_gps[2];
              Vec3 C_gps_res;
              C_gps_res(0) = C_gps[camId][0];
              C_gps_res(1) = C_gps[camId][1];
              C_gps_res(2) = C_gps[camId][2];
              Mat3 R_imu_trans;
              ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
              Mat3 R_imu_mat;
              ceres::AngleAxisToRotationMatrix(&R_imu[camId][0], R_imu_mat.data());

              Mat3 R_res;
              R_res = Rr_trans * R_imu_trans * R_imu_mat;
              Vec3 t_res;
              t_res = Rr_trans * gps_trans - Rr_trans * R_imu_trans * R_imu_mat * C_gps_res + tr_trans;

              Pose3 & pose = itPose->second;
              pose = Pose3(R_res, -R_res.transpose() * t_res);
            }
          }

      }


      // Update camera intrinsics with refined data
      if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
      {
        for (Intrinsics::iterator itIntrinsic = sfm_data.intrinsics.begin();
          itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
        {
          const IndexT indexCam = itIntrinsic->first;

          const std::vector<double> & vec_params = map_intrinsics[indexCam];
          itIntrinsic->second.get()->updateFromParams(vec_params);
        }
      }
      // Structure is already updated directly if needed (no data wrapping)
      return true;
    }



    }


bool Bundle_Adjustment_Ceres::Adjust_down_gps_cail_second_change
(SfM_Data & sfm_data,     // the SfM scene to refine
 const Optimize_Options options,
 std::vector<double>& in_transformation_R_imu,
 std::vector<double>& in_tramsformation_x_gps
)
{
    //----------
    // Add camera parameters
    // - intrinsics
    // - poses [R|t]

    // Create residuals for each observation in the bundle adjustment problem. The
    // parameters for cameras and points are added automatically.
    //----------

    ceres::Problem problem;

    //set parameters about gps imu
    Hash_Map<IndexT, std::vector<double> > R_imu;
    Hash_Map<IndexT, std::vector<double> > C_gps;
    std::vector<double> transformation_R_imu = in_transformation_R_imu;//{0.0,0.0,0.0};//{-0.0120417, 0.0201986, 0.0258046};//{0.0,0.0,0.0}; //{-0.02074, 0.0199567, 0.0271603}; //{-0.0177615, 0.0177593, 0.0196154};  //{0.0,0.0,0.0}; //{-0.000427019, 0.0125482, 0.0134619};//{-0.00594065, 0.0111603, 0.0170654}; //{0.0,0.0,0.0};//{-0.00986379, 0.00976848, 0.0253123};
    std::vector<double> tramsformation_x_gps = in_tramsformation_x_gps;//{0.0,0.0,0.0};//{-2.01425, 2.23307, -28.3916};


    double * parameter_transformation_R_imu = &transformation_R_imu[0];
    problem.AddParameterBlock(parameter_transformation_R_imu, 3);

    double * parameter_tramsformation_x_gps = &tramsformation_x_gps[0];
    problem.AddParameterBlock(parameter_tramsformation_x_gps, 3);


    //problem.SetParameterBlockConstant(parameter_transformation_R_imu);
    //problem.SetParameterBlockConstant(parameter_tramsformation_x_gps);


    //get gps and imu data
  {
        {
            //read file and initialise external
            std::vector<double> c(3);
            std::string tFilePath = sfm_data.s_root_path + "/latlon_c.res";  //!!!
            std::ifstream in;
            in.open(tFilePath);
            if(!in)
            {
                std::cout << "open tFilePath file error!" << std::endl;
                return false;
            }
            std::size_t line = 0;
            const int count = sfm_data.views.size();
            while((!in.eof()) && (line < count))
            {
                 in >> c[0] >> c[1] >> c[2];
                 C_gps[line] = {c[0], c[1], c[2]};
                 double * parameter_block = &C_gps[line][0];
                 problem.AddParameterBlock(parameter_block, C_gps[line].size());
                 problem.SetParameterBlockConstant(parameter_block);  ///needn't change
                ++line;
            }
            in.close();

        }

        //read imu
        {
                std::ifstream in;
                in.open(sfm_data.s_root_path + "/IMU.res"); //imu file

                int line = 0;
                const int count = sfm_data.views.size();

        //        std::ofstream out_show_camera;
        //        out_show_camera.open(sfm_data_.s_root_path + "/show_imu_R.txt");

                while((!in.eof()) && (line < count))
                {
                    double roll, pitch, yaw;
                    in >> roll >> pitch >> yaw;  //z x y

                    Mat3 Rz, Ry, Rx;
                    Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
                          sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
                          0.0, 0.0, 1.0;
                    Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
                          0.0, 1.0, 0.0,
                          -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
                    Rx << 1.0, 0.0, 0.0,
                          0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
                          0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);
                    Mat3 R_imu_in =  Rz*Rx*Ry;
                    Mat3 changeAxis_M;
                    changeAxis_M << 1.0, 0.0, 0.0,
                                    0.0, -1.0, 0.0,
                                    0.0, 0.0, -1.0;
                    R_imu_in = changeAxis_M * R_imu_in;


                    Vec3 temp_rd;
                    ceres::RotationMatrixToAngleAxis(R_imu_in.data(), temp_rd.data());
                    R_imu[line] = {temp_rd(0), temp_rd(1), temp_rd(2)};

                    double * parameter_block = &R_imu[line][0];
                    problem.AddParameterBlock(parameter_block, R_imu[line].size());

                    ++line;
                }
                in.close();


            }



  }


//    Hash_Map<IndexT, std::vector<double> > single_poses;
//    for (Poses::const_iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose)
//    {
//      const IndexT indexPose = itPose->first;

//      double * parameter_block;

//      {
//          //using in second BA
//          const Pose3 & pose = itPose->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 t = pose.translation();

//          Vec3 C = pose.center();
//          Vec3 temp_rd;
//          ceres::RotationMatrixToAngleAxis(R.data(), temp_rd.data());
//          R_imu[indexPose] = {temp_rd(0), temp_rd(1), temp_rd(2)};

//          double * parameter_block = &R_imu[indexPose][0];
//          problem.AddParameterBlock(parameter_block, R_imu[indexPose].size());

//      }

//      if(parameter_block == NULL)
//      {
//          continue;
//      }


//    }

///-----------------------------------------
///
///

        //set reference
        double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};
        Vec3 t_dd;
        t_dd << 0.0, 0.0, 0.0;



        std::vector<double> transformation_dr(6);

        transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                             t_dd(0), t_dd(1), t_dd(2)};


    std::cout << std::endl << std::endl << "initial value: " << std::endl;

        std::cout << "transformation_dr : " << transformation_dr[0] << " "
                                            << transformation_dr[1] << " "
                                            << transformation_dr[2] << " "
                                            << transformation_dr[3] << " "
                                            << transformation_dr[4] << " "
                                            << transformation_dr[5] << " " << std::endl;

        std::cout << "transformation_R_imu : " << transformation_R_imu[0] << " "
                                               << transformation_R_imu[1] << " "
                                               << transformation_R_imu[2] << std::endl;

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;


        double * parameter_transformation_dr = &transformation_dr[0];
        problem.AddParameterBlock(parameter_transformation_dr, 6);

        problem.SetParameterBlockConstant(parameter_transformation_dr);  ///needn't change


        if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_transformation_dr);
        }
        else  // Subset parametrization
        {
          std::vector<int> vec_constant_extrinsic_br, vec_constant_extrinsic_fr,vec_constant_extrinsic_dr;
          // If we adjust only the translation, we must set ROTATION as constant
          if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
          {
            // Subset rotation parametrization
            vec_constant_extrinsic_dr.push_back(0);
            vec_constant_extrinsic_dr.push_back(1);
            vec_constant_extrinsic_dr.push_back(2);
          }
          // If we adjust only the rotation, we must set TRANSLATION as constant
          if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
          {
            // Subset translation parametrization
              vec_constant_extrinsic_dr.push_back(3);
              vec_constant_extrinsic_dr.push_back(4);
              vec_constant_extrinsic_dr.push_back(5);
          }
          if (!vec_constant_extrinsic_br.empty())
          {
            ceres::SubsetParameterization *subset_parameterization_dr =
              new ceres::SubsetParameterization(6, vec_constant_extrinsic_dr);
            problem.SetParameterization(parameter_transformation_dr, subset_parameterization_dr);
          }
        }


    // Data wrapper for refinement:
    Hash_Map<IndexT, std::vector<double> > map_intrinsics;

   // Setup Intrinsics data & subparametrization
    for (Intrinsics::const_iterator itIntrinsic = sfm_data.intrinsics.begin();
      itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
    {
      const IndexT indexCam = itIntrinsic->first;

      if (isValid(itIntrinsic->second->getType()))
      {
        map_intrinsics[indexCam] = itIntrinsic->second->getParams();

        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            itIntrinsic->second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
      else
      {
        std::cerr << "Unsupported camera type." << std::endl;
      }
    }


    // Set a LossFunction to be less penalized by false measurements
    //  - set it to NULL if you don't want use a lossFunction.
    ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(Square(4.0));
    // TODO: make the LOSS function and the parameter an option

std::cout << "start." << std::endl;

    // For all visibility add reprojections errors:
    for (Landmarks::iterator iterTracks = sfm_data.structure.begin();
         iterTracks!= sfm_data.structure.end(); ++iterTracks)
    {
      const Observations & obs = iterTracks->second.obs;

      for (Observations::const_iterator itObs = obs.begin();
        itObs != obs.end(); ++itObs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(itObs->first).get();

        //decide camera pose through pose id
        double * parameter_transformation;
        IndexT poseId = view->id_pose;



        {  //

            {  //down

//                ceres::CostFunction* cost_function2 =
//                        new ceres::AutoDiffCostFunction
//                          <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
//                            new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

//                if (cost_function2)
//                  problem.AddResidualBlock(cost_function2,
//                    p_LossFunction,
//                    &map_intrinsics[view->id_intrinsic][0],
//                    &R_imu[poseId][0],
//                    &C_gps[poseId][0],
//                    &transformation_R_imu[0],
//                    &tramsformation_x_gps[0],
//                    &transformation_dr[0],
//                    iterTracks->second.X.data());  //

                ceres::CostFunction* cost_function2 =
                        new ceres::AutoDiffCostFunction
                          <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
                            new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

                if (cost_function2)
                  problem.AddResidualBlock(cost_function2,
                    p_LossFunction,
                    &map_intrinsics[view->id_intrinsic][0],
                    &R_imu[poseId][0],
                    &C_gps[poseId][0],
                    &transformation_R_imu[0],
                    &tramsformation_x_gps[0],
                    &transformation_dr[0],
                    iterTracks->second.X.data());  //

            }

        }

      }
            if (options.structure_opt == Structure_Parameter_Type::NONE)
              problem.SetParameterBlockConstant(iterTracks->second.X.data());
    }



    // Configure a BA engine and run it
    //  Make Ceres automatically detect the bundle structure.
    ceres::Solver::Options ceres_config_options;
    ceres_config_options.max_num_iterations = 500;
    ceres_config_options.preconditioner_type =
      static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
    ceres_config_options.linear_solver_type =
      static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
    ceres_config_options.sparse_linear_algebra_library_type =
      static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
    ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
    ceres_config_options.logging_type = ceres::SILENT;
    ceres_config_options.num_threads = ceres_options_.nb_threads_;
    ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
    ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;


    // Solve BA
    ceres::Solver::Summary summary;
    ceres::Solve(ceres_config_options, &problem, &summary);
    if (ceres_options_.bCeres_summary_)
      std::cout << summary.FullReport() << std::endl;


    // If no error, get back refined parameters
    if (!summary.IsSolutionUsable())
    {
      if (ceres_options_.bVerbose_)
        std::cout << "Bundle Adjustment failed." << std::endl;
      return false;
    }
    else // Solution is usable
    {
      if (ceres_options_.bVerbose_)
      {
        // Display statistics about the minimization
        std::cout << std::endl
          << "Bundle Adjustment statistics (approximated RMSE) Adjust_threeCameras:\n"
          << " #views: " << sfm_data.views.size() << "\n"
          << " #poses: " << sfm_data.poses.size() << "\n"
          << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
          << " #tracks: " << sfm_data.structure.size() << "\n"
          << " #residuals: " << summary.num_residuals << "\n"
          << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
          << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
          << " Time (s): " << summary.total_time_in_seconds << "\n"
          << std::endl;


        //show reslut
        std::cout << "show reslut : " << std::endl;
        std::cout << "transformation_dr : " << transformation_dr[0] << " "
                                            << transformation_dr[1] << " "
                                            << transformation_dr[2] << " "
                                            << transformation_dr[3] << " "
                                            << transformation_dr[4] << " "
                                            << transformation_dr[5] << std::endl;

        std::cout << "transformation_R_imu : " << transformation_R_imu[0] << " "
                                               << transformation_R_imu[1] << " "
                                               << transformation_R_imu[2] << std::endl;

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;

        in_transformation_R_imu[0] = transformation_R_imu[0];
        in_transformation_R_imu[1] = transformation_R_imu[1];
        in_transformation_R_imu[2] = transformation_R_imu[2];

        in_tramsformation_x_gps[0] = tramsformation_x_gps[0];
        in_tramsformation_x_gps[1] = tramsformation_x_gps[1];
        in_tramsformation_x_gps[2] = tramsformation_x_gps[2];
      }

      // Update camera poses with refined data
      if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
      {
        for (Poses::iterator itPose = sfm_data.poses.begin();
          itPose != sfm_data.poses.end(); ++itPose)
        {
          const IndexT indexPose = itPose->first;
          IndexT camId = indexPose;

          //single
//          if (hashMap.find(camId)->second == 0)
//          {
//              Mat3 R_refined;
//              ceres::AngleAxisToRotationMatrix(&single_poses[indexPose][0], R_refined.data());
//              Vec3 t_refined(single_poses[indexPose][3], single_poses[indexPose][4], single_poses[indexPose][5]);
//              // Update the pose
//              Pose3 & pose = itPose->second;
//              pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
//          }else
          {  //unity
              Mat3 Rr_trans;
              Vec3 tr_trans;

              {
                  ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rr_trans.data());

                  tr_trans(0) = transformation_dr[3];
                  tr_trans(1) = transformation_dr[4];
                  tr_trans(2) = transformation_dr[5];

              }

              Vec3 gps_trans;
              gps_trans(0) = tramsformation_x_gps[0];
              gps_trans(1) = tramsformation_x_gps[1];
              gps_trans(2) = tramsformation_x_gps[2];
              Vec3 C_gps_res;
              C_gps_res(0) = C_gps[camId][0];
              C_gps_res(1) = C_gps[camId][1];
              C_gps_res(2) = C_gps[camId][2];
              Mat3 R_imu_trans;
              ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
              Mat3 R_imu_mat;
              ceres::AngleAxisToRotationMatrix(&R_imu[camId][0], R_imu_mat.data());

              Mat3 R_res;
              R_res = Rr_trans * R_imu_trans * R_imu_mat;
              Vec3 t_res;
              t_res = Rr_trans * gps_trans - Rr_trans * R_imu_trans * R_imu_mat * C_gps_res + tr_trans;


              Pose3 & pose = itPose->second;
              pose = Pose3(R_res, -R_res.transpose() * t_res);
            }
          }

      }


      // Update camera intrinsics with refined data
      if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
      {
        for (Intrinsics::iterator itIntrinsic = sfm_data.intrinsics.begin();
          itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
        {
          const IndexT indexCam = itIntrinsic->first;

          const std::vector<double> & vec_params = map_intrinsics[indexCam];
          itIntrinsic->second.get()->updateFromParams(vec_params);
        }
      }
      // Structure is already updated directly if needed (no data wrapping)
      return true;
    }



    }


bool Bundle_Adjustment_Ceres::Adjust_threecamera_gps_cail_first_change
(SfM_Data & sfm_data,     // the SfM scene to refine
 const Optimize_Options options,
 std::vector<double>& out_transformation_R_imu,
 std::vector<double>& out_tramsformation_x_gps,
 std::vector<double>& out_transformation_br,
 std::vector<double>& out_transformation_fr
)
{
    //----------
    // Add camera parameters
    // - intrinsics
    // - poses [R|t]

    // Create residuals for each observation in the bundle adjustment problem. The
    // parameters for cameras and points are added automatically.
    //----------

    ceres::Problem problem;

    //set parameters about gps imu
    Hash_Map<IndexT, std::vector<double> > R_imu;
    Hash_Map<IndexT, std::vector<double> > C_gps;
    std::vector<double> transformation_R_imu = {0.0,0.0,0.0};//{-0.02074, 0.0199567, 0.0271603}; //{-0.0177615, 0.0177593, 0.0196154};  //{0.0,0.0,0.0}; //{-0.000427019, 0.0125482, 0.0134619};//{-0.00594065, 0.0111603, 0.0170654}; //{0.0,0.0,0.0};//{-0.00986379, 0.00976848, 0.0253123};
    std::vector<double> tramsformation_x_gps = {0.0,0.0,0.0}; //{-0.82937, 0.462316, -2.59147}; //{0.0,0.0,0.0};//{-0.281802, 0.574862, 0.064961};

    double * parameter_transformation_R_imu = &transformation_R_imu[0];
    problem.AddParameterBlock(parameter_transformation_R_imu, 3);

    double * parameter_tramsformation_x_gps = &tramsformation_x_gps[0];
    problem.AddParameterBlock(parameter_tramsformation_x_gps, 3);


    //problem.SetParameterBlockConstant(parameter_transformation_R_imu);
    //problem.SetParameterBlockConstant(parameter_tramsformation_x_gps);


    //get gps and imu data
  {
        {
            //read file and initialise external
            std::vector<double> c(3);
            std::string tFilePath = sfm_data.s_root_path + "/latlon_c.res";  //!!!
            std::ifstream in;
            in.open(tFilePath);
            if(!in)
            {
                std::cout << "open tFilePath file error!" << std::endl;
                return false;
            }
            std::size_t line = 0;
            const int count = sfm_data.views.size()/3;
            while((!in.eof()) && (line < count))
            {
                 in >> c[0] >> c[1] >> c[2];
                 C_gps[line + count] = {c[0], c[1], c[2]};
                 double * parameter_block = &C_gps[line + count][0];
                 problem.AddParameterBlock(parameter_block, C_gps[line + count].size());
                 problem.SetParameterBlockConstant(parameter_block);  ///needn't change
                ++line;
            }
            in.close();

        }
        //read imu
        {
        std::ifstream in;
        in.open(sfm_data.s_root_path + "/IMU.res"); //imu file

        int line = 0;
        const int count = sfm_data.views.size()/3;

        std::ofstream out_show_camera;
        out_show_camera.open(sfm_data.s_root_path + "/show_imu_R.txt");

        while((!in.eof()) && (line < count))
        {
            double roll, pitch, yaw;
            in >> roll >> pitch >> yaw;  //z x y

            Mat3 Rz, Ry, Rx, Rx_, Ry_;
            Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
                  sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
                  0.0, 0.0, 1.0;
            Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
                  0.0, 1.0, 0.0,
                  -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
            Rx << 1.0, 0.0, 0.0,
                  0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
                  0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);

            Mat3 R_d =  Rz*Rx*Ry;

            {

//                Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
//                Vec3 c = {C_gps[line + count][0], C_gps[line + count][1], C_gps[line + count][2]};
//                axis_x = R_d.transpose() * (axis_x + R_d*c);
//                axis_y = R_d.transpose() * (axis_y + R_d*c);
//                axis_z = R_d.transpose() * (axis_z + R_d*c);
//                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                       << axis_x(0) << " " << axis_x(1) << " " << axis_x(2)<< ";" << std::endl;
//                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                       << axis_y(0) << " " << axis_y(1) << " " << axis_y(2)<< ";" << std::endl;
//                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                       << axis_z(0) << " " << axis_z(1) << " " << axis_z(2)<< ";" << std::endl;

//                out_show_camera << R_d << std::endl << std::endl;
            }

//            //change axis
//            Vec3 changeAxis_A;
//            changeAxis_A << M_PI*1/sqrt(2.0), M_PI*1/sqrt(2.0), 0.0;
            Mat3 changeAxis_M;
//            ceres::AngleAxisToRotationMatrix(changeAxis_A.data(), changeAxis_M.data());

            changeAxis_M << 1.0, 0.0, 0.0,
                            0.0, -1.0, 0.0,
                            0.0, 0.0, -1.0;
            R_d = changeAxis_M * R_d;


            {

                Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
                Vec3 c = {C_gps[line + count][0], C_gps[line + count][1], C_gps[line + count][2]};
                axis_x = R_d.transpose() * (axis_x + R_d*c);
                axis_y = R_d.transpose() * (axis_y + R_d*c);
                axis_z = R_d.transpose() * (axis_z + R_d*c);
                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
                       << axis_x(0) << " " << axis_x(1) << " " << axis_x(2)<< ";" << std::endl;
                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
                       << axis_y(0) << " " << axis_y(1) << " " << axis_y(2)<< ";" << std::endl;
                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
                       << axis_z(0) << " " << axis_z(1) << " " << axis_z(2)<< ";" << std::endl;

                out_show_camera << R_d << std::endl << std::endl;
            }





            Vec3 temp_rd;
            ceres::RotationMatrixToAngleAxis(R_d.data(), temp_rd.data());
            R_imu[line + count] = {temp_rd(0), temp_rd(1), temp_rd(2)};

            double * parameter_block = &R_imu[line + count][0];
            problem.AddParameterBlock(parameter_block, R_imu[line + count].size());

            ++line;
        }
        in.close();

        out_show_camera.close();
    }
  }


    //get hash map
    std::map<int, int> hashMap;
    filiter_ununited_pose(sfm_data.poses, sfm_data.views, hashMap);
    //Hash_Map<IndexT, std::vector<double> > reference_poses;
    Hash_Map<IndexT, std::vector<double> > single_poses;
    for (Poses::const_iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose)
    {
      const IndexT indexPose = itPose->first;

      double * parameter_block;
      if(hashMap.find(indexPose)->second == 0)
      {
          std::cout << "signal " << indexPose << std::endl;
          //set single poses
          const Pose3 & pose = itPose->second;
          const Mat3 R = pose.rotation();
          const Vec3 t = pose.translation();

          double angleAxis[3];
          ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
          // angleAxis + translation
          single_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

          parameter_block = &single_poses[indexPose][0];
          problem.AddParameterBlock(parameter_block, single_poses[indexPose].size());
      }else{
          //set unity pose
//          if( (indexPose < sfm_data.views.size()/3) || (indexPose >= sfm_data.views.size()/3*2) )
//          {
//              continue;
//          }

//          const Pose3 & pose = itPose->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 t = pose.translation();
//          Vec3 C = pose.center();

//          double angleAxis[3];
//          ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
//          // angleAxis + translation
//          reference_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(0), C(1), C(2)};

//          parameter_block = &reference_poses[indexPose][0];
//          problem.AddParameterBlock(parameter_block, reference_poses[indexPose].size());


          //test if imu right (!!!)
//          if( (indexPose < sfm_data.views.size()/3) || (indexPose >= sfm_data.views.size()/3*2) )
//          {
//              continue;
//          }

//          const Pose3 & pose = itPose->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 t = pose.translation();

//          Vec3 C = pose.center();
//          Vec3 temp_rd;
//          ceres::RotationMatrixToAngleAxis(R.data(), temp_rd.data());
//          R_imu[indexPose] = {temp_rd(0), temp_rd(1), temp_rd(2)};

//          double * parameter_block = &R_imu[indexPose][0];
//          problem.AddParameterBlock(parameter_block, R_imu[line + count].size());
//          problem.SetParameterBlockConstant(parameter_block);  ///needn't change

      }


      if(parameter_block == NULL)
      {
          continue;
      }

//      if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
//      {
//        // set the whole parameter block as constant for best performance
//        problem.SetParameterBlockConstant(parameter_block);
//      }
//      else  // Subset parametrization
//      {
//        std::vector<int> vec_constant_extrinsic;
//        // If we adjust only the translation, we must set ROTATION as constant
//        if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
//        {
//          // Subset rotation parametrization
//          vec_constant_extrinsic.push_back(0);
//          vec_constant_extrinsic.push_back(1);
//          vec_constant_extrinsic.push_back(2);
//        }
//        // If we adjust only the rotation, we must set TRANSLATION as constant
//        if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
//        {
//          // Subset translation parametrization
//          vec_constant_extrinsic.push_back(3);
//          vec_constant_extrinsic.push_back(4);
//          vec_constant_extrinsic.push_back(5);
//        }
//        if (!vec_constant_extrinsic.empty())
//        {
//          ceres::SubsetParameterization *subset_parameterization =
//            new ceres::SubsetParameterization(6, vec_constant_extrinsic);
//          problem.SetParameterization(parameter_block, subset_parameterization);
//        }
//      }

    }

///-----------------------------------------
///
///

    //test distortion
    //test distortion
    Mat3 r_tb, r_tf;
    r_tb << 1.0, 0.0, 0.0,
            0.0, 0.787926, -0.61577,
            0.0, 0.61577, 0.787926;
    r_tf << 1.0, 0.0, 0.0,
            0.0, 0.78796,  0.615727,
            0.0, -0.615727, 0.78796;
    Vec3 RA_tb, RA_tf;
    ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
    ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

        //set reference
        double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), -RA_tb(2)};//{0.665224, 0.00448363, 0.00291885};  //{0.703355, 0.00263495, -0.00149232};//{0.0, 0.0, 0.0}; //{0.0212232, 0.0206259, 0.0420509};  //{0.0, 0.0, 0.0}; //{0.669661, 0.00354668, 0.00220508}; //  //{0.0, 0.0, 0.0}; //{0.662013, -0.00281273, 0.000638113 }; //{0.0, 0.0, 0.0}; //{0.671924, 0.00432461, 0.000528427};
        double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), -RA_tf(2)};//{-0.644463, 0.0113345, -0.00156344}; //{-0.682501, 0.0147873, 0.00150509};//{0.0, 0.0, 0.0}; //{-0.0232707, -0.00288052, 0.00035795};  //{0.0, 0.0, 0.0}; //{-0.650191, 0.0101123, -0.000844533}; //  //{0.0, 0.0, 0.0}; //{-0.658446, 0.00532591, 0.0010695}; //{-0.656461, 0.00508743, -0.00299628};
        double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};

        Vec3 c_bd, c_fd;//b-d, f-d
        c_bd << 0.0, 0.098, 0.027;
        c_fd << 0.0, -0.101, 0.028;

        Vec3 t_db, t_df, t_dd;
        t_db = -r_tb*c_bd;
        t_df = -r_tf*c_fd;


        //double S = 1470.531783220221;
//        t_db << 0.0, 0.098, 0.027;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
//        t_df << 0.0, -0.101, 0.028;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
        t_dd << 0.0, 0.0, 0.0;


        std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
        transformation_br = {R_br_angleAxis[0], R_br_angleAxis[1], R_br_angleAxis[2],
                             t_db(0), t_db(1), t_db(2)};
        transformation_fr = {R_fr_angleAxis[0], R_fr_angleAxis[1], R_fr_angleAxis[2],
                             t_df(0), t_df(1), t_df(2)};
        transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                             t_dd(0), t_dd(1), t_dd(2)};


    std::cout << std::endl << std::endl << "initial value: " << std::endl;
        std::cout << "transformation_br : " << transformation_br[0] << " "
                                            << transformation_br[1] << " "
                                            << transformation_br[2] << " "
                                            << transformation_br[3] << " "
                                            << transformation_br[4] << " "
                                            << transformation_br[5] << " " << std::endl;

        std::cout << "transformation_fr : " << transformation_fr[0] << " "
                                            << transformation_fr[1] << " "
                                            << transformation_fr[2] << " "
                                            << transformation_fr[3] << " "
                                            << transformation_fr[4] << " "
                                            << transformation_fr[5] << " " << std::endl;

        std::cout << "transformation_dr : " << transformation_dr[0] << " "
                                            << transformation_dr[1] << " "
                                            << transformation_dr[2] << " "
                                            << transformation_dr[3] << " "
                                            << transformation_dr[4] << " "
                                            << transformation_dr[5] << " " << std::endl;

        std::cout << "transformation_R_imu : " << transformation_R_imu[0] << " "
                                               << transformation_R_imu[1] << " "
                                               << transformation_R_imu[2] << std::endl;

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;


        double * parameter_transformation_br = &transformation_br[0];
        double * parameter_transformation_fr = &transformation_fr[0];
        double * parameter_transformation_dr = &transformation_dr[0];
        problem.AddParameterBlock(parameter_transformation_br, 6);
        problem.AddParameterBlock(parameter_transformation_fr, 6);
        problem.AddParameterBlock(parameter_transformation_dr, 6);

        problem.SetParameterBlockConstant(parameter_transformation_dr);  ///needn't change

        //problem.SetParameterBlockConstant(parameter_transformation_br);
        //problem.SetParameterBlockConstant(parameter_transformation_fr);


#ifndef CALIBRATION_THREE_CAMERAS
        //problem.SetParameterBlockConstant(parameter_transformation_br);
        //problem.SetParameterBlockConstant(parameter_transformation_fr);
#endif

        if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_transformation_br);
          problem.SetParameterBlockConstant(parameter_transformation_fr);
          problem.SetParameterBlockConstant(parameter_transformation_dr);
        }
        else  // Subset parametrization
        {
          std::vector<int> vec_constant_extrinsic_br, vec_constant_extrinsic_fr,vec_constant_extrinsic_dr;
          // If we adjust only the translation, we must set ROTATION as constant
          if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
          {
            // Subset rotation parametrization
            vec_constant_extrinsic_br.push_back(0);
            vec_constant_extrinsic_br.push_back(1);
            vec_constant_extrinsic_br.push_back(2);
            vec_constant_extrinsic_fr.push_back(0);
            vec_constant_extrinsic_fr.push_back(1);
            vec_constant_extrinsic_fr.push_back(2);
            vec_constant_extrinsic_dr.push_back(0);
            vec_constant_extrinsic_dr.push_back(1);
            vec_constant_extrinsic_dr.push_back(2);
          }
          // If we adjust only the rotation, we must set TRANSLATION as constant
          if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
          {
            // Subset translation parametrization
              vec_constant_extrinsic_br.push_back(3);
              vec_constant_extrinsic_br.push_back(4);
              vec_constant_extrinsic_br.push_back(5);
              vec_constant_extrinsic_fr.push_back(3);
              vec_constant_extrinsic_fr.push_back(4);
              vec_constant_extrinsic_fr.push_back(5);
              vec_constant_extrinsic_dr.push_back(3);
              vec_constant_extrinsic_dr.push_back(4);
              vec_constant_extrinsic_dr.push_back(5);
          }
          if (!vec_constant_extrinsic_br.empty())
          {
            ceres::SubsetParameterization *subset_parameterization_br =
              new ceres::SubsetParameterization(6, vec_constant_extrinsic_br);
            problem.SetParameterization(parameter_transformation_br, subset_parameterization_br);

            ceres::SubsetParameterization *subset_parameterization_fr =
              new ceres::SubsetParameterization(6, vec_constant_extrinsic_fr);
            problem.SetParameterization(parameter_transformation_fr, subset_parameterization_fr);

            ceres::SubsetParameterization *subset_parameterization_dr =
              new ceres::SubsetParameterization(6, vec_constant_extrinsic_dr);
            problem.SetParameterization(parameter_transformation_dr, subset_parameterization_dr);
          }
        }


    // Data wrapper for refinement:
    Hash_Map<IndexT, std::vector<double> > map_intrinsics;

   // Setup Intrinsics data & subparametrization
    for (Intrinsics::const_iterator itIntrinsic = sfm_data.intrinsics.begin();
      itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
    {
      const IndexT indexCam = itIntrinsic->first;

      if (isValid(itIntrinsic->second->getType()))
      {
        map_intrinsics[indexCam] = itIntrinsic->second->getParams();

        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            itIntrinsic->second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
      else
      {
        std::cerr << "Unsupported camera type." << std::endl;
      }
    }


    // Set a LossFunction to be less penalized by false measurements
    //  - set it to NULL if you don't want use a lossFunction.
    ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(Square(4.0));
    // TODO: make the LOSS function and the parameter an option

std::cout << "start." << std::endl;

    // For all visibility add reprojections errors:
    for (Landmarks::iterator iterTracks = sfm_data.structure.begin();
         iterTracks!= sfm_data.structure.end(); ++iterTracks)
    {
      const Observations & obs = iterTracks->second.obs;

      for (Observations::const_iterator itObs = obs.begin();
        itObs != obs.end(); ++itObs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(itObs->first).get();

        //decide camera pose through pose id
        double * parameter_transformation;
        IndexT poseId = view->id_pose;


        //single
        if (hashMap.find(poseId)->second == 0)
        {

            // Each Residual block takes a point and a camera as input and outputs a 2
            // dimensional residual. Internally, the cost function stores the observed
            // image location and compares the reprojection against the observation.
            ceres::CostFunction* cost_function =
              IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), itObs->second.x);

            if (cost_function)
              problem.AddResidualBlock(cost_function,
                p_LossFunction,
                &map_intrinsics[view->id_intrinsic][0],
                &single_poses[view->id_pose][0],
                iterTracks->second.X.data());

        }else{  //unity
            if(poseId < sfm_data.views.size()/3)  //back
            {
                poseId = poseId + sfm_data.views.size()/3;

//                  ceres::CostFunction* cost_function2 =
//                          new ceres::AutoDiffCostFunction
//                            <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
//                              new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

//                  if (cost_function2)
//                    problem.AddResidualBlock(cost_function2,
//                      p_LossFunction,
//                      &map_intrinsics[view->id_intrinsic][0],
//                      &R_imu[poseId][0],
//                      &C_gps[poseId][0],
//                      &transformation_R_imu[0],
//                      &tramsformation_x_gps[0],
//                      &transformation_br[0],
//                      iterTracks->second.X.data());  //


                  ceres::CostFunction* cost_function2 =
                          new ceres::AutoDiffCostFunction
                            <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
                              new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

                  if (cost_function2)
                    problem.AddResidualBlock(cost_function2,
                      p_LossFunction,
                      &map_intrinsics[view->id_intrinsic][0],
                      &R_imu[poseId][0],
                      &C_gps[poseId][0],
                      &transformation_R_imu[0],
                      &tramsformation_x_gps[0],
                      &transformation_br[0],
                      iterTracks->second.X.data());  //




            }else if(poseId >= sfm_data.views.size()/3*2)  //front
            {
                poseId = poseId - sfm_data.views.size()/3;

//                ceres::CostFunction* cost_function2 =
//                        new ceres::AutoDiffCostFunction
//                          <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
//                            new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

//                if (cost_function2)
//                  problem.AddResidualBlock(cost_function2,
//                    p_LossFunction,
//                    &map_intrinsics[view->id_intrinsic][0],
//                    &R_imu[poseId][0],
//                    &C_gps[poseId][0],
//                    &transformation_R_imu[0],
//                    &tramsformation_x_gps[0],
//                    &transformation_fr[0],
//                    iterTracks->second.X.data());  //


                ceres::CostFunction* cost_function2 =
                        new ceres::AutoDiffCostFunction
                          <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
                            new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

                if (cost_function2)
                  problem.AddResidualBlock(cost_function2,
                    p_LossFunction,
                    &map_intrinsics[view->id_intrinsic][0],
                    &R_imu[poseId][0],
                    &C_gps[poseId][0],
                    &transformation_R_imu[0],
                    &tramsformation_x_gps[0],
                    &transformation_fr[0],
                    iterTracks->second.X.data());  //


            }else{  //down

//                ceres::CostFunction* cost_function2 =
//                        new ceres::AutoDiffCostFunction
//                          <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
//                            new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

//                if (cost_function2)
//                  problem.AddResidualBlock(cost_function2,
//                    p_LossFunction,
//                    &map_intrinsics[view->id_intrinsic][0],
//                    &R_imu[poseId][0],
//                    &C_gps[poseId][0],
//                    &transformation_R_imu[0],
//                    &tramsformation_x_gps[0],
//                    &transformation_dr[0],
//                    iterTracks->second.X.data());  //



                ceres::CostFunction* cost_function2 =
                        new ceres::AutoDiffCostFunction
                          <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
                            new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

                if (cost_function2)
                  problem.AddResidualBlock(cost_function2,
                    p_LossFunction,
                    &map_intrinsics[view->id_intrinsic][0],
                    &R_imu[poseId][0],
                    &C_gps[poseId][0],
                    &transformation_R_imu[0],
                    &tramsformation_x_gps[0],
                    &transformation_dr[0],
                    iterTracks->second.X.data());  //


            }

        }

      }
            if (options.structure_opt == Structure_Parameter_Type::NONE)
              problem.SetParameterBlockConstant(iterTracks->second.X.data());
    }



    // Configure a BA engine and run it
    //  Make Ceres automatically detect the bundle structure.
    ceres::Solver::Options ceres_config_options;
    ceres_config_options.max_num_iterations = 500;
    ceres_config_options.preconditioner_type =
      static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
    ceres_config_options.linear_solver_type =
      static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
    ceres_config_options.sparse_linear_algebra_library_type =
      static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
    ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
    ceres_config_options.logging_type = ceres::SILENT;
    ceres_config_options.num_threads = ceres_options_.nb_threads_;
    ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
    ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;


    // Solve BA
    ceres::Solver::Summary summary;
    ceres::Solve(ceres_config_options, &problem, &summary);
    if (ceres_options_.bCeres_summary_)
      std::cout << summary.FullReport() << std::endl;


    // If no error, get back refined parameters
    if (!summary.IsSolutionUsable())
    {
      if (ceres_options_.bVerbose_)
        std::cout << "Bundle Adjustment failed." << std::endl;
      return false;
    }
    else // Solution is usable
    {
      if (ceres_options_.bVerbose_)
      {
        // Display statistics about the minimization
        std::cout << std::endl
          << "Bundle Adjustment statistics (approximated RMSE) Adjust_threeCameras:\n"
          << " #views: " << sfm_data.views.size() << "\n"
          << " #poses: " << sfm_data.poses.size() << "\n"
          << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
          << " #tracks: " << sfm_data.structure.size() << "\n"
          << " #residuals: " << summary.num_residuals << "\n"
          << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
          << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
          << " Time (s): " << summary.total_time_in_seconds << "\n"
          << std::endl;


        //show reslut
        std::cout << "show reslut : " << std::endl;
        std::cout << "transformation_br : " << transformation_br[0] << " "
                                            << transformation_br[1] << " "
                                            << transformation_br[2] << " "
                                            << transformation_br[3] << " "
                                            << transformation_br[4] << " "
                                            << transformation_br[5] << std::endl;

        std::cout << "transformation_fr : " << transformation_fr[0] << " "
                                            << transformation_fr[1] << " "
                                            << transformation_fr[2] << " "
                                            << transformation_fr[3] << " "
                                            << transformation_fr[4] << " "
                                            << transformation_fr[5] << std::endl;

        std::cout << "transformation_dr : " << transformation_dr[0] << " "
                                            << transformation_dr[1] << " "
                                            << transformation_dr[2] << " "
                                            << transformation_dr[3] << " "
                                            << transformation_dr[4] << " "
                                            << transformation_dr[5] << std::endl;

        std::cout << "transformation_R_imu : " << transformation_R_imu[0] << " "
                                               << transformation_R_imu[1] << " "
                                               << transformation_R_imu[2] << std::endl;

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;

        out_transformation_R_imu[0] = transformation_R_imu[0];
        out_transformation_R_imu[1] = transformation_R_imu[1];
        out_transformation_R_imu[2] = transformation_R_imu[2];

        out_tramsformation_x_gps[0] = tramsformation_x_gps[0];
        out_tramsformation_x_gps[1] = tramsformation_x_gps[1];
        out_tramsformation_x_gps[2] = tramsformation_x_gps[2];

        out_transformation_br = transformation_br;
        out_transformation_fr = transformation_fr;
      }

      // Update camera poses with refined data
      if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
      {
        for (Poses::iterator itPose = sfm_data.poses.begin();
          itPose != sfm_data.poses.end(); ++itPose)
        {
          const IndexT indexPose = itPose->first;
          IndexT camId = indexPose;

          //single
          if (hashMap.find(camId)->second == 0)
          {
              Mat3 R_refined;
              ceres::AngleAxisToRotationMatrix(&single_poses[indexPose][0], R_refined.data());
              Vec3 t_refined(single_poses[indexPose][3], single_poses[indexPose][4], single_poses[indexPose][5]);
              // Update the pose
              Pose3 & pose = itPose->second;
              pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
          }else{  //unity
              Mat3 Rr_trans;
              Vec3 tr_trans;

              //std::vector<double> Rr_refined = reference_poses[indexPose];
              if(camId < sfm_data.views.size()/3)  //back
              {
                  camId = camId + sfm_data.views.size()/3;
                  ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rr_trans.data());

                  tr_trans(0) = transformation_br[3];
                  tr_trans(1) = transformation_br[4];
                  tr_trans(2) = transformation_br[5];

              }else if(camId >= sfm_data.views.size()/3*2)  //front
              {
                  camId = camId - sfm_data.views.size()/3;
                  ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rr_trans.data());

                  tr_trans(0) = transformation_fr[3];
                  tr_trans(1) = transformation_fr[4];
                  tr_trans(2) = transformation_fr[5];

              }else
              {
                  ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rr_trans.data());

                  tr_trans(0) = transformation_dr[3];
                  tr_trans(1) = transformation_dr[4];
                  tr_trans(2) = transformation_dr[5];

              }

              Vec3 gps_trans;
              gps_trans(0) = tramsformation_x_gps[0];
              gps_trans(1) = tramsformation_x_gps[1];
              gps_trans(2) = tramsformation_x_gps[2];
              Vec3 C_gps_res;
              C_gps_res(0) = C_gps[camId][0];
              C_gps_res(1) = C_gps[camId][1];
              C_gps_res(2) = C_gps[camId][2];
              Mat3 R_imu_trans;
              ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
              Mat3 R_imu_mat;
              ceres::AngleAxisToRotationMatrix(&R_imu[camId][0], R_imu_mat.data());

              Mat3 R_res;
              R_res = Rr_trans * R_imu_trans * R_imu_mat;
              Vec3 t_res;
              t_res = Rr_trans * gps_trans - Rr_trans * R_imu_trans * R_imu_mat * C_gps_res + tr_trans;


              //use svd





              Pose3 & pose = itPose->second;
              pose = Pose3(R_res, -R_res.transpose() * t_res);
            }
          }

      }


      // Update camera intrinsics with refined data
      if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
      {
        for (Intrinsics::iterator itIntrinsic = sfm_data.intrinsics.begin();
          itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
        {
          const IndexT indexCam = itIntrinsic->first;

          const std::vector<double> & vec_params = map_intrinsics[indexCam];
          itIntrinsic->second.get()->updateFromParams(vec_params);
        }
      }
      // Structure is already updated directly if needed (no data wrapping)
      return true;
    }



    }


bool Bundle_Adjustment_Ceres::Adjust_threecamera_gps_cail_second_change
(SfM_Data & sfm_data,     // the SfM scene to refine
 const Optimize_Options options,
 std::vector<double>& in_transformation_R_imu,
 std::vector<double>& in_tramsformation_x_gps,
 std::vector<double>& in_transformation_br,
 std::vector<double>& in_transformation_fr
)
{
    //----------
    // Add camera parameters
    // - intrinsics
    // - poses [R|t]

    // Create residuals for each observation in the bundle adjustment problem. The
    // parameters for cameras and points are added automatically.
    //----------

    ceres::Problem problem;

    //set parameters about gps imu
    Hash_Map<IndexT, std::vector<double> > R_imu;
    Hash_Map<IndexT, std::vector<double> > C_gps;
    std::vector<double> transformation_R_imu = in_transformation_R_imu;//{0.0,0.0,0.0};//{-0.0120417, 0.0201986, 0.0258046};//{0.0,0.0,0.0}; //{-0.02074, 0.0199567, 0.0271603}; //{-0.0177615, 0.0177593, 0.0196154};  //{0.0,0.0,0.0}; //{-0.000427019, 0.0125482, 0.0134619};//{-0.00594065, 0.0111603, 0.0170654}; //{0.0,0.0,0.0};//{-0.00986379, 0.00976848, 0.0253123};
    std::vector<double> tramsformation_x_gps = in_tramsformation_x_gps;//{0.0,0.0,0.0};//{-2.01425, 2.23307, -28.3916};

    double * parameter_transformation_R_imu = &transformation_R_imu[0];
    problem.AddParameterBlock(parameter_transformation_R_imu, 3);

    double * parameter_tramsformation_x_gps = &tramsformation_x_gps[0];
    problem.AddParameterBlock(parameter_tramsformation_x_gps, 3);


    //problem.SetParameterBlockConstant(parameter_transformation_R_imu);
    //problem.SetParameterBlockConstant(parameter_tramsformation_x_gps);


    //get gps and imu data
  {
        {
            //read file and initialise external
            std::vector<double> c(3);
            std::string tFilePath = sfm_data.s_root_path + "/latlon_c.res";  //!!!
            std::ifstream in;
            in.open(tFilePath);
            if(!in)
            {
                std::cout << "open tFilePath file error!" << std::endl;
                return false;
            }
            std::size_t line = 0;
            const int count = sfm_data.views.size()/3;
            while((!in.eof()) && (line < count))
            {
                 in >> c[0] >> c[1] >> c[2];
                 C_gps[line + count] = {c[0], c[1], c[2]};
                 double * parameter_block = &C_gps[line + count][0];
                 problem.AddParameterBlock(parameter_block, C_gps[line + count].size());
                 problem.SetParameterBlockConstant(parameter_block);  ///needn't change
                ++line;
            }
            in.close();

        }
        //read imu
        {
        std::ifstream in;
        in.open(sfm_data.s_root_path + "/IMU.res"); //imu file

        int line = 0;
        const int count = sfm_data.views.size()/3;

        std::ofstream out_show_camera;
        out_show_camera.open(sfm_data.s_root_path + "/show_imu_R.txt");

        while((!in.eof()) && (line < count))
        {
            double roll, pitch, yaw;
            in >> roll >> pitch >> yaw;  //z x y

            Mat3 Rz, Ry, Rx, Rx_, Ry_;
            Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
                  sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
                  0.0, 0.0, 1.0;
            Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
                  0.0, 1.0, 0.0,
                  -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
            Rx << 1.0, 0.0, 0.0,
                  0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
                  0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);

            Mat3 R_d =  Rz*Rx*Ry;

            {

//                Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
//                Vec3 c = {C_gps[line + count][0], C_gps[line + count][1], C_gps[line + count][2]};
//                axis_x = R_d.transpose() * (axis_x + R_d*c);
//                axis_y = R_d.transpose() * (axis_y + R_d*c);
//                axis_z = R_d.transpose() * (axis_z + R_d*c);
//                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                       << axis_x(0) << " " << axis_x(1) << " " << axis_x(2)<< ";" << std::endl;
//                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                       << axis_y(0) << " " << axis_y(1) << " " << axis_y(2)<< ";" << std::endl;
//                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                       << axis_z(0) << " " << axis_z(1) << " " << axis_z(2)<< ";" << std::endl;

//                out_show_camera << R_d << std::endl << std::endl;
            }

//            //change axis
//            Vec3 changeAxis_A;
//            changeAxis_A << M_PI*1/sqrt(2.0), M_PI*1/sqrt(2.0), 0.0;
            Mat3 changeAxis_M;
//            ceres::AngleAxisToRotationMatrix(changeAxis_A.data(), changeAxis_M.data());

            changeAxis_M << 1.0, 0.0, 0.0,
                            0.0, -1.0, 0.0,
                            0.0, 0.0, -1.0;
            R_d = changeAxis_M * R_d;


            {

                Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
                Vec3 c = {C_gps[line + count][0], C_gps[line + count][1], C_gps[line + count][2]};
                axis_x = R_d.transpose() * (axis_x + R_d*c);
                axis_y = R_d.transpose() * (axis_y + R_d*c);
                axis_z = R_d.transpose() * (axis_z + R_d*c);
                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
                       << axis_x(0) << " " << axis_x(1) << " " << axis_x(2)<< ";" << std::endl;
                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
                       << axis_y(0) << " " << axis_y(1) << " " << axis_y(2)<< ";" << std::endl;
                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
                       << axis_z(0) << " " << axis_z(1) << " " << axis_z(2)<< ";" << std::endl;

                out_show_camera << R_d << std::endl << std::endl;
            }





            Vec3 temp_rd;
            ceres::RotationMatrixToAngleAxis(R_d.data(), temp_rd.data());
            R_imu[line + count] = {temp_rd(0), temp_rd(1), temp_rd(2)};

            double * parameter_block = &R_imu[line + count][0];
            problem.AddParameterBlock(parameter_block, R_imu[line + count].size());

            ++line;
        }
        in.close();

        out_show_camera.close();
    }
  }


    //get hash map
    std::map<int, int> hashMap;
    filiter_ununited_pose(sfm_data.poses, sfm_data.views, hashMap);
    //Hash_Map<IndexT, std::vector<double> > reference_poses;
    Hash_Map<IndexT, std::vector<double> > single_poses;
    for (Poses::const_iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose)
    {
      const IndexT indexPose = itPose->first;

      double * parameter_block;
      if(hashMap.find(indexPose)->second == 0)
      {
          std::cout << "signal " << indexPose << std::endl;
          //set single poses
          const Pose3 & pose = itPose->second;
          const Mat3 R = pose.rotation();
          const Vec3 t = pose.translation();

          double angleAxis[3];
          ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
          // angleAxis + translation
          single_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

          parameter_block = &single_poses[indexPose][0];
          problem.AddParameterBlock(parameter_block, single_poses[indexPose].size());
      }else{
          //set unity pose
//          if( (indexPose < sfm_data.views.size()/3) || (indexPose >= sfm_data.views.size()/3*2) )
//          {
//              continue;
//          }

//          const Pose3 & pose = itPose->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 t = pose.translation();
//          Vec3 C = pose.center();

//          double angleAxis[3];
//          ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
//          // angleAxis + translation
//          reference_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(0), C(1), C(2)};

//          parameter_block = &reference_poses[indexPose][0];
//          problem.AddParameterBlock(parameter_block, reference_poses[indexPose].size());


          //test if imu right (!!!)
//          if( (indexPose < sfm_data.views.size()/3) || (indexPose >= sfm_data.views.size()/3*2) )
//          {
//              continue;
//          }

//          const Pose3 & pose = itPose->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 t = pose.translation();

//          Vec3 C = pose.center();
//          Vec3 temp_rd;
//          ceres::RotationMatrixToAngleAxis(R.data(), temp_rd.data());
//          R_imu[indexPose] = {temp_rd(0), temp_rd(1), temp_rd(2)};

//          double * parameter_block = &R_imu[indexPose][0];
//          problem.AddParameterBlock(parameter_block, R_imu[line + count].size());
//          problem.SetParameterBlockConstant(parameter_block);  ///needn't change

      }


      if(parameter_block == NULL)
      {
          continue;
      }

//      if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
//      {
//        // set the whole parameter block as constant for best performance
//        problem.SetParameterBlockConstant(parameter_block);
//      }
//      else  // Subset parametrization
//      {
//        std::vector<int> vec_constant_extrinsic;
//        // If we adjust only the translation, we must set ROTATION as constant
//        if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
//        {
//          // Subset rotation parametrization
//          vec_constant_extrinsic.push_back(0);
//          vec_constant_extrinsic.push_back(1);
//          vec_constant_extrinsic.push_back(2);
//        }
//        // If we adjust only the rotation, we must set TRANSLATION as constant
//        if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
//        {
//          // Subset translation parametrization
//          vec_constant_extrinsic.push_back(3);
//          vec_constant_extrinsic.push_back(4);
//          vec_constant_extrinsic.push_back(5);
//        }
//        if (!vec_constant_extrinsic.empty())
//        {
//          ceres::SubsetParameterization *subset_parameterization =
//            new ceres::SubsetParameterization(6, vec_constant_extrinsic);
//          problem.SetParameterization(parameter_block, subset_parameterization);
//        }
//      }

    }

///-----------------------------------------
///
///

    //test distortion
//    Mat3 r_tb, r_tf;
//    r_tb << 1.0, 0.0, 0.0,
//            0.0, 0.787926, -0.61577,
//            0.0, 0.61577, 0.787926;
//    r_tf << 1.0, 0.0, 0.0,
//            0.0, 0.78796,  0.615727,
//            0.0, -0.615727, 0.78796;
//    Vec3 RA_tb, RA_tf;
//    ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
//    ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

        //set reference
        double R_br_angleAxis[3]= {0.0, 0.0, 0.0};//{RA_tb(0), RA_tb(1), -RA_tb(2)};//{0.665224, 0.00448363, 0.00291885};  //{0.703355, 0.00263495, -0.00149232};//{0.0, 0.0, 0.0}; //{0.0212232, 0.0206259, 0.0420509};  //{0.0, 0.0, 0.0}; //{0.669661, 0.00354668, 0.00220508}; //  //{0.0, 0.0, 0.0}; //{0.662013, -0.00281273, 0.000638113 }; //{0.0, 0.0, 0.0}; //{0.671924, 0.00432461, 0.000528427};
        double R_fr_angleAxis[3]= {0.0, 0.0, 0.0};//{RA_tf(0), RA_tf(1), -RA_tf(2)};//{-0.644463, 0.0113345, -0.00156344}; //{-0.682501, 0.0147873, 0.00150509};//{0.0, 0.0, 0.0}; //{-0.0232707, -0.00288052, 0.00035795};  //{0.0, 0.0, 0.0}; //{-0.650191, 0.0101123, -0.000844533}; //  //{0.0, 0.0, 0.0}; //{-0.658446, 0.00532591, 0.0010695}; //{-0.656461, 0.00508743, -0.00299628};
        double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};
        Vec3 t_db, t_df, t_dd;
        //double S = 1470.531783220221;
        t_db << 0.0, 0.0, 0.0;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
        t_df << 0.0, 0.0, 0.0;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
        t_dd << 0.0, 0.0, 0.0;


        std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
        transformation_br = in_transformation_br;
        transformation_fr = in_transformation_fr;
        transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                             t_dd(0), t_dd(1), t_dd(2)};


    std::cout << std::endl << std::endl << "initial value: " << std::endl;
        std::cout << "transformation_br : " << transformation_br[0] << " "
                                            << transformation_br[1] << " "
                                            << transformation_br[2] << " "
                                            << transformation_br[3] << " "
                                            << transformation_br[4] << " "
                                            << transformation_br[5] << " " << std::endl;

        std::cout << "transformation_fr : " << transformation_fr[0] << " "
                                            << transformation_fr[1] << " "
                                            << transformation_fr[2] << " "
                                            << transformation_fr[3] << " "
                                            << transformation_fr[4] << " "
                                            << transformation_fr[5] << " " << std::endl;

        std::cout << "transformation_dr : " << transformation_dr[0] << " "
                                            << transformation_dr[1] << " "
                                            << transformation_dr[2] << " "
                                            << transformation_dr[3] << " "
                                            << transformation_dr[4] << " "
                                            << transformation_dr[5] << " " << std::endl;

        std::cout << "transformation_R_imu : " << transformation_R_imu[0] << " "
                                               << transformation_R_imu[1] << " "
                                               << transformation_R_imu[2] << std::endl;

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;


        double * parameter_transformation_br = &transformation_br[0];
        double * parameter_transformation_fr = &transformation_fr[0];
        double * parameter_transformation_dr = &transformation_dr[0];
        problem.AddParameterBlock(parameter_transformation_br, 6);
        problem.AddParameterBlock(parameter_transformation_fr, 6);
        problem.AddParameterBlock(parameter_transformation_dr, 6);

        problem.SetParameterBlockConstant(parameter_transformation_dr);  ///needn't change

        //problem.SetParameterBlockConstant(parameter_transformation_br);
        //problem.SetParameterBlockConstant(parameter_transformation_fr);


#ifndef CALIBRATION_THREE_CAMERAS
        //problem.SetParameterBlockConstant(parameter_transformation_br);
        //problem.SetParameterBlockConstant(parameter_transformation_fr);
#endif

        if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_transformation_br);
          problem.SetParameterBlockConstant(parameter_transformation_fr);
          problem.SetParameterBlockConstant(parameter_transformation_dr);
        }
        else  // Subset parametrization
        {
          std::vector<int> vec_constant_extrinsic_br, vec_constant_extrinsic_fr,vec_constant_extrinsic_dr;
          // If we adjust only the translation, we must set ROTATION as constant
          if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
          {
            // Subset rotation parametrization
            vec_constant_extrinsic_br.push_back(0);
            vec_constant_extrinsic_br.push_back(1);
            vec_constant_extrinsic_br.push_back(2);
            vec_constant_extrinsic_fr.push_back(0);
            vec_constant_extrinsic_fr.push_back(1);
            vec_constant_extrinsic_fr.push_back(2);
            vec_constant_extrinsic_dr.push_back(0);
            vec_constant_extrinsic_dr.push_back(1);
            vec_constant_extrinsic_dr.push_back(2);
          }
          // If we adjust only the rotation, we must set TRANSLATION as constant
          if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
          {
            // Subset translation parametrization
              vec_constant_extrinsic_br.push_back(3);
              vec_constant_extrinsic_br.push_back(4);
              vec_constant_extrinsic_br.push_back(5);
              vec_constant_extrinsic_fr.push_back(3);
              vec_constant_extrinsic_fr.push_back(4);
              vec_constant_extrinsic_fr.push_back(5);
              vec_constant_extrinsic_dr.push_back(3);
              vec_constant_extrinsic_dr.push_back(4);
              vec_constant_extrinsic_dr.push_back(5);
          }
          if (!vec_constant_extrinsic_br.empty())
          {
            ceres::SubsetParameterization *subset_parameterization_br =
              new ceres::SubsetParameterization(6, vec_constant_extrinsic_br);
            problem.SetParameterization(parameter_transformation_br, subset_parameterization_br);

            ceres::SubsetParameterization *subset_parameterization_fr =
              new ceres::SubsetParameterization(6, vec_constant_extrinsic_fr);
            problem.SetParameterization(parameter_transformation_fr, subset_parameterization_fr);

            ceres::SubsetParameterization *subset_parameterization_dr =
              new ceres::SubsetParameterization(6, vec_constant_extrinsic_dr);
            problem.SetParameterization(parameter_transformation_dr, subset_parameterization_dr);
          }
        }


    // Data wrapper for refinement:
    Hash_Map<IndexT, std::vector<double> > map_intrinsics;

   // Setup Intrinsics data & subparametrization
    for (Intrinsics::const_iterator itIntrinsic = sfm_data.intrinsics.begin();
      itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
    {
      const IndexT indexCam = itIntrinsic->first;

      if (isValid(itIntrinsic->second->getType()))
      {
        map_intrinsics[indexCam] = itIntrinsic->second->getParams();

        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            itIntrinsic->second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
      else
      {
        std::cerr << "Unsupported camera type." << std::endl;
      }
    }


    // Set a LossFunction to be less penalized by false measurements
    //  - set it to NULL if you don't want use a lossFunction.
    ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(Square(4.0));
    // TODO: make the LOSS function and the parameter an option

std::cout << "start." << std::endl;

    // For all visibility add reprojections errors:
    for (Landmarks::iterator iterTracks = sfm_data.structure.begin();
         iterTracks!= sfm_data.structure.end(); ++iterTracks)
    {
      const Observations & obs = iterTracks->second.obs;

      for (Observations::const_iterator itObs = obs.begin();
        itObs != obs.end(); ++itObs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(itObs->first).get();

        //decide camera pose through pose id
        double * parameter_transformation;
        IndexT poseId = view->id_pose;


        //single
        if (hashMap.find(poseId)->second == 0)
        {

            // Each Residual block takes a point and a camera as input and outputs a 2
            // dimensional residual. Internally, the cost function stores the observed
            // image location and compares the reprojection against the observation.
            ceres::CostFunction* cost_function =
              IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), itObs->second.x);

            if (cost_function)
              problem.AddResidualBlock(cost_function,
                p_LossFunction,
                &map_intrinsics[view->id_intrinsic][0],
                &single_poses[view->id_pose][0],
                iterTracks->second.X.data());

        }else{  //unity
            if(poseId < sfm_data.views.size()/3)  //back
            {
                poseId = poseId + sfm_data.views.size()/3;

//                  ceres::CostFunction* cost_function2 =
//                          new ceres::AutoDiffCostFunction
//                            <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
//                              new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

//                  if (cost_function2)
//                    problem.AddResidualBlock(cost_function2,
//                      p_LossFunction,
//                      &map_intrinsics[view->id_intrinsic][0],
//                      &R_imu[poseId][0],
//                      &C_gps[poseId][0],
//                      &transformation_R_imu[0],
//                      &tramsformation_x_gps[0],
//                      &transformation_br[0],
//                      iterTracks->second.X.data());  //


                  ceres::CostFunction* cost_function2 =
                          new ceres::AutoDiffCostFunction
                            <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
                              new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

                  if (cost_function2)
                    problem.AddResidualBlock(cost_function2,
                      p_LossFunction,
                      &map_intrinsics[view->id_intrinsic][0],
                      &R_imu[poseId][0],
                      &C_gps[poseId][0],
                      &transformation_R_imu[0],
                      &tramsformation_x_gps[0],
                      &transformation_br[0],
                      iterTracks->second.X.data());  //




            }else if(poseId >= sfm_data.views.size()/3*2)  //front
            {
                poseId = poseId - sfm_data.views.size()/3;

//                ceres::CostFunction* cost_function2 =
//                        new ceres::AutoDiffCostFunction
//                          <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
//                            new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

//                if (cost_function2)
//                  problem.AddResidualBlock(cost_function2,
//                    p_LossFunction,
//                    &map_intrinsics[view->id_intrinsic][0],
//                    &R_imu[poseId][0],
//                    &C_gps[poseId][0],
//                    &transformation_R_imu[0],
//                    &tramsformation_x_gps[0],
//                    &transformation_fr[0],
//                    iterTracks->second.X.data());  //


                ceres::CostFunction* cost_function2 =
                        new ceres::AutoDiffCostFunction
                          <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
                            new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

                if (cost_function2)
                  problem.AddResidualBlock(cost_function2,
                    p_LossFunction,
                    &map_intrinsics[view->id_intrinsic][0],
                    &R_imu[poseId][0],
                    &C_gps[poseId][0],
                    &transformation_R_imu[0],
                    &tramsformation_x_gps[0],
                    &transformation_fr[0],
                    iterTracks->second.X.data());  //


            }else{  //down

//                ceres::CostFunction* cost_function2 =
//                        new ceres::AutoDiffCostFunction
//                          <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
//                            new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

//                if (cost_function2)
//                  problem.AddResidualBlock(cost_function2,
//                    p_LossFunction,
//                    &map_intrinsics[view->id_intrinsic][0],
//                    &R_imu[poseId][0],
//                    &C_gps[poseId][0],
//                    &transformation_R_imu[0],
//                    &tramsformation_x_gps[0],
//                    &transformation_dr[0],
//                    iterTracks->second.X.data());  //



                ceres::CostFunction* cost_function2 =
                        new ceres::AutoDiffCostFunction
                          <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
                            new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

                if (cost_function2)
                  problem.AddResidualBlock(cost_function2,
                    p_LossFunction,
                    &map_intrinsics[view->id_intrinsic][0],
                    &R_imu[poseId][0],
                    &C_gps[poseId][0],
                    &transformation_R_imu[0],
                    &tramsformation_x_gps[0],
                    &transformation_dr[0],
                    iterTracks->second.X.data());  //


            }

        }

      }
            if (options.structure_opt == Structure_Parameter_Type::NONE)
              problem.SetParameterBlockConstant(iterTracks->second.X.data());
    }



    // Configure a BA engine and run it
    //  Make Ceres automatically detect the bundle structure.
    ceres::Solver::Options ceres_config_options;
    ceres_config_options.max_num_iterations = 500;
    ceres_config_options.preconditioner_type =
      static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
    ceres_config_options.linear_solver_type =
      static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
    ceres_config_options.sparse_linear_algebra_library_type =
      static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
    ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
    ceres_config_options.logging_type = ceres::SILENT;
    ceres_config_options.num_threads = ceres_options_.nb_threads_;
    ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
    ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;


    // Solve BA
    ceres::Solver::Summary summary;
    ceres::Solve(ceres_config_options, &problem, &summary);
    if (ceres_options_.bCeres_summary_)
      std::cout << summary.FullReport() << std::endl;


    // If no error, get back refined parameters
    if (!summary.IsSolutionUsable())
    {
      if (ceres_options_.bVerbose_)
        std::cout << "Bundle Adjustment failed." << std::endl;
      return false;
    }
    else // Solution is usable
    {
      if (ceres_options_.bVerbose_)
      {
        // Display statistics about the minimization
        std::cout << std::endl
          << "Bundle Adjustment statistics (approximated RMSE) Adjust_threeCameras:\n"
          << " #views: " << sfm_data.views.size() << "\n"
          << " #poses: " << sfm_data.poses.size() << "\n"
          << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
          << " #tracks: " << sfm_data.structure.size() << "\n"
          << " #residuals: " << summary.num_residuals << "\n"
          << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
          << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
          << " Time (s): " << summary.total_time_in_seconds << "\n"
          << std::endl;


        //show reslut
        std::cout << "show reslut : " << std::endl;
        std::cout << "transformation_br : " << transformation_br[0] << " "
                                            << transformation_br[1] << " "
                                            << transformation_br[2] << " "
                                            << transformation_br[3] << " "
                                            << transformation_br[4] << " "
                                            << transformation_br[5] << std::endl;

        std::cout << "transformation_fr : " << transformation_fr[0] << " "
                                            << transformation_fr[1] << " "
                                            << transformation_fr[2] << " "
                                            << transformation_fr[3] << " "
                                            << transformation_fr[4] << " "
                                            << transformation_fr[5] << std::endl;

        std::cout << "transformation_dr : " << transformation_dr[0] << " "
                                            << transformation_dr[1] << " "
                                            << transformation_dr[2] << " "
                                            << transformation_dr[3] << " "
                                            << transformation_dr[4] << " "
                                            << transformation_dr[5] << std::endl;

        std::cout << "transformation_R_imu : " << transformation_R_imu[0] << " "
                                               << transformation_R_imu[1] << " "
                                               << transformation_R_imu[2] << std::endl;

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;

        in_transformation_R_imu[0] = transformation_R_imu[0];
        in_transformation_R_imu[1] = transformation_R_imu[1];
        in_transformation_R_imu[2] = transformation_R_imu[2];

        in_tramsformation_x_gps[0] = tramsformation_x_gps[0];
        in_tramsformation_x_gps[1] = tramsformation_x_gps[1];
        in_tramsformation_x_gps[2] = tramsformation_x_gps[2];

        in_transformation_br = transformation_br;
        in_transformation_fr = transformation_fr;
      }

      // Update camera poses with refined data
      if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
      {
        for (Poses::iterator itPose = sfm_data.poses.begin();
          itPose != sfm_data.poses.end(); ++itPose)
        {
          const IndexT indexPose = itPose->first;
          IndexT camId = indexPose;

          //single
          if (hashMap.find(camId)->second == 0)
          {
              Mat3 R_refined;
              ceres::AngleAxisToRotationMatrix(&single_poses[indexPose][0], R_refined.data());
              Vec3 t_refined(single_poses[indexPose][3], single_poses[indexPose][4], single_poses[indexPose][5]);
              // Update the pose
              Pose3 & pose = itPose->second;
              pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
          }else{  //unity
              Mat3 Rr_trans;
              Vec3 tr_trans;

              //std::vector<double> Rr_refined = reference_poses[indexPose];
              if(camId < sfm_data.views.size()/3)  //back
              {
                  camId = camId + sfm_data.views.size()/3;
                  ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rr_trans.data());

                  tr_trans(0) = transformation_br[3];
                  tr_trans(1) = transformation_br[4];
                  tr_trans(2) = transformation_br[5];

              }else if(camId >= sfm_data.views.size()/3*2)  //front
              {
                  camId = camId - sfm_data.views.size()/3;
                  ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rr_trans.data());

                  tr_trans(0) = transformation_fr[3];
                  tr_trans(1) = transformation_fr[4];
                  tr_trans(2) = transformation_fr[5];

              }else
              {
                  ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rr_trans.data());

                  tr_trans(0) = transformation_dr[3];
                  tr_trans(1) = transformation_dr[4];
                  tr_trans(2) = transformation_dr[5];

              }

              Vec3 gps_trans;
              gps_trans(0) = tramsformation_x_gps[0];
              gps_trans(1) = tramsformation_x_gps[1];
              gps_trans(2) = tramsformation_x_gps[2];
              Vec3 C_gps_res;
              C_gps_res(0) = C_gps[camId][0];
              C_gps_res(1) = C_gps[camId][1];
              C_gps_res(2) = C_gps[camId][2];
              Mat3 R_imu_trans;
              ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
              Mat3 R_imu_mat;
              ceres::AngleAxisToRotationMatrix(&R_imu[camId][0], R_imu_mat.data());

              Mat3 R_res;
              R_res = Rr_trans * R_imu_trans * R_imu_mat;
              Vec3 t_res;
              t_res = Rr_trans * gps_trans - Rr_trans * R_imu_trans * R_imu_mat * C_gps_res + tr_trans;


              //use svd





              Pose3 & pose = itPose->second;
              pose = Pose3(R_res, -R_res.transpose() * t_res);
            }
          }

      }


      // Update camera intrinsics with refined data
      if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
      {
        for (Intrinsics::iterator itIntrinsic = sfm_data.intrinsics.begin();
          itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
        {
          const IndexT indexCam = itIntrinsic->first;

          const std::vector<double> & vec_params = map_intrinsics[indexCam];
          itIntrinsic->second.get()->updateFromParams(vec_params);
        }
      }
      // Structure is already updated directly if needed (no data wrapping)
      return true;
    }



    }

bool Bundle_Adjustment_Ceres::Adjust_threecamera_gps_cail_first_change_c4
(SfM_Data & sfm_data,     // the SfM scene to refine
 const Optimize_Options options,
 std::vector<double>& out_transformation_R_imu,
 std::vector<double>& out_tramsformation_x_gps,
 std::vector<double>& out_transformation_br,
 std::vector<double>& out_transformation_fr
)
{
    //----------
    // Add camera parameters
    // - intrinsics
    // - poses [R|t]

    // Create residuals for each observation in the bundle adjustment problem. The
    // parameters for cameras and points are added automatically.
    //----------

    ceres::Problem problem;

    //set parameters about gps imu
    Hash_Map<IndexT, std::vector<double> > R_imu;
    Hash_Map<IndexT, std::vector<double> > C_gps;
    std::vector<double> transformation_R_imu = {0.0,0.0,0.0};//{-0.02074, 0.0199567, 0.0271603}; //{-0.0177615, 0.0177593, 0.0196154};  //{0.0,0.0,0.0}; //{-0.000427019, 0.0125482, 0.0134619};//{-0.00594065, 0.0111603, 0.0170654}; //{0.0,0.0,0.0};//{-0.00986379, 0.00976848, 0.0253123};
    std::vector<double> tramsformation_x_gps = {0.0,0.0,0.0}; //{-0.82937, 0.462316, -2.59147}; //{0.0,0.0,0.0};//{-0.281802, 0.574862, 0.064961};

    double * parameter_transformation_R_imu = &transformation_R_imu[0];
    problem.AddParameterBlock(parameter_transformation_R_imu, 3);

    double * parameter_tramsformation_x_gps = &tramsformation_x_gps[0];
    problem.AddParameterBlock(parameter_tramsformation_x_gps, 3);


    //problem.SetParameterBlockConstant(parameter_transformation_R_imu);
    problem.SetParameterBlockConstant(parameter_tramsformation_x_gps);


    //get gps and imu data
  {
        {
            //read file and initialise external
            std::vector<double> c(3);
            std::string tFilePath = sfm_data.s_root_path + "/latlon_c.res";  //!!!
            std::ifstream in;
            in.open(tFilePath);
            if(!in)
            {
                std::cout << "open tFilePath file error!" << std::endl;
                return false;
            }
            std::size_t line = 0;
            const int count = sfm_data.views.size()/3;
            while((!in.eof()) && (line < count))
            {
                 in >> c[0] >> c[1] >> c[2];
                 C_gps[line + count] = {c[0], c[1], c[2]};
                 double * parameter_block = &C_gps[line + count][0];
                 problem.AddParameterBlock(parameter_block, C_gps[line + count].size());
                 problem.SetParameterBlockConstant(parameter_block);  ///needn't change
                ++line;
            }
            in.close();

        }
        //read imu
        {
        std::ifstream in;
        in.open(sfm_data.s_root_path + "/IMU.res"); //imu file

        int line = 0;
        const int count = sfm_data.views.size()/3;

        std::ofstream out_show_camera;
        out_show_camera.open(sfm_data.s_root_path + "/show_imu_R.txt");

        while((!in.eof()) && (line < count))
        {
            double roll, pitch, yaw;
            in >> roll >> pitch >> yaw;  //z x y

            Mat3 Rz, Ry, Rx, Rx_, Ry_;
            Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
                  sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
                  0.0, 0.0, 1.0;
            Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
                  0.0, 1.0, 0.0,
                  -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
            Rx << 1.0, 0.0, 0.0,
                  0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
                  0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);

            Mat3 R_d =  Rz*Rx*Ry;

            {

//                Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
//                Vec3 c = {C_gps[line + count][0], C_gps[line + count][1], C_gps[line + count][2]};
//                axis_x = R_d.transpose() * (axis_x + R_d*c);
//                axis_y = R_d.transpose() * (axis_y + R_d*c);
//                axis_z = R_d.transpose() * (axis_z + R_d*c);
//                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                       << axis_x(0) << " " << axis_x(1) << " " << axis_x(2)<< ";" << std::endl;
//                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                       << axis_y(0) << " " << axis_y(1) << " " << axis_y(2)<< ";" << std::endl;
//                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                       << axis_z(0) << " " << axis_z(1) << " " << axis_z(2)<< ";" << std::endl;

//                out_show_camera << R_d << std::endl << std::endl;
            }

//            //change axis
//            Vec3 changeAxis_A;
//            changeAxis_A << M_PI*1/sqrt(2.0), M_PI*1/sqrt(2.0), 0.0;
            Mat3 changeAxis_M;
//            ceres::AngleAxisToRotationMatrix(changeAxis_A.data(), changeAxis_M.data());

            changeAxis_M << 1.0, 0.0, 0.0,
                            0.0, -1.0, 0.0,
                            0.0, 0.0, -1.0;
            R_d = changeAxis_M * R_d;


            {

                Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
                Vec3 c = {C_gps[line + count][0], C_gps[line + count][1], C_gps[line + count][2]};
                axis_x = R_d.transpose() * (axis_x + R_d*c);
                axis_y = R_d.transpose() * (axis_y + R_d*c);
                axis_z = R_d.transpose() * (axis_z + R_d*c);
                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
                       << axis_x(0) << " " << axis_x(1) << " " << axis_x(2)<< ";" << std::endl;
                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
                       << axis_y(0) << " " << axis_y(1) << " " << axis_y(2)<< ";" << std::endl;
                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
                       << axis_z(0) << " " << axis_z(1) << " " << axis_z(2)<< ";" << std::endl;

                out_show_camera << R_d << std::endl << std::endl;
            }





            Vec3 temp_rd;
            ceres::RotationMatrixToAngleAxis(R_d.data(), temp_rd.data());
            R_imu[line + count] = {temp_rd(0), temp_rd(1), temp_rd(2)};

            double * parameter_block = &R_imu[line + count][0];
            problem.AddParameterBlock(parameter_block, R_imu[line + count].size());

            ++line;
        }
        in.close();

        out_show_camera.close();
    }
  }


    //get hash map
    std::map<int, int> hashMap;
    filiter_ununited_pose(sfm_data.poses, sfm_data.views, hashMap);
    //Hash_Map<IndexT, std::vector<double> > reference_poses;
    Hash_Map<IndexT, std::vector<double> > single_poses;
    for (Poses::const_iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose)
    {
      const IndexT indexPose = itPose->first;

      double * parameter_block;
      if(hashMap.find(indexPose)->second == 0)
      {
          std::cout << "signal " << indexPose << std::endl;
          //set single poses
          const Pose3 & pose = itPose->second;
          const Mat3 R = pose.rotation();
          const Vec3 t = pose.translation();

          double angleAxis[3];
          ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
          // angleAxis + translation
          single_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

          parameter_block = &single_poses[indexPose][0];
          problem.AddParameterBlock(parameter_block, single_poses[indexPose].size());
      }else{
          //set unity pose
//          if( (indexPose < sfm_data.views.size()/3) || (indexPose >= sfm_data.views.size()/3*2) )
//          {
//              continue;
//          }

//          const Pose3 & pose = itPose->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 t = pose.translation();
//          Vec3 C = pose.center();

//          double angleAxis[3];
//          ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
//          // angleAxis + translation
//          reference_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(0), C(1), C(2)};

//          parameter_block = &reference_poses[indexPose][0];
//          problem.AddParameterBlock(parameter_block, reference_poses[indexPose].size());


          //test if imu right (!!!)
//          if( (indexPose < sfm_data.views.size()/3) || (indexPose >= sfm_data.views.size()/3*2) )
//          {
//              continue;
//          }

//          const Pose3 & pose = itPose->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 t = pose.translation();

//          Vec3 C = pose.center();
//          Vec3 temp_rd;
//          ceres::RotationMatrixToAngleAxis(R.data(), temp_rd.data());
//          R_imu[indexPose] = {temp_rd(0), temp_rd(1), temp_rd(2)};

//          double * parameter_block = &R_imu[indexPose][0];
//          problem.AddParameterBlock(parameter_block, R_imu[line + count].size());
//          problem.SetParameterBlockConstant(parameter_block);  ///needn't change

      }


      if(parameter_block == NULL)
      {
          continue;
      }

//      if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
//      {
//        // set the whole parameter block as constant for best performance
//        problem.SetParameterBlockConstant(parameter_block);
//      }
//      else  // Subset parametrization
//      {
//        std::vector<int> vec_constant_extrinsic;
//        // If we adjust only the translation, we must set ROTATION as constant
//        if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
//        {
//          // Subset rotation parametrization
//          vec_constant_extrinsic.push_back(0);
//          vec_constant_extrinsic.push_back(1);
//          vec_constant_extrinsic.push_back(2);
//        }
//        // If we adjust only the rotation, we must set TRANSLATION as constant
//        if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
//        {
//          // Subset translation parametrization
//          vec_constant_extrinsic.push_back(3);
//          vec_constant_extrinsic.push_back(4);
//          vec_constant_extrinsic.push_back(5);
//        }
//        if (!vec_constant_extrinsic.empty())
//        {
//          ceres::SubsetParameterization *subset_parameterization =
//            new ceres::SubsetParameterization(6, vec_constant_extrinsic);
//          problem.SetParameterization(parameter_block, subset_parameterization);
//        }
//      }

    }

///-----------------------------------------
///
///

    //test distortion
    Mat3 r_tb, r_tf;
    r_tb << 1.0, 0.0, 0.0,
            0.0, 0.787926, -0.61577,
            0.0, 0.61577, 0.787926;
    r_tf << 1.0, 0.0, 0.0,
            0.0, 0.78796,  0.615727,
            0.0, -0.615727, 0.78796;
    Vec3 RA_tb, RA_tf;
    ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
    ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

        //set reference
        double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), -RA_tb(2)};//{0.665224, 0.00448363, 0.00291885};  //{0.703355, 0.00263495, -0.00149232};//{0.0, 0.0, 0.0}; //{0.0212232, 0.0206259, 0.0420509};  //{0.0, 0.0, 0.0}; //{0.669661, 0.00354668, 0.00220508}; //  //{0.0, 0.0, 0.0}; //{0.662013, -0.00281273, 0.000638113 }; //{0.0, 0.0, 0.0}; //{0.671924, 0.00432461, 0.000528427};
        double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), -RA_tf(2)};//{-0.644463, 0.0113345, -0.00156344}; //{-0.682501, 0.0147873, 0.00150509};//{0.0, 0.0, 0.0}; //{-0.0232707, -0.00288052, 0.00035795};  //{0.0, 0.0, 0.0}; //{-0.650191, 0.0101123, -0.000844533}; //  //{0.0, 0.0, 0.0}; //{-0.658446, 0.00532591, 0.0010695}; //{-0.656461, 0.00508743, -0.00299628};
        double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};
        Vec3 t_db, t_df, t_dd;
        //double S = 1470.531783220221;
        t_db << 0.0, 0.0, 0.0;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
        t_df << 0.0, 0.0, 0.0;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
        t_dd << 0.0, 0.0, 0.0;


        std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
        transformation_br = {R_br_angleAxis[0], R_br_angleAxis[1], R_br_angleAxis[2],
                             t_db(0), t_db(1), t_db(2)};
        transformation_fr = {R_fr_angleAxis[0], R_fr_angleAxis[1], R_fr_angleAxis[2],
                             t_df(0), t_df(1), t_df(2)};
        transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                             t_dd(0), t_dd(1), t_dd(2)};


    std::cout << std::endl << std::endl << "initial value: " << std::endl;
        std::cout << "transformation_br : " << transformation_br[0] << " "
                                            << transformation_br[1] << " "
                                            << transformation_br[2] << " "
                                            << transformation_br[3] << " "
                                            << transformation_br[4] << " "
                                            << transformation_br[5] << " " << std::endl;

        std::cout << "transformation_fr : " << transformation_fr[0] << " "
                                            << transformation_fr[1] << " "
                                            << transformation_fr[2] << " "
                                            << transformation_fr[3] << " "
                                            << transformation_fr[4] << " "
                                            << transformation_fr[5] << " " << std::endl;

        std::cout << "transformation_dr : " << transformation_dr[0] << " "
                                            << transformation_dr[1] << " "
                                            << transformation_dr[2] << " "
                                            << transformation_dr[3] << " "
                                            << transformation_dr[4] << " "
                                            << transformation_dr[5] << " " << std::endl;

        std::cout << "transformation_R_imu : " << transformation_R_imu[0] << " "
                                               << transformation_R_imu[1] << " "
                                               << transformation_R_imu[2] << std::endl;

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;


        double * parameter_transformation_br = &transformation_br[0];
        double * parameter_transformation_fr = &transformation_fr[0];
        double * parameter_transformation_dr = &transformation_dr[0];
        problem.AddParameterBlock(parameter_transformation_br, 6);
        problem.AddParameterBlock(parameter_transformation_fr, 6);
        problem.AddParameterBlock(parameter_transformation_dr, 6);

        problem.SetParameterBlockConstant(parameter_transformation_dr);  ///needn't change

        //problem.SetParameterBlockConstant(parameter_transformation_br);
        //problem.SetParameterBlockConstant(parameter_transformation_fr);


#ifndef CALIBRATION_THREE_CAMERAS
        //problem.SetParameterBlockConstant(parameter_transformation_br);
        //problem.SetParameterBlockConstant(parameter_transformation_fr);
#endif

        if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_transformation_br);
          problem.SetParameterBlockConstant(parameter_transformation_fr);
          problem.SetParameterBlockConstant(parameter_transformation_dr);
        }
        else  // Subset parametrization
        {
          std::vector<int> vec_constant_extrinsic_br, vec_constant_extrinsic_fr,vec_constant_extrinsic_dr;
          // If we adjust only the translation, we must set ROTATION as constant
          if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
          {
            // Subset rotation parametrization
            vec_constant_extrinsic_br.push_back(0);
            vec_constant_extrinsic_br.push_back(1);
            vec_constant_extrinsic_br.push_back(2);
            vec_constant_extrinsic_fr.push_back(0);
            vec_constant_extrinsic_fr.push_back(1);
            vec_constant_extrinsic_fr.push_back(2);
            vec_constant_extrinsic_dr.push_back(0);
            vec_constant_extrinsic_dr.push_back(1);
            vec_constant_extrinsic_dr.push_back(2);
          }
          // If we adjust only the rotation, we must set TRANSLATION as constant
          if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
          {
            // Subset translation parametrization
              vec_constant_extrinsic_br.push_back(3);
              vec_constant_extrinsic_br.push_back(4);
              vec_constant_extrinsic_br.push_back(5);
              vec_constant_extrinsic_fr.push_back(3);
              vec_constant_extrinsic_fr.push_back(4);
              vec_constant_extrinsic_fr.push_back(5);
              vec_constant_extrinsic_dr.push_back(3);
              vec_constant_extrinsic_dr.push_back(4);
              vec_constant_extrinsic_dr.push_back(5);
          }
          if (!vec_constant_extrinsic_br.empty())
          {
            ceres::SubsetParameterization *subset_parameterization_br =
              new ceres::SubsetParameterization(6, vec_constant_extrinsic_br);
            problem.SetParameterization(parameter_transformation_br, subset_parameterization_br);

            ceres::SubsetParameterization *subset_parameterization_fr =
              new ceres::SubsetParameterization(6, vec_constant_extrinsic_fr);
            problem.SetParameterization(parameter_transformation_fr, subset_parameterization_fr);

            ceres::SubsetParameterization *subset_parameterization_dr =
              new ceres::SubsetParameterization(6, vec_constant_extrinsic_dr);
            problem.SetParameterization(parameter_transformation_dr, subset_parameterization_dr);
          }
        }


    // Data wrapper for refinement:
    Hash_Map<IndexT, std::vector<double> > map_intrinsics;

   // Setup Intrinsics data & subparametrization
    for (Intrinsics::const_iterator itIntrinsic = sfm_data.intrinsics.begin();
      itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
    {
      const IndexT indexCam = itIntrinsic->first;

      if (isValid(itIntrinsic->second->getType()))
      {
        map_intrinsics[indexCam] = itIntrinsic->second->getParams();

        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            itIntrinsic->second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
      else
      {
        std::cerr << "Unsupported camera type." << std::endl;
      }
    }


    // Set a LossFunction to be less penalized by false measurements
    //  - set it to NULL if you don't want use a lossFunction.
    ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(Square(4.0));
    // TODO: make the LOSS function and the parameter an option

std::cout << "start." << std::endl;

    // For all visibility add reprojections errors:
    for (Landmarks::iterator iterTracks = sfm_data.structure.begin();
         iterTracks!= sfm_data.structure.end(); ++iterTracks)
    {
      const Observations & obs = iterTracks->second.obs;

      for (Observations::const_iterator itObs = obs.begin();
        itObs != obs.end(); ++itObs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(itObs->first).get();

        //decide camera pose through pose id
        double * parameter_transformation;
        IndexT poseId = view->id_pose;


        //single
        if (hashMap.find(poseId)->second == 0)
        {

            // Each Residual block takes a point and a camera as input and outputs a 2
            // dimensional residual. Internally, the cost function stores the observed
            // image location and compares the reprojection against the observation.
            ceres::CostFunction* cost_function =
              IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), itObs->second.x);

            if (cost_function)
              problem.AddResidualBlock(cost_function,
                p_LossFunction,
                &map_intrinsics[view->id_intrinsic][0],
                &single_poses[view->id_pose][0],
                iterTracks->second.X.data());

        }else{  //unity
            if(poseId < sfm_data.views.size()/3)  //back
            {
                poseId = poseId + sfm_data.views.size()/3;

//                  ceres::CostFunction* cost_function2 =
//                          new ceres::AutoDiffCostFunction
//                            <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
//                              new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

//                  if (cost_function2)
//                    problem.AddResidualBlock(cost_function2,
//                      p_LossFunction,
//                      &map_intrinsics[view->id_intrinsic][0],
//                      &R_imu[poseId][0],
//                      &C_gps[poseId][0],
//                      &transformation_R_imu[0],
//                      &tramsformation_x_gps[0],
//                      &transformation_br[0],
//                      iterTracks->second.X.data());  //


                  ceres::CostFunction* cost_function2 =
                          new ceres::AutoDiffCostFunction
                            <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_ThreeCamera_IMU_GPS, 2, 8, 3, 3, 3, 3, 6, 3>(
                              new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_ThreeCamera_IMU_GPS(itObs->second.x.data()));

                  if (cost_function2)
                    problem.AddResidualBlock(cost_function2,
                      p_LossFunction,
                      &map_intrinsics[view->id_intrinsic][0],
                      &R_imu[poseId][0],
                      &C_gps[poseId][0],
                      &transformation_R_imu[0],
                      &tramsformation_x_gps[0],
                      &transformation_br[0],
                      iterTracks->second.X.data());  //




            }else if(poseId >= sfm_data.views.size()/3*2)  //front
            {
                poseId = poseId - sfm_data.views.size()/3;

//                ceres::CostFunction* cost_function2 =
//                        new ceres::AutoDiffCostFunction
//                          <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
//                            new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

//                if (cost_function2)
//                  problem.AddResidualBlock(cost_function2,
//                    p_LossFunction,
//                    &map_intrinsics[view->id_intrinsic][0],
//                    &R_imu[poseId][0],
//                    &C_gps[poseId][0],
//                    &transformation_R_imu[0],
//                    &tramsformation_x_gps[0],
//                    &transformation_fr[0],
//                    iterTracks->second.X.data());  //


                ceres::CostFunction* cost_function2 =
                        new ceres::AutoDiffCostFunction
                          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_ThreeCamera_IMU_GPS, 2, 8, 3, 3, 3, 3, 6, 3>(
                            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_ThreeCamera_IMU_GPS(itObs->second.x.data()));

                if (cost_function2)
                  problem.AddResidualBlock(cost_function2,
                    p_LossFunction,
                    &map_intrinsics[view->id_intrinsic][0],
                    &R_imu[poseId][0],
                    &C_gps[poseId][0],
                    &transformation_R_imu[0],
                    &tramsformation_x_gps[0],
                    &transformation_fr[0],
                    iterTracks->second.X.data());  //


            }else{  //down

//                ceres::CostFunction* cost_function2 =
//                        new ceres::AutoDiffCostFunction
//                          <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
//                            new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

//                if (cost_function2)
//                  problem.AddResidualBlock(cost_function2,
//                    p_LossFunction,
//                    &map_intrinsics[view->id_intrinsic][0],
//                    &R_imu[poseId][0],
//                    &C_gps[poseId][0],
//                    &transformation_R_imu[0],
//                    &tramsformation_x_gps[0],
//                    &transformation_dr[0],
//                    iterTracks->second.X.data());  //



                ceres::CostFunction* cost_function2 =
                        new ceres::AutoDiffCostFunction
                          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_ThreeCamera_IMU_GPS, 2, 8, 3, 3, 3, 3, 6, 3>(
                            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_ThreeCamera_IMU_GPS(itObs->second.x.data()));

                if (cost_function2)
                  problem.AddResidualBlock(cost_function2,
                    p_LossFunction,
                    &map_intrinsics[view->id_intrinsic][0],
                    &R_imu[poseId][0],
                    &C_gps[poseId][0],
                    &transformation_R_imu[0],
                    &tramsformation_x_gps[0],
                    &transformation_dr[0],
                    iterTracks->second.X.data());  //


            }

        }

      }
            if (options.structure_opt == Structure_Parameter_Type::NONE)
              problem.SetParameterBlockConstant(iterTracks->second.X.data());
    }



    // Configure a BA engine and run it
    //  Make Ceres automatically detect the bundle structure.
    ceres::Solver::Options ceres_config_options;
    ceres_config_options.max_num_iterations = 500;
    ceres_config_options.preconditioner_type =
      static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
    ceres_config_options.linear_solver_type =
      static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
    ceres_config_options.sparse_linear_algebra_library_type =
      static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
    ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
    ceres_config_options.logging_type = ceres::SILENT;
    ceres_config_options.num_threads = ceres_options_.nb_threads_;
    ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
    ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;


    // Solve BA
    ceres::Solver::Summary summary;
    ceres::Solve(ceres_config_options, &problem, &summary);
    if (ceres_options_.bCeres_summary_)
      std::cout << summary.FullReport() << std::endl;


    // If no error, get back refined parameters
    if (!summary.IsSolutionUsable())
    {
      if (ceres_options_.bVerbose_)
        std::cout << "Bundle Adjustment failed." << std::endl;
      return false;
    }
    else // Solution is usable
    {
      if (ceres_options_.bVerbose_)
      {
        // Display statistics about the minimization
        std::cout << std::endl
          << "Bundle Adjustment statistics (approximated RMSE) Adjust_threeCameras:\n"
          << " #views: " << sfm_data.views.size() << "\n"
          << " #poses: " << sfm_data.poses.size() << "\n"
          << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
          << " #tracks: " << sfm_data.structure.size() << "\n"
          << " #residuals: " << summary.num_residuals << "\n"
          << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
          << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
          << " Time (s): " << summary.total_time_in_seconds << "\n"
          << std::endl;


        //show reslut
        std::cout << "show reslut : " << std::endl;
        std::cout << "transformation_br : " << transformation_br[0] << " "
                                            << transformation_br[1] << " "
                                            << transformation_br[2] << " "
                                            << transformation_br[3] << " "
                                            << transformation_br[4] << " "
                                            << transformation_br[5] << std::endl;

        std::cout << "transformation_fr : " << transformation_fr[0] << " "
                                            << transformation_fr[1] << " "
                                            << transformation_fr[2] << " "
                                            << transformation_fr[3] << " "
                                            << transformation_fr[4] << " "
                                            << transformation_fr[5] << std::endl;

        std::cout << "transformation_dr : " << transformation_dr[0] << " "
                                            << transformation_dr[1] << " "
                                            << transformation_dr[2] << " "
                                            << transformation_dr[3] << " "
                                            << transformation_dr[4] << " "
                                            << transformation_dr[5] << std::endl;

        std::cout << "transformation_R_imu : " << transformation_R_imu[0] << " "
                                               << transformation_R_imu[1] << " "
                                               << transformation_R_imu[2] << std::endl;

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;

        out_transformation_R_imu[0] = transformation_R_imu[0];
        out_transformation_R_imu[1] = transformation_R_imu[1];
        out_transformation_R_imu[2] = transformation_R_imu[2];

        out_tramsformation_x_gps[0] = tramsformation_x_gps[0];
        out_tramsformation_x_gps[1] = tramsformation_x_gps[1];
        out_tramsformation_x_gps[2] = tramsformation_x_gps[2];

        out_transformation_br = transformation_br;
        out_transformation_fr = transformation_fr;
      }

      // Update camera poses with refined data
      if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
      {
        for (Poses::iterator itPose = sfm_data.poses.begin();
          itPose != sfm_data.poses.end(); ++itPose)
        {
          const IndexT indexPose = itPose->first;
          IndexT camId = indexPose;

          //single
          if (hashMap.find(camId)->second == 0)
          {
              Mat3 R_refined;
              ceres::AngleAxisToRotationMatrix(&single_poses[indexPose][0], R_refined.data());
              Vec3 t_refined(single_poses[indexPose][3], single_poses[indexPose][4], single_poses[indexPose][5]);
              // Update the pose
              Pose3 & pose = itPose->second;
              pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
          }else{  //unity
              Mat3 Rr_trans;
              Vec3 tr_trans;

              //std::vector<double> Rr_refined = reference_poses[indexPose];
              if(camId < sfm_data.views.size()/3)  //back
              {
                  camId = camId + sfm_data.views.size()/3;
                  ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rr_trans.data());

                  tr_trans(0) = transformation_br[3];
                  tr_trans(1) = transformation_br[4];
                  tr_trans(2) = transformation_br[5];

              }else if(camId >= sfm_data.views.size()/3*2)  //front
              {
                  camId = camId - sfm_data.views.size()/3;
                  ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rr_trans.data());

                  tr_trans(0) = transformation_fr[3];
                  tr_trans(1) = transformation_fr[4];
                  tr_trans(2) = transformation_fr[5];

              }else
              {
                  ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rr_trans.data());

                  tr_trans(0) = transformation_dr[3];
                  tr_trans(1) = transformation_dr[4];
                  tr_trans(2) = transformation_dr[5];

              }

              Vec3 gps_trans;
              gps_trans(0) = tramsformation_x_gps[0];
              gps_trans(1) = tramsformation_x_gps[1];
              gps_trans(2) = tramsformation_x_gps[2];
              Vec3 C_gps_res;
              C_gps_res(0) = C_gps[camId][0];
              C_gps_res(1) = C_gps[camId][1];
              C_gps_res(2) = C_gps[camId][2];
              Mat3 R_imu_trans;
              ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
              Mat3 R_imu_mat;
              ceres::AngleAxisToRotationMatrix(&R_imu[camId][0], R_imu_mat.data());

              Mat3 R_res;
              R_res = Rr_trans * R_imu_trans * R_imu_mat;
              Vec3 t_res;
              t_res = Rr_trans * gps_trans - Rr_trans * R_imu_trans * R_imu_mat * C_gps_res + tr_trans;


              //use svd





              Pose3 & pose = itPose->second;
              pose = Pose3(R_res, -R_res.transpose() * t_res);
            }
          }

      }


      // Update camera intrinsics with refined data
      if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
      {
        for (Intrinsics::iterator itIntrinsic = sfm_data.intrinsics.begin();
          itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
        {
          const IndexT indexCam = itIntrinsic->first;

          const std::vector<double> & vec_params = map_intrinsics[indexCam];
          itIntrinsic->second.get()->updateFromParams(vec_params);
        }
      }
      // Structure is already updated directly if needed (no data wrapping)
      return true;
    }



    }


bool Bundle_Adjustment_Ceres::Adjust_threecamera_gps_cail_second_change_c4
(SfM_Data & sfm_data,     // the SfM scene to refine
 const Optimize_Options options,
 std::vector<double>& in_transformation_R_imu,
 std::vector<double>& in_tramsformation_x_gps,
 std::vector<double>& in_transformation_br,
 std::vector<double>& in_transformation_fr
)
{
    //----------
    // Add camera parameters
    // - intrinsics
    // - poses [R|t]

    // Create residuals for each observation in the bundle adjustment problem. The
    // parameters for cameras and points are added automatically.
    //----------

    ceres::Problem problem;

    //set parameters about gps imu
    Hash_Map<IndexT, std::vector<double> > R_imu;
    Hash_Map<IndexT, std::vector<double> > C_gps;
    std::vector<double> transformation_R_imu = in_transformation_R_imu;//{0.0,0.0,0.0};//{-0.0120417, 0.0201986, 0.0258046};//{0.0,0.0,0.0}; //{-0.02074, 0.0199567, 0.0271603}; //{-0.0177615, 0.0177593, 0.0196154};  //{0.0,0.0,0.0}; //{-0.000427019, 0.0125482, 0.0134619};//{-0.00594065, 0.0111603, 0.0170654}; //{0.0,0.0,0.0};//{-0.00986379, 0.00976848, 0.0253123};
    std::vector<double> tramsformation_x_gps = in_tramsformation_x_gps;//{0.0,0.0,0.0};//{-2.01425, 2.23307, -28.3916};

    double * parameter_transformation_R_imu = &transformation_R_imu[0];
    problem.AddParameterBlock(parameter_transformation_R_imu, 3);

    double * parameter_tramsformation_x_gps = &tramsformation_x_gps[0];
    problem.AddParameterBlock(parameter_tramsformation_x_gps, 3);


    //problem.SetParameterBlockConstant(parameter_transformation_R_imu);
    //problem.SetParameterBlockConstant(parameter_tramsformation_x_gps);


    //get gps and imu data
  {
        {
            //read file and initialise external
            std::vector<double> c(3);
            std::string tFilePath = sfm_data.s_root_path + "/latlon_c.res";  //!!!
            std::ifstream in;
            in.open(tFilePath);
            if(!in)
            {
                std::cout << "open tFilePath file error!" << std::endl;
                return false;
            }
            std::size_t line = 0;
            const int count = sfm_data.views.size()/3;
            while((!in.eof()) && (line < count))
            {
                 in >> c[0] >> c[1] >> c[2];
                 C_gps[line + count] = {c[0], c[1], c[2]};
                 double * parameter_block = &C_gps[line + count][0];
                 problem.AddParameterBlock(parameter_block, C_gps[line + count].size());
                 problem.SetParameterBlockConstant(parameter_block);  ///needn't change
                ++line;
            }
            in.close();

        }
        //read imu
        {
        std::ifstream in;
        in.open(sfm_data.s_root_path + "/IMU.res"); //imu file

        int line = 0;
        const int count = sfm_data.views.size()/3;

        std::ofstream out_show_camera;
        out_show_camera.open(sfm_data.s_root_path + "/show_imu_R.txt");

        while((!in.eof()) && (line < count))
        {
            double roll, pitch, yaw;
            in >> roll >> pitch >> yaw;  //z x y

            Mat3 Rz, Ry, Rx, Rx_, Ry_;
            Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
                  sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
                  0.0, 0.0, 1.0;
            Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
                  0.0, 1.0, 0.0,
                  -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
            Rx << 1.0, 0.0, 0.0,
                  0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
                  0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);

            Mat3 R_d =  Rz*Rx*Ry;

            {

//                Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
//                Vec3 c = {C_gps[line + count][0], C_gps[line + count][1], C_gps[line + count][2]};
//                axis_x = R_d.transpose() * (axis_x + R_d*c);
//                axis_y = R_d.transpose() * (axis_y + R_d*c);
//                axis_z = R_d.transpose() * (axis_z + R_d*c);
//                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                       << axis_x(0) << " " << axis_x(1) << " " << axis_x(2)<< ";" << std::endl;
//                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                       << axis_y(0) << " " << axis_y(1) << " " << axis_y(2)<< ";" << std::endl;
//                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                       << axis_z(0) << " " << axis_z(1) << " " << axis_z(2)<< ";" << std::endl;

//                out_show_camera << R_d << std::endl << std::endl;
            }

//            //change axis
//            Vec3 changeAxis_A;
//            changeAxis_A << M_PI*1/sqrt(2.0), M_PI*1/sqrt(2.0), 0.0;
            Mat3 changeAxis_M;
//            ceres::AngleAxisToRotationMatrix(changeAxis_A.data(), changeAxis_M.data());

            changeAxis_M << 1.0, 0.0, 0.0,
                            0.0, -1.0, 0.0,
                            0.0, 0.0, -1.0;
            R_d = changeAxis_M * R_d;


            {

                Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
                Vec3 c = {C_gps[line + count][0], C_gps[line + count][1], C_gps[line + count][2]};
                axis_x = R_d.transpose() * (axis_x + R_d*c);
                axis_y = R_d.transpose() * (axis_y + R_d*c);
                axis_z = R_d.transpose() * (axis_z + R_d*c);
                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
                       << axis_x(0) << " " << axis_x(1) << " " << axis_x(2)<< ";" << std::endl;
                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
                       << axis_y(0) << " " << axis_y(1) << " " << axis_y(2)<< ";" << std::endl;
                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
                       << axis_z(0) << " " << axis_z(1) << " " << axis_z(2)<< ";" << std::endl;

                out_show_camera << R_d << std::endl << std::endl;
            }





            Vec3 temp_rd;
            ceres::RotationMatrixToAngleAxis(R_d.data(), temp_rd.data());
            R_imu[line + count] = {temp_rd(0), temp_rd(1), temp_rd(2)};

            double * parameter_block = &R_imu[line + count][0];
            problem.AddParameterBlock(parameter_block, R_imu[line + count].size());

            ++line;
        }
        in.close();

        out_show_camera.close();
    }
  }


    //get hash map
    std::map<int, int> hashMap;
    filiter_ununited_pose(sfm_data.poses, sfm_data.views, hashMap);
    //Hash_Map<IndexT, std::vector<double> > reference_poses;
    Hash_Map<IndexT, std::vector<double> > single_poses;
    for (Poses::const_iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose)
    {
      const IndexT indexPose = itPose->first;

      double * parameter_block;
      if(hashMap.find(indexPose)->second == 0)
      {
          std::cout << "signal " << indexPose << std::endl;
          //set single poses
          const Pose3 & pose = itPose->second;
          const Mat3 R = pose.rotation();
          const Vec3 t = pose.translation();

          double angleAxis[3];
          ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
          // angleAxis + translation
          single_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

          parameter_block = &single_poses[indexPose][0];
          problem.AddParameterBlock(parameter_block, single_poses[indexPose].size());
      }else{
          //set unity pose
//          if( (indexPose < sfm_data.views.size()/3) || (indexPose >= sfm_data.views.size()/3*2) )
//          {
//              continue;
//          }

//          const Pose3 & pose = itPose->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 t = pose.translation();
//          Vec3 C = pose.center();

//          double angleAxis[3];
//          ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
//          // angleAxis + translation
//          reference_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(0), C(1), C(2)};

//          parameter_block = &reference_poses[indexPose][0];
//          problem.AddParameterBlock(parameter_block, reference_poses[indexPose].size());


          //test if imu right (!!!)
//          if( (indexPose < sfm_data.views.size()/3) || (indexPose >= sfm_data.views.size()/3*2) )
//          {
//              continue;
//          }

//          const Pose3 & pose = itPose->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 t = pose.translation();

//          Vec3 C = pose.center();
//          Vec3 temp_rd;
//          ceres::RotationMatrixToAngleAxis(R.data(), temp_rd.data());
//          R_imu[indexPose] = {temp_rd(0), temp_rd(1), temp_rd(2)};

//          double * parameter_block = &R_imu[indexPose][0];
//          problem.AddParameterBlock(parameter_block, R_imu[line + count].size());
//          problem.SetParameterBlockConstant(parameter_block);  ///needn't change

      }


      if(parameter_block == NULL)
      {
          continue;
      }

//      if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
//      {
//        // set the whole parameter block as constant for best performance
//        problem.SetParameterBlockConstant(parameter_block);
//      }
//      else  // Subset parametrization
//      {
//        std::vector<int> vec_constant_extrinsic;
//        // If we adjust only the translation, we must set ROTATION as constant
//        if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
//        {
//          // Subset rotation parametrization
//          vec_constant_extrinsic.push_back(0);
//          vec_constant_extrinsic.push_back(1);
//          vec_constant_extrinsic.push_back(2);
//        }
//        // If we adjust only the rotation, we must set TRANSLATION as constant
//        if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
//        {
//          // Subset translation parametrization
//          vec_constant_extrinsic.push_back(3);
//          vec_constant_extrinsic.push_back(4);
//          vec_constant_extrinsic.push_back(5);
//        }
//        if (!vec_constant_extrinsic.empty())
//        {
//          ceres::SubsetParameterization *subset_parameterization =
//            new ceres::SubsetParameterization(6, vec_constant_extrinsic);
//          problem.SetParameterization(parameter_block, subset_parameterization);
//        }
//      }

    }

///-----------------------------------------
///
///

    //test distortion
//    Mat3 r_tb, r_tf;
//    r_tb << 1.0, 0.0, 0.0,
//            0.0, 0.787926, -0.61577,
//            0.0, 0.61577, 0.787926;
//    r_tf << 1.0, 0.0, 0.0,
//            0.0, 0.78796,  0.615727,
//            0.0, -0.615727, 0.78796;
//    Vec3 RA_tb, RA_tf;
//    ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
//    ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

        //set reference
        double R_br_angleAxis[3]= {0.0, 0.0, 0.0};//{RA_tb(0), RA_tb(1), -RA_tb(2)};//{0.665224, 0.00448363, 0.00291885};  //{0.703355, 0.00263495, -0.00149232};//{0.0, 0.0, 0.0}; //{0.0212232, 0.0206259, 0.0420509};  //{0.0, 0.0, 0.0}; //{0.669661, 0.00354668, 0.00220508}; //  //{0.0, 0.0, 0.0}; //{0.662013, -0.00281273, 0.000638113 }; //{0.0, 0.0, 0.0}; //{0.671924, 0.00432461, 0.000528427};
        double R_fr_angleAxis[3]= {0.0, 0.0, 0.0};//{RA_tf(0), RA_tf(1), -RA_tf(2)};//{-0.644463, 0.0113345, -0.00156344}; //{-0.682501, 0.0147873, 0.00150509};//{0.0, 0.0, 0.0}; //{-0.0232707, -0.00288052, 0.00035795};  //{0.0, 0.0, 0.0}; //{-0.650191, 0.0101123, -0.000844533}; //  //{0.0, 0.0, 0.0}; //{-0.658446, 0.00532591, 0.0010695}; //{-0.656461, 0.00508743, -0.00299628};
        double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};
        Vec3 t_db, t_df, t_dd;
        //double S = 1470.531783220221;
        t_db << 0.0, 0.0, 0.0;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
        t_df << 0.0, 0.0, 0.0;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
        t_dd << 0.0, 0.0, 0.0;


        std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
        transformation_br = in_transformation_br;
        transformation_fr = in_transformation_fr;
        transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                             t_dd(0), t_dd(1), t_dd(2)};


    std::cout << std::endl << std::endl << "initial value: " << std::endl;
        std::cout << "transformation_br : " << transformation_br[0] << " "
                                            << transformation_br[1] << " "
                                            << transformation_br[2] << " "
                                            << transformation_br[3] << " "
                                            << transformation_br[4] << " "
                                            << transformation_br[5] << " " << std::endl;

        std::cout << "transformation_fr : " << transformation_fr[0] << " "
                                            << transformation_fr[1] << " "
                                            << transformation_fr[2] << " "
                                            << transformation_fr[3] << " "
                                            << transformation_fr[4] << " "
                                            << transformation_fr[5] << " " << std::endl;

        std::cout << "transformation_dr : " << transformation_dr[0] << " "
                                            << transformation_dr[1] << " "
                                            << transformation_dr[2] << " "
                                            << transformation_dr[3] << " "
                                            << transformation_dr[4] << " "
                                            << transformation_dr[5] << " " << std::endl;

        std::cout << "transformation_R_imu : " << transformation_R_imu[0] << " "
                                               << transformation_R_imu[1] << " "
                                               << transformation_R_imu[2] << std::endl;

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;


        double * parameter_transformation_br = &transformation_br[0];
        double * parameter_transformation_fr = &transformation_fr[0];
        double * parameter_transformation_dr = &transformation_dr[0];
        problem.AddParameterBlock(parameter_transformation_br, 6);
        problem.AddParameterBlock(parameter_transformation_fr, 6);
        problem.AddParameterBlock(parameter_transformation_dr, 6);

        problem.SetParameterBlockConstant(parameter_transformation_dr);  ///needn't change

        //problem.SetParameterBlockConstant(parameter_transformation_br);
        //problem.SetParameterBlockConstant(parameter_transformation_fr);


#ifndef CALIBRATION_THREE_CAMERAS
        //problem.SetParameterBlockConstant(parameter_transformation_br);
        //problem.SetParameterBlockConstant(parameter_transformation_fr);
#endif

        if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_transformation_br);
          problem.SetParameterBlockConstant(parameter_transformation_fr);
          problem.SetParameterBlockConstant(parameter_transformation_dr);
        }
        else  // Subset parametrization
        {
          std::vector<int> vec_constant_extrinsic_br, vec_constant_extrinsic_fr,vec_constant_extrinsic_dr;
          // If we adjust only the translation, we must set ROTATION as constant
          if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
          {
            // Subset rotation parametrization
            vec_constant_extrinsic_br.push_back(0);
            vec_constant_extrinsic_br.push_back(1);
            vec_constant_extrinsic_br.push_back(2);
            vec_constant_extrinsic_fr.push_back(0);
            vec_constant_extrinsic_fr.push_back(1);
            vec_constant_extrinsic_fr.push_back(2);
            vec_constant_extrinsic_dr.push_back(0);
            vec_constant_extrinsic_dr.push_back(1);
            vec_constant_extrinsic_dr.push_back(2);
          }
          // If we adjust only the rotation, we must set TRANSLATION as constant
          if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
          {
            // Subset translation parametrization
              vec_constant_extrinsic_br.push_back(3);
              vec_constant_extrinsic_br.push_back(4);
              vec_constant_extrinsic_br.push_back(5);
              vec_constant_extrinsic_fr.push_back(3);
              vec_constant_extrinsic_fr.push_back(4);
              vec_constant_extrinsic_fr.push_back(5);
              vec_constant_extrinsic_dr.push_back(3);
              vec_constant_extrinsic_dr.push_back(4);
              vec_constant_extrinsic_dr.push_back(5);
          }
          if (!vec_constant_extrinsic_br.empty())
          {
            ceres::SubsetParameterization *subset_parameterization_br =
              new ceres::SubsetParameterization(6, vec_constant_extrinsic_br);
            problem.SetParameterization(parameter_transformation_br, subset_parameterization_br);

            ceres::SubsetParameterization *subset_parameterization_fr =
              new ceres::SubsetParameterization(6, vec_constant_extrinsic_fr);
            problem.SetParameterization(parameter_transformation_fr, subset_parameterization_fr);

            ceres::SubsetParameterization *subset_parameterization_dr =
              new ceres::SubsetParameterization(6, vec_constant_extrinsic_dr);
            problem.SetParameterization(parameter_transformation_dr, subset_parameterization_dr);
          }
        }


    // Data wrapper for refinement:
    Hash_Map<IndexT, std::vector<double> > map_intrinsics;

   // Setup Intrinsics data & subparametrization
    for (Intrinsics::const_iterator itIntrinsic = sfm_data.intrinsics.begin();
      itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
    {
      const IndexT indexCam = itIntrinsic->first;

      if (isValid(itIntrinsic->second->getType()))
      {
        map_intrinsics[indexCam] = itIntrinsic->second->getParams();

        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            itIntrinsic->second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
      else
      {
        std::cerr << "Unsupported camera type." << std::endl;
      }
    }


    // Set a LossFunction to be less penalized by false measurements
    //  - set it to NULL if you don't want use a lossFunction.
    ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(Square(4.0));
    // TODO: make the LOSS function and the parameter an option

std::cout << "start." << std::endl;

    // For all visibility add reprojections errors:
    for (Landmarks::iterator iterTracks = sfm_data.structure.begin();
         iterTracks!= sfm_data.structure.end(); ++iterTracks)
    {
      const Observations & obs = iterTracks->second.obs;

      for (Observations::const_iterator itObs = obs.begin();
        itObs != obs.end(); ++itObs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(itObs->first).get();

        //decide camera pose through pose id
        double * parameter_transformation;
        IndexT poseId = view->id_pose;


        //single
        if (hashMap.find(poseId)->second == 0)
        {

            // Each Residual block takes a point and a camera as input and outputs a 2
            // dimensional residual. Internally, the cost function stores the observed
            // image location and compares the reprojection against the observation.
            ceres::CostFunction* cost_function =
              IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), itObs->second.x);

            if (cost_function)
              problem.AddResidualBlock(cost_function,
                p_LossFunction,
                &map_intrinsics[view->id_intrinsic][0],
                &single_poses[view->id_pose][0],
                iterTracks->second.X.data());

        }else{  //unity
            if(poseId < sfm_data.views.size()/3)  //back
            {
                poseId = poseId + sfm_data.views.size()/3;

//                  ceres::CostFunction* cost_function2 =
//                          new ceres::AutoDiffCostFunction
//                            <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
//                              new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

//                  if (cost_function2)
//                    problem.AddResidualBlock(cost_function2,
//                      p_LossFunction,
//                      &map_intrinsics[view->id_intrinsic][0],
//                      &R_imu[poseId][0],
//                      &C_gps[poseId][0],
//                      &transformation_R_imu[0],
//                      &tramsformation_x_gps[0],
//                      &transformation_br[0],
//                      iterTracks->second.X.data());  //


                  ceres::CostFunction* cost_function2 =
                          new ceres::AutoDiffCostFunction
                            <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 8, 3, 3, 3, 3, 6, 3>(
                              new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

                  if (cost_function2)
                    problem.AddResidualBlock(cost_function2,
                      p_LossFunction,
                      &map_intrinsics[view->id_intrinsic][0],
                      &R_imu[poseId][0],
                      &C_gps[poseId][0],
                      &transformation_R_imu[0],
                      &tramsformation_x_gps[0],
                      &transformation_br[0],
                      iterTracks->second.X.data());  //




            }else if(poseId >= sfm_data.views.size()/3*2)  //front
            {
                poseId = poseId - sfm_data.views.size()/3;

//                ceres::CostFunction* cost_function2 =
//                        new ceres::AutoDiffCostFunction
//                          <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
//                            new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

//                if (cost_function2)
//                  problem.AddResidualBlock(cost_function2,
//                    p_LossFunction,
//                    &map_intrinsics[view->id_intrinsic][0],
//                    &R_imu[poseId][0],
//                    &C_gps[poseId][0],
//                    &transformation_R_imu[0],
//                    &tramsformation_x_gps[0],
//                    &transformation_fr[0],
//                    iterTracks->second.X.data());  //


                ceres::CostFunction* cost_function2 =
                        new ceres::AutoDiffCostFunction
                          <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 8, 3, 3, 3, 3, 6, 3>(
                            new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

                if (cost_function2)
                  problem.AddResidualBlock(cost_function2,
                    p_LossFunction,
                    &map_intrinsics[view->id_intrinsic][0],
                    &R_imu[poseId][0],
                    &C_gps[poseId][0],
                    &transformation_R_imu[0],
                    &tramsformation_x_gps[0],
                    &transformation_fr[0],
                    iterTracks->second.X.data());  //


            }else{  //down

//                ceres::CostFunction* cost_function2 =
//                        new ceres::AutoDiffCostFunction
//                          <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
//                            new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

//                if (cost_function2)
//                  problem.AddResidualBlock(cost_function2,
//                    p_LossFunction,
//                    &map_intrinsics[view->id_intrinsic][0],
//                    &R_imu[poseId][0],
//                    &C_gps[poseId][0],
//                    &transformation_R_imu[0],
//                    &tramsformation_x_gps[0],
//                    &transformation_dr[0],
//                    iterTracks->second.X.data());  //



                ceres::CostFunction* cost_function2 =
                        new ceres::AutoDiffCostFunction
                          <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 8, 3, 3, 3, 3, 6, 3>(
                            new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

                if (cost_function2)
                  problem.AddResidualBlock(cost_function2,
                    p_LossFunction,
                    &map_intrinsics[view->id_intrinsic][0],
                    &R_imu[poseId][0],
                    &C_gps[poseId][0],
                    &transformation_R_imu[0],
                    &tramsformation_x_gps[0],
                    &transformation_dr[0],
                    iterTracks->second.X.data());  //


            }

        }

      }
            if (options.structure_opt == Structure_Parameter_Type::NONE)
              problem.SetParameterBlockConstant(iterTracks->second.X.data());
    }



    // Configure a BA engine and run it
    //  Make Ceres automatically detect the bundle structure.
    ceres::Solver::Options ceres_config_options;
    ceres_config_options.max_num_iterations = 500;
    ceres_config_options.preconditioner_type =
      static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
    ceres_config_options.linear_solver_type =
      static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
    ceres_config_options.sparse_linear_algebra_library_type =
      static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
    ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
    ceres_config_options.logging_type = ceres::SILENT;
    ceres_config_options.num_threads = ceres_options_.nb_threads_;
    ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
    ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;


    // Solve BA
    ceres::Solver::Summary summary;
    ceres::Solve(ceres_config_options, &problem, &summary);
    if (ceres_options_.bCeres_summary_)
      std::cout << summary.FullReport() << std::endl;


    // If no error, get back refined parameters
    if (!summary.IsSolutionUsable())
    {
      if (ceres_options_.bVerbose_)
        std::cout << "Bundle Adjustment failed." << std::endl;
      return false;
    }
    else // Solution is usable
    {
      if (ceres_options_.bVerbose_)
      {
        // Display statistics about the minimization
        std::cout << std::endl
          << "Bundle Adjustment statistics (approximated RMSE) Adjust_threeCameras:\n"
          << " #views: " << sfm_data.views.size() << "\n"
          << " #poses: " << sfm_data.poses.size() << "\n"
          << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
          << " #tracks: " << sfm_data.structure.size() << "\n"
          << " #residuals: " << summary.num_residuals << "\n"
          << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
          << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
          << " Time (s): " << summary.total_time_in_seconds << "\n"
          << std::endl;


        //show reslut
        std::cout << "show reslut : " << std::endl;
        std::cout << "transformation_br : " << transformation_br[0] << " "
                                            << transformation_br[1] << " "
                                            << transformation_br[2] << " "
                                            << transformation_br[3] << " "
                                            << transformation_br[4] << " "
                                            << transformation_br[5] << std::endl;

        std::cout << "transformation_fr : " << transformation_fr[0] << " "
                                            << transformation_fr[1] << " "
                                            << transformation_fr[2] << " "
                                            << transformation_fr[3] << " "
                                            << transformation_fr[4] << " "
                                            << transformation_fr[5] << std::endl;

        std::cout << "transformation_dr : " << transformation_dr[0] << " "
                                            << transformation_dr[1] << " "
                                            << transformation_dr[2] << " "
                                            << transformation_dr[3] << " "
                                            << transformation_dr[4] << " "
                                            << transformation_dr[5] << std::endl;

        std::cout << "transformation_R_imu : " << transformation_R_imu[0] << " "
                                               << transformation_R_imu[1] << " "
                                               << transformation_R_imu[2] << std::endl;

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;

        in_transformation_R_imu[0] = transformation_R_imu[0];
        in_transformation_R_imu[1] = transformation_R_imu[1];
        in_transformation_R_imu[2] = transformation_R_imu[2];

        in_tramsformation_x_gps[0] = tramsformation_x_gps[0];
        in_tramsformation_x_gps[1] = tramsformation_x_gps[1];
        in_tramsformation_x_gps[2] = tramsformation_x_gps[2];

        in_transformation_br = transformation_br;
        in_transformation_fr = transformation_fr;
      }

      // Update camera poses with refined data
      if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
      {
        for (Poses::iterator itPose = sfm_data.poses.begin();
          itPose != sfm_data.poses.end(); ++itPose)
        {
          const IndexT indexPose = itPose->first;
          IndexT camId = indexPose;

          //single
          if (hashMap.find(camId)->second == 0)
          {
              Mat3 R_refined;
              ceres::AngleAxisToRotationMatrix(&single_poses[indexPose][0], R_refined.data());
              Vec3 t_refined(single_poses[indexPose][3], single_poses[indexPose][4], single_poses[indexPose][5]);
              // Update the pose
              Pose3 & pose = itPose->second;
              pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
          }else{  //unity
              Mat3 Rr_trans;
              Vec3 tr_trans;

              //std::vector<double> Rr_refined = reference_poses[indexPose];
              if(camId < sfm_data.views.size()/3)  //back
              {
                  camId = camId + sfm_data.views.size()/3;
                  ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rr_trans.data());

                  tr_trans(0) = transformation_br[3];
                  tr_trans(1) = transformation_br[4];
                  tr_trans(2) = transformation_br[5];

              }else if(camId >= sfm_data.views.size()/3*2)  //front
              {
                  camId = camId - sfm_data.views.size()/3;
                  ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rr_trans.data());

                  tr_trans(0) = transformation_fr[3];
                  tr_trans(1) = transformation_fr[4];
                  tr_trans(2) = transformation_fr[5];

              }else
              {
                  ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rr_trans.data());

                  tr_trans(0) = transformation_dr[3];
                  tr_trans(1) = transformation_dr[4];
                  tr_trans(2) = transformation_dr[5];

              }

              Vec3 gps_trans;
              gps_trans(0) = tramsformation_x_gps[0];
              gps_trans(1) = tramsformation_x_gps[1];
              gps_trans(2) = tramsformation_x_gps[2];
              Vec3 C_gps_res;
              C_gps_res(0) = C_gps[camId][0];
              C_gps_res(1) = C_gps[camId][1];
              C_gps_res(2) = C_gps[camId][2];
              Mat3 R_imu_trans;
              ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
              Mat3 R_imu_mat;
              ceres::AngleAxisToRotationMatrix(&R_imu[camId][0], R_imu_mat.data());

              Mat3 R_res;
              R_res = Rr_trans * R_imu_trans * R_imu_mat;
              Vec3 t_res;
              t_res = Rr_trans * gps_trans - Rr_trans * R_imu_trans * R_imu_mat * C_gps_res + tr_trans;


              //use svd





              Pose3 & pose = itPose->second;
              pose = Pose3(R_res, -R_res.transpose() * t_res);
            }
          }

      }


      // Update camera intrinsics with refined data
      if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
      {
        for (Intrinsics::iterator itIntrinsic = sfm_data.intrinsics.begin();
          itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
        {
          const IndexT indexCam = itIntrinsic->first;

          const std::vector<double> & vec_params = map_intrinsics[indexCam];
          itIntrinsic->second.get()->updateFromParams(vec_params);
        }
      }
      // Structure is already updated directly if needed (no data wrapping)
      return true;
    }



    }


bool Bundle_Adjustment_Ceres::Adjust_threecamera_gps_cail_first_cail_new_c1
(SfM_Data & sfm_data,     // the SfM scene to refine
 const Optimize_Options options,
 PoseStatusMap &poseStatusMap,
 std::vector<double>& out_transformation_R_imu,
 std::vector<double>& out_tramsformation_x_gps,
 std::vector<double>& out_transformation_br,
 std::vector<double>& out_transformation_fr,
 Hash_Map<IndexT, std::vector<double> >& C_gps_Map,
 Hash_Map<IndexT, std::vector<double> >& R_imu_Map
)
{
    //----------
    // Add camera parameters
    // - intrinsics
    // - poses [R|t]

    // Create residuals for each observation in the bundle adjustment problem. The
    // parameters for cameras and points are added automatically.
    //----------

    ceres::Problem problem;

    //set parameters about gps imu
//    Hash_Map<IndexT, std::vector<double> > R_imu;
//    Hash_Map<IndexT, std::vector<double> > C_gps;
    std::vector<double> transformation_R_imu = out_transformation_R_imu;//{0.0,0.0,0.0};//{-0.02074, 0.0199567, 0.0271603}; //{-0.0177615, 0.0177593, 0.0196154};  //{0.0,0.0,0.0}; //{-0.000427019, 0.0125482, 0.0134619};//{-0.00594065, 0.0111603, 0.0170654}; //{0.0,0.0,0.0};//{-0.00986379, 0.00976848, 0.0253123};
    std::vector<double> tramsformation_x_gps = out_tramsformation_x_gps;//{0.0,0.0,0.0}; //{-0.82937, 0.462316, -2.59147}; //{0.0,0.0,0.0};//{-0.281802, 0.574862, 0.064961};

    double * parameter_transformation_R_imu = &transformation_R_imu[0];
    problem.AddParameterBlock(parameter_transformation_R_imu, 3);

    double * parameter_tramsformation_x_gps = &tramsformation_x_gps[0];
    problem.AddParameterBlock(parameter_tramsformation_x_gps, 3);


    //problem.SetParameterBlockConstant(parameter_transformation_R_imu);
//    problem.SetParameterBlockConstant(parameter_tramsformation_x_gps);


    // Setup Poses data & subparametrization
    for (Hash_Map<IndexT, std::vector<double> >::iterator itC= C_gps_Map.begin();
         itC != C_gps_Map.end(); ++itC)
    {
      const IndexT indexPose = itC->first;
      double * parameter_block_C = &C_gps_Map[indexPose][0];
      problem.AddParameterBlock(parameter_block_C, C_gps_Map[indexPose].size());
      problem.SetParameterBlockConstant(parameter_block_C);  ///needn't change

      double * parameter_block_R = &R_imu_Map[indexPose][0];
      problem.AddParameterBlock(parameter_block_R, R_imu_Map[indexPose].size());

      if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
      {
        // set the whole parameter block as constant for best performance
        problem.SetParameterBlockConstant(parameter_block_R);
        problem.SetParameterBlockConstant(parameter_block_C);

        problem.SetParameterBlockConstant(parameter_tramsformation_x_gps);
      }
      else  // Subset parametrization
      {
        // If we adjust only the translation, we must set ROTATION as constant
        if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
        {
            problem.SetParameterBlockConstant(parameter_block_R);
        }
        // If we adjust only the rotation, we must set TRANSLATION as constant
        if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
        {
            problem.SetParameterBlockConstant(parameter_block_C);

            problem.SetParameterBlockConstant(parameter_tramsformation_x_gps);
        }
      }
    }


///-----------------------------------------
///
///
    //set reference
    double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};
    Vec3 t_dd;
    t_dd << 0.0, 0.0, 0.0;

    std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
    transformation_br = out_transformation_br;
    transformation_fr = out_transformation_fr;
    transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                         t_dd(0), t_dd(1), t_dd(2)};


    //show initial value
    {
        std::cout << std::endl << std::endl << "initial value: " << std::endl;
        std::cout << "transformation_br : " << transformation_br[0] << " "
                                            << transformation_br[1] << " "
                                            << transformation_br[2] << " "
                                            << transformation_br[3] << " "
                                            << transformation_br[4] << " "
                                            << transformation_br[5] << " " << std::endl;

        std::cout << "transformation_fr : " << transformation_fr[0] << " "
                                            << transformation_fr[1] << " "
                                            << transformation_fr[2] << " "
                                            << transformation_fr[3] << " "
                                            << transformation_fr[4] << " "
                                            << transformation_fr[5] << " " << std::endl;

        std::cout << "transformation_dr : " << transformation_dr[0] << " "
                                            << transformation_dr[1] << " "
                                            << transformation_dr[2] << " "
                                            << transformation_dr[3] << " "
                                            << transformation_dr[4] << " "
                                            << transformation_dr[5] << " " << std::endl;

        std::cout << "transformation_R_imu : " << transformation_R_imu[0] << " "
                                               << transformation_R_imu[1] << " "
                                               << transformation_R_imu[2] << std::endl;

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;
    }

    double * parameter_transformation_br = &transformation_br[0];
    double * parameter_transformation_fr = &transformation_fr[0];
    double * parameter_transformation_dr = &transformation_dr[0];
    problem.AddParameterBlock(parameter_transformation_br, 6);
    problem.AddParameterBlock(parameter_transformation_fr, 6);
    problem.AddParameterBlock(parameter_transformation_dr, 6);

    problem.SetParameterBlockConstant(parameter_transformation_dr);  ///needn't change

    //problem.SetParameterBlockConstant(parameter_transformation_br);
    //problem.SetParameterBlockConstant(parameter_transformation_fr);


//#ifndef CALIBRATION_THREE_CAMERAS
//        //problem.SetParameterBlockConstant(parameter_transformation_br);
//        //problem.SetParameterBlockConstant(parameter_transformation_fr);
//#endif

    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_transformation_br);
      problem.SetParameterBlockConstant(parameter_transformation_fr);
      problem.SetParameterBlockConstant(parameter_transformation_dr);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic_br, vec_constant_extrinsic_fr,vec_constant_extrinsic_dr;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic_br.push_back(0);
        vec_constant_extrinsic_br.push_back(1);
        vec_constant_extrinsic_br.push_back(2);
        vec_constant_extrinsic_fr.push_back(0);
        vec_constant_extrinsic_fr.push_back(1);
        vec_constant_extrinsic_fr.push_back(2);
        vec_constant_extrinsic_dr.push_back(0);
        vec_constant_extrinsic_dr.push_back(1);
        vec_constant_extrinsic_dr.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic_br.push_back(3);
        vec_constant_extrinsic_br.push_back(4);
        vec_constant_extrinsic_br.push_back(5);
        vec_constant_extrinsic_fr.push_back(3);
        vec_constant_extrinsic_fr.push_back(4);
        vec_constant_extrinsic_fr.push_back(5);
        vec_constant_extrinsic_dr.push_back(3);
        vec_constant_extrinsic_dr.push_back(4);
        vec_constant_extrinsic_dr.push_back(5);
      }
      if (!vec_constant_extrinsic_br.empty())
      {
        ceres::SubsetParameterization *subset_parameterization_br =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic_br);
          problem.SetParameterization(parameter_transformation_br, subset_parameterization_br);

        ceres::SubsetParameterization *subset_parameterization_fr =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic_fr);
          problem.SetParameterization(parameter_transformation_fr, subset_parameterization_fr);

        ceres::SubsetParameterization *subset_parameterization_dr =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic_dr);
          problem.SetParameterization(parameter_transformation_dr, subset_parameterization_dr);
       }
    }


    // Data wrapper for refinement:
    Hash_Map<IndexT, std::vector<double> > map_intrinsics;

   // Setup Intrinsics data & subparametrization
    for (Intrinsics::const_iterator itIntrinsic = sfm_data.intrinsics.begin();
      itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
    {
      const IndexT indexCam = itIntrinsic->first;

      if (isValid(itIntrinsic->second->getType()))
      {
        map_intrinsics[indexCam] = itIntrinsic->second->getParams();

        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
        {
          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);
        }
        else
        {
          const std::vector<int> vec_constant_intrinsic =
            itIntrinsic->second->subsetParameterization(options.intrinsics_opt);
          if (!vec_constant_intrinsic.empty())
          {
            ceres::SubsetParameterization *subset_parameterization =
              new ceres::SubsetParameterization(
                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
            problem.SetParameterization(parameter_block, subset_parameterization);
          }
        }
      }
      else
      {
        std::cerr << "Unsupported camera type." << std::endl;
      }
    }


    // Set a LossFunction to be less penalized by false measurements
    //  - set it to NULL if you don't want use a lossFunction.
    ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(Square(4.0));
    // TODO: make the LOSS function and the parameter an option

std::cout << "start." << std::endl;



    int intrinsicNum = map_intrinsics.at(0).size();
    // For all visibility add reprojections errors:
    for (Landmarks::iterator iterTracks = sfm_data.structure.begin();
         iterTracks!= sfm_data.structure.end(); ++iterTracks)
    {
      const Observations & obs = iterTracks->second.obs;

      for (Observations::const_iterator itObs = obs.begin();
        itObs != obs.end(); ++itObs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(itObs->first).get();

        //decide camera pose through pose id
        IndexT poseId = view->id_pose;

//        //single
//        if (hashMap.find(poseId)->second == 0)
//        {

//            // Each Residual block takes a point and a camera as input and outputs a 2
//            // dimensional residual. Internally, the cost function stores the observed
//            // image location and compares the reprojection against the observation.
//            ceres::CostFunction* cost_function =
//              IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), itObs->second.x);

//            if (cost_function)
//              problem.AddResidualBlock(cost_function,
//                p_LossFunction,
//                &map_intrinsics[view->id_intrinsic][0],
//                &single_poses[view->id_pose][0],
//                iterTracks->second.X.data());

//        }else
        {  //unity
            ceres::CostFunction* cost_function;
            if(intrinsicNum == 3)
            {
                cost_function =
                        new ceres::AutoDiffCostFunction
                          <ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE, 2, 3, 3, 3, 3, 3, 6, 3>(
                            new ResidualErrorFunctor_Pinhole_Intrinsic_IMU_GPS_USE(itObs->second.x.data()));

            }else if(intrinsicNum == 8)
            {
                cost_function =
                        new ceres::AutoDiffCostFunction
                          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_ThreeCamera_IMU_GPS, 2, 8, 3, 3, 3, 3, 6, 3>(
                            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_ThreeCamera_IMU_GPS(itObs->second.x.data()));

            }else
            {
                std::cout << "intrinsic error !" << std::endl;
            }
            int imgId;
            if(poseId / 100000 == 1) //back
            {
                imgId = poseId + 100000;
                if (cost_function)
                  problem.AddResidualBlock(cost_function,
                    p_LossFunction,
                    &map_intrinsics[view->id_intrinsic][0],
                    &R_imu_Map[imgId][0],
                    &C_gps_Map[imgId][0],
                    &transformation_R_imu[0],
                    &tramsformation_x_gps[0],
                    &transformation_br[0],
                    iterTracks->second.X.data());  //
            }else if(poseId / 100000 == 2) //down
            {
                imgId = poseId;
                if (cost_function)
                  problem.AddResidualBlock(cost_function,
                    p_LossFunction,
                    &map_intrinsics[view->id_intrinsic][0],
                    &R_imu_Map[imgId][0],
                    &C_gps_Map[imgId][0],
                    &transformation_R_imu[0],
                    &tramsformation_x_gps[0],
                    &transformation_dr[0],
                    iterTracks->second.X.data());  //
            }else if(poseId / 100000 == 3) //front
            {
                imgId = poseId - 100000;
                if (cost_function)
                  problem.AddResidualBlock(cost_function,
                    p_LossFunction,
                    &map_intrinsics[view->id_intrinsic][0],
                    &R_imu_Map[imgId][0],
                    &C_gps_Map[imgId][0],
                    &transformation_R_imu[0],
                    &tramsformation_x_gps[0],
                    &transformation_fr[0],
                    iterTracks->second.X.data());  //
            }





        }

      }
            if (options.structure_opt == Structure_Parameter_Type::NONE)
              problem.SetParameterBlockConstant(iterTracks->second.X.data());
    }



    // Configure a BA engine and run it
    //  Make Ceres automatically detect the bundle structure.
    ceres::Solver::Options ceres_config_options;
    ceres_config_options.max_num_iterations = 500;
    ceres_config_options.preconditioner_type =
      static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
    ceres_config_options.linear_solver_type =
      static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
    ceres_config_options.sparse_linear_algebra_library_type =
      static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
    ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
    ceres_config_options.logging_type = ceres::SILENT;
    ceres_config_options.num_threads = ceres_options_.nb_threads_;
    ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
    ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;


    // Solve BA
    ceres::Solver::Summary summary;
    ceres::Solve(ceres_config_options, &problem, &summary);
    if (ceres_options_.bCeres_summary_)
      std::cout << summary.FullReport() << std::endl;


    // If no error, get back refined parameters
    if (!summary.IsSolutionUsable())
    {
      if (ceres_options_.bVerbose_)
        std::cout << "Bundle Adjustment failed." << std::endl;
      return false;
    }
    else // Solution is usable
    {
      if (ceres_options_.bVerbose_)
      {
        // Display statistics about the minimization
        std::cout << std::endl
          << "Bundle Adjustment statistics (approximated RMSE) Adjust_threeCameras:\n"
          << " #views: " << sfm_data.views.size() << "\n"
          << " #poses: " << sfm_data.poses.size() << "\n"
          << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
          << " #tracks: " << sfm_data.structure.size() << "\n"
          << " #residuals: " << summary.num_residuals << "\n"
          << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
          << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
          << " Time (s): " << summary.total_time_in_seconds << "\n"
          << std::endl;


        //show reslut
        std::cout << "show reslut : " << std::endl;
        std::cout << "transformation_br : " << transformation_br[0] << " "
                                            << transformation_br[1] << " "
                                            << transformation_br[2] << " "
                                            << transformation_br[3] << " "
                                            << transformation_br[4] << " "
                                            << transformation_br[5] << std::endl;

        std::cout << "transformation_fr : " << transformation_fr[0] << " "
                                            << transformation_fr[1] << " "
                                            << transformation_fr[2] << " "
                                            << transformation_fr[3] << " "
                                            << transformation_fr[4] << " "
                                            << transformation_fr[5] << std::endl;

        std::cout << "transformation_dr : " << transformation_dr[0] << " "
                                            << transformation_dr[1] << " "
                                            << transformation_dr[2] << " "
                                            << transformation_dr[3] << " "
                                            << transformation_dr[4] << " "
                                            << transformation_dr[5] << std::endl;

        std::cout << "transformation_R_imu : " << transformation_R_imu[0] << " "
                                               << transformation_R_imu[1] << " "
                                               << transformation_R_imu[2] << std::endl;

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;

        out_transformation_R_imu[0] = transformation_R_imu[0];
        out_transformation_R_imu[1] = transformation_R_imu[1];
        out_transformation_R_imu[2] = transformation_R_imu[2];

        out_tramsformation_x_gps[0] = tramsformation_x_gps[0];
        out_tramsformation_x_gps[1] = tramsformation_x_gps[1];
        out_tramsformation_x_gps[2] = tramsformation_x_gps[2];

        out_transformation_br = transformation_br;
        out_transformation_fr = transformation_fr;
      }

      // Update camera poses with refined data
      if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
      {
        for (Poses::iterator itPose = sfm_data.poses.begin();
          itPose != sfm_data.poses.end(); ++itPose)
        {
          const IndexT indexPose = itPose->first;
          IndexT camId = indexPose;

//          //single
//          if (hashMap.find(camId)->second == 0)
//          {
//              Mat3 R_refined;
//              ceres::AngleAxisToRotationMatrix(&single_poses[indexPose][0], R_refined.data());
//              Vec3 t_refined(single_poses[indexPose][3], single_poses[indexPose][4], single_poses[indexPose][5]);
//              // Update the pose
//              Pose3 & pose = itPose->second;
//              pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
//          }else
          {  //unity
              Mat3 Rr_trans;
              Vec3 tr_trans;

              if(camId / 100000 == 1) //back
              {
                  camId = camId + 100000;
                  ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rr_trans.data());

                  tr_trans(0) = transformation_br[3];
                  tr_trans(1) = transformation_br[4];
                  tr_trans(2) = transformation_br[5];

              }else if(camId / 100000 == 2) //down
              {
                  ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rr_trans.data());

                  tr_trans(0) = transformation_dr[3];
                  tr_trans(1) = transformation_dr[4];
                  tr_trans(2) = transformation_dr[5];

              }else if(camId / 100000 == 3) //front
              {
                  camId = camId - 100000;
                  ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rr_trans.data());

                  tr_trans(0) = transformation_fr[3];
                  tr_trans(1) = transformation_fr[4];
                  tr_trans(2) = transformation_fr[5];

              }

              Vec3 gps_trans;
              gps_trans(0) = tramsformation_x_gps[0];
              gps_trans(1) = tramsformation_x_gps[1];
              gps_trans(2) = tramsformation_x_gps[2];
              Vec3 C_gps_res;
              C_gps_res(0) = C_gps_Map[camId][0];
              C_gps_res(1) = C_gps_Map[camId][1];
              C_gps_res(2) = C_gps_Map[camId][2];
              Mat3 R_imu_trans;
              ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
              Mat3 R_imu_mat;
              ceres::AngleAxisToRotationMatrix(&R_imu_Map[camId][0], R_imu_mat.data());

              Mat3 R_res;
              R_res = Rr_trans * R_imu_trans * R_imu_mat;
              Vec3 t_res;
              t_res = Rr_trans * gps_trans - Rr_trans * R_imu_trans * R_imu_mat * C_gps_res + tr_trans;

              Pose3 & pose = itPose->second;
              pose = Pose3(R_res, -R_res.transpose() * t_res);

            }
          }

      }


      // Update camera intrinsics with refined data
      if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
      {
        for (Intrinsics::iterator itIntrinsic = sfm_data.intrinsics.begin();
          itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
        {
          const IndexT indexCam = itIntrinsic->first;

          const std::vector<double> & vec_params = map_intrinsics[indexCam];
          itIntrinsic->second.get()->updateFromParams(vec_params);
        }
      }
      // Structure is already updated directly if needed (no data wrapping)
      return true;
    }



    }


} // namespace sfm
} // namespace openMVG
