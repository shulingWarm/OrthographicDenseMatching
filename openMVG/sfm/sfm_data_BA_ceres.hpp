// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_DATA_BA_CERES_HPP
#define OPENMVG_SFM_SFM_DATA_BA_CERES_HPP

#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/sfm/sfm_data_BA.hpp"

#include "openMVG/numeric/numeric.h"
#include "ceres/types.h"
#include "ceres/cost_function.h"

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/types.hpp"

namespace ceres { class CostFunction; }
namespace openMVG { namespace cameras { struct IntrinsicBase; } }
namespace openMVG { namespace sfm { struct SfM_Data; } }

namespace openMVG {
namespace sfm {

struct PoseStatus
{
    geometry::Pose3 pose;
//    int imgId;
    int status;
    int bdf;
};
using PoseStatusMap = Hash_Map<IndexT, PoseStatus>;  //<imageId, pose and status>

/// Create the appropriate cost functor according the provided input camera intrinsic model
/// Can be residual cost functor can be weighetd if desired (default 0.0 means no weight).
ceres::CostFunction * IntrinsicsToCostFunction
(
  cameras::IntrinsicBase * intrinsic,
  const Vec2 & observation,
  const double weight = 0.0
);

class Bundle_Adjustment_Ceres : public Bundle_Adjustment
{
  public:
  struct BA_Ceres_options
  {
    bool bVerbose_;
    unsigned int nb_threads_;
    bool bCeres_summary_;
    int linear_solver_type_;
    int preconditioner_type_;
    int sparse_linear_algebra_library_type_;
    double parameter_tolerance_;
    bool bUse_loss_function_;

    BA_Ceres_options(const bool bVerbose = true, bool bmultithreaded = true);
  };
  private:
    BA_Ceres_options ceres_options_;

  public:
  explicit Bundle_Adjustment_Ceres
  (
    Bundle_Adjustment_Ceres::BA_Ceres_options options =
    std::move(BA_Ceres_options())
  );

  BA_Ceres_options & ceres_options();

  bool Adjust
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  )override;

  //赵志豪 20201230 用来优化Z坐标的函数
  bool adjustWaterZ(sfm::SfM_Data&,const Optimize_Options &);

  bool Adjust_water_xy
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options,
          bool getPriorDiff=false //获取当前点云光心位置和先验光心位置的差别
  );

  bool Adjust_water_z
  (// the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    double &water_plane,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_water_z_fixC1
  (// the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    double &water_plane,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_water_min_z
  (// the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    double &water_plane,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_water_z_n
  (// the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    double &water_plane,
    double &ref_N,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_water_c1_D
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_water_c1_S
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_water_c4_D
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_water_c1_pi
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_water_c4_S
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_step23
  (// the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    const std::map<int, bool> &fixedImgIdList,
    const std::map<int, bool> &fixedXIdList,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_water_prepare
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_BAStep
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_water_c1
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_water2_c1
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_water3_c1
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_water5_c1
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_water6_c1
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_water4_c1
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_water4_c1_2
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_water4_c4
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_6
  (
    const std::map<int, bool> & controlImgMap,
    const std::map<int, bool> & controlStructureMap,
    const std::map<int, bool> & computeImageIdMap,
    std::map<int, bool> & computeStrctureIdMap,
    Landmarks & newStructure,
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_6_water
  (
    const std::map<int, bool> & controlImgMap,
    const std::map<int, bool> & controlStructureMap,
    const std::map<int, bool> & computeImageIdMap,
    std::map<int, bool> & computeStrctureIdMap,
    Landmarks & newStructure,
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_fix_some
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options,
    std::map<int, Pose3>& fixIdList
  );


  bool Adjust_down
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_front
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_back
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_down_front
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_fix2
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_old
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  );

  bool Adjust_single_fix_fuv
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options options
  );

  bool Adjust_down_gps_cail_first_change //c1
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options options,
    std::vector<double>& out_transformation_R_imu,
    std::vector<double>& out_tramsformation_x_gps
  );

  bool Adjust_down_gps_cail_second_change  //c1
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options options,
    std::vector<double>& in_transformation_R_imu,
    std::vector<double>& in_tramsformation_x_gps
  );

  bool Adjust_threecamera_gps_cail_first_change //c1
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options options,
    std::vector<double>& out_transformation_R_imu,
    std::vector<double>& out_tramsformation_x_gps,
    std::vector<double> &out_transformation_br,
    std::vector<double> &out_transformation_fr
  );

  bool Adjust_threecamera_gps_cail_second_change  //c1
  (// the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options options,
    std::vector<double>& in_transformation_R_imu,
    std::vector<double>& in_tramsformation_x_gps,
    std::vector<double> &out_transformation_br,
    std::vector<double> &out_transformation_fr);

  bool Adjust_threecamera_gps_cail_first_change_c4 //c4
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options options,
    std::vector<double>& out_transformation_R_imu,
    std::vector<double>& out_tramsformation_x_gps,
    std::vector<double> &out_transformation_br,
    std::vector<double> &out_transformation_fr
  );

  bool Adjust_threecamera_gps_cail_second_change_c4  //c4
  (// the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options options,
    std::vector<double>& in_transformation_R_imu,
    std::vector<double>& in_tramsformation_x_gps,
    std::vector<double> &out_transformation_br,
    std::vector<double> &out_transformation_fr);

  bool Adjust_threecamera_gps_cail_first_cail_new_c1 //c1
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options options,
    PoseStatusMap &poseStatusMap,
    std::vector<double>& out_transformation_R_imu,
    std::vector<double>& out_tramsformation_x_gps,
    std::vector<double> &out_transformation_br,
    std::vector<double> &out_transformation_fr,
    Hash_Map<IndexT, std::vector<double> >& C_gps_Map,
    Hash_Map<IndexT, std::vector<double> >& R_imu_Map
  );
private:
bool filiter_ununited_pose(sfm::Poses &poses, Views &views, std::map<int, int> &resMap);
bool filiter_ununited_pose_test(sfm::Poses &poses, Views &views, std::map<int, int> &resMap);

};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_DATA_BA_CERES_HPP
