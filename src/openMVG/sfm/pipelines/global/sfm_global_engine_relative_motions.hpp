// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_GLOBAL_ENGINE_RELATIVE_MOTIONS_HPP
#define OPENMVG_SFM_GLOBAL_ENGINE_RELATIVE_MOTIONS_HPP

#include <memory>
#include <string>

#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"
#include "openMVG/sfm/pipelines/sfm_engine.hpp"

#include "openMVG/sfm/sfm_data_BA_ceres.hpp"

namespace htmlDocument { class htmlDocumentStream; }
namespace openMVG { namespace matching { struct PairWiseMatches; } }
namespace openMVG { namespace sfm { struct Features_Provider; } }
namespace openMVG { namespace sfm { struct Matches_Provider; } }
namespace openMVG { namespace sfm { struct SfM_Data; } }

namespace openMVG{
namespace sfm{

/// Global SfM Pipeline Reconstruction Engine.
/// - Method: Global Fusion of Relative Motions.
class GlobalSfMReconstructionEngine_RelativeMotions : public ReconstructionEngine
{
public:

  GlobalSfMReconstructionEngine_RelativeMotions(
    const SfM_Data & sfm_data,
    const std::string & soutDirectory,
    const std::string & loggingFile = "");

  ~GlobalSfMReconstructionEngine_RelativeMotions() override;

  void SetFeaturesProvider(Features_Provider * provider);
  void SetMatchesProvider(Matches_Provider * provider);

  void SetRotationAveragingMethod(ERotationAveragingMethod eRotationAveragingMethod);
  void SetTranslationAveragingMethod(ETranslationAveragingMethod eTranslation_averaging_method_);

  bool Process() override;
  bool Process_water_xy_z();
  bool Process_water_xy_z_synthetic();
  bool Process_water_xy_z_n();
  bool Process_GCP();
  bool Process_GCP_GPS(std::string rtFilePath);
  bool Process_GCP_GPS_BAStep_save(std::string rtFilePath, std::string sOutDir, int bId);
  bool Process_GCP_GPS_BAStep_load(std::string rtFilePath, std::string sOutDir, int bId);
  bool Process_GCP_GPS_water(std::string rtFilePath);
  bool Process_GCP_GPS_water2(std::string rtFilePath);
  bool Process_GCP_GPS_water3(std::string rtFilePath);
  bool Process_GCP_GPS_water4(std::string rtFilePath);
  bool Process_GCP_GPS_water5(std::string rtFilePath);
  bool Process_GCP_GPS_water6(std::string rtFilePath);
  bool Process_GCP_GPS_6(std::string rtFilePath);
  bool Process_GCP_GPS_4_margeTest(std::string rtFilePath);
  bool Process_GCP_step2(std::string outPath);
  bool Process_single_fix_fuv(std::string sMatchesDir);
  bool Process_down_imu_gps(std::string sMatchesDir);
  bool Process_threeCameras_imu_gps(std::string sMatchesDir);
  bool Process_threeCameras_imu_gps_c4(std::string sMatchesDir);
  bool Process_threeCameras_imu_gps_c4_changeFileForm(std::string sMatchesDir);
  bool Process_threeCameras_imu_gps_cail_new();
  bool Process_threeCameras_imu_gps_init_new(std::string rtFilePath);
  bool Process_threeCameras_change_step(std::string rtFilePath, std::string refSfmDataIdList);
  bool Process_threeCameras_for_large(std::string rtFilePath);
  bool Process_threeCameras_imu_gps_init_new_second(std::string rtFilePath, Mat3 changeR, Vec3 changet, double changeS);
  bool Process_threeCameras_imu_gps_init_new_test();
  bool Process_test(std::string MatchesDir);
  bool Process_step21();
  bool Process_step22();
  bool Process_step23(const std::map<int, bool> &fixedImgIdList, const std::map<int, bool> &fixedXIdList);
  bool Process_forSimulation();
  bool Process_forSimulation_afterOptim();
  bool Process_water6();

  bool Refine_initValue();
protected:
  /// Compute from relative rotations the global rotations of the camera poses
  bool Compute_Global_Rotations
  (
    const openMVG::rotation_averaging::RelativeRotations & vec_relatives_R,
    Hash_Map<IndexT, Mat3> & map_globalR
  );

  /// Compute/refine relative translations and compute global translations
  bool Compute_Global_Translations
  (
    const Hash_Map<IndexT, Mat3> & global_rotations,
    matching::PairWiseMatches & tripletWise_matches
  );
  bool Compute_Global_Translations_gps
  (
    const Hash_Map<IndexT, Mat3> & global_rotations,
    matching::PairWiseMatches & tripletWise_matches
  );

  /// Compute the initial structure of the scene
  bool Compute_Initial_Structure
  (
    matching::PairWiseMatches & tripletWise_matches
  );

  // Adjust the scene (& remove outliers)
  bool Adjust();
  bool Adjust_nonremove();
  bool Adjust_water_xy_z();
  bool Adjust_water_xy_z_n();
  bool Adjust_GCP();
  bool Adjust_GCP_step2(std::string outPath);
  bool Adjust_single_fix_fuv();
  bool Adjust_down_gps_change();
  bool Adjust_threecamera_gps_change();
  bool Adjust_threecamera_gps_change_c4();
  bool Adjust_threecamera_gps_cail_new(PoseStatusMap& poseStatusMap,
          const std::vector<double> &transformation_R_imu,
          const std::vector<double> &tramsformation_x_gps,
          const std::vector<double> &transformation_br,
          const std::vector<double> &transformation_fr,
          Hash_Map<IndexT, std::vector<double> > &C_gps_Map,
          Hash_Map<IndexT, std::vector<double> > &R_imu_Map);
  bool Adjust_init();
  bool Adjust_init_BAStep();
  bool Adjust_init_6();
  bool Adjust_init_4_margeTest();
  bool Adjust_init_water_c1();
  bool Adjust_init_water_c4();
  bool Adjust_init_water2_c1();
  bool Adjust_init_water3_c1();
  bool Adjust_init_water4_c1();
  bool Adjust_init_water5_c1();
  bool Adjust_init_water5_c4();
  bool Adjust_init_water6_c1();
  bool Adjust_init_water4_c1_2();
  bool Adjust_init_water4_c4();
  bool Adjust_change_step(std::map<int, Pose3>& fixIdList);
  bool Adjust_for_large();
  bool Adjust_init_fix2();
  bool Adjust_init_old();
  bool Adjust_step21();
  bool Adjust_step22();
  bool Adjust_step23(const std::map<int, bool> &fixedImgIdList, const std::map<int, bool> &fixedXIdList);
private:
  /// Compute relative rotations
  void Compute_Relative_Rotations
  (
    openMVG::rotation_averaging::RelativeRotations & vec_relatives_R
  );

  //----
  //-- Data
  //----

  // HTML logger
  std::shared_ptr<htmlDocument::htmlDocumentStream> html_doc_stream_;
  std::string sLogging_file_;

  // Parameter
  ERotationAveragingMethod eRotation_averaging_method_;
  ETranslationAveragingMethod eTranslation_averaging_method_;

  //-- Data provider
  Features_Provider  * features_provider_;
  Matches_Provider  * matches_provider_;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_GLOBAL_ENGINE_RELATIVE_MOTIONS_HPP
