// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/cameras/Cameras_Common_command_line_helper.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"
#include "openMVG/sfm/pipelines/global/sfm_global_engine_relative_motions.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_report.hpp"
#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include "openMVG/multiview/triangulation_nview.hpp"

#include <cstdlib>
#include <memory>
#include <string>

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <ceres/problem.h>
#include <ceres/solver.h>
#include <ceres/types.h>
#include <ceres/cost_function.h>

using namespace openMVG;
using namespace openMVG::sfm;



struct WaterPlaneCostFunction
{
  const std::vector<double> CI_Inf;
  const std::vector<double> CJ_Inf;

  WaterPlaneCostFunction
  (
    const std::vector<double> & _CI_Inf,
    const std::vector<double> & _CJ_Inf
  ): CI_Inf(_CI_Inf), CJ_Inf(_CJ_Inf)
  {
  }

  template <typename T> bool
  operator()
  (
    const T* const Za,
    T* residuals
  )
  const
  {
      const T rpI = T(CI_Inf[0]);
      const T ZcI = T(CI_Inf[1]);
      const T ZpI_ = T(CI_Inf[2]);
      const T lI_xy = T(CI_Inf[3]);
      const T lI_z = T(CI_Inf[4]);
      const T lI = lI_xy / lI_z;

      const T rpJ = T(CJ_Inf[0]);
      const T ZcJ = T(CJ_Inf[1]);
      const T ZpJ_ = T(CJ_Inf[2]);
      const T lJ_xy = T(CJ_Inf[3]);
      const T lJ_z = T(CJ_Inf[4]);
      const T lJ = lJ_xy / lJ_z;

      const T water_n = T(1.34);
      const T resI = ((water_n*water_n-T(1.0))*(rpI-(*Za-ZcI)*lI)*(rpI-(*Za-ZcI)*lI) + water_n*water_n*ZpI_*ZpI_);
      const T resJ = ((water_n*water_n-T(1.0))*(rpJ-(*Za-ZcJ)*lJ)*(rpJ-(*Za-ZcJ)*lJ) + water_n*water_n*ZpJ_*ZpJ_);

      residuals[0] = abs(resI-resJ);//*(resI-resJ);

    return true;
  }
};
///
///
///
///
///
///
///
struct WaterPlaneCostFunction2
{
  std::vector<double> CI_Inf;
  std::vector<double> CJ_Inf;
  Vec3 CI, CJ;

  WaterPlaneCostFunction2
  (
    const std::vector<double> & _CI_Inf,
    const std::vector<double> & _CJ_Inf,
    const Vec3 & _CI,
    const Vec3 & _CJ
  ): CI_Inf(_CI_Inf), CJ_Inf(_CJ_Inf), CI(_CI), CJ(_CJ)
  {
  }

  template <typename T> bool
  operator()
  (
      const T* const Za,
      const T* const P_xy,
      T* residuals
  )
  const
  {
      const T rpI = sqrt((P_xy[0]-CI(0))*(P_xy[0]-CI(0)) + (P_xy[1]-CI(1))*(P_xy[1]-CI(1)));
//      const T rpI = T(CI_Inf[0]);
      const T ZcI = T(CI_Inf[1]);
      const T ZpI_ = T(CI_Inf[2]);
      const T lI_xy = T(CI_Inf[3]);
      const T lI_z = T(CI_Inf[4]);
      const T lI = lI_xy / lI_z;
      //
      const T rpJ = sqrt((P_xy[0]-CJ(0))*(P_xy[0]-CJ(0)) + (P_xy[1]-CJ(1))*(P_xy[1]-CJ(1)));
//      const T rpJ = T(CJ_Inf[0]);
      const T ZcJ = T(CJ_Inf[1]);
      const T ZpJ_ = T(CJ_Inf[2]);
      const T lJ_xy = T(CJ_Inf[3]);
      const T lJ_z = T(CJ_Inf[4]);
      const T lJ = lJ_xy / lJ_z;

      const T water_n = T(1.34);
      const T resI = ((water_n*water_n-T(1.0))*(rpI-(*Za-ZcI)*lI)*(rpI-(*Za-ZcI)*lI) + water_n*water_n*ZpI_*ZpI_);
      const T resJ = ((water_n*water_n-T(1.0))*(rpJ-(*Za-ZcJ)*lJ)*(rpJ-(*Za-ZcJ)*lJ) + water_n*water_n*ZpJ_*ZpJ_);

//      residuals[0] = abs(resI-resJ);//*(resI-resJ);
//      residuals[0] = resI-resJ;//*(resI-resJ);
      residuals[0] = abs(resI-resJ);//*(resI-resJ);


    return true;
  }
};
//
//
//
//
//
//
//
//
struct WaterPlaneCostFunction3
{
  std::vector<double> CI_Inf;
  std::vector<double> CJ_Inf;
  Vec3 imgIFeatC, imgJFeatC;

  WaterPlaneCostFunction3
  (
    const std::vector<double> & _CI_Inf,
    const std::vector<double> & _CJ_Inf,
    const Vec3 & _imgIFeatC,
    const Vec3 & _imgJFeatC
  ): CI_Inf(_CI_Inf), CJ_Inf(_CJ_Inf), imgIFeatC(_imgIFeatC), imgJFeatC(_imgJFeatC)
  {
  }

  template <typename T> bool
  operator()
  (
      const T* const Za,
      const T* const P_xy,
      const T* const cam_extrinsics_I, // R_t
      const T* const cam_extrinsics_J, // R_t
      T* residuals
  )
  const
  {
      // CI PI
      const T * cam_R_I = &cam_extrinsics_I[0];
      const T * cam_t_I = &cam_extrinsics_I[3];
      const T cam_R_transpose_I[3] = {-cam_R_I[0], -cam_R_I[1], -cam_R_I[2]};
      T pose_center_I[3];
      ceres::AngleAxisRotatePoint(cam_R_transpose_I, cam_t_I, pose_center_I);
      pose_center_I[0] *= T(-1);
      pose_center_I[1] *= T(-1);
      pose_center_I[2] *= T(-1);
      const T temp_CL_I[3] = {T(imgIFeatC(0))-cam_t_I[0], T(imgIFeatC(1))-cam_t_I[1], T(imgIFeatC(2))-cam_t_I[2]};
      T PCI[3];
      ceres::AngleAxisRotatePoint(cam_R_transpose_I, temp_CL_I, PCI);
      // CJ PI
      const T * cam_R_J = &cam_extrinsics_J[0];
      const T * cam_t_J = &cam_extrinsics_J[3];
      const T cam_R_transpose_J[3] = {-cam_R_J[0], -cam_R_J[1], -cam_R_J[2]};
      T pose_center_J[3];
      ceres::AngleAxisRotatePoint(cam_R_transpose_J, cam_t_J, pose_center_J);
      pose_center_J[0] *= T(-1);
      pose_center_J[1] *= T(-1);
      pose_center_J[2] *= T(-1);
      const T temp_CL_J[3] = {T(imgJFeatC(0))-cam_t_J[0], T(imgJFeatC(1))-cam_t_J[1], T(imgJFeatC(2))-cam_t_J[2]};
      T PCJ[3];
      ceres::AngleAxisRotatePoint(cam_R_transpose_J, temp_CL_J, PCJ);
      //
      // others
      const T rpI = sqrt((PCI[0]-pose_center_I[0])*(PCI[0]-pose_center_I[0]) + (PCI[1]-pose_center_I[1])*(PCI[1]-pose_center_I[1]));
//      const T rpI = T(CI_Inf[0]);
      const T ZcI = T(CI_Inf[1]);
      const T ZpI_ = T(CI_Inf[2]);
      const T lI_xy = T(CI_Inf[3]);
      const T lI_z = T(CI_Inf[4]);
      const T lI = lI_xy / lI_z;
      //
      const T rpJ = sqrt((PCJ[0]-pose_center_I[0])*(PCJ[0]-pose_center_I[0]) + (PCJ[1]-pose_center_I[1])*(PCJ[1]-pose_center_I[1]));
//      const T rpJ = T(CJ_Inf[0]);
      const T ZcJ = T(CJ_Inf[1]);
      const T ZpJ_ = T(CJ_Inf[2]);
      const T lJ_xy = T(CJ_Inf[3]);
      const T lJ_z = T(CJ_Inf[4]);
      const T lJ = lJ_xy / lJ_z;

      const T water_n = T(1.34);
      const T resI = ((water_n*water_n-T(1.0))*(rpI-(*Za-ZcI)*lI)*(rpI-(*Za-ZcI)*lI) + water_n*water_n*ZpI_*ZpI_);
      const T resJ = ((water_n*water_n-T(1.0))*(rpJ-(*Za-ZcJ)*lJ)*(rpJ-(*Za-ZcJ)*lJ) + water_n*water_n*ZpJ_*ZpJ_);

//      residuals[0] = abs(resI-resJ);//*(resI-resJ);
//      residuals[0] = resI-resJ;//*(resI-resJ);
      residuals[0] = abs(resI-resJ);//*(resI-resJ);


    return true;
  }
};


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



int main(int argc, char **argv)
{
  using namespace std;
  std::cout << std::endl
//    << "-----------------------------------------------------------\n"
//    << "Global Structure from Motion:\n"
//    << "-----------------------------------------------------------\n"
//    << "Open Source implementation of the paper:\n"
//    << "\"Global Fusion of Relative Motions for "
//    << "Robust, Accurate and Scalable Structure from Motion.\"\n"
//    << "Pierre Moulon, Pascal Monasse and Renaud Marlet. "
//    << " ICCV 2013." << std::endl
//    << "------------------------------------------------------------"
    << std::endl;


  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sMatchesDir;
  std::string sOutDir = "";
  int iRotationAveragingMethod = int (ROTATION_AVERAGING_L2);
  int iTranslationAveragingMethod = int (TRANSLATION_AVERAGING_SOFTL1);
  std::string sIntrinsic_refinement_options = "ADJUST_ALL";
  bool b_use_motion_priors = false;
  std::string rtFilePath = "";
  bool bNeedSixteenth = false;

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('m', sMatchesDir, "matchdir") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('r', iRotationAveragingMethod, "rotationAveraging") );
  cmd.add( make_option('t', iTranslationAveragingMethod, "translationAveraging") );
  cmd.add( make_option('f', sIntrinsic_refinement_options, "refineIntrinsics") );
  cmd.add( make_switch('P', "prior_usage") );
  cmd.add( make_option('n', rtFilePath, "file path of R t file, n_gps_imu_321_x.res") );
  cmd.add( make_option('g', bNeedSixteenth, "if need resize images into one-sixteenth") );


  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch (const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_file] path to a SfM_Data scene\n"
    << "[-m|--matchdir] path to the matches that corresponds to the provided SfM_Data scene\n"
    << "[-o|--outdir] path where the output data will be stored\n"
    << "\n[Optional]\n"
    << "[-r|--rotationAveraging]\n"
      << "\t 1 -> L1 minimization\n"
      << "\t 2 -> L2 minimization (default)\n"
    << "[-t|--translationAveraging]:\n"
      << "\t 1 -> L1 minimization\n"
      << "\t 2 -> L2 minimization of sum of squared Chordal distances\n"
      << "\t 3 -> SoftL1 minimization (default)\n"
    << "[-f|--refineIntrinsics] Intrinsic parameters refinement option\n"
      << "\t ADJUST_ALL -> refine all existing parameters (default) \n"
      << "\t NONE -> intrinsic parameters are held as constant\n"
      << "\t ADJUST_FOCAL_LENGTH -> refine only the focal length\n"
      << "\t ADJUST_PRINCIPAL_POINT -> refine only the principal point position\n"
      << "\t ADJUST_DISTORTION -> refine only the distortion coefficient(s) (if any)\n"
      << "\t -> NOTE: options can be combined thanks to '|'\n"
      << "\t ADJUST_FOCAL_LENGTH|ADJUST_PRINCIPAL_POINT\n"
      <<      "\t\t-> refine the focal length & the principal point position\n"
      << "\t ADJUST_FOCAL_LENGTH|ADJUST_DISTORTION\n"
      <<      "\t\t-> refine the focal length & the distortion coefficient(s) (if any)\n"
      << "\t ADJUST_PRINCIPAL_POINT|ADJUST_DISTORTION\n"
      <<      "\t\t-> refine the principal point position & the distortion coefficient(s) (if any)\n"
    << "[-P|--prior_usage] Enable usage of motion priors (i.e GPS positions)\n"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  if (iRotationAveragingMethod < ROTATION_AVERAGING_L1 ||
      iRotationAveragingMethod > ROTATION_AVERAGING_L2 )  {
    std::cerr << "\n Rotation averaging method is invalid" << std::endl;
    return EXIT_FAILURE;
  }

  const cameras::Intrinsic_Parameter_Type intrinsic_refinement_options =
    cameras::StringTo_Intrinsic_Parameter_Type(sIntrinsic_refinement_options);
  if (intrinsic_refinement_options == static_cast<cameras::Intrinsic_Parameter_Type>(0) )
  {
    std::cerr << "Invalid input for Bundle Adjusment Intrinsic parameter refinement option" << std::endl;
    return EXIT_FAILURE;
  }

  if (iTranslationAveragingMethod < TRANSLATION_AVERAGING_L1 ||
      iTranslationAveragingMethod > TRANSLATION_AVERAGING_SOFTL1 )  {
    std::cerr << "\n Translation averaging method is invalid" << std::endl;
    return EXIT_FAILURE;
  }

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }


  Save(sfm_data,
    stlplus::create_filespec(sOutDir, "cloud_and_poses_testInit", ".ply"),
    ESfM_Data(ALL));

  // Init the regions_type from the image describer file (used for image regions extraction)
  using namespace openMVG::features;
  const std::string sImage_describer = stlplus::create_filespec(sMatchesDir, "image_describer", "json");
  std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
  if (!regions_type)
  {
    std::cerr << "Invalid: "
      << sImage_describer << " regions type file." << std::endl;
    return EXIT_FAILURE;
  }

  // Features reading
  std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
  if (!feats_provider->load(sfm_data, sMatchesDir, regions_type)) {
    std::cerr << std::endl
      << "Invalid features." << std::endl;
    return EXIT_FAILURE;
  }
  // Matches reading
  std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
  if // Try to read the two matches file formats
  (
    !(matches_provider->load(sfm_data, stlplus::create_filespec(sMatchesDir, "matches.e.txt")) ||
      matches_provider->load(sfm_data, stlplus::create_filespec(sMatchesDir, "matches.e.bin")))
  )
  {
    std::cerr << std::endl
      << "Invalid matches file." << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutDir))
  {
    if (!stlplus::folder_create(sOutDir))
    {
      std::cerr << "\nCannot create the output directory" << std::endl;
    }
  }

  //---------------------------------------
  // Global SfM reconstruction process
  //---------------------------------------

  openMVG::system::Timer timer;
  GlobalSfMReconstructionEngine_RelativeMotions sfmEngine(
    sfm_data,
    sOutDir,
    stlplus::create_filespec(sOutDir, "Reconstruction_Report.html"));

  // Configure the features_provider & the matches_provider
  sfmEngine.SetFeaturesProvider(feats_provider.get());
  sfmEngine.SetMatchesProvider(matches_provider.get());

  // Configure reconstruction parameters
  sfmEngine.Set_Intrinsics_Refinement_Type(intrinsic_refinement_options);
  b_use_motion_priors = cmd.used('P');
  sfmEngine.Set_Use_Motion_Prior(b_use_motion_priors);

  // Configure motion averaging method
  sfmEngine.SetRotationAveragingMethod(
    ERotationAveragingMethod(iRotationAveragingMethod));
  sfmEngine.SetTranslationAveragingMethod(
    ETranslationAveragingMethod(iTranslationAveragingMethod));

  if(bNeedSixteenth)
  {
      for (auto itView = sfm_data.views.begin(); itView != sfm_data.views.end(); ++ itView)
      {
        const IndexT imgId = itView->first;
        #pragma omp parallel for schedule(dynamic)
        for(size_t i = 0; i < feats_provider->feats_per_view[imgId].size(); ++i)
        {
            PointFeatures::iterator itPF = feats_provider->feats_per_view[imgId].begin();
            std::advance(itPF, i);
            (*itPF).x() *= 4.0;
            (*itPF).y() *= 4.0;
        }
      }
  }

  SfM_Data newdata;
  if (sfmEngine.Process())//_GCP_GPS(rtFilePath))//Process_threeCameras_imu_gps_init_new())//Process_test(sMatchesDir))//_down_imu_gps(sMatchesDir))
  {

      {
          newdata = sfmEngine.Get_SfM_Data();
          newdata.intrinsics = sfmEngine.Get_SfM_Data().intrinsics;
          Mat3 changeR;
          Vec3 changet;
          float s;
          {
              //run2
//              s = 141.783;
//              changeR << 13.886957168579, 139.231246948242, -22.894561767578,
//                         140.898696899414, -12.451700210571, 9.739796638489,
//                         7.553865909576, -23.705772399902, -139.58268737793;
//              changet << Vec3(89.912673950195, -147.658294677734, 163.292388916016);
//              s = 182.257;
//              changeR << 42.258850097656, 175.210388183594, -27.074182510376,
//                         177.236633300781, -42.432579040527, 2.03838300705,
//                         -4.343770980835, -26.801080703735, -180.223022460938;
//              changet << Vec3(106.549041748047, 63.339290618896, 170.729095458984);
//              s = 307.739;
//              changeR << -307.702819824219, 1.30886900425, -4.563062667847,
//                         2.248091220856, 300.714599609375, -65.339447021484,
//                         4.181000709534, -65.365013122559, -300.688385009766;
//              changet << Vec3(39.385131835938, 147.40007019043, 366.163726806641);

              s = 0.985698;
              changeR << 0.985211133957, 0.001476573758, 0.030937716365,
                         0.004451564513, 0.967604577541, -0.187941178679,
                         -0.030651364475, 0.187988087535, 0.967120110989;
              changet <<  -0.039469640702,  0.225849807262,  -0.146059274673;



              Mat3 changeZ;
              changeZ << -1, 0, 0,
                         0, 1, 0,
                         0, 0, -1;


              changeR = changeZ * changeR;
//              s = 519.0267046517480;
                //folder 16 Z = 0
//              changeR << -37.933506011963, -512.732482910156, 71.099311828613,
//                         -514.708435058594, 44.935962677002, 49.443969726562,
//                         -54.999961853027, -66.894119262695, -511.750793457031;
//              changet = Vec3(-414.110778808594, -43.797954559326, 543.712890625);
//              changeR << 1.0, 0.0, 0.0,
//                         0.0, 1.0, 0.0,
//                         0.0, 0.0, 1.0;

              changeR = changeR/s;
              changet = changet/s;

              //s = 0.98445;//0.30558845;//1.0;//3.1908;
          }

          for (const auto & pose_it : sfmEngine.Get_SfM_Data().poses)
          {
            const IndexT indexPose = pose_it.first;

            const Pose3 & pose = pose_it.second;
            const Mat3 R = pose.rotation();
            const Vec3 C = pose.center();

            Mat3 Rnew = R * changeR.inverse();
//            Vec3 Cnew = changet + Rnew.inverse() * R * C;
            Vec3 Cnew = changeR * s * C + changet;

            newdata.poses[indexPose] = Pose3(Rnew, Cnew);

          }

          newdata.structure = sfmEngine.Get_SfM_Data().structure;
          //set all X
          for(Landmarks::iterator itX = newdata.structure.begin();
              itX != newdata.structure.end(); ++itX)
          {
              const Vec3 & X = itX->second.X;

              Vec3 newX = changeR * s * X + changet;
//                  Vec3 newX = R_gc * X + t_gc;
//                  Vec3 newX = R0.inverse() * X;
              itX->second.X = newX;
          }

          //-- Export to disk computed scene (data & visualizable results)
          std::cout << "...Export SfM_Data to disk." << std::endl;

          Save(newdata,
            stlplus::create_filespec(sOutDir, "sfm_data2", ".bin"),
            ESfM_Data(ALL));

          Save(newdata,
            stlplus::create_filespec(sOutDir, "cloud_and_poses2", ".ply"),
            ESfM_Data(ALL));

      }
    std::cout << std::endl << " Total Ac-Global-Sfm took (s): " << timer.elapsed() << std::endl;

    std::cout << "...Generating SfM_Report.html" << std::endl;
    Generate_SfM_Report(sfmEngine.Get_SfM_Data(),
      stlplus::create_filespec(sOutDir, "SfMReconstruction_Report.html"));

    //-- Export to disk computed scene (data & visualizable results)
    std::cout << "...Export SfM_Data to disk." << std::endl;
    Save(sfmEngine.Get_SfM_Data(),
      stlplus::create_filespec(sOutDir, "sfm_data", ".bin"),
      ESfM_Data(ALL));

    Save(sfmEngine.Get_SfM_Data(),
      stlplus::create_filespec(sOutDir, "cloud_and_poses", ".ply"),
      ESfM_Data(ALL));

//    return EXIT_SUCCESS;
  }


  std::map<std::pair<int, int>, Vec3> P_Map;
  {
      //preparing P
      for(Landmarks::const_iterator itX = newdata.structure.begin();
          itX != newdata.structure.end(); ++itX)
      {
          Vec3 thisP = itX->second.X;
          for(Observations::const_iterator itObs = itX->second.obs.begin();
              itObs != itX->second.obs.end(); ++itObs)
          {
              std::pair<int, int> imgFeatPair;
              imgFeatPair.first = itObs->first;
              imgFeatPair.second = itObs->second.id_feat;
              //
              P_Map.insert(std::pair<std::pair<int, int>, Vec3>(imgFeatPair, thisP));

          }
      }

  }


  {
      double n = 1.34;
      double N = (n*n - 1.0) / (n*n);

      matching::PairWiseMatches matchesPairs = matches_provider->pairWise_matches_;

      //std::map<std::pair<>, >
//      newdata.structure;
      for(std::size_t matchPairId = 0; matchPairId < matchesPairs.size(); ++matchPairId)
      {
          //public std::map< Pair, IndMatches > PairWiseMatches的具体形式是这样
          //map里面的第1个pair存储的是两个包含对应点的图片
          matching::PairWiseMatches::const_iterator matchPairIter = matchesPairs.begin();
          std::advance(matchPairIter, matchPairId);//获取当前id位置对应的迭代器
          //取出包含对应点的两个图片的标号
          int imgI = matchPairIter->first.first;
          int imgJ = matchPairIter->first.second;
//          std::cout << "imgI : " << imgI << " imgJ : " << imgJ << std::endl;
          //
          /// get features
          std::vector<double> WPlanList;
          std::vector<double> WPlanList_;
          std::vector<double> WPlanList__;
          //
          std::vector<std::vector<double> > imgIInfoList, imgJInfoList;
          std::vector<Vec3> imgICList, imgJCList;
          std::vector<std::vector<double>> PxyList;
          std::vector<std::vector<double>> EIList, EJList;
          std::vector<Vec3> imgPCIList, imgPCJList;

          //遍历两个图片之间的所有对应点
          for(std::size_t matchPointId = 0; matchPointId < matchPairIter->second.size(); ++matchPointId)
          {
              //根据对应点的标号,取出对应点信息的迭代器
              matching::IndMatches::const_iterator matchPointIter = matchPairIter->second.begin();
              std::advance(matchPointIter, matchPointId);
              //两个对应点分别在特征容器中的索引标号
              int featureIID = matchPointIter->i_;
              int featureJID = matchPointIter->j_;
              //根据标号取出两个对应点的特征
              PointFeature fI = feats_provider->feats_per_view[imgI][featureIID];
              PointFeature fJ = feats_provider->feats_per_view[imgJ][featureJID];
              //
              // I
              std::shared_ptr<cameras::IntrinsicBase> itIntrinsicI = newdata.intrinsics[newdata.views[imgI]->id_intrinsic];
              std::vector<double> intrinsicI = itIntrinsicI->getParams();
              Pose3 poseI = newdata.poses[imgI];
              Mat3 rMatI = poseI.rotation();//旋转矩阵
              Vec3 tVecI = poseI.translation();//平移向量
              Vec3 cameraCenterI = poseI.center();//光芯坐标
              Vec3 reprojI = rMatI.transpose()*(Vec3(fI.x()-intrinsicI[1], fI.x()-intrinsicI[2], intrinsicI[0]) / intrinsicI[0] - tVecI);
              Vec2 lI = Vec2(sqrt((reprojI(0)-cameraCenterI(0))*(reprojI(0)-cameraCenterI(0)) + (reprojI(1)-cameraCenterI(1))*(reprojI(1)-cameraCenterI(1))), reprojI(2)-cameraCenterI(2));
              Vec3 CLI = reprojI-cameraCenterI;
              // J
              std::shared_ptr<cameras::IntrinsicBase> itIntrinsicJ = newdata.intrinsics[newdata.views[imgJ]->id_intrinsic];
              std::vector<double> intrinsicJ = itIntrinsicJ->getParams();
              Pose3 poseJ = newdata.poses[imgJ];
              Mat3 RJ = poseJ.rotation();
              Vec3 tJ = poseJ.translation();
              Vec3 CJ = poseJ.center();
              Vec3 reprojJ = RJ.transpose()*(Vec3(fJ.x()-intrinsicJ[1], fJ.x()-intrinsicJ[2], intrinsicJ[0]) / intrinsicJ[0] - tJ);
              Vec2 lJ = Vec2(sqrt((reprojJ(0)-CJ(0))*(reprojJ(0)-CJ(0)) + (reprojJ(1)-CJ(1))*(reprojJ(1)-CJ(1))), reprojJ(2)-CJ(2));
              Vec3 CLJ = reprojJ-CJ;
              //
              ///get P_xy
              ///
                Eigen::MatrixXd functionA;
                functionA.setZero(4,3);
                functionA << CLI(1), -CLI(0), 0.0,
                             CLI(2), 0.0, -CLI(0),
                             CLJ(1), -CLJ(0), 0.0,
                             CLJ(2), 0.0, -CLJ(0);
                Eigen::Vector4d function_b;
//                Eigen::VectorXd function_b;
                function_b << CLI(0)*cameraCenterI(1)-CLI(1)*cameraCenterI(0),
                              CLI(0)*cameraCenterI(2)-CLI(2)*cameraCenterI(0),
                              CLJ(0)*CJ(1)-CLJ(1)*CJ(0),
                              CLJ(0)*CJ(2)-CLJ(2)*CJ(0);
                //重复的计算过程,好像没用
                //functionA.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(-function_b);
                Eigen::Vector3d testP = functionA.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(-function_b);
                const Vec2 P_xy = Vec2(testP(0), testP(1));
              //
//                openMVG::Triangulation trian_obj;
//                trian_obj.add(
//                  itIntrinsicI->get_projective_equivalent(poseI),
//                  Vec2(fI.x(), fI.y()));
//                trian_obj.add(
//                  itIntrinsicJ->get_projective_equivalent(poseJ),
//                  Vec2(fJ.x(), fJ.y()));
//                //
//                Vec3 X = trian_obj.compute();
//                const Vec2 P_xy = Vec2(X(0), X(1));
              //
              std::map<std::pair<int, int>, Vec3>::const_iterator itPI = P_Map.find(std::pair<int, int>(imgI, featureIID));
              std::map<std::pair<int, int>, Vec3>::const_iterator itPJ = P_Map.find(std::pair<int, int>(imgJ, featureJID));
              if(itPI == P_Map.end() || itPJ == P_Map.end())
              {
                  continue;
              }
//              const Vec2 P_xy = Vec2(itPI->second(0), itPI->second(1));
              //
              if(itPI->second(0) < -0.366938 || itPI->second(0) > 0.154272 || itPI->second(1) < -1.183584 || itPI->second(1) > -0.388861)
//                if(itP->second(2) < -1.5 || itP->second(2) > -1.2)
              {
                  continue;

              }
              if(P_xy(0) < -0.366938 || P_xy(0) > 0.154272 || P_xy(1) < -1.183584 || P_xy(1) > -0.388861)
//                if(itP->second(2) < -1.5 || itP->second(2) > -1.2)
              {
                  continue;

              }

              lI(0) = sqrt((P_xy(0)-cameraCenterI(0))*(P_xy(0)-cameraCenterI(0)) + (P_xy(1)-cameraCenterI(1))*(P_xy(1)-cameraCenterI(1)));
              lJ(0) = sqrt((P_xy(0)-CJ(0))*(P_xy(0)-CJ(0)) + (P_xy(1)-CJ(1))*(P_xy(1)-CJ(1)));
              // I
              double rpI = sqrt((P_xy(0)-cameraCenterI(0))*(P_xy(0)-cameraCenterI(0)) + (P_xy(1)-cameraCenterI(1))*(P_xy(1)-cameraCenterI(1)));
              double ZpI_ = cameraCenterI(2) + lI(1)*rpI/lI(0);
              // J
              double rpJ = sqrt((P_xy(0)-CJ(0))*(P_xy(0)-CJ(0)) + (P_xy(1)-CJ(1))*(P_xy(1)-CJ(1)));
              double ZpJ_ = CJ(2) + lJ(1)*rpJ/lJ(0);
              //
              {

                  std::vector<double> imgJInfo;
                  imgJInfo.push_back(rpJ);
                  imgJInfo.push_back(CJ(2));
                  imgJInfo.push_back(ZpJ_);
                  imgJInfo.push_back(lJ(0));
                  imgJInfo.push_back(lJ(1));
                  imgJInfoList.push_back(imgJInfo);
                  std::vector<double> imgIInfo;
                  imgIInfo.push_back(rpI);
                  imgIInfo.push_back(cameraCenterI(2));
                  imgIInfo.push_back(ZpI_);
                  imgIInfo.push_back(lI(0));
                  imgIInfo.push_back(lI(1));
                  imgIInfoList.push_back(imgIInfo);

                  imgJCList.push_back(CJ);
                  imgICList.push_back(cameraCenterI);

                  std::vector<double> tempPxy;
                  tempPxy.push_back(P_xy(0));
                  tempPxy.push_back(P_xy(1));
                  PxyList.push_back(tempPxy);

                  double angleAxisI[3];
                  ceres::RotationMatrixToAngleAxis((const double*)rMatI.data(), angleAxisI);
                  std::vector<double> EIdata = {angleAxisI[0], angleAxisI[1], angleAxisI[2], tVecI(0), tVecI(1), tVecI(2)};
                  EIList.push_back(EIdata);
                  double angleAxisJ[3];
                  ceres::RotationMatrixToAngleAxis((const double*)RJ.data(), angleAxisJ);
                  std::vector<double> EJdata = {angleAxisJ[0], angleAxisJ[1], angleAxisJ[2], tVecI(0), tVecI(1), tVecI(2)};
                  EJList.push_back(EJdata);

                  imgPCIList.push_back((Vec3(fI.x()-intrinsicI[1], fI.x()-intrinsicI[2], intrinsicI[0]) / intrinsicI[0] - tVecI));
                  imgPCJList.push_back((Vec3(fJ.x()-intrinsicJ[1], fJ.x()-intrinsicJ[2], intrinsicJ[0]) / intrinsicJ[0] - tJ));


              }
              //
              // compute a b c
              double LI = lI(0)/lI(1);
              double LJ = lJ(0)/lJ(1);
              double lIZCI = LI*cameraCenterI(2);
              double lJZCJ = LJ*CJ(2);
              double a = N*(LI*LI - LJ*LJ);
              double b = 2.0*N*(LJ*(lJZCJ+rpJ) - LI*(lIZCI+rpI)) - 2.0*(ZpI_-ZpJ_);
              double c = N*(lIZCI*(lIZCI+2.0*rpI) - lJZCJ*(lJZCJ+2.0*rpJ) + (rpI+rpJ)*(rpI-rpJ)) + ZpI_*ZpI_ - ZpJ_*ZpJ_;
              //
              double b2_4ac = b*b - 4.0*a*c;
              if(b2_4ac < 0)
              {
                  std::cout << "<0 featI " << featureIID << " featJ " << featureJID << std::endl;
              }else{
                  double res1 = (-b + sqrt(b2_4ac)) / (2.0*a);
                  double res2 = (-b - sqrt(b2_4ac)) / (2.0*a);
//                    if(res1 > testP(2) && res2 > testP(2))
//                    if(res1 < testP(2) && res2 < testP(2))
//                    if(res1 < X(2) && res2 < X(2))
//                    if(res1 < testP(2) && res2 < testP(2))
                  double maxZp = ZpI_ > ZpJ_ ? ZpI_ : ZpJ_;
                  double minZp = ZpI_ > ZpJ_ ? ZpJ_ : ZpI_;

                  double tempRes, tempRes_;
                  if ((res1<maxZp)&&(res1>minZp))
                  {
                      tempRes = res2;
                      tempRes_ = res1;
//                        WPlanList_.push_back(res1);
//                        WPlanList.push_back(res2);
                  }else if ((res2<maxZp)&&(res2>minZp))
                  {
                      tempRes = res1;
                      tempRes_ = res2;
//                        WPlanList_.push_back(res2);
//                        WPlanList.push_back(res1);
                  }else{
                      continue;
                  }
                  if(tempRes>itPI->second(2))
                  {
                      WPlanList_.push_back(tempRes);
                      WPlanList.push_back(tempRes_);

                  }
              }

          }
          //
          ///
          ceres::Problem problem;
          double resZa = -1.3;
          if(imgJInfoList.size() < 20)
          {
              continue;
          }
          std::cout << "imgI : " << imgI << " imgJ : " << imgJ << std::endl;
          std::cout << "size : " << imgJInfoList.size() << std::endl;
          for (std::size_t idInfo = 0; idInfo < imgJInfoList.size(); ++idInfo)
          {
            {
//              ceres::CostFunction * cost_function =
//                new ceres::AutoDiffCostFunction<WaterPlaneCostFunction, 1, 1>(
//                  new WaterPlaneCostFunction(imgIInfoList[idInfo], imgJInfoList[idInfo]));
//              problem.AddResidualBlock(cost_function, nullptr, &resZa);
            }

              {
                  ceres::CostFunction * cost_function =
                    new ceres::AutoDiffCostFunction<WaterPlaneCostFunction2, 1, 1, 2>(
                      new WaterPlaneCostFunction2(imgIInfoList[idInfo], imgJInfoList[idInfo], imgICList[idInfo], imgJCList[idInfo] ));
                  problem.AddResidualBlock(cost_function, nullptr, &resZa, &PxyList[idInfo][0]);
              }

              {
//                  ceres::CostFunction * cost_function =
//                    new ceres::AutoDiffCostFunction<WaterPlaneCostFunction3, 1, 1, 2, 6, 6>(
//                      new WaterPlaneCostFunction3(imgIInfoList[idInfo], imgJInfoList[idInfo], imgPCIList[idInfo], imgPCJList[idInfo] ));
//                  problem.AddResidualBlock(cost_function, nullptr, &resZa, &PxyList[idInfo][0], &EIList[idInfo][0], &EJList[idInfo][0]);

              }
          }
          ceres::Solver::Options ceres_config_options;
          ceres_config_options.max_num_iterations = 500;
          ceres_config_options.preconditioner_type = ceres::IDENTITY;
//            static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
          ceres_config_options.linear_solver_type = ceres::DENSE_SCHUR;
//            static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
          ceres_config_options.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
//            static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
//          ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
          ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
          ceres_config_options.num_threads = 60;//ceres_options_.nb_threads_;
//          ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
          ceres_config_options.parameter_tolerance = 0.000000001;//ceres_options_.parameter_tolerance_;

          ceres::Solver::Summary summary;
          ceres::Solve(ceres_config_options, &problem, &summary);

          // If no error, get back refined parameters
//          if (!summary.IsSolutionUsable())
//          {
//            std::cout << "Bundle Adjustment failed." << std::endl;
//            return false;
//          }
//          else // Solution is usable
          {
              std::cout << "resZa " << resZa << std::endl;

          }



      }



  }




  return EXIT_FAILURE;
}



