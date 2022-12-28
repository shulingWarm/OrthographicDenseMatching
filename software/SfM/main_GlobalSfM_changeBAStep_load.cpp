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

#include "openMVG/sfm/sfm_data_filters.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <memory>
#include <string>

#include "ceres/rotation.h"

#include <boost/algorithm/string.hpp>

using namespace openMVG;
using namespace openMVG::sfm;

int main(int argc, char **argv)
{
  using namespace std;

  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sResData_Folder;
  std::string sOutDir = "";
  int iRotationAveragingMethod = int (ROTATION_AVERAGING_L2);
  int iTranslationAveragingMethod = int (TRANSLATION_AVERAGING_SOFTL1);
  bool b_use_motion_priors = false;
  std::string rtFilePath = "";
  bool bNeedSixteenth = true;
  int iBId = 0;
//  std::string sIntrinsic_refinement_options = "";
  std::string subPieceIdListStr = "1_2_3_4";

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('r', sResData_Folder, "resData_Folder") );
  cmd.add( make_option('o', sOutDir, "outdir") );
//  cmd.add( make_option('r', iRotationAveragingMethod, "rotationAveraging") );
//  cmd.add( make_option('t', iTranslationAveragingMethod, "translationAveraging") );
//  cmd.add( make_switch('P', "prior_usage") );
//  cmd.add( make_option('n', rtFilePath, "file path of R t file, n_gps_imu_321_x.res") );
//  cmd.add( make_option('e', bNeedSixteenth, "needSixteenth") );
//  cmd.add( make_option('b', iBId, "bId") );
  cmd.add( make_option('s', subPieceIdListStr, "subPieceIdListStr") );


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
    << "[-P|--prior_usage] Enable usage of motion priors (i.e GPS positions)\n"
    << "[-e|--needSixteenth] if need resize images into one-sixteenth"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
  {
    if (!stlplus::folder_create(sOutDir))
    {
      std::cerr << "Cannot create output directory" << std::endl;
      return EXIT_FAILURE;
    }
  }
  //
  // Load input SfM_Data scene
  SfM_Data sfm_data_all;
  if (!Load(sfm_data_all, sSfM_Data_Filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }


//  std::stringstream ss;
//  ss << iBId;
//  std::string sBId;
//  ss >> sBId;
  //
  /// Load datum
  ///
  std::vector<SfM_Data> sfm_data_list;
  std::vector<std::string> subPieceIdList;
  boost::split(subPieceIdList, subPieceIdListStr, boost::is_any_of("_"));
  //
  for(std::size_t bListId = 0; bListId < subPieceIdList.size(); ++bListId)
  {
      std::string sBId = subPieceIdList[bListId];
      std::cout << "this b Id " << sBId << std::endl;
      //
      /// load extrinsic
      //
      /// load extrinsics
      ///
      //
      std::cout << "loading extrinsics" << std::endl;
      std::ifstream in_extrinsics;
      std::string in_extrinsics_path = sResData_Folder + "/extrinsics_finished_" + sBId + ".res";
      in_extrinsics.open(in_extrinsics_path);
      if(!in_extrinsics)
      {
          std::cout << "open " << in_extrinsics_path << " file failed!" << std::endl;
          return EXIT_FAILURE;
      }
      int count_extrinsics;
      in_extrinsics >> count_extrinsics;
      int line_extrinsic = 0;
      while(line_extrinsic < count_extrinsics && !in_extrinsics.eof())
      {
          int eId;
          double q0, q1, q2, q3, t0, t1, t2;
          in_extrinsics >> eId >> q0 >> q1 >> q2 >> q3 >> t0 >> t1 >> t2;
//          std::vector<double> q;//{q0, q1, q2, q3};
//          q.push_back(q1);
//          q.push_back(q2);
//          q.push_back(q3);
//          q.push_back(q0);
          Vec3 t;
          t << t0, t1, t2;
//          Mat3 R;
          Eigen::Quaterniond q = Eigen::Quaterniond(q0, q1, q2, q3);
          Eigen::Matrix3d R = q.toRotationMatrix();
//          ceres::QuaternionToRotation(q.data(), R.data());
          Pose3 thePose(R, -R.inverse()*t);
          //
          sfm_data_all.poses[eId] = thePose;
          //
          ++line_extrinsic;
      }
      in_extrinsics.close();
      //
      /// load landmarks
      ///
      Landmarks::const_iterator itLEnd = sfm_data_all.structure.end();
      int landmarksStartId;
      if(sfm_data_all.structure.size() == 0)
      {
          landmarksStartId = 0;

      }else{
          --itLEnd;
          landmarksStartId = itLEnd->first+1;
      }
      //
      std::ifstream in_landmark;
      std::string in_landmark_path = sResData_Folder + "/landmarks_finished_" + sBId + ".res";
      in_landmark.open(in_landmark_path);
      if(!in_landmark)
      {
          std::cout << "open " << in_landmark_path << " file failed!" << std::endl;
          return EXIT_FAILURE;
      }
      int count_landmarks;
      in_landmark >> count_landmarks;
      int line_landmark = 0;
      while(line_landmark < count_landmarks && !in_landmark.eof())
      {
          int tempLandmarkId;
          double tempX_x, tempX_y, tempX_z;
          in_landmark >> tempLandmarkId >> tempX_x >> tempX_y >> tempX_z;
          Vec3 tempX(tempX_x, tempX_y, tempX_z);
          sfm_data_all.structure[tempLandmarkId + landmarksStartId].X = tempX;
          {
              int count_obs;
              in_landmark >> count_obs;
              int line_obs = 0;
              while(line_obs < count_obs && !in_landmark.eof())
              {
                  int tempImgId;
                  double tempFeat_x, tempFeat_y;
                  in_landmark >> tempImgId >> tempFeat_x >> tempFeat_y;
                  Vec2 tempFeat(tempFeat_x, tempFeat_y);
                  Observation obsTemp(tempFeat, 0);
                  sfm_data_all.structure[tempLandmarkId + landmarksStartId].obs[tempImgId] = obsTemp;
                  //
                  ++line_obs;
              }
          }
          //
          ++line_landmark;
      }
      in_landmark.close();
      //
      /// load new X
      // Xes
//      std::cout << "load result of Xes" << std::endl;
//      std::ifstream in_Xres;
//      std::string in_Xres_path = sSfM_Data_Filename + "/Xes_finished_" + sBId + ".res";
//      in_Xres.open(in_Xres_path);
//      if (!in_Xres)
//      {
//          std::cout << "create " << in_Xres_path << " file failed!" << std::endl;
//          return EXIT_FAILURE;
//      }
//      int Xcount;
//      int line = 0;
//      in_Xres >> Xcount;
//      while(line < Xcount && !in_Xres.eof())
//      {
//          int XId;
//          double X0, X1, X2;
//          in_Xres >> XId >> X0 >> X1 >> X2;
//          Vec3 Xtemp;
//          Xtemp << X0, X1, X2;
//          sfm_data_all.structure[XId + landmarksStartId].X = Xtemp;
//          //
//          ++line;
//      }
//      in_Xres.close();

  }
  //
  std::cout << "save datum" << std::endl;
  std::string newS = "";
  for(std::size_t bListId = 0; bListId < subPieceIdList.size(); ++bListId)
  {
      newS +=  subPieceIdList[bListId];
  }
  std::cout << newS;
//  Save(sfm_data_all,
//       stlplus::create_filespec(sOutDir, "sfm_data_init_"+newS+".bin" ).c_str(),
//       ESfM_Data(ALL));
  Save(sfm_data_all,
       stlplus::create_filespec(sOutDir, "sfm_data_finished_"+newS+".ply").c_str(),
       ESfM_Data(EXTRINSICS | STRUCTURE));
//  Save(sfm_data_all,
//       stlplus::create_filespec(sOutDir, "sfm_data_EXTRINSICS.ply").c_str(),
//       ESfM_Data(EXTRINSICS));
  Save(sfm_data_all,
       stlplus::create_filespec(sOutDir, "sfm_data_finished_"+newS+".json").c_str(),
       ESfM_Data(EXTRINSICS | STRUCTURE));

  if(subPieceIdList.size() > 1)
  {
      Save(sfm_data_all,
           stlplus::create_filespec(sOutDir, "sfm_data_init_"+newS+".ply").c_str(),
           ESfM_Data(EXTRINSICS | STRUCTURE));
      Save(sfm_data_all,
           stlplus::create_filespec(sOutDir, "sfm_data_init_"+newS+".json").c_str(),
           ESfM_Data(EXTRINSICS | STRUCTURE));


//      // Remove outliers (max_angle, residual error)
//      const size_t pointcount_initial = sfm_data_all.structure.size();
//      RemoveOutliers_PixelResidualError(sfm_data_all, 4.0);
//      const size_t pointcount_pixelresidual_filter = sfm_data_all.structure.size();
//      RemoveOutliers_AngleError(sfm_data_all, 2.0);
//      const size_t pointcount_angular_filter = sfm_data_all.structure.size();
//      std::cout << "Outlier removal (remaining #points):\n"
//        << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
//        << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
//        << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

      Save(sfm_data_all,
           stlplus::create_filespec(sOutDir, "sfm_data_delete_1234.ply").c_str(),
           ESfM_Data(EXTRINSICS | STRUCTURE));

      std::ofstream out_landmark;
      std::string out_landmark_path = sOutDir + "/landmarks_" + newS + ".res";
      out_landmark.open(out_landmark_path);
      if(!out_landmark)
      {
          std::cout << "create " << out_landmark_path << " file failed!" << std::endl;
      }
//      std::ofstream out_X;
//      std::string out_X_path = sfm_data_.s_root_path + "/X_all.res";
//      out_X.open(out_X_path);
//      if(!out_X)
//      {
//          std::cout << "create " << out_X_path << " file failed!" << std::endl;
//      }
      //
      int countL = sfm_data_all.structure.size();
      out_landmark << countL << std::endl;
      //
      for(Landmarks::const_iterator itL = sfm_data_all.structure.begin(); itL != sfm_data_all.structure.end(); ++itL)
      {
          int Xid = itL->first;
          Vec3 X = itL->second.X;
          out_landmark << Xid << " " << X(0) << " " << X(1) << " " << X(2) << std::endl;
          //
          int countObs = itL->second.obs.size();
          out_landmark << countObs << " ";
          for(Observations::const_iterator itObs = itL->second.obs.begin(); itObs != itL->second.obs.end(); ++itObs)
          {
              int imgId = itObs->first;
              Vec2 featId = itObs->second.x;
              out_landmark << imgId << " " << featId(0) << " " << featId(1) << " ";
          }
          out_landmark << std::endl;
      }
      out_landmark.close();
      //
      /// save extrinsics
      std::ofstream out_extrinsic;
      std::string out_extrinsic_path = sOutDir + "/extrinsics_" + newS + ".res";
      out_extrinsic.open(out_extrinsic_path);
      if(!out_extrinsic)
      {
          std::cout << "create " << out_extrinsic_path << " file failed!" << std::endl;
      }
      out_extrinsic << sfm_data_all.poses.size() << std::endl;
      for(Poses::const_iterator itPose = sfm_data_all.poses.begin(); itPose != sfm_data_all.poses.end(); ++itPose)
      {
          Vec3 t = itPose->second.translation();
          Mat3 R = itPose->second.rotation();
//          Vec4 q;
//          std::vector<double> q;
//          ceres::RotationMatrixToQuaternion(R.data(), q.data());
          Eigen::Quaterniond q(R);
          //
//          out_extrinsic << itPose->first << " " << q(0) << " "<< q(1) << " " << q(2) << " " << q(3) << " " << t(0) << " " << t(1) << " " << t(2) << std::endl;
          out_extrinsic << itPose->first << " " << q.w() << " "<< q.x() << " " << q.y() << " " << q.z() << " " << t(0) << " " << t(1) << " " << t(2) << std::endl;
          //
          Vec3 C = itPose->second.center();
      }
      out_extrinsic.close();
      //
      /// save intrinsics
      std::ofstream out_intrinsic;
      std::string out_intrinsic_path = sOutDir + "/intrinsics_" + newS + ".res";
      out_intrinsic.open(out_intrinsic_path);
      if(!out_intrinsic)
      {
          std::cout << "create " << out_intrinsic_path << " file failed!" << std::endl;
      }
      for(Intrinsics::const_iterator itIntrinsic = sfm_data_all.intrinsics.begin(); itIntrinsic != sfm_data_all.intrinsics.end(); ++itIntrinsic)
      {
          std::vector<double> thisDatum = itIntrinsic->second->getParams();
          out_intrinsic << thisDatum[0] << " " << thisDatum[1] << " " << thisDatum[2] << std::endl;
      }
      out_intrinsic.close();
  }


{
//  const cameras::Intrinsic_Parameter_Type intrinsic_refinement_options =
//    cameras::StringTo_Intrinsic_Parameter_Type("NONE");
//  if (intrinsic_refinement_options == static_cast<cameras::Intrinsic_Parameter_Type>(0) )
//  {
//    std::cerr << "Invalid input for Bundle Adjusment Intrinsic parameter refinement option" << std::endl;
//    return EXIT_FAILURE;
//  }

//  GlobalSfMReconstructionEngine_RelativeMotions sfmEngine(
//    sfm_data_all,
//    sOutDir,
//    stlplus::create_filespec(sOutDir, "Reconstruction_Report.html"));

//  // Configure the features_provider & the matches_provider
////  sfmEngine.SetFeaturesProvider(feats_provider.get());
////  sfmEngine.SetMatchesProvider(matches_provider.get());

//  // Configure reconstruction parameters
//  sfmEngine.Set_Intrinsics_Refinement_Type(intrinsic_refinement_options);
//  sfmEngine.Set_Use_Motion_Prior(false);

//  // Configure motion averaging method
//  sfmEngine.SetRotationAveragingMethod(
//    ERotationAveragingMethod(iRotationAveragingMethod));
//  sfmEngine.SetTranslationAveragingMethod(
//    ETranslationAveragingMethod(iTranslationAveragingMethod));

//  if (sfmEngine.Process_GCP_GPS_water2("")) //Process_threeCameras_imu_gps_init_new())//Process_test(sMatchesDir))//_down_imu_gps(sMatchesDir))
//  {
//      Save(sfm_data_all,
//           stlplus::create_filespec(sOutDir, "sfm_data_testAdjust.ply").c_str(),
//           ESfM_Data(EXTRINSICS | STRUCTURE));

//  }
}

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_all.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_all, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_all.structure.size();
  RemoveOutliers_AngleError(sfm_data_all, 2.0);
  const size_t pointcount_angular_filter = sfm_data_all.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  Save(sfm_data_all,
       stlplus::create_filespec(sOutDir, "sfm_data_delete.ply").c_str(),
       ESfM_Data(EXTRINSICS | STRUCTURE));




  return EXIT_SUCCESS;
}
