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
#include <list>

#include <Eigen/Dense>

using namespace openMVG;
using namespace openMVG::sfm;

#define WATER_N 1.33

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
  int iRotationAveragingMethod = int (ROTATION_AVERAGING_L1);//(ROTATION_AVERAGING_L2);
  int iTranslationAveragingMethod = int (TRANSLATION_AVERAGING_SOFTL1);
  std::string sIntrinsic_refinement_options = "NONE";
  bool b_use_motion_priors = false;
  std::string rtFilePath = "";
  bool bNeedSixteenth = true;

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('m', sMatchesDir, "matchdir") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('r', iRotationAveragingMethod, "rotationAveraging") );
  cmd.add( make_option('t', iTranslationAveragingMethod, "translationAveraging") );
  cmd.add( make_option('f', sIntrinsic_refinement_options, "refineIntrinsics") );
  cmd.add( make_switch('P', "prior_usage") );
  cmd.add( make_option('n', rtFilePath, "file path of R t file, n_gps_imu_321_x.res") );
  cmd.add( make_option('e', bNeedSixteenth, "if need resize images into one-sixteenth") );


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


//  Save(sfm_data,
//    stlplus::create_filespec(sOutDir, "cloud_and_poses_testInit", ".ply"),
//    ESfM_Data(ALL));

  // Init the regions_type from the image describer file (used for image regions extraction)
//  using namespace openMVG::features;
//  const std::string sImage_describer = stlplus::create_filespec(sMatchesDir, "image_describer", "json");
//  std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
//  if (!regions_type)
//  {
//    std::cerr << "Invalid: "
//      << sImage_describer << " regions type file." << std::endl;
//    return EXIT_FAILURE;
//  }

  // Features reading
  std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
//  if (!feats_provider->load(sfm_data, sMatchesDir, regions_type)) {
//    std::cerr << std::endl
//      << "Invalid features." << std::endl;
//    return EXIT_FAILURE;
//  }
  // set feats matches extrinsics
  std::ifstream in_f, in_l, in_e;
  std::string f_path = sfm_data.s_root_path+"/x.txt";
  std::string l_path = sfm_data.s_root_path+"/landmarks_res.txt";
  std::string e_path = sfm_data.s_root_path+"/RC_res.txt";
//  in_f.open(f_path);
  in_l.open(l_path);
  in_e.open(e_path);
//  if(!in_f)
//  {
//      std::cout << "open file " << f_path << std::endl;
//      return EXIT_FAILURE;
//  }
  if(!in_l)
  {
      std::cout << "open file " << l_path << std::endl;
      return EXIT_FAILURE;
  }
  if(!in_e)
  {
      std::cout << "open file " << e_path << std::endl;
      return EXIT_FAILURE;
  }
//  //
  std::cout << "error 1 " << std::endl;
  const int imgSize = 9;//18;//18;
  int imgId = 0;
//  while(!in_f.eof() && imgId < imgSize)
//  {
//      openMVG::features::PointFeatures eachImg;
//      int thisImgId, featSize;
//      in_f >> thisImgId >> featSize;
//      for(std::size_t fId = 0; fId < featSize; ++fId)
//      {
//          double x, y;
//          in_f >> x >> y;
//          eachImg.push_back(openMVG::features::PointFeature(x,y));//[fId] = features::PointFeatures(x, y);
//      }
//      feats_provider->feats_per_view[thisImgId-1] = eachImg;
//      ++imgId;
//  }
//  in_f.close();
  imgId = 0;
  while(!in_e.eof() && imgId < imgSize)
  {
      double R[9], C[3];

      in_e >> R[0] >> R[1] >> R[2] >> C[0]
           >> R[3] >> R[4] >> R[5] >> C[1]
           >> R[6] >> R[7] >> R[8] >> C[2];
      Mat3 poseR;
      poseR << R[0], R[1], R[2],
               R[3], R[4], R[5],
               R[6], R[7], R[8];
      Vec3 poseC;
      poseC << C[0], C[1], C[2];
      Mat3 transR;
      transR << 0.0, 1.0, 0.0,
                1.0, 0.0, 0.0,
                0.0, 0.0, -1.0;
      Vec3 transC;
      transC << 0.0, 0.0, 1.0;
      //sfm_data.poses[imgId] = Pose3(poseR*transR, transR*poseC+transC);
      sfm_data.poses[imgId] = Pose3(poseR, poseC);
//      std::cout << "error 1" << std::endl;
////      sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(sfm_data.GetViews().at(imgId).get());
////      View& view_it = sfm_data.GetViews().begin();
//      View& view_it = sfm_data.GetViews();
//      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
////      sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(sfm_data.GetViews().at(imgId));
//      std::cout << "error 2" << std::endl;
////      prior->SetPoseCenterPrior(poseC, Vec3(1.0,1.0,1.0));
//      prior->center_weight_ = Vec3(1.0,1.0,1.0);
//      prior->pose_center_ = poseC;
//      std::cout << "error 3" << std::endl;
////      /sfm_data.views[imgId].SetPoseCenterPrior(poseC,Vec3(1.0,1.0,1.0));
      ++imgId;
  }
  in_e.close();

//  for()
//  {

//  }

  std::cout << "error 2 " << std::endl;
  //set landmarks
  Landmarks input_l;
  int l_count;
  in_l >> l_count;
  for (std::size_t lId = 0; lId < l_count; ++lId)
  {
      int obs_count;
      in_l >> obs_count;
      Landmark thisL;
      Observations thisObs;
      for (std::size_t oId = 0; oId < obs_count; ++oId)
      {
          int imgId, featId;
          float x1, x2;
          in_l >> imgId >> featId >> x1 >> x2;
          //
          Observation newObs = Observation(Vec2(x1,x2), featId-1);
          thisObs[imgId-1] = newObs;
      }
      float X1, X2, X3;
      in_l >> X1 >> X2 >> X3;
      //
      thisL.X = Vec3(X1, X2, X3);
      thisL.obs = thisObs;
      input_l.insert(std::pair<IndexT, Landmark>(lId,thisL));
  }
  in_l.close();
  sfm_data.structure = input_l;


  std::cout << "error 3" << std::endl;
//  //set match
//  matching::PairWiseMatches pairWise_matches_;
//  const int pairSize = 149;//153;//149;//
//  int pairId = 0;
//  while (!in_m.eof() && pairId < pairSize)
//  {

//      int imgIdI, imgIdJ, listSize;

//      in_m >> imgIdI >> imgIdJ >> listSize;
//      Pair imgPair = Pair(imgIdI-1, imgIdJ-1);
//      matching::IndMatches indM;
//      for (std::size_t fId = 0; fId < listSize; ++fId)
//      {
//          int fIdI, fIdJ;
//          in_m >> fIdI >> fIdJ;
//          indM.push_back(matching::IndMatch(fIdI-1, fIdJ-1));
//      }
//      pairWise_matches_.insert(std::pair<Pair, matching::IndMatches>(imgPair, indM));
//      //
//      ++pairId;
//  }
//  in_m.close();
//  std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
//  matches_provider->pairWise_matches_ = pairWise_matches_;


//  //set match
//  matching::PairWiseMatches pairWise_matches_;
//  for (std::size_t imgIdI = 0; imgIdI < imgSize-1; ++imgIdI)
//  {
//      for (std::size_t imgIdJ = imgIdI+1; imgIdJ < imgSize; ++imgIdJ)
//      {
//          Pair imgPair = Pair(imgIdI, imgIdJ);
//          matching::IndMatches indM;
//          //
//          for (std::size_t fId = 0; fId < featSize; ++fId)
//          {
//              indM.push_back(matching::IndMatch(fId, fId));
//          }
//          //
//          pairWise_matches_.insert(std::pair<Pair, matching::IndMatches>(imgPair, indM));
//      }
//  }
//  std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
//  matches_provider->pairWise_matches_ = pairWise_matches_;










  // Matches reading
//  std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
//  if // Try to read the two matches file formats
//  (
//    !(matches_provider->load(sfm_data, stlplus::create_filespec(sMatchesDir, "matches.e.txt")) ||
//      matches_provider->load(sfm_data, stlplus::create_filespec(sMatchesDir, "matches.e.bin")))
//  )
//  {
//    std::cerr << std::endl
//      << "Invalid matches file." << std::endl;
//    return EXIT_FAILURE;
//  }

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

  ///!!!
//  if(bNeedSixteenth)
//  {
//      for (auto itView = sfm_data.views.begin(); itView != sfm_data.views.end(); ++ itView)
//      {
//        const IndexT imgId = itView->first;
//        #pragma omp parallel for schedule(dynamic)
//        for(size_t i = 0; i < feats_provider->feats_per_view[imgId].size(); ++i)
//        {
//            openMVG::features::PointFeatures::iterator itPF = feats_provider->feats_per_view[imgId].begin();
//            std::advance(itPF, i);
//            (*itPF).x() *= 4.0;
//            (*itPF).y() *= 4.0;
//        }
//      }
//  }

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
  //sfmEngine.SetMatchesProvider(matches_provider.get());

  // Configure reconstruction parameters
  sfmEngine.Set_Intrinsics_Refinement_Type(intrinsic_refinement_options);
  b_use_motion_priors = false;//true;// cmd.used('P');
  sfmEngine.Set_Use_Motion_Prior(b_use_motion_priors);

  // Configure motion averaging method
  sfmEngine.SetRotationAveragingMethod(
    ERotationAveragingMethod(iRotationAveragingMethod));
  sfmEngine.SetTranslationAveragingMethod(
    ETranslationAveragingMethod(iTranslationAveragingMethod));

  std::cout << "error 4 " << std::endl;
  SfM_Data newdata, undiswaterdata;
  if (sfmEngine.Process_water_xy_z_synthetic())//_GCP_GPS_water(rtFilePath)) //Process_threeCameras_imu_gps_init_new())//Process_test(sMatchesDir))//_down_imu_gps(sMatchesDir))
  //if (sfmEngine.Process_water_xy_z())
  {
      newdata = sfmEngine.Get_SfM_Data();
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
      stlplus::create_filespec(sOutDir, "sfm_data", ".json"),
      ESfM_Data(ALL));

    Save(sfmEngine.Get_SfM_Data(),
      stlplus::create_filespec(sOutDir, "cloud_and_poses", ".ply"),
      ESfM_Data(ALL));

    {
//        SfM_Data newdata = sfmEngine.Get_SfM_Data();
        SfM_Data undiswaterdata = sfmEngine.Get_SfM_Data();
        std::ofstream out_e, out_k;
        std::string e_path = sfmEngine.Get_SfM_Data().s_root_path+"/RC_new_withWater.txt";//undisWater.txt";//
        std::string k_path = sfmEngine.Get_SfM_Data().s_root_path+"/K_new_withWater.txt";//undisWater.txt";//
        out_e.open(e_path);
        out_k.open(k_path);
        if(!out_e)
        {
            std::cout << "open file " << e_path << std::endl;
            return EXIT_FAILURE;
        }

        for(int imgId = 0; imgId < sfmEngine.Get_SfM_Data().poses.size(); ++imgId)
        {
            const Pose3 pose = sfmEngine.Get_SfM_Data().poses.at(imgId);
            Mat3 R = pose.rotation();
            Vec3 C = pose.center();

            out_e << R(0,0) << " " << R(0,1) << " " << R(0,2) << " " << C(0) << std::endl
                  << R(1,0) << " " << R(1,1) << " " << R(1,2) << " " << C(1) << std::endl
                  << R(2,0) << " " << R(2,1) << " " << R(2,2) << " " << C(2) << std::endl;

            std::vector<double> cam_intrinsics = sfm_data.intrinsics[sfm_data.views[imgId]->id_intrinsic]->getParams();
            out_k << cam_intrinsics[0] << " " << cam_intrinsics[1] << " " << cam_intrinsics[2] << std::endl;
        }
        out_e.close();
        out_k.close();

        //save for analyse error
        std::ofstream out_X, out_xm;
        out_X.open(sOutDir+"/X_withWater.txt");//undisWater.txt");//withWater.txt");//
        out_xm.open(sfmEngine.Get_SfM_Data().s_root_path+"/xm_withWater.txt");//undisWater.txt");//withWater.txt");
        for(Landmarks::const_iterator itX = undiswaterdata.structure.begin();
            itX != undiswaterdata.structure.end(); ++itX)
        {
            //for pot17
//            if(itX->second.X(0)<-0.706062 || itX->second.X(0) > 0.618697 || itX->second.X(1)<-0.691933 || itX->second.X(1)>0.29423 || itX->second.X(2) > 2.551696)
//            {

//                continue;
//            }

//            if(itX->second.X(2) > 1.114)
//            {
//                continue;
//            }

            out_xm << itX->second.obs.size() << " ";
            for(Observations::const_iterator itObs = itX->second.obs.begin();
                itObs != itX->second.obs.end(); ++itObs)
            {
                out_xm << itObs->first << " " << itObs->second.id_feat << " " << itObs->second.x(0) << " " << itObs->second.x(1) << " ";

            }
            out_xm << itX->second.X(0) << " "
                   << itX->second.X(1) << " "
                   << itX->second.X(2) << std::endl;
            out_X << itX->second.obs.begin()->first << " "
                  << itX->second.obs.begin()->second.id_feat << " "
                  << itX->second.X(0) << " "
                  << itX->second.X(1) << " "
                  << itX->second.X(2) << std::endl;

        }
        out_X.close();
        out_xm.close();
    }

    return EXIT_SUCCESS;
    //save for analyse error
//    std::ofstream out_X, out_xm;
//    out_X.open(sOutDir+"/X_withWater.txt");//undisWater.txt");//withWater.txt");//
//    out_xm.open(sOutDir+"/xm_withWater.txt");//undisWater.txt");//withWater.txt");
//    for(Landmarks::const_iterator itX = undiswaterdata.structure.begin();
//        itX != undiswaterdata.structure.end(); ++itX)
//    {
//        out_xm << itX->second.obs.size() << " ";
//        for(Observations::const_iterator itObs = itX->second.obs.begin();
//            itObs != itX->second.obs.end(); ++itObs)
//        {
//            out_xm << itObs->first << " " << itObs->second.id_feat << " " << itObs->second.x(0) << " " << itObs->second.x(1) << " ";

//        }
//        out_xm << itX->second.X(0) << " "
//               << itX->second.X(1) << " "
//               << itX->second.X(2) << std::endl;
//        out_X << itX->second.obs.begin()->first << " "
//              << itX->second.obs.begin()->second.id_feat << " "
//              << itX->second.X(0) << " "
//              << itX->second.X(1) << " "
//              << itX->second.X(2) << std::endl;

//    }
//    out_X.close();
//    out_xm.close();
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

  // underwater scene
  {

      const double n = 1.2796;
      const double N = (n*n - 1.0) / (n*n);
      const double waterPlane = -0.5075;//0.25;
      for(Landmarks::const_iterator itX = newdata.structure.begin();
          itX != newdata.structure.end(); ++itX)
      {
          Vec3 thisP = itX->second.X;
          if(thisP(2)>waterPlane)
          {
              continue;
          }

//          for(Observations::const_iterator itObs = itX->second.obs.begin();
//              itObs != itX->second.obs.end(); ++itObs)
//          {
//              std::pair<int, int> imgFeatPair;
//              imgFeatPair.first = itObs->first;
//              imgFeatPair.second = itObs->second.id_feat;
//              //
//              P_Map.insert(std::pair<std::pair<int, int>, Vec3>(imgFeatPair, thisP));

//          }

          //prepare line
//          std::vector<std::pair<Eigen::Vector4d, Eigen::Vector4d> > allLines(0); //lanmda*l + p ; <l, p>
          std::list<std::pair<Eigen::Vector4d, Eigen::Vector4d> > allLines;
          std::list<std::pair<Eigen::Vector3d, Eigen::Vector3d> > Llist;
//          std::cout << allLines.size() << std::endl;
          for(Observations::const_iterator itObs = itX->second.obs.begin();
              itObs != itX->second.obs.end(); ++itObs)
          {
              const Pose3 pose = newdata.poses[itObs->first];
              const Vec3 C = pose.center();
              Mat3 K;
              K << 3047.0, 0.0, 2000.0,
                   0.0, 3047.0, 1500.0,
                   0.0, 0.0, 1.0;
              const Vec3 X = pose.rotation().inverse()*(K.inverse() * Vec3(itObs->second.x(0), itObs->second.x(1), 1.0) - pose.translation());
              const Vec3 l_il = X - C;
              const double lanmda = (waterPlane - C(2))/l_il(2);
              const Vec3 c_il_w = lanmda * l_il + C;
              double l_ilx2 = l_il(0)*l_il(0);
              double l_ily2 = l_il(1)*l_il(1);
              double sin_in = sqrt((l_ilx2 + l_ily2)/(l_ilx2 + l_ily2 + l_il(2)*l_il(2)));
              double ref_l_z = -sqrt( (l_ilx2 + l_ily2)/(sin_in*sin_in / (WATER_N*WATER_N)) - l_ilx2 - l_ily2);
//              double ref_l_z = waterPlane-sqrt( (l_ilx2 + l_ily2)/(sin_in*sin_in / (WATER_N*WATER_N)) - l_ilx2 - l_ily2);
              // get reslut
              Eigen::Vector4d thisLinePlane1, thisLinePlane2;
              thisLinePlane1 << l_il(1), -l_il(0), 0.0, l_il(0)*c_il_w(1)-l_il(1)*c_il_w(0);
              thisLinePlane2 << 0.0, ref_l_z, -l_il(1), l_il(1)*c_il_w(2)-ref_l_z*c_il_w(1);
              std::pair<Eigen::Vector4d, Eigen::Vector4d> thisLine;
//              std::cout << sqrt(thisLinePlane1(0)*thisLinePlane1(0) + thisLinePlane1(1)*thisLinePlane1(1) + thisLinePlane1(2)*thisLinePlane1(2)) << std::endl;
//              std::cout << sqrt(thisLinePlane2(0)*thisLinePlane2(0) + thisLinePlane2(1)*thisLinePlane2(1) + thisLinePlane2(2)*thisLinePlane2(2)) << std::endl;
              thisLine.first = thisLinePlane1 / sqrt(thisLinePlane1(0)*thisLinePlane1(0) + thisLinePlane1(1)*thisLinePlane1(1) + thisLinePlane1(2)*thisLinePlane1(2));
              thisLine.second = thisLinePlane2 / sqrt(thisLinePlane2(0)*thisLinePlane2(0) + thisLinePlane2(1)*thisLinePlane2(1) + thisLinePlane2(2)*thisLinePlane2(2));

              allLines.push_back(thisLine);

              //
              std::pair<Eigen::Vector3d, Eigen::Vector3d> thisL;
              thisL.first = C;
              thisL.second = X;
              Llist.push_back(thisL);
          }
          // get resX
          {
//              Eigen::MatrixX4d A(allLines.size()*2, 4);
              Eigen::MatrixXd A(allLines.size()*2, 4);
              A = Eigen::MatrixXd::Zero(allLines.size()*2, 4);
              int iter = 0;
              for(std::list<std::pair<Eigen::Vector4d, Eigen::Vector4d> >::const_iterator itAllL = allLines.begin();
                  itAllL != allLines.end(); ++itAllL, ++iter)//std::size_t iter = 0; iter < allLines.size(); ++iter)
              {

                  A(iter*2, 0) = itAllL->first(0);
                  A(iter*2, 1) = itAllL->first(1);
                  A(iter*2, 2) = itAllL->first(2);
                  A(iter*2, 3) = itAllL->first(3);
                  //
                  A(iter*2+1, 0) = itAllL->second(0);
                  A(iter*2+1, 1) = itAllL->second(1);
                  A(iter*2+1, 2) = itAllL->second(2);
                  A(iter*2+1, 3) = itAllL->second(3);

              }
              // do SVD
              Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU|Eigen::ComputeThinV);
              Eigen::Matrix4Xd V = svd.matrixV();
    //                std::cout << " svd.matrixV() : " << svd.matrixV() << std::endl;
    //                std::cout << " V : " << std::endl << V << std::endl;
              int col = 3;//allLines.size()*2-1;
    //                std::cout << " col : " << col << std::endl;
              Vec3 resX = Vec3(V(0,col)/V(3,col), V(1,col)/V(3,col), V(2,col)/V(3,col));

              undiswaterdata.structure[itX->first].X = resX;
          }

      }


      Save(undiswaterdata,
        stlplus::create_filespec(sOutDir, "cloud_and_poses_undiswater", ".ply"),
        ESfM_Data(ALL));

      //save for analyse error
      std::ofstream out_X, out_xm;
      out_X.open(sOutDir+"/X_undisWater.txt");
      out_xm.open(sOutDir+"/xm_undisWater.txt");
      for(Landmarks::const_iterator itX = undiswaterdata.structure.begin();
          itX != undiswaterdata.structure.end(); ++itX)
      {
          out_xm << itX->second.obs.size() << " ";
          for(Observations::const_iterator itObs = itX->second.obs.begin();
              itObs != itX->second.obs.end(); ++itObs)
          {
              out_xm << itObs->first << " " << itObs->second.id_feat << " " << itObs->second.x(0) << " " << itObs->second.x(1) << " ";

          }
          out_xm << itX->second.X(0) << " "
                 << itX->second.X(1) << " "
                 << itX->second.X(2) << std::endl;
          out_X << itX->second.obs.begin()->first << " "
                << itX->second.obs.begin()->second.id_feat << " "
                << itX->second.X(0) << " "
                << itX->second.X(1) << " "
                << itX->second.X(2) << std::endl;

      }
      out_X.close();
      out_xm.close();



  }
  return EXIT_SUCCESS;
}
