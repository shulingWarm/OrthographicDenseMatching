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

#include <cstdlib>
#include <memory>
#include <string>
#include <sys/stat.h>
#include <list>

#include <boost/algorithm/string.hpp>

using namespace openMVG;
using namespace openMVG::sfm;

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
  bool needDivideRegion = false;
  bool needDivideForDensify = true;
  bool needDivideGroup = true;
  std::string rtFilePath = "";
  std::string subPieceIdListStr = "1_2_3_4";
  bool bNeedSixteenth = true;

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('m', sMatchesDir, "matchdir") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('r', iRotationAveragingMethod, "rotationAveraging") );
  cmd.add( make_option('t', iTranslationAveragingMethod, "translationAveraging") );
  cmd.add( make_option('f', sIntrinsic_refinement_options, "refineIntrinsics") );
  cmd.add( make_switch('P', "prior_usage") );
  cmd.add( make_option('d', needDivideRegion, "need do dividing region process") );
  cmd.add( make_option('n', rtFilePath, "file path of R t file, n_gps_imu_321_x.res") );
  cmd.add( make_option('g', needDivideGroup, "need do dividing for mutli-thread depth map process according groups") );
  cmd.add( make_option('s', subPieceIdListStr, "subPieceIdListStr") );
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
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS|CONTROL_POINTS))) { ///!!!
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

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
  //prepare input
  std::vector<std::string> subPieceIdList;
  boost::split(subPieceIdList, subPieceIdListStr, boost::is_any_of("_"));
//  std::cout << refSfmDataIdList.size() << std::endl;
//  std::cout << "split id: " << std::endl;
//  std::cout << refIdList.size()<< std::endl;

  std::shared_ptr<Matches_Provider> matches_provider_all = std::make_shared<Matches_Provider>();

  std::vector<std::shared_ptr<Matches_Provider> > sub_Matches_Provider_List;
  if(subPieceIdList.size() != 0)
  {
      /// load
      for(std::size_t sPId = 0; sPId < subPieceIdList.size(); ++sPId)
      {
          std::shared_ptr<Matches_Provider> matches_provider_sub = std::make_shared<Matches_Provider>();
          std::string sMatchesDir_sub = sMatchesDir;//sfm_data.s_root_path + "/matches_"+subPieceIdList[sPId];
          std::cout << "sub matches dir : " << sMatchesDir_sub << std::endl;

          if // Try to read the two matches file formats
          (
            !(matches_provider_sub->load(sfm_data, stlplus::create_filespec(sMatchesDir_sub, "matches_"+subPieceIdList[sPId]+".e.txt")) ||
              matches_provider_sub->load(sfm_data, stlplus::create_filespec(sMatchesDir_sub, "matches_"+subPieceIdList[sPId]+".e.bin")))
          )
          {
            std::cerr << std::endl
              << "Invalid matches file. " << sMatchesDir_sub << std::endl;
            return EXIT_FAILURE;
          }

          sub_Matches_Provider_List.push_back(matches_provider_sub);

      }

      ///marge
      matching::PairWiseMatches map_PutativesMatches_all;
      for(std::vector<std::shared_ptr<Matches_Provider> >::const_iterator itSMP = sub_Matches_Provider_List.begin();
          itSMP != sub_Matches_Provider_List.end(); ++itSMP)
      {
          for (matching::PairWiseMatches::const_iterator iterPM = (*(itSMP))->pairWise_matches_.begin();
            iterPM != (*(itSMP))->pairWise_matches_.end(); ++iterPM)
          {
              map_PutativesMatches_all.insert(*iterPM);
          }
      }
      matches_provider_all->pairWise_matches_ = map_PutativesMatches_all;



  }else{
      std::cout << "input sub-pieces id error" << std::endl;
      return EXIT_FAILURE;
  }



//  std::shared_ptr<Matches_Provider> matches_provider_all = std::make_shared<Matches_Provider>();
//  std::shared_ptr<Matches_Provider> matches_provider_1 = std::make_shared<Matches_Provider>();
//  std::shared_ptr<Matches_Provider> matches_provider_2 = std::make_shared<Matches_Provider>();
//  std::shared_ptr<Matches_Provider> matches_provider_3 = std::make_shared<Matches_Provider>();
//  std::shared_ptr<Matches_Provider> matches_provider_4 = std::make_shared<Matches_Provider>();
//  /// load datum
//  {
//      std::string sMatchesDir_1 = sfm_data.s_root_path + "/matches_1";
//      std::string sMatchesDir_2 = sfm_data.s_root_path + "/matches_2";
//      std::string sMatchesDir_3 = sfm_data.s_root_path + "/matches_3";
//      std::string sMatchesDir_4 = sfm_data.s_root_path + "/matches_4";
//      if // Try to read the two matches file formats
//      (
//        !(matches_provider_1->load(sfm_data, stlplus::create_filespec(sMatchesDir_1, "matches.e.txt")) ||
//          matches_provider_1->load(sfm_data, stlplus::create_filespec(sMatchesDir_1, "matches.e.bin")))
//      )
//      {
//        std::cerr << std::endl
//          << "Invalid matches file. " << sMatchesDir_1 << std::endl;
//        return EXIT_FAILURE;
//      }
//      if // Try to read the two matches file formats
//      (
//        !(matches_provider_2->load(sfm_data, stlplus::create_filespec(sMatchesDir_2, "matches.e.txt")) ||
//          matches_provider_2->load(sfm_data, stlplus::create_filespec(sMatchesDir_2, "matches.e.bin")))
//      )
//      {
//        std::cerr << std::endl
//          << "Invalid matches file. " << sMatchesDir_2 << std::endl;
//        return EXIT_FAILURE;
//      }
//      if // Try to read the two matches file formats
//      (
//        !(matches_provider_3->load(sfm_data, stlplus::create_filespec(sMatchesDir_3, "matches.e.txt")) ||
//          matches_provider_3->load(sfm_data, stlplus::create_filespec(sMatchesDir_3, "matches.e.bin")))
//      )
//      {
//        std::cerr << std::endl
//          << "Invalid matches file. " << sMatchesDir_3 << std::endl;
//        return EXIT_FAILURE;
//      }
//      if // Try to read the two matches file formats
//      (
//        !(matches_provider_4->load(sfm_data, stlplus::create_filespec(sMatchesDir_4, "matches.e.txt")) ||
//          matches_provider_4->load(sfm_data, stlplus::create_filespec(sMatchesDir_4, "matches.e.bin")))
//      )
//      {
//        std::cerr << std::endl
//          << "Invalid matches file. " << sMatchesDir_4 << std::endl;
//        return EXIT_FAILURE;
//      }
//  }
//  ///merge 1234 to all
//  {
//      matching::PairWiseMatches map_PutativesMatches_all;

//      for (matching::PairWiseMatches::const_iterator iter1 = matches_provider_1->pairWise_matches_.begin();
//        iter1 != matches_provider_1->pairWise_matches_.end(); ++iter1)
//      {
//          map_PutativesMatches_all.insert(*iter1);
//      }
//      for (matching::PairWiseMatches::const_iterator iter2 = matches_provider_2->pairWise_matches_.begin();
//        iter2 != matches_provider_2->pairWise_matches_.end(); ++iter2)
//      {
//          map_PutativesMatches_all.insert(*iter2);
//      }
//      for (matching::PairWiseMatches::const_iterator iter3 = matches_provider_3->pairWise_matches_.begin();
//        iter3 != matches_provider_3->pairWise_matches_.end(); ++iter3)
//      {
//          map_PutativesMatches_all.insert(*iter3);
//      }
//      for (matching::PairWiseMatches::const_iterator iter4 = matches_provider_4->pairWise_matches_.begin();
//        iter4 != matches_provider_4->pairWise_matches_.end(); ++iter4)
//      {
//          map_PutativesMatches_all.insert(*iter4);
//      }

//      matches_provider_all->pairWise_matches_ = map_PutativesMatches_all;
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

  //---------------------------------------
  // Global SfM reconstruction process
  //---------------------------------------

  openMVG::system::Timer timer;
  GlobalSfMReconstructionEngine_RelativeMotions sfmEngine(
    sfm_data,
    sOutDir,
    stlplus::create_filespec(sOutDir, "Reconstruction_Report.html"));


  ///!!!
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

  // Configure the features_provider & the matches_provider
  sfmEngine.SetFeaturesProvider(feats_provider.get());
  sfmEngine.SetMatchesProvider(matches_provider_all.get());

  // Configure reconstruction parameters
  sfmEngine.Set_Intrinsics_Refinement_Type(intrinsic_refinement_options);
  b_use_motion_priors = cmd.used('P');
  sfmEngine.Set_Use_Motion_Prior(b_use_motion_priors);

  // Configure motion averaging method
  sfmEngine.SetRotationAveragingMethod(
    ERotationAveragingMethod(iRotationAveragingMethod));
  sfmEngine.SetTranslationAveragingMethod(
    ETranslationAveragingMethod(iTranslationAveragingMethod));

  if (sfmEngine.Process())//_GCP_GPS_6(rtFilePath))//_GCP())//_threeCameras_imu_gps_init_new(rtFilePath))//Process())//Process_threeCameras_imu_gps(sMatchesDir))//Process())
  {
      {
//          SfM_Data newdata = sfmEngine.Get_SfM_Data();
//          newdata.intrinsics = sfmEngine.Get_SfM_Data().intrinsics;
//          Mat3 changeR;
//          Vec3 changet;
//          float s;
//          {
//              //run2
////              s = 141.783;
////              changeR << 13.886957168579, 139.231246948242, -22.894561767578,
////                         140.898696899414, -12.451700210571, 9.739796638489,
////                         7.553865909576, -23.705772399902, -139.58268737793;
////              changet << Vec3(89.912673950195, -147.658294677734, 163.292388916016);
////              s = 182.257;
////              changeR << 42.258850097656, 175.210388183594, -27.074182510376,
////                         177.236633300781, -42.432579040527, 2.03838300705,
////                         -4.343770980835, -26.801080703735, -180.223022460938;
////              changet << Vec3(106.549041748047, 63.339290618896, 170.729095458984);

//              s = 307.739;
//              changeR << -307.702819824219, 1.30886900425, -4.563062667847,
//                         2.248091220856, 300.714599609375, -65.339447021484,
//                         4.181000709534, -65.365013122559, -300.688385009766;
//              changet << Vec3(39.385131835938, 147.40007019043, 366.163726806641);




////              s = 519.0267046517480;
//                //folder 16 Z = 0
////              changeR << -37.933506011963, -512.732482910156, 71.099311828613,
////                         -514.708435058594, 44.935962677002, 49.443969726562,
////                         -54.999961853027, -66.894119262695, -511.750793457031;
////              changet = Vec3(-414.110778808594, -43.797954559326, 543.712890625);
////              changeR << 1.0, 0.0, 0.0,
////                         0.0, 1.0, 0.0,
////                         0.0, 0.0, 1.0;

//              changeR = changeR/s;
//              changet = changet/s;

//              //s = 0.98445;//0.30558845;//1.0;//3.1908;
//          }

//          for (const auto & pose_it : sfmEngine.Get_SfM_Data().poses)
//          {
//            const IndexT indexPose = pose_it.first;

//            const Pose3 & pose = pose_it.second;
//            const Mat3 R = pose.rotation();
//            const Vec3 C = pose.center();

//            Mat3 Rnew = R * changeR.inverse();
////            Vec3 Cnew = changet + Rnew.inverse() * R * C;
//            Vec3 Cnew = changeR * s * C + changet;

//            newdata.poses[indexPose] = Pose3(Rnew, Cnew);

//          }

//          newdata.structure = sfmEngine.Get_SfM_Data().structure;
//          //set all X
//          for(Landmarks::iterator itX = newdata.structure.begin();
//              itX != newdata.structure.end(); ++itX)
//          {
//              const Vec3 & X = itX->second.X;

//              Vec3 newX = changeR * s * X + changet;
////                  Vec3 newX = R_gc * X + t_gc;
////                  Vec3 newX = R0.inverse() * X;
//              itX->second.X = newX;
//          }

//          //-- Export to disk computed scene (data & visualizable results)
//          std::cout << "...Export SfM_Data to disk." << std::endl;

//          Save(newdata,
//            stlplus::create_filespec(sOutDir, "sfm_data2", ".bin"),
//            ESfM_Data(ALL));

//          Save(newdata,
//            stlplus::create_filespec(sOutDir, "cloud_and_poses2", ".ply"),
//            ESfM_Data(ALL));

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

    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
