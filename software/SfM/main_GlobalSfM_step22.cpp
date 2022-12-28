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

#include "Eigen/Dense"

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

  std::string sGridFilename;
  std::string sMatchesDir;
  std::string sOutDir = "";
  std::string sSubPiecesPath = "";
  int iRotationAveragingMethod = int (ROTATION_AVERAGING_L2);
  int iTranslationAveragingMethod = int (TRANSLATION_AVERAGING_SOFTL1);
  std::string sIntrinsic_refinement_options = "ADJUST_ALL";
  bool b_use_motion_priors = false;
  bool needDivideRegion = false;
  bool needDivideForDensify = true;
  bool needDivideGroup = true;
  std::string rtFilePath = "";
  int iBId = 0;
  std::string sSubPieceIdList = "11_21_31_41";
  std::string sSubMatchIdList = "1_2_3_4";

  cmd.add( make_option('i', sGridFilename, "input_file") );
  cmd.add( make_option('m', sMatchesDir, "matchdir") );
  cmd.add( make_option('M', sSubMatchIdList, "subMatchIdList") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('r', iRotationAveragingMethod, "rotationAveraging") );
  cmd.add( make_option('t', iTranslationAveragingMethod, "translationAveraging") );
  cmd.add( make_option('f', sIntrinsic_refinement_options, "refineIntrinsics") );
  cmd.add( make_switch('P', "prior_usage") );
  cmd.add( make_option('c', sSubPiecesPath, "subPiecesPath") );
  cmd.add( make_option('b', iBId, "bId") );
  cmd.add( make_option('s', sSubPieceIdList, "subPieceIdList") );

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
  SfM_Data sfm_data_grid;
  if (!Load(sfm_data_grid, sGridFilename, ESfM_Data(VIEWS|INTRINSICS|CONTROL_POINTS))) { ///!!!
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sGridFilename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }
  //
  /// Load and set views and poses
  ///
  std::vector<SfM_Data> sfm_data_list;
  std::vector<std::string> subPieceIdList;
  boost::split(subPieceIdList, sSubPieceIdList, boost::is_any_of("_"));
  //
  for(std::size_t subPListId = 0; subPListId < subPieceIdList.size(); ++subPListId)
  {
      SfM_Data sfm_data;
      if (!Load(sfm_data, sSubPiecesPath+"/sfm_data_"+subPieceIdList[subPListId]+".bin", ESfM_Data(ALL))) {
        std::cerr << std::endl
          << "The input SfM_Data file \""<< sSubPiecesPath+"/sfm_data_s_"+subPieceIdList[subPListId]+".bin" << "\" cannot be read." << std::endl;
        return EXIT_FAILURE;
      }
      sfm_data_list.push_back(sfm_data);
  }
  //
  /// set sfm_data_grid
  for(std::size_t subPListId = 0; subPListId < sfm_data_list.size(); ++subPListId)
  {
      SfM_Data thisData = sfm_data_list[subPListId];
      for(Views::const_iterator itSubView = thisData.views.begin();
          itSubView != thisData.views.end(); ++itSubView)
      {
          if(sfm_data_grid.views.find(itSubView->first) != thisData.views.end()) // find it
          {
              const View * viewData = itSubView->second.get();
              sfm_data_grid.poses[viewData->id_pose] = sfm_data_grid.poses[viewData->id_pose];
          }
      }
  }
  //
  /// check for test
  for (Views::const_iterator itGridView = sfm_data_grid.views.begin();
       itGridView != sfm_data_grid.views.end(); ++itGridView)
  {
      const View * viewData = itGridView->second.get();
      Poses::const_iterator itGridPose = sfm_data_grid.poses.find(viewData->id_pose);
      if(itGridPose == sfm_data_grid.poses.end())
      {
          std::cout << "!!! error view id : " << itGridView->first << std::endl;
      }
  }
  //
  //
  /// load matches
  ///
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
  if (!feats_provider->load(sfm_data_grid, sMatchesDir, regions_type)) {
    std::cerr << std::endl
      << "Invalid features." << std::endl;
    return EXIT_FAILURE;
  }
  // Matches reading
  /// load
  std::vector<std::shared_ptr<Matches_Provider> > sub_Matches_Provider_List;
  std::vector<std::string> subMatchIdList;
  boost::split(subMatchIdList, sSubMatchIdList, boost::is_any_of("_"));
  for(std::size_t sPId = 0; sPId < subMatchIdList.size(); ++sPId)
  {
      std::shared_ptr<Matches_Provider> matches_provider_sub = std::make_shared<Matches_Provider>();
      std::string sMatchesDir_sub = sMatchesDir;
      std::cout << "sub matches dir : " << sMatchesDir_sub << std::endl;

      if // Try to read the two matches file formats
      (
        !(matches_provider_sub->load(sfm_data_grid, stlplus::create_filespec(sMatchesDir_sub, "matches_"+subMatchIdList[sPId]+".e.txt")) ||
          matches_provider_sub->load(sfm_data_grid, stlplus::create_filespec(sMatchesDir_sub, "matches_"+subMatchIdList[sPId]+".e.bin")))
      )
      {
        std::cerr << std::endl
          << "Invalid matches file. " << sMatchesDir_sub << std::endl;
        return EXIT_FAILURE;
      }

      sub_Matches_Provider_List.push_back(matches_provider_sub);

  }
  ///merge matches
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
  std::shared_ptr<Matches_Provider> matches_provider_all = std::make_shared<Matches_Provider>();
  matches_provider_all->pairWise_matches_ = map_PutativesMatches_all;
  //
  //
  //
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
    sfm_data_grid,
    sOutDir,
    stlplus::create_filespec(sOutDir, "Reconstruction_Report.html"));


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

  if (sfmEngine.Process_step22())//_GCP_GPS(rtFilePath))//_GCP())//_threeCameras_imu_gps_init_new(rtFilePath))//Process())//Process_threeCameras_imu_gps(sMatchesDir))//Process())
  {
    std::cout << std::endl << " Total Ac-Global-Sfm took (s): " << timer.elapsed() << std::endl;

    std::cout << "...Generating SfM_Report.html" << std::endl;
    Generate_SfM_Report(sfmEngine.Get_SfM_Data(),
      stlplus::create_filespec(sOutDir, "SfMReconstruction_Report.html"));

    std::stringstream ss;
    ss << iBId;
    std::string s;
    ss >> s;

    //-- Export to disk computed scene (data & visualizable results)
    std::cout << "...Export SfM_Data to disk." << std::endl;
    Save(sfmEngine.Get_SfM_Data(),
      stlplus::create_filespec(sOutDir, "sfm_data_"+iBId, ".bin"),
      ESfM_Data(ALL));

    Save(sfmEngine.Get_SfM_Data(),
      stlplus::create_filespec(sOutDir, "cloud_and_poses", ".ply"),
      ESfM_Data(ALL));

    return EXIT_SUCCESS;

  }

  return EXIT_FAILURE;
}
