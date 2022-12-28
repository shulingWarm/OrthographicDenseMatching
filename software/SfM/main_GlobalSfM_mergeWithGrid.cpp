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

  std::string sRootPath;
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
  bool bNeedSixteenth = true;
  std::string sSubPieceIdList = "11_22_33_44";
  std::string sGridDataPath = "";
  std::string sPieceDataPath = "";

//  cmd.add( make_option('i', sRootPath, "input_folder") );
//  cmd.add( make_option('m', sMatchesDir, "matchdir") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('r', iRotationAveragingMethod, "rotationAveraging") );
  cmd.add( make_option('t', iTranslationAveragingMethod, "translationAveraging") );
  cmd.add( make_option('f', sIntrinsic_refinement_options, "refineIntrinsics") );
  cmd.add( make_switch('P', "prior_usage") );
  cmd.add( make_option('d', needDivideRegion, "need do dividing region process") );
  cmd.add( make_option('n', rtFilePath, "file path of R t file, n_gps_imu_321_x.res") );
  cmd.add( make_option('g', sGridDataPath, "GridDataPath") );
  cmd.add( make_option('p', sPieceDataPath, "PieceDataPath") );
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
  //
  /// Load datum
  ///
  SfM_Data sfm_data_grid;
  if (!Load(sfm_data_grid, sGridDataPath, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data sfm_data_grid file \""<< sGridDataPath << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }
  //
  std::vector<SfM_Data> sfm_data_piece_list;
  std::vector<std::string> subPieceIdList;
  boost::split(subPieceIdList, sSubPieceIdList, boost::is_any_of("_"));
  //
  for(std::size_t subPListId = 0; subPListId < subPieceIdList.size(); ++subPListId)
  {
      SfM_Data sfm_data_piece;
      if (!Load(sfm_data_piece, sPieceDataPath+"/sfm_data_"+subPieceIdList[subPListId]+".bin", ESfM_Data(ALL))) {
        std::cerr << std::endl
          << "The input SfM_Data file \""<< sPieceDataPath+"/sfm_data_s_"+subPieceIdList[subPListId]+".bin" << "\" cannot be read." << std::endl;
        return EXIT_FAILURE;
      }
      sfm_data_piece_list.push_back(sfm_data_piece);
  }
  //
  std::map<int, bool> fixedPoseIdList;
  for(Views::const_iterator itGridView = sfm_data_grid.views.begin();
      itGridView != sfm_data_grid.views.end(); ++itGridView)
  {
      fixedPoseIdList[itGridView->first] = true;
  }
  //
  /// merge pieces all
  ///
  SfM_Data sfm_data_piece_all = sfm_data_piece_list[0];
  for(std::size_t subPListId = 1; subPListId < subPieceIdList.size(); ++subPListId)
  {
      SfM_Data* sfm_data_this = &sfm_data_piece_list[subPListId];
      // views and poses
      for(Views::const_iterator itPieceView = sfm_data_this->views.begin();
          itPieceView != sfm_data_this->views.end(); ++itPieceView)
      {
          Views::const_iterator itPieceAllView = sfm_data_piece_all.views.find(itPieceView->first);
          if (itPieceAllView == sfm_data_piece_all.views.end())
          {
              const View * viewData = itPieceView->second.get();
              sfm_data_piece_all.views[itPieceView->first] = itPieceView->second;
              sfm_data_piece_all.poses[viewData->id_pose] = sfm_data_this->poses[viewData->id_pose];
          }
      }
      // landmarks
      Landmarks::const_iterator itSEnd = sfm_data_piece_all.structure.end();
      --itSEnd;
      const int landmark_startId = itSEnd->first + 1;
      for (Landmarks::const_iterator itPieceStruct = sfm_data_this->structure.begin();
           itPieceStruct != sfm_data_this->structure.end(); ++itPieceStruct)
      {
          sfm_data_piece_all.structure[landmark_startId + itPieceStruct->first] = itPieceStruct->second;
      }
  }
  //
  /// merge
  /// views and poses
  SfM_Data sfm_data_all = sfm_data_grid;
  for (Views::const_iterator itPieceView = sfm_data_piece_all.views.begin();
      itPieceView != sfm_data_piece_all.views.end(); ++itPieceView)
  {
      Views::const_iterator itGridView = sfm_data_all.views.find(itPieceView->first);
      if(itGridView == sfm_data_all.views.end())
      {
          const View * viewData = itPieceView->second.get();
          //
          sfm_data_all.views[itPieceView->first] = itPieceView->second;
          sfm_data_all.poses[viewData->id_pose] = sfm_data_piece_all.poses[viewData->id_pose];
      }
      //
  }
  //
  /// record fixed landmarks
  ///
  std::map<int, std::map<int, int> > fixedSList; //<imgId, <featId, structureId> >
  std::map<int, bool> fixedXIdList;
  for (Landmarks::const_iterator itGridStruct = sfm_data_grid.structure.begin();
       itGridStruct != sfm_data_grid.structure.end(); ++itGridStruct)
  {
      fixedXIdList[itGridStruct->first] = true;
      for(Observations::const_iterator itGridObs = itGridStruct->second.obs.begin();
          itGridObs != itGridStruct->second.obs.end(); ++itGridObs)
      {
          int imgId = itGridObs->first;
          int featId = itGridObs->second.id_feat;
          fixedSList[imgId][featId] = itGridStruct->first;
      }
  }
  //
  /// merge and record changed landmarks
  ///
  std::map<int, int> recordChangeIdList; //<large_struct_id, piece_struct_id>
  Landmarks::const_iterator itSEnd = sfm_data_grid.structure.end();
  --itSEnd;
  const int startSId = itSEnd->first + 1;
  for (Landmarks::const_iterator itPieceStruct = sfm_data_piece_all.structure.begin();
       itPieceStruct != sfm_data_piece_all.structure.end(); ++itPieceStruct)
  {
      for (Observations::const_iterator itPieceObs = itPieceStruct->second.obs.begin();
          itPieceObs != itPieceStruct->second.obs.end(); ++itPieceObs)
      {
          std::map<int, std::map<int, int> >::iterator findImgId = fixedSList.find(itPieceObs->first);
          if (findImgId != fixedSList.end()) // find it (imgId)
          {
              std::map<int, int>::const_iterator findFeatId = findImgId->second.find(itPieceObs->second.id_feat);
              if (findFeatId != findImgId->second.end()) // find it (featId)
              {
                  recordChangeIdList[findFeatId->second] = itPieceStruct->first;
                  continue;
              }
          }
          //
          sfm_data_all.structure[startSId+itPieceStruct->first] = itPieceStruct->second;
      }
  }
  //
  for (std::map<int, int>::const_iterator itRecord = recordChangeIdList.begin();
       itRecord != recordChangeIdList.end(); ++itRecord)
  {
      Landmarks::iterator itAllS = sfm_data_all.structure.find(itRecord->first);
      Landmarks::const_iterator itPieceS = sfm_data_piece_all.structure.find(itRecord->second);
      for(Observations::const_iterator itPieceObs = itPieceS->second.obs.begin();
          itPieceObs != itPieceS->second.obs.end(); ++itPieceObs)
      {
          Observations::iterator itAllObs = itAllS->second.obs.find(itPieceObs->first); // find the image id
          if(itAllObs == itAllS->second.obs.end())
          {
              itAllS->second.obs[itPieceObs->first] = itPieceObs->second;
          }else{
              // for dubug
              if(itAllObs->second.id_feat != itPieceObs->second.id_feat)
              {
                  std::cout << "landmarks id feat error " << std::endl;
              }
          }
      }
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
    sfm_data_all,
    sOutDir,
    stlplus::create_filespec(sOutDir, "Reconstruction_Report.html"));

  // Configure the features_provider & the matches_provider
//  sfmEngine.SetFeaturesProvider(feats_provider.get());
//  sfmEngine.SetMatchesProvider(matches_provider_all.get());

  // Configure reconstruction parameters
  sfmEngine.Set_Intrinsics_Refinement_Type(intrinsic_refinement_options);
  b_use_motion_priors = cmd.used('P');
  sfmEngine.Set_Use_Motion_Prior(b_use_motion_priors);

  // Configure motion averaging method
  sfmEngine.SetRotationAveragingMethod(
    ERotationAveragingMethod(iRotationAveragingMethod));
  sfmEngine.SetTranslationAveragingMethod(
    ETranslationAveragingMethod(iTranslationAveragingMethod));

  if (sfmEngine.Process_step23(fixedPoseIdList, fixedXIdList))//_GCP())//_threeCameras_imu_gps_init_new(rtFilePath))//Process())//Process_threeCameras_imu_gps(sMatchesDir))//Process())
  {
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
