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
  std::string subPieceIdListStr = "1_2_3_4";
  bool bNeedSixteenth = true;

  cmd.add( make_option('i', sRootPath, "input_folder") );
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

//  // Load input SfM_Data scene
//  SfM_Data sfm_data;
//  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS|CONTROL_POINTS))) { ///!!!
//    std::cerr << std::endl
//      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
//    return EXIT_FAILURE;
//  }

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

//  // Features reading
//  std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
//  if (!feats_provider->load(sfm_data, sMatchesDir, regions_type)) {
//    std::cerr << std::endl
//      << "Invalid features." << std::endl;
//    return EXIT_FAILURE;
//  }
  // Matches reading
  //prepare input
//  std::vector<std::string> subPieceIdList;
//  boost::split(subPieceIdList, subPieceIdListStr, boost::is_any_of("_"));
//  std::cout << refSfmDataIdList.size() << std::endl;
//  std::cout << "split id: " << std::endl;
//  std::cout << refIdList.size()<< std::endl;

//  std::shared_ptr<Matches_Provider> matches_provider_all = std::make_shared<Matches_Provider>();

//  std::vector<std::shared_ptr<Matches_Provider> > sub_Matches_Provider_List;
//  if(subPieceIdList.size() != 0)
//  {
//      /// load
//      for(std::size_t sPId = 0; sPId < subPieceIdList.size(); ++sPId)
//      {
//          std::shared_ptr<Matches_Provider> matches_provider_sub = std::make_shared<Matches_Provider>();
//          std::string sMatchesDir_sub = sfm_data.s_root_path + "/matches_"+subPieceIdList[sPId];
//          std::cout << "sub matches dir : " << sMatchesDir_sub << std::endl;

//          if // Try to read the two matches file formats
//          (
//            !(matches_provider_sub->load(sfm_data, stlplus::create_filespec(sMatchesDir_sub, "matches.e.txt")) ||
//              matches_provider_sub->load(sfm_data, stlplus::create_filespec(sMatchesDir_sub, "matches.e.bin")))
//          )
//          {
//            std::cerr << std::endl
//              << "Invalid matches file. " << sMatchesDir_sub << std::endl;
//            return EXIT_FAILURE;
//          }

//          sub_Matches_Provider_List.push_back(matches_provider_sub);

//      }

//      ///marge
//      matching::PairWiseMatches map_PutativesMatches_all;
//      for(std::vector<std::shared_ptr<Matches_Provider> >::const_iterator itSMP = sub_Matches_Provider_List.begin();
//          itSMP != sub_Matches_Provider_List.end(); ++itSMP)
//      {
//          for (matching::PairWiseMatches::const_iterator iterPM = (*(itSMP))->pairWise_matches_.begin();
//            iterPM != (*(itSMP))->pairWise_matches_.end(); ++iterPM)
//          {
//              map_PutativesMatches_all.insert(*iterPM);
//          }
//      }
//      matches_provider_all->pairWise_matches_ = map_PutativesMatches_all;



//  }else{
//      std::cout << "input sub-pieces id error" << std::endl;
//      return EXIT_FAILURE;
//  }



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

  /// Load datum
  ///
  std::vector<SfM_Data> sfm_data_list;
  std::vector<std::string> subPieceIdList;
  boost::split(subPieceIdList, subPieceIdListStr, boost::is_any_of("_"));
  //
  for(std::size_t subPListId = 0; subPListId < subPieceIdList.size(); ++subPListId)
  {
      SfM_Data sfm_data;
      if (!Load(sfm_data, sRootPath+"/matches_"+subPieceIdList[subPListId]+"/sfm_data.bin", ESfM_Data(ALL))) {
        std::cerr << std::endl
          << "The input SfM_Data file \""<< sRootPath+"/matches_"+subPieceIdList[subPListId]+"/sfm_data.bin" << "\" cannot be read." << std::endl;
        return EXIT_FAILURE;
      }
      // Rt
      {
          {
              SfM_Data newSfMData;
              newSfMData = sfm_data;
    //          std::stringstream bIdss;
    //          bIdss << bId;
    //          std::string bIds;
    //          bIdss >> bIds;

              std::map<int, Vec3> c_gps;
              std::ifstream in_gps;
              std::string inGpsPath = sRootPath + "/n_gps_imu_321.res";
              std::cout << "inGpsPath : " << inGpsPath << std::endl;
              in_gps.open(inGpsPath);
              if(in_gps)
              {
                  int count;
                  in_gps >> count;
                  int line = 0;
                  while(line < count && !in_gps.eof())
                  {
                      int id;
                      Vec3 xyz, temp;
                      int tempSatus;

                      in_gps >> id >> xyz(0) >> xyz(1) >> xyz(2) >> temp(0) >> temp(1) >> temp(2) >> tempSatus;

                      std::pair<int, Vec3> tempPair;
                      tempPair.first = id;
                      tempPair.second = xyz;
                      c_gps.insert(tempPair);

                      ++line;
                  }

              }else{
                  std::cout << "open n_gps_imu_321.res file failed, please check it !" << std::endl;
                  return EXIT_FAILURE;
              }
              in_gps.close();

              double length_gps = 0.0;
              double length_c = 0.0;
              sfm::Poses::iterator itEnd = sfm_data.poses.end();
              -- itEnd;
              int length_N = 0;
//              for(std::size_t poseListId1 = 0; poseListId1 != computeImgMap.size(); ++poseListId1)
              for(sfm::Poses::iterator itPose = sfm_data.poses.begin();
                  itPose != itEnd; ++itPose)
              {
                  int pId1 = itPose->first;//computeImgMap[poseListId1];
                  const Vec3 c1 = sfm_data.poses[pId1].center();

                  if (pId1/100000 != 2)
                  {
                      continue;
                  }
                  const Vec3 cg1 = c_gps.find(pId1)->second;

                  sfm::Poses::iterator itPose2 = itPose;
                  ++ itPose2;
                  for(;itPose2 != sfm_data.poses.end(); ++itPose2)
//                  for(std::size_t poseListId2 = poseListId1; poseListId2 != computeImgMap.size(); ++poseListId2)
                  {
                      int pId2 = itPose2->first;//computeImgMap[poseListId2];
                      const Vec3 c2 = sfm_data.poses[pId2].center();

                      if (pId2/100000 != 2)
                      {
                          continue;
                      }
                      const Vec3 cg2 = c_gps.find(pId2)->second;

                      double lg = (cg1-cg2).norm();
                      double lc = (c1-c2).norm();

                      length_gps += lg;
                      length_c += lc;

                      ++length_N;

                  }

              }
              length_gps /= (double)(length_N);
              length_c /= (double)(length_N);

              double s = length_gps / length_c;

              //compute H -> Xg = H * Xc
              std::vector<Vec3> Xc, Xg;
//              for(std::size_t poseListId = 0; poseListId != computeImgMap.size(); ++poseListId)
              for(sfm::Poses::iterator itPose = sfm_data.poses.begin();
                  itPose != itEnd; ++itPose)
              {
                  int pId = itPose->first;
                  if (pId/100000 != 2)
                  {
                      continue;
                  }
                  //prepare Xg
                  Vec3 Xg_temp;
                  Xg_temp = c_gps.find(pId)->second;
                  Xg.push_back(Xg_temp);

                  //prepare Xc
                  Vec3 Xc_temp;
                  Xc_temp = s * sfm_data.poses[pId].center();//itPose->second.center();
                  Xc.push_back(Xc_temp);

              }

              //prepare Xg
              Vec3 barycenter_g;
              barycenter_g << 0.0, 0.0, 0.0;
              for(std::size_t i = 0; i < Xg.size(); ++i)
              {
                  barycenter_g += Xg[i];
              }
              barycenter_g /= (double)(Xg.size());
              //
              std::vector<Vec3> Xg_bary;
              for(std::size_t i = 0; i < Xg.size(); ++i)
              {
                  Vec3 Xg_bary_temp;
                  Xg_bary_temp = Xg[i] - barycenter_g;
                  Xg_bary.push_back(Xg_bary_temp);
              }
              //prepare Xc
              Vec3 barycenter_c;
              barycenter_c << 0.0, 0.0, 0.0;
              for(std::size_t i = 0; i < Xc.size(); ++i)
              {
                  barycenter_c += Xc[i];
              }
              barycenter_c /= (double)(Xc.size());

              std::vector<Vec3> Xc_bary;
              for(std::size_t i = 0; i < Xc.size(); ++i)
              {
                  Vec3 Xc_bary_temp;
                  Xc_bary_temp = Xc[i] - barycenter_c;
                  Xc_bary.push_back(Xc_bary_temp);
              }
              //
              Mat3 H_gc;
              H_gc << 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0;
              for(std::size_t i = 0; i < Xc.size(); ++i)
              {
                  H_gc += Xg_bary[i] * Xc_bary[i].transpose();
              }
              //
              Eigen::JacobiSVD<Mat3> svd_gc(H_gc, Eigen::ComputeFullU | Eigen::ComputeFullV);
              Mat3 s_gc = svd_gc.matrixU();
              Mat3 d_gc = svd_gc.matrixV();
              //
              Mat3 R_gc;
              if(s_gc.determinant()*d_gc.determinant() >= 0)
              {
                  R_gc = s_gc * d_gc.transpose();
              }else{
                  Mat3 A;
                  A << 1.0, 0.0, 0.0,
                       0.0, 1.0, 0.0,
                       0.0, 0.0, -1.0;
                  R_gc = s_gc*A*d_gc.transpose();
              }
              //
              Vec3 t_gc;
              t_gc << 0.0, 0.0, 0.0;
              for(std::size_t i = 0; i < Xc.size(); ++i)
              {
                  Vec3 t_gc_temp;
                  t_gc_temp = Xg[i] - R_gc * Xc[i];
                  t_gc += t_gc_temp;
              }
              t_gc /= Xc.size();
              //
              std::cout << "s: " << s << std::endl;
              std::cout << "R: " << R_gc << std::endl;
              std::cout << "t: " << t_gc << std::endl;
              for (Poses::const_iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose)
//              for(std::size_t poseListId = 0; poseListId != computeImgMap.size(); ++poseListId)
              {
                  const Pose3 & pose = itPose->second;//newData_step4.poses[computeImgMap[poseListId]];
                  const Mat3 R = pose.rotation();
                  const Vec3 t = pose.translation();
                  const Vec3 C = pose.center();
                  Mat3 newR = R * R_gc.inverse();
                  Vec3 newC = R_gc * s * C + t_gc;
                  newSfMData.poses[itPose->first] = Pose3(newR, newC);
              }
              //set all X
              for(Landmarks::const_iterator itX = sfm_data.structure.begin();
                  itX != sfm_data.structure.end(); ++itX)
              {
                  const Vec3 & X = itX->second.X;
                  Vec3 newX = R_gc * s* X + t_gc;
                  newSfMData.structure[itX->first].X = newX;
              }
              Save(newSfMData,
                stlplus::create_filespec(newSfMData.s_root_path, "new_cloud_and_poses"+subPieceIdList[subPListId], ".ply"),
                ESfM_Data(ALL));
//              Save(newSfMData,
//                   stlplus::create_filespec(newSfMData.s_root_path, "new_cloud_and_poses", ".bin"),
//                   ESfM_Data(ALL));

              std::cout << "save newSfMData finish" << std::endl;

              sfm_data = newSfMData;
          }


      }
      sfm_data_list.push_back(sfm_data);
  }
  //
  SfM_Data sfm_data_all;
  sfm_data_all = sfm_data_list[0];
  for(std::size_t subPListId = 1; subPListId < subPieceIdList.size(); ++subPListId)
  {
      SfM_Data thisSub = sfm_data_list[subPListId];
      ///view
      for(std::size_t viewListId = 0; viewListId < thisSub.views.size(); ++viewListId)
      {
          Views::const_iterator itV = thisSub.views.begin();
          std::advance(itV, viewListId);
          //
          sfm_data_all.views[itV->first] = itV->second;
          sfm_data_all.poses[itV->first] = thisSub.poses[itV->first];
      }
      //Landmarks
      Landmarks::const_iterator itS = sfm_data_all.structure.end();
      --itS;
      int lastLId = itS->first + 1;
      for(std::size_t structureListId = 0; structureListId < thisSub.structure.size(); ++structureListId)
      {
          Landmarks::const_iterator itSubS = thisSub.structure.begin();
          std::advance(itSubS, structureListId);
          //
          sfm_data_all.structure[lastLId + structureListId] = itSubS->second;
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

  if (sfmEngine.Process_GCP_GPS_4_margeTest(rtFilePath))//_GCP())//_threeCameras_imu_gps_init_new(rtFilePath))//Process())//Process_threeCameras_imu_gps(sMatchesDir))//Process())
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
