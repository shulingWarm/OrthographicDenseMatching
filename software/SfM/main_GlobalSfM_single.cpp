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


/////!!!
//for (auto itView = sfm_data.views.begin(); itView != sfm_data.views.end(); ++ itView)
//{
//const IndexT imgId = itView->first;
//#pragma omp parallel for schedule(dynamic)
//for(size_t i = 0; i < feats_provider->feats_per_view[imgId].size(); ++i)
//{
//    PointFeatures::iterator itPF = feats_provider->feats_per_view[imgId].begin();
//    std::advance(itPF, i);
//    (*itPF).x() *= 4.0;
//    (*itPF).y() *= 4.0;
//}
//}

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

  if (sfmEngine.Process())//_GCP_GPS(rtFilePath))//_GCP())//_threeCameras_imu_gps_init_new(rtFilePath))//Process())//Process_threeCameras_imu_gps(sMatchesDir))//Process())
  {
      SfM_Data newdata;
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
            stlplus::create_filespec(sOutDir, "sfm_data2", ".json"),
            ESfM_Data(EXTRINSICS));

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

    Save(sfmEngine.Get_SfM_Data(),
      stlplus::create_filespec(sOutDir, "intrinsics", ".json"),
      ESfM_Data(INTRINSICS));

    //////////////////////////////////////////////////////
    needDivideRegion = false;
    if(needDivideRegion)
    {
        std::cout << "need divide region." << std::endl;
        SfM_Data res_sfm_data = sfmEngine.Get_SfM_Data();

        std::vector<std::vector<Vec3>> divideLineList;
        std::vector<int> idList;
        //input divide line
        std::ifstream in;
        in.open(res_sfm_data.s_root_path + "/RegionPiece.res");
        if(in)
        {
            std::size_t pieceCount = 0;
            in >> pieceCount;
            std::size_t line = 0;
            while(line < pieceCount && !in.eof())
            {
                int id;
                double p0[2], p1[2], p2[2], p3[2];
                in >> id >> p0[0] >> p0[1] >> p1[0] >> p1[1] >> p2[0] >> p2[1] >> p3[0] >> p3[1];

                std::vector<Vec3> divide_lines(4);
                double k01 = (p0[1]-p1[1])/(p0[0]-p1[0]);
                double b01 = p0[1] - k01*p0[0];
                double k12 = (p1[1]-p2[1])/(p1[0]-p2[0]);
                double b12 = p1[1] - k12*p1[0];
                double k23 = (p2[1]-p3[1])/(p2[0]-p3[0]);
                double b23 = p2[1] - k23*p2[0];
                double k30 = (p3[1]-p0[1])/(p3[0]-p0[0]);
                double b30 = p3[1] - k30*p3[0];

                divide_lines[0] = Vec3(1.0, -1.0/k01, b01/k01);
                divide_lines[1] = Vec3(1.0, -1.0/k12, b12/k12);
                divide_lines[2] = Vec3(1.0, -1.0/k23, b23/k23);
                divide_lines[3] = Vec3(1.0, -1.0/k30, b30/k30);

                divideLineList.push_back(divide_lines);
                idList.push_back(id);
            }

        }else{
            std::cout << "open file RegionPeice.res failed" << std::endl;
        }

        for(std::size_t i = 0; i < divideLineList.size(); ++i)
        {
            std::cout << "===============" << std::endl;

            //divide region
            SfM_Data sub_sfm_data;
            sub_sfm_data.s_root_path = res_sfm_data.s_root_path;
            sub_sfm_data.intrinsics =  res_sfm_data.intrinsics;

            std::vector<bool> view_map(sfmEngine.Get_SfM_Data().views.size(),false);
            //traverse structure
            for (Landmarks::const_iterator it = res_sfm_data.structure.begin();
                 it != res_sfm_data.structure.end(); ++it)
            {
                Vec3 x = Vec3(it->second.X(0), it->second.X(1), 1.0);
                double r0 = (divideLineList[i][0].transpose()*x);
                double r2 = (divideLineList[i][2].transpose()*x);
                double r_1 = r0 * r2;
                double r1 = (divideLineList[i][1].transpose()*x);
                double r3 = (divideLineList[i][3].transpose()*x);
                double r_2 = r1 * r3;
                if( (r_1 <= 0.0)&&(r_2 <= 0.0) )
                {
                    sub_sfm_data.structure[it->first] = it->second;

                    for(Observations::const_iterator itObs = it->second.obs.begin();
                        itObs != it->second.obs.end(); ++itObs )
                    {
                        //add charge to avoid evaluation-repetition
                        if(view_map[itObs->first] == false) //
                        {
//                            sub_sfm_data.views[itObs->first] = res_sfm_data.views[itObs->first];
//                            sub_sfm_data.poses[itObs->first] = res_sfm_data.poses[itObs->first];
                            view_map[itObs->first] = true;
                            std::cout << itObs->first << std::endl;
                        }
                    }

                }

            }

            Matches_Provider* provider = matches_provider.get();
            const Pair_Set pairs = provider->getPairs();
            Pair_Set res_pair;
            for(Pair_Set::const_iterator it = pairs.begin();
                it != pairs.end(); ++it)
            {
                Pair temp;
                temp.first = it->first;
                temp.second = it->second;
                res_pair.insert(*it);
                res_pair.insert(temp);
            }

            for(std::size_t j = 0; j < view_map.size(); ++j)
            {
                if(view_map[j] == false)
                    continue;

                Pair_Set::const_iterator iter = pairs.begin();
                while(iter != pairs.end())
                {
                    if(iter->first == j)
                        view_map[iter->second] = true;
                    ++iter;
                }
            }

            for(std::size_t k = 0; k < view_map.size(); ++k)
            {
                if(view_map[k] == true)
                {
                    sub_sfm_data.views[k] = res_sfm_data.views[k];
                    sub_sfm_data.poses[k] = res_sfm_data.poses[k];
                }
            }


            //Export
            std::cout << "...Export SfM_Data to disk." << std::endl;
            std::stringstream ss;
            ss << idList[i];
            std::string s;
            ss >> s;

            std::string fileName_bin = "sfm_data_" + s;
            Save(sub_sfm_data,
                 stlplus::create_filespec(sOutDir, fileName_bin, ".bin"),
                 ESfM_Data(ALL));

            std::string fileName_ply = "cloud_and_poses_" + s;
            Save(sub_sfm_data,
                 stlplus::create_filespec(sOutDir, fileName_ply, ".ply"),
                 ESfM_Data(ALL));

        }
    }


    ////////////////////////////////////////////////////////
    //test divide image for densify
    needDivideForDensify = false;
    if(needDivideForDensify)
    {
        std::cout << "need divide for densify." << std::endl;
        int posesCount = 0;
        SfM_Data res_sfm_data = sfmEngine.Get_SfM_Data();
        std::vector<bool> poses_map(sfmEngine.Get_SfM_Data().views.size(),false);
        for(Poses::const_iterator it_pose = res_sfm_data.poses.begin();
            it_pose != res_sfm_data.poses.end(); ++it_pose)
        {
            poses_map[it_pose->first] = true;
            ++posesCount;
        }

        //output file EachCameraInfo.txt
        //for trans and depth-map
        std::ofstream out;
        out.open(res_sfm_data.s_root_path + "/EachCameraInfo.txt");
        out << posesCount << std::endl;

        std::string dir_path = res_sfm_data.s_root_path + "/EachCam";
        ::mkdir(dir_path.c_str(), S_IRWXU | S_IRGRP | S_IXGRP);
    #pragma omp parallel for
        for(std::size_t poseId = 0; poseId < poses_map.size(); ++poseId)
        {
            if(poses_map[poseId] == false)
                continue;

            std::cout << "===============" << std::endl;
            std::vector<bool> view_map(sfmEngine.Get_SfM_Data().views.size(),false);
            std::vector<bool> intrinsic_map(3,false);
            SfM_Data sub_sfm_data;
            sub_sfm_data.s_root_path = res_sfm_data.s_root_path;
            //sub_sfm_data.intrinsics =  res_sfm_data.intrinsics;

            const int resId = poseId;  //aim pic id
            view_map[resId] = true;

            //prepare
            Matches_Provider* provider = matches_provider.get();
            const Pair_Set pairs = provider->getPairs();
            Pair_Set res_pair;
            for(Pair_Set::const_iterator it = pairs.begin();
                it != pairs.end(); ++it)
            {
                Pair temp;
                temp.first = it->second;
                temp.second = it->first;
                res_pair.insert(*it);
                res_pair.insert(temp);
            }

            for(std::size_t j = 0; j < view_map.size(); ++j)
            {
                Pair_Set::const_iterator iter = res_pair.begin();
                while(iter != res_pair.end())
                {
                    if(iter->first == resId)
                        view_map[iter->second] = true;
                    ++iter;
                }
            }

            //traverse structure
            for (Landmarks::const_iterator it = res_sfm_data.structure.begin();
                 it != res_sfm_data.structure.end(); ++it)
            {
                for(std::size_t j = 0; j < view_map.size(); ++j)
                {
                    if(view_map[j] == false)
                        continue;

                    if(it->second.obs.find(j) != it->second.obs.end())
                    {
                        sub_sfm_data.structure[it->first] = it->second;
                    }

                }

            }

            //add camera and traverse camera for instrinsic
            for(std::size_t k = 0; k < view_map.size(); ++k)
            {
                if(view_map[k] == true)
                {
                    //std::cout <<k<< std::endl;
                    sub_sfm_data.views[k] = res_sfm_data.views[k];
                    sub_sfm_data.poses[k] = res_sfm_data.poses[k];
                    const View * view = res_sfm_data.views.at(k).get();
                    intrinsic_map[view->id_intrinsic] = true;
                }
            }

            //add intrinsic
            for(std::size_t k = 0; k < intrinsic_map.size(); ++k)
            {
                if(intrinsic_map[k] == true)
                {
                    //std::cout <<k<< std::endl;
                    sub_sfm_data.intrinsics[k] = res_sfm_data.intrinsics[k];
                }
            }

            //get new id
            int newId = 0;
            for(std::size_t k = 0; k < view_map.size(); ++k)
            {
                if(k == resId)
                    break;
                if(view_map[k] == true)
                    ++newId;
            }


            //Export
            std::stringstream ss;
            ss << resId;
            std::string s;
            ss >> s;
            std::cout << "...Export SfM_Data " << s <<" to disk." << std::endl;

            std::string fileName_bin = "/depth_map_" + s + "_";
            Save(sub_sfm_data,
                 stlplus::create_filespec(dir_path, fileName_bin, ".bin"),
                 ESfM_Data(ALL));

            std::string fileName_ply = "cloud_and_poses_density_" + s;
            Save(sub_sfm_data,
                 stlplus::create_filespec(sOutDir, fileName_ply, ".ply"),
                 ESfM_Data(ALL));

    #pragma omp critical
            {
                out << poseId << " " << newId << " " << dir_path + fileName_bin << std::endl;
            }
        }

        out << "0 0 " << sOutDir + "/sfm_data" << std::endl;
        out.close();
    }


    ////////////////////////////////////////////////////////
    //change in version 20170907
    needDivideGroup = false;
    if(needDivideGroup)
    {
        std::cout << "need divide in to group." << std::endl;
        std::ifstream in_group;
        SfM_Data sfm_data = sfmEngine.Get_SfM_Data();
        in_group.open(sfm_data.s_root_path + "/group.res");
        if(!in_group)
        {
            std::cout << "do not load file group.res, please check it !" << std::endl;
            return EXIT_FAILURE;
        }

        //
        std::list<std::pair<int, std::vector<int> > > groupList; //pair::<downId, <lastId> >

        int groupCount;
        in_group >> groupCount;
        int groupLine = 0;
        while(groupLine < groupCount && !in_group.eof())
        {
            int downId;
            in_group >> downId;
            std::vector<int> restIdList;
            for(std::size_t i = 0; i < 24; ++i)
            {
                int restId;
                in_group >> restId;
                restIdList.push_back(restId);
            }
            std::pair<int, std::vector<int> > groupPair;
            groupPair.first = downId;
            groupPair.second = restIdList;
            groupList.push_back(groupPair);

            ++groupLine;
        }
        in_group.close();


        //prepare data
        int posesCount = 0;
        std::list<bool> posesList(sfm_data.views.size(),false);

        #pragma omp parallel for
        for(std::size_t k = 0; k < sfm_data.poses.size(); ++k)
  //            it_pose != sfm_data.poses.end(); ++it_pose)
        {
            Poses::const_iterator it_pose = sfm_data.poses.begin();
            std::advance(it_pose, k);
            std::list<bool>::iterator it_pList = posesList.begin();
            std::advance(it_pList, it_pose->first);
            *it_pList = true;
            #pragma omp critical
            {
                ++posesCount;
            }

        }

        //prepare folder
        std::string dir_path = sfm_data.s_root_path + "/EachCam";
        ::mkdir(dir_path.c_str(), S_IRWXU | S_IRGRP | S_IXGRP);

        std::ofstream out_groupBoundingBox;
        out_groupBoundingBox.open(sfm_data.s_root_path + "/groupBoundingBox.res");
        if(!out_groupBoundingBox)
        {
            std::cout << "create file groupBoundingBox.res failed !" << std::endl;
            return EXIT_FAILURE;
        }
        out_groupBoundingBox << groupCount << std::endl;

        std::ofstream out_dmapReference;
        out_dmapReference.open(sfm_data.s_root_path + "/dmapReference.res");
        if(!out_dmapReference)
        {
            std::cout << "create file dmapReference.res failed !" << std::endl;
            return EXIT_FAILURE;
        }
        out_dmapReference << groupCount*5 << std::endl;

        for(std::size_t i = 0; i < groupList.size(); ++i)
        {
            std::list<std::pair<int, std::vector<int> > >::iterator it_group = groupList.begin();
            std::advance(it_group, i);

            int downId = it_group->first;

            std::list<bool>::iterator it_pos = posesList.begin();
            std::advance(it_pos, downId);
            if(*it_pos == false)
            {
                continue;
            }

            SfM_Data sub_sfm_data;
            sub_sfm_data.s_root_path = sfm_data.s_root_path;


            double x_min, x_max, y_min, y_max;//bounding box
            bool findFirst = false;

            //traverse structure
            for (Landmarks::const_iterator it = sfm_data.structure.begin();
                 it != sfm_data.structure.end(); ++it)
            {
                Observations::const_iterator it_obs = it->second.obs.find(downId);
                if(it_obs != it->second.obs.end())
                {
  //                  sub_sfm_data.structure[it->first] = it->second;
  //                  sub_sfm_data.structure[it->first].obs[downId] = it->second.obs.find(downId)->second;
                    sub_sfm_data.structure[it->first].obs[downId] = it_obs->second;
                    sub_sfm_data.structure[it->first].X = it->second.X;
                    if(findFirst == false)
                    {
                        findFirst = true;
                        x_min = it->second.X[0];
                        x_max = it->second.X[0];
                        y_min = it->second.X[1];
                        y_max = it->second.X[1];
                    }else{
                        if(it->second.X[0] > x_max)
                        {
                            x_max = it->second.X[0];
                        }else if(it->second.X[0] < x_min)
                        {
                            x_min = it->second.X[0];
                        }

                        if(it->second.X[1] > y_max)
                        {
                            y_max = it->second.X[1];
                        }else if(it->second.X[1] < y_min)
                        {
                            y_min = it->second.X[1];
                        }
                    }
                }

                #pragma omp parallel for
                for(std::size_t j = 0; j < it_group->second.size(); ++j)
                {
                    std::vector<int>::iterator it_restId = it_group->second.begin();
                    std::advance(it_restId, j);
                    Observations::const_iterator it_obs = it->second.obs.find(*it_restId);
                    if(it_obs != it->second.obs.end())
                    {
                        #pragma omp critical
                        {
                            sub_sfm_data.structure[it->first].obs[*it_restId] = it_obs->second;
                            sub_sfm_data.structure[it->first].X = it->second.X;

                        }
  //                      sub_sfm_data.structure[it->first] = it->second;
  //                      sub_sfm_data.structure[it->first].obs[*it_restId] = it->second.obs.find(*it_restId)->second;

                    }
                }


            }

            //save bounding box
            std::stringstream ss_downId, ss_xmin, ss_xmax, ss_ymin, ss_ymax;
            ss_downId << downId;
            ss_xmin << x_min;
            ss_xmax << x_max;
            ss_ymin << y_min;
            ss_ymax << y_max;
            std::string s_downId, s_xmin, s_xmax, s_ymin, s_ymax;
            ss_downId >> s_downId;
            ss_xmin >> s_xmin;
            ss_xmax >> s_xmax;
            ss_ymin >> s_ymin;
            ss_ymax >> s_ymax;

            std::vector<int> dmapIds;//a group id in which need compute dmap
            dmapIds.push_back(2);
            dmapIds.push_back(9);
            dmapIds.push_back(14);
            dmapIds.push_back(21);
            #pragma omp critical
            {
                out_groupBoundingBox << s_downId << " ";
                for(std::size_t j = 0; j < dmapIds.size(); ++j)
                {
                    std::stringstream ss;
                    ss << it_group->second[dmapIds[j]];
                    std::string s;
                    ss >> s;
                    out_groupBoundingBox << s << " ";
                }
                out_groupBoundingBox << std::endl;
                out_groupBoundingBox << s_xmin << " " << s_xmax << " " << s_ymin << " " << s_ymax << std::endl;
            }

            //add camera and traverse camera for instrinsic
            std::vector<bool> intrinsic_map(3,false);
            sub_sfm_data.views[downId] = sfm_data.views[downId];
            sub_sfm_data.poses[downId] = sfm_data.poses[downId];
            const View * view = sfm_data.views.at(downId).get();
            intrinsic_map[view->id_intrinsic] = true;
            //#pragma omp parallel for
            for(std::size_t j = 0; j < it_group->second.size(); ++j)
            {
                std::vector<int>::iterator it_restId = it_group->second.begin();
                std::advance(it_restId, j);
                std::list<bool>::iterator it_pList = posesList.begin();
                std::advance(it_pList, *it_restId);
                if(*it_pList == true)
                {
                    sub_sfm_data.views[*it_restId] = sfm_data.views[*it_restId];
                    sub_sfm_data.poses[*it_restId] = sfm_data.poses[*it_restId];
                    const View * view = sfm_data.views.at(*it_restId).get();
                    intrinsic_map[view->id_intrinsic] = true;
                }
            }

            //add intrinsic
            for(std::size_t k = 0; k < intrinsic_map.size(); ++k)
            {
                if(intrinsic_map[k] == true)
                {
                    sub_sfm_data.intrinsics[k] = sfm_data.intrinsics[k];
                }
            }




            //sorting
            std::vector<int> sortingList = it_group->second;
            sortingList.push_back(it_group->first);
            std::sort(sortingList.begin(),sortingList.end());

            //get new id
            std::vector<int> newIdList; //new id for down 2 9 14 21
            for(std::size_t k = 0; k < sortingList.size(); ++k)  //for down
            {
                if(downId == sortingList[k])
                {
                    newIdList.push_back(k);
                }
            }
            for(std::size_t j = 0; j < 4; ++j)  //for rest
            {
                for(std::size_t k = 0; k < sortingList.size(); ++k)
                {
                    if(it_group->second[dmapIds[j]] == sortingList[k])
                    {
                        newIdList.push_back(k);
                        break;
                    }
                }
            }

            std::vector<std::vector<int> > newRefIdMap; //reference image new id, in order : down, parallel front, back, vertical front, back
            std::vector<int> listId_down; //a group list id in which need to use as reference id for down image
            listId_down.push_back(4);
            listId_down.push_back(5);
            listId_down.push_back(6);
            listId_down.push_back(7);
            std::vector<int> newRefIdList_down;
            for(std::size_t j = 0; j < 4; ++j)
            {
                for(std::size_t k = 0; k < sortingList.size(); ++k)
                {
                    if(it_group->second[listId_down[j]] == sortingList[k])
                    {
                        newRefIdList_down.push_back(k);
                        break;
                    }
                }
            }
            newRefIdMap.push_back(newRefIdList_down);

            std::vector<std::vector<int> > listId_rest(4);//in order : parallel front,back, vetical front,back
            listId_rest[0].push_back(0);
            listId_rest[0].push_back(1);
            listId_rest[0].push_back(3);
            listId_rest[1].push_back(8);
            listId_rest[1].push_back(10);
            listId_rest[1].push_back(11);
            listId_rest[2].push_back(12);
            listId_rest[2].push_back(13);
            listId_rest[2].push_back(15);
            listId_rest[3].push_back(20);
            listId_rest[3].push_back(22);
            listId_rest[3].push_back(23);
            for(std::size_t p = 0; p < 4; p++)
            {
                std::vector<int> newRefIdList_rest;
                for(std::size_t j = 0; j < 3; ++j)
                {
                    for(std::size_t k = 0; k < sortingList.size(); ++k)
                    {
                        if(it_group->second[listId_rest[p][j]] == sortingList[k])
                        {
                            newRefIdList_rest.push_back(k);
                        }
                    }
                }
                newRefIdMap.push_back(newRefIdList_rest);

            }



            //Export
            std::stringstream ss_down;
            ss_down << downId;
            std::string s_down;
            ss_down >> s_down;
            std::string fileName_down = "/depth_map_" + s_down + "_";
            std::string filePath_down = dir_path + "/depth_map_" + s_down + "_";
            std::cout << "...Export SfM_Data " << s_down <<" to disk." << std::endl;

            Save(sub_sfm_data,
                 stlplus::create_filespec(dir_path, fileName_down, ".bin"),
                 ESfM_Data(ALL));
            std::string fileName_ply = "cloud_and_poses_density_" + s_down;
            Save(sub_sfm_data,
                 stlplus::create_filespec(dir_path, fileName_ply, ".ply"),
                 ESfM_Data(ALL));

            //save dmapReference.res
            #pragma omp critical
            {
                std::stringstream ss_newD;
                ss_newD << newIdList[0];
                std::string s_newD;
                ss_newD >> s_newD;

                out_dmapReference << s_down << " " << s_newD << " " << filePath_down << std::endl;
                std::stringstream ss_count;
                ss_count << newRefIdMap[0].size();
                std::string s_count;
                ss_count >> s_count;
                out_dmapReference << s_count << " ";
                for(std::size_t j = 0; j < newRefIdMap[0].size(); ++j)
                {
                    std::stringstream ss;
                    ss << newRefIdMap[0][j];
                    std::string s;
                    ss >> s;
                    out_dmapReference << s << " ";
                }
                out_dmapReference << std::endl;

                for(std::size_t p = 0; p < 4; ++p)
                {
                    std::stringstream ss_old, ss_new;
                    ss_old << it_group->second[dmapIds[p]];
                    ss_new << newIdList[p+1];
                    std::string s_old, s_new;
                    ss_old >> s_old;
                    ss_new >> s_new;
                    out_dmapReference << s_old << " " << s_new << " " << filePath_down << std::endl;
                    std::stringstream ss_count;
                    ss_count << newRefIdMap[p+1].size();
                    std::string s_count;
                    ss_count >> s_count;
                    out_dmapReference << s_count << " ";
                    for(std::size_t k = 0; k < newRefIdMap[p+1].size(); ++k)
                    {
                        std::stringstream ss;
                        ss << newRefIdMap[p+1][k];
                        std::string s;
                        ss >> s;
                        out_dmapReference << s << " ";
                    }
                    out_dmapReference << std::endl;
                }
            }

        }

        out_dmapReference << "0 0 " << sOutDir + "/sfm_data" << std::endl;

        out_groupBoundingBox.close();
        out_dmapReference.close();

    }


    //test water plane
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
        for(std::size_t idMP = 0; idMP < matchesPairs.size(); ++idMP)
        {
            matching::PairWiseMatches::const_iterator itMP = matchesPairs.begin();
            std::advance(itMP, idMP);
            int imgI = itMP->first.first;
            int imgJ = itMP->first.second;
//            std::cout << "imgI : " << imgI << " imgJ : " << imgJ << std::endl;
            //
            /// get features
            std::vector<double> WPlanList;
            std::vector<double> WPlanList_;
            std::vector<double> WPlanList__;
            for(std::size_t idMF = 0; idMF < itMP->second.size(); ++idMF)
            {
                matching::IndMatches::const_iterator itMF = itMP->second.begin();
                std::advance(itMF, idMF);
                //
                int featI = itMF->i_;
                int featJ = itMF->j_;
                //
                PointFeature fI = feats_provider->feats_per_view[imgI][featI];
                PointFeature fJ = feats_provider->feats_per_view[imgJ][featJ];
                //
                // I
                std::shared_ptr<cameras::IntrinsicBase> itIntrinsicI = newdata.intrinsics[newdata.views[imgI]->id_intrinsic];
                std::vector<double> intrinsicI = itIntrinsicI->getParams();
                Pose3 poseI = newdata.poses[imgI];
                Mat3 RI = poseI.rotation();
                Vec3 tI = poseI.translation();
                Vec3 CI = poseI.center();
                Vec3 reprojI = RI.transpose()*(Vec3(fI.x()-intrinsicI[1], fI.x()-intrinsicI[2], intrinsicI[0]) / intrinsicI[0] - tI);
                Vec2 lI = Vec2(sqrt((reprojI(0)-CI(0))*(reprojI(0)-CI(0)) + (reprojI(1)-CI(1))*(reprojI(1)-CI(1))), reprojI(2)-CI(2));
                Vec3 CLI = reprojI-CI;
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
//                Eigen::MatrixXd functionA;
//                functionA.setZero(4,3);
//                functionA << CLI(1), -CLI(0), 0.0,
//                             CLI(2), 0.0, -CLI(0),
//                             CLJ(1), -CLJ(0), 0.0,
//                             CLJ(2), 0.0, -CLJ(0);
//                Eigen::Vector4d function_b;
////                Eigen::VectorXd function_b;
//                function_b << CLI(0)*CI(1)-CLI(1)*CI(0),
//                              CLI(0)*CI(2)-CLI(2)*CI(0),
//                              CLJ(0)*CJ(1)-CLJ(1)*CJ(0),
//                              CLJ(0)*CJ(2)-CLJ(2)*CJ(0);
//                functionA.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(-function_b);
//                Eigen::Vector3d testP = functionA.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(-function_b);
//                const Vec2 P_xy = Vec2(testP(0), testP(1));
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
                std::map<std::pair<int, int>, Vec3>::const_iterator itPI = P_Map.find(std::pair<int, int>(imgI, featI));
                std::map<std::pair<int, int>, Vec3>::const_iterator itPJ = P_Map.find(std::pair<int, int>(imgJ, featJ));
                if(itPI == P_Map.end() || itPJ == P_Map.end())
                {
                    continue;
                }
                const Vec2 P_xy = Vec2(itPI->second(0), itPI->second(1));
                //
                if(itPI->second(0) < -0.366938 || itPI->second(0) > 0.154272 || itPI->second(1) < -1.183584 || itPI->second(1) > -0.388861)
//                if(itP->second(2) < -1.5 || itP->second(2) > -1.2)
                {
                    continue;

                }

                lI(0) = sqrt((P_xy(0)-CI(0))*(P_xy(0)-CI(0)) + (P_xy(1)-CI(1))*(P_xy(1)-CI(1)));
                lJ(0) = sqrt((P_xy(0)-CJ(0))*(P_xy(0)-CJ(0)) + (P_xy(1)-CJ(1))*(P_xy(1)-CJ(1)));





                // I
                double rpI = sqrt((P_xy(0)-CI(0))*(P_xy(0)-CI(0)) + (P_xy(1)-CI(1))*(P_xy(1)-CI(1)));
                double ZpI_ = CI(2) + lI(1)*rpI/lI(0);
                // J
                double rpJ = sqrt((P_xy(0)-CJ(0))*(P_xy(0)-CJ(0)) + (P_xy(1)-CJ(1))*(P_xy(1)-CJ(1)));
                double ZpJ_ = CJ(2) + lJ(1)*rpJ/lJ(0);
                // compute a b c
                double LI = lI(0)/lI(1);
                double LJ = lJ(0)/lJ(1);
                double lIZCI = LI*CI(2);
                double lJZCJ = LJ*CJ(2);
                double a = N*(LI*LI - LJ*LJ);
                double b = 2.0*N*(LJ*(lJZCJ+rpJ) - LI*(lIZCI+rpI)) - 2.0*(ZpI_-ZpJ_);
                double c = N*(lIZCI*(lIZCI+2.0*rpI) - lJZCJ*(lJZCJ+2.0*rpJ) + (rpI+rpJ)*(rpI-rpJ)) + ZpI_*ZpI_ - ZpJ_*ZpJ_;
                //
                double b2_4ac = b*b - 4.0*a*c;
                if(b2_4ac < 0)
                {
                    std::cout << "<0 featI " << featI << " featJ " << featJ << std::endl;
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


//                    if(res1 < itPI->second(2) && res2 < itPI->second(2))
//                    {
//                        continue;
//                    }
//                    if((res1 > res2 ? res2 : res1) < itPI->second(2))// && res2 < itP->second(2))
//                    {
//                        continue;
//                    }
//                    WPlanList.push_back((res1 > res2 ? res1 : res2));
//                    WPlanList_.push_back((res1 > res2 ? res2 : res1));





                }





//                cameras::IntrinsicBase intrinsicI = newdata.intrinsics[newdata.views[imgI].id_intrinsic];


//                cameras::IntrinsicBase intrinsicJ = newdata.intrinsics[newdata.views[imgJ].id_intrinsic];










//                std::cout << "featI : " << featI << " featJ : " << featJ << std::endl;


  //                feats_provider->feats_per_view[]

            }
            //
            WPlanList;
            std::cout << WPlanList.size() << std::endl;
            double mean1 = 0;
            int size1 = 0;
            for(std::size_t wplId = 0; wplId < WPlanList.size(); ++wplId)
            {
//                if(WPlanList[wplId] > 5.0)
//                {
//                    continue;
//                }
                ++size1;
                mean1 += WPlanList[wplId];
            }
            mean1 = mean1 / size1;

            //
            double mean2 = 0;
            int size2 = 0;
            for(std::size_t wplId = 0; wplId < WPlanList_.size(); ++wplId)
            {
//                if(WPlanList_[wplId] < -5.0)
//                {
//                    continue;
//                }
                std::cout << "WPlanList_[wplId] : " << WPlanList_[wplId] << std::endl;
                ++size2;
                mean2 += WPlanList_[wplId];
            }
            mean2 = mean2 / size2;

            //
            std::cout << "mean1 " << mean1 << std::endl;
            std::cout << "mean2 " << mean2 << std::endl;

             std::cout << WPlanList_.size() << std::endl;
            WPlanList_;

        }



    }

    return EXIT_SUCCESS;



  }

  return EXIT_FAILURE;
}
