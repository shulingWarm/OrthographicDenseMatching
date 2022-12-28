// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 cDc <cdc.seacave@gmail.com>, Pierre MOULON

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//change line: 519

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/cameras/Camera_undistort_image.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"

#define _USE_EIGEN
#include "InterfaceMVS.h"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress_display.hpp"

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::image;
using namespace openMVG::sfm;

#include <cstdlib>
#include <string>
#include <sys/stat.h>
#include <list>


#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <ceres/problem.h>
#include <ceres/solver.h>
#include <ceres/types.h>
#include <ceres/cost_function.h>

#include <iomanip>

//ADD_EXECUTABLE(UAVP_main_GetFile main_GetFile.cpp)
//TARGET_LINK_LIBRARIES(UAVP_main_GetFile
//  openMVG_system
//  openMVG_image
//  openMVG_features
//  openMVG_sfm
//  stlplus
//  ${CERES_LIBRARIES}
//  )
//target_include_directories(
//  UAVP_main_GetFile
//  PRIVATE
//  ${CERES_INCLUDE_DIRS})



struct ResidualErrorFunctor
{
  ResidualErrorFunctor(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum : uint8_t {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2,
  };


  template <typename T>
  bool operator()(
    const T* const s,
    const T* const R,
    const T* const t,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
      //--
      // Apply external parameters (Pose)
      //--

      const T * cam_R = R;
      const T * cam_t = t;

      //double a = (double)(angle[0]);
      //std::vector<double> RTrans = getR_angle((angle[0]), (angle[1]), (angle[2]));

      const T& S = s[0];
      T transPFirst[3];
      transPFirst[0] = S * pos_3dpoint[0];
      transPFirst[1] = S*pos_3dpoint[1];
      transPFirst[2] = S*pos_3dpoint[2];

      T pos_proj[3];
      // Rotate the point according the camera rotation
      ceres::AngleAxisRotatePoint(cam_R, transPFirst, pos_proj);

      // Apply the camera translation
      pos_proj[0] += cam_t[0];
      pos_proj[1] += cam_t[1];
      pos_proj[2] += cam_t[2];


      // Compute and return the error is the difference between the predicted
      //  and observed position
      out_residuals[0] = pos_proj[0] - m_pos_2dpoint[0];
      out_residuals[1] = pos_proj[1] - m_pos_2dpoint[1];
      out_residuals[2] = pos_proj[2] - m_pos_2dpoint[2];

      return true;

  }

  static int num_residuals() { return 3; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const std::vector<double> & observation
  )
  {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor, 3, 1, 3, 3, 3>(
            new ResidualErrorFunctor(observation.data())));

  }

  const double * m_pos_2dpoint; // The 2D observation
};


bool settingAndRuning_LS(double& s, Mat3& R, Vec3& t, Poses poses, std::map<int, std::vector<double> > c_gps)
{
    ///prepare optimization computation
    ///
    ///
    ceres::Problem problem;

    std::vector<double> R_a;
    ceres::RotationMatrixToAngleAxis(R.data(), R_a.data());
    double * parameter_block_Ra = R_a.data();
    problem.AddParameterBlock(parameter_block_Ra, 3);


    double * parameter_block_t = t.data();
    problem.AddParameterBlock(parameter_block_t, 3);

    double * parameter_block_s = &s;
    problem.AddParameterBlock(parameter_block_s, 1);


    ///set cost function
    ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(16.0);
    for(Poses::iterator it_pose = poses.begin();
        it_pose != poses.end(); ++it_pose)
    {
        if(it_pose->first / 100000 != 2)
        {
            continue;
        }

        ceres::CostFunction* cost_function =
                ResidualErrorFunctor::Create((c_gps.find(it_pose->first)->second));

        problem.AddResidualBlock(cost_function,
          p_LossFunction,
          &R_a[0],
          &t[0],
          &s,
          &it_pose->second.center()[0]);

        problem.SetParameterBlockConstant(&it_pose->second.center()[0]);

    }

    /// Configure a BA engine and run it
    ///
    {

        //  Make Ceres automatically detect the bundle structure.

        ceres::Solver::Options ceres_config_options;
        ceres_config_options.linear_solver_type = ceres::DENSE_SCHUR;
        ceres_config_options.max_num_iterations = 500;
        ceres_config_options.minimizer_progress_to_stdout = true;
        ceres_config_options.logging_type = ceres::SILENT;


        // Solve BA
        ceres::Solver::Summary summary;
        ceres::Solve(ceres_config_options, &problem, &summary);

        std::cout << summary.BriefReport() << "\n";


        if (!summary.IsSolutionUsable())
        {
          //if (ceres_options_.bVerbose_)
          std::cout << "Bundle Adjustment failed." << std::endl;
          return false;
        }
        else // Solution is usable
        {
          //if (ceres_options_.bVerbose_)
          {
            // Display statistics about the minimization
            std::cout << std::endl
              << "Bundle Adjustment statistics (approximated RMSE):\n"
              << " #residuals: " << summary.num_residuals << "\n"
              << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
              << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
              << " Time (s): " << summary.total_time_in_seconds << "\n"
              << std::endl;

          }
        }
    }

//    std::vector<double> R_a;
    ceres::AngleAxisToRotationMatrix(R_a.data(), R.data());
    return true;

}




struct ImgGroup
{
    std::vector<std::pair<int, std::vector<int> > > groupList;
    std::vector<int> allId;
    int id;
};
typedef std::vector<ImgGroup> ImgGroupList;


int main(int argc, char *argv[])
{
  CmdLine cmd;
  std::string sSfM_Data_Filename;
//  std::string sOutFile = "scene.mvs";
  std::string sOutDir = "undistorted_images";
  std::string sInputFileName = "";
  bool undistortImage = false;
  std::string needDivideForDensify = ""; //damp_x.res
  std::string needDivideGroup = ""; //group_x.res
  bool needRt = false;
  bool needChangeName2Id = false;
  bool needChangeAxie = false;
//  bool demoSituation = false;
  bool needAutoDivideMesh = false;
  int timesAirway = 1;
  std::string changeRootPath = "";
  bool changeSubBlock = false;
  bool newFolder = false;

  int bId = -1;
  cmd.add( make_option('i', sSfM_Data_Filename, "sfm_data.bin") );
//  cmd.add( make_option('o', sOutFile, "outfile") );
  cmd.add( make_option('d', sOutDir, "outdir") );
  cmd.add( make_option('f', sInputFileName, "input file name") );
  cmd.add( make_option('r', needRt, "need out put rt file") );
  cmd.add( make_option('u', undistortImage, "undistort image") );
  cmd.add( make_option('n', needDivideForDensify, "need do dividing for mutli-thread depth map process") );
  cmd.add( make_option('g', needDivideGroup, "need do dividing for mutli-thread depth map process according groups") );
  cmd.add( make_option('c', needChangeName2Id, "need prepare for dividing for mutli-thread depth map process according groups to change image file name to image id") );
  cmd.add( make_option('a', needAutoDivideMesh, "need get divide mesh file automaticlly") );
  cmd.add( make_option('x', needChangeAxie, "need change aeix-z") );
  cmd.add( make_option('b', bId, "block id, same as root_path/block_path") );
  cmd.add( make_option('t', timesAirway, "times of airway step in divide mesh process") );
  cmd.add( make_option('e', changeRootPath, "need change sfm_data_GCP root_path") );
  cmd.add( make_option('j', changeSubBlock, "for sub blocks") );
  cmd.add( make_option('N', newFolder, "need GroupNew and EachCamNew folder") );

//  cmd.add( make_option('s', demoSituation, "need prepare for dividing for mutli-thread depth map process according groups to change image file name to image id, just in demo situation") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch (const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--sfmdata] filename, the SfM_Data file to convert\n"
      << "[-o|--outfile] OpenMVS scene file\n"
      << "[-d|--outdir] undistorted images path\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  // Read the input SfM scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  if(changeRootPath != "")
  {
      sfm_data.s_root_path = changeRootPath;
  }

  for(Poses::const_iterator it = sfm_data.poses.begin();
      it != sfm_data.poses.end(); ++it)
  {
      std::cout << "id : " << it->first << std::endl;
  }



  Save(
      sfm_data,
      "/media/guang/Data/0706_bin/sfm_data_intrinsic.json",
      ESfM_Data(INTRINSICS));
  std::cout << sfm_data.poses.begin()->second.rotation() << std::endl;
  return 0;


//  for(Poses::const_iterator it = sfm_data.poses.begin();
//      it != sfm_data.poses.end(); ++it)
//  {
//      const Pose3 & pose = it->second;
//      const Mat3 R = pose.rotation();
//      const Vec3 t = pose.translation();

//      std::cout << "pose :" << it->first <<std::endl << R <<std::endl << t << std::endl;
//      std::cout << "P :" << it->first <<std::endl << "[" << R(0,0) << "," << R(0,1) << ","<< R(0,2) << "," << t(0) << ";" << std::endl
//                                                  << R(1,0) << "," << R(1,1) << ","<< R(1,2) << "," << t(1) << ";" << std::endl
//                                                  << R(2,0) << "," << R(2,1) << ","<< R(2,2) << "," << t(2) << "]"<<std::endl;

//  }

//  return 0;

//  std::string sOut = "/home/guang/data/20170831/m2/";
//  Save(sfm_data,
//    stlplus::create_filespec(sOut, "sfm_data", ".bin"),
//    ESfM_Data(ALL));

//  Save(sfm_data,
//    stlplus::create_filespec(sOut, "cloud_and_poses", ".ply"),
//    ESfM_Data(ALL));

  //out
  if(needRt)
  {
          std::ofstream out_show_camera;
          out_show_camera.open(sfm_data.s_root_path + "/out_show_camera_points.txt");
          for(Poses::const_iterator it = sfm_data.poses.begin();
              it != sfm_data.poses.end(); ++it)
          {
              const Pose3 & pose = it->second;
              const Mat3 R = pose.rotation();
              const Vec3 c = pose.center();

              Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
              axis_x = R.transpose() * (axis_x + R*c);
              axis_y = R.transpose() * (axis_y + R*c);
              axis_z = R.transpose() * (axis_z + R*c);
              out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
                     << axis_x(0) << " " << axis_x(1) << " " << axis_x(2)<< ";" << std::endl;
              out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
                     << axis_y(0) << " " << axis_y(1) << " " << axis_y(2)<< ";" << std::endl;
              out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
                     << axis_z(0) << " " << axis_z(1) << " " << axis_z(2)<< ";" << std::endl;
          }
          out_show_camera.close();

          std::ofstream out_show_camera_Rt;
          out_show_camera_Rt.open(sfm_data.s_root_path + "/out_show_camera_Rc.txt");
          for(Poses::const_iterator it = sfm_data.poses.begin();
              it != sfm_data.poses.end(); ++it)
          {
              const Pose3 & pose = it->second;
              const Mat3 R = pose.rotation();
              const Vec3 c = pose.center();

              Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
              axis_x = R.transpose() * (axis_x + R*c);
              axis_y = R.transpose() * (axis_y + R*c);
              axis_z = R.transpose() * (axis_z + R*c);
              out_show_camera_Rt <<R(0)<<" "<< R(1)<< " "<<R(2)<< " "
                                 <<R(3)<<" "<< R(4)<< " "<<R(5)<< " "
                                 <<R(6)<<" "<< R(7)<< " "<<R(8)<< " "
                                 <<c(0)<<" "<< c(1)<< " "<<c(2) << std::endl;
          }
          out_show_camera_Rt.close();

  }

  {
//      std::ifstream in_group;//group.res
//      in_group.open(needDivideGroup);
//      if(!in_group)
//      {
//          std::cout << "do not load file group.res, please check it !" << std::endl;
//          return EXIT_FAILURE;
//      }

//      //
//      //std::list<std::pair<int, std::vector<int> > > groupList; //pair::<downId, <lastId> >

//      ImgGroupList imgGroupList;

//      int groupCount;
//      in_group >> groupCount;
//      int groupLine = 0;
//      while(groupLine < groupCount && !in_group.eof())
//      {
//          int image_loacl = 0;
//          ImgGroup imgGroup;
//          std::vector<int> allId;

//          int groupId;
//          in_group >> groupId;
//          imgGroup.id = groupId;

//          int row_count;
//          in_group >> row_count;
//          while(image_loacl < row_count && !in_group.eof())
//          {
//              int imgCount;
//              in_group >> imgCount;

//              int img_local = 1;

//              int thisImageId;
//              in_group >> thisImageId;
//              allId.push_back(thisImageId);
//              std::vector<int> refList;

//              while(img_local < imgCount && !in_group.eof())
//              {
//                  int imgId;
//                  in_group >> imgId;
//                  refList.push_back(imgId);
//                  allId.push_back(imgId);

//                  ++img_local;
//              }

//              std::pair<int, std::vector<int> > tempPair;
//              tempPair.first = thisImageId;
//              tempPair.second = refList;

//              imgGroup.groupList.push_back(tempPair);

//              ++image_loacl;
//          }
//          imgGroup.allId = allId;

//          imgGroupList.push_back(imgGroup);

//          ++groupLine;
//      }
//      in_group.close();


//      SfM_Data sub_sfm_data;
//      sub_sfm_data.s_root_path = sfm_data.s_root_path;
//      std::vector<int> useGroupIdList;
//      useGroupIdList.push_back(3); useGroupIdList.push_back(4); useGroupIdList.push_back(5); useGroupIdList.push_back(6);useGroupIdList.push_back(7);
//      useGroupIdList.push_back(12); useGroupIdList.push_back(13); useGroupIdList.push_back(10); useGroupIdList.push_back(9);useGroupIdList.push_back(8);
//      useGroupIdList.push_back(19); useGroupIdList.push_back(20); useGroupIdList.push_back(21); useGroupIdList.push_back(22);useGroupIdList.push_back(23);
//      useGroupIdList.push_back(28); useGroupIdList.push_back(27); useGroupIdList.push_back(26); useGroupIdList.push_back(25);useGroupIdList.push_back(24);
//      useGroupIdList.push_back(35); useGroupIdList.push_back(36); useGroupIdList.push_back(37); useGroupIdList.push_back(38);useGroupIdList.push_back(39);
////      useGroupIdList.push_back(32); useGroupIdList.push_back(33); useGroupIdList.push_back(34); useGroupIdList.push_back(35);

//      for(std::size_t i = 0; i < useGroupIdList.size(); ++i)
//      {
//          ImgGroupList::iterator it_group = imgGroupList.begin();
//          std::advance(it_group, useGroupIdList[i]);

//          std::vector<int> dmapIds;//a group id in which need compute dmap
//          for(std::size_t j = 0; j < it_group->groupList.size(); ++j)
//          {
//              dmapIds.push_back(it_group->groupList[j].first);
//          }

//          //add camera and traverse camera for instrinsic
//          std::vector<bool> intrinsic_map(3,false);

//          for(std::size_t j = 0; j < it_group->allId.size(); ++j)
//          {
//              std::vector<int>::iterator it_restId = it_group->allId.begin();
//              std::advance(it_restId, j);

//              int imgId = *it_restId;

//              sub_sfm_data.s_root_path = sfm_data.s_root_path;
//              sub_sfm_data.views[imgId] = sfm_data.views[imgId];
//              sub_sfm_data.poses[imgId] = sfm_data.poses[imgId];
//              const View * view = sfm_data.views.at(imgId).get();
//              intrinsic_map[view->id_intrinsic] = true;
//           }

//           //add intrinsic
//           for(std::size_t k = 0; k < intrinsic_map.size(); ++k)
//           {
//               if(intrinsic_map[k] == true)
//               {
//                   sub_sfm_data.intrinsics[k] = sfm_data.intrinsics[k];
//               }
//           }


//          //sorting
//          std::vector<int> sortingList = it_group->allId;
//          std::sort(sortingList.begin(),sortingList.end());

//          sortingList.erase(std::unique(sortingList.begin(), sortingList.end()), sortingList.end());

//      }

//      Save(
//          sub_sfm_data,
//          stlplus::create_filespec( sfm_data.s_root_path, "sfm_data_sub1.json" ).c_str(),
//          ESfM_Data(VIEWS|INTRINSICS));



//      std::cout << "get sub file finished." << std::endl;
//      return 0;

  }

  {

      {
////          sfm_data = sfmEngine.Get_SfM_Data();

//          std::stringstream bIdss;
//          bIdss << bId;
//          std::string bIds;
//          bIdss >> bIds;

//          std::map<int, Vec3> c_gps;
//          std::ifstream in_gps;
//          std::string inGpsPath = sfm_data.s_root_path + "/n_gps_imu_321_" + bIds + ".res";
//          std::cout << "inGpsPath: " << inGpsPath << std::endl;
//          in_gps.open(inGpsPath);
//          if(in_gps)
//          {
//              int count;
//              in_gps >> count;
//              int line = 0;
//              while(line < count && !in_gps.eof())
//              {
//                  int id;
//                  Vec3 xyz, temp;
//                  int tempSatus;

//                  in_gps >> id >> xyz(0) >> xyz(1) >> xyz(2) >> temp(0) >> temp(1) >> temp(2) >> tempSatus;

//                  std::pair<int, Vec3> tempPair;
//                  tempPair.first = id;
//                  tempPair.second = xyz;
//                  c_gps.insert(tempPair);

//                  ++line;
//              }

//          }else{
//              std::cout << "open n_gps_imu_321.res file failed, please check it !" << std::endl;
//              return EXIT_FAILURE;
//          }
//          in_gps.close();

//          double length_gps = 0.0;
//          double length_c = 0.0;
//          sfm::Poses::iterator itEnd = sfm_data.poses.end();
//          -- itEnd;
//          int length_N = 0;
//          for(sfm::Poses::iterator itPose = sfm_data.poses.begin();
//              itPose != itEnd; ++itPose)
//          {
//              int pId1 = itPose->first;
//              if (pId1/100000 != 2)
//              {
//                  continue;
//              }
//              sfm::Poses::iterator itPose2 = itPose;
//              ++ itPose2;
//              for(;itPose2 != sfm_data.poses.end(); ++itPose2)
//              {
//                  int pId2 = itPose2->first;
//                  if (pId2/100000 != 2)
//                  {
//                      continue;
//                  }
//                  ++length_N;
//              }
//          }
//          for(sfm::Poses::iterator itPose = sfm_data.poses.begin();
//              itPose != itEnd; ++itPose)
//          {
//              const Vec3 c1 = itPose->second.center();
//              int pId1 = itPose->first;
//              if (pId1/100000 != 2)
//              {
//                  continue;
//              }
//              const Vec3 cg1 = c_gps.find(pId1)->second;

//              sfm::Poses::iterator itPose2 = itPose;
//              ++ itPose2;
//              for(;itPose2 != sfm_data.poses.end(); ++itPose2)
//              {
//                  const Vec3 c2 = itPose2->second.center();
//                  int pId2 = itPose2->first;
//                  if (pId2/100000 != 2)
//                  {
//                      continue;
//                  }
//                  const Vec3 cg2 = c_gps.find(pId2)->second;

//                  double lg = (cg1-cg2).norm();
//                  double lc = (c1-c2).norm();

//                  length_gps += (lg/(double)(length_N));
//                  length_c += (lc/(double)(length_N));

//              }
//          }

//          double  changeS = length_gps / length_c;

//          //compute H -> Xg = H * Xc
//          std::vector<Vec3> Xc, Xg;
//          for(sfm::Poses::iterator itPose = sfm_data.poses.begin();
//              itPose != sfm_data.poses.end(); ++itPose)
//          {
//              int pId = itPose->first;
//              if (pId/100000 != 2)
//              {
//                  continue;
//              }
//              //prepare Xg
//              Vec3 Xg_temp;
//              Xg_temp = c_gps.find(itPose->first)->second;
//              Xg.push_back(Xg_temp);

//              //prepare Xc
//              Vec3 Xc_temp;
//              Xc_temp = changeS * itPose->second.center();
//              Xc.push_back(Xc_temp);

//          }

//          //prepare Xg
//          Vec3 barycenter_g;
//          barycenter_g << 0.0, 0.0, 0.0;
//          for(std::size_t i = 0; i < Xg.size(); ++i)
//          {
//              barycenter_g += Xg[i];
//          }
//          barycenter_g /= (double)(Xg.size());

//          std::vector<Vec3> Xg_bary;
//          for(std::size_t i = 0; i < Xg.size(); ++i)
//          {
//              Vec3 Xg_bary_temp;
//              Xg_bary_temp = Xg[i] - barycenter_g;
//              Xg_bary.push_back(Xg_bary_temp);
//          }

//          //prepare Xc
//          Vec3 barycenter_c;
//          barycenter_c << 0.0, 0.0, 0.0;
//          for(std::size_t i = 0; i < Xc.size(); ++i)
//          {
//              barycenter_c += Xc[i];
//          }
//          barycenter_c /= (double)(Xc.size());

//          std::vector<Vec3> Xc_bary;
//          for(std::size_t i = 0; i < Xc.size(); ++i)
//          {
//              Vec3 Xc_bary_temp;
//              Xc_bary_temp = Xc[i] - barycenter_c;
//              Xc_bary.push_back(Xc_bary_temp);
//          }

//          Mat3 H_gc;
//          H_gc << 0.0, 0.0, 0.0,
//                  0.0, 0.0, 0.0,
//                  0.0, 0.0, 0.0;
//          for(std::size_t i = 0; i < Xc.size(); ++i)
//          {
//              H_gc += Xg_bary[i] * Xc_bary[i].transpose();
//          }

//          Eigen::JacobiSVD<Mat3> svd_gc(H_gc, Eigen::ComputeFullU | Eigen::ComputeFullV);
//          Mat3 s_gc = svd_gc.matrixU();
//          Mat3 d_gc = svd_gc.matrixV();

//          Mat3 R_gc;
//          if(s_gc.determinant()*d_gc.determinant() >= 0)
//          {
//              R_gc = s_gc * d_gc.transpose();
//          }else{
//              Mat3 A;
//              A << 1.0, 0.0, 0.0,
//                   0.0, 1.0, 0.0,
//                   0.0, 0.0, -1.0;
//              R_gc = s_gc*A*d_gc.transpose();
//          }

//          Vec3 t_gc;
//          t_gc << 0.0, 0.0, 0.0;
//          for(std::size_t i = 0; i < Xc.size(); ++i)
//          {
//              Vec3 t_gc_temp;
//              t_gc_temp = Xg[i] - R_gc * Xc[i];
//              t_gc += t_gc_temp;

//          }
//          t_gc /= Xc.size();

//          //get and set all pose
//    //              Poses::const_iterator itPose = sfm_data.poses.begin();
//    //              std::advance(itPose,)
//          //const Mat3 R0 = sfm_data.poses[imgGroupList[14].allId[0]].rotation();
//          std::cout << "s: " << changeS << std::endl;
//          std::cout << "R: " << R_gc << std::endl;
//          std::cout << "t: " << t_gc << std::endl;

//          //save
//          std::string outGpsPath = sfm_data.s_root_path + "/sRt_" + bIds + ".res";

//          std::ofstream outGps;
//          outGps.open(outGpsPath);
//          if(outGps)
//          {
//              outGps << std::setprecision(15) << changeS << std::endl;
//              outGps << std::setprecision(15) << R_gc << std::endl;
//              outGps << std::setprecision(15) << t_gc << std::endl;

//          }else{
//              std::cout << "create sRt.res file failed." << std::endl;
//              return EXIT_FAILURE;
//          }
//          outGps.close();
//          std::cout << "save sRt file finish" << std::endl;



//          //test
//          {
//              SfM_Data newSfMData = sfm_data;
//              for (Poses::const_iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose)
//              {
//                  const Pose3 & pose = itPose->second;
//                  const Mat3 R = pose.rotation();
//                  const Vec3 t = pose.translation();

//                  Mat3 newR = (1/changeS)* R*R_gc.inverse();
//                  Vec3 newt = t - newR * t_gc;

//                  newSfMData.poses[itPose->first] = Pose3(newR, -newR.inverse() * newt);
//              }

//              //set all X
//              for(Landmarks::const_iterator itX = sfm_data.structure.begin();
//                  itX != sfm_data.structure.end(); ++itX)
//              {
//                  const Vec3 & X = itX->second.X;
//                  Vec3 newX = R_gc * changeS* X + t_gc;

//                  newSfMData.structure[itX->first].X = newX;
//              }


//              Save(newSfMData,
//                stlplus::create_filespec(newSfMData.s_root_path, "/matches/test_sRt_cloud_and_poses", ".ply"),
//                ESfM_Data(ALL));

//              Save(newSfMData,
//                   stlplus::create_filespec(newSfMData.s_root_path, "/matches/test_sRt_cloud_and_poses", ".bin"),
//                   ESfM_Data(ALL));

//          }

      }


  }
  //////////////////////////////////////////////////////
  if(changeSubBlock)
  {

      SfM_Data sub_sfm_data;

      SfM_Data ref_sfm_data;
      std::stringstream ss;
      ss << bId;
      std::string sbId;
      ss >> sbId;
      std::string refPath = changeRootPath + "/sfm_data_"+ sbId + ".json";
      if (!Load(ref_sfm_data, refPath, ESfM_Data(ALL))) {
        std::cerr << std::endl
          << "The input SfM_Data file \""<< refPath << "\" cannot be read." << std::endl;
        return EXIT_FAILURE;
      }
      /// set poses
      sub_sfm_data = ref_sfm_data;
      for(Views::const_iterator itView = ref_sfm_data.views.begin(); itView != ref_sfm_data.views.end(); ++itView)
      {
          sub_sfm_data.poses[itView->first] = sfm_data.poses[itView->first];
          std::cout << "view id : " << itView->first << std::endl;
      }
      /// set structure
      sub_sfm_data.structure = sfm_data.structure;
      /// set intrinsics
      sub_sfm_data.intrinsics = sfm_data.intrinsics;
      //      for(Landmarks::const_iterator itStrecut = ref_sfm_data.views.begin(); itView != ref_sfm_data.views.end(); ++itView)


      sfm_data = sub_sfm_data;
      //sub_sfm_data.poses

      Save(
            sfm_data,
            stlplus::create_filespec( sfm_data.s_root_path+"/matches_"+sbId, "sfm_data_"+sbId+".bin" ).c_str(),
            ESfM_Data(ALL));




  }

  //////////////////////////////////////////////////////
  if(needAutoDivideMesh)
  {

      SfM_Data newSfMData = sfm_data;
      {//break up

//              if(needChangeAxie)
//              {
//    //              std::map<int, Vec3> c_gps;
//    //              std::ifstream in_gps;
//    //              in_gps.open(sfm_data.s_root_path + "/n_gps_imu_321_2.res");
//    //              if(in_gps)
//    //              {
//    //                  int count;
//    //                  in_gps >> count;
//    //                  int line = 0;
//    //                  while(line < count && !in_gps.eof())
//    //                  {
//    //                      int id;
//    //                      Vec3 xyz, temp;
//    //                      int tempSatus;

//    //                      in_gps >> id >> xyz(0) >> xyz(1) >> xyz(2) >> temp(0) >> temp(1) >> temp(2) >> tempSatus;

//    //                      std::pair<int, Vec3> tempPair;
//    //                      tempPair.first = id;
//    //                      tempPair.second = xyz;
//    //                      c_gps.insert(tempPair);

//    //                      ++line;
//    //                  }

//    //              }else{
//    //                  std::cout << "open n_gps_imu_321.res file failed, please check it !" << std::endl;
//    //                  return EXIT_FAILURE;
//    //              }
//    //              in_gps.close();

//    //              Vec3 c_gps0 = c_gps.begin()->second;
//    //              std::cout << "gps begin : " << c_gps0(0) << " "<< c_gps0(1) << " "<< c_gps0(2) << std::endl;
//    //              for(std::map<int, Vec3>::iterator itC = c_gps.begin();
//    //                  itC != c_gps.end(); ++itC)
//    //              {
//    //                  itC->second -= c_gps0;
//    //              }
//    //              std::cout <<  "C gps end " << std::endl;

//    //              double length_gps = 0.0;
//    //              double length_c = 0.0;
//    //              for(std::size_t i = 0; i < imgGroupList.size()-1; ++i)
//    //              {
//    //                  const Vec3 cg1 = c_gps.find(imgGroupList[i].allId[0])->second;//
//    //                  const Vec3 cg2 = c_gps.find(imgGroupList[i+1].allId[0])->second;
//    //                  double lg = cg1.transpose() * cg2;

//    //                  const Vec3 c1 = sfm_data.poses[imgGroupList[i].allId[0]].center();
//    //                  const Vec3 c2 = sfm_data.poses[imgGroupList[i+1].allId[0]].center();
//    //                  double lc = c1.transpose() * c2;

//    //                  length_gps += lg;
//    //                  length_c += lc;

//    //              }
//    //              length_gps /= (double)(imgGroupList.size()-1);
//    //              length_c /= (double)(imgGroupList.size()-1);

//                  double s = 1.0;//0.993516;//length_gps / length_c;

//                  //compute H -> Xg = H * Xc
//                  //prepare Xc
//    //              Vec3 Xc1, Xc2, Xc3;
//    //              Xc1 = s* sfm_data.poses[imgGroupList[0].allId[0]].center();
//    //              Xc2 = s* sfm_data.poses[imgGroupList[14].allId[0]].center();
//    //              Xc3 = s* sfm_data.poses[imgGroupList[21].allId[0]].center();
//    //              Vec3 barycenter_c;
//    //              barycenter_c = (Xc1 + Xc2 + Xc3) / 3.0;
//    //              Vec3 Xc1_bary, Xc2_bary, Xc3_bary;
//    //              Xc1_bary = Xc1 - barycenter_c;
//    //              Xc2_bary = Xc2 - barycenter_c;
//    //              Xc3_bary = Xc3 - barycenter_c;

//    //              //prepare Xg
//    //              Vec3 Xg1, Xg2, Xg3;
//    //              Xg1 = c_gps.find(imgGroupList[0].allId[0])->second;
//    //              Xg2 = c_gps.find(imgGroupList[14].allId[0])->second;
//    //              Xg3 = c_gps.find(imgGroupList[21].allId[0])->second;

//    ////              std::cout << "Xg1 : " << Xg1(0) << " "<< Xg1(1) << " "<< Xg1(2) << std::endl;
//    ////              std::cout << "Xg2 : " << Xg2(0) << " "<< Xg2(1) << " "<< Xg2(2) << std::endl;
//    ////              std::cout << "Xg3 : " << Xg3(0) << " "<< Xg3(1) << " "<< Xg3(2) << std::endl;

//    //              Vec3 barycenter_g;
//    //              barycenter_g = (Xg1 + Xg2 + Xg3) / 3.0;
//    //              Vec3 Xg1_bary, Xg2_bary, Xg3_bary;
//    //              Xg1_bary = Xg1 - barycenter_g;
//    //              Xg2_bary = Xg2 - barycenter_g;
//    //              Xg3_bary = Xg3 - barycenter_g;

//    //              //Xg = H * Xc
//    //              Mat3 H_gc;
//    //              H_gc = Xg1_bary * Xc1_bary.transpose() + Xg2_bary * Xc2_bary.transpose() + Xg3_bary * Xc3_bary.transpose();

//    //              Eigen::JacobiSVD<Mat3> svd_gc(H_gc, Eigen::ComputeFullU | Eigen::ComputeFullV);
//    //              Mat3 s_gc = svd_gc.matrixU();
//    //              Mat3 d_gc = svd_gc.matrixV();

//                  Mat3 R_gc;
//                  R_gc << 1.007046818733, -0.106725797057, 0.030234582722,
//                          0.106586441398, 1.007495999336, 0.006227229256,
//                          -0.030722213909, -0.003008983564, 1.012667179108;
//    //              if(s_gc.determinant()*d_gc.determinant() >= 0)
//    //              {
//    //                  R_gc = s_gc * d_gc.transpose();
//    //              }else{
//    //                  Mat3 A;
//    //                  A << 1.0, 0.0, 0.0,
//    //                       0.0, 1.0, 0.0,
//    //                       0.0, 0.0, -1.0;
//    //                  R_gc = s_gc*A*d_gc.transpose();
//    //              }

//                  Vec3 t_gc;
//    //              t_gc = Xg1 - R_gc * Xc1;
//                  t_gc << 973.55322265625, 2700.741943359375, 1877.64501953125;

//    //              std::cout << "R_gc : " << std::endl << R_gc << std::endl << "t_gc : " << t_gc(0) << " " << t_gc(1) << " " << t_gc(2) << std::endl;

//                  //use least square
//                  {
//    //                  std::map<int, std::vector<double> >c_gps_temp;
//    //                  for(std::map<int, Vec3>::iterator itC = c_gps.begin();
//    //                      itC != c_gps.end(); ++itC)
//    //                  {
//    //                      std::pair<int, std::vector<double> >tempPair;
//    //                      tempPair.first = itC->first;
//    //                      tempPair.second.push_back(itC->second[0]);
//    //                      tempPair.second.push_back(itC->second[1]);
//    //                      tempPair.second.push_back(itC->second[2]);
//    //                      c_gps_temp.insert(tempPair);

//    //                  }
//    //                  settingAndRuning_LS(s, R_gc, t_gc, sfm_data.poses, c_gps_temp);

//    //                  for (const auto& view : sfm_data.GetViews())
//    //                  {
//    //                      if(sfm_data.poses.find(view->id_pose) == poses.end() && view->id_pose/100000 ==2)
//    //                      {
//    //                          sfm_data.poses[view->id_pose] = SfM::Pose (c_gps.find(view->id_pose))



//    //                      }
//    //                  }

//                  }


//                  //get and set all pose
//    //              Poses::const_iterator itPose = sfm_data.poses.begin();
//    //              std::advance(itPose,)
//                  const Mat3 R0 = sfm_data.poses[imgGroupList[14].allId[0]].rotation();
//                  for (Poses::const_iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose)
//                  {
//                      const Pose3 & pose = itPose->second;
//                      const Mat3 R = pose.rotation();
//                      const Vec3 t = pose.translation();
//                      const Vec3 C = pose.center();

//                      Mat3 newR = (1.0/s)*R*R_gc.inverse();
//    //                  Mat3 newR = R*R_gc.inverse();
//    //                  Mat3 newR = R*(s*R_gc).transpose();
//    //                  Vec3 newt = t - (1.0/s)*newR * t_gc;
//                      Vec3 newt = t - newR * t_gc;

//                      Vec3 newC = - newR.transpose() * newt;
//                      Vec3 testNewC = s * R_gc * C + t_gc;

//    //                  newSfMData.poses[itPose->first] = Pose3(newR, -newR.inverse() * newt);
//                      newSfMData.poses[itPose->first] = Pose3(newR, -newR.inverse() * newt);
//    //                  newSfMData.poses[itPose->first] = Pose3(newR, -newR.inverse() * newt);
//    //                  newSfMData.poses[itPose->first] = Pose3(newR, testNewC);
//    //                  newSfMData.poses[itPose->first] = Pose3(R*R0.inverse(), R0.inverse()*C);
//    //                  newSfMData.poses[itPose->first] = Pose3(R, testNewC);

//                      Vec3 testNewC2 = s * R_gc * C;
//    //                  Vec3 testNewC2 = s * C;
//    //                  std::cout <<  testNewC2(0) << " " << testNewC2(1) << " " << testNewC2(2) << " 0 255 0 255" << std::endl;

//    //                  std::cout << itPose->first << std::endl;
//    //                  std::cout << newR << std::endl << "7" << std::endl;
//    //                  std::cout << newSfMData.poses[itPose->first].rotation() << std::endl<<"6" << std::endl;
//    //                  std::cout << C(0) << " " << C(1) << " " << C(2) << std::endl;
//                  }

//                  //set all X
//                  for(Landmarks::const_iterator itX = sfm_data.structure.begin();
//                      itX != sfm_data.structure.end(); ++itX)
//                  {
//                      const Vec3 & X = itX->second.X;

//                      Vec3 newX = s * R_gc * X + t_gc;
//    //                  Vec3 newX = R_gc * X + t_gc;
//    //                  Vec3 newX = R0.inverse() * X;
//                      newSfMData.structure[itX->first].X = newX;
//                  }

//                  Save(newSfMData,
//                    stlplus::create_filespec(newSfMData.s_root_path, "new_cloud_and_poses", ".ply"),
//                    ESfM_Data(ALL));


//    //              sfm_data = newSfMData;

//                  std::cout << "save newSfMData finish" << std::endl;





//              }


      }

      std::ifstream in_group;//group.res
      in_group.open(needDivideGroup);
      if(!in_group)
      {
          std::cout << "do not load file group.res, please check it !" << std::endl;
          return EXIT_FAILURE;
      }

      ImgGroupList imgGroupList;

      int groupCount;
      in_group >> groupCount;
      int groupLine = 0;
      while(groupLine < groupCount && !in_group.eof())
      {
          int image_loacl = 0;
          ImgGroup imgGroup;
          std::vector<int> allId;

          int groupId;
          in_group >> groupId;
          imgGroup.id = groupId;

          int row_count;
          in_group >> row_count;
          while(image_loacl < row_count && !in_group.eof())
          {
              int imgCount;
              in_group >> imgCount;

              int img_local = 1;

              int thisImageId;
              in_group >> thisImageId;
              allId.push_back(thisImageId);
              std::vector<int> refList;

              while(img_local < imgCount && !in_group.eof())
              {
                  int imgId;
                  in_group >> imgId;
                  refList.push_back(imgId);
                  allId.push_back(imgId);

                  ++img_local;
              }

              std::pair<int, std::vector<int> > tempPair;
              tempPair.first = thisImageId;
              tempPair.second = refList;

              imgGroup.groupList.push_back(tempPair);

              ++image_loacl;
          }
          imgGroup.allId = allId;

          imgGroupList.push_back(imgGroup);

          ++groupLine;
      }
      in_group.close();

      if(needChangeAxie)
      {
          std::stringstream bIdss;
          bIdss << bId;
          std::string bIds;
          bIdss >> bIds;

          std::map<int, Vec3> c_gps;
          std::ifstream in_gps;
          std::string inGpsPath = sfm_data.s_root_path + "/n_gps_imu_321_" + bIds + ".res";
          std::cout << "inGpsPath: " << inGpsPath << std::endl;
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

          Vec3 c_gps0 = c_gps.begin()->second;
          std::ofstream out_origC;
          out_origC.open(sfm_data.s_root_path + "/origC.txt");
          std::cout << "gps begin : " << c_gps0(0) << " "<< c_gps0(1) << " "<< c_gps0(2) << std::endl;
          for(std::map<int, Vec3>::iterator itC = c_gps.begin();
              itC != c_gps.end(); ++itC)
          {
              itC->second -= c_gps0;
//                  /!!!!!!!!!!!!!!!!!
              out_origC << itC->second(0) << " " << itC->second(1) << " " << itC->second(2) << std::endl;

          }
          out_origC.close();

          std::cout <<  "C gps end " << std::endl;

          double length_gps = 0.0;
          double length_c = 0.0;
          sfm::Poses::iterator itEnd = sfm_data.poses.end();
          -- itEnd;
          int length_N = 0;
          for(sfm::Poses::iterator itPose = sfm_data.poses.begin();
              itPose != itEnd; ++itPose)
          {
              const Vec3 c1 = itPose->second.center();
              int pId1 = itPose->first;
              if (pId1/100000 != 2)
              {
                  continue;
              }
              const Vec3 cg1 = c_gps.find(pId1)->second;

              sfm::Poses::iterator itPose2 = itPose;
              ++ itPose2;
              for(;itPose2 != sfm_data.poses.end(); ++itPose2)
              {
                  const Vec3 c2 = itPose2->second.center();
                  int pId2 = itPose2->first;
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
//              double s = 1.0;//0.993516;//length_gps / length_c;


          //compute H -> Xg = H * Xc
          std::vector<Vec3> Xc, Xg;
          for(sfm::Poses::iterator itPose = sfm_data.poses.begin();
              itPose != sfm_data.poses.end(); ++itPose)
          {
              int pId = itPose->first;
              if (pId/100000 != 2)
              {
                  continue;
              }
              //prepare Xg
              Vec3 Xg_temp;
              Xg_temp = c_gps.find(itPose->first)->second;
              Xg.push_back(Xg_temp);

              //prepare Xc
              Vec3 Xc_temp;
              Xc_temp = s * itPose->second.center();
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

          Mat3 H_gc;
          H_gc << 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0;
          for(std::size_t i = 0; i < Xc.size(); ++i)
          {
              H_gc += Xg_bary[i] * Xc_bary[i].transpose();
          }

          Eigen::JacobiSVD<Mat3> svd_gc(H_gc, Eigen::ComputeFullU | Eigen::ComputeFullV);
          Mat3 s_gc = svd_gc.matrixU();
          Mat3 d_gc = svd_gc.matrixV();

          Mat3 R_gc;
//              R_gc << 1.007046818733, -0.106725797057, 0.030234582722,
//                      0.106586441398, 1.007495999336, 0.006227229256,
//                      -0.030722213909, -0.003008983564, 1.012667179108;
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

          Vec3 t_gc;
          t_gc << 0.0, 0.0, 0.0;
          for(std::size_t i = 0; i < Xc.size(); ++i)
          {
              Vec3 t_gc_temp;
              t_gc_temp = Xg[i] - R_gc * Xc[i];
              t_gc += t_gc_temp;

          }
          t_gc /= Xc.size();

//              t_gc << 973.55322265625, 2700.741943359375, 1877.64501953125;

//              std::cout << "R_gc : " << std::endl << R_gc << std::endl << "t_gc : " << t_gc(0) << " " << t_gc(1) << " " << t_gc(2) << std::endl;


          //get and set all pose
//              Poses::const_iterator itPose = sfm_data.poses.begin();
//              std::advance(itPose,)
          //const Mat3 R0 = sfm_data.poses[imgGroupList[14].allId[0]].rotation();
          std::cout << "s: " << s << std::endl;
          std::cout << "R: " << R_gc << std::endl;
          std::cout << "t: " << t_gc << std::endl;
          for (Poses::const_iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose)
          {
              const Pose3 & pose = itPose->second;
              const Mat3 R = pose.rotation();
              const Vec3 t = pose.translation();
              const Vec3 C = pose.center();

              Mat3 newR = (1/s)* R*R_gc.inverse();
//                  Mat3 newR = R*R_gc.inverse();
//                  Mat3 newR = R*(s*R_gc).transpose();
//                  Vec3 newt = t - (1.0/s)*newR * t_gc;
              Vec3 newt = t - newR * t_gc;

              Vec3 newC = - newR.transpose() * newt;
              Vec3 testNewC = R_gc * C + t_gc;

//                  newSfMData.poses[itPose->first] = Pose3(newR, -newR.inverse() * newt);
              newSfMData.poses[itPose->first] = Pose3(newR, -newR.inverse() * newt);
//                  newSfMData.poses[itPose->first] = Pose3(newR, -newR.inverse() * newt);
//                  newSfMData.poses[itPose->first] = Pose3(newR, testNewC);
//                  newSfMData.poses[itPose->first] = Pose3(R*R0.inverse(), R0.inverse()*C);
//                  newSfMData.poses[itPose->first] = Pose3(R, testNewC);

              Vec3 testNewC2 = s * R_gc * C;
//                  Vec3 testNewC2 = s * C;
//                  std::cout <<  testNewC2(0) << " " << testNewC2(1) << " " << testNewC2(2) << " 0 255 0 255" << std::endl;

//                  std::cout << itPose->first << std::endl;
//                  std::cout << newR << std::endl << "7" << std::endl;
//                  std::cout << newSfMData.poses[itPose->first].rotation() << std::endl<<"6" << std::endl;
//                  std::cout << C(0) << " " << C(1) << " " << C(2) << std::endl;
          }

          //set all X
          for(Landmarks::const_iterator itX = sfm_data.structure.begin();
              itX != sfm_data.structure.end(); ++itX)
          {
              const Vec3 & X = itX->second.X;

              Vec3 newX = R_gc * s* X + t_gc;
//                  Vec3 newX = R_gc * X + t_gc;
//                  Vec3 newX = R0.inverse() * X;
              newSfMData.structure[itX->first].X = newX;
          }

          Save(newSfMData,
            stlplus::create_filespec(newSfMData.s_root_path, "new_cloud_and_poses", ".ply"),
            ESfM_Data(ALL));


          Save(newSfMData,
               stlplus::create_filespec(newSfMData.s_root_path, "new_cloud_and_poses", ".bin"),
               ESfM_Data(ALL));


//              sfm_data = newSfMData;

          std::cout << "save newSfMData finish" << std::endl;





      }
      sfm_data = newSfMData;
      //////////////////////////////////////////////////////


      //set reference direction
      const Vec3 c0 = sfm_data.poses[imgGroupList[0].allId[0]].center();
      const Vec3 c1 = sfm_data.poses[imgGroupList[1].allId[0]].center();
      double ref_dir_x = c1(0) - c0(0);
      double ref_dir_y = c1(1) - c0(1);

      //find and compute length when divide region
      std::vector<int> downIdList;
      //std::vector<std::pair<int, int> > anotherLineList; //prepare another id in  airway to get region length
      std::vector<std::pair<double, double> > lengthList; //prepare another id in  airway to get region length

      std::vector<double> tempLenList;
      std::vector<int> tempLineNumList;
      int tempNum = 1;
      downIdList.push_back(imgGroupList[0].allId[0]);
      bool changeNow = false;
      for(std::size_t i = 1; i < imgGroupList.size(); ++i)
//              for(std::size_t i = 1; i < 8; ++i)
      {
          downIdList.push_back(imgGroupList[i].allId[0]);

          const Vec3 c0 = sfm_data.poses[downIdList[i-1]].center();
          const Vec3 c1 = sfm_data.poses[downIdList[i]].center();

          double dir_x = c1(0) - c0(0);
          double dir_y = c1(1) - c0(1);

          double angle = acos((ref_dir_x*dir_x + ref_dir_y*dir_y) /
                              (sqrt(ref_dir_x*ref_dir_x + ref_dir_y*ref_dir_y) * sqrt(dir_x*dir_x + dir_y*dir_y)));
//              std::cout << "angle : " << angle << std::endl;

          double a_ = angle - 1.5707963;//0.785;
          if(a_ < 0.0)// || angle > 2.355)
          {
              ++ tempNum;

              //change reference dir
              ref_dir_x = dir_x;
              ref_dir_y = dir_y;

          }else{
              //change to next line
              tempLineNumList.push_back(tempNum);

              std::vector<double> forMinL;
              for(std::size_t j = 0; j < tempNum; ++j)
              {
                  const Vec3 cl = sfm_data.poses[downIdList[i-j-1]].center();
                  double tempLen = sqrt((c1(0)-cl(0))*(c1(0)-cl(0)) + (c1(1)-cl(1))*(c1(1)-cl(1)));

                  forMinL.push_back(tempLen);
              }

              double findtempLen = *(std::min_element(std::begin(forMinL), std::end(forMinL)) );//sqrt((c1(0)-c0(0))*(c1(0)-c0(0)) + (c1(1)-c0(1))*(c1(1)-c0(1)));
              tempLenList.push_back(findtempLen);

              tempNum = 1;
              changeNow = true;

              //change reference dir
              //Vec3 c_now = c1;
              if(i+2 < imgGroupList.size())
              {
                  Vec3 c_next = sfm_data.poses[imgGroupList[i+1].allId[0]].center();

                  ref_dir_x = c_next(0) - c1(0);
                  ref_dir_y = c_next(1) - c1(1);
              }

          }

      }
      {
          //for last line
          tempLineNumList.push_back(tempNum);
          tempLenList.push_back(tempLenList[tempLenList.size()-1]);
      }

      std::cout << tempLenList.size() << std::endl<< std::endl;


      std::vector<bool> tempWidList;
      {
          double lastLen = tempLenList[0];
          for(std::size_t j = 0; j < tempLineNumList.size(); ++j)
          {
              for(std::size_t k = 0; k < tempLineNumList[j]; ++k)
              {
                  std::pair<double, double> tempLengthPair;
                  tempLengthPair.first = lastLen;
                  tempLengthPair.second = tempLenList[j];
                  lengthList.push_back(tempLengthPair);

                  tempWidList.push_back(false);

              }

              lastLen = tempLenList[j];
              tempWidList[tempWidList.size() - 1] = true;
          }
      }

      //////find length
      ///
      bool findL = false;
      for(std::size_t i = 1; i < downIdList.size(); ++i)
//              for(std::size_t i = 1; i < 8; ++i)
      {
          const Vec3 c0 = sfm_data.poses[downIdList[i-1]].center();
          const Vec3 c1 = sfm_data.poses[downIdList[i]].center();
          const double dir_x = c1(0) - c0(0);
          const double dir_y = c1(1) - c0(1);

          double angle = acos((ref_dir_x*dir_x + ref_dir_y*dir_y) /
                              (sqrt(ref_dir_x*ref_dir_x + ref_dir_y*ref_dir_y) * sqrt(dir_x*dir_x + dir_y*dir_y)));

          //if(angle < 0.785)// || angle > 2.355)// angle < 45 or angle > 135 -> compute with last
          double a_ = angle - 1.5707963;//0.785;
          if(a_ > 0.0)
          {
              findL = true;
              const Vec3 cl = sfm_data.poses[downIdList[0]].center();
              double tempLen = sqrt((c1(0)-cl(0))*(c1(0)-cl(0)) + (c1(1)-cl(1))*(c1(1)-cl(1)));
              lengthList[0].first = tempLen;
              break;


          }



      }
      std::cout << "lengthList : " << lengthList[0].first << std::endl;

      ref_dir_x = c1(0) - c0(0);
      ref_dir_y = c1(1) - c0(1);

      std::vector<std::vector<double> > resList;
      std::vector<std::vector<double> > resOverlappingList;
      //for c0
      {
          double l0_A, l0_B, l0_C, l1_A, l1_B, l1_C, l2_A, l2_B, l2_C, l3_A, l3_B, l3_C;
          //const Vec3 cW1 = sfm_data.poses[anotherLineList[0].first].center();
          //const Vec3 cW2 = sfm_data.poses[anotherLineList[0].second].center();
          double paral_k = ref_dir_y/ref_dir_x;
          double width = sqrt((c1(0)-c0(0))*(c1(0)-c0(0)) + (c1(1)-c0(1))*(c1(1)-c0(1)));
          double length = std::max(lengthList[0].first, lengthList[0].second);//sqrt((cW2(0)-cW1(0))*(cW2(0)-cW1(0)) + (cW2(1)-cW1(1))*(cW2(1)-cW1(1)));
          double angleK_paral = atan(paral_k);
          double mid_paralB = c0(1) - c0(0)*paral_k;
          l0_A = 1.0;
          l0_B = -1.0/paral_k;
          l0_C = (mid_paralB + length*5.0/8.0 / cos(angleK_paral)) / paral_k;
          l2_A = 1.0;
          l2_B = l0_B;
          l2_C = (mid_paralB - length*5.0/8.0 / cos(angleK_paral)) / paral_k;
          double ver_k = -1.0/paral_k;
          double angleK_ver = atan(ver_k);
          double mid_verB = c0(1) - c0(0)*ver_k;
          l1_A = 1.0;
          l1_B = -1.0/ver_k;
          l1_C = (mid_verB + width*5.0/8.0 / cos(angleK_ver)) / ver_k;
          l3_A = 1.0;
          l3_B = l1_B;
          l3_C = (mid_verB - width*5.0/8.0 / cos(angleK_ver)) / ver_k;

//              std::cout << "length : " << length << std::endl;

          //record
          std::vector<double> res(12);
          res[0] = l0_A; res[1] = l0_B; res[2] = l0_C;
          res[3] = l1_A; res[4] = l1_B; res[5] = l1_C;
          res[6] = l2_A; res[7] = l2_B; res[8] = l2_C;
          res[9] = l3_A; res[10] = l3_B; res[11] = l3_C;
          resList.push_back(res);

          double l0_A_over = 1.0, l0_B_over = -1.0/paral_k, l0_C_over = (mid_paralB + length / cos(angleK_paral)) / paral_k;
          double l1_A_over = 1.0, l1_B_over = -1.0/ver_k, l1_C_over = (mid_verB + width / cos(angleK_ver)) / ver_k;
          double l2_A_over = 1.0, l2_B_over = l0_B_over, l2_C_over = (mid_paralB - length / cos(angleK_paral)) / paral_k;
          double l3_A_over = 1.0, l3_B_over = l1_B_over, l3_C_over = (mid_verB - width / cos(angleK_ver)) / ver_k;

          std::vector<double> res_over(12);
          res_over[0] = l0_A_over; res_over[1] = l0_B_over; res_over[2] = l0_C_over;
          res_over[3] = l1_A_over; res_over[4] = l1_B_over; res_over[5] = l1_C_over;
          res_over[6] = l2_A_over; res_over[7] = l2_B_over; res_over[8] = l2_C_over;
          res_over[9] = l3_A_over; res_over[10] = l3_B_over; res_over[11] = l3_C_over;
          resOverlappingList.push_back(res_over);

      }


      //record result
      //0-----1  0-----1
      //  c0   ->  c1
      //2-----3  2-----3
      //
      //---l0---    ---l0---
      //l3 c0 l1 -> l3  c1 l1   |
      //---l2---    ---l2---
      const Vec3 C0 = sfm_data.poses[downIdList[0]].center();
      for(std::size_t i = 1; i < downIdList.size(); ++i)
//              for(std::size_t i = 1; i < 8; ++i)
      {
          const Vec3 c0 = sfm_data.poses[downIdList[i-1]].center();
          const Vec3 c1 = sfm_data.poses[downIdList[i]].center();
          const double dir_x = c1(0) - c0(0);
          const double dir_y = c1(1) - c0(1);

          double angle = acos((ref_dir_x*dir_x + ref_dir_y*dir_y) /
                              (sqrt(ref_dir_x*ref_dir_x + ref_dir_y*ref_dir_y) * sqrt(dir_x*dir_x + dir_y*dir_y)));
//              std::cout << "angle : " << angle << std::endl;

          double paral_k;
          double width;
          //if(angle < 0.785)// || angle > 2.355)// angle < 45 or angle > 135 -> compute with last
          double a_ = angle - 1.5707963;//0.785;
          if(a_ < 0.0)// || angle > 2.355)
          {

              if(i+1 < downIdList.size() && !tempWidList[i])
              {
                  const Vec3 c2 = sfm_data.poses[downIdList[i+1]].center();
                  double width10 = sqrt((c1(0)-c0(0))*(c1(0)-c0(0)) + (c1(1)-c0(1))*(c1(1)-c0(1)));
                  double width21 = sqrt((c2(0)-c1(0))*(c2(0)-c1(0)) + (c2(1)-c1(1))*(c2(1)-c1(1)));
                  width = std::max(width10, width21);

              }else{
                  width = sqrt((c1(0)-c0(0))*(c1(0)-c0(0)) + (c1(1)-c0(1))*(c1(1)-c0(1)));

              }

              paral_k = dir_y/dir_x;
              //change reference dir
              ref_dir_x = dir_x;
              ref_dir_y = dir_y;
          }
          else{ //-> compute with next
              const Vec3 c2 = sfm_data.poses[downIdList[i+1]].center();
              const double dir_x_21 = c2(0) - c1(0);
              const double dir_y_21 = c2(1) - c1(1);

              paral_k = dir_y_21 / dir_x_21;
              width = sqrt((c2(0)-c1(0))*(c2(0)-c1(0)) + (c2(1)-c1(1))*(c2(1)-c1(1)));

              //change reference dir
              Vec3 c_now = c1;
              if(i+1 < downIdList.size())
              {
//                      const Vec3 c2 = sfm_data.poses[downIdList[i+1]].center();
//                      double width10 = sqrt((c1(0)-c0(0))*(c1(0)-c0(0)) + (c1(1)-c0(1))*(c1(1)-c0(1)));
//                      double width21 = sqrt((c2(0)-c1(0))*(c2(0)-c1(0)) + (c2(1)-c1(1))*(c2(1)-c1(1)));
//                      width = std::max(width10, width21);

                  Vec3 c_next = sfm_data.poses[downIdList[i+1]].center();
                  ref_dir_x = c_next(0) - c_now(0);
                  ref_dir_y = c_next(1) - c_now(1);

              }
          }
          double l0_A, l0_B, l0_C, l1_A, l1_B, l1_C, l2_A, l2_B, l2_C, l3_A, l3_B, l3_C;
          double ver_k = -1.0/paral_k;
          double angleK_ver = atan(ver_k);
          double mid_verB = c1(1) - c1(0)*ver_k;
          l1_A = 1.0;
          l1_B = -1.0/ver_k;
          l1_C = (mid_verB + width*5.0/8.0 / cos(angleK_ver)) / ver_k;
          l3_A = 1.0;
          l3_B = l1_B;
          l3_C = (mid_verB - width*5.0/8.0 / cos(angleK_ver)) / ver_k;


          //const Vec3 cW1 = sfm_data.poses[anotherLineList[i].first].center();
          //const Vec3 cW2 = sfm_data.poses[anotherLineList[i].second].center();
          double length = std::max(lengthList[0].first, lengthList[0].second);//std::max(lengthList[i].first, lengthList[i].second );//sqrt((cW2(0)-cW1(0))*(cW2(0)-cW1(0)) + (cW2(1)-cW1(1))*(cW2(1)-cW1(1)));
          double angleK_paral = atan(paral_k);
          double mid_paralB = c1(1) - c1(0)*paral_k;

          {
//                  //l2-second
//                  double length_s = lengthList[i].second;
//                  double l2_A_t1 = 1.0;
//                  double l2_B_t1 = -1.0/paral_k;
//                  double l2_C_t1 = (mid_paralB + length_s*9.0/16.0 / cos(angleK_paral)) / paral_k;
//                  //l2-second
//                  double l2_A_t2 = 1.0;
//                  double l2_B_t2 = -1.0/paral_k;
//                  double l2_C_t2 = (mid_paralB - length_s*9.0/16.0 / cos(angleK_paral)) / paral_k;

//                  double d_c0_l2_t1 = abs(l2_A_t1*C0(0) + l2_B_t1*C0(1) + l2_C_t1) / sqrt(l2_A_t1*l2_A_t1 + l2_B_t1*l2_B_t1);
//                  double d_c0_l2_t2 = abs(l2_A_t2*C0(0) + l2_B_t2*C0(1) + l2_C_t2) / sqrt(l2_A_t2*l2_A_t2 + l2_B_t2*l2_B_t2);;
//                  if(d_c0_l2_t1 > d_c0_l2_t2)
//                  {
//                      l2_A = l2_A_t1;
//                      l2_B = l2_B_t1;
//                      l2_C = l2_C_t1;
//                  }else{
//                      l2_A = l2_A_t2;
//                      l2_B = l2_B_t2;
//                      l2_C = l2_C_t2;
//                  }

//                  //l0-first
//                  double length_f = lengthList[i].first;
//                  double l0_A_t1 = 1.0;
//                  double l0_B_t1 = -1.0/paral_k;
//                  double l0_C_t1 = (mid_paralB + length_f*9.0/16.0 / cos(angleK_paral)) / paral_k;
//                  //l0-first
//                  double l0_A_t2 = 1.0;
//                  double l0_B_t2 = -1.0/paral_k;
//                  double l0_C_t2 = (mid_paralB - length_f*9.0/16.0 / cos(angleK_paral)) / paral_k;

//                  double d_c0_l0_t1 = abs(l0_A_t1*C0(0) + l0_B_t1*C0(1) + l0_C_t1) / sqrt(l0_A_t1*l0_A_t1 + l0_B_t1*l0_B_t1);
//                  double d_c0_l0_t2 = abs(l0_A_t2*C0(0) + l0_B_t2*C0(1) + l0_C_t2) / sqrt(l0_A_t2*l0_A_t2 + l0_B_t2*l0_B_t2);;
//                  if(d_c0_l0_t1 < d_c0_l0_t2)
//                  {
//                      l0_A = l0_A_t1;
//                      l0_B = l0_B_t1;
//                      l0_C = l0_C_t1;
//                  }else{
//                      l0_A = l0_A_t2;
//                      l0_B = l0_B_t2;
//                      l0_C = l0_C_t2;
//                  }


          }



          l0_A = 1.0;
          l0_B = -1.0/paral_k;
          l0_C = (mid_paralB + length*5.0/8.0 / cos(angleK_paral)) / paral_k;
          l2_A = 1.0;
          l2_B = l0_B;
          l2_C = (mid_paralB - length*5.0/8.0 / cos(angleK_paral)) / paral_k;

          //record
          std::vector<double> res(12);
          res[0] = l0_A; res[1] = l0_B; res[2] = l0_C;
          res[3] = l1_A; res[4] = l1_B; res[5] = l1_C;
          res[6] = l2_A; res[7] = l2_B; res[8] = l2_C;
          res[9] = l3_A; res[10] = l3_B; res[11] = l3_C;
          resList.push_back(res);

          double l0_A_over = 1.0, l0_B_over = -1.0/paral_k, l0_C_over = (mid_paralB + length / cos(angleK_paral)) / paral_k;
          double l1_A_over = 1.0, l1_B_over = -1.0/ver_k, l1_C_over = (mid_verB + width / cos(angleK_ver)) / ver_k;
          double l2_A_over = 1.0, l2_B_over = l0_B_over, l2_C_over = (mid_paralB - length / cos(angleK_paral)) / paral_k;
          double l3_A_over = 1.0, l3_B_over = l1_B_over, l3_C_over = (mid_verB - width / cos(angleK_ver)) / ver_k;

          std::vector<double> res_over(12);
          res_over[0] = l0_A_over; res_over[1] = l0_B_over; res_over[2] = l0_C_over;
          res_over[3] = l1_A_over; res_over[4] = l1_B_over; res_over[5] = l1_C_over;
          res_over[6] = l2_A_over; res_over[7] = l2_B_over; res_over[8] = l2_C_over;
          res_over[9] = l3_A_over; res_over[10] = l3_B_over; res_over[11] = l3_C_over;
          resOverlappingList.push_back(res_over);

      }
      //for test
      std::ofstream out_c;
      out_c.open(sfm_data.s_root_path + "/OutC.res");
      Vec3 C0_out = sfm_data.poses[downIdList[0]].center();
      out_c << C0_out(0) << " " << C0_out(1) << std::endl;
      for(std::size_t i = 1; i < downIdList.size()+1; ++i)
      {
//              Vec3 c0 = sfm_data.poses[downIdList[i-1]].center();
          Vec3 c1 = sfm_data.poses[downIdList[i]].center();
          out_c << c1(0) << " " << c1(1) << std::endl;
      }
      out_c.close();

//          //for test
//          std::ofstream out_cW;
//          out_cW.open(sfm_data.s_root_path + "/OutCW.res");
//          Vec3 CW0 = sfm_data.poses[anotherLineList[0]].center();
//          out_cW << CW0(0) << " " << CW0(1) << std::endl;
//          for(std::size_t i = 1; i < anotherLineList.size()+1; ++i)
//          {
////              Vec3 c0 = sfm_data.poses[downIdList[i-1]].center();
//              Vec3 cW = sfm_data.poses[anotherLineList[i]].center();
//              out_cW << cW(0) << " " << cW(1) << std::endl;
//          }
//          out_cW.close();

      //output
      std::ofstream out_res;
      out_res.open(sfm_data.s_root_path + "/RegionPiece_line.res");
      if(!out_res)
      {
          std::cout << "create file RegionPiece_line.res failed, please check it!" << std::endl;
      }
      out_res << resList.size() << std::endl;
      int firstId = imgGroupList[0].id;
      for(std::size_t j = 0; j < resList.size(); ++j)
      {
          std::stringstream ssId;
          ssId << (j+firstId);
          std::string sId;
          ssId >> sId;
          out_res << sId << " ";
          for(std::size_t k = 0; k < 12; ++k)
          {
              std::stringstream ss;
              ss << resList[j][k];
              std::string s;
              ss >> s;
              out_res << s << " ";
          }
          out_res << std::endl;
      }
      out_res.close();

      std::ofstream out_resOver;
      out_resOver.open(sfm_data.s_root_path + "/RegionPiece_overlapping_line.res");
      if(!out_resOver)
      {
          std::cout << "create file RegionPiece_overlapping_line.res failed, please check it!" << std::endl;
      }
      out_resOver << resOverlappingList.size() << std::endl;
      for(std::size_t j = 0; j < resOverlappingList.size(); ++j)
      {
          std::stringstream ssId;
          ssId << j;
          std::string sId;
          ssId >> sId;
          out_resOver << sId << " ";
          for(std::size_t k = 0; k < 12; ++k)
          {
              std::stringstream ss;
              ss << resOverlappingList[j][k];
              std::string s;
              ss >> s;
              out_resOver << s << " ";
          }
          out_resOver << std::endl;

      }
      out_resOver.close();


  }

///////////////////////////////////////////////////////////////////////////////
  std::map<int, std::vector<std::pair<int, Landmark> > > recordLandmarkMap;
  if(needDivideGroup != "")
  {
      std::cout << "need divide in to group." << std::endl;

      if(needChangeName2Id)
      {

          std::ofstream out_group;
          std::ifstream in_groupName;

          in_groupName.open(sfm_data.s_root_path + "/groupName.res");
          if(!in_groupName)
          {
              std::cout << "open file " << sfm_data.s_root_path + "/groupName.res failed, please check it!" << std::endl;
              return EXIT_FAILURE;
          }

          out_group.open(needDivideGroup);
          if(!out_group)
          {
              std::cout << "create file " << sfm_data.s_root_path + "/group.res failed, please check it!" << std::endl;
              return EXIT_FAILURE;
          }

          int groupCount;
          in_groupName >> groupCount;
          out_group << groupCount << std::endl;
          int groupLine = 0;
          while(groupLine < groupCount && !in_groupName.eof())
          {
              int groupId;
              in_groupName >> groupId;
              out_group << groupId << std::endl;

              int image_loacl = 0;
              ImgGroup imgGroup;
              std::vector<int> allId;

              int row_count;
              in_groupName >> row_count;
              out_group << row_count << std::endl;
              while(image_loacl < row_count && !in_groupName.eof())
              {
                  int imgCount;
                  in_groupName >> imgCount;
                  std::stringstream ss;
                  ss << imgCount;
                  std::string s;
                  ss >> s;
                  out_group << s << " ";

                  int img_local = 0;

                  while(img_local < imgCount && !in_groupName.eof())
                  {

                      std::string thisImageName;
                      in_groupName >> thisImageName;

                      //find id
                      for(Views::iterator it_View = sfm_data.views.begin();
                          it_View != sfm_data.views.end(); ++ it_View)
                      {
                          std::cout << "it_View->s_Img_path : " << it_View->second->s_Img_path << std::endl;

//                          int f = it_View->second->s_Img_path.find_last_of('/');
//                          std::string sub_str;
//                          sub_str = it_View->second->s_Img_path.substr(f);
//                          std::cout << "sub_str : " << sub_str << std::endl;
                          if(it_View->second->s_Img_path == thisImageName)
                          {
                              std::stringstream ss;
                              ss << it_View->second->id_view;
                              std::string s;
                              ss >> s;
                              out_group << s <<" ";
                          }

                      }

                      ++img_local;
                  }

                  out_group << std::endl;
                  ++image_loacl;

              }

              ++groupLine;
          }

          out_group.close();
          in_groupName.close();

      }

      std::ifstream in_group;//group.res
      in_group.open(needDivideGroup);
      if(!in_group)
      {
          std::cout << "do not load file group.res, please check it !" << std::endl;
          return EXIT_FAILURE;
      }

      //
      //std::list<std::pair<int, std::vector<int> > > groupList; //pair::<downId, <lastId> >

      ImgGroupList imgGroupList;

      int groupCount;
      in_group >> groupCount;
      int groupLine = 0;
      while(groupLine < groupCount && !in_group.eof())
      {
          int image_loacl = 0;
          ImgGroup imgGroup;
          std::vector<int> allId;

          int groupId;
          in_group >> groupId;
          imgGroup.id = groupId;

          int row_count;
          in_group >> row_count;
          while(image_loacl < row_count && !in_group.eof())
          {
              int imgCount;
              in_group >> imgCount;

              int img_local = 1;

              int thisImageId;
              in_group >> thisImageId;
              allId.push_back(thisImageId);
              std::vector<int> refList;

              while(img_local < imgCount && !in_group.eof())
              {
                  int imgId;
                  in_group >> imgId;
                  refList.push_back(imgId);
                  allId.push_back(imgId);

                  ++img_local;
              }

              std::pair<int, std::vector<int> > tempPair;
              tempPair.first = thisImageId;
              tempPair.second = refList;

              imgGroup.groupList.push_back(tempPair);

              ++image_loacl;
          }
          imgGroup.allId = allId;

          imgGroupList.push_back(imgGroup);

          ++groupLine;
      }
      in_group.close();


      //prepare data
      //int posesCount = 0;
      std::list<bool> posesList(sfm_data.views.size(), true);//false);

//      #pragma omp parallel for
//      for(std::size_t k = 0; k < sfm_data.poses.size(); ++k)
////            it_pose != sfm_data.poses.end(); ++it_pose)
//      {
//          Poses::const_iterator it_pose = sfm_data.poses.begin();
//          std::advance(it_pose, k);
//          std::list<bool>::iterator it_pList = posesList.begin();
//          std::advance(it_pList, it_pose->first);
//          *it_pList = true;
//          #pragma omp critical
//          {
//              ++posesCount;
//          }

//      }

      //prepare folder
      std::string dir_path = sfm_data.s_root_path + "/Group";
      if(newFolder)
      {
          dir_path = sfm_data.s_root_path + "/GroupNew";
      }
      std::cout << "group folder : " << dir_path << std::endl;
      ::mkdir(dir_path.c_str(), S_IRWXU | S_IRGRP | S_IXGRP);

      std::ofstream out_dmapReference;
      if(newFolder)
      {
          out_dmapReference.open(sfm_data.s_root_path + "/dmapReferenceNew.res");
      }else{
          out_dmapReference.open(sfm_data.s_root_path + "/dmapReference.res");
      }

      if(!out_dmapReference)
      {
          std::cout << "create file dmapReference.res failed !" << std::endl;
          return EXIT_FAILURE;
      }
      groupCount = 0;
      for(std::size_t i = 0; i < imgGroupList.size(); ++i)
      {
          groupCount += imgGroupList[i].groupList.size();

      }
      out_dmapReference << groupCount << std::endl;



      //////////////////////////////////////////////////////






      //traverse structure
      //if(!needDivideForDensify)
      std::map<int, std::vector<std::pair<int, Landmark> > > recordLandmarkMap; //<imgId, <structureId, data> >
      {
          //prepare map data
          for(std::size_t i = 0; i < imgGroupList.size(); ++i)
          {
              for(std::size_t j = 0; j < imgGroupList[i].allId.size(); ++j)
              {
                  int imgId = imgGroupList[i].allId[j];
                  if(recordLandmarkMap.find(imgId) == recordLandmarkMap.end())
                  {
                      std::pair<int, std::vector<std::pair<int, Landmark> > > pair;
                      pair.first = imgId;

                      recordLandmarkMap.insert(pair);
                  }
              }
          }
        #pragma omp parallel for
        for (std::size_t structureId = 0; structureId < sfm_data.structure.size(); ++structureId)
        {
            Landmark it_Landmark = sfm_data.structure[structureId];
            for(Observations::iterator it_obs = it_Landmark.obs.begin();
                it_obs != it_Landmark.obs.end(); ++it_obs)
            {
                int imgId = it_obs->first;
                {
                    std::map<int, std::vector<std::pair<int, Landmark> > >::iterator it_record = recordLandmarkMap.find(imgId);

                    if(it_record == recordLandmarkMap.end())
                    {
                        continue;
                    }

                    #pragma omp critical
                    {
                        std::pair<int, Landmark> tempPair;
                        tempPair.first = structureId;
                        tempPair.second.obs[imgId] = it_Landmark.obs[imgId];
                        tempPair.second.X = it_Landmark.X;

                        it_record->second.push_back(tempPair);

                    }
                }
            }
        }

      }

      #pragma omp parallel for
      //cannot use parallel
      for(std::size_t i = 0; i < imgGroupList.size(); ++i)
      {
          ImgGroupList::iterator it_group = imgGroupList.begin();
          std::advance(it_group, i);

          int downId = it_group->groupList[0].first;

          std::list<bool>::iterator it_pos = posesList.begin();
          std::advance(it_pos, downId);
          if(*it_pos == false)
          {
              continue;
          }

          SfM_Data sub_sfm_data;
          sub_sfm_data.s_root_path = sfm_data.s_root_path;

          std::vector<int> dmapIds;//a group id in which need compute dmap
          for(std::size_t j = 0; j < it_group->groupList.size(); ++j)
          {
              dmapIds.push_back(it_group->groupList[j].first);
          }

          //add camera and traverse camera for instrinsic
          std::vector<bool> intrinsic_map(3,false);

          for(std::size_t j = 0; j < it_group->allId.size(); ++j)
          {
              std::vector<int>::iterator it_restId = it_group->allId.begin();
              std::advance(it_restId, j);

              int imgId = *it_restId;
              sub_sfm_data.s_root_path = sfm_data.s_root_path;
              sub_sfm_data.views[imgId] = sfm_data.views[imgId];
              sub_sfm_data.poses[imgId] = sfm_data.poses[imgId];
              const View * view = sfm_data.views.at(imgId).get();
              intrinsic_map[view->id_intrinsic] = true;

              //structure
              std::map<int, std::vector<std::pair<int, Landmark> > >::iterator it_record = recordLandmarkMap.find(imgId);
              for(std::vector<std::pair<int, Landmark> >::iterator it_map =  it_record->second.begin();
                  it_map != it_record->second.end(); ++it_map)
              {
                  sub_sfm_data.structure[it_map->first].obs[imgId] = it_map->second.obs[imgId];
                  sub_sfm_data.structure[it_map->first].X = it_map->second.X;
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
          std::vector<int> sortingList = it_group->allId;
          std::sort(sortingList.begin(),sortingList.end());

          sortingList.erase(std::unique(sortingList.begin(), sortingList.end()), sortingList.end());

          //get new id
          std::vector<int> newIdList(it_group->groupList.size()); //new id for dmap image (6)
          std::vector<std::vector<int> > newRefIdList(it_group->groupList.size());
          for(std::size_t k = 0; k < sortingList.size(); ++k)  //for down
          {
              for(std::size_t n =0; n < newIdList.size(); ++n)
              {
                  if(sortingList[k] == it_group->groupList[n].first)
                  {
                      newIdList[n] = k;
                  }
              }

              for(std::size_t n =0; n < newIdList.size(); ++n)
              {
                  for(std::vector<int>::iterator it = it_group->groupList[n].second.begin();
                      it != it_group->groupList[n].second.end(); ++it)
                  {
                      int thisId = *it;
                      if(sortingList[k] == thisId)
                      {
                          newRefIdList[n].push_back(k);
                      }
                  }

              }
          }



          //save dmapReference.res
          #pragma omp critical
          {

              //Export
              std::stringstream ss_down;
              ss_down << downId;
              std::string s_down;
              ss_down >> s_down;
              std::string fileName_down = "/group_" + s_down + "_";
              std::string filePath_down = dir_path + "/group_" + s_down + "_";
              std::cout << "...Export SfM_Data " << s_down <<" to disk." << std::endl;

              Save(sub_sfm_data,
                   stlplus::create_filespec(dir_path, fileName_down, ".bin"),
                   ESfM_Data(ALL));
              std::string fileName_ply = "cloud_and_poses_group_" + s_down;
              Save(sub_sfm_data,
                   stlplus::create_filespec(dir_path, fileName_ply, ".ply"),
                   ESfM_Data(ALL));

              for(std::size_t p = 0; p < it_group->groupList.size(); ++p)
              {
                  std::stringstream ss_old, ss_new;
                  ss_old << it_group->groupList[p].first;
//                  ss_old << it_group->second[dmapIds[p]];
                  ss_new << newIdList[p];
                  std::string s_old, s_new;
                  ss_old >> s_old;
                  ss_new >> s_new;
                  out_dmapReference << s_old << " " << s_new << " " << filePath_down << std::endl;
                  std::stringstream ss_count;
                  ss_count << newRefIdList[p].size();
                  std::string s_count;
                  ss_count >> s_count;
                  out_dmapReference << s_count << " ";
                  for(std::size_t k = 0; k < newRefIdList[p].size(); ++k)
                  {
                      std::stringstream ss;
                      ss << newRefIdList[p][k];
                      std::string s;
                      ss >> s;
                      out_dmapReference << s << " ";
                  }
                  out_dmapReference << std::endl;
              }
          }

      }

      std::string subFileName;
      int lastNum = sSfM_Data_Filename.find_last_of('.');
      subFileName = sSfM_Data_Filename.substr(0,lastNum);
      out_dmapReference << "0 0 " << subFileName << std::endl;
      out_dmapReference.close();


      /////////////////////////////////////////////////////
      //get file L0SceneName.txt
      //prepare file path
      std::string outFilePath = sfm_data.s_root_path + "/L0";//!!!need par
      ::mkdir(outFilePath.c_str(), S_IRWXU | S_IRGRP | S_IXGRP);
      std::string fileName = outFilePath + "/L0SceneName.txt";
      std::cout << "get L0SceneName.txt file" << std::endl;
      std::ofstream out_recordPath;
      out_recordPath.open(fileName);
      if(!out_recordPath)
      {
          std::cout << "create out_recordPath " << fileName << " failed !";
          return EXIT_FAILURE;
          //continue;
      }
      out_recordPath << "0 " << outFilePath << std::endl;
      out_recordPath << imgGroupList.size() << std::endl;

      imgGroupList;
      std::stringstream ssBid;
      ssBid << bId;
      std::string sBid;
      ssBid >> sBid;
      for(std::size_t i = 0; i < imgGroupList.size(); ++i)
      {
          int id = imgGroupList[i].id;
          std::stringstream ssRid;
          ssRid << id;
          std::string sRid;
          ssRid >> sRid;
          const std::string finalScenePath = outFilePath + "/R" + sRid + "/Scene_L0_B" + sBid + "_R" + sRid;

          out_recordPath << sRid << " " << finalScenePath << std::endl;

      }
      out_recordPath.close();


//      std::pair<int, std::string> regionPathPair;
//      regionPathPair.first = regionId;
//      regionPathPair.second = finalScenePath;
//      record_regionPath.insert(regionPathPair);



//      for(std::map<int, std::string>::iterator it = record_regionPath.begin();
//          it != record_regionPath.end(); ++it)
//      {
//          std::stringstream ss;
//          ss << it->first;
//          std::string s;
//          ss >> s;
//          out_recordPath << s << " " << it->second << std::endl;
//      }

//          out_recordPath.close();

//          std::cout << "!!!!!!!!!!!!!! finish" << std::endl;


  }

///////////////////////////////////////////////////////////////////////////////
  if(needDivideForDensify != "")
  {
      std::cout << "need divide for densify." << std::endl;

      Save(sfm_data,
        stlplus::create_filespec(sfm_data.s_root_path, "test_cloud_and_poses", ".ply"),
        ESfM_Data(ALL));

      std::string dir_path = sfm_data.s_root_path + "/EachCam";
      if(newFolder)
      {
          dir_path = sfm_data.s_root_path + "/EachCamNew";
      }
      std::cout << "EachCam folder : " << dir_path << std::endl;
      ::mkdir(dir_path.c_str(), S_IRWXU | S_IRGRP | S_IXGRP);

      int posesCount = 0;
//      std::vector<bool> poses_map(sfm_data.views.size(),false);
      for(Poses::const_iterator it_pose = sfm_data.poses.begin();
          it_pose != sfm_data.poses.end(); ++it_pose)
      {
//          poses_map[it_pose->first] = true;
          ++posesCount;
      }
      ///-------------------------------------------------///

      std::ifstream in_groupDmap;
      in_groupDmap.open(needDivideForDensify);
      if(!in_groupDmap)
      {
          std::cout << "do not load file dmap_x.res, please check it !" << std::endl;
          return EXIT_FAILURE;
      }

      std::vector<std::vector<int> > groupImageIdList;

      int dmapCount;
      in_groupDmap >> dmapCount;
      int groupLine = 0;
      while(groupLine < dmapCount && !in_groupDmap.eof())
      {
          int imageCount;
          in_groupDmap >> imageCount;

          std::vector<int> imgIdList;
          int row = 0;
          while(row < imageCount && !in_groupDmap.eof())
          {
              int imgId;
              in_groupDmap >> imgId;
              imgIdList.push_back(imgId);

              ++row;
          }

          groupImageIdList.push_back(imgIdList);

          ++groupLine;
      }
      in_groupDmap.close();
      ///-------------------------------------------------///

//      for (Poses::const_iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose)
//      {
//          const Pose3 & pose = itPose->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 t = pose.translation();
//          const Vec3 C = pose.center();
//          std::cout << R << std::endl << std::endl;

//      }

      //output file EachCameraInfo.txt
      //for trans and depth-map
      std::ofstream out_EachCameraInfo;
      if(newFolder)
      {
          out_EachCameraInfo.open(sfm_data.s_root_path + "/EachCameraInfoNew.txt");
      }else{
          out_EachCameraInfo.open(sfm_data.s_root_path + "/EachCameraInfo.txt");
      }
      std::stringstream ss_dmapSize;
      ss_dmapSize << groupImageIdList.size();
      std::string s_dmapSize;
      ss_dmapSize >> s_dmapSize;
      out_EachCameraInfo << s_dmapSize << std::endl;

      std::map<int, std::vector<std::pair<int, Landmark> > > recordLandmarkMap; //<imgId, <structureId, data> >
      //prepare map data
      for(std::size_t i = 0; i < groupImageIdList.size(); ++i)
      {
          for(std::size_t j = 0; j < groupImageIdList[i].size(); ++j)
          {
              int imgId = groupImageIdList[i][j];
              if(recordLandmarkMap.find(imgId) == recordLandmarkMap.end())
              {
                  std::pair<int, std::vector<std::pair<int, Landmark> > > pair;
                  pair.first = imgId;

                  recordLandmarkMap.insert(pair);
              }
          }
      }

      //traverse structure
      #pragma omp parallel for
      for (std::size_t structureId = 0; structureId < sfm_data.structure.size(); ++structureId)
      {
          Landmark it_Landmark = sfm_data.structure[structureId];
          for(Observations::iterator it_obs = it_Landmark.obs.begin();
              it_obs != it_Landmark.obs.end(); ++it_obs)
          {
              int imgId = it_obs->first;
              {
                  std::map<int, std::vector<std::pair<int, Landmark> > >::iterator it_record = recordLandmarkMap.find(imgId);

                  if(it_record == recordLandmarkMap.end())
                  {
                      continue;
                  }

                  #pragma omp critical
                  {
                      std::pair<int, Landmark> tempPair;
                      tempPair.first = structureId;
                      tempPair.second.obs[imgId] = it_Landmark.obs[imgId];
                      tempPair.second.X = it_Landmark.X;

                      it_record->second.push_back(tempPair);

                  }
              }
          }
      }
      std::cout << "traverse structure finished." << std::endl;

      //for every dmap group
      #pragma omp parallel for
      for(std::size_t i = 0; i < groupImageIdList.size(); ++i)
      {
          SfM_Data sub_sfmData;
          std::vector<bool> intrinsic_map(3,false);
          for(std::size_t j = 0; j < groupImageIdList[i].size(); ++j)
          {
              int imgId = groupImageIdList[i][j];

              sub_sfmData.s_root_path = sfm_data.s_root_path;
              sub_sfmData.views[imgId] = sfm_data.views[imgId];
              sub_sfmData.poses[imgId] = sfm_data.poses[imgId];
              const View * view = sfm_data.views.at(imgId).get();
              intrinsic_map[view->id_intrinsic] = true;

              //structure
              std::map<int, std::vector<std::pair<int, Landmark> > >::iterator it_record = recordLandmarkMap.find(imgId);
              for(std::vector<std::pair<int, Landmark> >::iterator it_map =  it_record->second.begin();
                  it_map != it_record->second.end(); ++it_map)
              {
                  sub_sfmData.structure[it_map->first].obs[imgId] = it_map->second.obs[imgId];
                  sub_sfmData.structure[it_map->first].X = it_map->second.X;
              }
          }

          //add intrinsic
          for(std::size_t k = 0; k < intrinsic_map.size(); ++k)
          {
              if(intrinsic_map[k] == true)
              {
                  sub_sfmData.intrinsics[k] = sfm_data.intrinsics[k];
              }
          }

          //sorting
          std::vector<int> sortingList = groupImageIdList[i];
          std::sort(sortingList.begin(),sortingList.end());

          sortingList.erase(std::unique(sortingList.begin(), sortingList.end()), sortingList.end());

          //get new id
          std::vector<int> newIdList(groupImageIdList[i].size()); //new id for dmap image (6)
          for(std::size_t k = 0; k < sortingList.size(); ++k)  //for down
          {
              for(std::size_t n =0; n < newIdList.size(); ++n)
              {
                  if(sortingList[k] == groupImageIdList[i][n])
                  {
                      newIdList[n] = k;
                  }
              }
          }

          //Export
          std::stringstream ss_oldId, ss_newId;
          ss_oldId << groupImageIdList[i][0];
          ss_newId << newIdList[0];
          std::string s_oldId, s_newId;
          ss_oldId >> s_oldId;
          ss_newId >> s_newId;





          #pragma omp critical
          {
              std::cout << "...Export Dmap_Data " << s_oldId <<" to disk." << std::endl;

              std::string fileName_bin = "/depth_map_" + s_oldId + "_";
              Save(sub_sfmData,
                   stlplus::create_filespec(dir_path, fileName_bin, ".bin"),
                   ESfM_Data(ALL));

              std::string fileName_ply = "cloud_and_poses_sparse_" + s_oldId;
              Save(sub_sfmData,
                   stlplus::create_filespec(dir_path, fileName_ply, ".ply"),
                   ESfM_Data(ALL));


              out_EachCameraInfo << s_oldId << " " << s_newId << " " << dir_path + fileName_bin << std::endl;

              std::stringstream ss_size;
              ss_size << groupImageIdList[i].size() - 1;
              std::string s_size;
              ss_size >> s_size;

              out_EachCameraInfo << s_size << " ";
              for(std::size_t k = 1; k < newIdList.size(); ++k)
              {
                  std::stringstream ss_newId;
                  ss_newId << newIdList[k];
                  std::string s_newId;
                  ss_newId >> s_newId;

                  out_EachCameraInfo << s_newId << " ";

              }
              out_EachCameraInfo << std::endl;
          }
      }

      if(changeSubBlock)
      {
          std::stringstream ss;
          ss << bId;
          std::string sbId;
          ss >> sbId;
          std::string subFileName;
          subFileName = changeRootPath + "/matches_"+sbId+"/sfm_data_"+ sbId;
          out_EachCameraInfo << "0 0 " << subFileName << std::endl;

      }else{
          std::string subFileName;
          int lastNum = sSfM_Data_Filename.find_last_of('.');
          subFileName = sSfM_Data_Filename.substr(0,lastNum);
          out_EachCameraInfo << "0 0 " << subFileName << std::endl;
      }
      out_EachCameraInfo.close();


  }

///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////

  return EXIT_SUCCESS;
}
