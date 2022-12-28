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
  std::string checkPointsList_Filename;
////  std::string sOutFile = "scene.mvs";
//  std::string sOutDir = "undistorted_images";
//  std::string sInputFileName = "";
//  bool undistortImage = false;
//  std::string needDivideForDensify = ""; //damp_x.res
//  std::string needDivideGroup = ""; //group_x.res
//  bool needRt = false;
//  bool needChangeName2Id = false;
//  bool needChangeAxie = false;
////  bool demoSituation = false;
//  bool needAutoDivideMesh = false;
//  int timesAirway = 1;

//  int bId = -1;
  cmd.add( make_option('i', sSfM_Data_Filename, "sfm_data.bin") );
  cmd.add( make_option('l', checkPointsList_Filename, "checkPointsList") );


//  cmd.add( make_option('d', sOutDir, "outdir") );
//  cmd.add( make_option('f', sInputFileName, "input file name") );
//  cmd.add( make_option('r', needRt, "need out put rt file") );
//  cmd.add( make_option('u', undistortImage, "undistort image") );
//  cmd.add( make_option('n', needDivideForDensify, "need do dividing for mutli-thread depth map process") );
//  cmd.add( make_option('g', needDivideGroup, "need do dividing for mutli-thread depth map process according groups") );
//  cmd.add( make_option('c', needChangeName2Id, "need prepare for dividing for mutli-thread depth map process according groups to change image file name to image id") );
//  cmd.add( make_option('a', needAutoDivideMesh, "need get divide mesh file automaticlly") );
//  cmd.add( make_option('x', needChangeAxie, "need change aeix-z") );
//  cmd.add( make_option('b', bId, "block id, same as root_path/block_path") );
//  cmd.add( make_option('t', timesAirway, "times of airway step in divide mesh process") );
////  cmd.add( make_option('s', demoSituation, "need prepare for dividing for mutli-thread depth map process according groups to change image file name to image id, just in demo situation") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch (const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--sfmdata] filename, the SfM_Data file to convert\n"
      << "[-l|--checkPointsList] a file of check points list\n"
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

  std::ifstream in;
  in.open(checkPointsList_Filename);
  std::map<int, Vec3> initPMap;//<GCP_Id, location>
  if(in)
  {
      int line = 0;
      in >> line;
      int l = 0;
      while(!in.eof() && l < line)
      {
          Vec3 tempP;
          int id;
          in >> id >> tempP(0) >> tempP(1) >> tempP(2);
          std::pair<int, Vec3> tempPair;
          tempPair.first = id;
          tempPair.second = tempP;
          initPMap.insert(tempPair);

          ++l;
      }

  }else{
      std::cout << "open file " << checkPointsList_Filename << " failed !" << std::endl;
      return EXIT_FAILURE;
  }
  in.close();


  std::map<int, std::vector<std::pair<int, Vec2> > > resMap;//<GCP id, <imgId, P_xy> >
  for(std::size_t t = 0; t < initPMap.size(); ++t)
  {
      std::map<int, Vec3>::iterator itP = initPMap.begin();
      std::advance(itP, t);
      Vec3 thisP = itP->second;
      std::vector<std::pair<int, Vec2> > thisResList;
      for(auto itPose : sfm_data.poses)
      {
          const Mat3 cam_R = itPose.second.rotation();
          const Vec3 cam_t = itPose.second.translation();

          Vec3 pos_proj;
          pos_proj = cam_R * thisP + cam_t;

          const double x_u = pos_proj(0) / pos_proj(2);
          const double y_u = pos_proj(1) / pos_proj(2);

          std::vector<double> cam_intrinsics = sfm_data.intrinsics[sfm_data.views[itPose.first]->id_intrinsic]->getParams();
          const double focal = cam_intrinsics[0];
          const double principal_point_x = cam_intrinsics[1];
          const double principal_point_y = cam_intrinsics[2];

          const double projected_x = principal_point_x + focal * x_u;
          const double projected_y = principal_point_y + focal * y_u;

          const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(sfm_data.views[itPose.first]->id_intrinsic).get();
          Vec2 projected_d = cam->get_d_pixel( Vec2(projected_x, projected_y) );

          int h = sfm_data.intrinsics[sfm_data.views[itPose.first]->id_intrinsic].get()->h();
          int w = sfm_data.intrinsics[sfm_data.views[itPose.first]->id_intrinsic].get()->w();

          if(projected_d(0) > (double)w|| projected_d(0) < 0.0 || projected_d(1) > (double)h || projected_d(1) < 0.0)
          {
              continue;
          }else{
              std::pair<int, Vec2> tempPxy;
              tempPxy.first = itPose.first;
              tempPxy.second = Vec2(projected_d(0), projected_d(1));
              thisResList.push_back(tempPxy);
          }

      }

      std::pair<int, std::vector<std::pair<int, Vec2> > > resPair;
      resPair.first = itP->first;
      resPair.second = thisResList;

      resMap.insert(resPair);

  }

  std::string sOutFile = sfm_data.s_root_path+"/checkPointsPreparing.txt";
  std::cout << sOutFile << std::endl;
  std::ofstream out;
  out.open(sOutFile);
  if(out)
  {
      for(std::size_t pId = 0; pId < resMap.size(); ++pId)
      {
          std::map<int, std::vector<std::pair<int, Vec2> > >::iterator itThis = resMap.begin();
          std::advance(itThis, pId);
          for(std::size_t pairId = 0; pairId < itThis->second.size(); ++pairId)
          {
              out << itThis->first << " "
                  << itThis->second[pairId].first << " "
                  << itThis->second[pairId].second(0) << " " << itThis->second[pairId].second(1) << std::endl;
              std::cout << itThis->first << " "
                        << itThis->second[pairId].first << " "
                        << itThis->second[pairId].second(0) << " " << itThis->second[pairId].second(1) << std::endl;
          }
      }

  }else{
      std::cout << "create file " << sOutFile << " failed !" << std::endl;
      return EXIT_FAILURE;
  }
  out.close();

  return EXIT_SUCCESS;
}
