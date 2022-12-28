// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 cDc <cdc.seacave@gmail.com>, Pierre MOULON

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/cameras/Camera_undistort_image.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"

#define _USE_EIGEN
#include "InterfaceMVS.h"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress_display.hpp"

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/cameras/Cameras_Common_command_line_helper.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"
#include "openMVG/sfm/pipelines/global/sfm_global_engine_relative_motions.hpp"

#include "openMVG/sfm/sfm_report.hpp"
#include "openMVG/system/timer.hpp"



using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::image;
using namespace openMVG::sfm;

#include <cstdlib>
#include <string>
#include <sys/stat.h>
#include <list>

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
  bool needDivideForDensify = false;
  bool needDivideGroup = false;
  bool needRt = false;
  bool needChangeName2Id = false;
  bool needChangeAxie = false;
//  bool demoSituation = false;
  bool needAutoDivideMesh = false;
  int timesAirway = 1;

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
  cmd.add( make_option('t', timesAirway, "times of airway step in divide mesh process") );
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

//  std::string sOut = "/home/guang/data/20170831/m2/";
//  Save(sfm_data,
//    stlplus::create_filespec(sOut, "sfm_data", ".bin"),
//    ESfM_Data(ALL));

//  Save(sfm_data,
//    stlplus::create_filespec(sOut, "cloud_and_poses", ".ply"),
//    ESfM_Data(ALL));

  bool needPoses = false;//strue;
  bool needTracks = true;

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

///////////////////////////////////////////////////////////////////////////////
  if(needPoses)
  {
      std::ofstream out_show_camera;
      out_show_camera.open(sfm_data.s_root_path + "/out_poses.txt");
      for(Poses::const_iterator it = sfm_data.poses.begin();
          it != sfm_data.poses.end(); ++it)
      {
          const Pose3 & pose = it->second;
          const Mat3 R = pose.rotation();
          const Vec3 c = pose.center();

          out_show_camera <<R(0)<<" "<< R(1)<< " "<<R(2)<< " "
                         <<R(3)<<" "<< R(4)<< " "<<R(5)<< " "
                         <<R(6)<<" "<< R(7)<< " "<<R(8)<< " "
                         <<c(0)<<" "<< c(1)<< " "<<c(2) << std::endl;

      }
      out_show_camera.close();
  }

///////////////////////////////////////////////////////////////////////////////
  if(needTracks)
  {
      std::string sMatchesDir = sfm_data.s_root_path+"/matches_grid";
      const std::string sImage_describer = stlplus::create_filespec(sMatchesDir, "image_describer", "json");
      std::unique_ptr<openMVG::features::Regions> regions_type = openMVG::features::Init_region_type_from_file(sImage_describer);
      if (!regions_type)
      {
        std::cerr << "Invalid: "
          << sImage_describer << " regions type file." << std::endl;
        return EXIT_FAILURE;
      }
      // Features reading
      std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
      feats_provider->load(sfm_data, sMatchesDir, regions_type);

      std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
      matches_provider->load(sfm_data, stlplus::create_filespec(sMatchesDir, "matches.e.bin"));

      const matching::PairWiseMatches& pairMatch = matches_provider->pairWise_matches_;
      Pair_Set pairSet = matches_provider->getPairs();

      using IndMatches = std::vector<matching::IndMatch>;

      for(std::size_t i = 0; i < pairMatch.size(); ++i)
      {
          std::map< Pair, IndMatches >::iterator itMP = matches_provider->pairWise_matches_.begin();
          std::advance(itMP, i);
          IndMatches& idM = itMP->second;//->second;

          int imgI = itMP->first.first;
          int imgJ = itMP->first.second;
          std::stringstream ssImgI, ssImgJ;
          ssImgI << imgI;ssImgJ << imgJ;
          std::string sImgI, sImgJ;
          ssImgI >> sImgI; ssImgJ >> sImgJ;

          std::ofstream out_feature;
          out_feature.open(sfm_data.s_root_path + "/tracks/out_tracks_"+sImgI+"_"+sImgJ+".txt");
          out_feature << pairMatch.size() << std::endl;

          for(std::size_t j = 0; j < idM.size(); ++j)
          {
              int i_ = idM[j].i_;
              int j_ = idM[j].j_;

              float imgIx = feats_provider->feats_per_view[imgI][i_].x();
              float imgIy = feats_provider->feats_per_view[imgI][i_].y();
              float imgJx = feats_provider->feats_per_view[imgJ][j_].x();
              float imgJy = feats_provider->feats_per_view[imgJ][j_].y();

              out_feature << imgIx << " " << imgIy << " " << imgJx << " " << imgJy << std::endl;
              //std::cout << imgIx << " " << imgIy << " " << imgJx << " " << imgJy << std::endl;

          }

          out_feature.close();


      }


//      pairs = getPairs(matches);





  }

///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////

  return EXIT_SUCCESS;
}
