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
  std::string sSfM_Data_Path;
  int blockCount = 0;


////  std::string sOutFile = "scene.mvs";
//  std::string sOutDir = "undistorted_images";
//  std::string sInputFileName = "";
//  bool undistortImage = false;
//  bool needDivideForDensify = false;
//  bool needDivideGroup = false;
//  bool needRt = false;
//  bool needChangeName2Id = false;
//  bool needChangeAxie = false;
////  bool demoSituation = false;
//  bool needAutoDivideMesh = false;
//  int timesAirway = 1;

  cmd.add( make_option('i', sSfM_Data_Path, "a path with all small blocks") );
//  cmd.add( make_option('o', sOutFile, "outfile") );
  cmd.add( make_option('c', blockCount, "count of small blocks") );


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


  //out path
  const std::string dir_path = sSfM_Data_Path + "/whole";
  ::mkdir(dir_path.c_str(), S_IRWXU | S_IRGRP | S_IXGRP);
  const std::string out_path = sSfM_Data_Path + "/whole/matches";
  ::mkdir(out_path.c_str(), S_IRWXU | S_IRGRP | S_IXGRP);

  SfM_Data finalData;
  finalData.s_root_path = out_path;
  for(std::size_t idPiece = 0; idPiece < blockCount; ++idPiece)
  {

      std::stringstream ss;
      ss << idPiece;
      std::string s;
      ss >> s;
      std::string thisPath = sSfM_Data_Path + "/" + s + "/matches/sfm_data.bin";

      std::cout << "load : " << thisPath << std::endl;

      SfM_Data thisData;
      if (!Load(thisData, thisPath, ESfM_Data(ALL))) {
        std::cerr << std::endl
          << "The input SfM_Data file \""<< thisPath << "\" cannot be read." << std::endl;
        return EXIT_FAILURE;
      }

      for(Views::const_iterator itViews = thisData.views.begin();
          itViews != thisData.views.end(); ++itViews)
      {
          finalData.views[itViews->first] = itViews->second;
          finalData.views[itViews->first]->s_Img_path = thisPath + "/" + itViews->second->s_Img_path;
      }
      for(Poses::const_iterator itPoses = thisData.poses.begin();
          itPoses != thisData.poses.end(); ++itPoses)
      {
          finalData.poses[itPoses->first] = itPoses->second;
      }
      for(Intrinsics::const_iterator itIntrinsics = thisData.intrinsics.begin();
          itIntrinsics != thisData.intrinsics.end(); ++itIntrinsics)
      {
          finalData.intrinsics[itIntrinsics->first] = itIntrinsics->second;
      }
      const int oldSize = finalData.structure.size();
      #pragma omp parallel for
      for(std::size_t i=0; i< thisData.structure.size(); ++i)
              //(Landmarks::const_iterator itLandmarks = thisData.landmarks.begin();
          //itLandmarks != thisData.landmarks.end(); ++itLandmarks)
      {
          Landmarks::const_iterator itLandmarks = thisData.structure.begin();
          std::advance(itLandmarks, i);
          int newId = itLandmarks->first + oldSize;
          #pragma omp critical
          {
              finalData.structure[newId] = itLandmarks->second;
          }

      }


  }


  std::cout << "...Export SfM_Data to disk." << std::endl;
  Save(finalData,
    stlplus::create_filespec(out_path, "sfm_data", ".bin"),
    ESfM_Data(ALL));

  Save(finalData,
    stlplus::create_filespec(out_path, "cloud_and_poses", ".ply"),
    ESfM_Data(ALL));

//  cmd.add( make_option('f', sInputFileName, "input file name") );
//  cmd.add( make_option('r', needRt, "need out put rt file") );
//  cmd.add( make_option('u', undistortImage, "undistort image") );
//  cmd.add( make_option('n', num, "need do dividing for mutli-thread depth map process") );
//  cmd.add( make_option('g', needDivideGroup, "need do dividing for mutli-thread depth map process according groups") );
//  cmd.add( make_option('c', needChangeName2Id, "need prepare for dividing for mutli-thread depth map process according groups to change image file name to image id") );
//  cmd.add( make_option('a', needAutoDivideMesh, "need get divide mesh file automaticlly") );
//  cmd.add( make_option('x', needChangeAxie, "need change aeix-z") );
//  cmd.add( make_option('t', timesAirway, "times of airway step in divide mesh process") );
//  cmd.add( make_option('s', demoSituation, "need prepare for dividing for mutli-thread depth map process according groups to change image file name to image id, just in demo situation") );





///////////////////////////////////////////////////////////////////////////////

  return EXIT_SUCCESS;
}
