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


int main(int argc, char *argv[])
{
  CmdLine cmd;
  std::string sSfM_Data_Filename;
  std::string sOutFile = "scene.mvs";
  std::string sOutDir = "undistorted_images";
  std::string sInputFileName_eachCamInfo = "";
  std::string sInputFileName_group = "";
  bool undistortImage = true;
  bool needTransWholeScene = true;
  bool forSparseScene = false;

  cmd.add( make_option('i', sSfM_Data_Filename, "sfmdata") );
  cmd.add( make_option('o', sOutDir, "outfile") );
//  cmd.add( make_option('d', sOutDir, "outdir") );
  cmd.add( make_option('f', sInputFileName_eachCamInfo, "input file name, file EachCameraInfo.txt" ) );
  cmd.add( make_option('g', sInputFileName_group, "input file name, file dmapReference.res" ) );
  cmd.add( make_option('u', undistortImage, "undistort image") );
  cmd.add( make_option('w', needTransWholeScene, "need trans whole Scene from sfm_data.bin to Scene.mvs") );
  cmd.add( make_option('t', forSparseScene, "for doing sparse scene mesh and texture.") );

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

  {
      std::cout << "out dir : " << sOutDir << std::endl;

      // Read the input SfM scene
      SfM_Data sfm_data;
      if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
        std::cerr << std::endl
          << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
        return EXIT_FAILURE;
      }

      //
      double d_l = 1.0, d_s = 1.0, d_f = 100.0;
      std::ofstream out_cExtrinsics;
      std::string sOutCEFile = sOutDir + "/CameraExtrinsics.txt";
      out_cExtrinsics.open(sOutCEFile.c_str());

      Mat3 K;
      K << 3047, 0, 2000,
              0, 3047, 1500,
              0, 0, 1;
      double c_scale = 0.05;
      for(const auto pose : sfm_data.poses)
      {
          const Vec3 C = pose.second.center();
          const Mat3 R = pose.second.rotation();
          const Vec3 t = -R*C;

//          Vec3 point0 = Vec3(-d_l/2.0, d_s/2.0, d_f);
//          Vec3 point1 = Vec3(d_l/2.0, d_s/2.0, d_f);
//          Vec3 point2 = Vec3(d_l/2.0, -d_s/2.0, d_f);
//          Vec3 point3 = Vec3(-d_l/2.0, -d_s/2.0, d_f);

          Vec3 point0 = Vec3(0,0,1);
          Vec3 point1 = Vec3(4000,0,1);
          Vec3 point2 = Vec3(4000,3000,1);
          Vec3 point3 = Vec3(0,3000,1);

//          Vec3 newP0 = R*point0 + t;
//          Vec3 newP1 = R*point1 + t;
//          Vec3 newP2 = R*point2 + t;
//          Vec3 newP3 = R*point3 + t;

          Vec3 newP0 = R.transpose()*(K.inverse()*point0*c_scale+R*C);
          Vec3 newP1 = R.transpose()*(K.inverse()*point1*c_scale+R*C);
          Vec3 newP2 = R.transpose()*(K.inverse()*point2*c_scale+R*C);
          Vec3 newP3 = R.transpose()*(K.inverse()*point3*c_scale+R*C);

          out_cExtrinsics << C(0) << " " << C(1) << " " << C(2) << std::endl << newP0(0) << " " << newP0(1) << " " << newP0(2) << std::endl << std::endl;
          out_cExtrinsics << C(0) << " " << C(1) << " " << C(2) << std::endl << newP1(0) << " " << newP1(1) << " " << newP1(2) << std::endl << std::endl;
          out_cExtrinsics << C(0) << " " << C(1) << " " << C(2) << std::endl << newP2(0) << " " << newP2(1) << " " << newP2(2) << std::endl << std::endl;
          out_cExtrinsics << C(0) << " " << C(1) << " " << C(2) << std::endl << newP3(0) << " " << newP3(1) << " " << newP3(2) << std::endl << std::endl;

          out_cExtrinsics << newP0(0) << " " << newP0(1) << " " << newP0(2) << std::endl << newP1(0) << " " << newP1(1) << " " << newP1(2) << std::endl << std::endl;
          out_cExtrinsics << newP1(0) << " " << newP1(1) << " " << newP1(2) << std::endl << newP2(0) << " " << newP2(1) << " " << newP2(2) << std::endl << std::endl;
          out_cExtrinsics << newP2(0) << " " << newP2(1) << " " << newP2(2) << std::endl << newP3(0) << " " << newP3(1) << " " << newP3(2) << std::endl << std::endl;
          out_cExtrinsics << newP3(0) << " " << newP3(1) << " " << newP3(2) << std::endl << newP0(0) << " " << newP0(1) << " " << newP0(2) << std::endl << std::endl;

      }
      out_cExtrinsics.close();

      std::ofstream out_points;

      std::string sOutPFile = sOutDir + "/RemoveWaterPoints.txt";
      out_points.open(sOutPFile.c_str());
      for(const auto thisP : sfm_data.structure)
      {
          const Vec3 X = thisP.second.X;
          out_points << X(0) << " " << X(1) << " " << X(2) << std::endl;
      }
      out_points.close();



  }


  return EXIT_SUCCESS;
}
