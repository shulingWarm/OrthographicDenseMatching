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

//#include <opencv2/core/core.hpp>
//#include <opencv2/highgui/highgui.hpp>
//#include <opencv2/imgproc/imgproc.hpp>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::image;
using namespace openMVG::sfm;

#include <cstdlib>
#include <string>

bool exportToOpenMVS(
  const SfM_Data & sfm_data,
  const std::string & sOutFile,
  const std::string & sOutDir,
  const bool undistortImage,
  const bool single_mvs
  )
{
  // Create undistorted images directory structure
  if (!stlplus::is_folder(sOutDir))
  {
    stlplus::folder_create(sOutDir);
    if (!stlplus::is_folder(sOutDir))
    {
      std::cerr << "Cannot access to one of the desired output directory" << std::endl;
      return false;
    }
  }

  // Export data :
  MVS::Interface scene;
  size_t nPoses(0);
  const uint32_t nViews((uint32_t)sfm_data.GetViews().size());

  C_Progress_display my_progress_bar(nViews);

  // OpenMVG can have not contiguous index, use a map to create the required OpenMVS contiguous ID index
  std::map<openMVG::IndexT, uint32_t> map_intrinsic, map_view;

  // define a platform with all the intrinsic group
  for (const auto& intrinsic: sfm_data.GetIntrinsics())
  {
    if (isPinhole(intrinsic.second->getType()))
    {
      const Pinhole_Intrinsic * cam = dynamic_cast<const Pinhole_Intrinsic*>(intrinsic.second.get());
      if (map_intrinsic.count(intrinsic.first) == 0)
        map_intrinsic.insert(std::make_pair(intrinsic.first, scene.platforms.size()));
      MVS::Interface::Platform platform;
      // add the camera
      MVS::Interface::Platform::Camera camera;
      camera.K = cam->K();
      // sub-pose
      camera.R = Mat3::Identity();
      camera.C = Vec3::Zero();
      platform.cameras.push_back(camera);
      scene.platforms.push_back(platform);
    }
  }

  /// prepare undistort step's remap file
  ///
  std::map<int, std::pair<std::vector<std::vector<double> >, std::vector<std::vector<double> > > > map_d_list;
  if(undistortImage)
  {

      const uint32_t images_w = sfm_data.views.begin()->second->ui_width;
      const uint32_t images_h = sfm_data.views.begin()->second->ui_height;
      std::cout << "images_w : " << images_w << " images_h : " << images_h << std::endl;

      for(Intrinsics::const_iterator itI = sfm_data.intrinsics.begin();
          itI != sfm_data.intrinsics.end(); ++itI)
      {
          std::vector<std::vector<double> > tempMapD_x(images_w, std::vector<double>(images_h));
          std::vector<std::vector<double> > tempMapD_y(images_w, std::vector<double>(images_h));

          #pragma omp parallel for
          for ( uint32_t j = 0; j < images_h; ++j )
            for ( uint32_t i = 0; i < images_w; ++i )
            {
                const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(itI->first).get();
                Vec2 d = cam->get_d_pixel( Vec2(i,j) );
                tempMapD_x[i][j] = d(0);
                tempMapD_y[i][j] = d(1);
    //            #pragma omp critical
    //            std::cout << "i : " << i << " j : " << j << " d0 : " << d(0) << " d1 : " << d(1)<< std::endl;
            }
          std::pair<int, std::pair<std::vector<std::vector<double> >, std::vector<std::vector<double> > > > tempPair;
          tempPair.first = itI->first;
          tempPair.second = std::pair<std::vector<std::vector<double> >, std::vector<std::vector<double> > >(tempMapD_x, tempMapD_y);
          //
          map_d_list.insert(tempPair);

//          std::cout << "finish 1" << std::endl;

      }
  }


  //  cv::Mat map_d[3], map_u[3];
//  map_d[2].create(cv::Size(sfm_data.GetViews().begin()->second.ui_width, sfm_data.GetViews().begin()->second.ui_width), CV_32FC1);
//  map_u[2].create(cv::Size(sfm_data.GetViews().begin()->second.ui_width, sfm_data.GetViews().begin()->second.ui_width), CV_32FC1);
//  {
//      for()
//      {
//          map_x[1].at<float>(y, x) = sx;
//          map_y[1].at<float>(y, x) = sy;
//      }
//      cv::remap(srcImage,dstImage,map_x,map_y,CV_INTER_LINEAR,BORDER_CONSTANT,Scalar(255,0,0));
//      std::string configfileD = undistortFilesPath+ "/datd.conf";//stlplus::create_filespec(undistortFilesPath, "datd.conf");
//      std::ifstream conf;
//      conf.open(configfileD.c_str());//, std::ios::binary);
//      if(!conf)
//      {
//          std::cerr << "open file " << configfileD << " failed!" << std::endl;
//          return EXIT_FAILURE;
//      }
//      std::string lx_d;
//      while(getline(conf, lx_d))
//      {
//          std::stringstream line(lx_d);
//          int x, y;
//          double sx, sy;
//          line >> x >> y >> sx >> sy;
//          map_x[1].at<float>(y, x) = sx;
//          map_y[1].at<float>(y, x) = sy;
//      }
//      conf.close();
//      std::cout << "load undistort file datd.conf finished" << std::endl;
//  }

//  for()
//  {

//  }
  // define images & poses
  scene.images.reserve(nViews);
  for (const auto& view : sfm_data.GetViews())
  {
    map_view[view.first] = scene.images.size();
    MVS::Interface::Image image;
    const std::string srcImage = stlplus::create_filespec(sfm_data.s_root_path, view.second->s_Img_path);
    image.name = stlplus::create_filespec(sOutDir, view.second->s_Img_path);
    image.platformID = map_intrinsic.at(view.second->id_intrinsic);
    MVS::Interface::Platform& platform = scene.platforms[image.platformID];
    image.cameraID = 0;

//    bool b1 = view.second.get()->id_intrinsic != UndefinedIndexT;
//    bool b2 = view.second.get()->id_pose != UndefinedIndexT;
//    bool b3 = sfm_data.intrinsics.find(view.second.get()->id_intrinsic) != sfm_data.intrinsics.end();
//    bool b4 = sfm_data.poses.find(view.second.get()->id_pose) != sfm_data.poses.end();
//    std::cout << b1 << b2 << b3 << b4 << std::endl;

    if (sfm_data.IsPoseAndIntrinsicDefined(view.second.get()) )//&& stlplus::is_file(srcImage))
    {
      MVS::Interface::Platform::Pose pose;
      image.poseID = platform.poses.size();
      const openMVG::geometry::Pose3 poseMVG(sfm_data.GetPoseOrDie(view.second.get()));
      pose.R = poseMVG.rotation();
      pose.C = poseMVG.center();

      if(undistortImage)
      {
          // export undistorted images
          const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view.second->id_intrinsic).get();
          if (cam->have_disto())
          {
            // undistort image and save it
            Image<openMVG::image::RGBColor> imageRGB, imageRGB_ud;
            ReadImage(srcImage.c_str(), &imageRGB);
            UndistortImage(imageRGB, cam, imageRGB_ud, view.second->id_intrinsic, map_d_list);
            WriteImage(image.name.c_str(), imageRGB_ud);
//            std::cout << "error 1" << std::endl;
          }
          else
          {
//              std::cout << "error 2" << std::endl;
            // just copy image
            stlplus::file_copy(srcImage, image.name);
          }
      }
      else
      {
          if(!single_mvs)
          {
              // just copy image
              //stlplus::file_copy(srcImage, image.name);
          }
      }

      platform.poses.push_back(pose);
      ++nPoses;
    }
    else
    {
        std::cout << "image " << view.first << " have not valid pose" << std::endl;
        // image have not valid pose, so set an undefined pose
        image.poseID = NO_ID;
        if(!single_mvs)
        {
            // just copy the image
            stlplus::file_copy(srcImage, image.name);
        }
    }
    scene.images.emplace_back(image);
    ++my_progress_bar;
  }

  // define structure
  scene.vertices.reserve(sfm_data.GetLandmarks().size());
  for (const auto& vertex: sfm_data.GetLandmarks())
  {
    const Landmark & landmark = vertex.second;
    MVS::Interface::Vertex vert;
    MVS::Interface::Vertex::ViewArr& views = vert.views;
    for (const auto& observation: landmark.obs)
    {
      const auto it(map_view.find(observation.first));
      if (it != map_view.end()) {
        MVS::Interface::Vertex::View view;
        view.imageID = it->second;
        view.confidence = 0;
        views.push_back(view);
      }
    }
    if (views.size() < 2)
      continue;
    std::sort(
      views.begin(), views.end(),
      [] (const MVS::Interface::Vertex::View& view0, const MVS::Interface::Vertex::View& view1)
      {
        return view0.imageID < view1.imageID;
      }
    );
    vert.X = landmark.X.cast<float>();
    scene.vertices.push_back(vert);
  }

  // normalize camera intrinsics
  for (size_t p=0; p<scene.platforms.size(); ++p)
  {
    MVS::Interface::Platform& platform = scene.platforms[p];
    for (size_t c=0; c<platform.cameras.size(); ++c) {
      MVS::Interface::Platform::Camera& camera = platform.cameras[c];
      // find one image using this camera
      MVS::Interface::Image* pImage(nullptr);
      for (MVS::Interface::Image& image: scene.images)
      {
        if (image.platformID == p && image.cameraID == c && image.poseID != NO_ID)
        {
          pImage = &image;
          break;
        }
      }
      if (pImage == nullptr)
      {
        std::cerr << "error: no image using camera " << c << " of platform " << p << std::endl;
        continue;
      }
      // read image meta-data
      ImageHeader imageHeader;
      ReadImageHeader(pImage->name.c_str(), &imageHeader);
      const double fScale(1.0/std::max(imageHeader.width, imageHeader.height));
      camera.K(0, 0) *= fScale;
      camera.K(1, 1) *= fScale;
      camera.K(0, 2) *= fScale;
      camera.K(1, 2) *= fScale;
    }
  }

  // write OpenMVS data
  if (!MVS::ARCHIVE::SerializeSave(scene, sOutFile))
    return false;

  std::cout
    << "Scene saved to OpenMVS interface format:\n"
    << "\t" << scene.images.size() << " images (" << nPoses << " calibrated)\n"
    << "\t" << scene.vertices.size() << " Landmarks\n";
  return true;
}

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
  cmd.add( make_option('o', sOutFile, "outfile") );
  cmd.add( make_option('d', sOutDir, "outdir") );
  cmd.add( make_option('f', sInputFileName_eachCamInfo, "input file name, file EachCameraInfo.txt" ) );
  cmd.add( make_option('g', sInputFileName_group, "input file name, file dmapReference.res" ) );
  cmd.add( make_option('u', undistortImage, "undistort image") );
  cmd.add( make_option('w', needTransWholeScene, "need trans whole Scene from sfm_data.bin to Scene.mvs") );
  cmd.add( make_option('t', forSparseScene, "for doing sparse scene mesh and texture.") );

  std::cout << "error 0 " << std::endl;
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

std::cout << "error 1 " << std::endl;
  if ((sInputFileName_group != "") || (sInputFileName_eachCamInfo != ""))
  {
      std::cout << "error 2 " << std::endl;
      std::cout<<"sInputFileName_group : " << sInputFileName_group << std::endl;

      std::vector<std::string> allFilePath;
      std::string wholeScenePath;

      if(sInputFileName_group != "")
      {
          std::ifstream in_group;  //dmapReference.res
          in_group.open(sInputFileName_group);
          if(in_group)
          {
              int count;
              in_group >> count;
              int line = 0;
              while(line<count && !in_group.eof())
              {
                  std::string path;
                  int a,b;
                  in_group >> a >> b >> path;

                  int refCount;
                  in_group >> refCount;
                  int num = 0;
                  while(num < refCount && !in_group.eof())
                  {
                      int tempRef;
                      in_group >> tempRef;

                      ++num;
                  }

                  if((allFilePath.size() == 0) || (allFilePath[allFilePath.size()-1] != path))
                  {
                      std::cout << path << std::endl;
                      allFilePath.push_back(path);
                  }

                  ++line;
              }
          }
          int temp1, temp2;
          in_group >> temp1 >> temp2 >> wholeScenePath;
          std::cout << wholeScenePath << std::endl;
          in_group.close();

      }



      if(sInputFileName_eachCamInfo != "")
      {
          std::ifstream in_eachCamInfo;  //EachCameraInfo.txt
          in_eachCamInfo.open(sInputFileName_eachCamInfo);
          if(in_eachCamInfo)
          {
              int count;
              in_eachCamInfo >> count;
              int line = 0;
              while(line<count && !in_eachCamInfo.eof())
              {
                  std::string path;
                  int a,b;
                  in_eachCamInfo >> a >> b >> path;
                  allFilePath.push_back(path);

                  int refCount;
                  in_eachCamInfo >> refCount;
                  int num = 0;
                  while(num < refCount && !in_eachCamInfo.eof())
                  {
                      int tempRef;
                      in_eachCamInfo >> tempRef;

                      ++num;
                  }

                  ++line;
              }
          }
          int temp1, temp2;
          in_eachCamInfo >> temp1 >> temp2 >> wholeScenePath;
          std::cout << wholeScenePath << std::endl;
          in_eachCamInfo.close();
      }

      //sorting and delete
      std::sort(allFilePath.begin(),allFilePath.end());
      allFilePath.erase(std::unique(allFilePath.begin(), allFilePath.end()), allFilePath.end());

      ///transform whole sfm_data.bin
      if(needTransWholeScene)
      {
          SfM_Data whole_sfm_data;
          if (!Load(whole_sfm_data, wholeScenePath+".bin", ESfM_Data(ALL))) {
            std::cerr << std::endl
              << "The input SfM_Data file \""<< wholeScenePath << "\" cannot be read." << std::endl;
            return EXIT_FAILURE;
          }


          Save(
              whole_sfm_data,
              stlplus::create_filespec( whole_sfm_data.s_root_path, "sfm_data_wholeScenePath.json" ).c_str(),
              ESfM_Data(INTRINSICS));

          std::cout << "whole_sfm_data.s_root_path+/Dense/Scene.mvs  " << whole_sfm_data.s_root_path+"/Dense/Scene.mvs" << std::endl;
          exportToOpenMVS(whole_sfm_data, sOutDir+"/Scene.mvs", sOutDir, undistortImage, false);
      }


#pragma omp parallel for
      for (std::size_t i = 0; i < allFilePath.size(); ++i)
      {
          std::cout<<allFilePath[i]+".bin" << std::endl;
          // Read
          SfM_Data sfm_data;
          if (!Load(sfm_data, allFilePath[i]+".bin", ESfM_Data(ALL))) {
            std::cerr << std::endl
              << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
          }

//          int p = allFilePath[i].find_last_of("/");
//          std::string ptr = allFilePath[i].substr(p);//;,allFilePath[i].size());
          //std::cout << ptr << std::endl;

          exportToOpenMVS(sfm_data, allFilePath[i]+".mvs", sOutDir, false, true);
      }

  }else if(forSparseScene)
  {
      std::cout << "error 3 " << std::endl;
      std::cout << "for doing Sparse Scene" << std::endl;

      // Read the input SfM scene
      SfM_Data sfm_data;
      if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
        std::cerr << std::endl
          << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
        return EXIT_FAILURE;
      }

      SfM_Data new_sfm_data;
      new_sfm_data.intrinsics[1] = sfm_data.intrinsics[1];
      new_sfm_data.s_root_path = sfm_data.s_root_path;
      new_sfm_data.structure = sfm_data.structure;
      //
      for(Views::const_iterator itView = sfm_data.views.begin(); itView != sfm_data.views.end(); ++itView)
      {
          IndexT id = itView->first;
          std::stringstream ss;
          ss << id;
          std::string s;
          ss >> s;
          std::cout << "id view : " << s << std::endl;
          if(s[0] == '2')
          {
              new_sfm_data.views[id] = itView->second;
          }
      }
      //
      for(Poses::const_iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose)
      {
          IndexT id = itPose->first;
          std::stringstream ss;
          ss << id;
          std::string s;
          ss >> s;
          if(s[0] == '2')
          {
              std::cout << "id : " << s << std::endl;
              new_sfm_data.poses[id] = itPose->second;
          }
      }
      //
//      for(Landmarks::const_iterator itLand = sfm_data.structure.begin(); itLand != sfm_data.structure.end(); ++itLand)
//      {

//          for(Observations::const_iterator itObs = itLand->second.obs.begin(); itObs != itLand->second.obs.end(); ++itObs)
//          {
//              std::cout << "id obs : " << itObs->first << std::endl;
//              IndexT id = itObs->first;
//              std::stringstream ss;
//              ss << id;
//              std::string s;
//              ss >> s;
//              if(s[0] == '2')
//              {
//                    new_sfm_data.structure[itLand->first].obs[itObs->first] = itObs->second;
//              }

//          }
//          if(new_sfm_data.structure[itLand->first].obs.size() != 0)
//          {
//              new_sfm_data.structure[itLand->first].X = itLand->second.X;
//          }
//      }



      if (exportToOpenMVS(new_sfm_data, sOutFile, sOutDir, undistortImage, false))
        return( EXIT_SUCCESS );

  }else{

      std::cout << "error 4 " << std::endl;
      std::cout << "without sInputFileName" << std::endl;
      std::cout << undistortImage << std::endl;
      std::cout<<"sSfM_Data_Filename : " << sSfM_Data_Filename << std::endl;
      std::cout << "sOutFile : " << sOutFile << std::endl;
      std::cout << "sOutDir : " << sOutDir << std::endl;

      // Read the input SfM scene
      SfM_Data sfm_data;
      if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
        std::cerr << std::endl
          << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
        return EXIT_FAILURE;
      }

      if (exportToOpenMVS(sfm_data, sOutFile, sOutDir, undistortImage, false))
        return( EXIT_SUCCESS );


  }


  return EXIT_SUCCESS;
}
