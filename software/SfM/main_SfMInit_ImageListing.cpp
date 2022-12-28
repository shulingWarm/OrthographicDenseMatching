// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/cameras.hpp"
#include "openMVG/exif/exif_IO_EasyExif.hpp"
#include "openMVG/exif/sensor_width_database/ParseDatabase.hpp"
#include "openMVG/geodesy/geodesy.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_utils.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/sfm/sfm_view_priors.hpp"
#include "openMVG/types.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress_display.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <fstream>
#include <memory>
#include <string>
#include <utility>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::exif;
using namespace openMVG::geodesy;
using namespace openMVG::image;
using namespace openMVG::sfm;

/// Check that Kmatrix is a string like "f;0;ppx;0;f;ppy;0;0;1"
/// With f,ppx,ppy as valid numerical value
bool checkIntrinsicStringValidity(const std::string & Kmatrix, double & focal, double & ppx, double & ppy)
{
  std::vector<std::string> vec_str;
  stl::split(Kmatrix, ';', vec_str);
  if (vec_str.size() != 9)  {
    std::cerr << "\n Missing ';' character" << std::endl;
    return false;
  }
  // Check that all K matrix value are valid numbers
  for (size_t i = 0; i < vec_str.size(); ++i) {
    double readvalue = 0.0;
    std::stringstream ss;
    ss.str(vec_str[i]);
    if (! (ss >> readvalue) )  {
      std::cerr << "\n Used an invalid not a number character" << std::endl;
      return false;
    }
    if (i==0) focal = readvalue;
    if (i==2) ppx = readvalue;
    if (i==5) ppy = readvalue;
  }
  return true;
}

std::pair<bool, Vec3> checkGPS
(
  const std::string & filename,
  const int & GPS_to_XYZ_method = 0
)
{
  std::pair<bool, Vec3> val(false, Vec3::Zero());
  std::unique_ptr<Exif_IO> exifReader(new Exif_IO_EasyExif);
  if (exifReader)
  {
    // Try to parse EXIF metada & check existence of EXIF data
    if ( exifReader->open( filename ) && exifReader->doesHaveExifInfo() )
    {
      // Check existence of GPS coordinates
      double latitude, longitude, altitude;
      if ( exifReader->GPSLatitude( &latitude ) &&
           exifReader->GPSLongitude( &longitude ) &&
           exifReader->GPSAltitude( &altitude ) )
      {
        // Add ECEF or UTM XYZ position to the GPS position array
        val.first = true;
        switch (GPS_to_XYZ_method)
        {
          case 1:
            val.second = lla_to_utm( latitude, longitude, altitude );
            break;
          case 0:
          default:
            val.second = lla_to_ecef( latitude, longitude, altitude );
            break;
        }
      }
    }
  }
  return val;
}


/// Check string of prior weights
std::pair<bool, Vec3> checkPriorWeightsString
(
  const std::string &sWeights
)
{
  std::pair<bool, Vec3> val(true, Vec3::Zero());
  std::vector<std::string> vec_str;
  stl::split(sWeights, ';', vec_str);
  if (vec_str.size() != 3)
  {
    std::cerr << "\n Missing ';' character in prior weights" << std::endl;
    val.first = false;
  }
  // Check that all weight values are valid numbers
  for (size_t i = 0; i < vec_str.size(); ++i)
  {
    double readvalue = 0.0;
    std::stringstream ss;
    ss.str(vec_str[i]);
    if (! (ss >> readvalue) )  {
      std::cerr << "\n Used an invalid not a number character in local frame origin" << std::endl;
      val.first = false;
    }
    val.second[i] = readvalue;
  }
  return val;
}

//赵志豪 2021-1-29
//从一组二进制文件中读取光芯坐标,每3个为一组,用于读取先验坐标
void readPriorInfo(
        std::string fileName,//文件在硬盘中的路径
        std::vector<Eigen::Vector3d> &poseVec //读取后得到的位姿会被存储在此
        )
{
    //读取目标文本
    std::ifstream binHandle;
    binHandle.open(fileName,std::ios::binary);
    //读取总共有多少数据
    int dataNum;
    binHandle.read((char*)(&dataNum),sizeof(int));
    //向量预先开辟空间
    poseVec.reserve(dataNum);
    //遍历所有需要读取的数据
    for(int dataCount=0;dataCount<dataNum;++dataCount)
    {
        //新建最后的向量结果
        Eigen::Vector3d tempVec;
        //遍历向量的每个坐标
        for(int dimCount=0;dimCount<3;++dimCount)
        {
            //读取当前位置的数据
            double dimData;
            binHandle.read((char*)(&dimData),sizeof(double));
            //记录当前位置的数据
            tempVec(dimCount)=dimData;
        }
        //记录当前位置的光心坐标
        poseVec.push_back(tempVec);
    }
    //关闭文本
    binHandle.close();
}


//
// Create the description of an input image dataset for OpenMVG toolsuite
// - Export a SfM_Data file with View & Intrinsic data
//
int main(int argc, char **argv)
{
    //是否使用标定后得到的外参,如果使用的话,数据会从二进制文件里面读取
    bool useStandardExtrinsic=true;
    //二进制形式的外参数据的路径
    std::string extPosePath="/home/cvlab/workSpace/mainProject/topViewConstruct/extrinsicStandard/noWaterCenter.bin";


  CmdLine cmd;

  std::string sImageDir,
    sfileDatabase = "",
    sOutputDir = "",
    sKmatrix;

  std::string sPriorWeights;
  std::pair<bool, Vec3> prior_w_info(false, Vec3(1.0,1.0,1.0));

  int i_User_camera_model = PINHOLE_CAMERA_RADIAL3;

  bool b_Group_camera_model = true;

  int i_GPS_XYZ_method = 0;

  double focal_pixels = -1.0;

  cmd.add( make_option('i', sImageDir, "imageDirectory") );
  cmd.add( make_option('d', sfileDatabase, "sensorWidthDatabase") );
  cmd.add( make_option('o', sOutputDir, "outputDirectory") );
  cmd.add( make_option('f', focal_pixels, "focal") );
  cmd.add( make_option('k', sKmatrix, "intrinsics") );
  cmd.add( make_option('c', i_User_camera_model, "camera_model") );
  cmd.add( make_option('g', b_Group_camera_model, "group_camera_model") );
  cmd.add( make_switch('P', "use_pose_prior") );
  cmd.add( make_option('W', sPriorWeights, "prior_weigths"));
  cmd.add( make_option('m', i_GPS_XYZ_method, "gps_to_xyz_method") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch (const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--imageDirectory]\n"
      << "[-d|--sensorWidthDatabase]\n"
      << "[-o|--outputDirectory]\n"
      << "[-f|--focal] (pixels)\n"
      << "[-k|--intrinsics] Kmatrix: \"f;0;ppx;0;f;ppy;0;0;1\"\n"
      << "[-c|--camera_model] Camera model type:\n"
      << "\t 1: Pinhole\n"
      << "\t 2: Pinhole radial 1\n"
      << "\t 3: Pinhole radial 3 (default)\n"
      << "\t 4: Pinhole brown 2\n"
      << "\t 5: Pinhole with a simple Fish-eye distortion\n"
      << "[-g|--group_camera_model]\n"
      << "\t 0-> each view have it's own camera intrinsic parameters,\n"
      << "\t 1-> (default) view can share some camera intrinsic parameters\n"
      << "\n"
      << "[-P|--use_pose_prior] Use pose prior if GPS EXIF pose is available"
      << "[-W|--prior_weigths] \"x;y;z;\" of weights for each dimension of the prior (default: 1.0)\n"
      << "[-m|--gps_to_xyz_method] XZY Coordinate system:\n"
      << "\t 0: ECEF (default)\n"
      << "\t 1: UTM\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--imageDirectory " << sImageDir << std::endl
            << "--sensorWidthDatabase " << sfileDatabase << std::endl
            << "--outputDirectory " << sOutputDir << std::endl
            << "--focal " << focal_pixels << std::endl
            << "--intrinsics " << sKmatrix << std::endl
            << "--camera_model " << i_User_camera_model << std::endl
            << "--group_camera_model " << b_Group_camera_model << std::endl;

  // Expected properties for each image
  double width = -1, height = -1, focal = -1, ppx = -1,  ppy = -1;

  const EINTRINSIC e_User_camera_model = EINTRINSIC(i_User_camera_model);

  if ( !stlplus::folder_exists( sImageDir ) )
  {
    std::cerr << "\nThe input directory doesn't exist" << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutputDir.empty())
  {
    std::cerr << "\nInvalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if ( !stlplus::folder_exists( sOutputDir ) )
  {
    if ( !stlplus::folder_create( sOutputDir ))
    {
      std::cerr << "\nCannot create output directory" << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (sKmatrix.size() > 0 &&
    !checkIntrinsicStringValidity(sKmatrix, focal, ppx, ppy) )
  {
    std::cerr << "\nInvalid K matrix input" << std::endl;
    return EXIT_FAILURE;
  }

  if (sKmatrix.size() > 0 && focal_pixels != -1.0)
  {
    std::cerr << "\nCannot combine -f and -k options" << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<Datasheet> vec_database;
  if (!sfileDatabase.empty())
  {
    if ( !parseDatabase( sfileDatabase, vec_database ) )
    {
      std::cerr
       << "\nInvalid input database: " << sfileDatabase
       << ", please specify a valid file." << std::endl;
      return EXIT_FAILURE;
    }
  }

  // Check if prior weights are given
  if (cmd.used('P') && !sPriorWeights.empty())
  {
    prior_w_info = checkPriorWeightsString(sPriorWeights);
  }
  else if (cmd.used('P'))
  {
    prior_w_info.first = true;
  }

  std::vector<std::string> vec_image = stlplus::folder_files( sImageDir );
  std::sort(vec_image.begin(), vec_image.end());

  // Configure an empty scene with Views and their corresponding cameras
  SfM_Data sfm_data;
  sfm_data.s_root_path = sImageDir; // Setup main image root_path
  Views & views = sfm_data.views;
  Intrinsics & intrinsics = sfm_data.intrinsics;

  C_Progress_display my_progress_bar( vec_image.size(),
      std::cout, "\n- Image listing -\n" );
  std::ostringstream error_report_stream;
  //如果需要使用二进制的外参数据,需要提前在这里读取到内存里
  std::vector<Eigen::Vector3d> poseVector;
  if(useStandardExtrinsic)
  {
      readPriorInfo(extPosePath,poseVector);
  }
  for ( std::vector<std::string>::const_iterator iter_image = vec_image.begin();
    iter_image != vec_image.end();
    ++iter_image, ++my_progress_bar )
  {
    // Read meta data to fill camera parameter (w,h,focal,ppx,ppy) fields.
    width = height = ppx = ppy = focal = -1.0;

    //判断是否需要参考来自外部的外参
    if(useStandardExtrinsic)
    {
        //判断是否还有可用的外参
        if(views.size()>=poseVector.size())
            useStandardExtrinsic=false;
    }

    const std::string sImageFilename = stlplus::create_filespec( sImageDir, *iter_image );
    const std::string sImFilenamePart = stlplus::filename_part(sImageFilename);

    // Test if the image format is supported:
    if (openMVG::image::GetFormat(sImageFilename.c_str()) == openMVG::image::Unknown)
    {
      error_report_stream
          << sImFilenamePart << ": Unkown image file format." << "\n";
      continue; // image cannot be opened
    }

    if (sImFilenamePart.find("mask.png") != std::string::npos
       || sImFilenamePart.find("_mask.png") != std::string::npos)
    {
      error_report_stream
          << sImFilenamePart << " is a mask image" << "\n";
      continue;
    }

    ImageHeader imgHeader;
    if (!openMVG::image::ReadImageHeader(sImageFilename.c_str(), &imgHeader))
      continue; // image cannot be read

    width = imgHeader.width;
    height = imgHeader.height;
    ppx = width / 2.0;
    ppy = height / 2.0;

    std::unique_ptr<Exif_IO> exifReader(new Exif_IO_EasyExif);
    exifReader->open( sImageFilename );

    const bool bHaveValidExifMetadata =
      exifReader->doesHaveExifInfo()
      && !exifReader->getModel().empty();

    // Consider the case where the focal is provided manually
    if ( !bHaveValidExifMetadata || focal_pixels != -1)
    {
      if (sKmatrix.size() > 0) // Known user calibration K matrix
      {
        if (!checkIntrinsicStringValidity(sKmatrix, focal, ppx, ppy))
          focal = -1.0;
      }
      else // User provided focal length value
        if (focal_pixels != -1 )
          focal = focal_pixels;
    }
    else // If image contains meta data
    {
      const std::string sCamModel = exifReader->getModel();

      // Handle case where focal length is equal to 0
      if (exifReader->getFocal() == 0.0f)
      {
        error_report_stream
          << stlplus::basename_part(sImageFilename) << ": Focal length is missing." << "\n";
        focal = -1.0;
      }
      else
      // Create the image entry in the list file
      {
        Datasheet datasheet;
        if ( getInfo( sCamModel, vec_database, datasheet ))
        {
          // The camera model was found in the database so we can compute it's approximated focal length
          const double ccdw = datasheet.sensorSize_;
          focal = std::max ( width, height ) * exifReader->getFocal() / ccdw;
        }
        else
        {
          error_report_stream
            << stlplus::basename_part(sImageFilename)
            << "\" model \"" << sCamModel << "\" doesn't exist in the database" << "\n"
            << "Please consider add your camera model and sensor width in the database." << "\n";
        }
      }
    }

    // Build intrinsic parameter related to the view
    std::shared_ptr<IntrinsicBase> intrinsic;

    if (focal > 0 && ppx > 0 && ppy > 0 && width > 0 && height > 0)
    {
      // Create the desired camera type
      switch (e_User_camera_model)
      {
        case PINHOLE_CAMERA:
          intrinsic = std::make_shared<Pinhole_Intrinsic>
            (width, height, focal, ppx, ppy);
        break;
        case PINHOLE_CAMERA_RADIAL1:
          intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K1>
            (width, height, focal, ppx, ppy, 0.0); // setup no distortion as initial guess
        break;
        case PINHOLE_CAMERA_RADIAL3:
          intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K3>
            (width, height, focal, ppx, ppy, 0.0, 0.0, 0.0);  // setup no distortion as initial guess
        break;
        case PINHOLE_CAMERA_BROWN:
          intrinsic =std::make_shared<Pinhole_Intrinsic_Brown_T2>
            (width, height, focal, ppx, ppy, 0.0, 0.0, 0.0, 0.0, 0.0); // setup no distortion as initial guess
        break;
        case PINHOLE_CAMERA_FISHEYE:
          intrinsic =std::make_shared<Pinhole_Intrinsic_Fisheye>
            (width, height, focal, ppx, ppy, 0.0, 0.0, 0.0, 0.0); // setup no distortion as initial guess
        break;
        default:
          std::cerr << "Error: unknown camera model: " << (int) e_User_camera_model << std::endl;
          return EXIT_FAILURE;
      }
    }

    // Build the view corresponding to the image
    const std::pair<bool, Vec3> gps_info = checkGPS(sImageFilename, i_GPS_XYZ_method);
    if (gps_info.first && cmd.used('P'))
    {
      ViewPriors v(*iter_image, views.size(), views.size(), views.size(), width, height);

      // Add intrinsic related to the image (if any)
      if (intrinsic == nullptr)
      {
        //Since the view have invalid intrinsic data
        // (export the view, with an invalid intrinsic field value)
        v.id_intrinsic = UndefinedIndexT;
      }
      else
      {
        // Add the defined intrinsic to the sfm_container
        intrinsics[v.id_intrinsic] = intrinsic;
      }

      v.b_use_pose_center_ = true;
      v.pose_center_ = gps_info.second;
      // prior weights
      if (prior_w_info.first == true)
      {
        v.center_weight_ = prior_w_info.second;
      }

      // Add the view to the sfm_container
      views[v.id_view] = std::make_shared<ViewPriors>(v);
    }
    else if(useStandardExtrinsic)//使用标定外参的情况
    {
        ViewPriors v(*iter_image, views.size(), views.size(), views.size(), width, height);

        // Add intrinsic related to the image (if any)
        if (intrinsic == nullptr)
        {
          //Since the view have invalid intrinsic data
          // (export the view, with an invalid intrinsic field value)
          v.id_intrinsic = UndefinedIndexT;
        }
        else
        {
          // Add the defined intrinsic to the sfm_container
          intrinsics[v.id_intrinsic] = intrinsic;
        }

        v.b_use_pose_center_ = true;
        v.pose_center_ = poseVector[v.id_view];
        //权值默认都是1
        v.center_weight_=Eigen::Vector3d(1,1,1);

        // Add the view to the sfm_container
        views[v.id_view] = std::make_shared<ViewPriors>(v);
    }
    else
    {
      View v(*iter_image, views.size(), views.size(), views.size(), width, height);

      // Add intrinsic related to the image (if any)
      if (intrinsic == nullptr)
      {
        //Since the view have invalid intrinsic data
        // (export the view, with an invalid intrinsic field value)
        v.id_intrinsic = UndefinedIndexT;
      }
      else
      {
        // Add the defined intrinsic to the sfm_container
        intrinsics[v.id_intrinsic] = intrinsic;
      }

      // Add the view to the sfm_container
      views[v.id_view] = std::make_shared<View>(v);
    }
  }

  // Display saved warning & error messages if any.
  if (!error_report_stream.str().empty())
  {
    std::cerr
      << "\nWarning & Error messages:" << std::endl
      << error_report_stream.str() << std::endl;
  }

  // Group camera that share common properties if desired (leads to more faster & stable BA).
  if (b_Group_camera_model)
  {
    GroupSharedIntrinsics(sfm_data);
  }
  std::cout << "error 1" << std::endl;


  ////add cameras for model c 4
  {
          // Build intrinsic parameter related to the view
          //Hash_Map<IndexT, std::vector<double> > mid_intrinsics;
          //mid_intrinsics[0] = sfm_data.intrinsics[0]->getParams();



//          //
          std::shared_ptr<IntrinsicBase> intrinsic0 (NULL);
          intrinsic0 = std::make_shared<Pinhole_Intrinsic_Brown_T2>
                  (7360.0, 4912.0, 10268.536000, 3657.861207, 2458.938399, 0.0, 0.0, 0.0, 0.0, 0.0);
//          (1280.0, 1024.0, 930.0, mid_intrinsics[0][1], mid_intrinsics[0][2], 0.0, 0.0, 0.0, 0.0, 0.0);
           sfm_data.intrinsics[0] = intrinsic0;

           std::cout << "error 1 1" << std::endl;
           //7175.0
          std::shared_ptr<IntrinsicBase> intrinsic1 (NULL);
          intrinsic1 = std::make_shared<Pinhole_Intrinsic_Brown_T2>
                  (7360.0, 4912.0, 7274.47, 3687.08, 2456.77, 0.0, 0.0, 0.0, 0.0, 0.0);
//          (1280.0, 1024.0, 930.0, mid_intrinsics[0][1], mid_intrinsics[0][2], 0.0, 0.0, 0.0, 0.0, 0.0);
           sfm_data.intrinsics[1] = intrinsic1;
           std::cout << "error 1 2" << std::endl;

           //10304
           std::shared_ptr<IntrinsicBase> intrinsic2 (NULL);
           intrinsic2 = std::make_shared<Pinhole_Intrinsic_Brown_T2>
                   (7360.0, 4912.0, 10287.45, 3722.57, 2539.4, 0.0, 0.0, 0.0, 0.0, 0.0);
//           (1280.0, 1024.0, 930.0, mid_intrinsics[0][1], mid_intrinsics[0][2], 0.0, 0.0, 0.0, 0.0, 0.0);
            sfm_data.intrinsics[2] = intrinsic2;
            std::cout << "error 1 3" << std::endl;

//            std::shared_ptr<IntrinsicBase> intrinsic0 (NULL);
//            intrinsic0 = std::make_shared<Pinhole_Intrinsic>
//                    (7360.0, 4912.0, 10307.76, 3680.0, 2456.0);
//             sfm_data.intrinsics[0] = intrinsic0;

//             //7175.0
//            std::shared_ptr<IntrinsicBase> intrinsic1 (NULL);
//            intrinsic1 = std::make_shared<Pinhole_Intrinsic>
//                    (7360.0, 4912.0, 7233.94, 3680.0, 2456.0);
//             sfm_data.intrinsics[1] = intrinsic1;

//             //10304
//             std::shared_ptr<IntrinsicBase> intrinsic2 (NULL);
//             intrinsic2 = std::make_shared<Pinhole_Intrinsic>
//                     (7360.0, 4912.0, 10279.8, 3680.0, 2456.0);
//              sfm_data.intrinsics[2] = intrinsic2;



           for(Views::iterator it = sfm_data.views.begin(); it != sfm_data.views.end(); ++it)
           {

              View * view = it->second.get();
              IndexT poseId = view->id_pose;
              if((poseId >= sfm_data.views.size()/3) && (poseId < sfm_data.views.size()/3*2))
              {
                  view->id_intrinsic = 1;
              }else if(poseId >= sfm_data.views.size()/3*2)
              {
                  view->id_intrinsic = 2;
              }else if(poseId < sfm_data.views.size()/3)
              {
                  view->id_intrinsic = 0;
              }
           }

  }
  std::cout << "error 2" << std::endl;

  // Store SfM_Data views & intrinsic data
  if (!Save(
    sfm_data,
    stlplus::create_filespec( sOutputDir, "sfm_data.json" ).c_str(),
    ESfM_Data(VIEWS|INTRINSICS)))
  {
    return EXIT_FAILURE;
  }

  std::cout << std::endl
    << "SfMInit_ImageListing report:\n"
    << "listed #File(s): " << vec_image.size() << "\n"
    << "usable #File(s) listed in sfm_data: " << sfm_data.GetViews().size() << "\n"
    << "usable #Intrinsic(s) listed in sfm_data: " << sfm_data.GetIntrinsics().size() << std::endl;

  return EXIT_SUCCESS;
}
