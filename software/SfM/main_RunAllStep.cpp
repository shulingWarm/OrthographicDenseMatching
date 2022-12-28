// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// The <cereal/archives> headers are special and must be included first.
#include <cereal/archives/json.hpp>

#include "openMVG/features/image_describer_akaze_io.hpp"

#include "openMVG/features/sift/SIFT_Anatomy_Image_Describer_io.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/features/regions_factory_io.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress_display.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include "nonFree/sift/SIFT_describer_io.hpp"

#include <cereal/details/helpers.hpp>


#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"
#include "openMVG/matching/indMatch_io.hpp"

#include "openMVG/graph/graph.hpp"
#include "openMVG/features/akaze/image_describer_akaze.hpp"
#include "openMVG/features/descriptor.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/matching_image_collection/Matcher_Regions.hpp"
#include "openMVG/matching_image_collection/Cascade_Hashing_Matcher_Regions.hpp"
#include "openMVG/matching_image_collection/GeometricFilter.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider_cache.hpp"
#include "openMVG/matching_image_collection/F_ACRobust.hpp"
#include "openMVG/matching_image_collection/E_ACRobust.hpp"
#include "openMVG/matching_image_collection/H_ACRobust.hpp"
#include "openMVG/matching_image_collection/Pair_Builder.hpp"
#include "openMVG/matching/pairwiseAdjacencyDisplay.hpp"
#include "openMVG/stl/stl.hpp"

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/cameras/Cameras_Common_command_line_helper.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"
#include "openMVG/sfm/pipelines/global/sfm_global_engine_relative_motions.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_report.hpp"

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/cameras/Camera_undistort_image.hpp"

//#include "openMVG/cameras/cameras.hpp"
//#include "openMVG/geodesy/geodesy.hpp"
//#include "openMVG/image/image_io.hpp"
//#include "openMVG/numeric/eigen_alias_definition.hpp"
//#include "openMVG/sfm/sfm_data_utils.hpp"
//#include "openMVG/sfm/sfm_view.hpp"
//#include "openMVG/sfm/sfm_view_priors.hpp"
//#include "openMVG/types.hpp"

#include <cstdlib>
#include <fstream>
#include <string>
#include <iostream>

#include <sys/time.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif

using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::features;
using namespace openMVG::sfm;
using namespace std;

features::EDESCRIBER_PRESET stringToEnum(const std::string & sPreset)
{
  features::EDESCRIBER_PRESET preset;
  if (sPreset == "NORMAL")
    preset = features::NORMAL_PRESET;
  else
  if (sPreset == "HIGH")
    preset = features::HIGH_PRESET;
  else
  if (sPreset == "ULTRA")
    preset = features::ULTRA_PRESET;
  else
    preset = features::EDESCRIBER_PRESET(-1);
  return preset;
}

enum EGeometricModel
{
  FUNDAMENTAL_MATRIX = 0,
  ESSENTIAL_MATRIX   = 1,
  HOMOGRAPHY_MATRIX  = 2
};

enum EPairMode
{
  PAIR_EXHAUSTIVE = 0,
  PAIR_CONTIGUOUS = 1,
  PAIR_FROM_FILE  = 2
};

double get_wall_time()
{
    struct timeval time ;
    if (gettimeofday(&time,NULL)){
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

/// - Compute view image description (feature & descriptor extraction)
/// - Export computed data
int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sOutDir = "";
  bool bUpRight = false;
  std::string sImage_Describer_Method = "SIFT";
  bool bForce = false;
  std::string sFeaturePreset = "";
  std::string sPredefinedPairList = "";
  std::string sMatchesDirectory = "";
  float fDistRatio = 0.8f;
  std::string sIntrinsic_refinement_options = "NONE";//"ADJUST_ALL";
  std::string rcFilePath = "";
  std::string undistortFilesPath = "";
  std::string sNeedundistort = "";
  int iFeatureLimit = 1500;
  bool bNeedGenerateMask = true;
  bool bNeedOneSixteenth = true;

  const int iRotationAveragingMethod = int (ROTATION_AVERAGING_L2);
  const int iTranslationAveragingMethod = int (TRANSLATION_AVERAGING_SOFTL1);
#ifdef OPENMVG_USE_OPENMP
  int iNumThreads = 60;
#endif

  // required
  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('l', sPredefinedPairList, "pair_list") );
  cmd.add( make_option('r', fDistRatio, "ratio") );
  cmd.add( make_option('f', sIntrinsic_refinement_options, "refineIntrinsics") );
  cmd.add( make_option('c', rcFilePath, "file path of R t file, n_gps_imu_321.res") );
  cmd.add( make_option('d', sNeedundistort, "needundistort") );
  cmd.add( make_option('s', iFeatureLimit, "feature points limit") );
  cmd.add( make_option('e', bNeedGenerateMask, "needMask") );
  cmd.add( make_option('g', bNeedOneSixteenth, "needOneSixteenth") );

  // Optional
  cmd.add( make_option('m', sImage_Describer_Method, "describerMethod") );
  cmd.add( make_option('u', bUpRight, "upright") );
//  cmd.add( make_option('f', bForce, "force") );
  cmd.add( make_option('p', sFeaturePreset, "describerPreset") );

#ifdef OPENMVG_USE_OPENMP
  cmd.add( make_option('n', iNumThreads, "numThreads") );
#endif

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch (const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--input_file] a SfM_Data file \n"
      << "[-o|--outdir path] \n"
      << "\n[Optional]\n"
      << "[-f|--force] Force to recompute data\n"
      << "[-m|--describerMethod]\n"
      << "  (method to use to describe an image):\n"
      << "   SIFT (default),\n"
      << "   SIFT_ANATOMY,\n"
      << "   AKAZE_FLOAT: AKAZE with floating point descriptors,\n"
      << "   AKAZE_MLDB:  AKAZE with binary descriptors\n"
      << "[-u|--upright] Use Upright feature 0 or 1\n"
      << "[-p|--describerPreset]\n"
      << "  (used to control the Image_describer configuration):\n"
      << "   NORMAL (default),\n"
      << "   HIGH,\n"
      << "   ULTRA: !!Can take long time!!\n"
      << "[-e|--needMask] If need generate a mask\n"
      << "[-g|--needOneSixteenth] If need resize images into one-sixteenth\n"
#ifdef OPENMVG_USE_OPENMP
      << "[-n|--numThreads] number of parallel computations\n"
#endif
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--input_file " << sSfM_Data_Filename << std::endl
            << "--outdir " << sOutDir << std::endl
            << "--describerMethod " << sImage_Describer_Method << std::endl
            << "--upright " << bUpRight << std::endl
            << "--describerPreset " << (sFeaturePreset.empty() ? "NORMAL" : sFeaturePreset) << std::endl
            << "--force " << bForce << std::endl
            << "--needMask " << bNeedGenerateMask << std::endl
            << "--needOneSixteenth " << bNeedOneSixteenth << std::endl
#ifdef OPENMVG_USE_OPENMP
            << "--numThreads " << iNumThreads << std::endl
#endif
            << std::endl;


  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
  {
    if (!stlplus::folder_create(sOutDir))
    {
      std::cerr << "Cannot create output directory" << std::endl;
      return EXIT_FAILURE;
    }
  }


  std::shared_ptr<Regions_Provider> regions_provider; ///for step compute-matches
  const std::string sImage_describer = stlplus::create_filespec(sOutDir, "image_describer", "json");
  //---------------------------------------
  // a. Load input scene
  //---------------------------------------
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS))) {
    std::cerr << std::endl
      << "The input file \""<< sSfM_Data_Filename << "\" cannot be read" << std::endl;
    return false;
  }

  std::cout << "load " << sSfM_Data_Filename << " file finish" << std::endl;

  /// undistort images
  /// prepare step
//  cv::Mat map_x[3], map_y[3];
  std::map<int, std::pair<std::vector<std::vector<double> >, std::vector<std::vector<double> > > > map_d_list;
  {

      if(sNeedundistort != "")
      {
          SfM_Data intrinsic_data;
          if (!Load(intrinsic_data, sNeedundistort, ESfM_Data(INTRINSICS))) {
            std::cerr << std::endl
              << "The input file \""<< sSfM_Data_Filename << "\" cannot be read" << std::endl;
            std::cout << "load " << sNeedundistort << " file failed!" << std::endl;
            return false;
          }

          const uint32_t images_w = sfm_data.views.begin()->second->ui_width;
          const uint32_t images_h = sfm_data.views.begin()->second->ui_height;
          std::cout << "images_w : " << images_w << " images_h : " << images_h << std::endl;

          for(Intrinsics::const_iterator itI = intrinsic_data.intrinsics.begin();
              itI != intrinsic_data.intrinsics.end(); ++itI)
          {
              std::vector<std::vector<double> > tempMapD_x(images_w, std::vector<double>(images_h));
              std::vector<std::vector<double> > tempMapD_y(images_w, std::vector<double>(images_h));

              #pragma omp parallel for
              for ( uint32_t j = 0; j < images_h; ++j )
                for ( uint32_t i = 0; i < images_w; ++i )
                {
                    const openMVG::cameras::IntrinsicBase * cam = intrinsic_data.GetIntrinsics().at(itI->first).get();
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
          std::cout << "preparing undistort map list finished" << std::endl;
      }

//               const int width = sfm_data.views.begin()->second.get()->ui_width;
//               const int height = sfm_data.views.begin()->second.get()->ui_height;


//               map_x[0].create(cv::Size(width, height), CV_32FC1);
//               map_y[0].create(cv::Size(width, height), CV_32FC1);
//               {
//                   std::string configfileB = undistortFilesPath+"/datb.conf";//stlplus::create_filespec(undistortFilesPath, "datb.conf");
//                   std::ifstream conf;
//                   conf.open(configfileB.c_str());//, std::ios::binary);
//                   if(!conf)
//                   {
//                       std::cerr << "open file " << configfileB << " failed!" << std::endl;
//                       return EXIT_FAILURE;
//                   }
//                   std::string lx_b;
//                   while(getline(conf, lx_b))
//                   {
//                       std::stringstream line(lx_b);
//                       int x, y;
//                       double sx, sy;
//                       line >> x >> y >> sx >> sy;
//                       map_x[0].at<float>(y, x) = sx;
//                       map_y[0].at<float>(y, x) = sy;
//                   }
//                   conf.close();
//                   std::cout << "load undistort file datb.conf finished" << std::endl;
//               }
//               //
//               map_x[1].create(cv::Size(width, height), CV_32FC1);
//               map_y[1].create(cv::Size(width, height), CV_32FC1);
//               {
//                   std::string configfileD = undistortFilesPath+ "/datd.conf";//stlplus::create_filespec(undistortFilesPath, "datd.conf");
//                   std::ifstream conf;
//                   conf.open(configfileD.c_str());//, std::ios::binary);
//                   if(!conf)
//                   {
//                       std::cerr << "open file " << configfileD << " failed!" << std::endl;
//                       return EXIT_FAILURE;
//                   }
//                   std::string lx_d;
//                   while(getline(conf, lx_d))
//                   {
//                       std::stringstream line(lx_d);
//                       int x, y;
//                       double sx, sy;
//                       line >> x >> y >> sx >> sy;
//                       map_x[1].at<float>(y, x) = sx;
//                       map_y[1].at<float>(y, x) = sy;
//                   }
//                   conf.close();
//                   std::cout << "load undistort file datd.conf finished" << std::endl;
//               }

//               //
//               map_x[2].create(cv::Size(width, height), CV_32FC1);
//               map_y[2].create(cv::Size(width, height), CV_32FC1);
//               {
//                   std::string configfileF = undistortFilesPath+"datf.conf";//stlplus::create_filespec(undistortFilesPath, "datf.conf");
//                   std::ifstream conf;//, std::ios::binary);
//                   conf.open(configfileF.c_str(), std::ios::binary);
//                   if(!conf)
//                   {
//                       std::cerr << "open file " << configfileF << " failed!" << std::endl;
//                       return EXIT_FAILURE;
//                   }
//                   std::string lx_f;
//                   while(getline(conf, lx_f))
//                   {
//                       std::stringstream line(lx_f);
//                       int x, y;
//                       double sx, sy;
//                       line >> x >> y >> sx >> sy;
//                       map_x[2].at<float>(y, x) = sx;
//                       map_y[2].at<float>(y, x) = sy;
//                   }
//                   conf.close();
//                   std::cout << "load undistort file datf.conf finished" << std::endl;
//               }

  }


  ///Compute Feature ***
  ///
  ///
  ///
  ///
  /// ******************
  {
      std::vector<int> reFeatureIdList;
      {
          // b. Init the image_describer
          // - retrieve the used one in case of pre-computed features
          // - else create the desired one
          using namespace openMVG::features;
          std::unique_ptr<Image_describer> image_describer;
          image_describer.reset(new SIFT_Image_describer
            (SIFT_Image_describer::Params(), !bUpRight));
          //

          // Export the used Image_describer and region type for:
          // - dynamic future regions computation and/or loading
          {
            std::ofstream stream(sImage_describer.c_str());
            if (!stream.is_open())
              return false;

            cereal::JSONOutputArchive archive(stream);
            archive(cereal::make_nvp("image_describer", image_describer));
            auto regionsType = image_describer->Allocate();
            archive(cereal::make_nvp("regions_type", regionsType));
          }
          // Feature extraction routines
          // For each View of the SfM_Data container:
          // - if regions file exists continue,
          // - if no file, compute features

          regions_provider = std::make_shared<Regions_Provider>();




          {
              system::Timer timer;
              Image<unsigned char> globalMask;
              Image<unsigned char> imageGray;

              const std::string sGlobalMask_filename = stlplus::create_filespec(sOutDir, "mask.png");
              if (stlplus::file_exists(sGlobalMask_filename))
              {
                if (ReadImage(sGlobalMask_filename.c_str(), &globalMask))
                {
                  std::cout
                    << "Feature extraction will use a GLOBAL MASK:\n"
                    << sGlobalMask_filename << std::endl;
                }
              }

            C_Progress_display my_progress_bar( sfm_data.GetViews().size(),
              std::cout, "\n- EXTRACT FEATURES -\n" );

            if (!stlplus::folder_exists(sfm_data.s_root_path+"/Dense"))
            {
              if (!stlplus::folder_create(sfm_data.s_root_path+"/Dense"))
              {
                std::cerr << "Cannot create Dense directory" << std::endl;
                return EXIT_FAILURE;
              }
            }

            #ifdef OPENMVG_USE_OPENMP
            const unsigned int nb_max_thread = omp_get_max_threads();

            if (iNumThreads > 0) {
                omp_set_num_threads(iNumThreads);
            } else {
                omp_set_num_threads(nb_max_thread);
            }

//            std::ofstream out;
//            out.open(sfm_data.s_root_path+"/feat_record.txt");

            #pragma omp parallel for if (iNumThreads > 0) private(imageGray)
            #endif
            for (int i = 0; i < static_cast<int>(sfm_data.views.size()); ++i)
            {
              Views::const_iterator iterViews = sfm_data.views.begin();
              std::advance(iterViews, i);
              const View * view = iterViews->second.get();
              const std::string
                sView_filename = stlplus::create_filespec(sfm_data.s_root_path, view->s_Img_path),////!!!!
                sFeat = stlplus::create_filespec(sOutDir, stlplus::basename_part(sView_filename), "feat"),
                sDesc = stlplus::create_filespec(sOutDir, stlplus::basename_part(sView_filename), "desc");

//              #pragma omp critical
//              {
//                    out << view->s_Img_path << " in error1" << std::endl;
//              }

              //If features or descriptors file are missing, compute them
        //      if (!stlplus::file_exists(sFeat) || !stlplus::file_exists(sDesc))
              {

                std::string
                    out_filename = stlplus::create_filespec(sfm_data.s_root_path+"/Dense", view->s_Img_path);////!!!!
                if (!stlplus::folder_exists(sfm_data.s_root_path+"/Dense"))
                    stlplus::folder_create(sfm_data.s_root_path+"/Dense");
                cv::Mat image_in, uimage, image_resize, image_gray;


                if(sNeedundistort != "")
                {
                    const std::string sView_undisFileName = stlplus::create_filespec(sfm_data.s_root_path + "/Dense", view->s_Img_path);

                    const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
                    // undistort image and save it
                    Image<openMVG::image::RGBColor> imageRGB, imageRGB_ud;
                    ReadImage(sView_filename.c_str(), &imageRGB);
                    cameras::UndistortImage(imageRGB, cam, imageRGB_ud, view->id_intrinsic, map_d_list);
//                    #pragma omp critical
                    WriteImage(sView_undisFileName.c_str(), imageRGB_ud);
                    //
                    image_in.create(cvSize(imageRGB_ud.Width(), imageRGB_ud.Height()), CV_8UC3);
                    //
//                    {
//                        //save for test
//                        const std::string sView_undisFileName = stlplus::create_filespec(sfm_data.s_root_path+"/testIn1/", view->s_Img_path);
////                        WriteImage(sView_undisFileName.c_str(), image_in);
//                        cv::imwrite(sView_undisFileName.c_str(), image_in);
//                    }
                    //convert RGBA to gray
                    Image<unsigned char> image_in_transTemp;
//                    ConvertPixelType( imageRGB_ud, &image_in_transTemp);
                    #pragma omp parallel for
                    for(int i = 0; i < image_in_transTemp.Height(); i++)
                    {
                        #pragma omp parallel for
                        for(int j = 0; j < image_in_transTemp.Width(); j++)
                        {
//                            image_in.at<unsigned char>(i, j) = image_in_transTemp(i, j);
                            image_in.at<unsigned char>(i, j) = imageRGB_ud(i, j);
//                            imageGray(i, j)=image_gray.at<unsigned char>(i, j);
                        }
                    }
                    //
//                    image_in=cv::imread(sView_undisFileName.c_str());
                    //
                }else{
//                    #pragma omp critical
                    {
//                        out << view->s_Img_path << " read error2" << std::endl;
//                        std::cout << view->s_Img_path << " read error2" << std::endl;
                        image_in = cv::imread(sView_filename.c_str());
//                        out <<  view->s_Img_path << " read finished error3"  << std::endl;
//                        std::cout <<  view->s_Img_path << " read finished error3"  << std::endl;
                    }
                    //
                }
                if(image_in.empty())
                {
//                    out << "error0 can not load imag : "<<sView_filename<<std::endl;
                    std::cout<<"can not load imag : "<<sView_filename<<std::endl;
                    continue;
                }
                uimage.create(image_in.size(), image_in.type());
                uimage = image_in.clone();
                if(bNeedOneSixteenth)
                {
                    cv::resize(image_in, image_resize, cv::Size(0,0), 0.25, 0.25);
                }else{
                    image_resize = image_in.clone();
                }
                cv::cvtColor(image_resize, image_gray, CV_BGR2GRAY);

                imageGray.resize(image_gray.cols, image_gray.rows);
//                std::cout << image_gray.cols << " " << image_gray.rows << " " << imageGray.Height() << " " << imageGray.Width() << std::endl;
                #pragma omp parallel for
                for(int i = 0; i < imageGray.Height(); i++)
                    #pragma omp parallel for
                    for(int j = 0; j < imageGray.Width(); j++)
                    {
                        imageGray(i, j)=image_gray.at<unsigned char>(i, j);
                    }

//                #pragma omp critical
//                {
////                    std::cout << view->s_Img_path << " imageGray finished error4" << std::endl;
//                    out << view->s_Img_path << " imageGray finished error4" << std::endl;
//                }
        //        if (!ReadImage(sView_filename.c_str(), &imageGray))
        //          continue;

                Image<unsigned char> * mask = nullptr; // The mask is null by default
                const std::string sImageMask_filename =
                  stlplus::create_filespec(sfm_data.s_root_path,
                    stlplus::basename_part(sView_filename) + "_mask", "png");
                Image<unsigned char> imageMask;
//                if (stlplus::file_exists(sImageMask_filename))
//                  ReadImage(sImageMask_filename.c_str(), &imageMask);

                /// preparing mask
                ///   -needgenerateMask
                ///
                if(bNeedGenerateMask)
                {
                    /////////////////////////////////////
                    cv::Mat _imageMask = cv::Mat(image_gray.size(), CV_8UC1, cv::Scalar::all(0));

//                    std::cout << image_in.cols << " " << image_in.rows << std::endl;
//                    std::cout << _imageMask.cols << " " << _imageMask.rows << std::endl;

                    cv::threshold(image_gray, _imageMask, 240, 255, CV_THRESH_BINARY);
//                    std::cout << _imageMask.cols << " " << _imageMask.rows << std::endl;


                    cv::Mat Kernel = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(9, 9));
                    cv::dilate(_imageMask, _imageMask, Kernel);
//                    std::cout << _imageMask.cols << " " << _imageMask.rows << std::endl;

                    cv::threshold(_imageMask, _imageMask, 100, 255, CV_THRESH_BINARY);
                    cv::bitwise_not(_imageMask, _imageMask);
                    //cv::resize(_imageMask, _imageMask, cv::Size(0,0), 0.25, 0.25); ///!
//                    #pragma omp critical
//                    cv::imwrite(sImageMask_filename, _imageMask);

                    imageMask.resize(_imageMask.cols, _imageMask.rows);
//                    std::cout << _imageMask.cols << " " << _imageMask.rows << std::endl;
//                    #pragma omp critical
//                    {
////                        std::cout << view->s_Img_path << " imageMask finished error5" << std::endl;
//                        out << view->s_Img_path << " imageMask finished error5" << std::endl;
//                    }
                    //
                    #pragma omp parallel for
                    for(int i=0; i<imageMask.Height(); i++)
                        #pragma omp parallel for
                        for(int j=0; j<imageMask.Width(); j++)
                        {
                            imageMask(i, j) = _imageMask.at<unsigned char>(i, j);
                        }
                }

//                // The mask point to the globalMask, if a valid one exists for the current image
//                if (globalMask.Width() == imageGray.Width() && globalMask.Height() == imageGray.Height())
//                  mask = &globalMask;
//                // The mask point to the imageMask (individual mask) if a valid one exists for the current image
//                if (imageMask.Width() == imageGray.Width() && imageMask.Height() == imageGray.Height())
                  mask = &imageMask;

//                {
//                    const std::string sImageMask_filename =
//                      stlplus::create_filespec(sfm_data.s_root_path,
//                        stlplus::basename_part(sView_filename) + "_testmask", "png");
//                    WriteImage(sImageMask_filename.c_str(), imageGray);
//                }
                // Compute features and descriptors and export them to files
                auto regions = image_describer->Describe(imageGray, mask);

//                #pragma omp critical
//                auto regions = image_describer->Describe(imageGray, mask);
//                {
////                    auto regions = image_describer->Describe(imageGray, mask);
////                    std::cout << view->s_Img_path << " Describe finished error6" << std::endl;
//                    out << view->s_Img_path << " Describe finished error6" << std::endl;
//                }


                /// record for the feature less then threshold
                ///   -iFeatureLimit
                ///
                if(regions->RegionCount() < iFeatureLimit)
                {
//                    std::cout << "error 1 1" << std::endl;
                    #pragma omp critical
                    {
                        reFeatureIdList.push_back(i);
//                        std::cout << view->s_Img_path << " record finished error7.1" << std::endl;
//                        out << view->s_Img_path << " record finished error7.1" << std::endl;
                    }
                }else
                {
//                    #pragma omp critical
                    {
                        image_describer->Save(regions.get(), sFeat, sDesc);
                        regions_provider->set(regions, iterViews->first);
//                        std::cout << view->s_Img_path << " saveset finished error7.2" << std::endl;
//                        out << view->s_Img_path << " saveset finished error7.2" << std::endl;
                    }
                }

                }
                ///for matching step
        //        #pragma omp critical
                {
        //            regions_provider->set(regions, iterViews->first);
                }


              ++my_progress_bar;
            }
//            out.close();

            std::cout << "Task done in (s): " << timer.elapsed() << std::endl;
          }


          std::cout << "reFeatureIdList " << reFeatureIdList.size() << std::endl;

      }

      /// compute again for features less than threshold
      ///
      {
          Image<unsigned char> globalMask;
          const std::string sGlobalMask_filename = stlplus::create_filespec(sOutDir, "mask.png");
          if (stlplus::file_exists(sGlobalMask_filename))
          {
            if (ReadImage(sGlobalMask_filename.c_str(), &globalMask))
            {
              std::cout
                << "Feature extraction will use a GLOBAL MASK:\n"
                << sGlobalMask_filename << std::endl;
            }
          }
          C_Progress_display my_progress_bar( reFeatureIdList.size(),
            std::cout, "\n- EXTRACT FEATURES -\n" );

          //
          std::unique_ptr<Image_describer> image_describer_new;
          SIFT_Image_describer::Params params_new;
          params_new._peak_threshold /= 5.0;
          image_describer_new.reset(new SIFT_Image_describer(params_new, !bUpRight));

          #pragma omp parallel for schedule(dynamic)
          for (int i = 0; i < reFeatureIdList.size(); ++i)
          {

            Views::const_iterator iterViews = sfm_data.views.begin();
            std::advance(iterViews, reFeatureIdList[i]);
            const View * view = iterViews->second.get();
            const std::string
              sView_filename = stlplus::create_filespec(sfm_data.s_root_path, view->s_Img_path),////!!!!
              sFeat = stlplus::create_filespec(sOutDir, stlplus::basename_part(sView_filename), "feat"),
              sDesc = stlplus::create_filespec(sOutDir, stlplus::basename_part(sView_filename), "desc");

            {
              std::string
                  out_filename = stlplus::create_filespec(sfm_data.s_root_path+"/Dense", view->s_Img_path);////!!!!
              if (!stlplus::folder_exists(sfm_data.s_root_path+"/Dense"))
                  stlplus::folder_create(sfm_data.s_root_path+"/Dense");
              cv::Mat image_in, uimage, image_resize, image_gray;

              if(sNeedundistort!="")
              {
                  const std::string sView_undisFileName = stlplus::create_filespec(sfm_data.s_root_path + "/Dense", view->s_Img_path);
                  //
                  image_in=cv::imread(sView_undisFileName.c_str());
                  //
              }else{
                  image_in=cv::imread(sView_filename.c_str());
                  //
              }
              if(image_in.empty())
              {
                  std::cout<<"can not load imag :"<<sView_filename<<std::endl;
                  continue;
              }
              uimage.create(image_in.size(), image_in.type());
              uimage = image_in.clone();
              if(bNeedOneSixteenth)
              {
                  cv::resize(image_in, image_resize, cv::Size(0,0), 0.25, 0.25);
              }else{
                  image_resize = image_in.clone();
              }
              cv::cvtColor(image_resize, image_gray, CV_BGR2GRAY);

              Image<unsigned char> imageGray_new; //= imageGray;
              imageGray_new.resize(image_gray.cols, image_gray.rows);
              #pragma omp parallel for
              for(int i = 0; i < imageGray_new.Height(); i++)
                  #pragma omp parallel for
                  for(int j = 0; j < imageGray_new.Width(); j++)
                  {
                      imageGray_new(i, j)=image_gray.at<unsigned char>(i, j);
        //              imageGray(i, j)=gray.at<unsigned char>(i, j);
                  }

              Image<unsigned char> * mask_new = nullptr; // The mask is null by default
              const std::string sImageMask_filename =
                stlplus::create_filespec(sfm_data.s_root_path,
                  stlplus::basename_part(sView_filename) + "_mask", "png");
              Image<unsigned char> imageMask_new;
//              if (stlplus::file_exists(sImageMask_filename))
//                ReadImage(sImageMask_filename.c_str(), &imageMask_new);

              /// preparing mask
              ///   -needgenerateMask
              ///
              if(bNeedGenerateMask)
              {
                  /////////////////////////////////////
                  cv::Mat _imageMask = cv::Mat(image_gray.size(), CV_8UC1, cv::Scalar::all(0));

//                  std::cout << image_in.cols << " " << image_in.rows << std::endl;
//                  std::cout << _imageMask.cols << " " << _imageMask.rows << std::endl;

                  cv::threshold(image_gray, _imageMask, 240, 255, CV_THRESH_BINARY);
//                  std::cout << _imageMask.cols << " " << _imageMask.rows << std::endl;


                  cv::Mat Kernel = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(9, 9));
                  cv::dilate(_imageMask, _imageMask, Kernel);
//                  std::cout << _imageMask.cols << " " << _imageMask.rows << std::endl;

                  cv::threshold(_imageMask, _imageMask, 100, 255, CV_THRESH_BINARY);
                  cv::bitwise_not(_imageMask, _imageMask);
                  //cv::resize(_imageMask, _imageMask, cv::Size(0,0), 0.25, 0.25); ///!
                  cv::imwrite(sImageMask_filename, _imageMask);

                  imageMask_new.resize(_imageMask.cols, _imageMask.rows);
//                  std::cout << _imageMask.cols << " " << _imageMask.rows << std::endl;
                  //
                  #pragma omp parallel for
                  for(int i=0; i<imageMask_new.Height(); i++)
                      #pragma omp parallel for
                      for(int j=0; j<imageMask_new.Width(); j++)
                      {
                          imageMask_new(i, j) = _imageMask.at<unsigned char>(i, j);
                      }
              }

              // The mask point to the globalMask, if a valid one exists for the current image
//              if (globalMask.Width() == imageGray_new.Width() && globalMask.Height() == imageGray_new.Height())
//                mask_new = &globalMask;
              // The mask point to the imageMask (individual mask) if a valid one exists for the current image
//              if (imageMask_new.Width() == imageGray_new.Width() && imageMask_new.Height() == imageGray_new.Height())
//                mask_new = &imageMask_new;

              // Compute features and descriptors and export them to files
              auto regions_new = image_describer_new->Describe(imageGray_new, mask_new);
//              #pragma omp critical
              {
                  image_describer_new->Save(regions_new.get(), sFeat, sDesc);
                  regions_provider->set(regions_new, iterViews->first);
              }


            }
            ++my_progress_bar;
          }
      }
  }

  ///release
  ///
  {
//      map_x[0].release();
//      map_x[1].release();
//      map_x[2].release();
//      //delete [] map_x;
//      map_y[0].release();
//      map_y[1].release();
//      map_y[2].release();
//     // delete [] map_y;
  }

  // Create sift finished flag
  std::string finishedFile = sOutDir+ "/sift_finished";
  if (!stlplus::folder_exists(finishedFile))
  {
    if (!stlplus::folder_create(finishedFile))
    {
      std::cerr << "Cannot create finishedFile directory" << std::endl;
      return EXIT_FAILURE;
    }
  }

  /// Compute Matches ***
  ///   --need regions_provider holonomic
  ///   --map_PutativesMatches result
  ///
  ///
  /// ******************
  matching::PairWiseMatches map_PutativesMatches;
  {
      std::cout << "Compute Matches" << std::endl;
      //
      if(bNeedOneSixteenth)
      {
          // change Intrinsics
          for(Intrinsics::iterator itI = sfm_data.intrinsics.begin();
              itI != sfm_data.intrinsics.end(); ++itI)
          {
              std::cout << "intrinsic id : " << itI->first << std::endl;
              std::vector<double> thisDatum = itI->second->getParams();
              std::shared_ptr<cameras::IntrinsicBase> tempIntr(NULL);
              tempIntr = std::make_shared<cameras::Pinhole_Intrinsic>
                      (itI->second->w()/4, itI->second->h()/4, thisDatum[0]/4.0, thisDatum[1]/4.0, thisDatum[2]/4.0);
               sfm_data.intrinsics[itI->first] = tempIntr;
          }
      }
      //
      {
          std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);

          EPairMode ePairmode = PAIR_EXHAUSTIVE;
          if (sPredefinedPairList.length()) {
            ePairmode = PAIR_FROM_FILE;
          }
          sMatchesDirectory = sOutDir;
          if (sMatchesDirectory.empty() || !stlplus::is_folder(sMatchesDirectory))  {
            std::cerr << "\nIt is an invalid output directory" << std::endl;
            return EXIT_FAILURE;
          }
          //
          EGeometricModel eGeometricModelToCompute;
          eGeometricModelToCompute = ESSENTIAL_MATRIX;
          std::string sGeometricMatchesFilename = "";
          sGeometricMatchesFilename = "matches.e.bin";
          //
          // -----------------------------
          // - Load SfM_Data Views & intrinsics data
          // a. Compute putative descriptor matches
          // b. Geometric filtering of putative matches
          // + Export some statistics
          // -----------------------------
          //
          //---------------------------------------
          // Read SfM Scene (image view & intrinsics data)
          //---------------------------------------
          //...
          //
          //---------------------------------------
          // Load SfM Scene regions
          //---------------------------------------
          // Init the regions_type from the image describer file (used for image regions extraction)

          if (!regions_type)
          {
            std::cerr << "Invalid: "
              << sImage_describer << " regions type file." << std::endl;
            return EXIT_FAILURE;
          }
          std::cout << "error 2 " << std::endl;
          //
          //---------------------------------------
          // a. Compute putative descriptor matches
          //    - Descriptor matching (according user method choice)
          //    - Keep correspondences only if NearestNeighbor ratio is ok
          //---------------------------------------
          //..
          // Show the progress on the command line:
          C_Progress_display progress;
          //
          // Build some alias from SfM_Data Views data:
          // - List views as a vector of filenames & image sizes
          std::vector<std::string> vec_fileNames;
          std::vector<std::pair<size_t, size_t> > vec_imagesSize;
          {
            vec_fileNames.reserve(sfm_data.GetViews().size());
            vec_imagesSize.reserve(sfm_data.GetViews().size());
            for (Views::const_iterator iter = sfm_data.GetViews().begin();
              iter != sfm_data.GetViews().end();
              ++iter)
            {
              const View * v = iter->second.get();
              vec_fileNames.push_back(stlplus::create_filespec(sfm_data.s_root_path,
                  v->s_Img_path));
              vec_imagesSize.push_back( std::make_pair( v->ui_width, v->ui_height) );
            }
          }

          std::cout << std::endl << " - PUTATIVE MATCHES - " << std::endl;
          std::cout << "Using ANN_L2 matcher" << std::endl;
          std::unique_ptr<matching_image_collection::Matcher> collectionMatcher;
          collectionMatcher.reset(new matching_image_collection::Matcher_Regions(fDistRatio, matching::ANN_L2));
          // Perform the matching
          system::Timer timer;
          {
            // From matching mode compute the pair list that have to be matched:
            Pair_Set pairs;
            if (!loadPairs(sfm_data.GetViews().size(), sPredefinedPairList, pairs))
            {
                std::cout << "loadPairs failed!" << std::endl;
                return EXIT_FAILURE;
            }
            // Photometric matching of putative pairs
            collectionMatcher->Match(sfm_data, regions_provider, pairs, map_PutativesMatches, &progress);
            //---------------------------------------
            //-- Export putative matches
            //---------------------------------------
            //..
          }
          std::cout << "Task (Regions Matching) done in (s): " << timer.elapsed() << std::endl;
          //
          //-- export putative matches Adjacency matrix
          PairWiseMatchingToAdjacencyMatrixSVG(vec_fileNames.size(),
            map_PutativesMatches,
            stlplus::create_filespec(sMatchesDirectory, "PutativeAdjacencyMatrix", "svg"));
          //-- export view pair graph once putative graph matches have been computed
          {
//            std::set<IndexT> set_ViewIds;
//            std::transform(sfm_data.GetViews().begin(), sfm_data.GetViews().end(),
//              std::inserter(set_ViewIds, set_ViewIds.begin()), stl::RetrieveKey());
//            graph::indexedGraph putativeGraph(set_ViewIds, getPairs(map_PutativesMatches));
//            graph::exportToGraphvizData(
//              stlplus::create_filespec(sMatchesDirectory, "putative_matches"),
//              putativeGraph);
          }
          //
          //---------------------------------------
          // b. Geometric filtering of putative matches
          //    - AContrario Estimation of the desired geometric model
          //    - Use an upper bound for the a contrario estimated threshold
          //---------------------------------------
          int imax_iteration = 2048;
//          //
          std::unique_ptr<matching_image_collection::ImageCollectionGeometricFilter> filter_ptr(
            new matching_image_collection::ImageCollectionGeometricFilter(&sfm_data, regions_provider));
          //
          if (filter_ptr)
          {
            system::Timer timer;
            const double d_distance_ratio = 0.6;
            bool bGuided_matching = false;

            matching::PairWiseMatches map_GeometricMatches;
            filter_ptr->Robust_model_estimation(matching_image_collection::GeometricFilter_EMatrix_AC(4.0, imax_iteration),
              map_PutativesMatches, bGuided_matching, d_distance_ratio, &progress);
            map_GeometricMatches = filter_ptr->Get_geometric_matches();

            //-- Perform an additional check to remove pairs with poor overlap
            std::vector<matching::PairWiseMatches::key_type> vec_toRemove;
            for (const auto & pairwisematches_it : map_GeometricMatches)
            {
              const size_t putativePhotometricCount = map_PutativesMatches.find(pairwisematches_it.first)->second.size();
              const size_t putativeGeometricCount = pairwisematches_it.second.size();
              const float ratio = putativeGeometricCount / static_cast<float>(putativePhotometricCount);
              if (putativeGeometricCount < 1){// || ratio < .3f)  {
                // the pair will be removed
                vec_toRemove.push_back(pairwisematches_it.first);
              }
            }
            //-- remove discarded pairs
            for (const auto & pair_to_remove_it : vec_toRemove)
            {
              map_GeometricMatches.erase(pair_to_remove_it);
            }
            //
            //---------------------------------------
            //-- Export geometric filtered matches
            //---------------------------------------
            //..
//            map_GeometricMatches = map_PutativesMatches;//!!!
            if (!Save(map_GeometricMatches,
              std::string(sMatchesDirectory + "/" + sGeometricMatchesFilename)))
            {
              std::cerr
                  << "Cannot save computed matches in: "
                  << std::string(sMatchesDirectory + "/" + sGeometricMatchesFilename);
              return EXIT_FAILURE;
            }
            std::cout << "Task done in (s): " << timer.elapsed() << std::endl;
            //-- export Adjacency matrix
            std::cout << "\n Export Adjacency Matrix of the pairwise's geometric matches"
              << std::endl;
            matching::PairWiseMatchingToAdjacencyMatrixSVG(vec_fileNames.size(),
              map_GeometricMatches,
              stlplus::create_filespec(sMatchesDirectory, "GeometricAdjacencyMatrix", "svg"));

//            std::cout << "error 2 4 " << std::endl;

            //-- export view pair graph once geometric filter have been done
            {
//              std::set<IndexT> set_ViewIds;
//              std::transform(sfm_data.GetViews().begin(), sfm_data.GetViews().end(),
//                std::inserter(set_ViewIds, set_ViewIds.begin()), stl::RetrieveKey());
//              graph::indexedGraph putativeGraph(set_ViewIds, getPairs(map_GeometricMatches));
//              graph::exportToGraphvizData(
//                stlplus::create_filespec(sMatchesDirectory, "geometric_matches"),
//                putativeGraph);
            }
//            std::cout << "error 2 5 " << std::endl;
          }
      }
  }

  /// preparing for global step
  ///
//  std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
//  for (auto itView = sfm_data.views.begin(); itView != sfm_data.views.end(); ++ itView)
//  {
////      get(const IndexT x)
//      const IndexT imgId = itView->first;
//      feats_provider->feats_per_view[imgId] = regions_provider->get(imgId)->GetRegionsPositions();
//      #pragma omp parallel for schedule(dynamic)
//      for(size_t i = 0; i < feats_provider->feats_per_view[imgId].size(); ++i)
//      {
//          PointFeatures::iterator itPF = feats_provider->feats_per_view[imgId].begin();
//          std::advance(itPF, i);
//          (*itPF).x() *= 4.0;
//          (*itPF).y() *= 4.0;
//      }
//  }

//  //regions_provider.reset();

//  ///Compute global ***
//  ///
//  ///
//  ///
//  /// ******************
//  std::cout << "compute global" << std::endl;
//  {
//      const cameras::Intrinsic_Parameter_Type intrinsic_refinement_options =
//        cameras::StringTo_Intrinsic_Parameter_Type(sIntrinsic_refinement_options);
//      if (intrinsic_refinement_options == static_cast<cameras::Intrinsic_Parameter_Type>(0) )
//      {
//        std::cerr << "Invalid input for Bundle Adjusment Intrinsic parameter refinement option" << std::endl;
//        return EXIT_FAILURE;
//      }
//      std::cout << "error 3 2 " << std::endl;
//      // Matches reading
//      std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
//      matches_provider->pairWise_matches_ = map_PutativesMatches;

//      std::cout << "error 3 3 " << std::endl;
//      openMVG::system::Timer timer;
//      GlobalSfMReconstructionEngine_RelativeMotions sfmEngine(
//        sfm_data,
//        sOutDir,
//        stlplus::create_filespec(sOutDir, "Reconstruction_Report.html"));

//      std::cout << "error 3 4 " << std::endl;
//      // Configure the features_provider & the matches_provider
//      sfmEngine.SetFeaturesProvider(feats_provider.get());
//      sfmEngine.SetMatchesProvider(matches_provider.get());

//      // Configure reconstruction parameters
//      sfmEngine.Set_Intrinsics_Refinement_Type(intrinsic_refinement_options);
//      sfmEngine.Set_Use_Motion_Prior(false);

//      // Configure motion averaging method
//      sfmEngine.SetRotationAveragingMethod(
//        ERotationAveragingMethod(iRotationAveragingMethod));
//      sfmEngine.SetTranslationAveragingMethod(
//        ETranslationAveragingMethod(iTranslationAveragingMethod));

//      std::cout << "error 3 5 " << std::endl;
//      if (sfmEngine.Process_GCP_GPS(rcFilePath))//threeCameras_change_step(rtFilePath, refSfmDataIdList))//(rtFilePath, refSfmDataIdList))//Process_threeCameras_imu_gps_init_new())//Process_test(sMatchesDir))//_down_imu_gps(sMatchesDir))
//      {
//          std::cout << std::endl << " Total Ac-Global-Sfm took (s): " << timer.elapsed() << std::endl;

//          Save(
//              sfmEngine.Get_SfM_Data(),
//              stlplus::create_filespec( sfm_data.s_root_path, "sfm_data_extrinsics.json" ).c_str(),
//              ESfM_Data(VIEWS|INTRINSICS|EXTRINSICS));

//          std::cout << "...Generating SfM_Report.html" << std::endl;
//          Generate_SfM_Report(sfmEngine.Get_SfM_Data(),
//            stlplus::create_filespec(sOutDir, "SfMReconstruction_Report.html"));

//          //-- Export to disk computed scene (data & visualizable results)
//          std::cout << "...Export SfM_Data to disk." << std::endl;
//          std::cout << "error 3 1" << std::endl;
//          Save(sfmEngine.Get_SfM_Data(),
//            stlplus::create_filespec(sOutDir, "sfm_data", ".bin"),
//            ESfM_Data(ALL));

//          std::cout << "error 3 2" << std::endl;
//          Save(sfmEngine.Get_SfM_Data(),
//            stlplus::create_filespec(sOutDir, "cloud_and_poses", ".ply"),
//            ESfM_Data(ALL));

//          std::cout << "error 3 3" << std::endl;
//          return EXIT_SUCCESS;
//          std::cout << "error 3 4" << std::endl;

//      }else{
//          std::cout << "error 3 5" << std::endl;
//          std::cout << "failed when doing global-step!" << std::endl;
//          return EXIT_FAILURE;

//      }

//      std::cout << "error 3 6" << std::endl;




//  }

  return EXIT_SUCCESS;
}
