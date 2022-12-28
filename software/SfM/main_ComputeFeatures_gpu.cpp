
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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

#include <cstdlib>
#include <fstream>
#include <string>

//#include "openMVG/image/image.hpp"
//#include "openMVG/sfm/sfm.hpp"

///// Feature/Regions & Image describer interfaces
//#include "openMVG/features/features.hpp"
//#include "nonFree/sift/SIFT_describer.hpp"
//#include <cereal/archives/json.hpp>
//#include "openMVG/system/timer.hpp"
//#include "openMVG/features/feature.hpp"
//#include "openMVG/features/descriptor.hpp"

//#include "third_party/cmdLine/cmdLine.h"
//#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
//#include "third_party/progress/progress.hpp"

//#include <cstdlib>
//#include <fstream>

#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif


//for match
//#include "openMVG/sfm/sfm_data.hpp"
//#include "openMVG/sfm/sfm_data_io.hpp"
//#include "openMVG/sfm/pipelines/sfm_engine.hpp"
//#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
//#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"

/// Generic Image Collection image matching
//#include "openMVG/matching_image_collection/Matcher_Regions_AllInMemory.hpp"
//#include "openMVG/matching_image_collection/Cascade_Hashing_Matcher_Regions_AllInMemory.hpp"
//#include "openMVG/matching_image_collection/GeometricFilter.hpp"
//#include "openMVG/matching_image_collection/F_ACRobust.hpp"
//#include "openMVG/matching_image_collection/E_ACRobust.hpp"
//#include "openMVG/matching_image_collection/H_ACRobust.hpp"
//#include "openMVG/matching/pairwiseAdjacencyDisplay.hpp"
//#include "openMVG/matching/indMatch_utils.hpp"


//#include "openMVG/graph/graph.hpp"
//#include "openMVG/stl/stl.hpp"


//#include
//#include <third_party/Sift_GPU/SiftGPU.h>

#include "third_party/SIFT_GPU/SiftGPU/SiftGPU.h"
#include <iostream>
#include <vector>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <boost/timer.hpp>

#include <GL/gl.h>





using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::features;
using namespace openMVG::sfm;
using namespace std;

features::EDESCRIBER_PRESET stringToEnum(const std::string & sPreset)
{
  features::EDESCRIBER_PRESET preset;
  if(sPreset == "NORMAL")
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
#ifdef OPENMVG_USE_OPENMP
  int iNumThreads = 0;
#endif

  // required
  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  // Optional
  cmd.add( make_option('m', sImage_Describer_Method, "describerMethod") );
  cmd.add( make_option('u', bUpRight, "upright") );
  cmd.add( make_option('f', bForce, "force") );
  cmd.add( make_option('p', sFeaturePreset, "describerPreset") );

#ifdef OPENMVG_USE_OPENMP
  cmd.add( make_option('n', iNumThreads, "numThreads") );
#endif

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--input_file] a SfM_Data file \n"
      << "[-o|--outdir path] \n"
      << "\n[Optional]\n"
      << "[-f|--force] Force to recompute data\n"
      << "[-m|--describerMethod]\n"
      << "  (method to use to describe an image):\n"
      << "   SIFT (default),\n"
      << "   AKAZE_FLOAT: AKAZE with floating point descriptors,\n"
      << "   AKAZE_MLDB:  AKAZE with binary descriptors\n"
      << "[-u|--upright] Use Upright feature 0 or 1\n"
      << "[-p|--describerPreset]\n"
      << "  (used to control the Image_describer configuration):\n"
      << "   NORMAL (default),\n"
      << "   HIGH,\n"
      << "   ULTRA: !!Can take long time!!\n"
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

  //---------------------------------------
  // a. Load input scene
  //---------------------------------------
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS))) {
    std::cerr << std::endl
      << "The input file \""<< sSfM_Data_Filename << "\" cannot be read" << std::endl;
    return false;
  }

  // b. Init the image_describer
  // - retrieve the used one in case of pre-computed features
  // - else create the desired one

  using namespace openMVG::features;
  std::unique_ptr<Image_describer> image_describer;

  const std::string sImage_describer = stlplus::create_filespec(sOutDir, "image_describer", "json");
  if (!bForce && stlplus::is_file(sImage_describer))
  {
    // Dynamically load the image_describer from the file (will restore old used settings)
    std::ifstream stream(sImage_describer.c_str());
    if (!stream.is_open())
      return false;

    try
    {
      cereal::JSONInputArchive archive(stream);
      archive(cereal::make_nvp("image_describer", image_describer));
    }
    catch (const cereal::Exception & e)
    {
      std::cerr << e.what() << std::endl
        << "Cannot dynamically allocate the Image_describer interface." << std::endl;
      return EXIT_FAILURE;
    }
  }
  else
  {
      // Create the desired Image_describer method.
      // Don't use a factory, perform direct allocation
      if (sImage_Describer_Method == "SIFT")
      {
        image_describer.reset(new SIFT_Image_describer
          (SIFT_Image_describer::Params(), !bUpRight));
      }
      else
      if (sImage_Describer_Method == "SIFT_ANATOMY")
      {
        image_describer.reset(
          new SIFT_Anatomy_Image_describer(SIFT_Anatomy_Image_describer::Params()));
      }
      else
      if (sImage_Describer_Method == "AKAZE_FLOAT")
      {
        image_describer = AKAZE_Image_describer::create
          (AKAZE_Image_describer::Params(AKAZE::Params(), AKAZE_MSURF), !bUpRight);
      }
      else
      if (sImage_Describer_Method == "AKAZE_MLDB")
      {
        image_describer = AKAZE_Image_describer::create
          (AKAZE_Image_describer::Params(AKAZE::Params(), AKAZE_MLDB), !bUpRight);
      }
      if (!image_describer)
      {
        std::cerr << "Cannot create the designed Image_describer:"
          << sImage_Describer_Method << "." << std::endl;
        return EXIT_FAILURE;
      }
      else
      {
        if (!sFeaturePreset.empty())
        if (!image_describer->Set_configuration_preset(stringToEnum(sFeaturePreset)))
        {
          std::cerr << "Preset configuration failed." << std::endl;
          return EXIT_FAILURE;
        }
      }

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
  }

  // Feature extraction routines
  // For each View of the SfM_Data container:
  // - if regions file exists continue,
  // - if no file, compute features

  vector<float> descriptors1;//(128*num);
  vector<float> descriptors2;//(128*num);
  int num1, num2;
  //gpu
  {
      //initial SiftGPU
      SiftGPU sift;
      char* myargv[4] ={ "-fo", "-1", "-v", "1"};
      sift.ParseParam(4, myargv);

      //check SiftGPU supported
      int support = sift.CreateContextGL();
      if ( support != SiftGPU::SIFTGPU_FULL_SUPPORTED )
      {
          cerr<<"SiftGPU is not supported!"<<endl;
          return 2;
      }


    for(int i = 0; i < sfm_data.views.size(); ++i)
    {

      Views::const_iterator iterViews = sfm_data.views.begin();
      std::advance(iterViews, i);
      const View * view = iterViews->second.get();

      //prepare image path
      std::string s_res = sfm_data.s_root_path + view->s_Img_path;
      char * sView_filename = (char*) malloc((s_res.length())*sizeof(char));
      strcpy(sView_filename, s_res.c_str());

      //read image
      cout<<"running sift " << sView_filename <<endl;
      boost::timer timer;

      sift.RunSIFT(sView_filename);//000000-100006.jpg" );"/media/add7/E/run_scene/cali_three/pic/all/1DSC00225.jpg");//"/media/add7/E/run_scene/cali_three/pic/all/1DSC00225.JPG");//

      cout<<"siftgpu::runsift() cost time= "<<timer.elapsed()<<endl;

      //get key points and sift discriptors
      int num = sift.GetFeatureNum();

      cout<<"Feature number="<<num<<endl;
      vector<float> descriptors(128*num);
      vector<SiftGPU::SiftKeypoint> keys(num);
      sift.GetFeatureVector(&keys[0], &descriptors[0]);

      //transforms and save
      std::string sFeat = stlplus::create_filespec(sOutDir, stlplus::basename_part(sView_filename), "feat");
      std::string sDesc = stlplus::create_filespec(sOutDir, stlplus::basename_part(sView_filename), "desc");

      //transform formation

      std::unique_ptr<Regions> regions;
      regions.reset( new SIFT_Regions );
      for (std::size_t points_id = 0; points_id < num; ++points_id )
      {
          //each point
          //prepare one point data for transform format
//          vl_sift_pix descr[128];
          //Descriptor<float, 128> descr;///!!!
          vl_sift_pix descr[128];
          Descriptor<unsigned char, 128> descriptor;

          for (std::size_t p_vec = 0; p_vec < 128; ++p_vec)
          {
              descr[p_vec] = descriptors[points_id*128 + p_vec];
          }

          const float sum = accumulate(descr, descr+128, 0.0f);
          for (int k=0;k<128;++k)
            descriptor[k] = static_cast<unsigned char>(512.f*sqrt(descr[k]/sum));

          const SIOPointFeature fp(keys[points_id].x, keys[points_id].y, keys[points_id].s, keys[points_id].o);

          //prepare for save
          SIFT_Regions * regionsCasted = dynamic_cast<SIFT_Regions*>(regions.get());
          // reserve some memory for faster keypoint saving
          regionsCasted->Features().reserve(2000);
          regionsCasted->Descriptors().reserve(2000);

//          regionsCasted->Descriptors().push_back(descr);///!!!
          regionsCasted->Descriptors().push_back(descriptor);
          regionsCasted->Features().push_back(fp);

      }
      //call save function direct
      image_describer->Save(regions.get(), sFeat, sDesc);

      //delete regions
      regions.reset();

      std::stringstream ss;
      ss << i;
      std::string s;
      ss >> s;

      std::cout << "finish pic " << s << " using " << timer.elapsed() << std::endl;
      std::cout << "--------------------------------" << std::endl;


    }


//    ///--------------------------------
//    ///process for compute pair matches
//    PairWiseMatches map_PutativesMatches;
//    //load need image feature test
//    for(std::size_t i = 0; i < 1; ++i)
//    {
//        //load signal
//        std::unique_ptr<features::Regions> region_const_ptr;
//        region_const_ptr.reset( new SIFT_Regions );
//        //prepare path
//        Views::const_iterator iter = sfm_data.GetViews().begin();
//        std::advance(iter, i);
//        const std::string sImageName = stlplus::create_filespec(sfm_data.s_root_path, iter->second.get()->s_Img_path);
//        const std::string basename = stlplus::basename_part(sImageName);
//        const std::string featFile = stlplus::create_filespec(sOutDir, basename, ".feat");
//        const std::string descFile = stlplus::create_filespec(sOutDir, basename, ".desc");

//        if (!region_const_ptr->Load(featFile, descFile))
//        {
//          std::cerr << "Invalid regions files for the view: " << sImageName << std::endl;
//          //bContinue = false;
//        }
////        region_const_ptr = std::move(regions_ptr);

//        SIFT_Regions * region_const_Casted = dynamic_cast<SIFT_Regions*>(region_const_ptr.get());

//        region_const_Casted->Descriptors();
//        int num1 = region_const_Casted->Descriptors().size();
//        vector<float> descriptors_const(128*num1);
//        //all point of one image
//        for (std::size_t p_id = 0; p_id < region_const_Casted->Descriptors().size(); ++p_id)
//        {
//            Descriptor<float, 128> descriptor;
//            descriptor = region_const_Casted->Descriptors()[p_id];

//    //        const float sum = accumulate(descriptor, descriptor+128, 0.0f);

//            //float descr[128];
//            for (std::size_t p_vec = 0; p_vec < 128; ++p_vec)
//            {
//                descriptors_const[p_id*128 + p_vec] = descriptor[p_vec];
//            }

//        }

//        //------
//        //other image
//        for(std::size_t j = 1; j < 2; ++j)
//        {
//            //load and prepare
//            //load signal
//            std::unique_ptr<features::Regions> region_other_ptr;
//            region_other_ptr.reset( new SIFT_Regions );
//            //prepare path
//            Views::const_iterator iter = sfm_data.GetViews().begin();
//            std::advance(iter, j);
//            const std::string sImageName = stlplus::create_filespec(sfm_data.s_root_path, iter->second.get()->s_Img_path);
//            const std::string basename = stlplus::basename_part(sImageName);
//            const std::string featFile = stlplus::create_filespec(sOutDir, basename, ".feat");
//            const std::string descFile = stlplus::create_filespec(sOutDir, basename, ".desc");

//            if (!region_other_ptr->Load(featFile, descFile))
//            {
//              std::cerr << "Invalid regions files for the view: " << sImageName << std::endl;
//              //bContinue = false;
//            }
////            region_other_ptr = std::move(regions_ptr);

//            SIFT_Regions * region_other_Casted = dynamic_cast<SIFT_Regions*>(region_other_ptr.get());

//            region_other_Casted->Descriptors();
//            int num2 = region_other_Casted->Descriptors().size();
//            vector<float> descriptors_other(128*num2);
//            //all point of one image
//            for (std::size_t p_id = 0; p_id < region_other_Casted->Descriptors().size(); ++p_id)
//            {
//                Descriptor<float, 128> descriptor;
//                descriptor = region_other_Casted->Descriptors()[p_id];

//        //        const float sum = accumulate(descriptor, descriptor+128, 0.0f);

//                //float descr[128];
//                for (std::size_t p_vec = 0; p_vec < 128; ++p_vec)
//                {
//                    descriptors_other[p_id*128 + p_vec] = descriptor[p_vec];
//                }

//            }

//            //------
//            //compute matcher
//            {
//                boost::timer timer;

//                SiftMatchGPU matcher;// = pCreateNewSiftMatchGPU(169968);//4096);
//                matcher.VerifyContextGL(); //must call once

//                matcher.SetMaxSift(2*(num1+num2));
//                matcher.SetDescriptors(0, num1, &descriptors_const[0]); //image 1
//                matcher.SetDescriptors(1, num2, &descriptors_other[0]); //image 2

//                std::cout << "num1 : " << num1 << std::endl;
//                std::cout << "num2 : " << num2 << std::endl;

//                //match and get result.
//                int (*match_buf)[2] = new int[num1][2];
//                //use the default thresholds. Check the declaration in SiftGPU.h
//                int num_match = matcher.GetSiftMatch(num1, match_buf);
//                std::cout << num_match << " sift matches were found;\n";
//                std::cout<<"match img cost time= "<<timer.elapsed()<<std::endl;


//            }
//            //delete region_other_ptr ;




//        }
//    }

    std::cout << "=========" << std::endl;



  }












  return EXIT_SUCCESS;
}




