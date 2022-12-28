
//// Copyright (c) 2012, 2013 Pierre MOULON.

//// This Source Code Form is subject to the terms of the Mozilla Public
//// License, v. 2.0. If a copy of the MPL was not distributed with this
//// file, You can obtain one at http://mozilla.org/MPL/2.0/.


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

//#ifdef OPENMVG_USE_OPENMP
//#include <omp.h>
//#endif


////for match
//#include "openMVG/sfm/sfm_data.hpp"
//#include "openMVG/sfm/sfm_data_io.hpp"
//#include "openMVG/sfm/pipelines/sfm_engine.hpp"
//#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
//#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"

///// Generic Image Collection image matching
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



//#include <SiftGPU.h>

//#include <iostream>
//#include <vector>

//#include <opencv2/core/core.hpp>
//#include <opencv2/highgui/highgui.hpp>

//#include <boost/timer.hpp>

//#include <GL/gl.h>





//using namespace openMVG;
//using namespace openMVG::image;
//using namespace openMVG::features;
//using namespace openMVG::sfm;
//using namespace std;

//features::EDESCRIBER_PRESET stringToEnum(const std::string & sPreset)
//{
//  features::EDESCRIBER_PRESET preset;
//  if(sPreset == "NORMAL")
//    preset = features::NORMAL_PRESET;
//  else
//  if (sPreset == "HIGH")
//    preset = features::HIGH_PRESET;
//  else
//  if (sPreset == "ULTRA")
//    preset = features::ULTRA_PRESET;
//  else
//    preset = features::EDESCRIBER_PRESET(-1);
//  return preset;
//}

///// - Compute view image description (feature & descriptor extraction)
///// - Export computed data
int main(int argc, char **argv)
{}
//  CmdLine cmd;

//  std::string sSfM_Data_Filename;
//  std::string sOutDir = "";
//  bool bUpRight = false;
//  std::string sImage_Describer_Method = "SIFT";
//  bool bForce = false;
//  std::string sFeaturePreset = "";
//#ifdef OPENMVG_USE_OPENMP
//  int iNumThreads = 0;
//#endif

//  // required
//  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
//  cmd.add( make_option('o', sOutDir, "outdir") );
//  // Optional
//  cmd.add( make_option('m', sImage_Describer_Method, "describerMethod") );
//  cmd.add( make_option('u', bUpRight, "upright") );
//  cmd.add( make_option('f', bForce, "force") );
//  cmd.add( make_option('p', sFeaturePreset, "describerPreset") );

//#ifdef OPENMVG_USE_OPENMP
//  cmd.add( make_option('n', iNumThreads, "numThreads") );
//#endif

//  try {
//      if (argc == 1) throw std::string("Invalid command line parameter.");
//      cmd.process(argc, argv);
//  } catch(const std::string& s) {
//      std::cerr << "Usage: " << argv[0] << '\n'
//      << "[-i|--input_file] a SfM_Data file \n"
//      << "[-o|--outdir path] \n"
//      << "\n[Optional]\n"
//      << "[-f|--force] Force to recompute data\n"
//      << "[-m|--describerMethod]\n"
//      << "  (method to use to describe an image):\n"
//      << "   SIFT (default),\n"
//      << "   AKAZE_FLOAT: AKAZE with floating point descriptors,\n"
//      << "   AKAZE_MLDB:  AKAZE with binary descriptors\n"
//      << "[-u|--upright] Use Upright feature 0 or 1\n"
//      << "[-p|--describerPreset]\n"
//      << "  (used to control the Image_describer configuration):\n"
//      << "   NORMAL (default),\n"
//      << "   HIGH,\n"
//      << "   ULTRA: !!Can take long time!!\n"
//#ifdef OPENMVG_USE_OPENMP
//      << "[-n|--numThreads] number of parallel computations\n"
//#endif
//      << std::endl;

//      std::cerr << s << std::endl;
//      return EXIT_FAILURE;
//  }

//  std::cout << " You called : " <<std::endl
//            << argv[0] << std::endl
//            << "--input_file " << sSfM_Data_Filename << std::endl
//            << "--outdir " << sOutDir << std::endl
//            << "--describerMethod " << sImage_Describer_Method << std::endl
//            << "--upright " << bUpRight << std::endl
//            << "--describerPreset " << (sFeaturePreset.empty() ? "NORMAL" : sFeaturePreset) << std::endl
//            << "--force " << bForce << std::endl
//#ifdef OPENMVG_USE_OPENMP
//            << "--numThreads " << iNumThreads << std::endl
//#endif
//            << std::endl;


//  if (sOutDir.empty())  {
//    std::cerr << "\nIt is an invalid output directory" << std::endl;
//    return EXIT_FAILURE;
//  }

//  // Create output dir
//  if (!stlplus::folder_exists(sOutDir))
//  {
//    if (!stlplus::folder_create(sOutDir))
//    {
//      std::cerr << "Cannot create output directory" << std::endl;
//      return EXIT_FAILURE;
//    }
//  }

//  //---------------------------------------
//  // a. Load input scene
//  //---------------------------------------
//  SfM_Data sfm_data;
//  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS))) {
//    std::cerr << std::endl
//      << "The input file \""<< sSfM_Data_Filename << "\" cannot be read" << std::endl;
//    return false;
//  }

//  //---------------------------------------
//  // b. Compute matches in GPU
//  //---------------------------------------


//  // Feature extraction routines
//  // For each View of the SfM_Data container:
//  // - if regions file exists continue,
//  // - if no file, compute features

//  vector<float> descriptors1;//(128*num);//for debug
//  vector<float> descriptors2;//(128*num);//for debug
//  int num1, num2;//for debug

//  //initial SiftGPU
//  SiftGPU sift;
//  char* myargv[4] ={ "-fo", "-1", "-v", "1"};
//  sift.ParseParam(4, myargv);

//  //check SiftGPU supported
//  int support = sift.CreateContextGL();
//  if ( support != SiftGPU::SIFTGPU_FULL_SUPPORTED )
//  {
//      cerr<<"SiftGPU is not supported!"<<endl;
//      return 2;
//  }

//  //gpu
//  {

//    ///--------------------------------
//    ///process for compute pair matches
//    PairWiseMatches map_PutativesMatches;
//    //load need image feature test
//    for(std::size_t i = 0; i < sfm_data.GetViews().size() - 1; ++i)
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
//        for(std::size_t j = i + 1; j < sfm_data.GetViews().size(); ++j)
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

//                SiftMatchGPU matcher(4096);// = pCreateNewSiftMatchGPU(169968);//4096);

//                //SiftMatchGPU *matcher = new SiftMatchGPU(4096);
//                matcher.SetLanguage(SiftMatchGPU::SIFTMATCH_CUDA + 1);
//                //matcher->CreateContextGL();
//                matcher.VerifyContextGL(); //must call once

//                matcher.SetMaxSift(2*(num1+num2));
//                //descriptors_other = descriptors_const;
//                matcher.SetDescriptors(0, num1, &descriptors_const[0]); //image 1
//                matcher.SetDescriptors(1, num2, &descriptors_other[0]); //image 2

//                std::cout << "num1 : " << num1 << std::endl;
//                std::cout << "num2 : " << num2 << std::endl;

//                //match and get result.
////                int (*match_buf)[2] = new int[num1][2];
//                int match_buf[4096][2];
//                //use the default thresholds. Check the declaration in SiftGPU.h
//                int num_match = matcher.GetSiftMatch(4096, match_buf);

//                //recoding result




//                std::cout << num_match << " sift matches were found;\n";
//                std::cout<<"match img cost time= "<<timer.elapsed()<<std::endl;


////std::map< Pair, IndMatches >
//                IndMatches temp_matches;
//                //get matches id
//                for (std::size_t id_count = 0; id_count < num_match; ++id_count )
//                {
////                    std::stringstream ss0,ss1;
////                    std::string s0,s1;
////                    ss0 << match_buf[id_count][0];
////                    ss0 >> s0;
////                    ss1 << match_buf[id_count][1];
////                    ss1 >> s1;
////                    std::cout << s0 << " " << s1 << std::endl;
//                    IndMatch temp_match(match_buf[id_count][0], match_buf[id_count][1]);

//                    temp_matches.push_back(temp_match);

//                }

//                std::pair<Pair, IndMatches> temp_PairWise;
//                temp_PairWise.first = std::pair<IndexT, IndexT>(i,j);
//                temp_PairWise.second = temp_matches;

//                map_PutativesMatches.insert(temp_PairWise);



//            }


//        }
//    }

//    std::cout << "map_PutativesMatches size : " << map_PutativesMatches.size() << std::endl;

//    //about matches result
//    //---------------------------------------
//    //--c. Export putative matches
//    //---------------------------------------
//    if (!Save(map_PutativesMatches, std::string(sOutDir + "/matches.putative.bin")))
//    {
//      std::cerr
//        << "Cannot save computed matches in: "
//        << std::string(sOutDir + "/matches.putative.bin");
//      return EXIT_FAILURE;
//    }

//    //-- export putative matches Adjacency matrix
//    std::vector<std::string> vec_fileNames;
//    std::vector<std::pair<size_t, size_t> > vec_imagesSize;
//    {
//      vec_fileNames.reserve(sfm_data.GetViews().size());
//      vec_imagesSize.reserve(sfm_data.GetViews().size());
//      for (Views::const_iterator iter = sfm_data.GetViews().begin();
//        iter != sfm_data.GetViews().end();
//        ++iter)
//      {
//        const View * v = iter->second.get();
//        vec_fileNames.push_back(stlplus::create_filespec(sfm_data.s_root_path,
//            v->s_Img_path));
//        vec_imagesSize.push_back( std::make_pair( v->ui_width, v->ui_height) );
//      }
//    }

//    PairWiseMatchingToAdjacencyMatrixSVG(vec_fileNames.size(),
//      map_PutativesMatches,
//      stlplus::create_filespec(sOutDir, "PutativeAdjacencyMatrix", "svg"));
//    //-- export view pair graph once putative graph matches have been computed
//    {
//      std::set<IndexT> set_ViewIds;
//      std::transform(sfm_data.GetViews().begin(), sfm_data.GetViews().end(),
//        std::inserter(set_ViewIds, set_ViewIds.begin()), stl::RetrieveKey());
//      graph::indexedGraph putativeGraph(set_ViewIds, getPairs(map_PutativesMatches));
//      graph::exportToGraphvizData(
//        stlplus::create_filespec(sOutDir, "putative_matches"),
//        putativeGraph.g);
//    }

//    //---------------------------------------
//    // Geometric filtering of putative matches
//    //    - AContrario Estimation of the desired geometric model
//    //    - Use an upper bound for the a contrario estimated threshold
//    //---------------------------------------
//    // Load the corresponding view regions
//    const std::string sImage_describer = stlplus::create_filespec(sOutDir, "image_describer", "json");
//    std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
//    std::shared_ptr<Regions_Provider> regions_provider = std::make_shared<Regions_Provider>();
//    if (!regions_provider->load(sfm_data, sOutDir, regions_type)) {
//      std::cerr << std::endl << "Invalid regions." << std::endl;
//      return EXIT_FAILURE;
//    }
//    std::unique_ptr<matching_image_collection::ImageCollectionGeometricFilter> filter_ptr(
//      new matching_image_collection::ImageCollectionGeometricFilter(&sfm_data, regions_provider));



////    if (filter_ptr) ///????
////    {
//      system::Timer timer;
//      std::cout << std::endl << " - Geometric filtering - " << std::endl;

//      PairWiseMatches map_GeometricMatches;

//      bool bGuided_matching = false;
//      int imax_iteration = 2048;
//      filter_ptr->Robust_model_estimation(matching_image_collection::GeometricFilter_EMatrix_AC(400.0, imax_iteration),//4.0, imax_iteration),
//        map_PutativesMatches, bGuided_matching);
//      map_GeometricMatches = filter_ptr->Get_geometric_matches();

//      std::cout << "map_GeometricMatches size : " << map_GeometricMatches.size() << std::endl;

//      //-- Perform an additional check to remove pairs with poor overlap
//      std::vector<PairWiseMatches::key_type> vec_toRemove;
//      for (PairWiseMatches::const_iterator iterMap = map_GeometricMatches.begin();
//        iterMap != map_GeometricMatches.end(); ++iterMap)
//      {
//        const size_t putativePhotometricCount = map_PutativesMatches.find(iterMap->first)->second.size();
//        const size_t putativeGeometricCount = iterMap->second.size();
//        const float ratio = putativeGeometricCount / (float)putativePhotometricCount;

//        std::cout << "putativePhotometricCount   : " << putativePhotometricCount << std::endl;
//        std::cout << "putativeGeometricCount   : " << putativeGeometricCount << std::endl;
//        std::cout << "ratio   : " << ratio << std::endl;


//        if (putativeGeometricCount < 50 || ratio < .3f)
//        //if (putativeGeometricCount < 10 || ratio < .8f)
//        {
//          // the pair will be removed
//          vec_toRemove.push_back(iterMap->first);
//        }
//      }
//      //-- remove discarded pairs
//      for (std::vector<PairWiseMatches::key_type>::const_iterator
//        iter =  vec_toRemove.begin(); iter != vec_toRemove.end(); ++iter)
//      {
//        map_GeometricMatches.erase(*iter);
//      }

//      std::cout << "map_GeometricMatches size : " << map_GeometricMatches.size() << std::endl;


//      //---------------------------------------
//      //-- Export geometric filtered matches
//      //---------------------------------------
//      std::string sGeometricMatchesFilename = "matches.e.bin";
//      if (!Save(map_GeometricMatches,
//        std::string(sOutDir + "/" + sGeometricMatchesFilename)))
//      {
//        std::cerr
//            << "Cannot save computed matches in: "
//            << std::string(sOutDir + "/" + sGeometricMatchesFilename);
//        return EXIT_FAILURE;
//      }

//      std::cout << "Task done in (s): " << timer.elapsed() << std::endl;

//      //-- export Adjacency matrix
//      std::cout << "\n Export Adjacency Matrix of the pairwise's geometric matches"
//        << std::endl;

//      // Build some alias from SfM_Data Views data:
//      // - List views as a vector of filenames & image sizes
//      PairWiseMatchingToAdjacencyMatrixSVG(vec_fileNames.size(),
//        map_GeometricMatches,
//        stlplus::create_filespec(sOutDir, "GeometricAdjacencyMatrix", "svg"));

//      //-- export view pair graph once geometric filter have been done
//      {
//        std::set<IndexT> set_ViewIds;
//        std::transform(sfm_data.GetViews().begin(), sfm_data.GetViews().end(),
//          std::inserter(set_ViewIds, set_ViewIds.begin()), stl::RetrieveKey());
//        graph::indexedGraph putativeGraph(set_ViewIds, getPairs(map_GeometricMatches));
//        graph::exportToGraphvizData(
//          stlplus::create_filespec(sOutDir, "geometric_matches"),
//          putativeGraph.g);
//      }








//  }












////    /* export matches res for Bk */
////    std::vector<std::vector<int>> matches_res(sfm_data.views.size());
////    for(std::size_t i = 0; i < sfm_data.views.size(); ++i)
////    {
////        std::vector<int> mid(sfm_data.views.size());
////        for(std::size_t j = 0; j < sfm_data.views.size(); ++j)
////        {
////            mid[j] = 0;
////        }
////        matches_res[i] = mid;
////    }
////    //  for(PairWiseMatches::iterator it = map_PutativesMatches.begin(); it != map_PutativesMatches.end(); ++it)
////      for(PairWiseMatches::iterator it = map_GeometricMatches.begin(); it != map_GeometricMatches.end(); ++it)
////      {
////          it->first.second;
////          matches_res[it->first.first][it->first.second] = it->second.size();
////          matches_res[it->first.second][it->first.first] = it->second.size();
////      }
////      std::ofstream out;
////      out.open("/media/add7/E/testData/match_martix_filter.res");
////      if(!out)
////      {
////          std::cout << "create matches res file error"<<std::endl;
////      }
////      for(std::size_t i = 0; i < matches_res.size(); ++i)
////      {
////          for(std::size_t j = 0; j < matches_res.size(); ++j)
////          {
////              out << matches_res[i][j] << " ";
////          }
////          out << std::endl;
////      }
////      out.close();












//  return EXIT_SUCCESS;
//}






