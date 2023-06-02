// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/global/sfm_global_engine_relative_motions.hpp"

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/graph/graph.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp"
#include "openMVG/sfm/sfm_data_BA.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/sfm/sfm_data_triangulation.hpp"
#include "openMVG/sfm/sfm_filters.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/tracks/tracks.hpp"
#include "openMVG/types.hpp"

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/cameras/Camera_undistort_image.hpp"

#include "third_party/histogram/histogram.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"

#include <ceres/types.h>

#include "ceres/rotation.h"

#include <iostream>
#include<string>

#include "third_party/progress/progress.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/graph/connectedComponent.hpp"


//test struct from pose
#include "openMVG/matching_image_collection/Pair_Builder.hpp"
#include "openMVG/sfm/pipelines/structure_from_known_poses/structure_estimator.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/features/svg_features.hpp"
#include "openMVG/geometry/frustum.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider_cache.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_filters_frustum.hpp"
#include "openMVG/sfm/sfm_report.hpp"

#include "openMVG/geometry/rigid_transformation3D_srt.hpp"
#include "openMVG/geometry/Similarity3.hpp"
#include "openMVG/sfm/sfm_data_transform.hpp"

#include <boost/algorithm/string.hpp>

#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

namespace openMVG{
namespace sfm{

using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::features;

//struct PoseStatus
//{
//    geometry::Pose3 pose;
////    int imgId;
//    int status;
//    int bdf;//123
//};
//using PoseStatusMap = Hash_Map<IndexT, PoseStatus>;  //<imageId, pose and status>

GlobalSfMReconstructionEngine_RelativeMotions::GlobalSfMReconstructionEngine_RelativeMotions(
  const SfM_Data & sfm_data,
  const std::string & soutDirectory,
  const std::string & sloggingFile)
  : ReconstructionEngine(sfm_data, soutDirectory), sLogging_file_(sloggingFile)
{

  if (!sLogging_file_.empty())
  {
    // setup HTML logger
    html_doc_stream_ = std::make_shared<htmlDocument::htmlDocumentStream>("GlobalReconstructionEngine SFM report.");
    html_doc_stream_->pushInfo(
      htmlDocument::htmlMarkup("h1", std::string("GlobalSfMReconstructionEngine_RelativeMotions")));
    html_doc_stream_->pushInfo("<hr>");

    html_doc_stream_->pushInfo( "Dataset info:");
    html_doc_stream_->pushInfo( "Views count: " +
      htmlDocument::toString( sfm_data.GetViews().size()) + "<br>");
  }

  // Set default motion Averaging methods
  eRotation_averaging_method_ = ROTATION_AVERAGING_L2;
  eTranslation_averaging_method_ = TRANSLATION_AVERAGING_L1;
}

GlobalSfMReconstructionEngine_RelativeMotions::~GlobalSfMReconstructionEngine_RelativeMotions()
{
  if (!sLogging_file_.empty())
  {
    // Save the reconstruction Log
    std::ofstream htmlFileStream(sLogging_file_.c_str());
    htmlFileStream << html_doc_stream_->getDoc();
  }
}

void GlobalSfMReconstructionEngine_RelativeMotions::SetFeaturesProvider(Features_Provider * provider)
{
  features_provider_ = provider;
}

void GlobalSfMReconstructionEngine_RelativeMotions::SetMatchesProvider(Matches_Provider * provider)
{
  matches_provider_ = provider;
}

void GlobalSfMReconstructionEngine_RelativeMotions::SetRotationAveragingMethod
(
  ERotationAveragingMethod eRotationAveragingMethod
)
{
  eRotation_averaging_method_ = eRotationAveragingMethod;
}

void GlobalSfMReconstructionEngine_RelativeMotions::SetTranslationAveragingMethod
(
  ETranslationAveragingMethod eTranslationAveragingMethod
)
{
  eTranslation_averaging_method_ = eTranslationAveragingMethod;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process() {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if (set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

  openMVG::rotation_averaging::RelativeRotations relatives_R;
  Compute_Relative_Rotations(relatives_R);

  Hash_Map<IndexT, Mat3> global_rotations;
  if (!Compute_Global_Rotations(relatives_R, global_rotations))
  {
    std::cerr << "GlobalSfM:: Rotation Averaging failure!" << std::endl;
    return false;
  }

  matching::PairWiseMatches  tripletWise_matches;
  if (!Compute_Global_Translations(global_rotations, tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
    return false;
  }
  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }
  if (!Adjust())
//  if (!Adjust_init_water_c1())
//  if (!Adjust_init_water_c4())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process_water_xy_z() {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
    std::cout << "error 5 0" << std::endl;
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if (set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

    std::cout << "error 5 " << std::endl;
  openMVG::rotation_averaging::RelativeRotations relatives_R;
  Compute_Relative_Rotations(relatives_R);

  std::cout << "error 6 " << std::endl;
  Hash_Map<IndexT, Mat3> global_rotations;
  if (!Compute_Global_Rotations(relatives_R, global_rotations))
  {
    std::cerr << "GlobalSfM:: Rotation Averaging failure!" << std::endl;
    return false;
  }

  std::cout << "error 7 " << std::endl;
  matching::PairWiseMatches  tripletWise_matches;
  if (!Compute_Global_Translations(global_rotations, tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
    return false;
  }
  std::cout << "error 8 " << std::endl;
  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }
  std::cout << "error 9 " << std::endl;
  if (!Adjust_nonremove())
//  if (!Adjust_init_water_c1())
//  if (!Adjust_init_water_c4())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  if (!Adjust_water_xy_z())
//  if (!Adjust_init_water_c1())
//  if (!Adjust_init_water_c4())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process_water_xy_z_synthetic() {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
//    std::cout << "error 5 0" << std::endl;
//  {
////    const Pair_Set pairs = matches_provider_->getPairs();
////    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
////    if (set_remainingIds.empty())
////    {
////      std::cout << "Invalid input image graph for global SfM" << std::endl;
////      return false;
////    }
////    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
//  }

//    std::cout << "error 5 " << std::endl;
//  openMVG::rotation_averaging::RelativeRotations relatives_R;
//  Compute_Relative_Rotations(relatives_R);

//  std::cout << "error 6 " << std::endl;
//  Hash_Map<IndexT, Mat3> global_rotations;
//  if (!Compute_Global_Rotations(relatives_R, global_rotations))
//  {
//    std::cerr << "GlobalSfM:: Rotation Averaging failure!" << std::endl;
//    return false;
//  }

//  std::cout << "error 7 " << std::endl;
//  matching::PairWiseMatches  tripletWise_matches;
//  if (!Compute_Global_Translations(global_rotations, tripletWise_matches))
//  {
//    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
//    return false;
//  }
//  std::cout << "error 8 " << std::endl;
//  if (!Compute_Initial_Structure(tripletWise_matches))
//  {
//    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
//    return false;
//  }
//  std::cout << "error 9 " << std::endl;
  if (!Adjust_nonremove())
//  if (!Adjust_init_water_c1())
//  if (!Adjust_init_water_c4())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  if (!Adjust_water_xy_z_n())
//  if (!Adjust_init_water_c1())
//  if (!Adjust_init_water_c4())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Process_water_xy_z_n() {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if (set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }


    //使用现成的sfm_bin的时候,把下面的东西注释掉 20210118
    //下面运行的是传统算法的相关内容
//  openMVG::rotation_averaging::RelativeRotations relatives_R;
//  Compute_Relative_Rotations(relatives_R);


//  Hash_Map<IndexT, Mat3> global_rotations;
//  if (!Compute_Global_Rotations(relatives_R, global_rotations))
//  {
//    std::cerr << "GlobalSfM:: Rotation Averaging failure!" << std::endl;
//    return false;
//  }

//  matching::PairWiseMatches  tripletWise_matches;
//  if (!Compute_Global_Translations(global_rotations, tripletWise_matches))
//  {
//    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
//    return false;
//  }

//  if (!Compute_Initial_Structure(tripletWise_matches))
//  {
//    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
//    return false;
//  }

//  std::cout << "error 9 " << std::endl;

//  if (!Adjust())
//  {
//    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
//    return false;
//  }
//使用现成的sfm_bin文件的时候，把上面的部分注释掉

  if (!Adjust_water_xy_z_n())
//  if (!Adjust_init_water_c1())
//  if (!Adjust_init_water_c4())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }


  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}



bool GlobalSfMReconstructionEngine_RelativeMotions::Process_step21() {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if (set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

  openMVG::rotation_averaging::RelativeRotations relatives_R;
  std::cout <<"start Compute_Relative_Rotations" << std::endl;
  Compute_Relative_Rotations(relatives_R);

  Hash_Map<IndexT, Mat3> global_rotations;
  std::cout <<"start Compute_Global_Rotations" << std::endl;
  if (!Compute_Global_Rotations(relatives_R, global_rotations))
  {
    std::cerr << "GlobalSfM:: Rotation Averaging failure!" << std::endl;
    return false;
  }

  matching::PairWiseMatches  tripletWise_matches;
  std::cout <<"start Compute_Global_Translations" << std::endl;
  if (!Compute_Global_Translations(global_rotations, tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
    return false;
  }
  std::cout <<"start Compute_Initial_Structure" << std::endl;
  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }
  std::cout <<"start Adjust_step21" << std::endl;
  if (!Adjust_step21())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process_step22() {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if (set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }
    matching::PairWiseMatches  tripletWise_matches;
    std::cout <<"start Compute_Initial_Structure" << std::endl;
    if (!Compute_Initial_Structure(tripletWise_matches))
    {
      std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
      return false;
    }

  if (!Adjust_step22())
//  if (!Adjust_init_water_c1())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process_forSimulation() {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if (set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }



    std::cout << "v size: " << sfm_data_.views.size() << std::endl;

    std::cout << "error 5" <<std::endl;
    openMVG::rotation_averaging::RelativeRotations relatives_R;
    Compute_Relative_Rotations(relatives_R);

    Hash_Map<IndexT,Mat3> global_rotation;
    if(!Compute_Global_Rotations(relatives_R, global_rotation))
    {
        return false;
    }
    matching::PairWiseMatches tripletWise_matches;
    if(!Compute_Global_Translations(global_rotation, tripletWise_matches))
    {
      return false;
    }

//    for (auto & view_it : sfm_data_.GetViews())
//    {
//        Vec3 C = sfm_data_.poses.at(view_it.first).center();
//        Vec3 w = Vec3(1.0,1.0,1.0);
//        dynamic_cast<sfm::ViewPriors*>(view_it.second.get())->SetPoseCenterPrior(C, w);

//    }

      std::cout <<"start Compute_Initial_Structure" << std::endl;
      if (!Compute_Initial_Structure(tripletWise_matches))
      {
        std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
        return false;
      }

      return true;



//    for (auto & iterView : sfm_data_.views)
//    {
//        sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(iterView.second.get());
//        Vec3 v1 = Vec3(2.0,2.0,2.0);
//        Vec3 v2 = Vec3(2.0,2.0,2.0);
//        std::cout << prior->pose_center_<< std::endl;
////        prior->SetPoseCenterPrior(v1, v2);
////        prior->pose_center_ = Vec3(2.0,2.0,2.0);
//        sfm_data_.poses.at(iterView.first).center();

//    }
//  if (!Adjust_step22())
//  {}
  if (!Adjust())
//  {}
//  if (!Adjust_init_water_c1())
//  if (!Adjust_init_water_c4())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process_forSimulation_afterOptim() {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
//    const Pair_Set pairs = matches_provider_->getPairs();
//    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
//    if (set_remainingIds.empty())
//    {
//      std::cout << "Invalid input image graph for global SfM" << std::endl;
//      return false;
//    }
//    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }
//    std::cout << "v size: " << sfm_data_.views.size() << std::endl;

//    std::cout << "error 5" <<std::endl;
//    openMVG::rotation_averaging::RelativeRotations relatives_R;
//    Compute_Relative_Rotations(relatives_R);

//    Hash_Map<IndexT,Mat3> global_rotation;
//    if(!Compute_Global_Rotations(relatives_R, global_rotation))
//    {
//        return false;
//    }
    matching::PairWiseMatches tripletWise_matches;
//    if(!Compute_Global_Translations(global_rotation, tripletWise_matches))
//    {
//      return false;
//    }

//    for (auto & view_it : sfm_data_.GetViews())
//    {
//        Vec3 C = sfm_data_.poses.at(view_it.first).center();
//        Vec3 w = Vec3(1.0,1.0,1.0);
//        dynamic_cast<sfm::ViewPriors*>(view_it.second.get())->SetPoseCenterPrior(C, w);

//    }

//      std::cout <<"start Compute_Initial_Structure" << std::endl;
//      if (!Compute_Initial_Structure(tripletWise_matches))
//      {
//        std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
//        return false;
//      }

//      return true;



//    for (auto & iterView : sfm_data_.views)
//    {
//        sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(iterView.second.get());
//        Vec3 v1 = Vec3(2.0,2.0,2.0);
//        Vec3 v2 = Vec3(2.0,2.0,2.0);
//        std::cout << prior->pose_center_<< std::endl;
////        prior->SetPoseCenterPrior(v1, v2);
////        prior->pose_center_ = Vec3(2.0,2.0,2.0);
//        sfm_data_.poses.at(iterView.first).center();

//    }
//  if (!Adjust_step22())
//  {}
  if (!Adjust())
//  {}
//  if (!Adjust_init_water_c1())
//  if (!Adjust_init_water_c4())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Process_water6() {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if (set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

    openMVG::rotation_averaging::RelativeRotations relatives_R;
    Compute_Relative_Rotations(relatives_R);

    Hash_Map<IndexT, Mat3> global_rotations;
    if (!Compute_Global_Rotations(relatives_R, global_rotations))
    {
      std::cerr << "GlobalSfM:: Rotation Averaging failure!" << std::endl;
      return false;
    }

    matching::PairWiseMatches  tripletWise_matches;
    if (!Compute_Global_Translations(global_rotations, tripletWise_matches))
    {
      std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
      return false;
    }


    std::cout << "error 2 " << std::endl;

    std::cout << "v size: " << sfm_data_.views.size() << std::endl;

    std::cout << "error 5" <<std::endl;
//    openMVG::rotation_averaging::RelativeRotations relatives_R;
//    Compute_Relative_Rotations(relatives_R);

//    Hash_Map<IndexT,Mat3> global_rotation;
//    if(!Compute_Global_Rotations(relatives_R, global_rotation))
//    {
//        return false;
//    }
//    if(!Compute_Global_Translations(global_rotation, tripletWise_matches))
//    {
//      return false;
//    }

      std::cout <<"start Compute_Initial_Structure" << std::endl;
      if (!Compute_Initial_Structure(tripletWise_matches))
      {
        std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
        return false;
      }

//      return true;

      for (const auto & pose_it : sfm_data_.poses)
      {
          std::cout << pose_it.first << " ";

      }
      std::cout << std::endl;
      {
//          std::vector<int> useIdList;
//          useIdList.push_back(0);//0
//          useIdList.push_back(2);//1
//          useIdList.push_back(3);//2
//          useIdList.push_back(5);//3
//          useIdList.push_back(7);//4
//          useIdList.push_back(9);//5
//          useIdList.push_back(11);//6
//          useIdList.push_back(13);//7
//          useIdList.push_back(14);//8
//          useIdList.push_back(15);//9
//          useIdList.push_back(17);//10
//          useIdList.push_back(18);//11
//          useIdList.push_back(20);//12
//          useIdList.push_back(21);//13
//          useIdList.push_back(23);//14
//          useIdList.push_back(25);//15
//          useIdList.push_back(26);//16
//          useIdList.push_back(29);//17
//          useIdList.push_back(32);//18
//          SfM_Data newdata;
//          newdata.intrinsics = sfm_data_.intrinsics;
//          newdata.s_root_path = sfm_data_.s_root_path;

//          for (const auto & pose_it : sfm_data_.poses)
//          {
//            const IndexT indexPose = pose_it.first;
//            bool find = false;
//            std::size_t newId;
//            for(std::size_t tUILId = 0; tUILId < useIdList.size(); ++tUILId)
//            {
//                if(useIdList[tUILId] == indexPose)
//                {
//                    find = true;
//                    newId = tUILId;
//                    break;
//                }
//            }
//            if (find)
//            {
//                newdata.poses[newId] = pose_it.second;
//            }
//          }

//          //
//          for (const auto & view_it : sfm_data_.views)
//          {
//            const IndexT indexView = view_it.first;
//            bool find = false;
//            std::size_t newId;
//            for(std::size_t tUILId = 0; tUILId < useIdList.size(); ++tUILId)
//            {
//                if(useIdList[tUILId] == indexView)
//                {
//                    find = true;
//                    newId = tUILId;
//                    break;
//                }
//            }
//            if (find)
//            {
////                std::shared_ptr<View> v;
//                std::shared_ptr<View> v = std::make_shared<View>(
//                            view_it.second.get()->s_Img_path,
//                            newId,
//                            view_it.second.get()->id_intrinsic,
//                            newId,
//                            view_it.second.get()->ui_width,
//                            view_it.second.get()->ui_height);

//                newdata.views[newId] = v;//view_it.second;
////                newdata.views[newId]->id_view = newId;
////                newdata.views[newId]->id_pose = newId;
//            }
//          }
//          //
//          for(Landmarks::iterator itX = sfm_data_.structure.begin();
//              itX != sfm_data_.structure.end(); ++itX)
//          {
//              Observations tempObs;
//              for(Observations::iterator itObs = itX->second.obs.begin();
//                  itObs != itX->second.obs.end(); ++itObs)
//              {
//                  const View * view = sfm_data_.views.at(itObs->first).get();
//                  View * view2 = sfm_data_.views[itObs->first].get();
//                  std::cout << "view : " << view->id_view << " view2 : " << view2->id_view << " itObs->first : " << itObs->first <<std::endl;

//                  bool find = false;
//                  std::size_t newId;
//                  for(std::size_t tUILId = 0; tUILId < useIdList.size(); ++tUILId)
//                  {
//                      if(useIdList[tUILId] == view->id_view)
//                      {
//                          find = true;
//                          newId = tUILId;
//                          break;
//                      }
//                  }
//                  if(find)
//                  {
//                      tempObs[newId] = itObs->second;
//                  }
//              }
//              //
//              if(tempObs.size() >= 3)
//              {
//                  Landmark tempL;
//                  tempL.obs = tempObs;
//                  tempL.X = itX->second.X;
//                  newdata.structure[itX->first] = tempL;
//              }
//          }
//          //
//          sfm_data_ = newdata;

      }
//      return true;


//      std::cout << "error 7" <<std::endl;
//      for (auto & view_it : sfm_data_.GetViews())
//      {
//          std::cout << "error 7 " << view_it.first <<std::endl;
//  //        sfm_data_.poses.at(iterView.first).center();
//          Vec3 C = sfm_data_.poses.at(view_it.first).center();
//          Vec3 w = Vec3(1.0,1.0,1.0);
//  //        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
//          dynamic_cast<sfm::ViewPriors*>(view_it.second.get())->SetPoseCenterPrior(C, w);
//  //      sfm::ViewPriors prior = (view_it.second.get());

//      }

      std::cout << "error 8" <<std::endl;
//    for (auto & iterView : sfm_data_.views)
//    {
//        sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(iterView.second.get());
//        Vec3 v1 = Vec3(2.0,2.0,2.0);
//        Vec3 v2 = Vec3(2.0,2.0,2.0);
//        std::cout << prior->pose_center_<< std::endl;
////        prior->SetPoseCenterPrior(v1, v2);
////        prior->pose_center_ = Vec3(2.0,2.0,2.0);
//        sfm_data_.poses.at(iterView.first).center();

//    }
//  if (!Adjust_step22())
//  {}
//  if (!Adjust())
//  {}
  if (!Adjust_init_water_c4())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}



bool GlobalSfMReconstructionEngine_RelativeMotions::Process_GCP() {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if (set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

  openMVG::rotation_averaging::RelativeRotations relatives_R;
  Compute_Relative_Rotations(relatives_R);

  Hash_Map<IndexT, Mat3> global_rotations;
  if (!Compute_Global_Rotations(relatives_R, global_rotations))
  {
    std::cerr << "GlobalSfM:: Rotation Averaging failure!" << std::endl;
    return false;
  }

  matching::PairWiseMatches  tripletWise_matches;
  if (!Compute_Global_Translations(global_rotations, tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
    return false;
  }
  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }
  if (!Adjust_GCP())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process_GCP_GPS(std::string rtFilePath) {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if (set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

    Hash_Map<IndexT, Mat3> global_rotations;

    //init inter-constraint
    Mat3 Rb_trans, Rf_trans, Rd_trans;
    Vec3 tb_trans, tf_trans, td_trans;
    Mat3 R_imu_trans;
    Vec3 gps_trans;
    std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
    const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
    const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};
    {
        //test distortion
        Mat3 r_tb, r_tf;
        r_tb << 1.0, 0.0, 0.0,
                0.0, 0.787926, -0.61577,
                0.0, 0.61577, 0.787926;
        r_tf << 1.0, 0.0, 0.0,
                0.0, 0.78796,  0.615727,
                0.0, -0.615727, 0.78796;

  //      Rx << 1.0, 0.0, 0.0,
  //             0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //             0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


  //      r_tb << 1.0, 0.0, 0.0,
  //              0.0, cos(38.0/180.0*M_PI), -sin(38.0/180.0*M_PI),
  //              0.0, sin(38.0/180.0*M_PI), cos(38.0/180.0*M_PI);

  //      r_tf << 1.0, 0.0, 0.0,
  //              0.0, cos(-38.0/180.0*M_PI), -sin(-38.0/180.0*M_PI),
  //              0.0, sin(-38.0/180.0*M_PI), cos(-38.0/180.0*M_PI);

  //      r_tb << cos(-38.0/180.0*M_PI), 0.0, sin(-38.0/180.0*M_PI),//-38
  //              0.0, 1.0, 0.0,
  //              -sin(-38.0/180.0*M_PI), 0.0, cos(-38.0/180.0*M_PI);
  //      r_tf << cos(38.0/180.0*M_PI), 0.0, sin(38.0/180.0*M_PI),//38
  //              0.0, 1.0, 0.0,
  //              -sin(38.0/180.0*M_PI), 0.0, cos(38.0/180.0*M_PI);
        Vec3 RA_tb, RA_tf;
        ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
        ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

        //set reference
        double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), -RA_tb(2)};
        double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), -RA_tf(2)};
  //      double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), RA_tb(2)};
  //      double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), RA_tf(2)};
        double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};

        Vec3 c_bd, c_fd;//b-d, f-d
        c_bd << 0.0, 0.098, 0.027;
        c_fd << 0.0, -0.101, 0.028;
  //      c_bd << -0.098, 0.0, -0.027;
  //      c_fd << 0.101, 0.0, -0.028;

  //      c_bd << 0.098, 0.0, -0.027;
  //      c_fd << -0.101, 0.0, -0.028;
  //      c_bd << 0.0, 1.0, 0.027;
  //      c_fd << 0.0, -0.101, 0.028;

        Vec3 t_db, t_df, t_dd;
        t_db = -r_tb*c_bd;
        t_df = -r_tf*c_fd;

        //double S = 1470.531783220221;
      //        t_db << 0.0, 0.098, 0.027;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
      //        t_df << 0.0, -0.101, 0.028;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
        t_dd << 0.0, 0.0, 0.0;

        transformation_br = {R_br_angleAxis[0], R_br_angleAxis[1], R_br_angleAxis[2],
                             t_db(0), t_db(1), t_db(2)};
        transformation_fr = {R_fr_angleAxis[0], R_fr_angleAxis[1], R_fr_angleAxis[2],
                             t_df(0), t_df(1), t_df(2)};
        transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                             t_dd(0), t_dd(1), t_dd(2)};

  //      const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
  //      const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;

        //prepare for setting pose
        ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rb_trans.data());
        tb_trans = Vec3(transformation_br[3], transformation_br[4], transformation_br[5]);
        ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rf_trans.data());
        tf_trans = Vec3(transformation_fr[3], transformation_fr[4], transformation_fr[5]);
        ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rd_trans.data());
        td_trans = Vec3(transformation_dr[3], transformation_dr[4], transformation_dr[5]);


        ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
        gps_trans = Vec3(tramsformation_x_gps[0], tramsformation_x_gps[1], tramsformation_x_gps[2]);

    }


    //get gps and imu data
    PoseStatusMap poseStatusMap;
    Hash_Map<IndexT, std::vector<double> > C_gps_Map;
    Hash_Map<IndexT, std::vector<double> > R_imu_Map;
    std::vector<double> c0(3);
    {
      //        /Hash_Map<IndexT, std::vector<double> > R_imu;
  //            Hash_Map<IndexT, std::vector<double> > C_gps;
        std::string tFilePath = rtFilePath;  //!!!
        std::ifstream in;
        in.open(tFilePath);
        if(!in)
        {
            std::cout << "open tFilePath file error!" << std::endl;
            return false;
        }
        std::size_t line = 0;
        int count = -1;
        in >> count;  ///!!!

        bool isC0 = false;//true;

        while((!in.eof()) && (line < count))
        {
            // input
            int imgId, imgStatus;
            std::vector<double> c(3);
            //double pitch, roll, yaw;
            double omega, phi, kappa;
            in >> imgId
               >> c[0] >> c[1] >> c[2]
  //             >> c[1] >> c[0] >> c[2]
               >> omega >> phi >> kappa
               >> imgStatus;

  //          omega = omega + 2.4;
  //          phi = phi + 0.0713;
  //          kappa = kappa - 0.5805;
  //          kappa = kappa + 2.4;
  //          phi = phi + 0.0713;
  //          omega = omega - 0.5805;

            if(isC0 == false)
            {
                c0 = c;
                isC0 = true;
            }

//            c[0] -= c0[0];
//            c[1] -= c0[1];
//            c[2] -= c0[2];

            ++line;



  //          double roll, pitch, yaw;

  //          roll = omega;
  //          pitch = phi;
  //          yaw = kappa;

  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
  //                sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
  //                0.0, 1.0, 0.0,
  //                -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
  //          Rx << 1.0, 0.0, 0.0,
  //                0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
  //                0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);
  //          Mat3 R_imu =  Rz*Rx*Ry;
  //          Mat3 changeAxis_M;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, -1.0, 0.0,
  //                          0.0, 0.0, -1.0;
  //          R_imu = changeAxis_M * R_imu;




  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
  //                0.0, 1.0, 0.0,
  //                -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);
  //          Rx << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

  ////          Mat3 R_imu =  Rx*Ry*Rz;
  ////          Mat3 R_imu =  Rx*Rz*Ry;//
  ////          Mat3 R_imu =  Ry*Rz*Rx;
  ////          Mat3 R_imu =  Ry*Rx*Rz;
  //          Mat3 R_imu =  Rz*Rx*Ry;//
  ////          Mat3 R_imu =  Rz*Ry*Rx;


            //use now
            Mat3 Rz, Ry, Rx;
            Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
                  sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
                  0.0, 0.0, 1.0;
            Ry << cos(phi/180.0*M_PI), 0.0, sin(phi/180.0*M_PI),
                  0.0, 1.0, 0.0,
                  -sin(phi/180.0*M_PI), 0.0, cos(phi/180.0*M_PI);
            Rx << 1.0, 0.0, 0.0,
                  0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
                  0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
  //          Mat3 R_imu =  Rz*Ry*Rx;
            Mat3 R_imu =  Rx*Ry*Rz;



  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, -cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), -cos(phi/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;

  //          Rz << cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI), 0.0,
  //                sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, -cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
  //                0.0, sin(omega/180.0*M_PI), -cos(omega/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;

  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;          Mat3 Rz, Ry, Rx;
  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(phi/180.0*M_PI), sin(phi/180.0*M_PI), 0.0,
  //                -sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
  //                0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;



  //          Mat3 R_imu =  Ry*Rz*Rx;
  //          Mat3 R_imu =  Rz*Ry*Rx;
  //          Mat3 R_imu =  Rx*Ry*Rz;


  //          Mat3 Rk, Ra, RA;
  //          Rk << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          RA << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ra << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


  //          Mat3 R_imu =  Ry*Rx*Rz;

  //          Mat3 changeAxis_M;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, -1.0, 0.0,
  //                          0.0, 0.0, -1.0;  //M1
  //          changeAxis_M << 0.0, 1.0, 0.0,
  //                          1.0, 0.0, 0.0,
  //                          0.0, 0.0, -1.0;  //M2

  //          R_imu = changeAxis_M * R_imu;
  //          R_imu = R_imu * changeAxis_M;





  //              tempGet = getListXYZ(-phi3, kappa3, -omega3);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //               tempGet = getListXYZ(phi3, kappa3, omega3);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //          tempGet = getListXYZ(phi1, kappa1, omega1);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //          Mat3 Rz, Ry, Rx;
  //          Rx << 1.0, 0.0, 0.0,
  //                 0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                 0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

  //          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
  //                 0.0, 1.0, 0.0,
  //                 -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);

  //          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                 sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                 0.0, 0.0, 1.0;


  //          Mat3 R_imu =  Ry*Rz*Rx;
  //          Mat3 R_imu =  Rz*Rx*Ry;
  //          Mat3 R_imu =  Rx*Rz*Ry;
  //          Mat3 R_imu =  Rx*Ry*Rz;

            Mat3 changeAxis_M;
            changeAxis_M << 1.0, 0.0, 0.0,
                            0.0, -1.0, 0.0,
                            0.0, 0.0, -1.0;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, 1.0, 0.0,
  //                          0.0, 0.0, -1.0;

  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, 1.0, 0.0,
  //                          0.0, 0.0, -1.0;
            Mat3 R_;
            R_ << 1.0, 0.0, 0.0,
                  0.0, 1.0, 0.0,
                  0.0, 0.0, -1.0;
            R_imu = changeAxis_M * R_imu;//* R_;


            //Mat3 R2;
            //ceres::AngleAxisToRotationMatrix(R2_a.data(), R2.data());

            //R_imu = changeAxis_M * R_imu * R2;//* R_;


            std::vector<double> R_imu_a(3);
            ceres::RotationMatrixToAngleAxis(R_imu.data(), R_imu_a.data());

            Vec3 C_gps_res = Vec3(c[0], c[1], c[2]);


            {
                bool find = false;
                for(Views::const_iterator itV = sfm_data_.views.begin(); itV != sfm_data_.views.end(); ++itV)
                {
                    if(itV->first == imgId)
                    {
                        find = true;
                        break;
                    }
                }
                if(find == false)
                    continue;
            }
            // check status
            if(imgId/100000 == 1)//back
            {
                if(imgStatus != 321)
                {
                    C_gps_Map[imgId + 100000] = c;
                    R_imu_Map[imgId + 100000] = R_imu_a;
                }

                //set back-camera pose
                Mat3 Rb_res = Rb_trans * R_imu_trans * R_imu;

                Vec3 tb_res;
                tb_res = Rb_trans * gps_trans - Rb_trans * R_imu_trans * R_imu * C_gps_res + tb_trans;
                global_rotations[imgId] = Rb_res;

                sfm_data_.poses[imgId] = Pose3(Rb_res, -Rb_res.transpose() * tb_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 1;
                poseStatusMap[imgId] = poseStatus;

            }else if(imgId/100000 == 2)//down
            {
                C_gps_Map[imgId] = c;
                R_imu_Map[imgId] = R_imu_a;

                //set down-camera pose
                Mat3 Rd_res = Rd_trans * R_imu_trans * R_imu;

                Vec3 td_res;
                td_res = Rd_trans * gps_trans - Rd_trans * R_imu_trans * R_imu * C_gps_res + td_trans;
                global_rotations[imgId] = Rd_res;

                sfm_data_.poses[imgId] = Pose3(Rd_res, -Rd_res.transpose() * td_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 2;
                poseStatusMap[imgId] = poseStatus;

            }else if(imgId/100000 == 3)//front
            {
                if(imgStatus != 321)
                {
                    C_gps_Map[imgId - 100000] = c;
                    R_imu_Map[imgId - 100000] = R_imu_a;
                }

                //set front-camera pose
                Mat3 Rf_res = Rf_trans * R_imu_trans * R_imu;

                Vec3 tf_res;
                tf_res = Rf_trans * gps_trans - Rf_trans * R_imu_trans * R_imu * C_gps_res + tf_trans;
                global_rotations[imgId] = Rf_res;

                sfm_data_.poses[imgId] = Pose3(Rf_res, -Rf_res.transpose() * tf_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 3;
                poseStatusMap[imgId] = poseStatus;

            }

        }
        in.close();

    }


    matching::PairWiseMatches tripletWise_matches;
    if (!Compute_Global_Translations_gps(global_rotations, tripletWise_matches))
    {
      std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
      return false;
    }

  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }
  if (!Adjust_init())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process_GCP_GPS_BAStep_save(std::string rtFilePath, std::string sOutDir, int bId) {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if (set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

    Hash_Map<IndexT, Mat3> global_rotations;

    //init inter-constraint
    Mat3 Rb_trans, Rf_trans, Rd_trans;
    Vec3 tb_trans, tf_trans, td_trans;
    Mat3 R_imu_trans;
    Vec3 gps_trans;
    std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
    const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
    const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};
    {
        //test distortion
        Mat3 r_tb, r_tf;
        r_tb << 1.0, 0.0, 0.0,
                0.0, 0.787926, -0.61577,
                0.0, 0.61577, 0.787926;
        r_tf << 1.0, 0.0, 0.0,
                0.0, 0.78796,  0.615727,
                0.0, -0.615727, 0.78796;

//        r_tb << 1.0, 0.0, 0.0,
//                0.0, 0.787926, -0.61577,
//                0.0, 0.61577, 0.787926;
//        r_tf << 1.0, 0.0, 0.0,
//                0.0, 0.78796,  0.615727,
//                0.0, -0.615727, 0.78796;

  //      Rx << 1.0, 0.0, 0.0,
  //             0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //             0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


  //      r_tb << 1.0, 0.0, 0.0,
  //              0.0, cos(38.0/180.0*M_PI), -sin(38.0/180.0*M_PI),
  //              0.0, sin(38.0/180.0*M_PI), cos(38.0/180.0*M_PI);

  //      r_tf << 1.0, 0.0, 0.0,
  //              0.0, cos(-38.0/180.0*M_PI), -sin(-38.0/180.0*M_PI),
  //              0.0, sin(-38.0/180.0*M_PI), cos(-38.0/180.0*M_PI);

  //      r_tb << cos(-38.0/180.0*M_PI), 0.0, sin(-38.0/180.0*M_PI),//-38
  //              0.0, 1.0, 0.0,
  //              -sin(-38.0/180.0*M_PI), 0.0, cos(-38.0/180.0*M_PI);
  //      r_tf << cos(38.0/180.0*M_PI), 0.0, sin(38.0/180.0*M_PI),//38
  //              0.0, 1.0, 0.0,
  //              -sin(38.0/180.0*M_PI), 0.0, cos(38.0/180.0*M_PI);
        Vec3 RA_tb, RA_tf;
        ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
        ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

        //set reference
        double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), -RA_tb(2)};
        double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), -RA_tf(2)};
  //      double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), RA_tb(2)};
  //      double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), RA_tf(2)};
        double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};

        Vec3 c_bd, c_fd;//b-d, f-d
//        c_bd << 0.0, 0.0, 0.0;
//        c_fd << 0.0, 0.0, 0.0;

        c_bd << 0.0, 0.098, 0.027;
        c_fd << 0.0, -0.101, 0.028;

  //      c_bd << -0.098, 0.0, -0.027;
  //      c_fd << 0.101, 0.0, -0.028;

  //      c_bd << 0.098, 0.0, -0.027;
  //      c_fd << -0.101, 0.0, -0.028;
  //      c_bd << 0.0, 1.0, 0.027;
  //      c_fd << 0.0, -0.101, 0.028;

        Vec3 t_db, t_df, t_dd;
        t_db = -r_tb*c_bd;
        t_df = -r_tf*c_fd;

        //double S = 1470.531783220221;
      //        t_db << 0.0, 0.098, 0.027;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
      //        t_df << 0.0, -0.101, 0.028;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
        t_dd << 0.0, 0.0, 0.0;

        transformation_br = {R_br_angleAxis[0], R_br_angleAxis[1], R_br_angleAxis[2],
                             t_db(0), t_db(1), t_db(2)};
        transformation_fr = {R_fr_angleAxis[0], R_fr_angleAxis[1], R_fr_angleAxis[2],
                             t_df(0), t_df(1), t_df(2)};
        transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                             t_dd(0), t_dd(1), t_dd(2)};

  //      const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
  //      const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;

        //prepare for setting pose
        ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rb_trans.data());
        tb_trans = Vec3(transformation_br[3], transformation_br[4], transformation_br[5]);
        ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rf_trans.data());
        tf_trans = Vec3(transformation_fr[3], transformation_fr[4], transformation_fr[5]);
        ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rd_trans.data());
        td_trans = Vec3(transformation_dr[3], transformation_dr[4], transformation_dr[5]);


        ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
        gps_trans = Vec3(tramsformation_x_gps[0], tramsformation_x_gps[1], tramsformation_x_gps[2]);

    }


    //get gps and imu data
    PoseStatusMap poseStatusMap;
    Hash_Map<IndexT, std::vector<double> > C_gps_Map;
    Hash_Map<IndexT, std::vector<double> > R_imu_Map;
    std::vector<double> c0(3);
    {
      //        /Hash_Map<IndexT, std::vector<double> > R_imu;
  //            Hash_Map<IndexT, std::vector<double> > C_gps;
        std::string tFilePath = rtFilePath;  //!!!
        std::ifstream in;
        in.open(tFilePath);
        if(!in)
        {
            std::cout << "open tFilePath file error!" << std::endl;
            return false;
        }
        std::size_t line = 0;
        int count = -1;
        in >> count;  ///!!!

        bool isC0 = false;//true;

        while((!in.eof()) && (line < count))
        {
            // input
            int imgId, imgStatus;
            std::vector<double> c(3);
            //double pitch, roll, yaw;
            double omega, phi, kappa;
            in >> imgId
               >> c[0] >> c[1] >> c[2]
  //             >> c[1] >> c[0] >> c[2]
               >> omega >> phi >> kappa
               >> imgStatus;

  //          omega = omega + 2.4;
  //          phi = phi + 0.0713;
  //          kappa = kappa - 0.5805;
  //          kappa = kappa + 2.4;
  //          phi = phi + 0.0713;
  //          omega = omega - 0.5805;

            if(isC0 == false)
            {
                c0 = c;
                isC0 = true;
            }

//            c[0] -= c0[0];
//            c[1] -= c0[1];
//            c[2] -= c0[2];

            ++line;



  //          double roll, pitch, yaw;

  //          roll = omega;
  //          pitch = phi;
  //          yaw = kappa;

  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
  //                sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
  //                0.0, 1.0, 0.0,
  //                -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
  //          Rx << 1.0, 0.0, 0.0,
  //                0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
  //                0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);
  //          Mat3 R_imu =  Rz*Rx*Ry;
  //          Mat3 changeAxis_M;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, -1.0, 0.0,
  //                          0.0, 0.0, -1.0;
  //          R_imu = changeAxis_M * R_imu;




  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
  //                0.0, 1.0, 0.0,
  //                -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);
  //          Rx << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

  ////          Mat3 R_imu =  Rx*Ry*Rz;
  ////          Mat3 R_imu =  Rx*Rz*Ry;//
  ////          Mat3 R_imu =  Ry*Rz*Rx;
  ////          Mat3 R_imu =  Ry*Rx*Rz;
  //          Mat3 R_imu =  Rz*Rx*Ry;//
  ////          Mat3 R_imu =  Rz*Ry*Rx;


            //use now
            Mat3 Rz, Ry, Rx;
            Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
                  sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
                  0.0, 0.0, 1.0;
            Ry << cos(phi/180.0*M_PI), 0.0, sin(phi/180.0*M_PI),
                  0.0, 1.0, 0.0,
                  -sin(phi/180.0*M_PI), 0.0, cos(phi/180.0*M_PI);
            Rx << 1.0, 0.0, 0.0,
                  0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
                  0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
  //          Mat3 R_imu =  Rz*Ry*Rx;
            Mat3 R_imu =  Rx*Ry*Rz;



  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, -cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), -cos(phi/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;

  //          Rz << cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI), 0.0,
  //                sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, -cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
  //                0.0, sin(omega/180.0*M_PI), -cos(omega/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;

  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;          Mat3 Rz, Ry, Rx;
  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(phi/180.0*M_PI), sin(phi/180.0*M_PI), 0.0,
  //                -sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
  //                0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;



  //          Mat3 R_imu =  Ry*Rz*Rx;
  //          Mat3 R_imu =  Rz*Ry*Rx;
  //          Mat3 R_imu =  Rx*Ry*Rz;


  //          Mat3 Rk, Ra, RA;
  //          Rk << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          RA << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ra << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


  //          Mat3 R_imu =  Ry*Rx*Rz;

  //          Mat3 changeAxis_M;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, -1.0, 0.0,
  //                          0.0, 0.0, -1.0;  //M1
  //          changeAxis_M << 0.0, 1.0, 0.0,
  //                          1.0, 0.0, 0.0,
  //                          0.0, 0.0, -1.0;  //M2

  //          R_imu = changeAxis_M * R_imu;
  //          R_imu = R_imu * changeAxis_M;





  //              tempGet = getListXYZ(-phi3, kappa3, -omega3);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //               tempGet = getListXYZ(phi3, kappa3, omega3);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //          tempGet = getListXYZ(phi1, kappa1, omega1);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //          Mat3 Rz, Ry, Rx;
  //          Rx << 1.0, 0.0, 0.0,
  //                 0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                 0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

  //          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
  //                 0.0, 1.0, 0.0,
  //                 -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);

  //          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                 sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                 0.0, 0.0, 1.0;


  //          Mat3 R_imu =  Ry*Rz*Rx;
  //          Mat3 R_imu =  Rz*Rx*Ry;
  //          Mat3 R_imu =  Rx*Rz*Ry;
  //          Mat3 R_imu =  Rx*Ry*Rz;

            Mat3 changeAxis_M;
            changeAxis_M << 1.0, 0.0, 0.0,
                            0.0, -1.0, 0.0,
                            0.0, 0.0, -1.0;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, 1.0, 0.0,
  //                          0.0, 0.0, -1.0;

  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, 1.0, 0.0,
  //                          0.0, 0.0, -1.0;
            Mat3 R_;
            R_ << 1.0, 0.0, 0.0,
                  0.0, 1.0, 0.0,
                  0.0, 0.0, -1.0;
            R_imu = changeAxis_M * R_imu;//* R_;


            //Mat3 R2;
            //ceres::AngleAxisToRotationMatrix(R2_a.data(), R2.data());

            //R_imu = changeAxis_M * R_imu * R2;//* R_;


            std::vector<double> R_imu_a(3);
            ceres::RotationMatrixToAngleAxis(R_imu.data(), R_imu_a.data());

            Vec3 C_gps_res = Vec3(c[0], c[1], c[2]);


            {
                bool find = false;
                for(Views::const_iterator itV = sfm_data_.views.begin(); itV != sfm_data_.views.end(); ++itV)
                {
                    if(itV->first == imgId)
                    {
                        find = true;
                        break;
                    }
                }
                if(find == false)
                    continue;
            }
            // check status
            if(imgId/100000 == 1)//back
            {
                if(imgStatus != 321)
                {
                    C_gps_Map[imgId + 100000] = c;
                    R_imu_Map[imgId + 100000] = R_imu_a;
                }

                //set back-camera pose
                Mat3 Rb_res = Rb_trans * R_imu_trans * R_imu;

                Vec3 tb_res;
                tb_res = Rb_trans * gps_trans - Rb_trans * R_imu_trans * R_imu * C_gps_res + tb_trans;
                global_rotations[imgId] = Rb_res;

                sfm_data_.poses[imgId] = Pose3(Rb_res, -Rb_res.transpose() * tb_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 1;
                poseStatusMap[imgId] = poseStatus;

            }else if(imgId/100000 == 2)//down
            {
                C_gps_Map[imgId] = c;
                R_imu_Map[imgId] = R_imu_a;

                //set down-camera pose
                Mat3 Rd_res = Rd_trans * R_imu_trans * R_imu;

                Vec3 td_res;
                td_res = Rd_trans * gps_trans - Rd_trans * R_imu_trans * R_imu * C_gps_res + td_trans;
                global_rotations[imgId] = Rd_res;

                sfm_data_.poses[imgId] = Pose3(Rd_res, -Rd_res.transpose() * td_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 2;
                poseStatusMap[imgId] = poseStatus;

            }else if(imgId/100000 == 3)//front
            {
                if(imgStatus != 321)
                {
                    C_gps_Map[imgId - 100000] = c;
                    R_imu_Map[imgId - 100000] = R_imu_a;
                }

                //set front-camera pose
                Mat3 Rf_res = Rf_trans * R_imu_trans * R_imu;

                Vec3 tf_res;
                tf_res = Rf_trans * gps_trans - Rf_trans * R_imu_trans * R_imu * C_gps_res + tf_trans;
                global_rotations[imgId] = Rf_res;

                sfm_data_.poses[imgId] = Pose3(Rf_res, -Rf_res.transpose() * tf_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 3;
                poseStatusMap[imgId] = poseStatus;

            }

        }
        in.close();

    }


    matching::PairWiseMatches tripletWise_matches;
    if (!Compute_Global_Translations_gps(global_rotations, tripletWise_matches))
    {
      std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
      return false;
    }

  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }
  std::ofstream out_testPose, out_testView;
  out_testPose.open(sOutDir+"/out_testPose.res");
  out_testView.open(sOutDir+"/out_testView.res");

  {
      std::stringstream ss;
      ss << bId;
      std::string sBId;
      ss >> sBId;
      std::ofstream out_landmark;
      std::string out_landmark_path = sOutDir + "/landmarks_" + sBId + ".res";
      out_landmark.open(out_landmark_path);
      if(!out_landmark)
      {
          std::cout << "create " << out_landmark_path << " file failed!" << std::endl;
      }
//      std::ofstream out_X;
//      std::string out_X_path = sfm_data_.s_root_path + "/X_all.res";
//      out_X.open(out_X_path);
//      if(!out_X)
//      {
//          std::cout << "create " << out_X_path << " file failed!" << std::endl;
//      }
      //
      int countL = sfm_data_.structure.size();
      out_landmark << countL << std::endl;
      //
      for(Landmarks::const_iterator itL = sfm_data_.structure.begin(); itL != sfm_data_.structure.end(); ++itL)
      {
          int Xid = itL->first;
          Vec3 X = itL->second.X;
          out_landmark << Xid << " " << X(0) << " " << X(1) << " " << X(2) << std::endl;
          //
          int countObs = itL->second.obs.size();
          out_landmark << countObs << " ";
          for(Observations::const_iterator itObs = itL->second.obs.begin(); itObs != itL->second.obs.end(); ++itObs)
          {
              int imgId = itObs->first;
              Vec2 featId = itObs->second.x;
              out_landmark << imgId;
              out_landmark << std::setprecision(8) << " " << featId(0) << " " << featId(1) << " ";
          }
          out_landmark << std::endl;
      }
      out_landmark.close();
      //
      /// save extrinsics
      std::ofstream out_extrinsic;
      std::string out_extrinsic_path = sOutDir + "/extrinsics_" + sBId + ".res";
      out_extrinsic.open(out_extrinsic_path);
      if(!out_extrinsic)
      {
          std::cout << "create " << out_extrinsic_path << " file failed!" << std::endl;
      }
      out_extrinsic << sfm_data_.poses.size() << std::endl;
      for(Poses::const_iterator itPose = sfm_data_.poses.begin(); itPose != sfm_data_.poses.end(); ++itPose)
      {
          Vec3 t = itPose->second.translation();
          Mat3 R = itPose->second.rotation();
//          Vec4 q;
//          std::vector<double> q;
//          ceres::RotationMatrixToQuaternion(R.data(), q.data());
          Eigen::Quaterniond q(R);
          //
//          out_extrinsic << itPose->first << " " << q(0) << " "<< q(1) << " " << q(2) << " " << q(3) << " " << t(0) << " " << t(1) << " " << t(2) << std::endl;
          out_extrinsic << itPose->first;
          out_extrinsic << std::setprecision(8) << " " << q.w() << " "<< q.x() << " " << q.y() << " " << q.z() << " " << t(0) << " " << t(1) << " " << t(2) << std::endl;

          //
          Vec3 C = itPose->second.center();
          out_testPose << itPose->first;
          out_testPose << std::setprecision(8) << " " << C(0) << " " << C(1) << " " << C(2) << std::endl;


      }
      out_extrinsic.close();
      //
      for (const auto & view_it : sfm_data_.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        out_testView << prior->id_view << " " << prior->pose_center_(0) << " " << prior->pose_center_(1) << " " << prior->pose_center_(2) << std::endl;
      }
      //
      /// save intrinsics
      std::ofstream out_intrinsic;
      std::string out_intrinsic_path = sOutDir + "/intrinsics_" + sBId + ".res";
      out_intrinsic.open(out_intrinsic_path);
      if(!out_intrinsic)
      {
          std::cout << "create " << out_intrinsic_path << " file failed!" << std::endl;
      }
      for(Intrinsics::const_iterator itIntrinsic = sfm_data_.intrinsics.begin(); itIntrinsic != sfm_data_.intrinsics.end(); ++itIntrinsic)
      {
          std::vector<double> thisDatum = itIntrinsic->second->getParams();
          out_intrinsic << std::setprecision(18) << thisDatum[0] << " " << thisDatum[1] << " " << thisDatum[2] << std::endl;
      }
      out_intrinsic.close();
  }

  return true;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process_GCP_GPS_BAStep_load(std::string rtFilePath, std::string sOutDir, int bId) {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if (set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

    Hash_Map<IndexT, Mat3> global_rotations;

    //init inter-constraint
    Mat3 Rb_trans, Rf_trans, Rd_trans;
    Vec3 tb_trans, tf_trans, td_trans;
    Mat3 R_imu_trans;
    Vec3 gps_trans;
    std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
    const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
    const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};
    {
        //test distortion
        Mat3 r_tb, r_tf;
        r_tb << 1.0, 0.0, 0.0,
                0.0, 0.787926, -0.61577,
                0.0, 0.61577, 0.787926;
        r_tf << 1.0, 0.0, 0.0,
                0.0, 0.78796,  0.615727,
                0.0, -0.615727, 0.78796;

  //      Rx << 1.0, 0.0, 0.0,
  //             0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //             0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


  //      r_tb << 1.0, 0.0, 0.0,
  //              0.0, cos(38.0/180.0*M_PI), -sin(38.0/180.0*M_PI),
  //              0.0, sin(38.0/180.0*M_PI), cos(38.0/180.0*M_PI);

  //      r_tf << 1.0, 0.0, 0.0,
  //              0.0, cos(-38.0/180.0*M_PI), -sin(-38.0/180.0*M_PI),
  //              0.0, sin(-38.0/180.0*M_PI), cos(-38.0/180.0*M_PI);

  //      r_tb << cos(-38.0/180.0*M_PI), 0.0, sin(-38.0/180.0*M_PI),//-38
  //              0.0, 1.0, 0.0,
  //              -sin(-38.0/180.0*M_PI), 0.0, cos(-38.0/180.0*M_PI);
  //      r_tf << cos(38.0/180.0*M_PI), 0.0, sin(38.0/180.0*M_PI),//38
  //              0.0, 1.0, 0.0,
  //              -sin(38.0/180.0*M_PI), 0.0, cos(38.0/180.0*M_PI);
        Vec3 RA_tb, RA_tf;
        ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
        ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

        //set reference
        double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), -RA_tb(2)};
        double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), -RA_tf(2)};
  //      double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), RA_tb(2)};
  //      double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), RA_tf(2)};
        double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};

        Vec3 c_bd, c_fd;//b-d, f-d
        c_bd << 0.0, 0.098, 0.027;
        c_fd << 0.0, -0.101, 0.028;
  //      c_bd << -0.098, 0.0, -0.027;
  //      c_fd << 0.101, 0.0, -0.028;

  //      c_bd << 0.098, 0.0, -0.027;
  //      c_fd << -0.101, 0.0, -0.028;
  //      c_bd << 0.0, 1.0, 0.027;
  //      c_fd << 0.0, -0.101, 0.028;

        Vec3 t_db, t_df, t_dd;
        t_db = -r_tb*c_bd;
        t_df = -r_tf*c_fd;

        //double S = 1470.531783220221;
      //        t_db << 0.0, 0.098, 0.027;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
      //        t_df << 0.0, -0.101, 0.028;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
        t_dd << 0.0, 0.0, 0.0;

        transformation_br = {R_br_angleAxis[0], R_br_angleAxis[1], R_br_angleAxis[2],
                             t_db(0), t_db(1), t_db(2)};
        transformation_fr = {R_fr_angleAxis[0], R_fr_angleAxis[1], R_fr_angleAxis[2],
                             t_df(0), t_df(1), t_df(2)};
        transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                             t_dd(0), t_dd(1), t_dd(2)};

  //      const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
  //      const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;

        //prepare for setting pose
        ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rb_trans.data());
        tb_trans = Vec3(transformation_br[3], transformation_br[4], transformation_br[5]);
        ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rf_trans.data());
        tf_trans = Vec3(transformation_fr[3], transformation_fr[4], transformation_fr[5]);
        ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rd_trans.data());
        td_trans = Vec3(transformation_dr[3], transformation_dr[4], transformation_dr[5]);


        ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
        gps_trans = Vec3(tramsformation_x_gps[0], tramsformation_x_gps[1], tramsformation_x_gps[2]);

    }


    //get gps and imu data
    PoseStatusMap poseStatusMap;
    Hash_Map<IndexT, std::vector<double> > C_gps_Map;
    Hash_Map<IndexT, std::vector<double> > R_imu_Map;
    std::vector<double> c0(3);
    {
      //        /Hash_Map<IndexT, std::vector<double> > R_imu;
  //            Hash_Map<IndexT, std::vector<double> > C_gps;
        std::string tFilePath = rtFilePath;  //!!!
        std::ifstream in;
        in.open(tFilePath);
        if(!in)
        {
            std::cout << "open tFilePath file error!" << std::endl;
            return false;
        }
        std::size_t line = 0;
        int count = -1;
        in >> count;  ///!!!

        bool isC0 = false;//true;

        while((!in.eof()) && (line < count))
        {
            // input
            int imgId, imgStatus;
            std::vector<double> c(3);
            //double pitch, roll, yaw;
            double omega, phi, kappa;
            in >> imgId
               >> c[0] >> c[1] >> c[2]
  //             >> c[1] >> c[0] >> c[2]
               >> omega >> phi >> kappa
               >> imgStatus;

  //          omega = omega + 2.4;
  //          phi = phi + 0.0713;
  //          kappa = kappa - 0.5805;
  //          kappa = kappa + 2.4;
  //          phi = phi + 0.0713;
  //          omega = omega - 0.5805;

            if(isC0 == false)
            {
                c0 = c;
                isC0 = true;
            }

//            c[0] -= c0[0];
//            c[1] -= c0[1];
//            c[2] -= c0[2];

            ++line;



  //          double roll, pitch, yaw;

  //          roll = omega;
  //          pitch = phi;
  //          yaw = kappa;

  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
  //                sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
  //                0.0, 1.0, 0.0,
  //                -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
  //          Rx << 1.0, 0.0, 0.0,
  //                0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
  //                0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);
  //          Mat3 R_imu =  Rz*Rx*Ry;
  //          Mat3 changeAxis_M;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, -1.0, 0.0,
  //                          0.0, 0.0, -1.0;
  //          R_imu = changeAxis_M * R_imu;




  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
  //                0.0, 1.0, 0.0,
  //                -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);
  //          Rx << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

  ////          Mat3 R_imu =  Rx*Ry*Rz;
  ////          Mat3 R_imu =  Rx*Rz*Ry;//
  ////          Mat3 R_imu =  Ry*Rz*Rx;
  ////          Mat3 R_imu =  Ry*Rx*Rz;
  //          Mat3 R_imu =  Rz*Rx*Ry;//
  ////          Mat3 R_imu =  Rz*Ry*Rx;


            //use now
            Mat3 Rz, Ry, Rx;
            Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
                  sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
                  0.0, 0.0, 1.0;
            Ry << cos(phi/180.0*M_PI), 0.0, sin(phi/180.0*M_PI),
                  0.0, 1.0, 0.0,
                  -sin(phi/180.0*M_PI), 0.0, cos(phi/180.0*M_PI);
            Rx << 1.0, 0.0, 0.0,
                  0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
                  0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
  //          Mat3 R_imu =  Rz*Ry*Rx;
            Mat3 R_imu =  Rx*Ry*Rz;



  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, -cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), -cos(phi/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;

  //          Rz << cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI), 0.0,
  //                sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, -cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
  //                0.0, sin(omega/180.0*M_PI), -cos(omega/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;

  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;          Mat3 Rz, Ry, Rx;
  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(phi/180.0*M_PI), sin(phi/180.0*M_PI), 0.0,
  //                -sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
  //                0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;



  //          Mat3 R_imu =  Ry*Rz*Rx;
  //          Mat3 R_imu =  Rz*Ry*Rx;
  //          Mat3 R_imu =  Rx*Ry*Rz;


  //          Mat3 Rk, Ra, RA;
  //          Rk << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          RA << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ra << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


  //          Mat3 R_imu =  Ry*Rx*Rz;

  //          Mat3 changeAxis_M;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, -1.0, 0.0,
  //                          0.0, 0.0, -1.0;  //M1
  //          changeAxis_M << 0.0, 1.0, 0.0,
  //                          1.0, 0.0, 0.0,
  //                          0.0, 0.0, -1.0;  //M2

  //          R_imu = changeAxis_M * R_imu;
  //          R_imu = R_imu * changeAxis_M;





  //              tempGet = getListXYZ(-phi3, kappa3, -omega3);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //               tempGet = getListXYZ(phi3, kappa3, omega3);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //          tempGet = getListXYZ(phi1, kappa1, omega1);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //          Mat3 Rz, Ry, Rx;
  //          Rx << 1.0, 0.0, 0.0,
  //                 0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                 0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

  //          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
  //                 0.0, 1.0, 0.0,
  //                 -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);

  //          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                 sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                 0.0, 0.0, 1.0;


  //          Mat3 R_imu =  Ry*Rz*Rx;
  //          Mat3 R_imu =  Rz*Rx*Ry;
  //          Mat3 R_imu =  Rx*Rz*Ry;
  //          Mat3 R_imu =  Rx*Ry*Rz;

            Mat3 changeAxis_M;
            changeAxis_M << 1.0, 0.0, 0.0,
                            0.0, -1.0, 0.0,
                            0.0, 0.0, -1.0;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, 1.0, 0.0,
  //                          0.0, 0.0, -1.0;

  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, 1.0, 0.0,
  //                          0.0, 0.0, -1.0;
            Mat3 R_;
            R_ << 1.0, 0.0, 0.0,
                  0.0, 1.0, 0.0,
                  0.0, 0.0, -1.0;
            R_imu = changeAxis_M * R_imu;//* R_;


            //Mat3 R2;
            //ceres::AngleAxisToRotationMatrix(R2_a.data(), R2.data());

            //R_imu = changeAxis_M * R_imu * R2;//* R_;


            std::vector<double> R_imu_a(3);
            ceres::RotationMatrixToAngleAxis(R_imu.data(), R_imu_a.data());

            Vec3 C_gps_res = Vec3(c[0], c[1], c[2]);


            {
                bool find = false;
                for(Views::const_iterator itV = sfm_data_.views.begin(); itV != sfm_data_.views.end(); ++itV)
                {
                    if(itV->first == imgId)
                    {
                        find = true;
                        break;
                    }
                }
                if(find == false)
                    continue;
            }
            // check status
            if(imgId/100000 == 1)//back
            {
                if(imgStatus != 321)
                {
                    C_gps_Map[imgId + 100000] = c;
                    R_imu_Map[imgId + 100000] = R_imu_a;
                }

                //set back-camera pose
                Mat3 Rb_res = Rb_trans * R_imu_trans * R_imu;

                Vec3 tb_res;
                tb_res = Rb_trans * gps_trans - Rb_trans * R_imu_trans * R_imu * C_gps_res + tb_trans;
                global_rotations[imgId] = Rb_res;

                sfm_data_.poses[imgId] = Pose3(Rb_res, -Rb_res.transpose() * tb_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 1;
                poseStatusMap[imgId] = poseStatus;

            }else if(imgId/100000 == 2)//down
            {
                C_gps_Map[imgId] = c;
                R_imu_Map[imgId] = R_imu_a;

                //set down-camera pose
                Mat3 Rd_res = Rd_trans * R_imu_trans * R_imu;

                Vec3 td_res;
                td_res = Rd_trans * gps_trans - Rd_trans * R_imu_trans * R_imu * C_gps_res + td_trans;
                global_rotations[imgId] = Rd_res;

                sfm_data_.poses[imgId] = Pose3(Rd_res, -Rd_res.transpose() * td_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 2;
                poseStatusMap[imgId] = poseStatus;

            }else if(imgId/100000 == 3)//front
            {
                if(imgStatus != 321)
                {
                    C_gps_Map[imgId - 100000] = c;
                    R_imu_Map[imgId - 100000] = R_imu_a;
                }

                //set front-camera pose
                Mat3 Rf_res = Rf_trans * R_imu_trans * R_imu;

                Vec3 tf_res;
                tf_res = Rf_trans * gps_trans - Rf_trans * R_imu_trans * R_imu * C_gps_res + tf_trans;
                global_rotations[imgId] = Rf_res;

                sfm_data_.poses[imgId] = Pose3(Rf_res, -Rf_res.transpose() * tf_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 3;
                poseStatusMap[imgId] = poseStatus;

            }

        }
        in.close();

    }


    matching::PairWiseMatches tripletWise_matches;
    if (!Compute_Global_Translations_gps(global_rotations, tripletWise_matches))
    {
      std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
      return false;
    }

  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }
  {
      std::stringstream ss;
      ss << bId;
      std::string sBId;
      ss >> sBId;
      std::ofstream out_landmark;
      std::string out_landmark_path = sOutDir + "/landmarks_" + sBId + ".res";
      out_landmark.open(out_landmark_path);
      if(!out_landmark)
      {
          std::cout << "create " << out_landmark_path << " file failed!" << std::endl;
      }
//      std::ofstream out_X;
//      std::string out_X_path = sfm_data_.s_root_path + "/X_all.res";
//      out_X.open(out_X_path);
//      if(!out_X)
//      {
//          std::cout << "create " << out_X_path << " file failed!" << std::endl;
//      }
      //
      int countL = sfm_data_.structure.size();
      out_landmark << countL << std::endl;
      //
      for(Landmarks::const_iterator itL = sfm_data_.structure.begin(); itL != sfm_data_.structure.end(); ++itL)
      {
          int Xid = itL->first;
          Vec3 X = itL->second.X;
          out_landmark << Xid << " " << X(0) << " " << X(1) << " " << X(2) << std::endl;
          //
          int countObs = itL->second.obs.size();
          out_landmark << countObs << " ";
          for(Observations::const_iterator itObs = itL->second.obs.begin(); itObs != itL->second.obs.end(); ++itObs)
          {
              int imgId = itObs->first;
              Vec2 featId = itObs->second.x;
              out_landmark << imgId << " " << featId(0) << " " << featId(1) << " ";
          }
          out_landmark << std::endl;
      }
      out_landmark.close();
      //
      /// save extrinsics
      std::ofstream out_extrinsic;
      std::string out_extrinsic_path = sOutDir + "/extrinsics_" + sBId + ".res";
      out_extrinsic.open(out_extrinsic_path);
      if(!out_extrinsic)
      {
          std::cout << "create " << out_extrinsic_path << " file failed!" << std::endl;
      }
      out_extrinsic << sfm_data_.poses.size() << std::endl;
      for(Poses::const_iterator itPose = sfm_data_.poses.begin(); itPose != sfm_data_.poses.end(); ++itPose)
      {
          Vec3 t = itPose->second.translation();
          std::cout << "t" << std::endl;
          Mat3 R = itPose->second.rotation();
          std::cout << "R" << std::endl;
          Vec4 q;
//          std::vector<double> q;
          ceres::RotationMatrixToQuaternion(R.data(), q.data());
          std::cout << "q" << std::endl;
          //
//          out_extrinsic << itPose->first << " " << q[0] << " "<< q[1] << " " << q[2] << " " << q[3] << " " << t(0) << " " << t(1) << " " << t(2) << std::endl;
          out_extrinsic << itPose->first << " " << q(0) << " "<< q(1) << " " << q(2) << " " << q(3) << " " << t(0) << " " << t(1) << " " << t(2) << std::endl;
      }
      out_extrinsic.close();
  }

  return true;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Process_GCP_GPS_water(std::string rtFilePath) {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if (set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

    Hash_Map<IndexT, Mat3> global_rotations;

    //init inter-constraint
    Mat3 Rb_trans, Rf_trans, Rd_trans;
    Vec3 tb_trans, tf_trans, td_trans;
    Mat3 R_imu_trans;
    Vec3 gps_trans;
    std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
    const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
    const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};
    {
        //test distortion
        Mat3 r_tb, r_tf;
        r_tb << 1.0, 0.0, 0.0,
                0.0, 0.787926, -0.61577,
                0.0, 0.61577, 0.787926;
        r_tf << 1.0, 0.0, 0.0,
                0.0, 0.78796,  0.615727,
                0.0, -0.615727, 0.78796;

  //      Rx << 1.0, 0.0, 0.0,
  //             0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //             0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


  //      r_tb << 1.0, 0.0, 0.0,
  //              0.0, cos(38.0/180.0*M_PI), -sin(38.0/180.0*M_PI),
  //              0.0, sin(38.0/180.0*M_PI), cos(38.0/180.0*M_PI);

  //      r_tf << 1.0, 0.0, 0.0,
  //              0.0, cos(-38.0/180.0*M_PI), -sin(-38.0/180.0*M_PI),
  //              0.0, sin(-38.0/180.0*M_PI), cos(-38.0/180.0*M_PI);

  //      r_tb << cos(-38.0/180.0*M_PI), 0.0, sin(-38.0/180.0*M_PI),//-38
  //              0.0, 1.0, 0.0,
  //              -sin(-38.0/180.0*M_PI), 0.0, cos(-38.0/180.0*M_PI);
  //      r_tf << cos(38.0/180.0*M_PI), 0.0, sin(38.0/180.0*M_PI),//38
  //              0.0, 1.0, 0.0,
  //              -sin(38.0/180.0*M_PI), 0.0, cos(38.0/180.0*M_PI);
        Vec3 RA_tb, RA_tf;
        ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
        ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

        //set reference
        double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), -RA_tb(2)};
        double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), -RA_tf(2)};
  //      double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), RA_tb(2)};
  //      double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), RA_tf(2)};
        double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};

        Vec3 c_bd, c_fd;//b-d, f-d
        c_bd << 0.0, 0.098, 0.027;
        c_fd << 0.0, -0.101, 0.028;
  //      c_bd << -0.098, 0.0, -0.027;
  //      c_fd << 0.101, 0.0, -0.028;

  //      c_bd << 0.098, 0.0, -0.027;
  //      c_fd << -0.101, 0.0, -0.028;
  //      c_bd << 0.0, 1.0, 0.027;
  //      c_fd << 0.0, -0.101, 0.028;

        Vec3 t_db, t_df, t_dd;
        t_db = -r_tb*c_bd;
        t_df = -r_tf*c_fd;

        //double S = 1470.531783220221;
      //        t_db << 0.0, 0.098, 0.027;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
      //        t_df << 0.0, -0.101, 0.028;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
        t_dd << 0.0, 0.0, 0.0;

        transformation_br = {R_br_angleAxis[0], R_br_angleAxis[1], R_br_angleAxis[2],
                             t_db(0), t_db(1), t_db(2)};
        transformation_fr = {R_fr_angleAxis[0], R_fr_angleAxis[1], R_fr_angleAxis[2],
                             t_df(0), t_df(1), t_df(2)};
        transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                             t_dd(0), t_dd(1), t_dd(2)};

  //      const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
  //      const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;

        //prepare for setting pose
        ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rb_trans.data());
        tb_trans = Vec3(transformation_br[3], transformation_br[4], transformation_br[5]);
        ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rf_trans.data());
        tf_trans = Vec3(transformation_fr[3], transformation_fr[4], transformation_fr[5]);
        ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rd_trans.data());
        td_trans = Vec3(transformation_dr[3], transformation_dr[4], transformation_dr[5]);


        ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
        gps_trans = Vec3(tramsformation_x_gps[0], tramsformation_x_gps[1], tramsformation_x_gps[2]);

    }


    //get gps and imu data
    PoseStatusMap poseStatusMap;
    Hash_Map<IndexT, std::vector<double> > C_gps_Map;
    Hash_Map<IndexT, std::vector<double> > R_imu_Map;
    std::vector<double> c0(3);
    {
      //        /Hash_Map<IndexT, std::vector<double> > R_imu;
  //            Hash_Map<IndexT, std::vector<double> > C_gps;
        std::string tFilePath = rtFilePath;  //!!!
        std::ifstream in;
        in.open(tFilePath);
        if(!in)
        {
            std::cout << "open tFilePath file error!" << std::endl;
            return false;
        }
        std::size_t line = 0;
        int count = -1;
        in >> count;  ///!!!

        bool isC0 = false;//true;

        while((!in.eof()) && (line < count))
        {
            // input
            int imgId, imgStatus;
            std::vector<double> c(3);
            //double pitch, roll, yaw;
            double omega, phi, kappa;
            in >> imgId
               >> c[0] >> c[1] >> c[2]
  //             >> c[1] >> c[0] >> c[2]
               >> omega >> phi >> kappa
               >> imgStatus;

  //          omega = omega + 2.4;
  //          phi = phi + 0.0713;
  //          kappa = kappa - 0.5805;
  //          kappa = kappa + 2.4;
  //          phi = phi + 0.0713;
  //          omega = omega - 0.5805;

            if(isC0 == false)
            {
                c0 = c;
                isC0 = true;
            }

//            c[0] -= c0[0];
//            c[1] -= c0[1];
//            c[2] -= c0[2];

            ++line;



  //          double roll, pitch, yaw;

  //          roll = omega;
  //          pitch = phi;
  //          yaw = kappa;

  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
  //                sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
  //                0.0, 1.0, 0.0,
  //                -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
  //          Rx << 1.0, 0.0, 0.0,
  //                0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
  //                0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);
  //          Mat3 R_imu =  Rz*Rx*Ry;
  //          Mat3 changeAxis_M;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, -1.0, 0.0,
  //                          0.0, 0.0, -1.0;
  //          R_imu = changeAxis_M * R_imu;




  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
  //                0.0, 1.0, 0.0,
  //                -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);
  //          Rx << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

  ////          Mat3 R_imu =  Rx*Ry*Rz;
  ////          Mat3 R_imu =  Rx*Rz*Ry;//
  ////          Mat3 R_imu =  Ry*Rz*Rx;
  ////          Mat3 R_imu =  Ry*Rx*Rz;
  //          Mat3 R_imu =  Rz*Rx*Ry;//
  ////          Mat3 R_imu =  Rz*Ry*Rx;


            //use now
            Mat3 Rz, Ry, Rx;
            Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
                  sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
                  0.0, 0.0, 1.0;
            Ry << cos(phi/180.0*M_PI), 0.0, sin(phi/180.0*M_PI),
                  0.0, 1.0, 0.0,
                  -sin(phi/180.0*M_PI), 0.0, cos(phi/180.0*M_PI);
            Rx << 1.0, 0.0, 0.0,
                  0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
                  0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
  //          Mat3 R_imu =  Rz*Ry*Rx;
            Mat3 R_imu =  Rx*Ry*Rz;



  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, -cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), -cos(phi/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;

  //          Rz << cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI), 0.0,
  //                sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, -cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
  //                0.0, sin(omega/180.0*M_PI), -cos(omega/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;

  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;          Mat3 Rz, Ry, Rx;
  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(phi/180.0*M_PI), sin(phi/180.0*M_PI), 0.0,
  //                -sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
  //                0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;



  //          Mat3 R_imu =  Ry*Rz*Rx;
  //          Mat3 R_imu =  Rz*Ry*Rx;
  //          Mat3 R_imu =  Rx*Ry*Rz;


  //          Mat3 Rk, Ra, RA;
  //          Rk << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          RA << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ra << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


  //          Mat3 R_imu =  Ry*Rx*Rz;

  //          Mat3 changeAxis_M;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, -1.0, 0.0,
  //                          0.0, 0.0, -1.0;  //M1
  //          changeAxis_M << 0.0, 1.0, 0.0,
  //                          1.0, 0.0, 0.0,
  //                          0.0, 0.0, -1.0;  //M2

  //          R_imu = changeAxis_M * R_imu;
  //          R_imu = R_imu * changeAxis_M;





  //              tempGet = getListXYZ(-phi3, kappa3, -omega3);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //               tempGet = getListXYZ(phi3, kappa3, omega3);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //          tempGet = getListXYZ(phi1, kappa1, omega1);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //          Mat3 Rz, Ry, Rx;
  //          Rx << 1.0, 0.0, 0.0,
  //                 0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                 0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

  //          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
  //                 0.0, 1.0, 0.0,
  //                 -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);

  //          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                 sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                 0.0, 0.0, 1.0;


  //          Mat3 R_imu =  Ry*Rz*Rx;
  //          Mat3 R_imu =  Rz*Rx*Ry;
  //          Mat3 R_imu =  Rx*Rz*Ry;
  //          Mat3 R_imu =  Rx*Ry*Rz;

            Mat3 changeAxis_M;
            changeAxis_M << 1.0, 0.0, 0.0,
                            0.0, -1.0, 0.0,
                            0.0, 0.0, -1.0;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, 1.0, 0.0,
  //                          0.0, 0.0, -1.0;

  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, 1.0, 0.0,
  //                          0.0, 0.0, -1.0;
            Mat3 R_;
            R_ << 1.0, 0.0, 0.0,
                  0.0, 1.0, 0.0,
                  0.0, 0.0, -1.0;
            R_imu = changeAxis_M * R_imu;//* R_;


            //Mat3 R2;
            //ceres::AngleAxisToRotationMatrix(R2_a.data(), R2.data());

            //R_imu = changeAxis_M * R_imu * R2;//* R_;


            std::vector<double> R_imu_a(3);
            ceres::RotationMatrixToAngleAxis(R_imu.data(), R_imu_a.data());

            Vec3 C_gps_res = Vec3(c[0], c[1], c[2]);


            {
                bool find = false;
                for(Views::const_iterator itV = sfm_data_.views.begin(); itV != sfm_data_.views.end(); ++itV)
                {
                    if(itV->first == imgId)
                    {
                        find = true;
                        break;
                    }
                }
                if(find == false)
                    continue;
            }
            // check status
            if(imgId/100000 == 1)//back
            {
                if(imgStatus != 321)
                {
                    C_gps_Map[imgId + 100000] = c;
                    R_imu_Map[imgId + 100000] = R_imu_a;
                }

                //set back-camera pose
                Mat3 Rb_res = Rb_trans * R_imu_trans * R_imu;

                Vec3 tb_res;
                tb_res = Rb_trans * gps_trans - Rb_trans * R_imu_trans * R_imu * C_gps_res + tb_trans;
                global_rotations[imgId] = Rb_res;

                sfm_data_.poses[imgId] = Pose3(Rb_res, -Rb_res.transpose() * tb_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 1;
                poseStatusMap[imgId] = poseStatus;

            }else if(imgId/100000 == 2)//down
            {
                C_gps_Map[imgId] = c;
                R_imu_Map[imgId] = R_imu_a;

                //set down-camera pose
                Mat3 Rd_res = Rd_trans * R_imu_trans * R_imu;

                Vec3 td_res;
                td_res = Rd_trans * gps_trans - Rd_trans * R_imu_trans * R_imu * C_gps_res + td_trans;
                global_rotations[imgId] = Rd_res;

                sfm_data_.poses[imgId] = Pose3(Rd_res, -Rd_res.transpose() * td_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 2;
                poseStatusMap[imgId] = poseStatus;

            }else if(imgId/100000 == 3)//front
            {
                if(imgStatus != 321)
                {
                    C_gps_Map[imgId - 100000] = c;
                    R_imu_Map[imgId - 100000] = R_imu_a;
                }

                //set front-camera pose
                Mat3 Rf_res = Rf_trans * R_imu_trans * R_imu;

                Vec3 tf_res;
                tf_res = Rf_trans * gps_trans - Rf_trans * R_imu_trans * R_imu * C_gps_res + tf_trans;
                global_rotations[imgId] = Rf_res;

                sfm_data_.poses[imgId] = Pose3(Rf_res, -Rf_res.transpose() * tf_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 3;
                poseStatusMap[imgId] = poseStatus;

            }

        }
        in.close();

    }


    matching::PairWiseMatches tripletWise_matches;
    if (!Compute_Global_Translations_gps(global_rotations, tripletWise_matches))
    {
      std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
      return false;
    }

  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }
  ///save for g2o
  {
//      std::ofstream out_landmark;
//      std::string out_landmark_path = sfm_data_.s_root_path + "/landmarks.res";
//      out_landmark.open(out_landmark_path);
//      if(!out_landmark)
//      {
//          std::cout << "create " << out_landmark_path << " file failed!" << std::endl;
//      }
////      std::ofstream out_X;
////      std::string out_X_path = sfm_data_.s_root_path + "/X_all.res";
////      out_X.open(out_X_path);
////      if(!out_X)
////      {
////          std::cout << "create " << out_X_path << " file failed!" << std::endl;
////      }
//      //
//      int countL = sfm_data_.structure.size();
//      out_landmark << countL << std::endl;
//      //
//      for(Landmarks::const_iterator itL = sfm_data_.structure.begin(); itL != sfm_data_.structure.end(); ++itL)
//      {
//          int Xid = itL->first;
//          Vec3 X = itL->second.X;
//          out_landmark << Xid << " " << X(0) << " " << X(1) << " " << X(2) << std::endl;
//          //
//          int countObs = itL->second.obs.size();
//          out_landmark << countObs << " ";
//          for(Observations::const_iterator itObs = itL->second.obs.begin(); itObs != itL->second.obs.end(); ++itObs)
//          {
//              int imgId = itObs->first;
//              Vec2 featId = itObs->second.x;
//              out_landmark << imgId << " " << featId(0) << " " << featId(1) << " ";
//          }
//          out_landmark << std::endl;
//      }
//      out_landmark.close();
//      //
//      /// save extrinsics
//      std::ofstream out_extrinsic;
//      std::string out_extrinsic_path = sfm_data_.s_root_path + "/extrinsics.res";
//      out_extrinsic.open(out_extrinsic_path);
//      if(!out_extrinsic)
//      {
//          std::cout << "create " << out_extrinsic_path << " file failed!" << std::endl;
//      }
//      out_extrinsic << sfm_data_.poses.size() << std::endl;
//      for(Poses::const_iterator itPose = sfm_data_.poses.begin(); itPose != sfm_data_.poses.end(); ++itPose)
//      {
//          Vec3 t = itPose->second.translation();
//          std::cout << "t" << std::endl;
//          Mat3 R = itPose->second.rotation();
//          std::cout << "R" << std::endl;
//          Vec4 q;
////          std::vector<double> q;
//          ceres::RotationMatrixToQuaternion(R.data(), q.data());
//          std::cout << "q" << std::endl;
//          //
////          out_extrinsic << itPose->first << " " << q[0] << " "<< q[1] << " " << q[2] << " " << q[3] << " " << t(0) << " " << t(1) << " " << t(2) << std::endl;
//          out_extrinsic << itPose->first << " " << q(0) << " "<< q(1) << " " << q(2) << " " << q(3) << " " << t(0) << " " << t(1) << " " << t(2) << std::endl;
//      }
//      out_extrinsic.close();



  }
  if (!Adjust_init_water_c1())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Process_GCP_GPS_water2(std::string rtFilePath) {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
//  {
//    const Pair_Set pairs = matches_provider_->getPairs();
//    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
//    if (set_remainingIds.empty())
//    {
//      std::cout << "Invalid input image graph for global SfM" << std::endl;
//      return false;
//    }
//    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
//  }

//    Hash_Map<IndexT, Mat3> global_rotations;

//    //init inter-constraint
//    Mat3 Rb_trans, Rf_trans, Rd_trans;
//    Vec3 tb_trans, tf_trans, td_trans;
//    Mat3 R_imu_trans;
//    Vec3 gps_trans;
//    std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
//    const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
//    const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};
//    {
//        //test distortion
//        Mat3 r_tb, r_tf;
//        r_tb << 1.0, 0.0, 0.0,
//                0.0, 0.787926, -0.61577,
//                0.0, 0.61577, 0.787926;
//        r_tf << 1.0, 0.0, 0.0,
//                0.0, 0.78796,  0.615727,
//                0.0, -0.615727, 0.78796;

//  //      Rx << 1.0, 0.0, 0.0,
//  //             0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//  //             0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


//  //      r_tb << 1.0, 0.0, 0.0,
//  //              0.0, cos(38.0/180.0*M_PI), -sin(38.0/180.0*M_PI),
//  //              0.0, sin(38.0/180.0*M_PI), cos(38.0/180.0*M_PI);

//  //      r_tf << 1.0, 0.0, 0.0,
//  //              0.0, cos(-38.0/180.0*M_PI), -sin(-38.0/180.0*M_PI),
//  //              0.0, sin(-38.0/180.0*M_PI), cos(-38.0/180.0*M_PI);

//  //      r_tb << cos(-38.0/180.0*M_PI), 0.0, sin(-38.0/180.0*M_PI),//-38
//  //              0.0, 1.0, 0.0,
//  //              -sin(-38.0/180.0*M_PI), 0.0, cos(-38.0/180.0*M_PI);
//  //      r_tf << cos(38.0/180.0*M_PI), 0.0, sin(38.0/180.0*M_PI),//38
//  //              0.0, 1.0, 0.0,
//  //              -sin(38.0/180.0*M_PI), 0.0, cos(38.0/180.0*M_PI);
//        Vec3 RA_tb, RA_tf;
//        ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
//        ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

//        //set reference
//        double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), -RA_tb(2)};
//        double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), -RA_tf(2)};
//  //      double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), RA_tb(2)};
//  //      double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), RA_tf(2)};
//        double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};

//        Vec3 c_bd, c_fd;//b-d, f-d
//        c_bd << 0.0, 0.098, 0.027;
//        c_fd << 0.0, -0.101, 0.028;
//  //      c_bd << -0.098, 0.0, -0.027;
//  //      c_fd << 0.101, 0.0, -0.028;

//  //      c_bd << 0.098, 0.0, -0.027;
//  //      c_fd << -0.101, 0.0, -0.028;
//  //      c_bd << 0.0, 1.0, 0.027;
//  //      c_fd << 0.0, -0.101, 0.028;

//        Vec3 t_db, t_df, t_dd;
//        t_db = -r_tb*c_bd;
//        t_df = -r_tf*c_fd;

//        //double S = 1470.531783220221;
//      //        t_db << 0.0, 0.098, 0.027;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
//      //        t_df << 0.0, -0.101, 0.028;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
//        t_dd << 0.0, 0.0, 0.0;

//        transformation_br = {R_br_angleAxis[0], R_br_angleAxis[1], R_br_angleAxis[2],
//                             t_db(0), t_db(1), t_db(2)};
//        transformation_fr = {R_fr_angleAxis[0], R_fr_angleAxis[1], R_fr_angleAxis[2],
//                             t_df(0), t_df(1), t_df(2)};
//        transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
//                             t_dd(0), t_dd(1), t_dd(2)};

//  //      const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
//  //      const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};

//        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
//                                               << tramsformation_x_gps[1] << " "
//                                               << tramsformation_x_gps[2] << std::endl;

//        //prepare for setting pose
//        ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rb_trans.data());
//        tb_trans = Vec3(transformation_br[3], transformation_br[4], transformation_br[5]);
//        ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rf_trans.data());
//        tf_trans = Vec3(transformation_fr[3], transformation_fr[4], transformation_fr[5]);
//        ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rd_trans.data());
//        td_trans = Vec3(transformation_dr[3], transformation_dr[4], transformation_dr[5]);


//        ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
//        gps_trans = Vec3(tramsformation_x_gps[0], tramsformation_x_gps[1], tramsformation_x_gps[2]);

//    }


//    //get gps and imu data
//    PoseStatusMap poseStatusMap;
//    Hash_Map<IndexT, std::vector<double> > C_gps_Map;
//    Hash_Map<IndexT, std::vector<double> > R_imu_Map;
//    std::vector<double> c0(3);
//    {
//      //        /Hash_Map<IndexT, std::vector<double> > R_imu;
//  //            Hash_Map<IndexT, std::vector<double> > C_gps;
//        std::string tFilePath = rtFilePath;  //!!!
//        std::ifstream in;
//        in.open(tFilePath);
//        if(!in)
//        {
//            std::cout << "open tFilePath file error!" << std::endl;
//            return false;
//        }
//        std::size_t line = 0;
//        int count = -1;
//        in >> count;  ///!!!

//        bool isC0 = false;//true;

//        while((!in.eof()) && (line < count))
//        {
//            // input
//            int imgId, imgStatus;
//            std::vector<double> c(3);
//            //double pitch, roll, yaw;
//            double omega, phi, kappa;
//            in >> imgId
//               >> c[0] >> c[1] >> c[2]
//  //             >> c[1] >> c[0] >> c[2]
//               >> omega >> phi >> kappa
//               >> imgStatus;

//  //          omega = omega + 2.4;
//  //          phi = phi + 0.0713;
//  //          kappa = kappa - 0.5805;
//  //          kappa = kappa + 2.4;
//  //          phi = phi + 0.0713;
//  //          omega = omega - 0.5805;

//            if(isC0 == false)
//            {
//                c0 = c;
//                isC0 = true;
//            }

//            c[0] -= c0[0];
//            c[1] -= c0[1];
////            c[2] -= c0[2];

//            ++line;



//  //          double roll, pitch, yaw;

//  //          roll = omega;
//  //          pitch = phi;
//  //          yaw = kappa;

//  //          Mat3 Rz, Ry, Rx;
//  //          Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
//  //                sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
//  //                0.0, 0.0, 1.0;
//  //          Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
//  //                0.0, 1.0, 0.0,
//  //                -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
//  //          Rx << 1.0, 0.0, 0.0,
//  //                0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
//  //                0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);
//  //          Mat3 R_imu =  Rz*Rx*Ry;
//  //          Mat3 changeAxis_M;
//  //          changeAxis_M << 1.0, 0.0, 0.0,
//  //                          0.0, -1.0, 0.0,
//  //                          0.0, 0.0, -1.0;
//  //          R_imu = changeAxis_M * R_imu;




//  //          Mat3 Rz, Ry, Rx;
//  //          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//  //                0.0, 0.0, 1.0;
//  //          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
//  //                0.0, 1.0, 0.0,
//  //                -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);
//  //          Rx << 1.0, 0.0, 0.0,
//  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

//  ////          Mat3 R_imu =  Rx*Ry*Rz;
//  ////          Mat3 R_imu =  Rx*Rz*Ry;//
//  ////          Mat3 R_imu =  Ry*Rz*Rx;
//  ////          Mat3 R_imu =  Ry*Rx*Rz;
//  //          Mat3 R_imu =  Rz*Rx*Ry;//
//  ////          Mat3 R_imu =  Rz*Ry*Rx;


//            //use now
//            Mat3 Rz, Ry, Rx;
//            Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                  sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                  0.0, 0.0, 1.0;
//            Ry << cos(phi/180.0*M_PI), 0.0, sin(phi/180.0*M_PI),
//                  0.0, 1.0, 0.0,
//                  -sin(phi/180.0*M_PI), 0.0, cos(phi/180.0*M_PI);
//            Rx << 1.0, 0.0, 0.0,
//                  0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
//                  0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
//  //          Mat3 R_imu =  Rz*Ry*Rx;
//            Mat3 R_imu =  Rx*Ry*Rz;



//  //          Mat3 Rz, Ry, Rx;
//  //          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
//  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
//  //                0.0, 0.0, 1.0;
//  //          Ry << 1.0, 0.0, 0.0,
//  //                0.0, -cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//  //                0.0, sin(phi/180.0*M_PI), -cos(phi/180.0*M_PI);
//  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//  //                0.0, 0.0, 1.0;

//  //          Rz << cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI), 0.0,
//  //                sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
//  //                0.0, 0.0, 1.0;
//  //          Ry << 1.0, 0.0, 0.0,
//  //                0.0, -cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
//  //                0.0, sin(omega/180.0*M_PI), -cos(omega/180.0*M_PI);
//  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//  //                0.0, 0.0, 1.0;

//  //          Mat3 Rz, Ry, Rx;
//  //          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
//  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
//  //                0.0, 0.0, 1.0;
//  //          Ry << 1.0, 0.0, 0.0,
//  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);
//  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//  //                0.0, 0.0, 1.0;          Mat3 Rz, Ry, Rx;
//  //          Mat3 Rz, Ry, Rx;
//  //          Rz << cos(phi/180.0*M_PI), sin(phi/180.0*M_PI), 0.0,
//  //                -sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
//  //                0.0, 0.0, 1.0;
//  //          Ry << 1.0, 0.0, 0.0,
//  //                0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
//  //                0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
//  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//  //                0.0, 0.0, 1.0;



//  //          Mat3 R_imu =  Ry*Rz*Rx;
//  //          Mat3 R_imu =  Rz*Ry*Rx;
//  //          Mat3 R_imu =  Rx*Ry*Rz;


//  //          Mat3 Rk, Ra, RA;
//  //          Rk << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//  //                0.0, 0.0, 1.0;
//  //          RA << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
//  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
//  //                0.0, 0.0, 1.0;
//  //          Ra << 1.0, 0.0, 0.0,
//  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


//  //          Mat3 R_imu =  Ry*Rx*Rz;

//  //          Mat3 changeAxis_M;
//  //          changeAxis_M << 1.0, 0.0, 0.0,
//  //                          0.0, -1.0, 0.0,
//  //                          0.0, 0.0, -1.0;  //M1
//  //          changeAxis_M << 0.0, 1.0, 0.0,
//  //                          1.0, 0.0, 0.0,
//  //                          0.0, 0.0, -1.0;  //M2

//  //          R_imu = changeAxis_M * R_imu;
//  //          R_imu = R_imu * changeAxis_M;





//  //              tempGet = getListXYZ(-phi3, kappa3, -omega3);
//  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

//  //               tempGet = getListXYZ(phi3, kappa3, omega3);
//  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

//  //          tempGet = getListXYZ(phi1, kappa1, omega1);
//  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

//  //          Mat3 Rz, Ry, Rx;
//  //          Rx << 1.0, 0.0, 0.0,
//  //                 0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//  //                 0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

//  //          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
//  //                 0.0, 1.0, 0.0,
//  //                 -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);

//  //          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//  //                 sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//  //                 0.0, 0.0, 1.0;


//  //          Mat3 R_imu =  Ry*Rz*Rx;
//  //          Mat3 R_imu =  Rz*Rx*Ry;
//  //          Mat3 R_imu =  Rx*Rz*Ry;
//  //          Mat3 R_imu =  Rx*Ry*Rz;

//            Mat3 changeAxis_M;
//            changeAxis_M << 1.0, 0.0, 0.0,
//                            0.0, -1.0, 0.0,
//                            0.0, 0.0, -1.0;
//  //          changeAxis_M << 1.0, 0.0, 0.0,
//  //                          0.0, 1.0, 0.0,
//  //                          0.0, 0.0, -1.0;

//  //          changeAxis_M << 1.0, 0.0, 0.0,
//  //                          0.0, 1.0, 0.0,
//  //                          0.0, 0.0, -1.0;
//            Mat3 R_;
//            R_ << 1.0, 0.0, 0.0,
//                  0.0, 1.0, 0.0,
//                  0.0, 0.0, -1.0;
//            R_imu = changeAxis_M * R_imu;//* R_;


//            //Mat3 R2;
//            //ceres::AngleAxisToRotationMatrix(R2_a.data(), R2.data());

//            //R_imu = changeAxis_M * R_imu * R2;//* R_;


//            std::vector<double> R_imu_a(3);
//            ceres::RotationMatrixToAngleAxis(R_imu.data(), R_imu_a.data());

//            Vec3 C_gps_res = Vec3(c[0], c[1], c[2]);


//            {
//                bool find = false;
//                for(Views::const_iterator itV = sfm_data_.views.begin(); itV != sfm_data_.views.end(); ++itV)
//                {
//                    if(itV->first == imgId)
//                    {
//                        find = true;
//                        break;
//                    }
//                }
//                if(find == false)
//                    continue;
//            }
//            // check status
//            if(imgId/100000 == 1)//back
//            {
//                if(imgStatus != 321)
//                {
//                    C_gps_Map[imgId + 100000] = c;
//                    R_imu_Map[imgId + 100000] = R_imu_a;
//                }

//                //set back-camera pose
//                Mat3 Rb_res = Rb_trans * R_imu_trans * R_imu;

//                Vec3 tb_res;
//                tb_res = Rb_trans * gps_trans - Rb_trans * R_imu_trans * R_imu * C_gps_res + tb_trans;
//                global_rotations[imgId] = Rb_res;

//                sfm_data_.poses[imgId] = Pose3(Rb_res, -Rb_res.transpose() * tb_res);

//                PoseStatus poseStatus;
//                poseStatus.pose = sfm_data_.poses[imgId];
//                poseStatus.status = imgStatus;
//                poseStatus.bdf = 1;
//                poseStatusMap[imgId] = poseStatus;

//            }else if(imgId/100000 == 2)//down
//            {
//                C_gps_Map[imgId] = c;
//                R_imu_Map[imgId] = R_imu_a;

//                //set down-camera pose
//                Mat3 Rd_res = Rd_trans * R_imu_trans * R_imu;

//                Vec3 td_res;
//                td_res = Rd_trans * gps_trans - Rd_trans * R_imu_trans * R_imu * C_gps_res + td_trans;
//                global_rotations[imgId] = Rd_res;

//                sfm_data_.poses[imgId] = Pose3(Rd_res, -Rd_res.transpose() * td_res);

//                PoseStatus poseStatus;
//                poseStatus.pose = sfm_data_.poses[imgId];
//                poseStatus.status = imgStatus;
//                poseStatus.bdf = 2;
//                poseStatusMap[imgId] = poseStatus;

//            }else if(imgId/100000 == 3)//front
//            {
//                if(imgStatus != 321)
//                {
//                    C_gps_Map[imgId - 100000] = c;
//                    R_imu_Map[imgId - 100000] = R_imu_a;
//                }

//                //set front-camera pose
//                Mat3 Rf_res = Rf_trans * R_imu_trans * R_imu;

//                Vec3 tf_res;
//                tf_res = Rf_trans * gps_trans - Rf_trans * R_imu_trans * R_imu * C_gps_res + tf_trans;
//                global_rotations[imgId] = Rf_res;

//                sfm_data_.poses[imgId] = Pose3(Rf_res, -Rf_res.transpose() * tf_res);

//                PoseStatus poseStatus;
//                poseStatus.pose = sfm_data_.poses[imgId];
//                poseStatus.status = imgStatus;
//                poseStatus.bdf = 3;
//                poseStatusMap[imgId] = poseStatus;

//            }

//        }
//        in.close();

//    }


//    matching::PairWiseMatches tripletWise_matches;
//    if (!Compute_Global_Translations_gps(global_rotations, tripletWise_matches))
//    {
//      std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
//      return false;
//    }

//  if (!Compute_Initial_Structure(tripletWise_matches))
//  {
//    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
//    return false;
//  }
  if (!Adjust())//_init_water2_c1())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process_GCP_GPS_water3(std::string rtFilePath) {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if (set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

    Hash_Map<IndexT, Mat3> global_rotations;

    //init inter-constraint
    Mat3 Rb_trans, Rf_trans, Rd_trans;
    Vec3 tb_trans, tf_trans, td_trans;
    Mat3 R_imu_trans;
    Vec3 gps_trans;
    std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
    const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
    const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};
    {
        //test distortion
        Mat3 r_tb, r_tf;
        r_tb << 1.0, 0.0, 0.0,
                0.0, 0.787926, -0.61577,
                0.0, 0.61577, 0.787926;
        r_tf << 1.0, 0.0, 0.0,
                0.0, 0.78796,  0.615727,
                0.0, -0.615727, 0.78796;

  //      Rx << 1.0, 0.0, 0.0,
  //             0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //             0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


  //      r_tb << 1.0, 0.0, 0.0,
  //              0.0, cos(38.0/180.0*M_PI), -sin(38.0/180.0*M_PI),
  //              0.0, sin(38.0/180.0*M_PI), cos(38.0/180.0*M_PI);

  //      r_tf << 1.0, 0.0, 0.0,
  //              0.0, cos(-38.0/180.0*M_PI), -sin(-38.0/180.0*M_PI),
  //              0.0, sin(-38.0/180.0*M_PI), cos(-38.0/180.0*M_PI);

  //      r_tb << cos(-38.0/180.0*M_PI), 0.0, sin(-38.0/180.0*M_PI),//-38
  //              0.0, 1.0, 0.0,
  //              -sin(-38.0/180.0*M_PI), 0.0, cos(-38.0/180.0*M_PI);
  //      r_tf << cos(38.0/180.0*M_PI), 0.0, sin(38.0/180.0*M_PI),//38
  //              0.0, 1.0, 0.0,
  //              -sin(38.0/180.0*M_PI), 0.0, cos(38.0/180.0*M_PI);
        Vec3 RA_tb, RA_tf;
        ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
        ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

        //set reference
        double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), -RA_tb(2)};
        double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), -RA_tf(2)};
  //      double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), RA_tb(2)};
  //      double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), RA_tf(2)};
        double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};

        Vec3 c_bd, c_fd;//b-d, f-d
        c_bd << 0.0, 0.098, 0.027;
        c_fd << 0.0, -0.101, 0.028;
  //      c_bd << -0.098, 0.0, -0.027;
  //      c_fd << 0.101, 0.0, -0.028;

  //      c_bd << 0.098, 0.0, -0.027;
  //      c_fd << -0.101, 0.0, -0.028;
  //      c_bd << 0.0, 1.0, 0.027;
  //      c_fd << 0.0, -0.101, 0.028;

        Vec3 t_db, t_df, t_dd;
        t_db = -r_tb*c_bd;
        t_df = -r_tf*c_fd;

        //double S = 1470.531783220221;
      //        t_db << 0.0, 0.098, 0.027;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
      //        t_df << 0.0, -0.101, 0.028;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
        t_dd << 0.0, 0.0, 0.0;

        transformation_br = {R_br_angleAxis[0], R_br_angleAxis[1], R_br_angleAxis[2],
                             t_db(0), t_db(1), t_db(2)};
        transformation_fr = {R_fr_angleAxis[0], R_fr_angleAxis[1], R_fr_angleAxis[2],
                             t_df(0), t_df(1), t_df(2)};
        transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                             t_dd(0), t_dd(1), t_dd(2)};

  //      const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
  //      const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;

        //prepare for setting pose
        ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rb_trans.data());
        tb_trans = Vec3(transformation_br[3], transformation_br[4], transformation_br[5]);
        ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rf_trans.data());
        tf_trans = Vec3(transformation_fr[3], transformation_fr[4], transformation_fr[5]);
        ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rd_trans.data());
        td_trans = Vec3(transformation_dr[3], transformation_dr[4], transformation_dr[5]);


        ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
        gps_trans = Vec3(tramsformation_x_gps[0], tramsformation_x_gps[1], tramsformation_x_gps[2]);

    }


    //get gps and imu data
    PoseStatusMap poseStatusMap;
    Hash_Map<IndexT, std::vector<double> > C_gps_Map;
    Hash_Map<IndexT, std::vector<double> > R_imu_Map;
    std::vector<double> c0(3);
    {
      //        /Hash_Map<IndexT, std::vector<double> > R_imu;
  //            Hash_Map<IndexT, std::vector<double> > C_gps;
        std::string tFilePath = rtFilePath;  //!!!
        std::ifstream in;
        in.open(tFilePath);
        if(!in)
        {
            std::cout << "open tFilePath file error!" << std::endl;
            return false;
        }
        std::size_t line = 0;
        int count = -1;
        in >> count;  ///!!!

        bool isC0 = false;//true;

        while((!in.eof()) && (line < count))
        {
            // input
            int imgId, imgStatus;
            std::vector<double> c(3);
            //double pitch, roll, yaw;
            double omega, phi, kappa;
            in >> imgId
               >> c[0] >> c[1] >> c[2]
  //             >> c[1] >> c[0] >> c[2]
               >> omega >> phi >> kappa
               >> imgStatus;

  //          omega = omega + 2.4;
  //          phi = phi + 0.0713;
  //          kappa = kappa - 0.5805;
  //          kappa = kappa + 2.4;
  //          phi = phi + 0.0713;
  //          omega = omega - 0.5805;

            if(isC0 == false)
            {
                c0 = c;
                isC0 = true;
            }

//            c[0] -= c0[0];
//            c[1] -= c0[1];
//            c[2] -= c0[2];

            ++line;



  //          double roll, pitch, yaw;

  //          roll = omega;
  //          pitch = phi;
  //          yaw = kappa;

  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
  //                sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
  //                0.0, 1.0, 0.0,
  //                -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
  //          Rx << 1.0, 0.0, 0.0,
  //                0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
  //                0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);
  //          Mat3 R_imu =  Rz*Rx*Ry;
  //          Mat3 changeAxis_M;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, -1.0, 0.0,
  //                          0.0, 0.0, -1.0;
  //          R_imu = changeAxis_M * R_imu;




  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
  //                0.0, 1.0, 0.0,
  //                -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);
  //          Rx << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

  ////          Mat3 R_imu =  Rx*Ry*Rz;
  ////          Mat3 R_imu =  Rx*Rz*Ry;//
  ////          Mat3 R_imu =  Ry*Rz*Rx;
  ////          Mat3 R_imu =  Ry*Rx*Rz;
  //          Mat3 R_imu =  Rz*Rx*Ry;//
  ////          Mat3 R_imu =  Rz*Ry*Rx;


            //use now
            Mat3 Rz, Ry, Rx;
            Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
                  sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
                  0.0, 0.0, 1.0;
            Ry << cos(phi/180.0*M_PI), 0.0, sin(phi/180.0*M_PI),
                  0.0, 1.0, 0.0,
                  -sin(phi/180.0*M_PI), 0.0, cos(phi/180.0*M_PI);
            Rx << 1.0, 0.0, 0.0,
                  0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
                  0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
  //          Mat3 R_imu =  Rz*Ry*Rx;
            Mat3 R_imu =  Rx*Ry*Rz;



  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, -cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), -cos(phi/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;

  //          Rz << cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI), 0.0,
  //                sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, -cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
  //                0.0, sin(omega/180.0*M_PI), -cos(omega/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;

  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;          Mat3 Rz, Ry, Rx;
  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(phi/180.0*M_PI), sin(phi/180.0*M_PI), 0.0,
  //                -sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
  //                0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;



  //          Mat3 R_imu =  Ry*Rz*Rx;
  //          Mat3 R_imu =  Rz*Ry*Rx;
  //          Mat3 R_imu =  Rx*Ry*Rz;


  //          Mat3 Rk, Ra, RA;
  //          Rk << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          RA << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ra << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


  //          Mat3 R_imu =  Ry*Rx*Rz;

  //          Mat3 changeAxis_M;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, -1.0, 0.0,
  //                          0.0, 0.0, -1.0;  //M1
  //          changeAxis_M << 0.0, 1.0, 0.0,
  //                          1.0, 0.0, 0.0,
  //                          0.0, 0.0, -1.0;  //M2

  //          R_imu = changeAxis_M * R_imu;
  //          R_imu = R_imu * changeAxis_M;





  //              tempGet = getListXYZ(-phi3, kappa3, -omega3);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //               tempGet = getListXYZ(phi3, kappa3, omega3);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //          tempGet = getListXYZ(phi1, kappa1, omega1);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //          Mat3 Rz, Ry, Rx;
  //          Rx << 1.0, 0.0, 0.0,
  //                 0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                 0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

  //          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
  //                 0.0, 1.0, 0.0,
  //                 -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);

  //          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                 sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                 0.0, 0.0, 1.0;


  //          Mat3 R_imu =  Ry*Rz*Rx;
  //          Mat3 R_imu =  Rz*Rx*Ry;
  //          Mat3 R_imu =  Rx*Rz*Ry;
  //          Mat3 R_imu =  Rx*Ry*Rz;

            Mat3 changeAxis_M;
            changeAxis_M << 1.0, 0.0, 0.0,
                            0.0, -1.0, 0.0,
                            0.0, 0.0, -1.0;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, 1.0, 0.0,
  //                          0.0, 0.0, -1.0;

  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, 1.0, 0.0,
  //                          0.0, 0.0, -1.0;
            Mat3 R_;
            R_ << 1.0, 0.0, 0.0,
                  0.0, 1.0, 0.0,
                  0.0, 0.0, -1.0;
            R_imu = changeAxis_M * R_imu;//* R_;


            //Mat3 R2;
            //ceres::AngleAxisToRotationMatrix(R2_a.data(), R2.data());

            //R_imu = changeAxis_M * R_imu * R2;//* R_;


            std::vector<double> R_imu_a(3);
            ceres::RotationMatrixToAngleAxis(R_imu.data(), R_imu_a.data());

            Vec3 C_gps_res = Vec3(c[0], c[1], c[2]);


            {
                bool find = false;
                for(Views::const_iterator itV = sfm_data_.views.begin(); itV != sfm_data_.views.end(); ++itV)
                {
                    if(itV->first == imgId)
                    {
                        find = true;
                        break;
                    }
                }
                if(find == false)
                    continue;
            }
            // check status
            if(imgId/100000 == 1)//back
            {
                if(imgStatus != 321)
                {
                    C_gps_Map[imgId + 100000] = c;
                    R_imu_Map[imgId + 100000] = R_imu_a;
                }

                //set back-camera pose
                Mat3 Rb_res = Rb_trans * R_imu_trans * R_imu;

                Vec3 tb_res;
                tb_res = Rb_trans * gps_trans - Rb_trans * R_imu_trans * R_imu * C_gps_res + tb_trans;
                global_rotations[imgId] = Rb_res;

                sfm_data_.poses[imgId] = Pose3(Rb_res, -Rb_res.transpose() * tb_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 1;
                poseStatusMap[imgId] = poseStatus;

            }else if(imgId/100000 == 2)//down
            {
                C_gps_Map[imgId] = c;
                R_imu_Map[imgId] = R_imu_a;

                //set down-camera pose
                Mat3 Rd_res = Rd_trans * R_imu_trans * R_imu;

                Vec3 td_res;
                td_res = Rd_trans * gps_trans - Rd_trans * R_imu_trans * R_imu * C_gps_res + td_trans;
                global_rotations[imgId] = Rd_res;

                sfm_data_.poses[imgId] = Pose3(Rd_res, -Rd_res.transpose() * td_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 2;
                poseStatusMap[imgId] = poseStatus;

            }else if(imgId/100000 == 3)//front
            {
                if(imgStatus != 321)
                {
                    C_gps_Map[imgId - 100000] = c;
                    R_imu_Map[imgId - 100000] = R_imu_a;
                }

                //set front-camera pose
                Mat3 Rf_res = Rf_trans * R_imu_trans * R_imu;

                Vec3 tf_res;
                tf_res = Rf_trans * gps_trans - Rf_trans * R_imu_trans * R_imu * C_gps_res + tf_trans;
                global_rotations[imgId] = Rf_res;

                sfm_data_.poses[imgId] = Pose3(Rf_res, -Rf_res.transpose() * tf_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 3;
                poseStatusMap[imgId] = poseStatus;

            }

        }
        in.close();

    }


    matching::PairWiseMatches tripletWise_matches;
    if (!Compute_Global_Translations_gps(global_rotations, tripletWise_matches))
    {
      std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
      return false;
    }

  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }
  if (!Adjust_init_water3_c1())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process_GCP_GPS_water4(std::string rtFilePath) {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if (set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

    Hash_Map<IndexT, Mat3> global_rotations;

    //init inter-constraint
    Mat3 Rb_trans, Rf_trans, Rd_trans;
    Vec3 tb_trans, tf_trans, td_trans;
    Mat3 R_imu_trans;
    Vec3 gps_trans;
    std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
    const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
    const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};
    {
        //test distortion
        Mat3 r_tb, r_tf;
        r_tb << 1.0, 0.0, 0.0,
                0.0, 0.787926, -0.61577,
                0.0, 0.61577, 0.787926;
        r_tf << 1.0, 0.0, 0.0,
                0.0, 0.78796,  0.615727,
                0.0, -0.615727, 0.78796;

  //      Rx << 1.0, 0.0, 0.0,
  //             0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //             0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


  //      r_tb << 1.0, 0.0, 0.0,
  //              0.0, cos(38.0/180.0*M_PI), -sin(38.0/180.0*M_PI),
  //              0.0, sin(38.0/180.0*M_PI), cos(38.0/180.0*M_PI);

  //      r_tf << 1.0, 0.0, 0.0,
  //              0.0, cos(-38.0/180.0*M_PI), -sin(-38.0/180.0*M_PI),
  //              0.0, sin(-38.0/180.0*M_PI), cos(-38.0/180.0*M_PI);

  //      r_tb << cos(-38.0/180.0*M_PI), 0.0, sin(-38.0/180.0*M_PI),//-38
  //              0.0, 1.0, 0.0,
  //              -sin(-38.0/180.0*M_PI), 0.0, cos(-38.0/180.0*M_PI);
  //      r_tf << cos(38.0/180.0*M_PI), 0.0, sin(38.0/180.0*M_PI),//38
  //              0.0, 1.0, 0.0,
  //              -sin(38.0/180.0*M_PI), 0.0, cos(38.0/180.0*M_PI);
        Vec3 RA_tb, RA_tf;
        ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
        ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

        //set reference
        double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), -RA_tb(2)};
        double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), -RA_tf(2)};
  //      double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), RA_tb(2)};
  //      double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), RA_tf(2)};
        double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};

        Vec3 c_bd, c_fd;//b-d, f-d
        c_bd << 0.0, 0.098, 0.027;
        c_fd << 0.0, -0.101, 0.028;
  //      c_bd << -0.098, 0.0, -0.027;
  //      c_fd << 0.101, 0.0, -0.028;

  //      c_bd << 0.098, 0.0, -0.027;
  //      c_fd << -0.101, 0.0, -0.028;
  //      c_bd << 0.0, 1.0, 0.027;
  //      c_fd << 0.0, -0.101, 0.028;

        Vec3 t_db, t_df, t_dd;
        t_db = -r_tb*c_bd;
        t_df = -r_tf*c_fd;

        //double S = 1470.531783220221;
      //        t_db << 0.0, 0.098, 0.027;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
      //        t_df << 0.0, -0.101, 0.028;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
        t_dd << 0.0, 0.0, 0.0;

        transformation_br = {R_br_angleAxis[0], R_br_angleAxis[1], R_br_angleAxis[2],
                             t_db(0), t_db(1), t_db(2)};
        transformation_fr = {R_fr_angleAxis[0], R_fr_angleAxis[1], R_fr_angleAxis[2],
                             t_df(0), t_df(1), t_df(2)};
        transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                             t_dd(0), t_dd(1), t_dd(2)};

  //      const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
  //      const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;

        //prepare for setting pose
        ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rb_trans.data());
        tb_trans = Vec3(transformation_br[3], transformation_br[4], transformation_br[5]);
        ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rf_trans.data());
        tf_trans = Vec3(transformation_fr[3], transformation_fr[4], transformation_fr[5]);
        ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rd_trans.data());
        td_trans = Vec3(transformation_dr[3], transformation_dr[4], transformation_dr[5]);


        ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
        gps_trans = Vec3(tramsformation_x_gps[0], tramsformation_x_gps[1], tramsformation_x_gps[2]);

    }


    //get gps and imu data
    PoseStatusMap poseStatusMap;
    Hash_Map<IndexT, std::vector<double> > C_gps_Map;
    Hash_Map<IndexT, std::vector<double> > R_imu_Map;
    std::vector<double> c0(3);
    {
      //        /Hash_Map<IndexT, std::vector<double> > R_imu;
  //            Hash_Map<IndexT, std::vector<double> > C_gps;
        std::string tFilePath = rtFilePath;  //!!!
        std::ifstream in;
        in.open(tFilePath);
        if(!in)
        {
            std::cout << "open tFilePath file error!" << std::endl;
            return false;
        }
        std::size_t line = 0;
        int count = -1;
        in >> count;  ///!!!

        bool isC0 = false;//true;

        while((!in.eof()) && (line < count))
        {
            // input
            int imgId, imgStatus;
            std::vector<double> c(3);
            //double pitch, roll, yaw;
            double omega, phi, kappa;
            in >> imgId
               >> c[0] >> c[1] >> c[2]
  //             >> c[1] >> c[0] >> c[2]
               >> omega >> phi >> kappa
               >> imgStatus;

  //          omega = omega + 2.4;
  //          phi = phi + 0.0713;
  //          kappa = kappa - 0.5805;
  //          kappa = kappa + 2.4;
  //          phi = phi + 0.0713;
  //          omega = omega - 0.5805;

            if(isC0 == false)
            {
                c0 = c;
                isC0 = true;
            }

//            c[0] -= c0[0];
//            c[1] -= c0[1];
//            c[2] -= c0[2];

            ++line;



  //          double roll, pitch, yaw;

  //          roll = omega;
  //          pitch = phi;
  //          yaw = kappa;

  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
  //                sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
  //                0.0, 1.0, 0.0,
  //                -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
  //          Rx << 1.0, 0.0, 0.0,
  //                0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
  //                0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);
  //          Mat3 R_imu =  Rz*Rx*Ry;
  //          Mat3 changeAxis_M;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, -1.0, 0.0,
  //                          0.0, 0.0, -1.0;
  //          R_imu = changeAxis_M * R_imu;




  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
  //                0.0, 1.0, 0.0,
  //                -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);
  //          Rx << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

  ////          Mat3 R_imu =  Rx*Ry*Rz;
  ////          Mat3 R_imu =  Rx*Rz*Ry;//
  ////          Mat3 R_imu =  Ry*Rz*Rx;
  ////          Mat3 R_imu =  Ry*Rx*Rz;
  //          Mat3 R_imu =  Rz*Rx*Ry;//
  ////          Mat3 R_imu =  Rz*Ry*Rx;


            //use now
            Mat3 Rz, Ry, Rx;
            Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
                  sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
                  0.0, 0.0, 1.0;
            Ry << cos(phi/180.0*M_PI), 0.0, sin(phi/180.0*M_PI),
                  0.0, 1.0, 0.0,
                  -sin(phi/180.0*M_PI), 0.0, cos(phi/180.0*M_PI);
            Rx << 1.0, 0.0, 0.0,
                  0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
                  0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
  //          Mat3 R_imu =  Rz*Ry*Rx;
            Mat3 R_imu =  Rx*Ry*Rz;



  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, -cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), -cos(phi/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;

  //          Rz << cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI), 0.0,
  //                sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, -cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
  //                0.0, sin(omega/180.0*M_PI), -cos(omega/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;

  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;          Mat3 Rz, Ry, Rx;
  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(phi/180.0*M_PI), sin(phi/180.0*M_PI), 0.0,
  //                -sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
  //                0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;



  //          Mat3 R_imu =  Ry*Rz*Rx;
  //          Mat3 R_imu =  Rz*Ry*Rx;
  //          Mat3 R_imu =  Rx*Ry*Rz;


  //          Mat3 Rk, Ra, RA;
  //          Rk << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          RA << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ra << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


  //          Mat3 R_imu =  Ry*Rx*Rz;

  //          Mat3 changeAxis_M;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, -1.0, 0.0,
  //                          0.0, 0.0, -1.0;  //M1
  //          changeAxis_M << 0.0, 1.0, 0.0,
  //                          1.0, 0.0, 0.0,
  //                          0.0, 0.0, -1.0;  //M2

  //          R_imu = changeAxis_M * R_imu;
  //          R_imu = R_imu * changeAxis_M;





  //              tempGet = getListXYZ(-phi3, kappa3, -omega3);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //               tempGet = getListXYZ(phi3, kappa3, omega3);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //          tempGet = getListXYZ(phi1, kappa1, omega1);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //          Mat3 Rz, Ry, Rx;
  //          Rx << 1.0, 0.0, 0.0,
  //                 0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                 0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

  //          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
  //                 0.0, 1.0, 0.0,
  //                 -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);

  //          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                 sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                 0.0, 0.0, 1.0;


  //          Mat3 R_imu =  Ry*Rz*Rx;
  //          Mat3 R_imu =  Rz*Rx*Ry;
  //          Mat3 R_imu =  Rx*Rz*Ry;
  //          Mat3 R_imu =  Rx*Ry*Rz;

            Mat3 changeAxis_M;
            changeAxis_M << 1.0, 0.0, 0.0,
                            0.0, -1.0, 0.0,
                            0.0, 0.0, -1.0;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, 1.0, 0.0,
  //                          0.0, 0.0, -1.0;

  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, 1.0, 0.0,
  //                          0.0, 0.0, -1.0;
            Mat3 R_;
            R_ << 1.0, 0.0, 0.0,
                  0.0, 1.0, 0.0,
                  0.0, 0.0, -1.0;
            R_imu = changeAxis_M * R_imu;//* R_;


            //Mat3 R2;
            //ceres::AngleAxisToRotationMatrix(R2_a.data(), R2.data());

            //R_imu = changeAxis_M * R_imu * R2;//* R_;


            std::vector<double> R_imu_a(3);
            ceres::RotationMatrixToAngleAxis(R_imu.data(), R_imu_a.data());

            Vec3 C_gps_res = Vec3(c[0], c[1], c[2]);


            {
                bool find = false;
                for(Views::const_iterator itV = sfm_data_.views.begin(); itV != sfm_data_.views.end(); ++itV)
                {
                    if(itV->first == imgId)
                    {
                        find = true;
                        break;
                    }
                }
                if(find == false)
                    continue;
            }
            // check status
            if(imgId/100000 == 1)//back
            {
                if(imgStatus != 321)
                {
                    C_gps_Map[imgId + 100000] = c;
                    R_imu_Map[imgId + 100000] = R_imu_a;
                }

                //set back-camera pose
                Mat3 Rb_res = Rb_trans * R_imu_trans * R_imu;

                Vec3 tb_res;
                tb_res = Rb_trans * gps_trans - Rb_trans * R_imu_trans * R_imu * C_gps_res + tb_trans;
                global_rotations[imgId] = Rb_res;

                sfm_data_.poses[imgId] = Pose3(Rb_res, -Rb_res.transpose() * tb_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 1;
                poseStatusMap[imgId] = poseStatus;

            }else if(imgId/100000 == 2)//down
            {
                C_gps_Map[imgId] = c;
                R_imu_Map[imgId] = R_imu_a;

                //set down-camera pose
                Mat3 Rd_res = Rd_trans * R_imu_trans * R_imu;

                Vec3 td_res;
                td_res = Rd_trans * gps_trans - Rd_trans * R_imu_trans * R_imu * C_gps_res + td_trans;
                global_rotations[imgId] = Rd_res;

                sfm_data_.poses[imgId] = Pose3(Rd_res, -Rd_res.transpose() * td_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 2;
                poseStatusMap[imgId] = poseStatus;

            }else if(imgId/100000 == 3)//front
            {
                if(imgStatus != 321)
                {
                    C_gps_Map[imgId - 100000] = c;
                    R_imu_Map[imgId - 100000] = R_imu_a;
                }

                //set front-camera pose
                Mat3 Rf_res = Rf_trans * R_imu_trans * R_imu;

                Vec3 tf_res;
                tf_res = Rf_trans * gps_trans - Rf_trans * R_imu_trans * R_imu * C_gps_res + tf_trans;
                global_rotations[imgId] = Rf_res;

                sfm_data_.poses[imgId] = Pose3(Rf_res, -Rf_res.transpose() * tf_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 3;
                poseStatusMap[imgId] = poseStatus;

            }

        }
        in.close();

    }


    matching::PairWiseMatches tripletWise_matches;
    if (!Compute_Global_Translations_gps(global_rotations, tripletWise_matches))
    {
      std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
      return false;
    }

  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }
  if (!Adjust_init_water4_c1_2())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process_GCP_GPS_water5(std::string rtFilePath) {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if (set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

    openMVG::rotation_averaging::RelativeRotations relatives_R;
    Compute_Relative_Rotations(relatives_R);

    Hash_Map<IndexT, Mat3> global_rotations;
    if (!Compute_Global_Rotations(relatives_R, global_rotations))
    {
      std::cerr << "GlobalSfM:: Rotation Averaging failure!" << std::endl;
      return false;
    }

    matching::PairWiseMatches  tripletWise_matches;
    if (!Compute_Global_Translations(global_rotations, tripletWise_matches))
    {
      std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
      return false;
    }
    if (!Compute_Initial_Structure(tripletWise_matches))
    {
      std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
      return false;
    }

    ///save for g2o
    {
//        std::ofstream out_landmark;
//        std::string out_landmark_path = sfm_data_.s_root_path + "/landmarks.res";
//        out_landmark.open(out_landmark_path);
//        if(!out_landmark)
//        {
//            std::cout << "create " << out_landmark_path << " file failed!" << std::endl;
//        }
//  //      std::ofstream out_X;
//  //      std::string out_X_path = sfm_data_.s_root_path + "/X_all.res";
//  //      out_X.open(out_X_path);
//  //      if(!out_X)
//  //      {
//  //          std::cout << "create " << out_X_path << " file failed!" << std::endl;
//  //      }
//        //
//        int countL = sfm_data_.structure.size();
//        out_landmark << countL << std::endl;
//        //
//        for(Landmarks::const_iterator itL = sfm_data_.structure.begin(); itL != sfm_data_.structure.end(); ++itL)
//        {
//            int Xid = itL->first;
//            Vec3 X = itL->second.X;
//            out_landmark << Xid << " " << X(0) << " " << X(1) << " " << X(2) << std::endl;
//            //
//            int countObs = itL->second.obs.size();
//            out_landmark << countObs << " ";
//            for(Observations::const_iterator itObs = itL->second.obs.begin(); itObs != itL->second.obs.end(); ++itObs)
//            {
//                int imgId = itObs->first;
//                Vec2 featId = itObs->second.x;
//                out_landmark << imgId << " " << featId(0) << " " << featId(1) << " ";
//            }
//            out_landmark << std::endl;
//        }
//        out_landmark.close();
//        //
//        /// save extrinsics
//        std::ofstream out_extrinsic;
//        std::string out_extrinsic_path = sfm_data_.s_root_path + "/extrinsics.res";
//        out_extrinsic.open(out_extrinsic_path);
//        if(!out_extrinsic)
//        {
//            std::cout << "create " << out_extrinsic_path << " file failed!" << std::endl;
//        }
//        out_extrinsic << sfm_data_.poses.size() << std::endl;
//        for(Poses::const_iterator itPose = sfm_data_.poses.begin(); itPose != sfm_data_.poses.end(); ++itPose)
//        {
//            Vec3 t = itPose->second.translation();
//            std::cout << "t" << std::endl;
//            Mat3 R = itPose->second.rotation();
//            std::cout << "R" << std::endl;
//            Vec4 q;
//  //          std::vector<double> q;
//            ceres::RotationMatrixToQuaternion(R.data(), q.data());
//            std::cout << "q" << std::endl;
//            //
//  //          out_extrinsic << itPose->first << " " << q[0] << " "<< q[1] << " " << q[2] << " " << q[3] << " " << t(0) << " " << t(1) << " " << t(2) << std::endl;
//            out_extrinsic << itPose->first << " " << q(0) << " "<< q(1) << " " << q(2) << " " << q(3) << " " << t(0) << " " << t(1) << " " << t(2) << std::endl;
//        }
//        out_extrinsic.close();



    }


//  if (!Adjust_init_water5_c1())
//  if (!Adjust_init_water5_c4())
  if (!Adjust_init_water_c4())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process_GCP_GPS_water6(std::string rtFilePath) {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if (set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

    Hash_Map<IndexT, Mat3> global_rotations;

    //init inter-constraint
    Mat3 Rb_trans, Rf_trans, Rd_trans;
    Vec3 tb_trans, tf_trans, td_trans;
    Mat3 R_imu_trans;
    Vec3 gps_trans;
    std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
    const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
    const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};
    {
        //test distortion
        Mat3 r_tb, r_tf;
        r_tb << 1.0, 0.0, 0.0,
                0.0, 0.787926, -0.61577,
                0.0, 0.61577, 0.787926;
        r_tf << 1.0, 0.0, 0.0,
                0.0, 0.78796,  0.615727,
                0.0, -0.615727, 0.78796;

  //      Rx << 1.0, 0.0, 0.0,
  //             0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //             0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


  //      r_tb << 1.0, 0.0, 0.0,
  //              0.0, cos(38.0/180.0*M_PI), -sin(38.0/180.0*M_PI),
  //              0.0, sin(38.0/180.0*M_PI), cos(38.0/180.0*M_PI);

  //      r_tf << 1.0, 0.0, 0.0,
  //              0.0, cos(-38.0/180.0*M_PI), -sin(-38.0/180.0*M_PI),
  //              0.0, sin(-38.0/180.0*M_PI), cos(-38.0/180.0*M_PI);

  //      r_tb << cos(-38.0/180.0*M_PI), 0.0, sin(-38.0/180.0*M_PI),//-38
  //              0.0, 1.0, 0.0,
  //              -sin(-38.0/180.0*M_PI), 0.0, cos(-38.0/180.0*M_PI);
  //      r_tf << cos(38.0/180.0*M_PI), 0.0, sin(38.0/180.0*M_PI),//38
  //              0.0, 1.0, 0.0,
  //              -sin(38.0/180.0*M_PI), 0.0, cos(38.0/180.0*M_PI);
        Vec3 RA_tb, RA_tf;
        ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
        ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

        //set reference
        double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), -RA_tb(2)};
        double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), -RA_tf(2)};
  //      double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), RA_tb(2)};
  //      double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), RA_tf(2)};
        double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};

        Vec3 c_bd, c_fd;//b-d, f-d
        c_bd << 0.0, 0.098, 0.027;
        c_fd << 0.0, -0.101, 0.028;
  //      c_bd << -0.098, 0.0, -0.027;
  //      c_fd << 0.101, 0.0, -0.028;

  //      c_bd << 0.098, 0.0, -0.027;
  //      c_fd << -0.101, 0.0, -0.028;
  //      c_bd << 0.0, 1.0, 0.027;
  //      c_fd << 0.0, -0.101, 0.028;

        Vec3 t_db, t_df, t_dd;
        t_db = -r_tb*c_bd;
        t_df = -r_tf*c_fd;

        //double S = 1470.531783220221;
      //        t_db << 0.0, 0.098, 0.027;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
      //        t_df << 0.0, -0.101, 0.028;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
        t_dd << 0.0, 0.0, 0.0;

        transformation_br = {R_br_angleAxis[0], R_br_angleAxis[1], R_br_angleAxis[2],
                             t_db(0), t_db(1), t_db(2)};
        transformation_fr = {R_fr_angleAxis[0], R_fr_angleAxis[1], R_fr_angleAxis[2],
                             t_df(0), t_df(1), t_df(2)};
        transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                             t_dd(0), t_dd(1), t_dd(2)};

  //      const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
  //      const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;

        //prepare for setting pose
        ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rb_trans.data());
        tb_trans = Vec3(transformation_br[3], transformation_br[4], transformation_br[5]);
        ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rf_trans.data());
        tf_trans = Vec3(transformation_fr[3], transformation_fr[4], transformation_fr[5]);
        ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rd_trans.data());
        td_trans = Vec3(transformation_dr[3], transformation_dr[4], transformation_dr[5]);


        ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
        gps_trans = Vec3(tramsformation_x_gps[0], tramsformation_x_gps[1], tramsformation_x_gps[2]);

    }


    //get gps and imu data
    PoseStatusMap poseStatusMap;
    Hash_Map<IndexT, std::vector<double> > C_gps_Map;
    Hash_Map<IndexT, std::vector<double> > R_imu_Map;
    std::vector<double> c0(3);
    {
      //        /Hash_Map<IndexT, std::vector<double> > R_imu;
  //            Hash_Map<IndexT, std::vector<double> > C_gps;
        std::string tFilePath = rtFilePath;  //!!!
        std::ifstream in;
        in.open(tFilePath);
        if(!in)
        {
            std::cout << "open tFilePath file error!" << std::endl;
            return false;
        }
        std::size_t line = 0;
        int count = -1;
        in >> count;  ///!!!

        bool isC0 = false;//true;

        while((!in.eof()) && (line < count))
        {
            // input
            int imgId, imgStatus;
            std::vector<double> c(3);
            //double pitch, roll, yaw;
            double omega, phi, kappa;
            in >> imgId
               >> c[0] >> c[1] >> c[2]
  //             >> c[1] >> c[0] >> c[2]
               >> omega >> phi >> kappa
               >> imgStatus;

  //          omega = omega + 2.4;
  //          phi = phi + 0.0713;
  //          kappa = kappa - 0.5805;
  //          kappa = kappa + 2.4;
  //          phi = phi + 0.0713;
  //          omega = omega - 0.5805;

            if(isC0 == false)
            {
                c0 = c;
                isC0 = true;
            }

//            c[0] -= c0[0];
//            c[1] -= c0[1];
//            c[2] -= 280.0;//c0[2];

            ++line;



  //          double roll, pitch, yaw;

  //          roll = omega;
  //          pitch = phi;
  //          yaw = kappa;

  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
  //                sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
  //                0.0, 1.0, 0.0,
  //                -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
  //          Rx << 1.0, 0.0, 0.0,
  //                0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
  //                0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);
  //          Mat3 R_imu =  Rz*Rx*Ry;
  //          Mat3 changeAxis_M;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, -1.0, 0.0,
  //                          0.0, 0.0, -1.0;
  //          R_imu = changeAxis_M * R_imu;




  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
  //                0.0, 1.0, 0.0,
  //                -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);
  //          Rx << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

  ////          Mat3 R_imu =  Rx*Ry*Rz;
  ////          Mat3 R_imu =  Rx*Rz*Ry;//
  ////          Mat3 R_imu =  Ry*Rz*Rx;
  ////          Mat3 R_imu =  Ry*Rx*Rz;
  //          Mat3 R_imu =  Rz*Rx*Ry;//
  ////          Mat3 R_imu =  Rz*Ry*Rx;


            //use now
            Mat3 Rz, Ry, Rx;
            Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
                  sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
                  0.0, 0.0, 1.0;
            Ry << cos(phi/180.0*M_PI), 0.0, sin(phi/180.0*M_PI),
                  0.0, 1.0, 0.0,
                  -sin(phi/180.0*M_PI), 0.0, cos(phi/180.0*M_PI);
            Rx << 1.0, 0.0, 0.0,
                  0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
                  0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
  //          Mat3 R_imu =  Rz*Ry*Rx;
            Mat3 R_imu =  Rx*Ry*Rz;



  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, -cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), -cos(phi/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;

  //          Rz << cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI), 0.0,
  //                sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, -cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
  //                0.0, sin(omega/180.0*M_PI), -cos(omega/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;

  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;          Mat3 Rz, Ry, Rx;
  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(phi/180.0*M_PI), sin(phi/180.0*M_PI), 0.0,
  //                -sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
  //                0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;



  //          Mat3 R_imu =  Ry*Rz*Rx;
  //          Mat3 R_imu =  Rz*Ry*Rx;
  //          Mat3 R_imu =  Rx*Ry*Rz;


  //          Mat3 Rk, Ra, RA;
  //          Rk << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          RA << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ra << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


  //          Mat3 R_imu =  Ry*Rx*Rz;

  //          Mat3 changeAxis_M;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, -1.0, 0.0,
  //                          0.0, 0.0, -1.0;  //M1
  //          changeAxis_M << 0.0, 1.0, 0.0,
  //                          1.0, 0.0, 0.0,
  //                          0.0, 0.0, -1.0;  //M2

  //          R_imu = changeAxis_M * R_imu;
  //          R_imu = R_imu * changeAxis_M;





  //              tempGet = getListXYZ(-phi3, kappa3, -omega3);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //               tempGet = getListXYZ(phi3, kappa3, omega3);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //          tempGet = getListXYZ(phi1, kappa1, omega1);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //          Mat3 Rz, Ry, Rx;
  //          Rx << 1.0, 0.0, 0.0,
  //                 0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                 0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

  //          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
  //                 0.0, 1.0, 0.0,
  //                 -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);

  //          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                 sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                 0.0, 0.0, 1.0;


  //          Mat3 R_imu =  Ry*Rz*Rx;
  //          Mat3 R_imu =  Rz*Rx*Ry;
  //          Mat3 R_imu =  Rx*Rz*Ry;
  //          Mat3 R_imu =  Rx*Ry*Rz;

            Mat3 changeAxis_M;
            changeAxis_M << 1.0, 0.0, 0.0,
                            0.0, -1.0, 0.0,
                            0.0, 0.0, -1.0;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, 1.0, 0.0,
  //                          0.0, 0.0, -1.0;

  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, 1.0, 0.0,
  //                          0.0, 0.0, -1.0;
            Mat3 R_;
            R_ << 1.0, 0.0, 0.0,
                  0.0, 1.0, 0.0,
                  0.0, 0.0, -1.0;
            R_imu = changeAxis_M * R_imu;//* R_;


            //Mat3 R2;
            //ceres::AngleAxisToRotationMatrix(R2_a.data(), R2.data());

            //R_imu = changeAxis_M * R_imu * R2;//* R_;


            std::vector<double> R_imu_a(3);
            ceres::RotationMatrixToAngleAxis(R_imu.data(), R_imu_a.data());

            Vec3 C_gps_res = Vec3(c[0], c[1], c[2]);


            {
                bool find = false;
                for(Views::const_iterator itV = sfm_data_.views.begin(); itV != sfm_data_.views.end(); ++itV)
                {
                    if(itV->first == imgId)
                    {
                        find = true;
                        break;
                    }
                }
                if(find == false)
                    continue;
            }
            // check status
            if(imgId/100000 == 1)//back
            {
                if(imgStatus != 321)
                {
                    C_gps_Map[imgId + 100000] = c;
                    R_imu_Map[imgId + 100000] = R_imu_a;
                }

                //set back-camera pose
                Mat3 Rb_res = Rb_trans * R_imu_trans * R_imu;

                Vec3 tb_res;
                tb_res = Rb_trans * gps_trans - Rb_trans * R_imu_trans * R_imu * C_gps_res + tb_trans;
                global_rotations[imgId] = Rb_res;

                sfm_data_.poses[imgId] = Pose3(Rb_res, -Rb_res.transpose() * tb_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 1;
                poseStatusMap[imgId] = poseStatus;

            }else if(imgId/100000 == 2)//down
            {
                C_gps_Map[imgId] = c;
                R_imu_Map[imgId] = R_imu_a;

                //set down-camera pose
                Mat3 Rd_res = Rd_trans * R_imu_trans * R_imu;

                Vec3 td_res;
                td_res = Rd_trans * gps_trans - Rd_trans * R_imu_trans * R_imu * C_gps_res + td_trans;
                global_rotations[imgId] = Rd_res;

                sfm_data_.poses[imgId] = Pose3(Rd_res, -Rd_res.transpose() * td_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 2;
                poseStatusMap[imgId] = poseStatus;

            }else if(imgId/100000 == 3)//front
            {
                if(imgStatus != 321)
                {
                    C_gps_Map[imgId - 100000] = c;
                    R_imu_Map[imgId - 100000] = R_imu_a;
                }

                //set front-camera pose
                Mat3 Rf_res = Rf_trans * R_imu_trans * R_imu;

                Vec3 tf_res;
                tf_res = Rf_trans * gps_trans - Rf_trans * R_imu_trans * R_imu * C_gps_res + tf_trans;
                global_rotations[imgId] = Rf_res;

                sfm_data_.poses[imgId] = Pose3(Rf_res, -Rf_res.transpose() * tf_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 3;
                poseStatusMap[imgId] = poseStatus;

            }

        }
        in.close();

    }


    matching::PairWiseMatches tripletWise_matches;
    if (!Compute_Global_Translations_gps(global_rotations, tripletWise_matches))
    {
      std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
      return false;
    }

  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }
  if (!Adjust_init_water6_c1())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}



bool GlobalSfMReconstructionEngine_RelativeMotions::Process_GCP_GPS_6(std::string rtFilePath) {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if (set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

    Hash_Map<IndexT, Mat3> global_rotations;

    //init inter-constraint
    Mat3 Rb_trans, Rf_trans, Rd_trans;
    Vec3 tb_trans, tf_trans, td_trans;
    Mat3 R_imu_trans;
    Vec3 gps_trans;
    std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
    const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
    const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};
    {
        //test distortion
        Mat3 r_tb, r_tf;
        r_tb << 1.0, 0.0, 0.0,
                0.0, 0.787926, -0.61577,
                0.0, 0.61577, 0.787926;
        r_tf << 1.0, 0.0, 0.0,
                0.0, 0.78796,  0.615727,
                0.0, -0.615727, 0.78796;

  //      Rx << 1.0, 0.0, 0.0,
  //             0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //             0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


  //      r_tb << 1.0, 0.0, 0.0,
  //              0.0, cos(38.0/180.0*M_PI), -sin(38.0/180.0*M_PI),
  //              0.0, sin(38.0/180.0*M_PI), cos(38.0/180.0*M_PI);

  //      r_tf << 1.0, 0.0, 0.0,
  //              0.0, cos(-38.0/180.0*M_PI), -sin(-38.0/180.0*M_PI),
  //              0.0, sin(-38.0/180.0*M_PI), cos(-38.0/180.0*M_PI);

  //      r_tb << cos(-38.0/180.0*M_PI), 0.0, sin(-38.0/180.0*M_PI),//-38
  //              0.0, 1.0, 0.0,
  //              -sin(-38.0/180.0*M_PI), 0.0, cos(-38.0/180.0*M_PI);
  //      r_tf << cos(38.0/180.0*M_PI), 0.0, sin(38.0/180.0*M_PI),//38
  //              0.0, 1.0, 0.0,
  //              -sin(38.0/180.0*M_PI), 0.0, cos(38.0/180.0*M_PI);
        Vec3 RA_tb, RA_tf;
        ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
        ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

        //set reference
        double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), -RA_tb(2)};
        double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), -RA_tf(2)};
  //      double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), RA_tb(2)};
  //      double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), RA_tf(2)};
        double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};

        Vec3 c_bd, c_fd;//b-d, f-d
        c_bd << 0.0, 0.098, 0.027;
        c_fd << 0.0, -0.101, 0.028;
  //      c_bd << -0.098, 0.0, -0.027;
  //      c_fd << 0.101, 0.0, -0.028;

  //      c_bd << 0.098, 0.0, -0.027;
  //      c_fd << -0.101, 0.0, -0.028;
  //      c_bd << 0.0, 1.0, 0.027;
  //      c_fd << 0.0, -0.101, 0.028;

        Vec3 t_db, t_df, t_dd;
        t_db = -r_tb*c_bd;
        t_df = -r_tf*c_fd;

        //double S = 1470.531783220221;
      //        t_db << 0.0, 0.098, 0.027;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
      //        t_df << 0.0, -0.101, 0.028;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
        t_dd << 0.0, 0.0, 0.0;

        transformation_br = {R_br_angleAxis[0], R_br_angleAxis[1], R_br_angleAxis[2],
                             t_db(0), t_db(1), t_db(2)};
        transformation_fr = {R_fr_angleAxis[0], R_fr_angleAxis[1], R_fr_angleAxis[2],
                             t_df(0), t_df(1), t_df(2)};
        transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                             t_dd(0), t_dd(1), t_dd(2)};

  //      const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
  //      const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;

        //prepare for setting pose
        ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rb_trans.data());
        tb_trans = Vec3(transformation_br[3], transformation_br[4], transformation_br[5]);
        ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rf_trans.data());
        tf_trans = Vec3(transformation_fr[3], transformation_fr[4], transformation_fr[5]);
        ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rd_trans.data());
        td_trans = Vec3(transformation_dr[3], transformation_dr[4], transformation_dr[5]);


        ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
        gps_trans = Vec3(tramsformation_x_gps[0], tramsformation_x_gps[1], tramsformation_x_gps[2]);

    }


    //get gps and imu data
    PoseStatusMap poseStatusMap;
    Hash_Map<IndexT, std::vector<double> > C_gps_Map;
    Hash_Map<IndexT, std::vector<double> > R_imu_Map;
    std::vector<double> c0(3);
    {
      //        /Hash_Map<IndexT, std::vector<double> > R_imu;
  //            Hash_Map<IndexT, std::vector<double> > C_gps;
        std::string tFilePath = rtFilePath;  //!!!
        std::ifstream in;
        in.open(tFilePath);
        if(!in)
        {
            std::cout << "open tFilePath file error!" << std::endl;
            return false;
        }
        std::size_t line = 0;
        int count = -1;
        in >> count;  ///!!!

        bool isC0 = false;//true;

        while((!in.eof()) && (line < count))
        {
            // input
            int imgId, imgStatus;
            std::vector<double> c(3);
            //double pitch, roll, yaw;
            double omega, phi, kappa;
            in >> imgId
               >> c[0] >> c[1] >> c[2]
  //             >> c[1] >> c[0] >> c[2]
               >> omega >> phi >> kappa
               >> imgStatus;

  //          omega = omega + 2.4;
  //          phi = phi + 0.0713;
  //          kappa = kappa - 0.5805;
  //          kappa = kappa + 2.4;
  //          phi = phi + 0.0713;
  //          omega = omega - 0.5805;

            if(isC0 == false)
            {
                c0 = c;
                isC0 = true;
            }

  //          c[0] -= c0[0];
  //          c[1] -= c0[1];
  //          c[2] -= c0[2];

            ++line;



  //          double roll, pitch, yaw;

  //          roll = omega;
  //          pitch = phi;
  //          yaw = kappa;

  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
  //                sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
  //                0.0, 1.0, 0.0,
  //                -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
  //          Rx << 1.0, 0.0, 0.0,
  //                0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
  //                0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);
  //          Mat3 R_imu =  Rz*Rx*Ry;
  //          Mat3 changeAxis_M;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, -1.0, 0.0,
  //                          0.0, 0.0, -1.0;
  //          R_imu = changeAxis_M * R_imu;




  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
  //                0.0, 1.0, 0.0,
  //                -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);
  //          Rx << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

  ////          Mat3 R_imu =  Rx*Ry*Rz;
  ////          Mat3 R_imu =  Rx*Rz*Ry;//
  ////          Mat3 R_imu =  Ry*Rz*Rx;
  ////          Mat3 R_imu =  Ry*Rx*Rz;
  //          Mat3 R_imu =  Rz*Rx*Ry;//
  ////          Mat3 R_imu =  Rz*Ry*Rx;


            //use now
            Mat3 Rz, Ry, Rx;
            Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
                  sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
                  0.0, 0.0, 1.0;
            Ry << cos(phi/180.0*M_PI), 0.0, sin(phi/180.0*M_PI),
                  0.0, 1.0, 0.0,
                  -sin(phi/180.0*M_PI), 0.0, cos(phi/180.0*M_PI);
            Rx << 1.0, 0.0, 0.0,
                  0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
                  0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
  //          Mat3 R_imu =  Rz*Ry*Rx;
            Mat3 R_imu =  Rx*Ry*Rz;



  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, -cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), -cos(phi/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;

  //          Rz << cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI), 0.0,
  //                sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, -cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
  //                0.0, sin(omega/180.0*M_PI), -cos(omega/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;

  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;          Mat3 Rz, Ry, Rx;
  //          Mat3 Rz, Ry, Rx;
  //          Rz << cos(phi/180.0*M_PI), sin(phi/180.0*M_PI), 0.0,
  //                -sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ry << 1.0, 0.0, 0.0,
  //                0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
  //                0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
  //          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;



  //          Mat3 R_imu =  Ry*Rz*Rx;
  //          Mat3 R_imu =  Rz*Ry*Rx;
  //          Mat3 R_imu =  Rx*Ry*Rz;


  //          Mat3 Rk, Ra, RA;
  //          Rk << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          RA << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
  //                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
  //                0.0, 0.0, 1.0;
  //          Ra << 1.0, 0.0, 0.0,
  //                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


  //          Mat3 R_imu =  Ry*Rx*Rz;

  //          Mat3 changeAxis_M;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, -1.0, 0.0,
  //                          0.0, 0.0, -1.0;  //M1
  //          changeAxis_M << 0.0, 1.0, 0.0,
  //                          1.0, 0.0, 0.0,
  //                          0.0, 0.0, -1.0;  //M2

  //          R_imu = changeAxis_M * R_imu;
  //          R_imu = R_imu * changeAxis_M;





  //              tempGet = getListXYZ(-phi3, kappa3, -omega3);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //               tempGet = getListXYZ(phi3, kappa3, omega3);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //          tempGet = getListXYZ(phi1, kappa1, omega1);
  //          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

  //          Mat3 Rz, Ry, Rx;
  //          Rx << 1.0, 0.0, 0.0,
  //                 0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
  //                 0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

  //          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
  //                 0.0, 1.0, 0.0,
  //                 -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);

  //          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
  //                 sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
  //                 0.0, 0.0, 1.0;


  //          Mat3 R_imu =  Ry*Rz*Rx;
  //          Mat3 R_imu =  Rz*Rx*Ry;
  //          Mat3 R_imu =  Rx*Rz*Ry;
  //          Mat3 R_imu =  Rx*Ry*Rz;

            Mat3 changeAxis_M;
            changeAxis_M << 1.0, 0.0, 0.0,
                            0.0, -1.0, 0.0,
                            0.0, 0.0, -1.0;
  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, 1.0, 0.0,
  //                          0.0, 0.0, -1.0;

  //          changeAxis_M << 1.0, 0.0, 0.0,
  //                          0.0, 1.0, 0.0,
  //                          0.0, 0.0, -1.0;
            Mat3 R_;
            R_ << 1.0, 0.0, 0.0,
                  0.0, 1.0, 0.0,
                  0.0, 0.0, -1.0;
            R_imu = changeAxis_M * R_imu;//* R_;


            //Mat3 R2;
            //ceres::AngleAxisToRotationMatrix(R2_a.data(), R2.data());

            //R_imu = changeAxis_M * R_imu * R2;//* R_;


            std::vector<double> R_imu_a(3);
            ceres::RotationMatrixToAngleAxis(R_imu.data(), R_imu_a.data());

            Vec3 C_gps_res = Vec3(c[0], c[1], c[2]);


            {
//                bool find = false;
//                for(Views::const_iterator itV = sfm_data_.views.begin(); itV != sfm_data_.views.end(); ++itV)
//                {
//                    if(itV->first == imgId)
//                    {
//                        find = true;
//                        break;
//                    }
//                }
//                if(find == false)
//                    continue;
            }
            // check status
            if(imgId/100000 == 1)//back
            {
                if(imgStatus != 321)
                {
                    C_gps_Map[imgId + 100000] = c;
                    R_imu_Map[imgId + 100000] = R_imu_a;
                }

                //set back-camera pose
                Mat3 Rb_res = Rb_trans * R_imu_trans * R_imu;

                Vec3 tb_res;
                tb_res = Rb_trans * gps_trans - Rb_trans * R_imu_trans * R_imu * C_gps_res + tb_trans;
                global_rotations[imgId] = Rb_res;

                sfm_data_.poses[imgId] = Pose3(Rb_res, -Rb_res.transpose() * tb_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 1;
                poseStatusMap[imgId] = poseStatus;

            }else if(imgId/100000 == 2)//down
            {
                C_gps_Map[imgId] = c;
                R_imu_Map[imgId] = R_imu_a;

                //set down-camera pose
                Mat3 Rd_res = Rd_trans * R_imu_trans * R_imu;

                Vec3 td_res;
                td_res = Rd_trans * gps_trans - Rd_trans * R_imu_trans * R_imu * C_gps_res + td_trans;
                global_rotations[imgId] = Rd_res;

                sfm_data_.poses[imgId] = Pose3(Rd_res, -Rd_res.transpose() * td_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 2;
                poseStatusMap[imgId] = poseStatus;

            }else if(imgId/100000 == 3)//front
            {
                if(imgStatus != 321)
                {
                    C_gps_Map[imgId - 100000] = c;
                    R_imu_Map[imgId - 100000] = R_imu_a;
                }

                //set front-camera pose
                Mat3 Rf_res = Rf_trans * R_imu_trans * R_imu;

                Vec3 tf_res;
                tf_res = Rf_trans * gps_trans - Rf_trans * R_imu_trans * R_imu * C_gps_res + tf_trans;
                global_rotations[imgId] = Rf_res;

                sfm_data_.poses[imgId] = Pose3(Rf_res, -Rf_res.transpose() * tf_res);

                PoseStatus poseStatus;
                poseStatus.pose = sfm_data_.poses[imgId];
                poseStatus.status = imgStatus;
                poseStatus.bdf = 3;
                poseStatusMap[imgId] = poseStatus;

            }

        }
        in.close();

    }


    matching::PairWiseMatches tripletWise_matches;
    if (!Compute_Global_Translations_gps(global_rotations, tripletWise_matches))
    {
      std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
      return false;
    }

  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }
  if (!Adjust_init_6())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process_GCP_GPS_4_margeTest(std::string rtFilePath) {

  if (!Adjust_init_4_margeTest())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process_step23(const std::map<int, bool>& fixedImgIdList, const std::map<int, bool>& fixedXIdList) {

  if (!Adjust_step23(fixedImgIdList, fixedXIdList))
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Process_GCP_step2(std::string outPath) {

  if (!Adjust_GCP_step2(outPath))
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Process_single_fix_fuv(std::string sMatchesDir) {
    //-------------------
    // Keep only the largest biedge connected subgraph
    //-------------------
    {
      const Pair_Set pairs = matches_provider_->getPairs();
      const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
      if(set_remainingIds.empty())
      {
        std::cout << "Invalid input image graph for global SfM" << std::endl;
        return false;
      }
      KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
    }

    openMVG::rotation_averaging::RelativeRotations relatives_R;
    Compute_Relative_Rotations(relatives_R);

    Hash_Map<IndexT, Mat3> global_rotations;
    if (!Compute_Global_Rotations(relatives_R, global_rotations))
    {
      std::cerr << "GlobalSfM:: Rotation Averaging failure!" << std::endl;
      return false;
    }
    matching::PairWiseMatches  tripletWise_matches;
    if (!Compute_Global_Translations(global_rotations, tripletWise_matches))
    {
      std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
      return false;
    }
    if (!Compute_Initial_Structure(tripletWise_matches))
    {
      std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
      return false;
    }
    if (!Adjust_single_fix_fuv())
    {
      std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
      return false;
    }

    //-- Export statistics about the SfM process
    if (!sLogging_file_.empty())
    {
      using namespace htmlDocument;
      std::ostringstream os;
      os << "Structure from Motion statistics.";
      html_doc_stream_->pushInfo("<hr>");
      html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

      os.str("");
      os << "-------------------------------" << "<br>"
        << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
        << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
        << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
        << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
        << "-------------------------------" << "<br>";
      html_doc_stream_->pushInfo(os.str());
    }

    return true;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process_down_imu_gps(std::string sMatchesDir) {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if(set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

//    openMVG::rotation_averaging::RelativeRotations relatives_R;
//    Compute_Relative_Rotations(relatives_R);

    Hash_Map<IndexT, Mat3> global_rotations;
//    if (!Compute_Global_Rotations(relatives_R, global_rotations))
//    {
//      std::cerr << "GlobalSfM:: Rotation Averaging failure!" << std::endl;
//      return false;
//    }

    const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};
    const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};

        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;


        //get gps and imu data
      {
    //        /Hash_Map<IndexT, std::vector<double> > R_imu;
            Hash_Map<IndexT, std::vector<double> > C_gps;
            {
                //read file and initialise external
                std::vector<double> c(3);
                std::string tFilePath = sfm_data_.s_root_path + "/latlon_c.res";  //!!!
                std::ifstream in;
                in.open(tFilePath);
                if(!in)
                {
                    std::cout << "open tFilePath file error!" << std::endl;
                    return false;
                }
                std::size_t line = 0;
                std::size_t count = sfm_data_.views.size();
                while((!in.eof()) && (line < count))
                {
                     in >> c[0] >> c[1] >> c[2];
                     C_gps[line] = {c[0], c[1], c[2]};
                    ++line;
                }
                in.close();


            }

        //read imu
        {
                std::ifstream in;
                in.open(sfm_data_.s_root_path + "/IMU.res"); //imu file

                int line = 0;
                const int count = sfm_data_.views.size();

        //        std::ofstream out_show_camera;
        //        out_show_camera.open(sfm_data_.s_root_path + "/show_imu_R.txt");

                while((!in.eof()) && (line < count))
                {
                    double roll, pitch, yaw;
                    in >> roll >> pitch >> yaw;  //z x y

                    Mat3 Rz, Ry, Rx;
                    Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
                          sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
                          0.0, 0.0, 1.0;
                    Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
                          0.0, 1.0, 0.0,
                          -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
                    Rx << 1.0, 0.0, 0.0,
                          0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
                          0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);
                    Mat3 R_imu =  Rz*Rx*Ry;
                    Mat3 changeAxis_M;
                    changeAxis_M << 1.0, 0.0, 0.0,
                                    0.0, -1.0, 0.0,
                                    0.0, 0.0, -1.0;
                    R_imu = changeAxis_M * R_imu;



//                    Mat3 R_imu_trans;
//                    ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
//                    Vec3 gps_trans = Vec3(tramsformation_x_gps[0], tramsformation_x_gps[1], tramsformation_x_gps[2]);

//                    //set camera down
//                    Mat3 Rd_res = R_imu_trans * R_imu;
                    Vec3 td_res, temp_C_gps;
                    temp_C_gps << C_gps[line][0], C_gps[line][1], C_gps[line][2];
//                    td_res = gps_trans - R_imu_trans * R_imu * temp_C_gps;
                    global_rotations[line] = R_imu;//Rd_res;
//                    sfm_data_.poses[line] = Pose3(Rd_res, -Rd_res.transpose() * td_res);

                    sfm_data_.poses[line] = Pose3(R_imu, temp_C_gps);

                    ++line;
                }
                in.close();


            }
//            //set pose
//            {

//                Poses::iterator itPose = sfm_data_.poses.begin();

//                for(std::size_t i = 0; i < sfm_data_.views.size(); ++i)
//                {
//                    std::advance(itPose,i);
//                    Vec3 c;
//                    c << C_gps[i][0], C_gps[i][1], C_gps[i][2];
//                    sfm_data_.poses[i] = Pose3(global_rotations[i], c);

//                }

//            }
      }


  matching::PairWiseMatches tripletWise_matches;
  if (!Compute_Global_Translations_gps(global_rotations, tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
    return false;
  }

    if (!Compute_Initial_Structure(tripletWise_matches))
    {
      std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
      return false;
    }

//    for( int i = 0; i < sfm_data_.poses.size(); ++i)
//    {
//        double R0_angleAxis[3];
//        Mat3 mid_R0 = sfm_data_.poses[i].rotation();
//        ceres::RotationMatrixToAngleAxis((const double*)mid_R0.data(), R0_angleAxis);

//        std::cout << "R0_angleAxis : " << R0_angleAxis[0] << " " << R0_angleAxis[1] << " " << R0_angleAxis[2] << std::endl;

//    }

    {
//            std::cout<< "initial intrinsics :"<<std::endl;
//                Hash_Map<IndexT, std::vector<double> > map_intrinsics;
//                map_intrinsics[0] = sfm_data_.intrinsics[0]->getParams();
//                std::cout<< "map_intrinsics 0 " << map_intrinsics[0][0]<<" "<<map_intrinsics[0][1]<<" "<<
//                         map_intrinsics[0][2]<<" "<<map_intrinsics[0][3]<<" "<<map_intrinsics[0][4]<<" "<<
//                         map_intrinsics[0][5]<<" "<<map_intrinsics[0][6]<<" "<<map_intrinsics[0][7]<<std::endl;

//                map_intrinsics[1] = sfm_data_.intrinsics[1]->getParams();
//                std::cout<< "map_intrinsics 1 " << map_intrinsics[1][0]<<" "<<map_intrinsics[1][1]<<" "<<
//                         map_intrinsics[1][2]<<" "<<map_intrinsics[1][3]<<" "<<map_intrinsics[1][4]<<" "<<
//                         map_intrinsics[1][5]<<" "<<map_intrinsics[1][6]<<" "<<map_intrinsics[1][7]<<std::endl;

//                map_intrinsics[2] = sfm_data_.intrinsics[2]->getParams();
//                std::cout<< "map_intrinsics 2 " << map_intrinsics[2][0]<<" "<<map_intrinsics[2][1]<<" "<<
//                         map_intrinsics[2][2]<<" "<<map_intrinsics[2][3]<<" "<<map_intrinsics[2][4]<<" "<<
//                         map_intrinsics[2][5]<<" "<<map_intrinsics[2][6]<<" "<<map_intrinsics[2][7]<<std::endl;
//                std::cout <<"-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() <<std::endl;
    }


  if (!Adjust_down_gps_change())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }




  {
//          std::cout<< "map_intrinsics 1"<<std::endl;
//              Hash_Map<IndexT, std::vector<double> > map_intrinsics;
//              map_intrinsics[0] = sfm_data_.intrinsics[0]->getParams();
//              std::cout<< "map_intrinsics 0 " << map_intrinsics[0][0]<<" "<<map_intrinsics[0][1]<<" "<<
//                       map_intrinsics[0][2]<<" "<<map_intrinsics[0][3]<<" "<<map_intrinsics[0][4]<<" "<<
//                       map_intrinsics[0][5]<<" "<<map_intrinsics[0][6]<<" "<<map_intrinsics[0][7]<<std::endl;

//              map_intrinsics[1] = sfm_data_.intrinsics[1]->getParams();
//              std::cout<< "map_intrinsics 1 " << map_intrinsics[1][0]<<" "<<map_intrinsics[1][1]<<" "<<
//                       map_intrinsics[1][2]<<" "<<map_intrinsics[1][3]<<" "<<map_intrinsics[1][4]<<" "<<
//                       map_intrinsics[1][5]<<" "<<map_intrinsics[1][6]<<" "<<map_intrinsics[1][7]<<std::endl;

//              map_intrinsics[2] = sfm_data_.intrinsics[2]->getParams();
//              std::cout<< "map_intrinsics 2 " << map_intrinsics[2][0]<<" "<<map_intrinsics[2][1]<<" "<<
//                       map_intrinsics[2][2]<<" "<<map_intrinsics[2][3]<<" "<<map_intrinsics[2][4]<<" "<<
//                       map_intrinsics[2][5]<<" "<<map_intrinsics[2][6]<<" "<<map_intrinsics[2][7]<<std::endl;
//              std::cout <<"-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() <<std::endl;
  }


  //out
  {
      std::ofstream out_show_camera;
      out_show_camera.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_reference.txt");
      std::ofstream out_Rc;
      out_Rc.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_Rc.txt");
      for(Poses::const_iterator it = sfm_data_.poses.begin();
          it != sfm_data_.poses.end(); ++it)
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

          out_Rc << c(0) << " " << c(1)<< " " << c(2) << " " << std::endl;

      }
      out_Rc.close();
      out_show_camera.close();

  }

  //out
  {
//      std::ofstream out_c;
//      out_c.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_c.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Vec3 c = pose.center();

//          out_c << c(0) << " " << c(1)<< " " << c(2) << std::endl;

//      }
//      out_c.close();

  }

  {
      std::ofstream out_show_camera;
      out_show_camera.open("/media/add7/E/run_scene/testMerge/29/out_camera_res_R.txt");
      for(Poses::const_iterator it = sfm_data_.poses.begin();
          it != sfm_data_.poses.end(); ++it)
      {
          const Pose3 & pose = it->second;
          const Mat3 R = pose.rotation();
          const Vec3 c = pose.center();

          Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
          axis_x = R.transpose() * (axis_x + R*c);
          axis_y = R.transpose() * (axis_y + R*c);
          axis_z = R.transpose() * (axis_z + R*c);
          out_show_camera <<R(0,0)<<" "<< R(0,1)<< " "<<R(0,2)<< " "
                          <<R(1,0)<<" "<< R(1,1)<< " "<<R(1,2)<< " "
                          <<R(2,0)<<" "<< R(2,1)<< " "<<R(2,2)<< std::endl;
      }
      out_show_camera.close();
  }





  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }





  return true;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process_threeCameras_imu_gps(std::string sMatchesDir) {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if(set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

//    openMVG::rotation_averaging::RelativeRotations relatives_R;
//    Compute_Relative_Rotations(relatives_R);

    Hash_Map<IndexT, Mat3> global_rotations;
//    if (!Compute_Global_Rotations(relatives_R, global_rotations))
//    {
//      std::cerr << "GlobalSfM:: Rotation Averaging failure!" << std::endl;
//      return false;
//    }


    //test distortion
    Mat3 r_tb, r_tf;
    r_tb << 1.0, 0.0, 0.0,
            0.0, 0.787926, -0.61577,
            0.0, 0.61577, 0.787926;
    r_tf << 1.0, 0.0, 0.0,
            0.0, 0.78796,  0.615727,
            0.0, -0.615727, 0.78796;
    Vec3 RA_tb, RA_tf;
    ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
    ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

        //set reference
        double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), -RA_tb(2)};//{0.665224, 0.00448363, 0.00291885};  //{0.703355, 0.00263495, -0.00149232};//{0.0, 0.0, 0.0}; //{0.0212232, 0.0206259, 0.0420509};  //{0.0, 0.0, 0.0}; //{0.669661, 0.00354668, 0.00220508}; //  //{0.0, 0.0, 0.0}; //{0.662013, -0.00281273, 0.000638113 }; //{0.0, 0.0, 0.0}; //{0.671924, 0.00432461, 0.000528427};
        double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), -RA_tf(2)};//{-0.644463, 0.0113345, -0.00156344}; //{-0.682501, 0.0147873, 0.00150509};//{0.0, 0.0, 0.0}; //{-0.0232707, -0.00288052, 0.00035795};  //{0.0, 0.0, 0.0}; //{-0.650191, 0.0101123, -0.000844533}; //  //{0.0, 0.0, 0.0}; //{-0.658446, 0.00532591, 0.0010695}; //{-0.656461, 0.00508743, -0.00299628};
        double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};

        Vec3 c_bd, c_fd;//b-d, f-d
        c_bd << 0.0, 0.098, 0.027;
        c_fd << 0.0, -0.101, 0.028;

        Vec3 t_db, t_df, t_dd;
        t_db = -r_tb*c_bd;
        t_df = -r_tf*c_fd;


        //double S = 1470.531783220221;
//        t_db << 0.0, 0.098, 0.027;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
//        t_df << 0.0, -0.101, 0.028;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
        t_dd << 0.0, 0.0, 0.0;


        std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
        transformation_br = {R_br_angleAxis[0], R_br_angleAxis[1], R_br_angleAxis[2],
                             t_db(0), t_db(1), t_db(2)};
        transformation_fr = {R_fr_angleAxis[0], R_fr_angleAxis[1], R_fr_angleAxis[2],
                             t_df(0), t_df(1), t_df(2)};
        transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                             t_dd(0), t_dd(1), t_dd(2)};

    const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
    const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};


        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;


        //get gps and imu data
      {
    //        /Hash_Map<IndexT, std::vector<double> > R_imu;
            Hash_Map<IndexT, std::vector<double> > C_gps;
            {
                //read file and initialise external
                std::vector<double> c(3);
                std::string tFilePath = sfm_data_.s_root_path + "/latlon_c.res";  //!!!
                std::ifstream in;
                in.open(tFilePath);
                if(!in)
                {
                    std::cout << "open tFilePath file error!" << std::endl;
                    return false;
                }
                std::size_t line = 0;
                const int count = sfm_data_.views.size()/3;
                while((!in.eof()) && (line < count))
                {
                     in >> c[0] >> c[1] >> c[2];
                     C_gps[line + count] = {c[0], c[1], c[2]};
                    ++line;
                }
                in.close();

            }
            //read imu
            {
            std::ifstream in;
            in.open(sfm_data_.s_root_path + "/IMU.res"); //imu file

            int line = 0;
            const int count = sfm_data_.views.size()/3;

    //        std::ofstream out_show_camera;
    //        out_show_camera.open(sfm_data_.s_root_path + "/show_imu_R.txt");

            while((!in.eof()) && (line < count))
            {
                double roll, pitch, yaw;
                in >> roll >> pitch >> yaw;  //z x y

                Mat3 Rz, Ry, Rx;
                Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
                      sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
                      0.0, 0.0, 1.0;
                Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
                      0.0, 1.0, 0.0,
                      -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
                Rx << 1.0, 0.0, 0.0,
                      0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
                      0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);
                Mat3 R_imu =  Rz*Rx*Ry;
                Mat3 changeAxis_M;
                changeAxis_M << 1.0, 0.0, 0.0,
                                0.0, -1.0, 0.0,
                                0.0, 0.0, -1.0;
                R_imu = changeAxis_M * R_imu;


                {

    //                Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
    //                Vec3 c = {C_gps[line + count][0], C_gps[line + count][1], C_gps[line + count][2]};
    //                axis_x = R_d.transpose() * (axis_x + R_d*c);
    //                axis_y = R_d.transpose() * (axis_y + R_d*c);
    //                axis_z = R_d.transpose() * (axis_z + R_d*c);
    //                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
    //                       << axis_x(0) << " " << axis_x(1) << " " << axis_x(2)<< ";" << std::endl;
    //                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
    //                       << axis_y(0) << " " << axis_y(1) << " " << axis_y(2)<< ";" << std::endl;
    //                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
    //                       << axis_z(0) << " " << axis_z(1) << " " << axis_z(2)<< ";" << std::endl;

    //                out_show_camera << R_d << std::endl << std::endl;
                }

    //            global_rotations[line + count] = R_d;
    //            sfm_data_.poses[line + count] = Pose3(R_d, C_gps[line + count]);

                Mat3 Rb_trans, Rf_trans, Rd_trans;
                Vec3 tb_trans, tf_trans, td_trans;

                ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rb_trans.data());
                tb_trans = Vec3(transformation_br[3], transformation_br[4], transformation_br[5]);
                ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rf_trans.data());
                tf_trans = Vec3(transformation_fr[3], transformation_fr[4], transformation_fr[5]);
                ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rd_trans.data());
                td_trans = Vec3(transformation_dr[3], transformation_dr[4], transformation_dr[5]);


                Mat3 R_imu_trans;
                ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
                Vec3 gps_trans = Vec3(tramsformation_x_gps[0], tramsformation_x_gps[1], tramsformation_x_gps[2]);
                Vec3 C_gps_res = Vec3(C_gps[line + count][0], C_gps[line + count][1], C_gps[line + count][2]);

                //set f
                Mat3 Rf_res = Rf_trans * R_imu_trans * R_imu;
                Vec3 tf_res;
                tf_res = Rf_trans * gps_trans - Rf_trans * R_imu_trans * R_imu * C_gps_res + tf_trans;
                global_rotations[line + count*2] = Rf_res;
                sfm_data_.poses[line + count*2] = Pose3(Rf_res, -Rf_res.transpose() * tf_res);

                //set d
                Mat3 Rd_res = Rd_trans * R_imu_trans * R_imu;
                Vec3 td_res;
                td_res = Rd_trans * gps_trans - Rd_trans * R_imu_trans * R_imu * C_gps_res + td_trans;
                global_rotations[line + count] = Rd_res;
                sfm_data_.poses[line + count] = Pose3(Rd_res, -Rd_res.transpose() * td_res);

                //set b
                Mat3 Rb_res = Rb_trans * R_imu_trans * R_imu;
                Vec3 tb_res;
                tb_res = Rb_trans * gps_trans - Rb_trans * R_imu_trans * R_imu * C_gps_res + tb_trans;
                global_rotations[line] = Rb_res;
                sfm_data_.poses[line] = Pose3(Rb_res, -Rb_res.transpose() * tb_res);

                ++line;
            }
            in.close();

    //        out_show_camera.close();
        }
      }


  matching::PairWiseMatches tripletWise_matches;
  if (!Compute_Global_Translations_gps(global_rotations, tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
    return false;
  }

    if (!Compute_Initial_Structure(tripletWise_matches))
    {
      std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
      return false;
    }

//    for( int i = 0; i < sfm_data_.poses.size(); ++i)
//    {
//        double R0_angleAxis[3];
//        Mat3 mid_R0 = sfm_data_.poses[i].rotation();
//        ceres::RotationMatrixToAngleAxis((const double*)mid_R0.data(), R0_angleAxis);

//        std::cout << "R0_angleAxis : " << R0_angleAxis[0] << " " << R0_angleAxis[1] << " " << R0_angleAxis[2] << std::endl;

//    }

    {
//            std::cout<< "initial intrinsics :"<<std::endl;
//                Hash_Map<IndexT, std::vector<double> > map_intrinsics;
//                map_intrinsics[0] = sfm_data_.intrinsics[0]->getParams();
//                std::cout<< "map_intrinsics 0 " << map_intrinsics[0][0]<<" "<<map_intrinsics[0][1]<<" "<<
//                         map_intrinsics[0][2]<<" "<<map_intrinsics[0][3]<<" "<<map_intrinsics[0][4]<<" "<<
//                         map_intrinsics[0][5]<<" "<<map_intrinsics[0][6]<<" "<<map_intrinsics[0][7]<<std::endl;

//                map_intrinsics[1] = sfm_data_.intrinsics[1]->getParams();
//                std::cout<< "map_intrinsics 1 " << map_intrinsics[1][0]<<" "<<map_intrinsics[1][1]<<" "<<
//                         map_intrinsics[1][2]<<" "<<map_intrinsics[1][3]<<" "<<map_intrinsics[1][4]<<" "<<
//                         map_intrinsics[1][5]<<" "<<map_intrinsics[1][6]<<" "<<map_intrinsics[1][7]<<std::endl;

//                map_intrinsics[2] = sfm_data_.intrinsics[2]->getParams();
//                std::cout<< "map_intrinsics 2 " << map_intrinsics[2][0]<<" "<<map_intrinsics[2][1]<<" "<<
//                         map_intrinsics[2][2]<<" "<<map_intrinsics[2][3]<<" "<<map_intrinsics[2][4]<<" "<<
//                         map_intrinsics[2][5]<<" "<<map_intrinsics[2][6]<<" "<<map_intrinsics[2][7]<<std::endl;
//                std::cout <<"-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() <<std::endl;
    }


  if (!Adjust_threecamera_gps_change())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }




  {
//          std::cout<< "map_intrinsics 1"<<std::endl;
//              Hash_Map<IndexT, std::vector<double> > map_intrinsics;
//              map_intrinsics[0] = sfm_data_.intrinsics[0]->getParams();
//              std::cout<< "map_intrinsics 0 " << map_intrinsics[0][0]<<" "<<map_intrinsics[0][1]<<" "<<
//                       map_intrinsics[0][2]<<" "<<map_intrinsics[0][3]<<" "<<map_intrinsics[0][4]<<" "<<
//                       map_intrinsics[0][5]<<" "<<map_intrinsics[0][6]<<" "<<map_intrinsics[0][7]<<std::endl;

//              map_intrinsics[1] = sfm_data_.intrinsics[1]->getParams();
//              std::cout<< "map_intrinsics 1 " << map_intrinsics[1][0]<<" "<<map_intrinsics[1][1]<<" "<<
//                       map_intrinsics[1][2]<<" "<<map_intrinsics[1][3]<<" "<<map_intrinsics[1][4]<<" "<<
//                       map_intrinsics[1][5]<<" "<<map_intrinsics[1][6]<<" "<<map_intrinsics[1][7]<<std::endl;

//              map_intrinsics[2] = sfm_data_.intrinsics[2]->getParams();
//              std::cout<< "map_intrinsics 2 " << map_intrinsics[2][0]<<" "<<map_intrinsics[2][1]<<" "<<
//                       map_intrinsics[2][2]<<" "<<map_intrinsics[2][3]<<" "<<map_intrinsics[2][4]<<" "<<
//                       map_intrinsics[2][5]<<" "<<map_intrinsics[2][6]<<" "<<map_intrinsics[2][7]<<std::endl;
//              std::cout <<"-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() <<std::endl;
  }


  //out
  {
//      std::ofstream out_show_camera;
//      out_show_camera.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_reference.txt");
//      std::ofstream out_Rc;
//      out_Rc.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_Rc.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 c = pose.center();

//          Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
//          axis_x = R.transpose() * (axis_x + R*c);
//          axis_y = R.transpose() * (axis_y + R*c);
//          axis_z = R.transpose() * (axis_z + R*c);
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_x(0) << " " << axis_x(1) << " " << axis_x(2)<< ";" << std::endl;
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_y(0) << " " << axis_y(1) << " " << axis_y(2)<< ";" << std::endl;
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_z(0) << " " << axis_z(1) << " " << axis_z(2)<< ";" << std::endl;

//          out_Rc << c(0) << " " << c(1)<< " " << c(2) << " " << std::endl;

//      }
//      out_Rc.close();
//      out_show_camera.close();

  }

  //out
  {
//      std::ofstream out_c;
//      out_c.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_c.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Vec3 c = pose.center();

//          out_c << c(0) << " " << c(1)<< " " << c(2) << std::endl;

//      }
//      out_c.close();

  }

  {
//      std::ofstream out_show_camera;
//      out_show_camera.open("/media/add7/E/run_scene/testMerge/29/out_camera_res_R.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 c = pose.center();

//          Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
//          axis_x = R.transpose() * (axis_x + R*c);
//          axis_y = R.transpose() * (axis_y + R*c);
//          axis_z = R.transpose() * (axis_z + R*c);
//          out_show_camera <<R(0,0)<<" "<< R(0,1)<< " "<<R(0,2)<< " "
//                          <<R(1,0)<<" "<< R(1,1)<< " "<<R(1,2)<< " "
//                          <<R(2,0)<<" "<< R(2,1)<< " "<<R(2,2)<< std::endl;
//      }
//      out_show_camera.close();
  }





  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }





  return true;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process_threeCameras_imu_gps_cail_new() {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if(set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

  Hash_Map<IndexT, Mat3> global_rotations;

  //init inter-constraint
  Mat3 Rb_trans, Rf_trans, Rd_trans;
  Vec3 tb_trans, tf_trans, td_trans;
  Mat3 R_imu_trans;
  Vec3 gps_trans;
  std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
  const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
  const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};
  {
      //test distortion
      Mat3 r_tb, r_tf;
      r_tb << 1.0, 0.0, 0.0,
              0.0, 0.787926, -0.61577,
              0.0, 0.61577, 0.787926;
      r_tf << 1.0, 0.0, 0.0,
              0.0, 0.78796,  0.615727,
              0.0, -0.615727, 0.78796;
      Vec3 RA_tb, RA_tf;
      ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
      ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

      //set reference
      double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), -RA_tb(2)};
      double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), -RA_tf(2)};
      double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};

      Vec3 c_bd, c_fd;//b-d, f-d
      c_bd << 0.0, 0.098, 0.027;
      c_fd << 0.0, -0.101, 0.028;

      Vec3 t_db, t_df, t_dd;
      t_db = -r_tb*c_bd;
      t_df = -r_tf*c_fd;

      //double S = 1470.531783220221;
    //        t_db << 0.0, 0.098, 0.027;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
    //        t_df << 0.0, -0.101, 0.028;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
      t_dd << 0.0, 0.0, 0.0;

      transformation_br = {R_br_angleAxis[0], R_br_angleAxis[1], R_br_angleAxis[2],
                           t_db(0), t_db(1), t_db(2)};
      transformation_fr = {R_fr_angleAxis[0], R_fr_angleAxis[1], R_fr_angleAxis[2],
                           t_df(0), t_df(1), t_df(2)};
      transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                           t_dd(0), t_dd(1), t_dd(2)};

      const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
      const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};

      std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                             << tramsformation_x_gps[1] << " "
                                             << tramsformation_x_gps[2] << std::endl;

      //prepare for setting pose
      ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rb_trans.data());
      tb_trans = Vec3(transformation_br[3], transformation_br[4], transformation_br[5]);
      ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rf_trans.data());
      tf_trans = Vec3(transformation_fr[3], transformation_fr[4], transformation_fr[5]);
      ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rd_trans.data());
      td_trans = Vec3(transformation_dr[3], transformation_dr[4], transformation_dr[5]);


      ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
      gps_trans = Vec3(tramsformation_x_gps[0], tramsformation_x_gps[1], tramsformation_x_gps[2]);

  }

  //get gps and imu data
  PoseStatusMap poseStatusMap;
  Hash_Map<IndexT, std::vector<double> > C_gps_Map;
  Hash_Map<IndexT, std::vector<double> > R_imu_Map;
  {
    //        /Hash_Map<IndexT, std::vector<double> > R_imu;
//            Hash_Map<IndexT, std::vector<double> > C_gps;
      std::string tFilePath = sfm_data_.s_root_path + "/n_gps_imu_321.res";  //!!!
      std::ifstream in;
      in.open(tFilePath);
      if(!in)
      {
          std::cout << "open tFilePath file error!" << std::endl;
          return false;
      }
      std::size_t line = 0;
      int count = -1;
      in >> count;  ///!!!

      bool isC0 = false;
      std::vector<double> c0(3);
      while((!in.eof()) && (line < count))
      {
          // input
          int imgId, imgStatus;
          std::vector<double> c(3);
          //double pitch, roll, yaw;
          double omega, phi, kappa;
          in >> imgId
             >> c[0] >> c[1] >> c[2]
             >> omega >> phi >> kappa
             >> imgStatus;

          if(isC0 == false)
          {
              c0 = c;
              isC0 = true;
          }

          c[0] -= c0[0];
          c[1] -= c0[1];
          c[2] -= c0[2];

          ++line;

          Mat3 Rz, Ry, Rx;
          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
                0.0, 0.0, 1.0;
          Ry << cos(phi/180.0*M_PI), 0.0, -sin(phi/180.0*M_PI),
                0.0, 1.0, 0.0,
                sin(phi/180.0*M_PI), 0.0, cos(phi/180.0*M_PI);
          Rx << 1.0, 0.0, 0.0,
                0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
                0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
          Mat3 R_imu =  Ry*Rx*Rz;

          Mat3 changeAxis_M;
          changeAxis_M << 1.0, 0.0, 0.0,
                          0.0, -1.0, 0.0,
                          0.0, 0.0, -1.0;
          R_imu = changeAxis_M * R_imu;

          std::vector<double> R_imu_a(3);
          ceres::RotationMatrixToAngleAxis(R_imu.data(), R_imu_a.data());

          Vec3 C_gps_res = Vec3(c[0], c[1], c[2]);

          // check status
          if(imgId/100000 == 1)//back
          {
              if(imgStatus != 321)
              {
                  C_gps_Map[imgId + 100000] = c;
                  R_imu_Map[imgId + 100000] = R_imu_a;
              }

              //set back-camera pose
              Mat3 Rb_res = Rb_trans * R_imu_trans * R_imu;
              Vec3 tb_res;
              tb_res = Rb_trans * gps_trans - Rb_trans * R_imu_trans * R_imu * C_gps_res + tb_trans;
              global_rotations[imgId] = Rb_res;
              sfm_data_.poses[imgId] = Pose3(Rb_res, -Rb_res.transpose() * tb_res);

              PoseStatus poseStatus;
              poseStatus.pose = sfm_data_.poses[imgId];
              poseStatus.status = imgStatus;
              poseStatus.bdf = 1;
              poseStatusMap[imgId] = poseStatus;

          }else if(imgId/100000 == 2)//down
          {
              C_gps_Map[imgId] = c;
              R_imu_Map[imgId] = R_imu_a;

              //set down-camera pose
              Mat3 Rd_res = Rd_trans * R_imu_trans * R_imu;
              Vec3 td_res;
              td_res = Rd_trans * gps_trans - Rd_trans * R_imu_trans * R_imu * C_gps_res + td_trans;
              global_rotations[imgId] = Rd_res;
              sfm_data_.poses[imgId] = Pose3(Rd_res, -Rd_res.transpose() * td_res);

              PoseStatus poseStatus;
              poseStatus.pose = sfm_data_.poses[imgId];
              poseStatus.status = imgStatus;
              poseStatus.bdf = 2;
              poseStatusMap[imgId] = poseStatus;

          }else if(imgId/100000 == 3)//front
          {
              if(imgStatus != 321)
              {
                  C_gps_Map[imgId - 100000] = c;
                  R_imu_Map[imgId - 100000] = R_imu_a;
              }

              //set front-camera pose
              Mat3 Rf_res = Rf_trans * R_imu_trans * R_imu;
              Vec3 tf_res;
              tf_res = Rf_trans * gps_trans - Rf_trans * R_imu_trans * R_imu * C_gps_res + tf_trans;
              global_rotations[imgId] = Rf_res;
              sfm_data_.poses[imgId] = Pose3(Rf_res, -Rf_res.transpose() * tf_res);

              PoseStatus poseStatus;
              poseStatus.pose = sfm_data_.poses[imgId];
              poseStatus.status = imgStatus;
              poseStatus.bdf = 3;
              poseStatusMap[imgId] = poseStatus;

          }

      }
      in.close();


            {
//                //read file and initialise external
//                std::vector<double> c(3);
//                std::string tFilePath = sfm_data_.s_root_path + "/latlon_c.res";  //!!!
//                std::ifstream in;
//                in.open(tFilePath);
//                if(!in)
//                {
//                    std::cout << "open tFilePath file error!" << std::endl;
//                    return false;
//                }
//                std::size_t line = 0;
//                const int count = sfm_data_.views.size()/3;
//                while((!in.eof()) && (line < count))
//                {
//                     in >> c[0] >> c[1] >> c[2];
//                     C_gps[line + count] = {c[0], c[1], c[2]};
//                    ++line;
//                }
//                in.close();

            }
            //read imu
            {
//            std::ifstream in;
//            in.open(sfm_data_.s_root_path + "/IMU.res"); //imu file

//            int line = 0;
//            const int count = sfm_data_.views.size()/3;

//    //        std::ofstream out_show_camera;
//    //        out_show_camera.open(sfm_data_.s_root_path + "/show_imu_R.txt");

//            while((!in.eof()) && (line < count))
//            {
//                double roll, pitch, yaw;
//                in >> roll >> pitch >> yaw;  //z x y

//                Mat3 Rz, Ry, Rx;
//                Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
//                      sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
//                      0.0, 0.0, 1.0;
//                Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
//                      0.0, 1.0, 0.0,
//                      -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
//                Rx << 1.0, 0.0, 0.0,
//                      0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
//                      0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);
//                Mat3 R_imu =  Rz*Rx*Ry;
//                Mat3 changeAxis_M;
//                changeAxis_M << 1.0, 0.0, 0.0,
//                                0.0, -1.0, 0.0,
//                                0.0, 0.0, -1.0;
//                R_imu = changeAxis_M * R_imu;


//                {

//    //                Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
//    //                Vec3 c = {C_gps[line + count][0], C_gps[line + count][1], C_gps[line + count][2]};
//    //                axis_x = R_d.transpose() * (axis_x + R_d*c);
//    //                axis_y = R_d.transpose() * (axis_y + R_d*c);
//    //                axis_z = R_d.transpose() * (axis_z + R_d*c);
//    //                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//    //                       << axis_x(0) << " " << axis_x(1) << " " << axis_x(2)<< ";" << std::endl;
//    //                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//    //                       << axis_y(0) << " " << axis_y(1) << " " << axis_y(2)<< ";" << std::endl;
//    //                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//    //                       << axis_z(0) << " " << axis_z(1) << " " << axis_z(2)<< ";" << std::endl;

//    //                out_show_camera << R_d << std::endl << std::endl;
//                }

//    //            global_rotations[line + count] = R_d;
//    //            sfm_data_.poses[line + count] = Pose3(R_d, C_gps[line + count]);

//                Mat3 Rb_trans, Rf_trans, Rd_trans;
//                Vec3 tb_trans, tf_trans, td_trans;

//                ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rb_trans.data());
//                tb_trans = Vec3(transformation_br[3], transformation_br[4], transformation_br[5]);
//                ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rf_trans.data());
//                tf_trans = Vec3(transformation_fr[3], transformation_fr[4], transformation_fr[5]);
//                ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rd_trans.data());
//                td_trans = Vec3(transformation_dr[3], transformation_dr[4], transformation_dr[5]);


//                Mat3 R_imu_trans;
//                ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
//                Vec3 gps_trans = Vec3(tramsformation_x_gps[0], tramsformation_x_gps[1], tramsformation_x_gps[2]);
//                Vec3 C_gps_res = Vec3(C_gps[line + count][0], C_gps[line + count][1], C_gps[line + count][2]);

//                //set f
//                Mat3 Rf_res = Rf_trans * R_imu_trans * R_imu;
//                Vec3 tf_res;
//                tf_res = Rf_trans * gps_trans - Rf_trans * R_imu_trans * R_imu * C_gps_res + tf_trans;
//                global_rotations[line + count*2] = Rf_res;
//                sfm_data_.poses[line + count*2] = Pose3(Rf_res, -Rf_res.transpose() * tf_res);

//                //set d
//                Mat3 Rd_res = Rd_trans * R_imu_trans * R_imu;
//                Vec3 td_res;
//                td_res = Rd_trans * gps_trans - Rd_trans * R_imu_trans * R_imu * C_gps_res + td_trans;
//                global_rotations[line + count] = Rd_res;
//                sfm_data_.poses[line + count] = Pose3(Rd_res, -Rd_res.transpose() * td_res);

//                //set b
//                Mat3 Rb_res = Rb_trans * R_imu_trans * R_imu;
//                Vec3 tb_res;
//                tb_res = Rb_trans * gps_trans - Rb_trans * R_imu_trans * R_imu * C_gps_res + tb_trans;
//                global_rotations[line] = Rb_res;
//                sfm_data_.poses[line] = Pose3(Rb_res, -Rb_res.transpose() * tb_res);

//                ++line;
//            }
//            in.close();

//    //        out_show_camera.close();
        }
  }

  matching::PairWiseMatches tripletWise_matches;
  if (!Compute_Global_Translations_gps(global_rotations, tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
    return false;
  }

  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }

//    for( int i = 0; i < sfm_data_.poses.size(); ++i)
//    {
//        double R0_angleAxis[3];
//        Mat3 mid_R0 = sfm_data_.poses[i].rotation();
//        ceres::RotationMatrixToAngleAxis((const double*)mid_R0.data(), R0_angleAxis);

//        std::cout << "R0_angleAxis : " << R0_angleAxis[0] << " " << R0_angleAxis[1] << " " << R0_angleAxis[2] << std::endl;

//    }

    {
//            std::cout<< "initial intrinsics :"<<std::endl;
//                Hash_Map<IndexT, std::vector<double> > map_intrinsics;
//                map_intrinsics[0] = sfm_data_.intrinsics[0]->getParams();
//                std::cout<< "map_intrinsics 0 " << map_intrinsics[0][0]<<" "<<map_intrinsics[0][1]<<" "<<
//                         map_intrinsics[0][2]<<" "<<map_intrinsics[0][3]<<" "<<map_intrinsics[0][4]<<" "<<
//                         map_intrinsics[0][5]<<" "<<map_intrinsics[0][6]<<" "<<map_intrinsics[0][7]<<std::endl;

//                map_intrinsics[1] = sfm_data_.intrinsics[1]->getParams();
//                std::cout<< "map_intrinsics 1 " << map_intrinsics[1][0]<<" "<<map_intrinsics[1][1]<<" "<<
//                         map_intrinsics[1][2]<<" "<<map_intrinsics[1][3]<<" "<<map_intrinsics[1][4]<<" "<<
//                         map_intrinsics[1][5]<<" "<<map_intrinsics[1][6]<<" "<<map_intrinsics[1][7]<<std::endl;

//                map_intrinsics[2] = sfm_data_.intrinsics[2]->getParams();
//                std::cout<< "map_intrinsics 2 " << map_intrinsics[2][0]<<" "<<map_intrinsics[2][1]<<" "<<
//                         map_intrinsics[2][2]<<" "<<map_intrinsics[2][3]<<" "<<map_intrinsics[2][4]<<" "<<
//                         map_intrinsics[2][5]<<" "<<map_intrinsics[2][6]<<" "<<map_intrinsics[2][7]<<std::endl;
//                std::cout <<"-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() <<std::endl;
    }


  if (!Adjust_threecamera_gps_cail_new(poseStatusMap, transformation_R_imu, tramsformation_x_gps, transformation_br, transformation_fr, C_gps_Map, R_imu_Map))
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }




  {
//          std::cout<< "map_intrinsics 1"<<std::endl;
//              Hash_Map<IndexT, std::vector<double> > map_intrinsics;
//              map_intrinsics[0] = sfm_data_.intrinsics[0]->getParams();
//              std::cout<< "map_intrinsics 0 " << map_intrinsics[0][0]<<" "<<map_intrinsics[0][1]<<" "<<
//                       map_intrinsics[0][2]<<" "<<map_intrinsics[0][3]<<" "<<map_intrinsics[0][4]<<" "<<
//                       map_intrinsics[0][5]<<" "<<map_intrinsics[0][6]<<" "<<map_intrinsics[0][7]<<std::endl;

//              map_intrinsics[1] = sfm_data_.intrinsics[1]->getParams();
//              std::cout<< "map_intrinsics 1 " << map_intrinsics[1][0]<<" "<<map_intrinsics[1][1]<<" "<<
//                       map_intrinsics[1][2]<<" "<<map_intrinsics[1][3]<<" "<<map_intrinsics[1][4]<<" "<<
//                       map_intrinsics[1][5]<<" "<<map_intrinsics[1][6]<<" "<<map_intrinsics[1][7]<<std::endl;

//              map_intrinsics[2] = sfm_data_.intrinsics[2]->getParams();
//              std::cout<< "map_intrinsics 2 " << map_intrinsics[2][0]<<" "<<map_intrinsics[2][1]<<" "<<
//                       map_intrinsics[2][2]<<" "<<map_intrinsics[2][3]<<" "<<map_intrinsics[2][4]<<" "<<
//                       map_intrinsics[2][5]<<" "<<map_intrinsics[2][6]<<" "<<map_intrinsics[2][7]<<std::endl;
//              std::cout <<"-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() <<std::endl;
  }


  //out
  {
//      std::ofstream out_show_camera;
//      out_show_camera.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_reference.txt");
//      std::ofstream out_Rc;
//      out_Rc.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_Rc.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 c = pose.center();

//          Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
//          axis_x = R.transpose() * (axis_x + R*c);
//          axis_y = R.transpose() * (axis_y + R*c);
//          axis_z = R.transpose() * (axis_z + R*c);
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_x(0) << " " << axis_x(1) << " " << axis_x(2)<< ";" << std::endl;
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_y(0) << " " << axis_y(1) << " " << axis_y(2)<< ";" << std::endl;
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_z(0) << " " << axis_z(1) << " " << axis_z(2)<< ";" << std::endl;

//          out_Rc << c(0) << " " << c(1)<< " " << c(2) << " " << std::endl;

//      }
//      out_Rc.close();
//      out_show_camera.close();

  }

  //out
  {
//      std::ofstream out_c;
//      out_c.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_c.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Vec3 c = pose.center();

//          out_c << c(0) << " " << c(1)<< " " << c(2) << std::endl;

//      }
//      out_c.close();

  }

  {
//      std::ofstream out_show_camera;
//      out_show_camera.open("/media/add7/E/run_scene/testMerge/29/out_camera_res_R.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 c = pose.center();

//          Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
//          axis_x = R.transpose() * (axis_x + R*c);
//          axis_y = R.transpose() * (axis_y + R*c);
//          axis_z = R.transpose() * (axis_z + R*c);
//          out_show_camera <<R(0,0)<<" "<< R(0,1)<< " "<<R(0,2)<< " "
//                          <<R(1,0)<<" "<< R(1,1)<< " "<<R(1,2)<< " "
//                          <<R(2,0)<<" "<< R(2,1)<< " "<<R(2,2)<< std::endl;
//      }
//      out_show_camera.close();
  }





  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }





  return true;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process_threeCameras_imu_gps_init_new(std::string rtFilePath) {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if(set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

  Hash_Map<IndexT, Mat3> global_rotations;

  //init inter-constraint
  Mat3 Rb_trans, Rf_trans, Rd_trans;
  Vec3 tb_trans, tf_trans, td_trans;
  Mat3 R_imu_trans;
  Vec3 gps_trans;
  std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
  const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
  const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};
  {
      //test distortion
      Mat3 r_tb, r_tf;
      r_tb << 1.0, 0.0, 0.0,
              0.0, 0.787926, -0.61577,
              0.0, 0.61577, 0.787926;
      r_tf << 1.0, 0.0, 0.0,
              0.0, 0.78796,  0.615727,
              0.0, -0.615727, 0.78796;

//      Rx << 1.0, 0.0, 0.0,
//             0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//             0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


//      r_tb << 1.0, 0.0, 0.0,
//              0.0, cos(38.0/180.0*M_PI), -sin(38.0/180.0*M_PI),
//              0.0, sin(38.0/180.0*M_PI), cos(38.0/180.0*M_PI);

//      r_tf << 1.0, 0.0, 0.0,
//              0.0, cos(-38.0/180.0*M_PI), -sin(-38.0/180.0*M_PI),
//              0.0, sin(-38.0/180.0*M_PI), cos(-38.0/180.0*M_PI);

//      r_tb << cos(-38.0/180.0*M_PI), 0.0, sin(-38.0/180.0*M_PI),//-38
//              0.0, 1.0, 0.0,
//              -sin(-38.0/180.0*M_PI), 0.0, cos(-38.0/180.0*M_PI);
//      r_tf << cos(38.0/180.0*M_PI), 0.0, sin(38.0/180.0*M_PI),//38
//              0.0, 1.0, 0.0,
//              -sin(38.0/180.0*M_PI), 0.0, cos(38.0/180.0*M_PI);
      Vec3 RA_tb, RA_tf;
      ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
      ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

      //set reference
      double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), -RA_tb(2)};
      double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), -RA_tf(2)};
//      double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), RA_tb(2)};
//      double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), RA_tf(2)};
      double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};

      Vec3 c_bd, c_fd;//b-d, f-d
      c_bd << 0.0, 0.098, 0.027;
      c_fd << 0.0, -0.101, 0.028;
//      c_bd << -0.098, 0.0, -0.027;
//      c_fd << 0.101, 0.0, -0.028;

//      c_bd << 0.098, 0.0, -0.027;
//      c_fd << -0.101, 0.0, -0.028;
//      c_bd << 0.0, 1.0, 0.027;
//      c_fd << 0.0, -0.101, 0.028;

      Vec3 t_db, t_df, t_dd;
      t_db = -r_tb*c_bd;
      t_df = -r_tf*c_fd;

      //double S = 1470.531783220221;
    //        t_db << 0.0, 0.098, 0.027;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
    //        t_df << 0.0, -0.101, 0.028;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
      t_dd << 0.0, 0.0, 0.0;

      transformation_br = {R_br_angleAxis[0], R_br_angleAxis[1], R_br_angleAxis[2],
                           t_db(0), t_db(1), t_db(2)};
      transformation_fr = {R_fr_angleAxis[0], R_fr_angleAxis[1], R_fr_angleAxis[2],
                           t_df(0), t_df(1), t_df(2)};
      transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                           t_dd(0), t_dd(1), t_dd(2)};

//      const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
//      const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};

      std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                             << tramsformation_x_gps[1] << " "
                                             << tramsformation_x_gps[2] << std::endl;

      //prepare for setting pose
      ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rb_trans.data());
      tb_trans = Vec3(transformation_br[3], transformation_br[4], transformation_br[5]);
      ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rf_trans.data());
      tf_trans = Vec3(transformation_fr[3], transformation_fr[4], transformation_fr[5]);
      ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rd_trans.data());
      td_trans = Vec3(transformation_dr[3], transformation_dr[4], transformation_dr[5]);


      ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
      gps_trans = Vec3(tramsformation_x_gps[0], tramsformation_x_gps[1], tramsformation_x_gps[2]);

  }


  //get gps and imu data
  PoseStatusMap poseStatusMap;
  Hash_Map<IndexT, std::vector<double> > C_gps_Map;
  Hash_Map<IndexT, std::vector<double> > R_imu_Map;
  std::vector<double> c0(3);
  {
    //        /Hash_Map<IndexT, std::vector<double> > R_imu;
//            Hash_Map<IndexT, std::vector<double> > C_gps;
      std::string tFilePath = rtFilePath;  //!!!
      std::ifstream in;
      in.open(tFilePath);
      if(!in)
      {
          std::cout << "open tFilePath file error!" << std::endl;
          return false;
      }
      std::size_t line = 0;
      int count = -1;
      in >> count;  ///!!!

      bool isC0 = false;//true;

      while((!in.eof()) && (line < count))
      {
          // input
          int imgId, imgStatus;
          std::vector<double> c(3);
          //double pitch, roll, yaw;
          double omega, phi, kappa;
          in >> imgId
             >> c[0] >> c[1] >> c[2]
//             >> c[1] >> c[0] >> c[2]
             >> omega >> phi >> kappa
             >> imgStatus;

//          omega = omega + 2.4;
//          phi = phi + 0.0713;
//          kappa = kappa - 0.5805;
//          kappa = kappa + 2.4;
//          phi = phi + 0.0713;
//          omega = omega - 0.5805;

          if(isC0 == false)
          {
              c0 = c;
              isC0 = true;
          }

//          c[0] -= c0[0];
//          c[1] -= c0[1];
//          c[2] -= c0[2];

          ++line;



//          double roll, pitch, yaw;

//          roll = omega;
//          pitch = phi;
//          yaw = kappa;

//          Mat3 Rz, Ry, Rx;
//          Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
//                sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
//                0.0, 1.0, 0.0,
//                -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
//          Rx << 1.0, 0.0, 0.0,
//                0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
//                0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);
//          Mat3 R_imu =  Rz*Rx*Ry;
//          Mat3 changeAxis_M;
//          changeAxis_M << 1.0, 0.0, 0.0,
//                          0.0, -1.0, 0.0,
//                          0.0, 0.0, -1.0;
//          R_imu = changeAxis_M * R_imu;




//          Mat3 Rz, Ry, Rx;
//          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
//                0.0, 1.0, 0.0,
//                -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);
//          Rx << 1.0, 0.0, 0.0,
//                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

////          Mat3 R_imu =  Rx*Ry*Rz;
////          Mat3 R_imu =  Rx*Rz*Ry;//
////          Mat3 R_imu =  Ry*Rz*Rx;
////          Mat3 R_imu =  Ry*Rx*Rz;
//          Mat3 R_imu =  Rz*Rx*Ry;//
////          Mat3 R_imu =  Rz*Ry*Rx;


          //use now
          Mat3 Rz, Ry, Rx;
          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
                0.0, 0.0, 1.0;
          Ry << cos(phi/180.0*M_PI), 0.0, sin(phi/180.0*M_PI),
                0.0, 1.0, 0.0,
                -sin(phi/180.0*M_PI), 0.0, cos(phi/180.0*M_PI);
          Rx << 1.0, 0.0, 0.0,
                0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
                0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
//          Mat3 R_imu =  Rz*Ry*Rx;
          Mat3 R_imu =  Rx*Ry*Rz;



//          Mat3 Rz, Ry, Rx;
//          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
//                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ry << 1.0, 0.0, 0.0,
//                0.0, -cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//                0.0, sin(phi/180.0*M_PI), -cos(phi/180.0*M_PI);
//          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;

//          Rz << cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI), 0.0,
//                sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ry << 1.0, 0.0, 0.0,
//                0.0, -cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
//                0.0, sin(omega/180.0*M_PI), -cos(omega/180.0*M_PI);
//          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;

//          Mat3 Rz, Ry, Rx;
//          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
//                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ry << 1.0, 0.0, 0.0,
//                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);
//          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;          Mat3 Rz, Ry, Rx;
//          Mat3 Rz, Ry, Rx;
//          Rz << cos(phi/180.0*M_PI), sin(phi/180.0*M_PI), 0.0,
//                -sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ry << 1.0, 0.0, 0.0,
//                0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
//                0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
//          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;



//          Mat3 R_imu =  Ry*Rz*Rx;
//          Mat3 R_imu =  Rz*Ry*Rx;
//          Mat3 R_imu =  Rx*Ry*Rz;


//          Mat3 Rk, Ra, RA;
//          Rk << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          RA << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
//                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ra << 1.0, 0.0, 0.0,
//                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


//          Mat3 R_imu =  Ry*Rx*Rz;

//          Mat3 changeAxis_M;
//          changeAxis_M << 1.0, 0.0, 0.0,
//                          0.0, -1.0, 0.0,
//                          0.0, 0.0, -1.0;  //M1
//          changeAxis_M << 0.0, 1.0, 0.0,
//                          1.0, 0.0, 0.0,
//                          0.0, 0.0, -1.0;  //M2

//          R_imu = changeAxis_M * R_imu;
//          R_imu = R_imu * changeAxis_M;





//              tempGet = getListXYZ(-phi3, kappa3, -omega3);
//          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

//               tempGet = getListXYZ(phi3, kappa3, omega3);
//          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

//          tempGet = getListXYZ(phi1, kappa1, omega1);
//          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

//          Mat3 Rz, Ry, Rx;
//          Rx << 1.0, 0.0, 0.0,
//                 0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//                 0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

//          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
//                 0.0, 1.0, 0.0,
//                 -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);

//          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                 sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                 0.0, 0.0, 1.0;


//          Mat3 R_imu =  Ry*Rz*Rx;
//          Mat3 R_imu =  Rz*Rx*Ry;
//          Mat3 R_imu =  Rx*Rz*Ry;
//          Mat3 R_imu =  Rx*Ry*Rz;

          Mat3 changeAxis_M;
          changeAxis_M << 1.0, 0.0, 0.0,
                          0.0, -1.0, 0.0,
                          0.0, 0.0, -1.0;
//          changeAxis_M << 1.0, 0.0, 0.0,
//                          0.0, 1.0, 0.0,
//                          0.0, 0.0, -1.0;

//          changeAxis_M << 1.0, 0.0, 0.0,
//                          0.0, 1.0, 0.0,
//                          0.0, 0.0, -1.0;
          Mat3 R_;
          R_ << 1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, -1.0;
          R_imu = changeAxis_M * R_imu;//* R_;


          //Mat3 R2;
          //ceres::AngleAxisToRotationMatrix(R2_a.data(), R2.data());

          //R_imu = changeAxis_M * R_imu * R2;//* R_;


          std::vector<double> R_imu_a(3);
          ceres::RotationMatrixToAngleAxis(R_imu.data(), R_imu_a.data());

          Vec3 C_gps_res = Vec3(c[0], c[1], c[2]);


          {
              bool find = false;
              for(Views::const_iterator itV = sfm_data_.views.begin(); itV != sfm_data_.views.end(); ++itV)
              {
                  if(itV->first == imgId)
                  {
                      find = true;
                      break;
                  }
              }
              if(find == false)
                  continue;
          }
          // check status
          if(imgId/100000 == 1)//back
          {
              if(imgStatus != 321)
              {
                  C_gps_Map[imgId + 100000] = c;
                  R_imu_Map[imgId + 100000] = R_imu_a;
              }

              //set back-camera pose
              Mat3 Rb_res = Rb_trans * R_imu_trans * R_imu;

              Vec3 tb_res;
              tb_res = Rb_trans * gps_trans - Rb_trans * R_imu_trans * R_imu * C_gps_res + tb_trans;
              global_rotations[imgId] = Rb_res;

              sfm_data_.poses[imgId] = Pose3(Rb_res, -Rb_res.transpose() * tb_res);

              PoseStatus poseStatus;
              poseStatus.pose = sfm_data_.poses[imgId];
              poseStatus.status = imgStatus;
              poseStatus.bdf = 1;
              poseStatusMap[imgId] = poseStatus;

          }else if(imgId/100000 == 2)//down
          {
              C_gps_Map[imgId] = c;
              R_imu_Map[imgId] = R_imu_a;

              //set down-camera pose
              Mat3 Rd_res = Rd_trans * R_imu_trans * R_imu;

              Vec3 td_res;
              td_res = Rd_trans * gps_trans - Rd_trans * R_imu_trans * R_imu * C_gps_res + td_trans;
              global_rotations[imgId] = Rd_res;

              sfm_data_.poses[imgId] = Pose3(Rd_res, -Rd_res.transpose() * td_res);

              PoseStatus poseStatus;
              poseStatus.pose = sfm_data_.poses[imgId];
              poseStatus.status = imgStatus;
              poseStatus.bdf = 2;
              poseStatusMap[imgId] = poseStatus;

          }else if(imgId/100000 == 3)//front
          {
              if(imgStatus != 321)
              {
                  C_gps_Map[imgId - 100000] = c;
                  R_imu_Map[imgId - 100000] = R_imu_a;
              }

              //set front-camera pose
              Mat3 Rf_res = Rf_trans * R_imu_trans * R_imu;

              Vec3 tf_res;
              tf_res = Rf_trans * gps_trans - Rf_trans * R_imu_trans * R_imu * C_gps_res + tf_trans;
              global_rotations[imgId] = Rf_res;

              sfm_data_.poses[imgId] = Pose3(Rf_res, -Rf_res.transpose() * tf_res);

              PoseStatus poseStatus;
              poseStatus.pose = sfm_data_.poses[imgId];
              poseStatus.status = imgStatus;
              poseStatus.bdf = 3;
              poseStatusMap[imgId] = poseStatus;

          }

      }
      in.close();

  }

  //Refine_initValue();
  matching::PairWiseMatches tripletWise_matches;
  if (!Compute_Global_Translations_gps(global_rotations, tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
    return false;
  }


  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }



    {
//            std::cout<< "initial intrinsics :"<<std::endl;
//                Hash_Map<IndexT, std::vector<double> > map_intrinsics;
//                map_intrinsics[0] = sfm_data_.intrinsics[0]->getParams();
//                std::cout<< "map_intrinsics 0 " << map_intrinsics[0][0]<<" "<<map_intrinsics[0][1]<<" "<<
//                         map_intrinsics[0][2]<<" "<<map_intrinsics[0][3]<<" "<<map_intrinsics[0][4]<<" "<<
//                         map_intrinsics[0][5]<<" "<<map_intrinsics[0][6]<<" "<<map_intrinsics[0][7]<<std::endl;

//                map_intrinsics[1] = sfm_data_.intrinsics[1]->getParams();
//                std::cout<< "map_intrinsics 1 " << map_intrinsics[1][0]<<" "<<map_intrinsics[1][1]<<" "<<
//                         map_intrinsics[1][2]<<" "<<map_intrinsics[1][3]<<" "<<map_intrinsics[1][4]<<" "<<
//                         map_intrinsics[1][5]<<" "<<map_intrinsics[1][6]<<" "<<map_intrinsics[1][7]<<std::endl;

//                map_intrinsics[2] = sfm_data_.intrinsics[2]->getParams();
//                std::cout<< "map_intrinsics 2 " << map_intrinsics[2][0]<<" "<<map_intrinsics[2][1]<<" "<<
//                         map_intrinsics[2][2]<<" "<<map_intrinsics[2][3]<<" "<<map_intrinsics[2][4]<<" "<<
//                         map_intrinsics[2][5]<<" "<<map_intrinsics[2][6]<<" "<<map_intrinsics[2][7]<<std::endl;
//                std::cout <<"-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() <<std::endl;
    }


  Save(sfm_data_,
    stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "cloud_and_poses_first_init", ".ply"),
    ESfM_Data(ALL));

  if (!Adjust_init())//_fix2())//(!Adjust_threecamera_gps_cail_new(poseStatusMap, transformation_R_imu, tramsformation_x_gps, transformation_br, transformation_fr, C_gps_Map, R_imu_Map))
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }



  {
//      double focal = 7265.87;
//      Vec3 refP = Vec3(18449796.870187, 3496035.996928, 504.744);
//      refP(0) -= c0[0];
//      refP(1) -= c0[1];
//      refP(2) -= c0[2];

//      const Pose3 & pose1 = sfm_data_.poses[210092];
//      const Pose3 & pose2 = sfm_data_.poses[210093];
//      const Pose3 & pose3 = sfm_data_.poses[210094];

//      const Vec3 t1 = pose1.translation();
//      const Mat3 R1 = pose1.rotation();
//      Vec3 x1 = R1 * refP + t1;
//      Vec3 point1 = Vec3(x1(0)/x1(2)*focal+3680, x1(1)/x1(2)*focal+2456, 1.0);

//      const Vec3 t2 = pose2.translation();
//      const Mat3 R2 = pose2.rotation();
//      Vec3 x2 = R2 * refP + t2;
//      Vec3 point2 = Vec3(x2(0)/x2(2)*focal+3680, x2(1)/x2(2)*focal+2456, 1.0);

//      const Vec3 t3 = pose3.translation();
//      const Mat3 R3 = pose3.rotation();
//      Vec3 x3 = R3 * refP + t3;
//      Vec3 point3 = Vec3(x3(0)/x3(2)*focal+3680, x3(1)/x3(2)*focal+2456, 1.0);



//      std::cout <<"point1 " << point1(0)<< " "<< point1(1)<< " "<< point1(2) << std::endl;
//      std::cout <<"point2 " << point2(0)<< " "<< point2(1)<< " "<< point2(2) << std::endl;
//      std::cout <<"point3 " << point3(0)<< " " << point3(1)<< " " << point3(2) << std::endl<< std::endl<< std::endl;
  }


  {
//      std::cout << "...Export SfM_Data to disk." << std::endl;
//      Save(sfm_data_,
//        stlplus::create_filespec("/home/guang/data/test_Blocks/block_4/matches", "sfm_data6", ".bin"),
//        ESfM_Data(ALL));
  }


  {
//          std::cout<< "map_intrinsics 1"<<std::endl;
//              Hash_Map<IndexT, std::vector<double> > map_intrinsics;
//              map_intrinsics[0] = sfm_data_.intrinsics[0]->getParams();
//              std::cout<< "map_intrinsics 0 " << map_intrinsics[0][0]<<" "<<map_intrinsics[0][1]<<" "<<
//                       map_intrinsics[0][2]<<" "<<map_intrinsics[0][3]<<" "<<map_intrinsics[0][4]<<" "<<
//                       map_intrinsics[0][5]<<" "<<map_intrinsics[0][6]<<" "<<map_intrinsics[0][7]<<std::endl;

//              map_intrinsics[1] = sfm_data_.intrinsics[1]->getParams();
//              std::cout<< "map_intrinsics 1 " << map_intrinsics[1][0]<<" "<<map_intrinsics[1][1]<<" "<<
//                       map_intrinsics[1][2]<<" "<<map_intrinsics[1][3]<<" "<<map_intrinsics[1][4]<<" "<<
//                       map_intrinsics[1][5]<<" "<<map_intrinsics[1][6]<<" "<<map_intrinsics[1][7]<<std::endl;

//              map_intrinsics[2] = sfm_data_.intrinsics[2]->getParams();
//              std::cout<< "map_intrinsics 2 " << map_intrinsics[2][0]<<" "<<map_intrinsics[2][1]<<" "<<
//                       map_intrinsics[2][2]<<" "<<map_intrinsics[2][3]<<" "<<map_intrinsics[2][4]<<" "<<
//                       map_intrinsics[2][5]<<" "<<map_intrinsics[2][6]<<" "<<map_intrinsics[2][7]<<std::endl;
//              std::cout <<"-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() <<std::endl;
  }


  //out
  {
//      std::ofstream out_show_camera;
//      out_show_camera.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_reference.txt");
//      std::ofstream out_Rc;
//      out_Rc.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_Rc.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 c = pose.center();

//          Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
//          axis_x = R.transpose() * (axis_x + R*c);
//          axis_y = R.transpose() * (axis_y + R*c);
//          axis_z = R.transpose() * (axis_z + R*c);
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_x(0) << " " << axis_x(1) << " " << axis_x(2)<< ";" << std::endl;
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_y(0) << " " << axis_y(1) << " " << axis_y(2)<< ";" << std::endl;
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_z(0) << " " << axis_z(1) << " " << axis_z(2)<< ";" << std::endl;

//          out_Rc << c(0) << " " << c(1)<< " " << c(2) << " " << std::endl;

//      }
//      out_Rc.close();
//      out_show_camera.close();

  }

  //out
  {

//      std::ofstream out_c;
//      out_c.open("/home/guang/data/beichuan_testBA4/sfm3/out_show_camera_c.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Vec3 c = pose.center();
//          const Mat3 R = pose.rotation();

//          out_c <<it->first << " "<< c(0) << " " << c(1)<< " " << c(2) << " "
//                <<R(0,0)<<" "<< R(0,1)<< " "<<R(0,2)<< " "
//                <<R(1,0)<<" "<< R(1,1)<< " "<<R(1,2)<< " "
//                <<R(2,0)<<" "<< R(2,1)<< " "<<R(2,2)<< std::endl;

//      }
//      out_c.close();

  }

  //move
  {
//      //std::ofstream out_c;
//      //out_c.open("/home/guang/data/beichuan_testBA4/sfm3/out_show_camera_c.txt");
//      for(Poses::iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          Pose3 & pose = it->second;
//          const Vec3 c = pose.center();
//          const Mat3 R = pose.rotation();

//          Vec3 c_new = Vec3(c(0)+c0[0], c(1)+c0[1], c(2)+c0[2]);

//          pose = Pose3(R, c_new);

//      }

//      for (Landmarks::iterator structure_landmark_it = sfm_data_.structure.begin();
//           structure_landmark_it != sfm_data_.structure.end(); ++structure_landmark_it)
//      {
//          Vec3 X_old = structure_landmark_it->second.X;
//          Vec3 X_new = Vec3(X_old(0)+c0[0], X_old(1)+c0[1], X_old(2)+c0[2]);
//          structure_landmark_it->second.X = X_new;

//      }

  }

  {
//      std::ofstream out_show_camera;
//      out_show_camera.open("/media/add7/E/run_scene/testMerge/29/out_camera_res_R.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 c = pose.center();

//          Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
//          axis_x = R.transpose() * (axis_x + R*c);
//          axis_y = R.transpose() * (axis_y + R*c);
//          axis_z = R.transpose() * (axis_z + R*c);
//          out_show_camera <<R(0,0)<<" "<< R(0,1)<< " "<<R(0,2)<< " "
//                          <<R(1,0)<<" "<< R(1,1)<< " "<<R(1,2)<< " "
//                          <<R(2,0)<<" "<< R(2,1)<< " "<<R(2,2)<< std::endl;
//      }
//      out_show_camera.close();
  }





  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }





  return true;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process_threeCameras_change_step(std::string rtFilePath, std::string refSfmDataIdList) {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if(set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

  Hash_Map<IndexT, Mat3> global_rotations;

  //init inter-constraint
  Mat3 Rb_trans, Rf_trans, Rd_trans;
  Vec3 tb_trans, tf_trans, td_trans;
  Mat3 R_imu_trans;
  Vec3 gps_trans;
  std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
  const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
  const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};
  {
      //test distortion
      Mat3 r_tb, r_tf;
      r_tb << 1.0, 0.0, 0.0,
              0.0, 0.787926, -0.61577,
              0.0, 0.61577, 0.787926;
      r_tf << 1.0, 0.0, 0.0,
              0.0, 0.78796,  0.615727,
              0.0, -0.615727, 0.78796;

//      Rx << 1.0, 0.0, 0.0,
//             0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//             0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


//      r_tb << 1.0, 0.0, 0.0,
//              0.0, cos(38.0/180.0*M_PI), -sin(38.0/180.0*M_PI),
//              0.0, sin(38.0/180.0*M_PI), cos(38.0/180.0*M_PI);

//      r_tf << 1.0, 0.0, 0.0,
//              0.0, cos(-38.0/180.0*M_PI), -sin(-38.0/180.0*M_PI),
//              0.0, sin(-38.0/180.0*M_PI), cos(-38.0/180.0*M_PI);

//      r_tb << cos(-38.0/180.0*M_PI), 0.0, sin(-38.0/180.0*M_PI),//-38
//              0.0, 1.0, 0.0,
//              -sin(-38.0/180.0*M_PI), 0.0, cos(-38.0/180.0*M_PI);
//      r_tf << cos(38.0/180.0*M_PI), 0.0, sin(38.0/180.0*M_PI),//38
//              0.0, 1.0, 0.0,
//              -sin(38.0/180.0*M_PI), 0.0, cos(38.0/180.0*M_PI);
      Vec3 RA_tb, RA_tf;
      ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
      ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

      //set reference
      double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), -RA_tb(2)};
      double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), -RA_tf(2)};
//      double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), RA_tb(2)};
//      double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), RA_tf(2)};
      double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};

      Vec3 c_bd, c_fd;//b-d, f-d
      c_bd << 0.0, 0.098, 0.027;
      c_fd << 0.0, -0.101, 0.028;
//      c_bd << -0.098, 0.0, -0.027;
//      c_fd << 0.101, 0.0, -0.028;

//      c_bd << 0.098, 0.0, -0.027;
//      c_fd << -0.101, 0.0, -0.028;
//      c_bd << 0.0, 1.0, 0.027;
//      c_fd << 0.0, -0.101, 0.028;

      Vec3 t_db, t_df, t_dd;
      t_db = -r_tb*c_bd;
      t_df = -r_tf*c_fd;

      //double S = 1470.531783220221;
    //        t_db << 0.0, 0.098, 0.027;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
    //        t_df << 0.0, -0.101, 0.028;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
      t_dd << 0.0, 0.0, 0.0;

      transformation_br = {R_br_angleAxis[0], R_br_angleAxis[1], R_br_angleAxis[2],
                           t_db(0), t_db(1), t_db(2)};
      transformation_fr = {R_fr_angleAxis[0], R_fr_angleAxis[1], R_fr_angleAxis[2],
                           t_df(0), t_df(1), t_df(2)};
      transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                           t_dd(0), t_dd(1), t_dd(2)};

//      const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
//      const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};

      std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                             << tramsformation_x_gps[1] << " "
                                             << tramsformation_x_gps[2] << std::endl;

      //prepare for setting pose
      ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rb_trans.data());
      tb_trans = Vec3(transformation_br[3], transformation_br[4], transformation_br[5]);
      ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rf_trans.data());
      tf_trans = Vec3(transformation_fr[3], transformation_fr[4], transformation_fr[5]);
      ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rd_trans.data());
      td_trans = Vec3(transformation_dr[3], transformation_dr[4], transformation_dr[5]);


      ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
      gps_trans = Vec3(tramsformation_x_gps[0], tramsformation_x_gps[1], tramsformation_x_gps[2]);

  }


  //get gps and imu data
  PoseStatusMap poseStatusMap;
  Hash_Map<IndexT, std::vector<double> > C_gps_Map;
  Hash_Map<IndexT, std::vector<double> > R_imu_Map;
  std::vector<double> c0(3);
  {
    //        /Hash_Map<IndexT, std::vector<double> > R_imu;
//            Hash_Map<IndexT, std::vector<double> > C_gps;
      std::string tFilePath = rtFilePath;  //!!!
      std::ifstream in;
      in.open(tFilePath);
      if(!in)
      {
          std::cout << "open tFilePath file error!" << std::endl;
          return false;
      }
      std::size_t line = 0;
      int count = -1;
      in >> count;  ///!!!

      bool isC0 = false;//true;

      while((!in.eof()) && (line < count))
      {
          // input
          int imgId, imgStatus;
          std::vector<double> c(3);
          //double pitch, roll, yaw;
          double omega, phi, kappa;
          in >> imgId
             >> c[0] >> c[1] >> c[2]
//             >> c[1] >> c[0] >> c[2]
             >> omega >> phi >> kappa
             >> imgStatus;

//          omega = omega + 2.4;
//          phi = phi + 0.0713;
//          kappa = kappa - 0.5805;
//          kappa = kappa + 2.4;
//          phi = phi + 0.0713;
//          omega = omega - 0.5805;

          if(isC0 == false)
          {
              c0 = c;
              isC0 = true;
          }

//          c[0] -= c0[0];
//          c[1] -= c0[1];
//          c[2] -= c0[2];

          ++line;



//          double roll, pitch, yaw;

//          roll = omega;
//          pitch = phi;
//          yaw = kappa;

//          Mat3 Rz, Ry, Rx;
//          Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
//                sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
//                0.0, 1.0, 0.0,
//                -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
//          Rx << 1.0, 0.0, 0.0,
//                0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
//                0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);
//          Mat3 R_imu =  Rz*Rx*Ry;
//          Mat3 changeAxis_M;
//          changeAxis_M << 1.0, 0.0, 0.0,
//                          0.0, -1.0, 0.0,
//                          0.0, 0.0, -1.0;
//          R_imu = changeAxis_M * R_imu;




//          Mat3 Rz, Ry, Rx;
//          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
//                0.0, 1.0, 0.0,
//                -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);
//          Rx << 1.0, 0.0, 0.0,
//                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

////          Mat3 R_imu =  Rx*Ry*Rz;
////          Mat3 R_imu =  Rx*Rz*Ry;//
////          Mat3 R_imu =  Ry*Rz*Rx;
////          Mat3 R_imu =  Ry*Rx*Rz;
//          Mat3 R_imu =  Rz*Rx*Ry;//
////          Mat3 R_imu =  Rz*Ry*Rx;


          //use now
          Mat3 Rz, Ry, Rx;
          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
                0.0, 0.0, 1.0;
          Ry << cos(phi/180.0*M_PI), 0.0, sin(phi/180.0*M_PI),
                0.0, 1.0, 0.0,
                -sin(phi/180.0*M_PI), 0.0, cos(phi/180.0*M_PI);
          Rx << 1.0, 0.0, 0.0,
                0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
                0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
//          Mat3 R_imu =  Rz*Ry*Rx;
          Mat3 R_imu =  Rx*Ry*Rz;



//          Mat3 Rz, Ry, Rx;
//          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
//                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ry << 1.0, 0.0, 0.0,
//                0.0, -cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//                0.0, sin(phi/180.0*M_PI), -cos(phi/180.0*M_PI);
//          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;

//          Rz << cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI), 0.0,
//                sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ry << 1.0, 0.0, 0.0,
//                0.0, -cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
//                0.0, sin(omega/180.0*M_PI), -cos(omega/180.0*M_PI);
//          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;

//          Mat3 Rz, Ry, Rx;
//          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
//                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ry << 1.0, 0.0, 0.0,
//                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);
//          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;          Mat3 Rz, Ry, Rx;
//          Mat3 Rz, Ry, Rx;
//          Rz << cos(phi/180.0*M_PI), sin(phi/180.0*M_PI), 0.0,
//                -sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ry << 1.0, 0.0, 0.0,
//                0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
//                0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
//          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;



//          Mat3 R_imu =  Ry*Rz*Rx;
//          Mat3 R_imu =  Rz*Ry*Rx;
//          Mat3 R_imu =  Rx*Ry*Rz;


//          Mat3 Rk, Ra, RA;
//          Rk << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          RA << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
//                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ra << 1.0, 0.0, 0.0,
//                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


//          Mat3 R_imu =  Ry*Rx*Rz;

//          Mat3 changeAxis_M;
//          changeAxis_M << 1.0, 0.0, 0.0,
//                          0.0, -1.0, 0.0,
//                          0.0, 0.0, -1.0;  //M1
//          changeAxis_M << 0.0, 1.0, 0.0,
//                          1.0, 0.0, 0.0,
//                          0.0, 0.0, -1.0;  //M2

//          R_imu = changeAxis_M * R_imu;
//          R_imu = R_imu * changeAxis_M;





//              tempGet = getListXYZ(-phi3, kappa3, -omega3);
//          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

//               tempGet = getListXYZ(phi3, kappa3, omega3);
//          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

//          tempGet = getListXYZ(phi1, kappa1, omega1);
//          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

//          Mat3 Rz, Ry, Rx;
//          Rx << 1.0, 0.0, 0.0,
//                 0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//                 0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

//          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
//                 0.0, 1.0, 0.0,
//                 -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);

//          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                 sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                 0.0, 0.0, 1.0;


//          Mat3 R_imu =  Ry*Rz*Rx;
//          Mat3 R_imu =  Rz*Rx*Ry;
//          Mat3 R_imu =  Rx*Rz*Ry;
//          Mat3 R_imu =  Rx*Ry*Rz;

          Mat3 changeAxis_M;
          changeAxis_M << 1.0, 0.0, 0.0,
                          0.0, -1.0, 0.0,
                          0.0, 0.0, -1.0;
//          changeAxis_M << 1.0, 0.0, 0.0,
//                          0.0, 1.0, 0.0,
//                          0.0, 0.0, -1.0;

//          changeAxis_M << 1.0, 0.0, 0.0,
//                          0.0, 1.0, 0.0,
//                          0.0, 0.0, -1.0;
          Mat3 R_;
          R_ << 1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, -1.0;
          R_imu = changeAxis_M * R_imu;//* R_;


          //Mat3 R2;
          //ceres::AngleAxisToRotationMatrix(R2_a.data(), R2.data());

          //R_imu = changeAxis_M * R_imu * R2;//* R_;


          std::vector<double> R_imu_a(3);
          ceres::RotationMatrixToAngleAxis(R_imu.data(), R_imu_a.data());

          Vec3 C_gps_res = Vec3(c[0], c[1], c[2]);


          {
              bool find = false;
              for(Views::const_iterator itV = sfm_data_.views.begin(); itV != sfm_data_.views.end(); ++itV)
              {
                  if(itV->first == imgId)
                  {
                      find = true;
                      break;
                  }
              }
              if(find == false)
                  continue;
          }
          // check status
          if(imgId/100000 == 1)//back
          {
              if(imgStatus != 321)
              {
                  C_gps_Map[imgId + 100000] = c;
                  R_imu_Map[imgId + 100000] = R_imu_a;
              }

              //set back-camera pose
              Mat3 Rb_res = Rb_trans * R_imu_trans * R_imu;

              Vec3 tb_res;
              tb_res = Rb_trans * gps_trans - Rb_trans * R_imu_trans * R_imu * C_gps_res + tb_trans;
              global_rotations[imgId] = Rb_res;

              sfm_data_.poses[imgId] = Pose3(Rb_res, -Rb_res.transpose() * tb_res);

              PoseStatus poseStatus;
              poseStatus.pose = sfm_data_.poses[imgId];
              poseStatus.status = imgStatus;
              poseStatus.bdf = 1;
              poseStatusMap[imgId] = poseStatus;

          }else if(imgId/100000 == 2)//down
          {
              C_gps_Map[imgId] = c;
              R_imu_Map[imgId] = R_imu_a;

              //set down-camera pose
              Mat3 Rd_res = Rd_trans * R_imu_trans * R_imu;

              Vec3 td_res;
              td_res = Rd_trans * gps_trans - Rd_trans * R_imu_trans * R_imu * C_gps_res + td_trans;
              global_rotations[imgId] = Rd_res;

              sfm_data_.poses[imgId] = Pose3(Rd_res, -Rd_res.transpose() * td_res);

              PoseStatus poseStatus;
              poseStatus.pose = sfm_data_.poses[imgId];
              poseStatus.status = imgStatus;
              poseStatus.bdf = 2;
              poseStatusMap[imgId] = poseStatus;

          }else if(imgId/100000 == 3)//front
          {
              if(imgStatus != 321)
              {
                  C_gps_Map[imgId - 100000] = c;
                  R_imu_Map[imgId - 100000] = R_imu_a;
              }

              //set front-camera pose
              Mat3 Rf_res = Rf_trans * R_imu_trans * R_imu;

              Vec3 tf_res;
              tf_res = Rf_trans * gps_trans - Rf_trans * R_imu_trans * R_imu * C_gps_res + tf_trans;
              global_rotations[imgId] = Rf_res;

              sfm_data_.poses[imgId] = Pose3(Rf_res, -Rf_res.transpose() * tf_res);

              PoseStatus poseStatus;
              poseStatus.pose = sfm_data_.poses[imgId];
              poseStatus.status = imgStatus;
              poseStatus.bdf = 3;
              poseStatusMap[imgId] = poseStatus;

          }

      }
      in.close();

  }



  ///record res
  ///
  std::map<int, Pose3> fixIdList;

  //prepare input
  std::vector<std::string> refIdList;
  boost::split(refIdList, refSfmDataIdList, boost::is_any_of("_"));
//  std::cout << refSfmDataIdList.size() << std::endl;
//  std::cout << "split id: " << std::endl;
//  std::cout << refIdList.size()<< std::endl;

  if(refSfmDataIdList.size() != 0)
  {
      for(std::size_t i = 0; i < refIdList.size(); ++i)
      {
          std::cout << refIdList[i] << std::endl;
          // Load input SfM_Data scene
          SfM_Data ref_sfm_data;
          if (!Load(ref_sfm_data, sfm_data_.s_root_path+"/sfm_data_"+refIdList[i]+"_extrinsics.json", ESfM_Data(VIEWS|INTRINSICS|EXTRINSICS))) {
            //std::cerr << std::endl
              //poseRef<< "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
            return EXIT_FAILURE;
          }

          for(Poses::iterator itPthis = sfm_data_.poses.begin(); itPthis != sfm_data_.poses.end(); ++itPthis)
          {
              for(Poses::const_iterator itPref = ref_sfm_data.poses.begin(); itPref != ref_sfm_data.poses.end(); ++itPref)
              {
                  if(itPthis->first == itPref->first)
                  {
                      const Pose3 & poseRef = itPref->second;

                      std::pair<int, Pose3> tempPair;
                      tempPair.first = itPthis->first;
                      tempPair.second = poseRef;
                      fixIdList.insert(tempPair);
                  }
              }
          }

          //change
          {
              int length_N = 0;

              double length_gps = 0;
              double length_c = 0;

              double changeS = 0;
              for(std::map<int, Pose3>::const_iterator itP = fixIdList.begin(); itP != fixIdList.end(); ++itP)
              {
                  const Vec3 c1 = itP->second.center();
                  int pId1 = itP->first;
                  const Vec3 cg1 = fixIdList.find(pId1)->second.center();

                  std::map<int, Pose3>::const_iterator itPose2 = itP;
                  ++ itPose2;
                  for(;itPose2 != fixIdList.end(); ++itPose2)
                  {
                      const Vec3 c2 = itPose2->second.center();
                      int pId2 = itPose2->first;
                      const Vec3 cg2 = fixIdList.find(pId2)->second.center();

                      double lg = (cg1-cg2).norm();
                      double lc = (c1-c2).norm();

                      length_gps += lg;
                      length_c += lc;

                      ++length_N;

                  }

              }
              length_gps /= (double)(length_N);
              length_c /= (double)(length_N);

              changeS = length_gps / length_c;


              //compute H -> Xg = H * Xc
              std::vector<Vec3> Xc, Xg;
              for(std::map<int, Pose3>::const_iterator itPose = fixIdList.begin();
                  itPose != fixIdList.end(); ++itPose)
              {
                  int pId = itPose->first;
                  //prepare Xg
                  Vec3 Xg_temp;
                  Xg_temp = fixIdList.find(pId)->second.center();
                  Xg.push_back(Xg_temp);

                  //prepare Xc
                  Vec3 Xc_temp;
                  Xc_temp = changeS * sfm_data_.poses[pId].center();
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

              //get and set all pose
              std::cout << "s: " << changeS << std::endl;
              std::cout << "R: " << R_gc << std::endl;
              std::cout << "t: " << t_gc << std::endl;
              for (Poses::iterator itPose = sfm_data_.poses.begin(); itPose != sfm_data_.poses.end(); ++itPose)
              {
                  Pose3 & pose = itPose->second;
                  const Mat3 R = pose.rotation();
                  const Vec3 t = pose.translation();

                  Mat3 newR = (1/changeS)* R*R_gc.inverse();
                  Vec3 newt = t - newR * t_gc;

                  global_rotations[itPose->first] = newR;
                  pose = Pose3(newR, -newR.inverse() * newt);
              }
          }



      }

  }



  matching::PairWiseMatches tripletWise_matches;
  if (!Compute_Global_Translations_gps(global_rotations, tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
    return false;
  }


  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }

  {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "cloud_and_poses_first_init", ".ply"),
        ESfM_Data(ALL));
  }

  if(refSfmDataIdList.size() == 0)
  {
      if (!Adjust_init())//_fix2())//(!Adjust_threecamera_gps_cail_new(poseStatusMap, transformation_R_imu, tramsformation_x_gps, transformation_br, transformation_fr, C_gps_Map, R_imu_Map))
      {
        std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
        return false;
      }

  }else{
      if (!Adjust_change_step(fixIdList))//_fix2())//(!Adjust_threecamera_gps_cail_new(poseStatusMap, transformation_R_imu, tramsformation_x_gps, transformation_br, transformation_fr, C_gps_Map, R_imu_Map))
      {
        std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
        return false;
      }
  }



  {
//      double focal = 7265.87;
//      Vec3 refP = Vec3(18449796.870187, 3496035.996928, 504.744);
//      refP(0) -= c0[0];
//      refP(1) -= c0[1];
//      refP(2) -= c0[2];

//      const Pose3 & pose1 = sfm_data_.poses[210092];
//      const Pose3 & pose2 = sfm_data_.poses[210093];
//      const Pose3 & pose3 = sfm_data_.poses[210094];

//      const Vec3 t1 = pose1.translation();
//      const Mat3 R1 = pose1.rotation();
//      Vec3 x1 = R1 * refP + t1;
//      Vec3 point1 = Vec3(x1(0)/x1(2)*focal+3680, x1(1)/x1(2)*focal+2456, 1.0);

//      const Vec3 t2 = pose2.translation();
//      const Mat3 R2 = pose2.rotation();
//      Vec3 x2 = R2 * refP + t2;
//      Vec3 point2 = Vec3(x2(0)/x2(2)*focal+3680, x2(1)/x2(2)*focal+2456, 1.0);

//      const Vec3 t3 = pose3.translation();
//      const Mat3 R3 = pose3.rotation();
//      Vec3 x3 = R3 * refP + t3;
//      Vec3 point3 = Vec3(x3(0)/x3(2)*focal+3680, x3(1)/x3(2)*focal+2456, 1.0);



//      std::cout <<"point1 " << point1(0)<< " "<< point1(1)<< " "<< point1(2) << std::endl;
//      std::cout <<"point2 " << point2(0)<< " "<< point2(1)<< " "<< point2(2) << std::endl;
//      std::cout <<"point3 " << point3(0)<< " " << point3(1)<< " " << point3(2) << std::endl<< std::endl<< std::endl;
  }


  {
//      std::cout << "...Export SfM_Data to disk." << std::endl;
//      Save(sfm_data_,
//        stlplus::create_filespec("/home/guang/data/test_Blocks/block_4/matches", "sfm_data6", ".bin"),
//        ESfM_Data(ALL));
  }


  {
//          std::cout<< "map_intrinsics 1"<<std::endl;
//              Hash_Map<IndexT, std::vector<double> > map_intrinsics;
//              map_intrinsics[0] = sfm_data_.intrinsics[0]->getParams();
//              std::cout<< "map_intrinsics 0 " << map_intrinsics[0][0]<<" "<<map_intrinsics[0][1]<<" "<<
//                       map_intrinsics[0][2]<<" "<<map_intrinsics[0][3]<<" "<<map_intrinsics[0][4]<<" "<<
//                       map_intrinsics[0][5]<<" "<<map_intrinsics[0][6]<<" "<<map_intrinsics[0][7]<<std::endl;

//              map_intrinsics[1] = sfm_data_.intrinsics[1]->getParams();
//              std::cout<< "map_intrinsics 1 " << map_intrinsics[1][0]<<" "<<map_intrinsics[1][1]<<" "<<
//                       map_intrinsics[1][2]<<" "<<map_intrinsics[1][3]<<" "<<map_intrinsics[1][4]<<" "<<
//                       map_intrinsics[1][5]<<" "<<map_intrinsics[1][6]<<" "<<map_intrinsics[1][7]<<std::endl;

//              map_intrinsics[2] = sfm_data_.intrinsics[2]->getParams();
//              std::cout<< "map_intrinsics 2 " << map_intrinsics[2][0]<<" "<<map_intrinsics[2][1]<<" "<<
//                       map_intrinsics[2][2]<<" "<<map_intrinsics[2][3]<<" "<<map_intrinsics[2][4]<<" "<<
//                       map_intrinsics[2][5]<<" "<<map_intrinsics[2][6]<<" "<<map_intrinsics[2][7]<<std::endl;
//              std::cout <<"-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() <<std::endl;
  }


  //out
  {
//      std::ofstream out_show_camera;
//      out_show_camera.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_reference.txt");
//      std::ofstream out_Rc;
//      out_Rc.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_Rc.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 c = pose.center();

//          Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
//          axis_x = R.transpose() * (axis_x + R*c);
//          axis_y = R.transpose() * (axis_y + R*c);
//          axis_z = R.transpose() * (axis_z + R*c);
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_x(0) << " " << axis_x(1) << " " << axis_x(2)<< ";" << std::endl;
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_y(0) << " " << axis_y(1) << " " << axis_y(2)<< ";" << std::endl;
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_z(0) << " " << axis_z(1) << " " << axis_z(2)<< ";" << std::endl;

//          out_Rc << c(0) << " " << c(1)<< " " << c(2) << " " << std::endl;

//      }
//      out_Rc.close();
//      out_show_camera.close();

  }

  //out
  {

//      std::ofstream out_c;
//      out_c.open("/home/guang/data/beichuan_testBA4/sfm3/out_show_camera_c.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Vec3 c = pose.center();
//          const Mat3 R = pose.rotation();

//          out_c <<it->first << " "<< c(0) << " " << c(1)<< " " << c(2) << " "
//                <<R(0,0)<<" "<< R(0,1)<< " "<<R(0,2)<< " "
//                <<R(1,0)<<" "<< R(1,1)<< " "<<R(1,2)<< " "
//                <<R(2,0)<<" "<< R(2,1)<< " "<<R(2,2)<< std::endl;

//      }
//      out_c.close();

  }

  //move
  {
//      //std::ofstream out_c;
//      //out_c.open("/home/guang/data/beichuan_testBA4/sfm3/out_show_camera_c.txt");
//      for(Poses::iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          Pose3 & pose = it->second;
//          const Vec3 c = pose.center();
//          const Mat3 R = pose.rotation();

//          Vec3 c_new = Vec3(c(0)+c0[0], c(1)+c0[1], c(2)+c0[2]);

//          pose = Pose3(R, c_new);

//      }

//      for (Landmarks::iterator structure_landmark_it = sfm_data_.structure.begin();
//           structure_landmark_it != sfm_data_.structure.end(); ++structure_landmark_it)
//      {
//          Vec3 X_old = structure_landmark_it->second.X;
//          Vec3 X_new = Vec3(X_old(0)+c0[0], X_old(1)+c0[1], X_old(2)+c0[2]);
//          structure_landmark_it->second.X = X_new;

//      }

  }

  {
//      std::ofstream out_show_camera;
//      out_show_camera.open("/media/add7/E/run_scene/testMerge/29/out_camera_res_R.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 c = pose.center();

//          Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
//          axis_x = R.transpose() * (axis_x + R*c);
//          axis_y = R.transpose() * (axis_y + R*c);
//          axis_z = R.transpose() * (axis_z + R*c);
//          out_show_camera <<R(0,0)<<" "<< R(0,1)<< " "<<R(0,2)<< " "
//                          <<R(1,0)<<" "<< R(1,1)<< " "<<R(1,2)<< " "
//                          <<R(2,0)<<" "<< R(2,1)<< " "<<R(2,2)<< std::endl;
//      }
//      out_show_camera.close();
  }





  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }





  return true;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Process_threeCameras_for_large(std::string rtFilePath) {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if(set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

  Hash_Map<IndexT, Mat3> global_rotations;

  //init inter-constraint
  Mat3 Rb_trans, Rf_trans, Rd_trans;
  Vec3 tb_trans, tf_trans, td_trans;
  Mat3 R_imu_trans;
  Vec3 gps_trans;
  std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
  const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
  const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};
  {
      //test distortion
      Mat3 r_tb, r_tf;
      r_tb << 1.0, 0.0, 0.0,
              0.0, 0.787926, -0.61577,
              0.0, 0.61577, 0.787926;
      r_tf << 1.0, 0.0, 0.0,
              0.0, 0.78796,  0.615727,
              0.0, -0.615727, 0.78796;


      Vec3 RA_tb, RA_tf;
      ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
      ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

      //set reference
      double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), -RA_tb(2)};
      double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), -RA_tf(2)};
//      double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), RA_tb(2)};
//      double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), RA_tf(2)};
      double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};

      Vec3 c_bd, c_fd;//b-d, f-d
      c_bd << 0.0, 0.098, 0.027;
      c_fd << 0.0, -0.101, 0.028;


      Vec3 t_db, t_df, t_dd;
      t_db = -r_tb*c_bd;
      t_df = -r_tf*c_fd;

      //double S = 1470.531783220221;
    //        t_db << 0.0, 0.098, 0.027;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
    //        t_df << 0.0, -0.101, 0.028;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
      t_dd << 0.0, 0.0, 0.0;

      transformation_br = {R_br_angleAxis[0], R_br_angleAxis[1], R_br_angleAxis[2],
                           t_db(0), t_db(1), t_db(2)};
      transformation_fr = {R_fr_angleAxis[0], R_fr_angleAxis[1], R_fr_angleAxis[2],
                           t_df(0), t_df(1), t_df(2)};
      transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                           t_dd(0), t_dd(1), t_dd(2)};

//      const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
//      const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};

      std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                             << tramsformation_x_gps[1] << " "
                                             << tramsformation_x_gps[2] << std::endl;

      //prepare for setting pose
      ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rb_trans.data());
      tb_trans = Vec3(transformation_br[3], transformation_br[4], transformation_br[5]);
      ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rf_trans.data());
      tf_trans = Vec3(transformation_fr[3], transformation_fr[4], transformation_fr[5]);
      ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rd_trans.data());
      td_trans = Vec3(transformation_dr[3], transformation_dr[4], transformation_dr[5]);


      ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
      gps_trans = Vec3(tramsformation_x_gps[0], tramsformation_x_gps[1], tramsformation_x_gps[2]);

  }


  //get gps and imu data
  PoseStatusMap poseStatusMap;
  Hash_Map<IndexT, std::vector<double> > C_gps_Map;
  Hash_Map<IndexT, std::vector<double> > R_imu_Map;
  std::vector<double> c0(3);
  {
    //        /Hash_Map<IndexT, std::vector<double> > R_imu;
//            Hash_Map<IndexT, std::vector<double> > C_gps;
      std::string tFilePath = rtFilePath;  //!!!
      std::ifstream in;
      in.open(tFilePath);
      if(!in)
      {
          std::cout << "open tFilePath file error!" << std::endl;
          return false;
      }
      std::size_t line = 0;
      int count = -1;
      in >> count;  ///!!!

      bool isC0 = false;//true;

      while((!in.eof()) && (line < count))
      {
          // input
          int imgId, imgStatus;
          std::vector<double> c(3);
          //double pitch, roll, yaw;
          double omega, phi, kappa;
          in >> imgId
             >> c[0] >> c[1] >> c[2]
//             >> c[1] >> c[0] >> c[2]
             >> omega >> phi >> kappa
             >> imgStatus;

//          omega = omega + 2.4;
//          phi = phi + 0.0713;
//          kappa = kappa - 0.5805;
//          kappa = kappa + 2.4;
//          phi = phi + 0.0713;
//          omega = omega - 0.5805;

          if(isC0 == false)
          {
              c0 = c;
              isC0 = true;
          }

//          c[0] -= c0[0];
//          c[1] -= c0[1];
//          c[2] -= c0[2];

          ++line;


          //use now
          Mat3 Rz, Ry, Rx;
          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
                0.0, 0.0, 1.0;
          Ry << cos(phi/180.0*M_PI), 0.0, sin(phi/180.0*M_PI),
                0.0, 1.0, 0.0,
                -sin(phi/180.0*M_PI), 0.0, cos(phi/180.0*M_PI);
          Rx << 1.0, 0.0, 0.0,
                0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
                0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
//          Mat3 R_imu =  Rz*Ry*Rx;
          Mat3 R_imu =  Rx*Ry*Rz;


          Mat3 changeAxis_M;
          changeAxis_M << 1.0, 0.0, 0.0,
                          0.0, -1.0, 0.0,
                          0.0, 0.0, -1.0;
//          changeAxis_M << 1.0, 0.0, 0.0,
//                          0.0, 1.0, 0.0,
//                          0.0, 0.0, -1.0;

//          changeAxis_M << 1.0, 0.0, 0.0,
//                          0.0, 1.0, 0.0,
//                          0.0, 0.0, -1.0;
          Mat3 R_;
          R_ << 1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, -1.0;
          R_imu = changeAxis_M * R_imu;//* R_;


          //Mat3 R2;
          //ceres::AngleAxisToRotationMatrix(R2_a.data(), R2.data());

          //R_imu = changeAxis_M * R_imu * R2;//* R_;


          std::vector<double> R_imu_a(3);
          ceres::RotationMatrixToAngleAxis(R_imu.data(), R_imu_a.data());

          Vec3 C_gps_res = Vec3(c[0], c[1], c[2]);

          // check status
          if(imgId/100000 == 1)//back
          {
              if(imgStatus != 321)
              {
                  C_gps_Map[imgId + 100000] = c;
                  R_imu_Map[imgId + 100000] = R_imu_a;
              }

              //set back-camera pose
              Mat3 Rb_res = Rb_trans * R_imu_trans * R_imu;

              Vec3 tb_res;
              tb_res = Rb_trans * gps_trans - Rb_trans * R_imu_trans * R_imu * C_gps_res + tb_trans;
              global_rotations[imgId] = Rb_res;

              sfm_data_.poses[imgId] = Pose3(Rb_res, -Rb_res.transpose() * tb_res);

              PoseStatus poseStatus;
              poseStatus.pose = sfm_data_.poses[imgId];
              poseStatus.status = imgStatus;
              poseStatus.bdf = 1;
              poseStatusMap[imgId] = poseStatus;

          }else if(imgId/100000 == 2)//down
          {
              C_gps_Map[imgId] = c;
              R_imu_Map[imgId] = R_imu_a;

              //set down-camera pose
              Mat3 Rd_res = Rd_trans * R_imu_trans * R_imu;

              Vec3 td_res;
              td_res = Rd_trans * gps_trans - Rd_trans * R_imu_trans * R_imu * C_gps_res + td_trans;
              global_rotations[imgId] = Rd_res;

              sfm_data_.poses[imgId] = Pose3(Rd_res, -Rd_res.transpose() * td_res);

              PoseStatus poseStatus;
              poseStatus.pose = sfm_data_.poses[imgId];
              poseStatus.status = imgStatus;
              poseStatus.bdf = 2;
              poseStatusMap[imgId] = poseStatus;

          }else if(imgId/100000 == 3)//front
          {
              if(imgStatus != 321)
              {
                  C_gps_Map[imgId - 100000] = c;
                  R_imu_Map[imgId - 100000] = R_imu_a;
              }

              //set front-camera pose
              Mat3 Rf_res = Rf_trans * R_imu_trans * R_imu;

              Vec3 tf_res;
              tf_res = Rf_trans * gps_trans - Rf_trans * R_imu_trans * R_imu * C_gps_res + tf_trans;
              global_rotations[imgId] = Rf_res;

              sfm_data_.poses[imgId] = Pose3(Rf_res, -Rf_res.transpose() * tf_res);

              PoseStatus poseStatus;
              poseStatus.pose = sfm_data_.poses[imgId];
              poseStatus.status = imgStatus;
              poseStatus.bdf = 3;
              poseStatusMap[imgId] = poseStatus;

          }

      }
      in.close();

  }

  //Refine_initValue();
  matching::PairWiseMatches tripletWise_matches;
  if (!Compute_Global_Translations_gps(global_rotations, tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
    return false;
  }


  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }



  Save(sfm_data_,
    stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "cloud_and_poses_first_init", ".ply"),
    ESfM_Data(ALL));

  if (!Adjust_for_large())//_fix2())//(!Adjust_threecamera_gps_cail_new(poseStatusMap, transformation_R_imu, tramsformation_x_gps, transformation_br, transformation_fr, C_gps_Map, R_imu_Map))
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }



  {
//      double focal = 7265.87;
//      Vec3 refP = Vec3(18449796.870187, 3496035.996928, 504.744);
//      refP(0) -= c0[0];
//      refP(1) -= c0[1];
//      refP(2) -= c0[2];

//      const Pose3 & pose1 = sfm_data_.poses[210092];
//      const Pose3 & pose2 = sfm_data_.poses[210093];
//      const Pose3 & pose3 = sfm_data_.poses[210094];

//      const Vec3 t1 = pose1.translation();
//      const Mat3 R1 = pose1.rotation();
//      Vec3 x1 = R1 * refP + t1;
//      Vec3 point1 = Vec3(x1(0)/x1(2)*focal+3680, x1(1)/x1(2)*focal+2456, 1.0);

//      const Vec3 t2 = pose2.translation();
//      const Mat3 R2 = pose2.rotation();
//      Vec3 x2 = R2 * refP + t2;
//      Vec3 point2 = Vec3(x2(0)/x2(2)*focal+3680, x2(1)/x2(2)*focal+2456, 1.0);

//      const Vec3 t3 = pose3.translation();
//      const Mat3 R3 = pose3.rotation();
//      Vec3 x3 = R3 * refP + t3;
//      Vec3 point3 = Vec3(x3(0)/x3(2)*focal+3680, x3(1)/x3(2)*focal+2456, 1.0);



//      std::cout <<"point1 " << point1(0)<< " "<< point1(1)<< " "<< point1(2) << std::endl;
//      std::cout <<"point2 " << point2(0)<< " "<< point2(1)<< " "<< point2(2) << std::endl;
//      std::cout <<"point3 " << point3(0)<< " " << point3(1)<< " " << point3(2) << std::endl<< std::endl<< std::endl;
  }


  {
//      std::cout << "...Export SfM_Data to disk." << std::endl;
//      Save(sfm_data_,
//        stlplus::create_filespec("/home/guang/data/test_Blocks/block_4/matches", "sfm_data6", ".bin"),
//        ESfM_Data(ALL));
  }


  {
//          std::cout<< "map_intrinsics 1"<<std::endl;
//              Hash_Map<IndexT, std::vector<double> > map_intrinsics;
//              map_intrinsics[0] = sfm_data_.intrinsics[0]->getParams();
//              std::cout<< "map_intrinsics 0 " << map_intrinsics[0][0]<<" "<<map_intrinsics[0][1]<<" "<<
//                       map_intrinsics[0][2]<<" "<<map_intrinsics[0][3]<<" "<<map_intrinsics[0][4]<<" "<<
//                       map_intrinsics[0][5]<<" "<<map_intrinsics[0][6]<<" "<<map_intrinsics[0][7]<<std::endl;

//              map_intrinsics[1] = sfm_data_.intrinsics[1]->getParams();
//              std::cout<< "map_intrinsics 1 " << map_intrinsics[1][0]<<" "<<map_intrinsics[1][1]<<" "<<
//                       map_intrinsics[1][2]<<" "<<map_intrinsics[1][3]<<" "<<map_intrinsics[1][4]<<" "<<
//                       map_intrinsics[1][5]<<" "<<map_intrinsics[1][6]<<" "<<map_intrinsics[1][7]<<std::endl;

//              map_intrinsics[2] = sfm_data_.intrinsics[2]->getParams();
//              std::cout<< "map_intrinsics 2 " << map_intrinsics[2][0]<<" "<<map_intrinsics[2][1]<<" "<<
//                       map_intrinsics[2][2]<<" "<<map_intrinsics[2][3]<<" "<<map_intrinsics[2][4]<<" "<<
//                       map_intrinsics[2][5]<<" "<<map_intrinsics[2][6]<<" "<<map_intrinsics[2][7]<<std::endl;
//              std::cout <<"-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() <<std::endl;
  }


  //out
  {
//      std::ofstream out_show_camera;
//      out_show_camera.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_reference.txt");
//      std::ofstream out_Rc;
//      out_Rc.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_Rc.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 c = pose.center();

//          Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
//          axis_x = R.transpose() * (axis_x + R*c);
//          axis_y = R.transpose() * (axis_y + R*c);
//          axis_z = R.transpose() * (axis_z + R*c);
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_x(0) << " " << axis_x(1) << " " << axis_x(2)<< ";" << std::endl;
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_y(0) << " " << axis_y(1) << " " << axis_y(2)<< ";" << std::endl;
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_z(0) << " " << axis_z(1) << " " << axis_z(2)<< ";" << std::endl;

//          out_Rc << c(0) << " " << c(1)<< " " << c(2) << " " << std::endl;

//      }
//      out_Rc.close();
//      out_show_camera.close();

  }

  //out
  {

//      std::ofstream out_c;
//      out_c.open("/home/guang/data/beichuan_testBA4/sfm3/out_show_camera_c.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Vec3 c = pose.center();
//          const Mat3 R = pose.rotation();

//          out_c <<it->first << " "<< c(0) << " " << c(1)<< " " << c(2) << " "
//                <<R(0,0)<<" "<< R(0,1)<< " "<<R(0,2)<< " "
//                <<R(1,0)<<" "<< R(1,1)<< " "<<R(1,2)<< " "
//                <<R(2,0)<<" "<< R(2,1)<< " "<<R(2,2)<< std::endl;

//      }
//      out_c.close();

  }

  //move
  {
//      //std::ofstream out_c;
//      //out_c.open("/home/guang/data/beichuan_testBA4/sfm3/out_show_camera_c.txt");
//      for(Poses::iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          Pose3 & pose = it->second;
//          const Vec3 c = pose.center();
//          const Mat3 R = pose.rotation();

//          Vec3 c_new = Vec3(c(0)+c0[0], c(1)+c0[1], c(2)+c0[2]);

//          pose = Pose3(R, c_new);

//      }

//      for (Landmarks::iterator structure_landmark_it = sfm_data_.structure.begin();
//           structure_landmark_it != sfm_data_.structure.end(); ++structure_landmark_it)
//      {
//          Vec3 X_old = structure_landmark_it->second.X;
//          Vec3 X_new = Vec3(X_old(0)+c0[0], X_old(1)+c0[1], X_old(2)+c0[2]);
//          structure_landmark_it->second.X = X_new;

//      }

  }

  {
//      std::ofstream out_show_camera;
//      out_show_camera.open("/media/add7/E/run_scene/testMerge/29/out_camera_res_R.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 c = pose.center();

//          Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
//          axis_x = R.transpose() * (axis_x + R*c);
//          axis_y = R.transpose() * (axis_y + R*c);
//          axis_z = R.transpose() * (axis_z + R*c);
//          out_show_camera <<R(0,0)<<" "<< R(0,1)<< " "<<R(0,2)<< " "
//                          <<R(1,0)<<" "<< R(1,1)<< " "<<R(1,2)<< " "
//                          <<R(2,0)<<" "<< R(2,1)<< " "<<R(2,2)<< std::endl;
//      }
//      out_show_camera.close();
  }





  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }





  return true;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Process_threeCameras_imu_gps_init_new_second(std::string rtFilePath, Mat3 changeR, Vec3 changet, double changeS) {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if(set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

  Hash_Map<IndexT, Mat3> global_rotations;

  //init inter-constraint
  Mat3 Rb_trans, Rf_trans, Rd_trans;
  Vec3 tb_trans, tf_trans, td_trans;
  Mat3 R_imu_trans;
  Vec3 gps_trans;
  std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
  const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
  const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};
  {
      //test distortion
      Mat3 r_tb, r_tf;
      r_tb << 1.0, 0.0, 0.0,
              0.0, 0.787926, -0.61577,
              0.0, 0.61577, 0.787926;
      r_tf << 1.0, 0.0, 0.0,
              0.0, 0.78796,  0.615727,
              0.0, -0.615727, 0.78796;

//      Rx << 1.0, 0.0, 0.0,
//             0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//             0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


//      r_tb << 1.0, 0.0, 0.0,
//              0.0, cos(38.0/180.0*M_PI), -sin(38.0/180.0*M_PI),
//              0.0, sin(38.0/180.0*M_PI), cos(38.0/180.0*M_PI);

//      r_tf << 1.0, 0.0, 0.0,
//              0.0, cos(-38.0/180.0*M_PI), -sin(-38.0/180.0*M_PI),
//              0.0, sin(-38.0/180.0*M_PI), cos(-38.0/180.0*M_PI);

//      r_tb << cos(-38.0/180.0*M_PI), 0.0, sin(-38.0/180.0*M_PI),//-38
//              0.0, 1.0, 0.0,
//              -sin(-38.0/180.0*M_PI), 0.0, cos(-38.0/180.0*M_PI);
//      r_tf << cos(38.0/180.0*M_PI), 0.0, sin(38.0/180.0*M_PI),//38
//              0.0, 1.0, 0.0,
//              -sin(38.0/180.0*M_PI), 0.0, cos(38.0/180.0*M_PI);
      Vec3 RA_tb, RA_tf;
      ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
      ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

      //set reference
      double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), -RA_tb(2)};
      double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), -RA_tf(2)};
//      double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), RA_tb(2)};
//      double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), RA_tf(2)};
      double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};

      Vec3 c_bd, c_fd;//b-d, f-d
      c_bd << 0.0, 0.098, 0.027;
      c_fd << 0.0, -0.101, 0.028;
//      c_bd << -0.098, 0.0, -0.027;
//      c_fd << 0.101, 0.0, -0.028;

//      c_bd << 0.098, 0.0, -0.027;
//      c_fd << -0.101, 0.0, -0.028;
//      c_bd << 0.0, 1.0, 0.027;
//      c_fd << 0.0, -0.101, 0.028;

      Vec3 t_db, t_df, t_dd;
      t_db = -r_tb*c_bd;
      t_df = -r_tf*c_fd;

      //double S = 1470.531783220221;
    //        t_db << 0.0, 0.098, 0.027;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
    //        t_df << 0.0, -0.101, 0.028;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
      t_dd << 0.0, 0.0, 0.0;

      transformation_br = {R_br_angleAxis[0], R_br_angleAxis[1], R_br_angleAxis[2],
                           t_db(0), t_db(1), t_db(2)};
      transformation_fr = {R_fr_angleAxis[0], R_fr_angleAxis[1], R_fr_angleAxis[2],
                           t_df(0), t_df(1), t_df(2)};
      transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                           t_dd(0), t_dd(1), t_dd(2)};

//      const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
//      const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};

      std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                             << tramsformation_x_gps[1] << " "
                                             << tramsformation_x_gps[2] << std::endl;

      //prepare for setting pose
      ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rb_trans.data());
      tb_trans = Vec3(transformation_br[3], transformation_br[4], transformation_br[5]);
      ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rf_trans.data());
      tf_trans = Vec3(transformation_fr[3], transformation_fr[4], transformation_fr[5]);
      ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rd_trans.data());
      td_trans = Vec3(transformation_dr[3], transformation_dr[4], transformation_dr[5]);


      ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
      gps_trans = Vec3(tramsformation_x_gps[0], tramsformation_x_gps[1], tramsformation_x_gps[2]);

  }


//  Mat3 changeR;
//  Vec3 changet;
//  double s;
//  {
//      //folder 16 Z = 0
//      changeR << 0.999116, 0.0420285, -0.0012437,
//           -0.0420389, 0.999065, -0.0101322,
//           0.000816696, 0.0101755, 0.999948;
//      changet = Vec3(-778.693, -79.3427, -780.618);
//      s = 0.98445;//0.30558845;//1.0;//3.1908;
//  }


  //get gps and imu data
  PoseStatusMap poseStatusMap;
  Hash_Map<IndexT, std::vector<double> > C_gps_Map;
  Hash_Map<IndexT, std::vector<double> > R_imu_Map;
  std::vector<double> c0(3);
  {
    //        /Hash_Map<IndexT, std::vector<double> > R_imu;
//            Hash_Map<IndexT, std::vector<double> > C_gps;
      std::string tFilePath = rtFilePath;  //!!!
      std::ifstream in;
      in.open(tFilePath);
      if(!in)
      {
          std::cout << "open tFilePath file error!" << std::endl;
          return false;
      }
      std::size_t line = 0;
      int count = -1;
      in >> count;  ///!!!

      bool isC0 = false;//true;

      while((!in.eof()) && (line < count))
      {
          // input
          int imgId, imgStatus;
          std::vector<double> c(3);
          //double pitch, roll, yaw;
          double omega, phi, kappa;
          in >> imgId
             >> c[0] >> c[1] >> c[2]
             >> omega >> phi >> kappa
             >> imgStatus;

          if(isC0 == false)
          {
              c0 = c;
              isC0 = true;
          }

//          c[0] -= c0[0];
//          c[1] -= c0[1];
//          c[2] -= c0[2];

          ++line;



//          double roll, pitch, yaw;

//          roll = omega;
//          pitch = phi;
//          yaw = kappa;

//          Mat3 Rz, Ry, Rx;
//          Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
//                sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
//                0.0, 1.0, 0.0,
//                -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
//          Rx << 1.0, 0.0, 0.0,
//                0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
//                0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);
//          Mat3 R_imu =  Rz*Rx*Ry;
//          Mat3 changeAxis_M;
//          changeAxis_M << 1.0, 0.0, 0.0,
//                          0.0, -1.0, 0.0,
//                          0.0, 0.0, -1.0;
//          R_imu = changeAxis_M * R_imu;




//          Mat3 Rz, Ry, Rx;
//          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
//                0.0, 1.0, 0.0,
//                -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);
//          Rx << 1.0, 0.0, 0.0,
//                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

////          Mat3 R_imu =  Rx*Ry*Rz;
////          Mat3 R_imu =  Rx*Rz*Ry;//
////          Mat3 R_imu =  Ry*Rz*Rx;
////          Mat3 R_imu =  Ry*Rx*Rz;
//          Mat3 R_imu =  Rz*Rx*Ry;//
////          Mat3 R_imu =  Rz*Ry*Rx;


          //use now
          Mat3 Rz, Ry, Rx;
          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
                0.0, 0.0, 1.0;
          Ry << cos(phi/180.0*M_PI), 0.0, sin(phi/180.0*M_PI),
                0.0, 1.0, 0.0,
                -sin(phi/180.0*M_PI), 0.0, cos(phi/180.0*M_PI);
          Rx << 1.0, 0.0, 0.0,
                0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
                0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
//          Mat3 R_imu =  Rz*Ry*Rx;
          Mat3 R_imu =  Rx*Ry*Rz;



//          Mat3 Rz, Ry, Rx;
//          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
//                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ry << 1.0, 0.0, 0.0,
//                0.0, -cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//                0.0, sin(phi/180.0*M_PI), -cos(phi/180.0*M_PI);
//          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;

//          Rz << cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI), 0.0,
//                sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ry << 1.0, 0.0, 0.0,
//                0.0, -cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
//                0.0, sin(omega/180.0*M_PI), -cos(omega/180.0*M_PI);
//          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;

//          Mat3 Rz, Ry, Rx;
//          Rz << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
//                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ry << 1.0, 0.0, 0.0,
//                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);
//          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;          Mat3 Rz, Ry, Rx;
//          Mat3 Rz, Ry, Rx;
//          Rz << cos(phi/180.0*M_PI), sin(phi/180.0*M_PI), 0.0,
//                -sin(phi/180.0*M_PI), cos(phi/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ry << 1.0, 0.0, 0.0,
//                0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
//                0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
//          Rx << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;



//          Mat3 R_imu =  Ry*Rz*Rx;
//          Mat3 R_imu =  Rz*Ry*Rx;
//          Mat3 R_imu =  Rx*Ry*Rz;


//          Mat3 Rk, Ra, RA;
//          Rk << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          RA << cos(omega/180.0*M_PI), sin(omega/180.0*M_PI), 0.0,
//                -sin(omega/180.0*M_PI), cos(omega/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ra << 1.0, 0.0, 0.0,
//                0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//                0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);


//          Mat3 R_imu =  Ry*Rx*Rz;

//          Mat3 changeAxis_M;
//          changeAxis_M << 1.0, 0.0, 0.0,
//                          0.0, -1.0, 0.0,
//                          0.0, 0.0, -1.0;  //M1
//          changeAxis_M << 0.0, 1.0, 0.0,
//                          1.0, 0.0, 0.0,
//                          0.0, 0.0, -1.0;  //M2

//          R_imu = changeAxis_M * R_imu;
//          R_imu = R_imu * changeAxis_M;





//              tempGet = getListXYZ(-phi3, kappa3, -omega3);
//          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

//               tempGet = getListXYZ(phi3, kappa3, omega3);
//          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

//          tempGet = getListXYZ(phi1, kappa1, omega1);
//          Eigen::Matrix3d Rres2 = getRx(a1)*getRz(a2)*getRy(a3);

//          Mat3 Rz, Ry, Rx;
//          Rx << 1.0, 0.0, 0.0,
//                 0.0, cos(phi/180.0*M_PI), -sin(phi/180.0*M_PI),
//                 0.0, sin(phi/180.0*M_PI), cos(phi/180.0*M_PI);

//          Ry << cos(omega/180.0*M_PI), 0.0, sin(omega/180.0*M_PI),
//                 0.0, 1.0, 0.0,
//                 -sin(omega/180.0*M_PI), 0.0, cos(omega/180.0*M_PI);

//          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                 sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                 0.0, 0.0, 1.0;


//          Mat3 R_imu =  Ry*Rz*Rx;
//          Mat3 R_imu =  Rz*Rx*Ry;
//          Mat3 R_imu =  Rx*Rz*Ry;
//          Mat3 R_imu =  Rx*Ry*Rz;

          Mat3 changeAxis_M;
          changeAxis_M << 1.0, 0.0, 0.0,
                          0.0, -1.0, 0.0,
                          0.0, 0.0, -1.0;
//          changeAxis_M << 1.0, 0.0, 0.0,
//                          0.0, 1.0, 0.0,
//                          0.0, 0.0, -1.0;

//          changeAxis_M << 1.0, 0.0, 0.0,
//                          0.0, 1.0, 0.0,
//                          0.0, 0.0, -1.0;
          Mat3 R_;
          R_ << 1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, -1.0;
          R_imu = changeAxis_M * R_imu;//* R_;


          //Mat3 R2;
          //ceres::AngleAxisToRotationMatrix(R2_a.data(), R2.data());

          //R_imu = changeAxis_M * R_imu * R2;//* R_;


          std::vector<double> R_imu_a(3);
          ceres::RotationMatrixToAngleAxis(R_imu.data(), R_imu_a.data());

          Vec3 C_gps_res = Vec3(c[0], c[1], c[2]);

          // check status
          if(imgId/100000 == 1)//back
          {
              if(imgStatus != 321)
              {
                  C_gps_Map[imgId + 100000] = c;
                  R_imu_Map[imgId + 100000] = R_imu_a;
              }

              //set back-camera pose
              Mat3 Rb_res = Rb_trans * R_imu_trans * R_imu;
              //!!
              Mat3 Rb_res2 = (1/changeS)*Rb_res * changeR.inverse();

              Vec3 tb_res;
//              tb_res = Rb_trans * gps_trans - Rb_trans * R_imu_trans * R_imu * C_gps_res + tb_trans;
              tb_res = Rb_trans * gps_trans - (1/changeS) * Rb_trans * R_imu_trans * R_imu  * changeR.inverse()* C_gps_res + tb_trans;
              global_rotations[imgId] = Rb_res2;

//              tb_res = tb_res - Rb_res*changet;

//              sfm_data_.poses[imgId] = Pose3(Rb_res, -Rb_res.transpose() * tb_res);
              sfm_data_.poses[imgId] = Pose3(Rb_res2, -Rb_res2.transpose() * tb_res);

              PoseStatus poseStatus;
              poseStatus.pose = sfm_data_.poses[imgId];
              poseStatus.status = imgStatus;
              poseStatus.bdf = 1;
              poseStatusMap[imgId] = poseStatus;

          }else if(imgId/100000 == 2)//down
          {
              C_gps_Map[imgId] = c;
              R_imu_Map[imgId] = R_imu_a;

              //set down-camera pose
              Mat3 Rd_res = Rd_trans * R_imu_trans * R_imu;
              //!!
              Mat3 Rd_res2 = (1/changeS)*Rd_res * changeR.inverse();

              Vec3 td_res;
//              td_res = Rd_trans * gps_trans - Rd_trans * R_imu_trans * R_imu * C_gps_res + td_trans;
              td_res = Rd_trans * gps_trans - (1/changeS)*Rd_trans * R_imu_trans * R_imu  * changeR.inverse()* C_gps_res + td_trans;
              global_rotations[imgId] = Rd_res2;
              //!!
//              td_res = td_res - Rd_res*changet;

              sfm_data_.poses[imgId] = Pose3(Rd_res2, -Rd_res2.transpose() * td_res);

              PoseStatus poseStatus;
              poseStatus.pose = sfm_data_.poses[imgId];
              poseStatus.status = imgStatus;
              poseStatus.bdf = 2;
              poseStatusMap[imgId] = poseStatus;


          }else if(imgId/100000 == 3)//front
          {
              if(imgStatus != 321)
              {
                  C_gps_Map[imgId - 100000] = c;
                  R_imu_Map[imgId - 100000] = R_imu_a;
              }

              //set front-camera pose
              Mat3 Rf_res = Rf_trans * R_imu_trans * R_imu;
              //!!
              Mat3 Rf_res2 = (1/changeS)*Rf_res * changeR.inverse();

              Vec3 tf_res;
//              tf_res = Rf_trans * gps_trans - Rf_trans * R_imu_trans * R_imu * C_gps_res + tf_trans;
              tf_res = Rf_trans * gps_trans - (1/changeS)*Rf_trans * R_imu_trans * R_imu  * changeR.inverse()* C_gps_res + tf_trans;
              global_rotations[imgId] = Rf_res2;
              //!!
//              tf_res = tf_res - Rf_res*changet;

              sfm_data_.poses[imgId] = Pose3(Rf_res2, -Rf_res2.transpose() * tf_res);

              PoseStatus poseStatus;
              poseStatus.pose = sfm_data_.poses[imgId];
              poseStatus.status = imgStatus;
              poseStatus.bdf = 3;
              poseStatusMap[imgId] = poseStatus;

          }

      }
      in.close();

  }

  //Refine_initValue();
  matching::PairWiseMatches tripletWise_matches;
  if (!Compute_Global_Translations_gps(global_rotations, tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
    return false;
  }


  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }


  if (!Adjust_init())//_fix2())//(!Adjust_threecamera_gps_cail_new(poseStatusMap, transformation_R_imu, tramsformation_x_gps, transformation_br, transformation_fr, C_gps_Map, R_imu_Map))
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }


  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }





  return true;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Process_threeCameras_imu_gps_init_new_test() {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if(set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

  Hash_Map<IndexT, Mat3> global_rotations;

  for (Poses::const_iterator itPose = sfm_data_.poses.begin(); itPose != sfm_data_.poses.end(); ++itPose)
  {
      const Mat3 Rthis = itPose->second.rotation();
      global_rotations[itPose->first] = Rthis;
  }

  //Refine_initValue();
  matching::PairWiseMatches tripletWise_matches;
  if (!Compute_Global_Translations_gps(global_rotations, tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
    return false;
  }


  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }


  if (!Adjust_init())//(!Adjust_threecamera_gps_cail_new(poseStatusMap, transformation_R_imu, tramsformation_x_gps, transformation_br, transformation_fr, C_gps_Map, R_imu_Map))
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }


  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }





  return true;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Process_test(std::string MatchesDir) {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if(set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

  Hash_Map<IndexT, Mat3> global_rotations;

  //init inter-constraint
  Mat3 Rb_trans, Rf_trans, Rd_trans;
  Vec3 tb_trans, tf_trans, td_trans;
  Mat3 R_imu_trans;
  Vec3 gps_trans;
  std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
  const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
  const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};
  {
      //test distortion
      Mat3 r_tb, r_tf;
      r_tb << 1.0, 0.0, 0.0,
              0.0, 0.787926, -0.61577,
              0.0, 0.61577, 0.787926;
      r_tf << 1.0, 0.0, 0.0,
              0.0, 0.78796,  0.615727,
              0.0, -0.615727, 0.78796;
      Vec3 RA_tb, RA_tf;
      ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
      ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

      //set reference
      double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), -RA_tb(2)};
      double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), -RA_tf(2)};
      double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};

      Vec3 c_bd, c_fd;//b-d, f-d
      c_bd << 0.0, 0.098, 0.027;
      c_fd << 0.0, -0.101, 0.028;

      Vec3 t_db, t_df, t_dd;
      t_db = -r_tb*c_bd;
      t_df = -r_tf*c_fd;

      //double S = 1470.531783220221;
    //        t_db << 0.0, 0.098, 0.027;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
    //        t_df << 0.0, -0.101, 0.028;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
      t_dd << 0.0, 0.0, 0.0;

      transformation_br = {R_br_angleAxis[0], R_br_angleAxis[1], R_br_angleAxis[2],
                           t_db(0), t_db(1), t_db(2)};
      transformation_fr = {R_fr_angleAxis[0], R_fr_angleAxis[1], R_fr_angleAxis[2],
                           t_df(0), t_df(1), t_df(2)};
      transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                           t_dd(0), t_dd(1), t_dd(2)};

//      const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
//      const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};

      std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                             << tramsformation_x_gps[1] << " "
                                             << tramsformation_x_gps[2] << std::endl;

      //prepare for setting pose
      ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rb_trans.data());
      tb_trans = Vec3(transformation_br[3], transformation_br[4], transformation_br[5]);
      ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rf_trans.data());
      tf_trans = Vec3(transformation_fr[3], transformation_fr[4], transformation_fr[5]);
      ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rd_trans.data());
      td_trans = Vec3(transformation_dr[3], transformation_dr[4], transformation_dr[5]);


      ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
      gps_trans = Vec3(tramsformation_x_gps[0], tramsformation_x_gps[1], tramsformation_x_gps[2]);

  }

  //get gps and imu data
  PoseStatusMap poseStatusMap;
  Hash_Map<IndexT, std::vector<double> > C_gps_Map;
  Hash_Map<IndexT, std::vector<double> > R_imu_Map;
  {
//      std::string tFilePath = sfm_data_.s_root_path + "/n_gps_imu_321.res";  //!!!
//      std::ifstream in;
//      in.open(tFilePath);
//      if(!in)
//      {
//          std::cout << "open tFilePath file error!" << std::endl;
//          return false;
//      }
//      std::size_t line = 0;
//      int count = -1;
//      in >> count;  ///!!!

//      bool isC0 = false;
//      std::vector<double> c0(3);
//      while((!in.eof()) && (line < count))
//      {
//          // input
//          int imgId, imgStatus;
//          std::vector<double> c(3);
//          //double pitch, roll, yaw;
//          double omega, phi, kappa;
//          in >> imgId
//             >> c[0] >> c[1] >> c[2]
//             >> omega >> phi >> kappa
//             >> imgStatus;

//          if(isC0 == false)
//          {
//              c0 = c;
//              isC0 = true;
//          }

//          c[0] -= c0[0];
//          c[1] -= c0[1];
//          c[2] -= c0[2];

//          ++line;

//          Mat3 Rz, Ry, Rx;
//          Rz << cos(kappa/180.0*M_PI), -sin(kappa/180.0*M_PI), 0.0,
//                sin(kappa/180.0*M_PI), cos(kappa/180.0*M_PI), 0.0,
//                0.0, 0.0, 1.0;
//          Ry << cos(phi/180.0*M_PI), 0.0, -sin(phi/180.0*M_PI),
//                0.0, 1.0, 0.0,
//                sin(phi/180.0*M_PI), 0.0, cos(phi/180.0*M_PI);
//          Rx << 1.0, 0.0, 0.0,
//                0.0, cos(omega/180.0*M_PI), -sin(omega/180.0*M_PI),
//                0.0, sin(omega/180.0*M_PI), cos(omega/180.0*M_PI);
//          Mat3 R_imu =  Ry*Rx*Rz;

//          Mat3 changeAxis_M;
//          changeAxis_M << 1.0, 0.0, 0.0,
//                          0.0, -1.0, 0.0,
//                          0.0, 0.0, -1.0;
//          R_imu = changeAxis_M * R_imu;

//          std::vector<double> R_imu_a(3);
//          ceres::RotationMatrixToAngleAxis(R_imu.data(), R_imu_a.data());

//          Vec3 C_gps_res = Vec3(c[0], c[1], c[2]);

//          // check status
//          if(imgId/100000 == 1)//back
//          {
//              if(imgStatus != 321)
//              {
//                  C_gps_Map[imgId + 100000] = c;
//                  R_imu_Map[imgId + 100000] = R_imu_a;
//              }

//              //set back-camera pose
//              Mat3 Rb_res = Rb_trans * R_imu_trans * R_imu;
//              Vec3 tb_res;
//              tb_res = Rb_trans * gps_trans - Rb_trans * R_imu_trans * R_imu * C_gps_res + tb_trans;
//              global_rotations[imgId] = Rb_res;
//              sfm_data_.poses[imgId] = Pose3(Rb_res, -Rb_res.transpose() * tb_res);

//              PoseStatus poseStatus;
//              poseStatus.pose = sfm_data_.poses[imgId];
//              poseStatus.status = imgStatus;
//              poseStatus.bdf = 1;
//              poseStatusMap[imgId] = poseStatus;

//          }else if(imgId/100000 == 2)//down
//          {
//              C_gps_Map[imgId] = c;
//              R_imu_Map[imgId] = R_imu_a;

//              //set down-camera pose
//              Mat3 Rd_res = Rd_trans * R_imu_trans * R_imu;
//              Vec3 td_res;
//              td_res = Rd_trans * gps_trans - Rd_trans * R_imu_trans * R_imu * C_gps_res + td_trans;
//              global_rotations[imgId] = Rd_res;
//              sfm_data_.poses[imgId] = Pose3(Rd_res, -Rd_res.transpose() * td_res);

//              PoseStatus poseStatus;
//              poseStatus.pose = sfm_data_.poses[imgId];
//              poseStatus.status = imgStatus;
//              poseStatus.bdf = 2;
//              poseStatusMap[imgId] = poseStatus;

//          }else if(imgId/100000 == 3)//front
//          {
//              if(imgStatus != 321)
//              {
//                  C_gps_Map[imgId - 100000] = c;
//                  R_imu_Map[imgId - 100000] = R_imu_a;
//              }

//              //set front-camera pose
//              Mat3 Rf_res = Rf_trans * R_imu_trans * R_imu;
//              Vec3 tf_res;
//              tf_res = Rf_trans * gps_trans - Rf_trans * R_imu_trans * R_imu * C_gps_res + tf_trans;
//              global_rotations[imgId] = Rf_res;
//              sfm_data_.poses[imgId] = Pose3(Rf_res, -Rf_res.transpose() * tf_res);

//              PoseStatus poseStatus;
//              poseStatus.pose = sfm_data_.poses[imgId];
//              poseStatus.status = imgStatus;
//              poseStatus.bdf = 3;
//              poseStatusMap[imgId] = poseStatus;

//          }

//      }
//      in.close();

  }


  {
      std::string tFilePath = "/home/guang/data/beichuan_testBA4/sfm3/out_show_camera_c.txt";  //!!!
      std::ifstream in;
      in.open(tFilePath);
      if(!in)
      {
          std::cout << "open tFilePath file error!" << std::endl;
          return false;
      }
      std::size_t line = 0;
      int count = -1;
      in >> count;  ///!!!

      bool isC0 = false;
      std::vector<double> c0(3);
      while((!in.eof()) && (line < count))
      {
          int id;
          double c0, c1, c2, r0, r1, r2, r3, r4, r5, r6, r7, r8;
          in >> id >>c0 >> c1 >>c2>>r0>> r1>> r2>> r3>> r4>> r5>> r6>> r7>> r8;

          Vec3 c;
          c<< c0, c1, c2;
          Mat3 R;
          R <<r0, r1, r2, r3, r4, r5, r6, r7, r8;
          sfm_data_.poses[id] = Pose3(R, c);

      }
      in.close();
  }

  {
//    std::cout
//      << "=============================================================\n"
//      << "Robust triangulation of the tracks\n"
//      << " - Triangulation of guided epipolar geometry matches\n"
//      << "============================================================="
//      << std::endl;
//    //--
//    //- Pair selection method:
//    //  - geometry guided -> camera frustum intersection,
//    //  - putative matches guided (photometric matches)
//    //     (keep pairs that have valid Intrinsic & Pose ids).
//    //--
//    Pair_Set pairs = matches_provider_->getPairs();
//    std::cout << "error 1" << std::endl;
//    //const std::set<IndexT> valid_viewIdx = Get_Valid_Views(sfm_data_);
//    //pairs = Pair_filter(pairs, valid_viewIdx);
//    std::cout << "error 2" << std::endl;


//    openMVG::system::Timer timer;

//    //------------------------------------------
//    // Compute Structure from known camera poses
//    //------------------------------------------
//    double dMax_reprojection_error = 4.0;
//    std::shared_ptr<Regions_Provider> regions_provider;
//    regions_provider = std::make_shared<Regions_Provider>();
//    C_Progress_display progress;
//    const std::string sImage_describer = stlplus::create_filespec(MatchesDir, "image_describer", "json");
//    std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
//    if (!regions_provider->load(sfm_data_, MatchesDir, regions_type, &progress)) {
//      std::cerr << std::endl
//        << "Invalid regions." << std::endl;
//      return EXIT_FAILURE;
//    }

//    std::cout << "error 3" << std::endl;
//    SfM_Data_Structure_Estimation_From_Known_Poses structure_estimator(dMax_reprojection_error);
//    std::cout << "error 4" << std::endl;
//    std::cout << pairs.size() << std::endl;
//    structure_estimator.run(sfm_data_, pairs, regions_provider);
//    std::cout << "error 5" << std::endl;
//    std::cout << "\nStructure estimation took (s): " << timer.elapsed() << "." << std::endl;

//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "from_pose", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));

//    regions_provider.reset(); // Regions are not longer needed.
//    RemoveOutliers_AngleError(sfm_data_, 2.0);

  }

  {
      //Refine_initValue();
      matching::PairWiseMatches tripletWise_matches;
      if (!Compute_Global_Translations_gps(global_rotations, tripletWise_matches))
      {
        std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
        return false;
      }

      if (!Compute_Initial_Structure(tripletWise_matches))
      {
        std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
        return false;
      }
  }


  Save(sfm_data_,
    stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "from_pose_remove", "ply"),
    ESfM_Data(EXTRINSICS | STRUCTURE));

  // Check that poses & intrinsic cover some measures (after outlier removal)
  const IndexT minPointPerPose = 12; // 6 min
  const IndexT minTrackLength = 3; // 2 min
  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
  {
    KeepLargestViewCCTracks(sfm_data_);
    eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength);

    const size_t pointcount_cleaning = sfm_data_.structure.size();
    std::cout << "Point_cloud cleaning:\n"
      << "\t #3DPoints: " << pointcount_cleaning << "\n";
  }

  Save(sfm_data_,
    stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "from_pose_remove2", "ply"),
    ESfM_Data(EXTRINSICS | STRUCTURE));

    {
//            std::cout<< "initial intrinsics :"<<std::endl;
//                Hash_Map<IndexT, std::vector<double> > map_intrinsics;
//                map_intrinsics[0] = sfm_data_.intrinsics[0]->getParams();
//                std::cout<< "map_intrinsics 0 " << map_intrinsics[0][0]<<" "<<map_intrinsics[0][1]<<" "<<
//                         map_intrinsics[0][2]<<" "<<map_intrinsics[0][3]<<" "<<map_intrinsics[0][4]<<" "<<
//                         map_intrinsics[0][5]<<" "<<map_intrinsics[0][6]<<" "<<map_intrinsics[0][7]<<std::endl;

//                map_intrinsics[1] = sfm_data_.intrinsics[1]->getParams();
//                std::cout<< "map_intrinsics 1 " << map_intrinsics[1][0]<<" "<<map_intrinsics[1][1]<<" "<<
//                         map_intrinsics[1][2]<<" "<<map_intrinsics[1][3]<<" "<<map_intrinsics[1][4]<<" "<<
//                         map_intrinsics[1][5]<<" "<<map_intrinsics[1][6]<<" "<<map_intrinsics[1][7]<<std::endl;

//                map_intrinsics[2] = sfm_data_.intrinsics[2]->getParams();
//                std::cout<< "map_intrinsics 2 " << map_intrinsics[2][0]<<" "<<map_intrinsics[2][1]<<" "<<
//                         map_intrinsics[2][2]<<" "<<map_intrinsics[2][3]<<" "<<map_intrinsics[2][4]<<" "<<
//                         map_intrinsics[2][5]<<" "<<map_intrinsics[2][6]<<" "<<map_intrinsics[2][7]<<std::endl;
//                std::cout <<"-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() <<std::endl;
    }

  if (!Adjust_init())//_init())//(!Adjust_threecamera_gps_cail_new(poseStatusMap, transformation_R_imu, tramsformation_x_gps, transformation_br, transformation_fr, C_gps_Map, R_imu_Map))
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }




  {
//          std::cout<< "map_intrinsics 1"<<std::endl;
//              Hash_Map<IndexT, std::vector<double> > map_intrinsics;
//              map_intrinsics[0] = sfm_data_.intrinsics[0]->getParams();
//              std::cout<< "map_intrinsics 0 " << map_intrinsics[0][0]<<" "<<map_intrinsics[0][1]<<" "<<
//                       map_intrinsics[0][2]<<" "<<map_intrinsics[0][3]<<" "<<map_intrinsics[0][4]<<" "<<
//                       map_intrinsics[0][5]<<" "<<map_intrinsics[0][6]<<" "<<map_intrinsics[0][7]<<std::endl;

//              map_intrinsics[1] = sfm_data_.intrinsics[1]->getParams();
//              std::cout<< "map_intrinsics 1 " << map_intrinsics[1][0]<<" "<<map_intrinsics[1][1]<<" "<<
//                       map_intrinsics[1][2]<<" "<<map_intrinsics[1][3]<<" "<<map_intrinsics[1][4]<<" "<<
//                       map_intrinsics[1][5]<<" "<<map_intrinsics[1][6]<<" "<<map_intrinsics[1][7]<<std::endl;

//              map_intrinsics[2] = sfm_data_.intrinsics[2]->getParams();
//              std::cout<< "map_intrinsics 2 " << map_intrinsics[2][0]<<" "<<map_intrinsics[2][1]<<" "<<
//                       map_intrinsics[2][2]<<" "<<map_intrinsics[2][3]<<" "<<map_intrinsics[2][4]<<" "<<
//                       map_intrinsics[2][5]<<" "<<map_intrinsics[2][6]<<" "<<map_intrinsics[2][7]<<std::endl;
//              std::cout <<"-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() <<std::endl;
  }


  //out
  {
//      std::ofstream out_show_camera;
//      out_show_camera.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_reference.txt");
//      std::ofstream out_Rc;
//      out_Rc.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_Rc.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 c = pose.center();

//          Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
//          axis_x = R.transpose() * (axis_x + R*c);
//          axis_y = R.transpose() * (axis_y + R*c);
//          axis_z = R.transpose() * (axis_z + R*c);
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_x(0) << " " << axis_x(1) << " " << axis_x(2)<< ";" << std::endl;
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_y(0) << " " << axis_y(1) << " " << axis_y(2)<< ";" << std::endl;
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_z(0) << " " << axis_z(1) << " " << axis_z(2)<< ";" << std::endl;

//          out_Rc << c(0) << " " << c(1)<< " " << c(2) << " " << std::endl;

//      }
//      out_Rc.close();
//      out_show_camera.close();

  }

  //out
  {
//      std::ofstream out_c;
//      out_c.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_c.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Vec3 c = pose.center();

//          out_c << c(0) << " " << c(1)<< " " << c(2) << std::endl;

//      }
//      out_c.close();

  }

  {
//      std::ofstream out_show_camera;
//      out_show_camera.open("/media/add7/E/run_scene/testMerge/29/out_camera_res_R.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 c = pose.center();

//          Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
//          axis_x = R.transpose() * (axis_x + R*c);
//          axis_y = R.transpose() * (axis_y + R*c);
//          axis_z = R.transpose() * (axis_z + R*c);
//          out_show_camera <<R(0,0)<<" "<< R(0,1)<< " "<<R(0,2)<< " "
//                          <<R(1,0)<<" "<< R(1,1)<< " "<<R(1,2)<< " "
//                          <<R(2,0)<<" "<< R(2,1)<< " "<<R(2,2)<< std::endl;
//      }
//      out_show_camera.close();
  }





  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }





  return true;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Process_threeCameras_imu_gps_c4(std::string sMatchesDir) {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if(set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

//    openMVG::rotation_averaging::RelativeRotations relatives_R;
//    Compute_Relative_Rotations(relatives_R);

    Hash_Map<IndexT, Mat3> global_rotations;

    //test distortion
    Mat3 r_tb, r_tf;
    r_tb << 1.0, 0.0, 0.0,
            0.0, 0.787926, -0.61577,
            0.0, 0.61577, 0.787926;
    r_tf << 1.0, 0.0, 0.0,
            0.0, 0.78796,  0.615727,
            0.0, -0.615727, 0.78796;
    Vec3 RA_tb, RA_tf;
    ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
    ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

        //set reference
        double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), -RA_tb(2)};//{0.665224, 0.00448363, 0.00291885};  //{0.703355, 0.00263495, -0.00149232};//{0.0, 0.0, 0.0}; //{0.0212232, 0.0206259, 0.0420509};  //{0.0, 0.0, 0.0}; //{0.669661, 0.00354668, 0.00220508}; //  //{0.0, 0.0, 0.0}; //{0.662013, -0.00281273, 0.000638113 }; //{0.0, 0.0, 0.0}; //{0.671924, 0.00432461, 0.000528427};
        double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), -RA_tf(2)};//{-0.644463, 0.0113345, -0.00156344}; //{-0.682501, 0.0147873, 0.00150509};//{0.0, 0.0, 0.0}; //{-0.0232707, -0.00288052, 0.00035795};  //{0.0, 0.0, 0.0}; //{-0.650191, 0.0101123, -0.000844533}; //  //{0.0, 0.0, 0.0}; //{-0.658446, 0.00532591, 0.0010695}; //{-0.656461, 0.00508743, -0.00299628};
        double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};
        Vec3 t_db, t_df, t_dd;
        //double S = 1470.531783220221;
        t_db << 0.0, 0.0, 0.0;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
        t_df << 0.0, 0.0, 0.0;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
        t_dd << 0.0, 0.0, 0.0;


        std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
        transformation_br = {R_br_angleAxis[0], R_br_angleAxis[1], R_br_angleAxis[2],
                             t_db(0), t_db(1), t_db(2)};
        transformation_fr = {R_fr_angleAxis[0], R_fr_angleAxis[1], R_fr_angleAxis[2],
                             t_df(0), t_df(1), t_df(2)};
        transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                             t_dd(0), t_dd(1), t_dd(2)};


    const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
    const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};


        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;


        //get gps and imu data
      {
    //        /Hash_Map<IndexT, std::vector<double> > R_imu;
            Hash_Map<IndexT, std::vector<double> > C_gps;
            {
                //read file and initialise external
                std::vector<double> c(3);
                std::string tFilePath = sfm_data_.s_root_path + "/latlon_c.res";  //!!!
                std::ifstream in;
                in.open(tFilePath);
                if(!in)
                {
                    std::cout << "open tFilePath file error!" << std::endl;
                    return false;
                }
                std::size_t line = 0;
                const int count = sfm_data_.views.size()/3;
                while((!in.eof()) && (line < count))
                {
                     in >> c[0] >> c[1] >> c[2];
                     C_gps[line + count] = {c[0], c[1], c[2]};
                    ++line;
                }
                in.close();

            }
            //read imu
            {
            std::ifstream in;
            in.open(sfm_data_.s_root_path + "/IMU.res"); //imu file

            int line = 0;
            const int count = sfm_data_.views.size()/3;

    //        std::ofstream out_show_camera;
    //        out_show_camera.open(sfm_data_.s_root_path + "/show_imu_R.txt");

            while((!in.eof()) && (line < count))
            {
                double roll, pitch, yaw;
                in >> roll >> pitch >> yaw;  //z x y

                Mat3 Rz, Ry, Rx;
                Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
                      sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
                      0.0, 0.0, 1.0;
                Ry << cos(pitch/180.0*M_PI), 0.0, sin(pitch/180.0*M_PI),
                      0.0, 1.0, 0.0,
                      -sin(pitch/180.0*M_PI), 0.0, cos(pitch/180.0*M_PI);
                Rx << 1.0, 0.0, 0.0,
                      0.0, cos(roll/180.0*M_PI), -sin(roll/180.0*M_PI),
                      0.0, sin(roll/180.0*M_PI), cos(roll/180.0*M_PI);
                Mat3 R_imu =  Rz*Rx*Ry;
                Mat3 changeAxis_M;
                changeAxis_M << 1.0, 0.0, 0.0,
                                0.0, -1.0, 0.0,
                                0.0, 0.0, -1.0;
                R_imu = changeAxis_M * R_imu;


                {

    //                Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
    //                Vec3 c = {C_gps[line + count][0], C_gps[line + count][1], C_gps[line + count][2]};
    //                axis_x = R_d.transpose() * (axis_x + R_d*c);
    //                axis_y = R_d.transpose() * (axis_y + R_d*c);
    //                axis_z = R_d.transpose() * (axis_z + R_d*c);
    //                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
    //                       << axis_x(0) << " " << axis_x(1) << " " << axis_x(2)<< ";" << std::endl;
    //                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
    //                       << axis_y(0) << " " << axis_y(1) << " " << axis_y(2)<< ";" << std::endl;
    //                out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
    //                       << axis_z(0) << " " << axis_z(1) << " " << axis_z(2)<< ";" << std::endl;

    //                out_show_camera << R_d << std::endl << std::endl;
                }

    //            global_rotations[line + count] = R_d;
    //            sfm_data_.poses[line + count] = Pose3(R_d, C_gps[line + count]);

                Mat3 Rb_trans, Rf_trans, Rd_trans;
                Vec3 tb_trans, tf_trans, td_trans;

                ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rb_trans.data());
                tb_trans = Vec3(transformation_br[3], transformation_br[4], transformation_br[5]);
                ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rf_trans.data());
                tf_trans = Vec3(transformation_fr[3], transformation_fr[4], transformation_fr[5]);
                ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rd_trans.data());
                td_trans = Vec3(transformation_dr[3], transformation_dr[4], transformation_dr[5]);


                Mat3 R_imu_trans;
                ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
                Vec3 gps_trans = Vec3(tramsformation_x_gps[0], tramsformation_x_gps[1], tramsformation_x_gps[2]);
                Vec3 C_gps_res = Vec3(C_gps[line + count][0], C_gps[line + count][1], C_gps[line + count][2]);

                //set f
                Mat3 Rf_res = Rf_trans * R_imu_trans * R_imu;
                Vec3 tf_res;
                tf_res = Rf_trans * gps_trans - Rf_trans * R_imu_trans * R_imu * C_gps_res + tf_trans;
                global_rotations[line + count*2] = Rf_res;
                sfm_data_.poses[line + count*2] = Pose3(Rf_res, -Rf_res.transpose() * tf_res);

                //set d
                Mat3 Rd_res = Rd_trans * R_imu_trans * R_imu;
                Vec3 td_res;
                td_res = Rd_trans * gps_trans - Rd_trans * R_imu_trans * R_imu * C_gps_res + td_trans;
                global_rotations[line + count] = Rd_res;
                sfm_data_.poses[line + count] = Pose3(Rd_res, -Rd_res.transpose() * td_res);

                //set b
                Mat3 Rb_res = Rb_trans * R_imu_trans * R_imu;
                Vec3 tb_res;
                tb_res = Rb_trans * gps_trans - Rb_trans * R_imu_trans * R_imu * C_gps_res + tb_trans;
                global_rotations[line] = Rb_res;
                sfm_data_.poses[line] = Pose3(Rb_res, -Rb_res.transpose() * tb_res);

                ++line;
            }
            in.close();

    //        out_show_camera.close();
        }
      }


  matching::PairWiseMatches tripletWise_matches;
  if (!Compute_Global_Translations_gps(global_rotations, tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
    return false;
  }

    if (!Compute_Initial_Structure(tripletWise_matches))
    {
      std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
      return false;
    }

//    for( int i = 0; i < sfm_data_.poses.size(); ++i)
//    {
//        double R0_angleAxis[3];
//        Mat3 mid_R0 = sfm_data_.poses[i].rotation();
//        ceres::RotationMatrixToAngleAxis((const double*)mid_R0.data(), R0_angleAxis);

//        std::cout << "R0_angleAxis : " << R0_angleAxis[0] << " " << R0_angleAxis[1] << " " << R0_angleAxis[2] << std::endl;

//    }

    {
//            std::cout<< "initial intrinsics :"<<std::endl;
//                Hash_Map<IndexT, std::vector<double> > map_intrinsics;
//                map_intrinsics[0] = sfm_data_.intrinsics[0]->getParams();
//                std::cout<< "map_intrinsics 0 " << map_intrinsics[0][0]<<" "<<map_intrinsics[0][1]<<" "<<
//                         map_intrinsics[0][2]<<" "<<map_intrinsics[0][3]<<" "<<map_intrinsics[0][4]<<" "<<
//                         map_intrinsics[0][5]<<" "<<map_intrinsics[0][6]<<" "<<map_intrinsics[0][7]<<std::endl;

//                map_intrinsics[1] = sfm_data_.intrinsics[1]->getParams();
//                std::cout<< "map_intrinsics 1 " << map_intrinsics[1][0]<<" "<<map_intrinsics[1][1]<<" "<<
//                         map_intrinsics[1][2]<<" "<<map_intrinsics[1][3]<<" "<<map_intrinsics[1][4]<<" "<<
//                         map_intrinsics[1][5]<<" "<<map_intrinsics[1][6]<<" "<<map_intrinsics[1][7]<<std::endl;

//                map_intrinsics[2] = sfm_data_.intrinsics[2]->getParams();
//                std::cout<< "map_intrinsics 2 " << map_intrinsics[2][0]<<" "<<map_intrinsics[2][1]<<" "<<
//                         map_intrinsics[2][2]<<" "<<map_intrinsics[2][3]<<" "<<map_intrinsics[2][4]<<" "<<
//                         map_intrinsics[2][5]<<" "<<map_intrinsics[2][6]<<" "<<map_intrinsics[2][7]<<std::endl;
//                std::cout <<"-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() <<std::endl;
    }


  if (!Adjust_threecamera_gps_change_c4())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }




  {
//          std::cout<< "map_intrinsics 1"<<std::endl;
//              Hash_Map<IndexT, std::vector<double> > map_intrinsics;
//              map_intrinsics[0] = sfm_data_.intrinsics[0]->getParams();
//              std::cout<< "map_intrinsics 0 " << map_intrinsics[0][0]<<" "<<map_intrinsics[0][1]<<" "<<
//                       map_intrinsics[0][2]<<" "<<map_intrinsics[0][3]<<" "<<map_intrinsics[0][4]<<" "<<
//                       map_intrinsics[0][5]<<" "<<map_intrinsics[0][6]<<" "<<map_intrinsics[0][7]<<std::endl;

//              map_intrinsics[1] = sfm_data_.intrinsics[1]->getParams();
//              std::cout<< "map_intrinsics 1 " << map_intrinsics[1][0]<<" "<<map_intrinsics[1][1]<<" "<<
//                       map_intrinsics[1][2]<<" "<<map_intrinsics[1][3]<<" "<<map_intrinsics[1][4]<<" "<<
//                       map_intrinsics[1][5]<<" "<<map_intrinsics[1][6]<<" "<<map_intrinsics[1][7]<<std::endl;

//              map_intrinsics[2] = sfm_data_.intrinsics[2]->getParams();
//              std::cout<< "map_intrinsics 2 " << map_intrinsics[2][0]<<" "<<map_intrinsics[2][1]<<" "<<
//                       map_intrinsics[2][2]<<" "<<map_intrinsics[2][3]<<" "<<map_intrinsics[2][4]<<" "<<
//                       map_intrinsics[2][5]<<" "<<map_intrinsics[2][6]<<" "<<map_intrinsics[2][7]<<std::endl;
//              std::cout <<"-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() <<std::endl;
  }


  //out
  {
//      std::ofstream out_show_camera;
//      out_show_camera.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_reference.txt");
//      std::ofstream out_Rc;
//      out_Rc.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_Rc.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 c = pose.center();

//          Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
//          axis_x = R.transpose() * (axis_x + R*c);
//          axis_y = R.transpose() * (axis_y + R*c);
//          axis_z = R.transpose() * (axis_z + R*c);
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_x(0) << " " << axis_x(1) << " " << axis_x(2)<< ";" << std::endl;
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_y(0) << " " << axis_y(1) << " " << axis_y(2)<< ";" << std::endl;
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_z(0) << " " << axis_z(1) << " " << axis_z(2)<< ";" << std::endl;

//          out_Rc << c(0) << " " << c(1)<< " " << c(2) << " " << std::endl;

//      }
//      out_Rc.close();
//      out_show_camera.close();

  }

  //out
  {
//      std::ofstream out_c;
//      out_c.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_c.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Vec3 c = pose.center();

//          out_c << c(0) << " " << c(1)<< " " << c(2) << std::endl;

//      }
//      out_c.close();

  }

  {
//      std::ofstream out_show_camera;
//      out_show_camera.open("/media/add7/E/run_scene/testMerge/29/out_camera_res_R.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 c = pose.center();

//          Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
//          axis_x = R.transpose() * (axis_x + R*c);
//          axis_y = R.transpose() * (axis_y + R*c);
//          axis_z = R.transpose() * (axis_z + R*c);
//          out_show_camera <<R(0,0)<<" "<< R(0,1)<< " "<<R(0,2)<< " "
//                          <<R(1,0)<<" "<< R(1,1)<< " "<<R(1,2)<< " "
//                          <<R(2,0)<<" "<< R(2,1)<< " "<<R(2,2)<< std::endl;
//      }
//      out_show_camera.close();
  }





  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }





  return true;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Refine_initValue()
{
    //
    // Build the Relative pose graph from matches:
    //
    /// pairwise view relation between poseIds
    using PoseWiseMatches = std::map< Pair, Pair_Set >;

    // List shared correspondences (pairs) between poses
    PoseWiseMatches poseWiseMatches;
    for (const auto & iterMatches : matches_provider_->pairWise_matches_)
    {
      const Pair pair = iterMatches.first;
      const View * v1 = sfm_data_.GetViews().at(pair.first).get();
      const View * v2 = sfm_data_.GetViews().at(pair.second).get();
      poseWiseMatches[Pair(v1->id_pose, v2->id_pose)].insert(pair);
    }

    C_Progress_display my_progress_bar( poseWiseMatches.size(),
        std::cout, "\n- Refine init value -\n" );

  #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
  #endif
    // Compute the relative pose from pairwise point matches:
    for (int i = 0; i < static_cast<int>(poseWiseMatches.size()); ++i)
    {
      ++my_progress_bar;
      {
        PoseWiseMatches::const_iterator iter (poseWiseMatches.begin());
        std::advance(iter, i);
        const auto & relative_pose_iterator(*iter);
        const Pair relative_pose_pair = relative_pose_iterator.first;
        const Pair_Set & match_pairs = relative_pose_iterator.second;

        // If a pair has the same ID, discard it
        if (relative_pose_pair.first == relative_pose_pair.second)
        {
          continue;
        }

        // Select common bearing vectors
        if (match_pairs.size() > 1)
        {
          std::cerr << "Compute relative pose between more than two view is not supported" << std::endl;
          continue;
        }

        const Pair pairIterator = *(match_pairs.begin());

        const IndexT I = pairIterator.first;
        const IndexT J = pairIterator.second;

        const View * view_I = sfm_data_.views[I].get();
        const View * view_J = sfm_data_.views[J].get();

        // Check that valid cameras are existing for the pair of view
        if (sfm_data_.GetIntrinsics().count(view_I->id_intrinsic) == 0 ||
          sfm_data_.GetIntrinsics().count(view_J->id_intrinsic) == 0)
          continue;


        const IntrinsicBase * cam_I = sfm_data_.GetIntrinsics().at(view_I->id_intrinsic).get();
        const IntrinsicBase * cam_J = sfm_data_.GetIntrinsics().at(view_J->id_intrinsic).get();

        // Setup corresponding bearing vector
        const matching::IndMatches & matches = matches_provider_->pairWise_matches_.at(pairIterator);
        size_t nbBearing = matches.size();
        Mat x1(2, nbBearing), x2(2, nbBearing);
        nbBearing = 0;
        for (const auto & match : matches)
        {
          x1.col(nbBearing) = ((*cam_I)(cam_I->get_ud_pixel(features_provider_->feats_per_view[I][match.i_].coords().cast<double>()))).hnormalized();
          x2.col(nbBearing++) = ((*cam_J)(cam_J->get_ud_pixel(features_provider_->feats_per_view[J][match.j_].coords().cast<double>()))).hnormalized();
        }

        RelativePose_Info relativePose_info;
        // Compute max authorized error as geometric mean of camera plane tolerated residual error
        relativePose_info.initial_residual_tolerance = std::pow(
          cam_I->imagePlane_toCameraPlaneError(2.5) *
          cam_J->imagePlane_toCameraPlaneError(2.5),
          1./2.);

        // Since we use normalized features, we will use unit image size and intrinsic matrix:
        const std::pair<size_t, size_t> imageSize(1., 1.);
        const Mat3 K  = Mat3::Identity();

//        if (!robustRelativePose(K, K, x1, x2, relativePose_info, imageSize, imageSize, 256))
//        {
//          continue;
//        }
        const bool bRefine_using_BA = true;
        if (bRefine_using_BA)
        {
          // Refine the defined scene
          SfM_Data tiny_scene;
          tiny_scene.views.insert(*sfm_data_.GetViews().find(view_I->id_view));
          tiny_scene.views.insert(*sfm_data_.GetViews().find(view_J->id_view));
          tiny_scene.intrinsics.insert(*sfm_data_.GetIntrinsics().find(view_I->id_intrinsic));
          tiny_scene.intrinsics.insert(*sfm_data_.GetIntrinsics().find(view_J->id_intrinsic));

          // Init poses
          const Pose3 & Pose_I = tiny_scene.poses[view_I->id_pose] = Pose3(Mat3::Identity(), Vec3::Zero());
          const Pose3 & Pose_J = tiny_scene.poses[view_J->id_pose] = relativePose_info.relativePose;

          // Init structure
          const Mat34 P1 = cam_I->get_projective_equivalent(Pose_I);
          const Mat34 P2 = cam_J->get_projective_equivalent(Pose_J);
          Landmarks & landmarks = tiny_scene.structure;
          for (Mat::Index k = 0; k < x1.cols(); ++k)
          {
            const Vec2 x1_ = features_provider_->feats_per_view[I][matches[k].i_].coords().cast<double>();
            const Vec2 x2_ = features_provider_->feats_per_view[J][matches[k].j_].coords().cast<double>();
            Vec3 X;
            TriangulateDLT(P1, x1_.homogeneous(), P2, x2_.homogeneous(), &X);
            Observations obs;
            obs[view_I->id_view] = Observation(x1_, matches[k].i_);
            obs[view_J->id_view] = Observation(x2_, matches[k].j_);
            landmarks[k].obs = obs;
            landmarks[k].X = X;
          }
          // - refine only Structure and Rotations & translations (keep intrinsic constant)
          Bundle_Adjustment_Ceres::BA_Ceres_options options(false, false);
          options.linear_solver_type_ = ceres::DENSE_SCHUR;
          Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
          const Optimize_Options ba_refine_options
            (Intrinsic_Parameter_Type::NONE, // -> Keep intrinsic constant
            Extrinsic_Parameter_Type::ADJUST_ALL, // adjust camera motion
            Structure_Parameter_Type::ADJUST_ALL);// adjust scene structure
          if (bundle_adjustment_obj.Adjust(tiny_scene, ba_refine_options))
          {
              std::cout << "this ba" << std::endl;
            // --> to debug: save relative pair geometry on disk
            // std::ostringstream os;
            // os << relative_pose_pair.first << "_" << relative_pose_pair.second << ".ply";
            // Save(tiny_scene, os.str(), ESfM_Data(STRUCTURE | EXTRINSICS));
            //
            const Mat3 R1 = tiny_scene.poses[view_I->id_pose].rotation();
            const Mat3 R2 = tiny_scene.poses[view_J->id_pose].rotation();
            const Vec3 t1 = tiny_scene.poses[view_I->id_pose].translation();
            const Vec3 t2 = tiny_scene.poses[view_J->id_pose].translation();
            // Compute relative motion and save it
            Mat3 Rrel;
            Vec3 trel;
            RelativeCameraMotion(R1, t1, R2, t2, &Rrel, &trel);
            // Update found relative pose
            relativePose_info.relativePose = Pose3(Rrel, -Rrel.transpose() * trel);
          }
        }
  #ifdef OPENMVG_USE_OPENMP
        #pragma omp critical
  #endif
        {
//          // Add the relative rotation to the relative 'rotation' pose graph
//          using namespace openMVG::rotation_averaging;
//            vec_relatives_R.emplace_back(
//              relative_pose_pair.first, relative_pose_pair.second,
//              relativePose_info.relativePose.rotation(),
//              1.f);
        }
      }
    } // for all relative pose

    // Log input graph to the HTML report
    if (!sLogging_file_.empty() && !sOut_directory_.empty())
    {
      // Log a relative view graph
      {
//        std::set<IndexT> set_ViewIds;
//        std::transform(sfm_data_.GetViews().begin(), sfm_data_.GetViews().end(),
//          std::inserter(set_ViewIds, set_ViewIds.begin()), stl::RetrieveKey());
//        graph::indexedGraph putativeGraph(set_ViewIds, getPairs(matches_provider_->pairWise_matches_));
//        graph::exportToGraphvizData(
//          stlplus::create_filespec(sOut_directory_, "global_relative_rotation_view_graph"),
//          putativeGraph);
      }

      // Log a relative pose graph
      {
//        std::set<IndexT> set_pose_ids;
//        Pair_Set relative_pose_pairs;
//        for (const auto & relative_R : vec_relatives_R)
//        {
//          const Pair relative_pose_indices(relative_R.i, relative_R.j);
//          relative_pose_pairs.insert(relative_pose_indices);
//          set_pose_ids.insert(relative_R.i);
//          set_pose_ids.insert(relative_R.j);
//        }
//        const std::string sGraph_name = "global_relative_rotation_pose_graph";
//        graph::indexedGraph putativeGraph(set_pose_ids, relative_pose_pairs);
//        graph::exportToGraphvizData(
//          stlplus::create_filespec(sOut_directory_, sGraph_name),
//          putativeGraph);
//        using namespace htmlDocument;
//        std::ostringstream os;

//        os << "<br>" << "global_relative_rotation_pose_graph" << "<br>"
//           << "<img src=\""
//           << stlplus::create_filespec(sOut_directory_, "global_relative_rotation_pose_graph", "svg")
//           << "\" height=\"600\">\n";

//        html_doc_stream_->pushInfo(os.str());
      }
    }



}


struct PoseInfo
{
    std::vector<double> C_gps;
    std::vector<double> R_imu;
    int status;
};
///Hash_Map<IndexT, PoseInfo> poseInfoList; ///use this

bool GlobalSfMReconstructionEngine_RelativeMotions::Process_threeCameras_imu_gps_c4_changeFileForm(std::string sMatchesDir) {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if(set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

//    openMVG::rotation_averaging::RelativeRotations relatives_R;
//    Compute_Relative_Rotations(relatives_R);

    Hash_Map<IndexT, Mat3> global_rotations;

    //test distortion
    Mat3 r_tb, r_tf;
    r_tb << 1.0, 0.0, 0.0,
            0.0, 0.787926, -0.61577,
            0.0, 0.61577, 0.787926;
    r_tf << 1.0, 0.0, 0.0,
            0.0, 0.78796,  0.615727,
            0.0, -0.615727, 0.78796;
    Vec3 RA_tb, RA_tf;
    ceres::RotationMatrixToAngleAxis(r_tb.data(), RA_tb.data());
    ceres::RotationMatrixToAngleAxis(r_tf.data(), RA_tf.data());

        //set reference
        double R_br_angleAxis[3]= {RA_tb(0), RA_tb(1), -RA_tb(2)};//{0.665224, 0.00448363, 0.00291885};  //{0.703355, 0.00263495, -0.00149232};//{0.0, 0.0, 0.0}; //{0.0212232, 0.0206259, 0.0420509};  //{0.0, 0.0, 0.0}; //{0.669661, 0.00354668, 0.00220508}; //  //{0.0, 0.0, 0.0}; //{0.662013, -0.00281273, 0.000638113 }; //{0.0, 0.0, 0.0}; //{0.671924, 0.00432461, 0.000528427};
        double R_fr_angleAxis[3]= {RA_tf(0), RA_tf(1), -RA_tf(2)};//{-0.644463, 0.0113345, -0.00156344}; //{-0.682501, 0.0147873, 0.00150509};//{0.0, 0.0, 0.0}; //{-0.0232707, -0.00288052, 0.00035795};  //{0.0, 0.0, 0.0}; //{-0.650191, 0.0101123, -0.000844533}; //  //{0.0, 0.0, 0.0}; //{-0.658446, 0.00532591, 0.0010695}; //{-0.656461, 0.00508743, -0.00299628};
        double R_dr_angleAxis[3]= {0.0, 0.0, 0.0};
        Vec3 t_db, t_df, t_dd;
        //double S = 1470.531783220221;
        t_db << 0.0, 0.0, 0.0;//0.0, 0.42862, -0.260222;//-1.14174, 0.37576, -1.73691;  //-0.640808, 15.8409, 0.7518;//0.0, 0.0, 0.0; //2.3162, -4.69223, 1.55132;  //0.0, 0.0, 0.0; //-0.722651, 1.72471, -0.664053; // //0.0, 0.0, 0.0; //0.000207573, -5.92667e-05, -0.00025752; //0.0, 0.0, 0.0; //0.422496, 3.39275, 2.62513;  //S*0.000948552, S*0.000923492, S*-0.00271615;
        t_df << 0.0, 0.0, 0.0;//0.0, -0.181266, 0.457972;//-0.457481, 5.35573, -3.39704;  //-2.75629, -10.0797, -0.169824;//0.0, 0.0, 0.0; //-0.709151, 4.71803, 2.80816;  //0.0, 0.0, 0.0; //-0.123532, 3.24503, -3.61083; // //0.0, 0.0, 0.0; //0.000217762, 7.85402e-05, -0.000246061; //0.0, 0.0, 0.0; //-0.303247, 0.647175, -2.00489;  //S*0.000209318, S*-0.0018472, S*0.00271677;
        t_dd << 0.0, 0.0, 0.0;


        std::vector<double> transformation_br(6), transformation_fr(6), transformation_dr(6);
        transformation_br = {R_br_angleAxis[0], R_br_angleAxis[1], R_br_angleAxis[2],
                             t_db(0), t_db(1), t_db(2)};
        transformation_fr = {R_fr_angleAxis[0], R_fr_angleAxis[1], R_fr_angleAxis[2],
                             t_df(0), t_df(1), t_df(2)};
        transformation_dr = {R_dr_angleAxis[0], R_dr_angleAxis[1], R_dr_angleAxis[2],
                             t_dd(0), t_dd(1), t_dd(2)};


    const std::vector<double> transformation_R_imu = {0.0, 0.0, 0.0}; //{-0.000427019, 0.0125482, 0.0134619};  //{-0.00986379, 0.00976848, 0.0253123};
    const std::vector<double> tramsformation_x_gps = {0.0, 0.0, 0.0};  //{-0.340012, 0.33451, 3.08947};


        std::cout << "tramsformation_x_gps : " << tramsformation_x_gps[0] << " "
                                               << tramsformation_x_gps[1] << " "
                                               << tramsformation_x_gps[2] << std::endl;


        //get gps and imu data
      {
    //        /Hash_Map<IndexT, std::vector<double> > R_imu;
            //Hash_Map<IndexT, std::vector<double> > C_gps;
            Hash_Map<IndexT, PoseInfo> poseInfoList; ///use this
            {
                //read file and initialise external
                //std::vector<double> c(3);
                std::string tFilePath = sfm_data_.s_root_path + "/n_gps_imu_321.res";  //!!!
                std::ifstream in;
                in.open(tFilePath);
                if(!in)
                {
                    std::cout << "open tFilePath file error!" << std::endl;
                    return false;
                }
                std::size_t line = 0;
                int poseCount;
                in >> poseCount;
                while((!in.eof()) && (line < poseCount))
                {
                    int poseId;
                    double x,y,z;
                    double roll, pitch, yaw;
                    int status;
                    in >> poseId >> x >> y >> z >> pitch >> roll >> yaw >> status;

                    poseInfoList[poseCount].C_gps[0] = x;
                    poseInfoList[poseCount].C_gps[1] = y;
                    poseInfoList[poseCount].C_gps[2] = z;

                    Mat3 Rz, Ry, Rx;
                    Rz << cos(yaw/180.0*M_PI), -sin(yaw/180.0*M_PI), 0.0,
                          sin(yaw/180.0*M_PI), cos(yaw/180.0*M_PI), 0.0,
                          0.0, 0.0, 1.0;
                    Ry << cos(roll/180.0*M_PI), 0.0, sin(roll/180.0*M_PI),
                          0.0, 1.0, 0.0,
                          -sin(roll/180.0*M_PI), 0.0, cos(roll/180.0*M_PI);
                    Rx << 1.0, 0.0, 0.0,
                          0.0, cos(pitch/180.0*M_PI), -sin(pitch/180.0*M_PI),
                          0.0, sin(pitch/180.0*M_PI), cos(pitch/180.0*M_PI);
                    Mat3 R_imu =  Rx*Ry*Rz;
                    Mat3 changeAxis_M;
                    changeAxis_M << 1.0, 0.0, 0.0,
                                    0.0, -1.0, 0.0,
                                    0.0, 0.0, -1.0;
                    R_imu = changeAxis_M * R_imu;
                    /////////////////////////////////////////////////


                    //set
                    Mat3 R_imu_trans;
                    ceres::AngleAxisToRotationMatrix(&transformation_R_imu[0], R_imu_trans.data());
                    Vec3 gps_trans = Vec3(tramsformation_x_gps[0], tramsformation_x_gps[1], tramsformation_x_gps[2]);
                    Vec3 C_gps_res = Vec3(x, y, z);
                    if(poseId < 199999) //back
                    {
                        Mat3 Rb_trans;
                        Vec3 tb_trans;
                        ceres::AngleAxisToRotationMatrix(&transformation_br[0], Rb_trans.data());
                        tb_trans = Vec3(transformation_br[3], transformation_br[4], transformation_br[5]);

                        //set b
                        Mat3 Rb_res = Rb_trans * R_imu_trans * R_imu;
                        Vec3 tb_res;
                        tb_res = Rb_trans * gps_trans - Rb_trans * R_imu_trans * R_imu * C_gps_res + tb_trans;
                        global_rotations[poseId] = Rb_res;
                        sfm_data_.poses[poseId] = Pose3(Rb_res, -Rb_res.transpose() * tb_res);

                    }else if(poseId < 299999) //down
                    {
                        Mat3 Rd_trans;
                        Vec3 td_trans;
                        ceres::AngleAxisToRotationMatrix(&transformation_dr[0], Rd_trans.data());
                        td_trans = Vec3(transformation_dr[3], transformation_dr[4], transformation_dr[5]);

                        //set d
                        Mat3 Rd_res = Rd_trans * R_imu_trans * R_imu;
                        Vec3 td_res;
                        td_res = Rd_trans * gps_trans - Rd_trans * R_imu_trans * R_imu * C_gps_res + td_trans;
                        global_rotations[poseId] = Rd_res;
                        sfm_data_.poses[poseId] = Pose3(Rd_res, -Rd_res.transpose() * td_res);


                    }else if(poseId < 399999) //front
                    {
                        Mat3 Rf_trans;
                        Vec3 tf_trans;
                        ceres::AngleAxisToRotationMatrix(&transformation_fr[0], Rf_trans.data());
                        tf_trans = Vec3(transformation_fr[3], transformation_fr[4], transformation_fr[5]);

                        //set f
                        Mat3 Rf_res = Rf_trans * R_imu_trans * R_imu;
                        Vec3 tf_res;
                        tf_res = Rf_trans * gps_trans - Rf_trans * R_imu_trans * R_imu * C_gps_res + tf_trans;
                        global_rotations[poseId] = Rf_res;
                        sfm_data_.poses[poseId] = Pose3(Rf_res, -Rf_res.transpose() * tf_res);

                    }

                    ++line;
                }
                in.close();

            }
      }


  matching::PairWiseMatches tripletWise_matches;
  if (!Compute_Global_Translations_gps(global_rotations, tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
    return false;
  }

    if (!Compute_Initial_Structure(tripletWise_matches))
    {
      std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
      return false;
    }

//    for( int i = 0; i < sfm_data_.poses.size(); ++i)
//    {
//        double R0_angleAxis[3];
//        Mat3 mid_R0 = sfm_data_.poses[i].rotation();
//        ceres::RotationMatrixToAngleAxis((const double*)mid_R0.data(), R0_angleAxis);

//        std::cout << "R0_angleAxis : " << R0_angleAxis[0] << " " << R0_angleAxis[1] << " " << R0_angleAxis[2] << std::endl;

//    }

    {
//            std::cout<< "initial intrinsics :"<<std::endl;
//                Hash_Map<IndexT, std::vector<double> > map_intrinsics;
//                map_intrinsics[0] = sfm_data_.intrinsics[0]->getParams();
//                std::cout<< "map_intrinsics 0 " << map_intrinsics[0][0]<<" "<<map_intrinsics[0][1]<<" "<<
//                         map_intrinsics[0][2]<<" "<<map_intrinsics[0][3]<<" "<<map_intrinsics[0][4]<<" "<<
//                         map_intrinsics[0][5]<<" "<<map_intrinsics[0][6]<<" "<<map_intrinsics[0][7]<<std::endl;

//                map_intrinsics[1] = sfm_data_.intrinsics[1]->getParams();
//                std::cout<< "map_intrinsics 1 " << map_intrinsics[1][0]<<" "<<map_intrinsics[1][1]<<" "<<
//                         map_intrinsics[1][2]<<" "<<map_intrinsics[1][3]<<" "<<map_intrinsics[1][4]<<" "<<
//                         map_intrinsics[1][5]<<" "<<map_intrinsics[1][6]<<" "<<map_intrinsics[1][7]<<std::endl;

//                map_intrinsics[2] = sfm_data_.intrinsics[2]->getParams();
//                std::cout<< "map_intrinsics 2 " << map_intrinsics[2][0]<<" "<<map_intrinsics[2][1]<<" "<<
//                         map_intrinsics[2][2]<<" "<<map_intrinsics[2][3]<<" "<<map_intrinsics[2][4]<<" "<<
//                         map_intrinsics[2][5]<<" "<<map_intrinsics[2][6]<<" "<<map_intrinsics[2][7]<<std::endl;
//                std::cout <<"-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() <<std::endl;
    }


  if (!Adjust_threecamera_gps_change_c4())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }




  {
//          std::cout<< "map_intrinsics 1"<<std::endl;
//              Hash_Map<IndexT, std::vector<double> > map_intrinsics;
//              map_intrinsics[0] = sfm_data_.intrinsics[0]->getParams();
//              std::cout<< "map_intrinsics 0 " << map_intrinsics[0][0]<<" "<<map_intrinsics[0][1]<<" "<<
//                       map_intrinsics[0][2]<<" "<<map_intrinsics[0][3]<<" "<<map_intrinsics[0][4]<<" "<<
//                       map_intrinsics[0][5]<<" "<<map_intrinsics[0][6]<<" "<<map_intrinsics[0][7]<<std::endl;

//              map_intrinsics[1] = sfm_data_.intrinsics[1]->getParams();
//              std::cout<< "map_intrinsics 1 " << map_intrinsics[1][0]<<" "<<map_intrinsics[1][1]<<" "<<
//                       map_intrinsics[1][2]<<" "<<map_intrinsics[1][3]<<" "<<map_intrinsics[1][4]<<" "<<
//                       map_intrinsics[1][5]<<" "<<map_intrinsics[1][6]<<" "<<map_intrinsics[1][7]<<std::endl;

//              map_intrinsics[2] = sfm_data_.intrinsics[2]->getParams();
//              std::cout<< "map_intrinsics 2 " << map_intrinsics[2][0]<<" "<<map_intrinsics[2][1]<<" "<<
//                       map_intrinsics[2][2]<<" "<<map_intrinsics[2][3]<<" "<<map_intrinsics[2][4]<<" "<<
//                       map_intrinsics[2][5]<<" "<<map_intrinsics[2][6]<<" "<<map_intrinsics[2][7]<<std::endl;
//              std::cout <<"-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() <<std::endl;
  }


  //out
  {
//      std::ofstream out_show_camera;
//      out_show_camera.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_reference.txt");
//      std::ofstream out_Rc;
//      out_Rc.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_Rc.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 c = pose.center();

//          Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
//          axis_x = R.transpose() * (axis_x + R*c);
//          axis_y = R.transpose() * (axis_y + R*c);
//          axis_z = R.transpose() * (axis_z + R*c);
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_x(0) << " " << axis_x(1) << " " << axis_x(2)<< ";" << std::endl;
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_y(0) << " " << axis_y(1) << " " << axis_y(2)<< ";" << std::endl;
//          out_show_camera <<c(0)<<" "<< c(1)<< " "<<c(2)<< " "
//                 << axis_z(0) << " " << axis_z(1) << " " << axis_z(2)<< ";" << std::endl;

//          out_Rc << c(0) << " " << c(1)<< " " << c(2) << " " << std::endl;

//      }
//      out_Rc.close();
//      out_show_camera.close();

  }

  //out
  {
//      std::ofstream out_c;
//      out_c.open("/media/add7/E/run_scene/testMerge/29/out_show_camera_c.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Vec3 c = pose.center();

//          out_c << c(0) << " " << c(1)<< " " << c(2) << std::endl;

//      }
//      out_c.close();

  }

  {
//      std::ofstream out_show_camera;
//      out_show_camera.open("/media/add7/E/run_scene/testMerge/29/out_camera_res_R.txt");
//      for(Poses::const_iterator it = sfm_data_.poses.begin();
//          it != sfm_data_.poses.end(); ++it)
//      {
//          const Pose3 & pose = it->second;
//          const Mat3 R = pose.rotation();
//          const Vec3 c = pose.center();

//          Vec3 axis_x={1.0,0.0,0.0}, axis_y={0.0,1.0,0.0}, axis_z={0.0,0.0,1.0};
//          axis_x = R.transpose() * (axis_x + R*c);
//          axis_y = R.transpose() * (axis_y + R*c);
//          axis_z = R.transpose() * (axis_z + R*c);
//          out_show_camera <<R(0,0)<<" "<< R(0,1)<< " "<<R(0,2)<< " "
//                          <<R(1,0)<<" "<< R(1,1)<< " "<<R(1,2)<< " "
//                          <<R(2,0)<<" "<< R(2,1)<< " "<<R(2,2)<< std::endl;
//      }
//      out_show_camera.close();
  }





  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }





  return true;
}



/// Compute from relative rotations the global rotations of the camera poses
bool GlobalSfMReconstructionEngine_RelativeMotions::Compute_Global_Rotations
(
  const rotation_averaging::RelativeRotations & relatives_R,
  Hash_Map<IndexT, Mat3> & global_rotations
)
{
  if (relatives_R.empty())
    return false;
  // Log statistics about the relative rotation graph
  {
    std::set<IndexT> set_pose_ids;
    for (const auto & relative_R : relatives_R)
    {
      set_pose_ids.insert(relative_R.i);
      set_pose_ids.insert(relative_R.j);
    }

    std::cout << "\n-------------------------------" << "\n"
      << " Global rotations computation: " << "\n"
      << "  #relative rotations: " << relatives_R.size() << "\n"
      << "  #global rotations: " << set_pose_ids.size() << std::endl;
  }

  // Global Rotation solver:
  const ERelativeRotationInferenceMethod eRelativeRotationInferenceMethod =
    TRIPLET_ROTATION_INFERENCE_COMPOSITION_ERROR;
    //TRIPLET_ROTATION_INFERENCE_NONE;

  system::Timer t;
  GlobalSfM_Rotation_AveragingSolver rotation_averaging_solver;
  const bool b_rotation_averaging = rotation_averaging_solver.Run(
    eRotation_averaging_method_, eRelativeRotationInferenceMethod,
    relatives_R, global_rotations);

  std::cout
    << "Found #global_rotations: " << global_rotations.size() << "\n"
    << "Timing: " << t.elapsed() << " seconds" << std::endl;


  if (b_rotation_averaging)
  {
    // Compute & display rotation fitting residual errors
    std::vector<float> vec_rotation_fitting_error;
    vec_rotation_fitting_error.reserve(relatives_R.size());
    for (const auto & relative_R : relatives_R)
    {
      const Mat3 & Rij = relative_R.Rij;
      const IndexT i = relative_R.i;
      const IndexT j = relative_R.j;
      if (global_rotations.count(i)==0 || global_rotations.count(j)==0)
        continue;
      const Mat3 & Ri = global_rotations[i];
      const Mat3 & Rj = global_rotations[j];
      const Mat3 eRij(Rj.transpose()*Rij*Ri);
      const double angularErrorDegree = R2D(getRotationMagnitude(eRij));
      vec_rotation_fitting_error.push_back(angularErrorDegree);
    }

    if (!vec_rotation_fitting_error.empty())
    {
      const float error_max = *max_element(vec_rotation_fitting_error.begin(), vec_rotation_fitting_error.end());
      Histogram<float> histo(0.0f,error_max, 20);
      histo.Add(vec_rotation_fitting_error.begin(), vec_rotation_fitting_error.end());
      std::cout
        << "\nRelative/Global degree rotations residual errors {0," << error_max<< "}:"
        << histo.ToString() << std::endl;
      {
        Histogram<float> histo(0.0f, 5.0f, 20);
        histo.Add(vec_rotation_fitting_error.begin(), vec_rotation_fitting_error.end());
        std::cout
          << "\nRelative/Global degree rotations residual errors {0,5}:"
          << histo.ToString() << std::endl;
      }
      std::cout << "\nStatistics about global rotation evaluation:" << std::endl;
      minMaxMeanMedian<float>(vec_rotation_fitting_error.begin(), vec_rotation_fitting_error.end());
    }

    // Log input graph to the HTML report
    if (!sLogging_file_.empty() && !sOut_directory_.empty())
    {
      // Log a relative pose graph
      {
        std::set<IndexT> set_pose_ids;
        Pair_Set relative_pose_pairs;
        for (const auto & view : sfm_data_.GetViews())
        {
          const IndexT pose_id = view.second->id_pose;
          set_pose_ids.insert(pose_id);
        }
        const std::string sGraph_name = "global_relative_rotation_pose_graph_final";
        graph::indexedGraph putativeGraph(set_pose_ids, rotation_averaging_solver.GetUsedPairs());
//        graph::exportToGraphvizData(
//          stlplus::create_filespec(sOut_directory_, sGraph_name),
//          putativeGraph);

        using namespace htmlDocument;
        std::ostringstream os;

        os << "<br>" << sGraph_name << "<br>"
           << "<img src=\""
           << stlplus::create_filespec(sOut_directory_, sGraph_name, "svg")
           << "\" height=\"600\">\n";

        html_doc_stream_->pushInfo(os.str());
      }
    }
  }
  return b_rotation_averaging;
}

/// Compute/refine relative translations and compute global translations
bool GlobalSfMReconstructionEngine_RelativeMotions::Compute_Global_Translations
(
  const Hash_Map<IndexT, Mat3> & global_rotations,
  matching::PairWiseMatches & tripletWise_matches
)
{
  // Translation averaging (compute translations & update them to a global common coordinates system)
  GlobalSfM_Translation_AveragingSolver translation_averaging_solver;
  const bool bTranslationAveraging = translation_averaging_solver.Run(
    eTranslation_averaging_method_,
    sfm_data_,
    features_provider_,
    matches_provider_,
    global_rotations,
    tripletWise_matches);

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "cameraPath_translation_averaging", "ply"),
      ESfM_Data(EXTRINSICS));
  }

  return bTranslationAveraging;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Compute_Global_Translations_gps
(
  const Hash_Map<IndexT, Mat3> & global_rotations,
  matching::PairWiseMatches & tripletWise_matches
)
{
  // Translation averaging (compute translations & update them to a global common coordinates system)
  GlobalSfM_Translation_AveragingSolver translation_averaging_solver;
  const bool bTranslationAveraging = translation_averaging_solver.Run_gps(
    eTranslation_averaging_method_,
    sfm_data_,
    features_provider_,
    matches_provider_,
    global_rotations,
    tripletWise_matches);

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "cameraPath_translation_averaging", "ply"),
      ESfM_Data(EXTRINSICS));
  }

  return bTranslationAveraging;
}

//#define USE_ALL_VALID_MATCHES
/// Compute the initial structure of the scene
bool GlobalSfMReconstructionEngine_RelativeMotions::Compute_Initial_Structure
(
  matching::PairWiseMatches & tripletWise_matches
)
{
  // Build tracks from selected triplets (Union of all the validated triplet tracks (_tripletWise_matches))
  {
    using namespace openMVG::tracks;
    TracksBuilder tracksBuilder;
//#if defined USE_ALL_VALID_MATCHES // not used by default
    matching::PairWiseMatches pose_supported_matches;

    //for (const std::pair< Pair, matching::IndMatches > & match_info :  matches_provider_->pairWise_matches_)
    for (const auto & match_info :  matches_provider_->pairWise_matches_)
    {
      const View * vI = sfm_data_.GetViews().at(match_info.first.first).get();
      const View * vJ = sfm_data_.GetViews().at(match_info.first.second).get();
      if (sfm_data_.IsPoseAndIntrinsicDefined(vI) && sfm_data_.IsPoseAndIntrinsicDefined(vJ))
      {
        pose_supported_matches.insert(match_info);
      }
    }
    tracksBuilder.Build(pose_supported_matches);

//#else
    // Use triplet validated matches
//    tracksBuilder.Build(tripletWise_matches);
//#endif
    tracksBuilder.Filter(3);
    STLMAPTracks map_selectedTracks; // reconstructed track (visibility per 3D point)
    tracksBuilder.ExportToSTL(map_selectedTracks);

    // Fill sfm_data with the computed tracks (no 3D yet)
    Landmarks & structure = sfm_data_.structure;
    IndexT idx(0);
    for (STLMAPTracks::const_iterator itTracks = map_selectedTracks.begin();
      itTracks != map_selectedTracks.end();
      ++itTracks, ++idx)
    {
      const submapTrack & track = itTracks->second;
      structure[idx] = Landmark();
      Observations & obs = structure.at(idx).obs;
      for (submapTrack::const_iterator it = track.begin(); it != track.end(); ++it)
      {
        const size_t imaIndex = it->first;
        const size_t featIndex = it->second;
        const PointFeature & pt = features_provider_->feats_per_view.at(imaIndex)[featIndex];
        obs[imaIndex] = Observation(pt.coords().cast<double>(), featIndex);
      }
    }

    std::cout << std::endl << "Track stats" << std::endl;
    {
      std::ostringstream osTrack;
      //-- Display stats:
      //    - number of images
      //    - number of tracks
      std::set<uint32_t> set_imagesId;
      TracksUtilsMap::ImageIdInTracks(map_selectedTracks, set_imagesId);
      osTrack << "------------------" << "\n"
        << "-- Tracks Stats --" << "\n"
        << " Tracks number: " << tracksBuilder.NbTracks() << "\n"
        << " Images Id: " << "\n";
      std::copy(set_imagesId.begin(),
        set_imagesId.end(),
        std::ostream_iterator<uint32_t>(osTrack, ", "));
      osTrack << "\n------------------" << "\n";

      std::map<uint32_t, uint32_t> map_Occurence_TrackLength;
      TracksUtilsMap::TracksLength(map_selectedTracks, map_Occurence_TrackLength);
      osTrack << "TrackLength, Occurrence" << "\n";
      for (const auto & iter : map_Occurence_TrackLength)  {
        osTrack << "\t" << iter.first << "\t" << iter.second << "\n";
      }
      osTrack << "\n";
      std::cout << osTrack.str();
    }
  }

  // Compute 3D position of the landmark of the structure by triangulation of the observations
  {
    openMVG::system::Timer timer;

    const IndexT trackCountBefore = sfm_data_.GetLandmarks().size();
    SfM_Data_Structure_Computation_Blind structure_estimator(true);
    structure_estimator.triangulate(sfm_data_);

    std::cout << "\n#removed tracks (invalid triangulation): " <<
      trackCountBefore - IndexT(sfm_data_.GetLandmarks().size()) << std::endl;
    std::cout << std::endl << "  Triangulation took (s): " << timer.elapsed() << std::endl;

    // Export initial structure
    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "initial_structure", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }
  return !sfm_data_.structure.empty();
}

// Adjust the scene (& remove outliers)
bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

    std::cout << "this->b_use_motion_prior_ : " << this->b_use_motion_prior_ << std::endl;

    std::cout << "with pose id: " << std::endl;
    for (const auto & pose_it : sfm_data_.poses)
    {
      const IndexT indexPose = pose_it.first;
      std::cout << indexPose << " ";
    }
    std::cout << std::endl;

    Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = bundle_adjustment_obj.Adjust
    (
      sfm_data_,
      Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::NONE,//ADJUST_TRANSLATION, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
        Structure_Parameter_Type::ADJUST_ALL,
        Control_Point_Parameter(),
        this->b_use_motion_prior_)
    );
  if (b_BA_Status)
  {
    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }
  // Check that poses & intrinsic cover some measures (after outlier removal)
  const IndexT minPointPerPose = 12; // 6 min
  const IndexT minTrackLength = 3; // 2 min
  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
  {
    // TODO: must ensure that track graph is producing a single connected component

    const size_t pointcount_cleaning = sfm_data_.structure.size();
    std::cout << "Point_cloud cleaning:\n"
      << "\t #3DPoints: " << pointcount_cleaning << "\n";
  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  const Optimize_Options ba_refine_options(
    ReconstructionEngine::intrinsic_refinement_options_,
    Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
    Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
    Control_Point_Parameter(),
    this->b_use_motion_prior_);

  b_BA_Status = bundle_adjustment_obj.Adjust(sfm_data_, ba_refine_options);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  std::cout << "with pose id: " << std::endl;
  for (const auto & pose_it : sfm_data_.poses)
  {
    const IndexT indexPose = pose_it.first;
    std::cout << indexPose << " ";
  }
  std::cout << std::endl;


  return b_BA_Status;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_nonremove()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

    std::cout << "this->b_use_motion_prior_ : " << this->b_use_motion_prior_ << std::endl;

    std::cout << "with pose id: " << std::endl;
    for (const auto & pose_it : sfm_data_.poses)
    {
      const IndexT indexPose = pose_it.first;
      std::cout << indexPose << " ";
    }
    std::cout << std::endl;

    Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status;
  b_BA_Status= bundle_adjustment_obj.Adjust
    (
      sfm_data_,
      Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_ALL,//ADJUST_TRANSLATION, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
        Structure_Parameter_Type::ADJUST_ALL,
        Control_Point_Parameter(),
        this->b_use_motion_prior_)
    );
//  if (b_BA_Status)
//  {
//    if (!sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

//    // - refine only Structure and Rotations & translations
//    b_BA_Status = bundle_adjustment_obj.Adjust
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    if (b_BA_Status && !sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }
//  }

//  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
//    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
//    b_BA_Status = bundle_adjustment_obj.Adjust
//      (
//        sfm_data_,
//        Optimize_Options(
//          ReconstructionEngine::intrinsic_refinement_options_,
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    if (b_BA_Status && !sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }
//  }

//  // Remove outliers (max_angle, residual error)
//  const size_t pointcount_initial = sfm_data_.structure.size();
//  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
//  RemoveOutliers_AngleError(sfm_data_, 2.0);
//  const size_t pointcount_angular_filter = sfm_data_.structure.size();
//  std::cout << "Outlier removal (remaining #points):\n"
//    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
//    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
//    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

//  if (!sLogging_file_.empty())
//  {
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
//  }
//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
//  const Optimize_Options ba_refine_options(
//    ReconstructionEngine::intrinsic_refinement_options_,
//    Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
//    Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
//    Control_Point_Parameter(),
//    this->b_use_motion_prior_);

//  b_BA_Status = bundle_adjustment_obj.Adjust(sfm_data_, ba_refine_options);
//  if (b_BA_Status && !sLogging_file_.empty())
//  {
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
//  }

//  std::cout << "with pose id: " << std::endl;
//  for (const auto & pose_it : sfm_data_.poses)
//  {
//    const IndexT indexPose = pose_it.first;
//    std::cout << indexPose << " ";
//  }
//  std::cout << std::endl;


  return b_BA_Status;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_water_xy_z()
{
    {
      // Remove outliers (max_angle, residual error)
      const size_t pointcount_initial = sfm_data_.structure.size();
      RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
      const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
      RemoveOutliers_AngleError(sfm_data_, 2.0);
      const size_t pointcount_angular_filter = sfm_data_.structure.size();
      std::cout << "Outlier removal (remaining #points):\n"
        << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
        << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
        << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;
    }

  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

    std::cout << "this->b_use_motion_prior_ : " << this->b_use_motion_prior_ << std::endl;

    std::cout << "with pose id: " << std::endl;
    for (const auto & pose_it : sfm_data_.poses)
    {
      const IndexT indexPose = pose_it.first;
      std::cout << indexPose << " ";
    }
    std::cout << std::endl;

{
    // if need trans2plane
    Mat3 transR;
//    transR << 0.972892880440, 0.003697630484, 0.231226608157,
//            -0.032148022205, 0.992326200008, 0.119395114481,
//            -0.229010745883, -0.123592138290, 0.965546011925; //6

//    transR <<0.697977840900, 0.458458542824, 0.550129711628,
//            -0.019125504419, 0.779869437218, -0.625649988651,
//            -0.715863943100, 0.426168292761, 0.553099811077; //9

//    transR << 0.997934401035, -0.007309580687, 0.063823908567,
//            0.000000000000, 0.993505537510, 0.113783515990,
//            -0.064241118729, -0.113548487425, 0.991453409195; //7


//    transR << 0.949449598789, 0.007278595120, 0.313835054636,
//              -0.050830382854, 0.9901028871540, 0.130814984441,
//              -0.309776842594, -0.140154585242, 0.940422773361; //6_

//    transR << 0.947387576103, -0.050057616085, 0.316150307655,
//              0.003470299998, 0.989244163036, 0.146232604980,
//              -0.320069879293, -0.137441813946, 0.937371313572; //6_



//    transR << 0.821539223194, -0.006155421957, 0.5701187849040,
//              0.314995497465, 0.838380157948, -0.444855630398,
//              -0.475238025188, 0.5450511574750, 0.6907011270520; // 9_


//    transR << 0, 1, 0,
//              1, 0, 0,
//              0, 0, -1; // synthetic

    transR << 1, 0, 0,
              0, 1, 0,
              0, 0, 1; // real data

//    transR << 0.863876760006, 0.304972529411, 0.400884926319,
//              0.001411679666, 0.794406175613, -0.607385218143,
//             -0.503701269627, 0.525271892548, 0.685838520527 ; // 9_



    std::cout << "transR : " << transR << std::endl;

    for (auto & pose_it : sfm_data_.poses)
    {
      Pose3 & pose = pose_it.second;
      Mat3 R = pose.rotation();
      Vec3 C = pose.center();

      Mat3 newR = R*transR.transpose();
      Vec3 newC = transR*C;

      pose = Pose3(newR, newC);

    }
    //
    for (auto & structure_landmark_it : sfm_data_.structure)
    {
        structure_landmark_it.second.X = transR * structure_landmark_it.second.X;
    }
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_water_init_plane", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
}

    Bundle_Adjustment_Ceres bundle_adjustment_obj;
    bool b_BA_Status;
  // - refine only Structure and translations
  b_BA_Status = bundle_adjustment_obj.Adjust_water_xy
    (
      sfm_data_,
      Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::NONE,//ADJUST_TRANSLATION, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
        Structure_Parameter_Type::ADJUST_ALL,
        Control_Point_Parameter(),
        this->b_use_motion_prior_)
    );
  Save(sfm_data_,
    stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_water_refine_XY", "ply"),
    ESfM_Data(EXTRINSICS | STRUCTURE));

  std::cout << "Adjust_water_xy finished" << std::endl;


//  double water_plane = 1.04;
  //double water_plane = 0.9612;//1.0049;
  //double water_plane = 0.966;//1.0049;
  //double water_plane = 1.059828;//1.0049;
//  double water_plane = 1;//-0.9;//-1;//1.0049;
//  double water_plane = -0.35;//-0.9;//-1;//1.0049;
  double water_plane = 1.248;//-0.9;//-1;//1.0049;

  b_BA_Status = bundle_adjustment_obj.Adjust_water_z
    (
      sfm_data_,
      water_plane,
      Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::NONE,//ADJUST_TRANSLATION, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
        Structure_Parameter_Type::ADJUST_ALL,
        Control_Point_Parameter(),
        this->b_use_motion_prior_)
    );
  Save(sfm_data_,
    stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_water_refine_Z", "ply"),
    ESfM_Data(EXTRINSICS | STRUCTURE));

  std::cout << "Adjust_water_z finished" << std::endl;

  {
      // transform Z-axis
      Mat3 transR;
//      transR << 0.0, 1.0, 0.0,
//                1.0, 0.0, 0.0,
//                0.0, 0.0, -1.0;

      transR << 1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0;

      for (auto & pose_it : sfm_data_.poses)
      {
        Pose3 & pose = pose_it.second;
        Mat3 R = pose.rotation();
        Vec3 C = pose.center();

        Mat3 newR = R*transR.transpose();
        Vec3 newC = transR*C;

        pose = Pose3(newR, newC);

      }
      //
      for (auto & structure_landmark_it : sfm_data_.structure)
      {
          structure_landmark_it.second.X = transR * structure_landmark_it.second.X;
      }
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_water_init_plane_Z", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
  }
  return true;




  if (b_BA_Status)
  {
    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }
  // Check that poses & intrinsic cover some measures (after outlier removal)
  const IndexT minPointPerPose = 12; // 6 min
  const IndexT minTrackLength = 3; // 2 min
  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
  {
    // TODO: must ensure that track graph is producing a single connected component

    const size_t pointcount_cleaning = sfm_data_.structure.size();
    std::cout << "Point_cloud cleaning:\n"
      << "\t #3DPoints: " << pointcount_cleaning << "\n";
  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  const Optimize_Options ba_refine_options(
    ReconstructionEngine::intrinsic_refinement_options_,
    Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
    Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
    Control_Point_Parameter(),
    this->b_use_motion_prior_);

  b_BA_Status = bundle_adjustment_obj.Adjust(sfm_data_, ba_refine_options);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  std::cout << "with pose id: " << std::endl;
  for (const auto & pose_it : sfm_data_.poses)
  {
    const IndexT indexPose = pose_it.first;
    std::cout << indexPose << " ";
  }
  std::cout << std::endl;


  return b_BA_Status;
}

//输出重投影误差信息
void printReprojectError(SfM_Data& sfm_data)
{
    //重投影误差的和
    double sumError[2];
    sumError[0]=0;sumError[1]=0;
    //重投影误差的计算次数
    int allNum=0;
    for (auto & structure_landmark_it : sfm_data.structure)
    {
      const Observations & obs = structure_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        //当前位置的旋转矩阵
        Eigen::Matrix3d rotMat=sfm_data.poses[view->id_pose].rotation();
        //当前位置的平移
        Eigen::Vector3d translateVec=sfm_data.poses[view->id_pose].translation();
        //当前位置的3D点坐标
        Eigen::Vector3d pos3d=structure_landmark_it.second.X;
        //点坐标投影到相机坐标系下
        Eigen::Vector3d camPos=rotMat*pos3d+translateVec;
        //相机内参
        std::vector<double> intrinsicInfo=sfm_data.intrinsics[view->id_intrinsic]->getParams();
        //投影到图片上
        double reprojectPos[2];
        reprojectPos[0]=camPos[0]/camPos[2]*intrinsicInfo[0]+intrinsicInfo[1];
        reprojectPos[1]=camPos[1]/camPos[2]*intrinsicInfo[0]+intrinsicInfo[2];
        //计算重投影误差
        reprojectPos[0]=abs(reprojectPos[0]-obs_it.second.x[0]);
        reprojectPos[1]=abs(reprojectPos[1]-obs_it.second.x[1]);
        //输出重投影误差
        //std::cout<<reprojectPos[0]<<"\t\t"<<reprojectPos[1]<<"\n";
        //重投影误差计数
        sumError[0]+=reprojectPos[0];
        sumError[1]+=reprojectPos[1];
        //计算次数计数
        allNum++;

      }
    }
      //输出重投影误差的平均值
    std::cout<<"mean error:\n";
    std::cout<<sumError[0]/allNum<<"\t\t"<<sumError[1]/allNum<<"\n";
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_water_xy_z_n()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

    std::cout << "this->b_use_motion_prior_ : " << this->b_use_motion_prior_ << std::endl;
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_water_init_plane", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));

//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "sfm_data", ".bin"),
//      ESfM_Data(ALL));
    Bundle_Adjustment_Ceres bundle_adjustment_obj;
    bool b_BA_Status=true;
    //输出重投影误差信息
    //printReprojectError(sfm_data_);
    //过滤掉有效观测图太少的情况
   //sfm_data_.filterLossObv();
  // - refine only Structure and translations
//      b_BA_Status = bundle_adjustment_obj.Adjust_water_xy
//        (
//          sfm_data_,
//          Optimize_Options(
//            Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//            Extrinsic_Parameter_Type::ADJUST_ALL,//ADJUST_TRANSLATION, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//            Structure_Parameter_Type::ADJUST_ALL,
//            Control_Point_Parameter(),
//            this->b_use_motion_prior_)
//        );
    //指定根目录
    std::cout<<sfm_data_.s_root_path<<std::endl;
    sfm_data_.s_root_path="/media/cvlab/data/workSpace/mainProject/topViewConstruct/workspace/hekouzhen/imgs";
    //共享投影特征点
 //   sfm_data_.shareViewFeature();
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "sfmInit", ".bin"),
//      ESfM_Data(ALL));
    //载入view的图片
    sfm_data_.loadViewImgs();
    //按照dom图片上每个像素的顺序依次填充
    //sfm_data_.denseDomByDomPixel();
    try {
        sfm_data_.denseDomLikeMvs();
    } catch (int errorFlag) {
        //这里只负责显示拿到的错误信息
        std::cout<<errorFlag<<std::endl;
        throw errorFlag;
    }

    //释放图片内存
    sfm_data_.releaseImgs();

    //根据dom的分数，把它弄成点云，高度表示分数
//    sfm_data_.getScoreAsCloud();
//        Save(sfm_data_,
//          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "score_cloud", "ply"),
//          ESfM_Data(STRUCTURE));
    //根据dom图里面的坐标重构点云
    sfm_data_.getZAsCloud();
    //按照图片顺序往DOM上面贴图
      //sfm_data_.denseDomImageByImage();
    //做德劳内三角化
    //sfm_data_.delaunayTriangulation3D();
    //sfm_data_.delaunayTriangulation2D();
    //把点云稠密化
    //sfm_data_.densifyPointcloud();
    //保存最后的DOM结果
    sfm_data_.saveDomResult(stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "domResult", ".bmp"));
      //输出重投影误差信息
        //printReprojectError(sfm_data_);
      //去除部分坏点
    //  RemoveOutliers_PixelResidualError(sfm_data_, 500.0);//越小越严格
    //  RemoveOutliers_AngleError(sfm_data_, 2.0);//越大越严格
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_water_refine_XY", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
    //printReprojectError(sfm_data_);
//  std::cout << "Adjust_water_xy finished" << std::endl;
  //优化Z
//  bundle_adjustment_obj.adjustWaterZ(sfm_data_,
//                                     Optimize_Options(
//                                       Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//                                       Extrinsic_Parameter_Type::ADJUST_ALL,
//                                       Structure_Parameter_Type::ADJUST_ALL,
//                                       Control_Point_Parameter(),
//                                       this->b_use_motion_prior_));
  //输出重投影误差信息
  //printReprojectError(sfm_data_);
    return true;
//  double water_plane = 1.04;
  //double water_plane = 0.9612;//1.0049;
  //double water_plane = 0.966;//1.0049;
  //double water_plane = 1.059828;//1.0049;
//  double water_plane = 1.038;//73
//  double water_plane = 2.37;//72
  //double water_plane = 1.0075;//73
  //double water_plane = 1.0369;//71_1
  //double water_plane = 0.9;//run_9
//  double water_plane = 1.31;//run_6
    //double water_plane = -0.625;//c9_wp0.6
//    double water_plane = -0.352;//c9_wp0.35
//    double water_plane = -0.352;//c6_wp0.35
//    double water_plane = -0.625;//c6_wp0.6
    double water_plane = 0;//run_6
//    double water_plane = -0.4;//run_6

//    b_BA_Status = bundle_adjustment_obj.Adjust_water_min_z
//      (
//        sfm_data_,
//        water_plane,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::NONE,//ADJUST_TRANSLATION, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );

//    std::cout << "water_plane : " << water_plane << std::endl;

//  b_BA_Status = bundle_adjustment_obj.Adjust_water_z_fixC1
//    (
//      sfm_data_,
//      water_plane,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::NONE,//ADJUST_TRANSLATION, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );
//  Save(sfm_data_,
//    stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_water_refine_Z", "ply"),
//    ESfM_Data(EXTRINSICS | STRUCTURE));

//  b_BA_Status = bundle_adjustment_obj.Adjust_water_min_z
//    (
//      sfm_data_,
//      water_plane,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::NONE,//ADJUST_TRANSLATION, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );
//  Save(sfm_data_,
//    stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_water_refine_Z2", "ply"),
//    ESfM_Data(EXTRINSICS | STRUCTURE));

  std::cout << "water_plane : " << water_plane << std::endl;

  b_BA_Status = bundle_adjustment_obj.Adjust_water_z_fixC1
    (
      sfm_data_,
      water_plane,
      Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::NONE,//ADJUST_TRANSLATION, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
        Structure_Parameter_Type::ADJUST_ALL,
        Control_Point_Parameter(),
        this->b_use_motion_prior_)
    );
  Save(sfm_data_,
    stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_water_refine_Z3", "ply"),
    ESfM_Data(EXTRINSICS | STRUCTURE));

  double init_ref_n = 1.33;

  b_BA_Status = bundle_adjustment_obj.Adjust_water_z_n
    (
      sfm_data_,
      water_plane,
      init_ref_n,
      Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::NONE,//ADJUST_TRANSLATION, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
        Structure_Parameter_Type::ADJUST_ALL,
        Control_Point_Parameter(),
        this->b_use_motion_prior_)
    );
  Save(sfm_data_,
    stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_water_refine_Zn", "ply"),
    ESfM_Data(EXTRINSICS | STRUCTURE));

  std::cout << "water_plane : " << water_plane << std::endl;
  std::cout << "Adjust_water_z_n n = " << init_ref_n  << std::endl;
  std::cout << "Adjust_water_z_n finished" << std::endl;

  {
      // transform Z-axis
      Mat3 transR;
      transR << 0.0, 1.0, 0.0,
                1.0, 0.0, 0.0,
                0.0, 0.0, -1.0;
      transR << 1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0;

      for (auto & pose_it : sfm_data_.poses)
      {
        Pose3 & pose = pose_it.second;
        Mat3 R = pose.rotation();
        Vec3 C = pose.center();

        Mat3 newR = R*transR.transpose();
        Vec3 newC = transR*C;

        pose = Pose3(newR, newC);

      }
      //
      for (auto & structure_landmark_it : sfm_data_.structure)
      {
          structure_landmark_it.second.X = transR * structure_landmark_it.second.X;
      }
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_water_init_plane_Z", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
  }
  return true;




  if (b_BA_Status)
  {
    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }
  // Check that poses & intrinsic cover some measures (after outlier removal)
  const IndexT minPointPerPose = 12; // 6 min
  const IndexT minTrackLength = 3; // 2 min
  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
  {
    // TODO: must ensure that track graph is producing a single connected component

    const size_t pointcount_cleaning = sfm_data_.structure.size();
    std::cout << "Point_cloud cleaning:\n"
      << "\t #3DPoints: " << pointcount_cleaning << "\n";
  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  const Optimize_Options ba_refine_options(
    ReconstructionEngine::intrinsic_refinement_options_,
    Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
    Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
    Control_Point_Parameter(),
    this->b_use_motion_prior_);

  b_BA_Status = bundle_adjustment_obj.Adjust(sfm_data_, ba_refine_options);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  std::cout << "with pose id: " << std::endl;
  for (const auto & pose_it : sfm_data_.poses)
  {
    const IndexT indexPose = pose_it.first;
    std::cout << indexPose << " ";
  }
  std::cout << std::endl;


  return b_BA_Status;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_step21()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):
  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = bundle_adjustment_obj.Adjust
    (
      sfm_data_,
      Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant// // Rotations are held as constant
        Structure_Parameter_Type::ADJUST_ALL,
        Control_Point_Parameter(),
        this->b_use_motion_prior_)
    );
  if (b_BA_Status)
  {
    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }
  // Check that poses & intrinsic cover some measures (after outlier removal)
  const IndexT minPointPerPose = 12; // 6 min
  const IndexT minTrackLength = 3; // 2 min
  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
  {
    // TODO: must ensure that track graph is producing a single connected component

    const size_t pointcount_cleaning = sfm_data_.structure.size();
    std::cout << "Point_cloud cleaning:\n"
      << "\t #3DPoints: " << pointcount_cleaning << "\n";
  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  const Optimize_Options ba_refine_options(
    ReconstructionEngine::intrinsic_refinement_options_,
    Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
    Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
    Control_Point_Parameter(),
    this->b_use_motion_prior_);

  b_BA_Status = bundle_adjustment_obj.Adjust(sfm_data_, ba_refine_options);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  return b_BA_Status;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_step22()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):
  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = bundle_adjustment_obj.Adjust
    (
      sfm_data_,
      Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_ALL, // Rotations are held as constant// // Rotations are held as constant
        Structure_Parameter_Type::ADJUST_ALL,
        Control_Point_Parameter(),
                  true)
//        this->b_use_motion_prior_)
    );

  return true;

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL, // Rotations are held as constant// // Rotations are held as constant
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
                  false)

//          this->b_use_motion_prior_)
      );

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_step22", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }
  return b_BA_Status;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_GCP()
{  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):
//    std::cout << "with pose id: " << std::endl;
//    for (const auto & pose_it : sfm_data_.poses)
//    {
//      const IndexT indexPose = pose_it.first;
//      std::cout << indexPose << " ";
//    }
//    std::cout << std::endl;

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = bundle_adjustment_obj.Adjust
    (
      sfm_data_,
      Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
        Structure_Parameter_Type::ADJUST_ALL,
        Control_Point_Parameter(),
        this->b_use_motion_prior_)
    );

  if (b_BA_Status)
  {
    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  // Check that poses & intrinsic cover some measures (after outlier removal)
  const IndexT minPointPerPose = 12; // 6 min
  const IndexT minTrackLength = 3; // 2 min
  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
  {
    // TODO: must ensure that track graph is producing a single connected component

    const size_t pointcount_cleaning = sfm_data_.structure.size();
    std::cout << "Point_cloud cleaning:\n"
      << "\t #3DPoints: " << pointcount_cleaning << "\n";
  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  const Optimize_Options ba_refine_options(
    ReconstructionEngine::intrinsic_refinement_options_,
    Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
    Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
    Control_Point_Parameter(),
    this->b_use_motion_prior_);

  b_BA_Status = bundle_adjustment_obj.Adjust(sfm_data_, ba_refine_options);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

//  std::cout << "with pose id: " << std::endl;
//  for (const auto & pose_it : sfm_data_.poses)
//  {
//    const IndexT indexPose = pose_it.first;
//    std::cout << indexPose << " ";
//  }
//  std::cout << std::endl;

  ///apply GCPs
  ///
  ///
  ///

  //---
  // registration (coarse):
  // - compute the 3D points corresponding to the control point observation for the SfM scene
  // - compute a coarse registration between the controls points & the triangulated point
  // - transform the scene according the found transformation
  //---
  std::cout << "error 1" << std::endl;
  Hash_Map<IndexT, Vec3> vec_control_points, vec_triangulated;
  Hash_Map<IndexT, double> vec_triangulation_errors;
//  std::map<IndexT, Vec3> vec_control_points, vec_triangulated;
//  std::map<IndexT, double> vec_triangulation_errors;
  for (Landmarks::iterator iterCP = sfm_data_.control_points.begin();
    iterCP != sfm_data_.control_points.end(); ++iterCP)
  {
    Landmark & landmark = iterCP->second;
    //Triangulate the point:
    Triangulation trianObj;
    const Observations & obs = landmark.obs;
    for (Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      const View * view = sfm_data_.views.at(itObs->first).get();
      if (!sfm_data_.IsPoseAndIntrinsicDefined(view))
        continue;
      const openMVG::cameras::IntrinsicBase * cam = sfm_data_.GetIntrinsics().at(view->id_intrinsic).get();
      const openMVG::geometry::Pose3 pose = sfm_data_.GetPoseOrDie(view);
      trianObj.add(
        cam->get_projective_equivalent(pose),
        cam->get_ud_pixel(itObs->second.x));
    }
    // Compute the 3D point
    const Vec3 X = trianObj.compute();
    if (trianObj.minDepth() > 0) // Keep the point only if it have a positive depth
    {
      vec_triangulated[iterCP->first] = X;
      vec_control_points[iterCP->first] = landmark.X;
      vec_triangulation_errors[iterCP->first] = trianObj.error()/(double)trianObj.size();
    }
    else
    {
        std::cout << "error 1.1" << std::endl;
//      std::cout << "Invalid triangulation" << std::endl;
//      return;
    }
  }

  // compute the similarity
  {
    // data conversion to appropriate container
    Mat x1(3, vec_control_points.size()),
        x2(3, vec_control_points.size());
    for (int i = 0; i < vec_control_points.size(); ++i)
    {
      Hash_Map<IndexT, Vec3>::iterator it_vcp = vec_control_points.begin();
      std::advance(it_vcp, i);
      Hash_Map<IndexT, Vec3>::iterator it_vt = vec_triangulated.begin();
      std::advance(it_vt, i);
      x1.col(i) = it_vt->second;//vec_triangulated[i];
      x2.col(i) = it_vcp->second;//vec_control_points[i];
    }

    std::cout
      << "Control points observation triangulations:\n"
      << x1 << std::endl << std::endl
      << "Control points coords:\n"
      << x2 << std::endl << std::endl;

    Vec3 t;
    Mat3 R;
    double S;
    if (openMVG::geometry::FindRTS(x1, x2, &S, &t, &R))
    {
      openMVG::geometry::Refine_RTS(x1,x2,&S,&t,&R);
      std::cout << "Found transform:\n"
        << " scale: " << S << "\n"
        << " rotation:\n" << R << "\n"
        << " translation: "<< t.transpose() << std::endl;


      //--
      // Apply the found transformation as a 3D Similarity transformation matrix // S * R * X + t
      //--

      const openMVG::geometry::Similarity3 sim(geometry::Pose3(R, -R.transpose() * t/S), S);
      openMVG::sfm::ApplySimilarity(sim, sfm_data_);

      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_05_Similarity", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));

      // Display some statistics:
      std::stringstream os;
      for (Landmarks::const_iterator iterL = sfm_data_.control_points.begin();
        iterL != sfm_data_.control_points.end(); ++iterL)
      {
        const IndexT CPIndex = iterL->first;
        // If the control point has not been used, continue...
        if (vec_triangulation_errors.find(CPIndex) == vec_triangulation_errors.end())
          continue;

        os
          << "CP index: " << CPIndex << "\n"
          << "CP triangulation error: " << vec_triangulation_errors[CPIndex] << " pixel(s)\n"
          << "CP registration error: "
          << (sim(vec_triangulated[CPIndex]) - vec_control_points[CPIndex]).norm() << " user unit(s)"<< "\n\n";
      }
      std::cout << os.str();

//      QMessageBox msgBox;
//      msgBox.setText(QString::fromStdString(string_pattern_replace(os.str(), "\n", "<br>")));
//      msgBox.exec();
    }
    else
    {
//      QMessageBox msgBox;
//      msgBox.setText("Registration failed. Please check your Control Points coordinates.");
//      msgBox.exec();
    }
  }

  //---
  // Bundle adjustment with GCP
  //---
  {
    using namespace openMVG::sfm;
    Bundle_Adjustment_Ceres::BA_Ceres_options options;
    Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
    Control_Point_Parameter control_point_opt(20.0, true);
    if (!bundle_adjustment_obj.Adjust(sfm_data_,
        Optimize_Options
        (
          cameras::Intrinsic_Parameter_Type::NONE, // Keep intrinsic constant
          Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
          Structure_Parameter_Type::ADJUST_ALL, // Adjust structure
          control_point_opt // Use GCP and weight more their observation residuals
          )
        )
      )
    {
//      QMessageBox msgBox;
//      msgBox.setText("BA with GCP failed.");
//      msgBox.exec();
    }

    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_06_GCP", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));

  }

  b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),//(20.0, true),
          this->b_use_motion_prior_)
      );
  Save(sfm_data_,
    stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_07_after_GCP", "ply"),
    ESfM_Data(EXTRINSICS | STRUCTURE));



  return true;//b_BA_Status;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_GCP_step2(std::string outPath)
{
  ///apply GCPs
  ///
  ///
  ///

  bool b_BA_Status = true;
  //---
  // registration (coarse):
  // - compute the 3D points corresponding to the control point observation for the SfM scene
  // - compute a coarse registration between the controls points & the triangulated point
  // - transform the scene according the found transformation
  //---
  Hash_Map<IndexT, Vec3> vec_control_points, vec_triangulated;
  Hash_Map<IndexT, double> vec_triangulation_errors;
//  std::map<IndexT, Vec3> vec_control_points, vec_triangulated;
//  std::map<IndexT, double> vec_triangulation_errors;
  for (Landmarks::iterator iterCP = sfm_data_.control_points.begin();
    iterCP != sfm_data_.control_points.end(); ++iterCP)
  {
    Landmark & landmark = iterCP->second;
    //Triangulate the point:
    Triangulation trianObj;
    const Observations & obs = landmark.obs;
    for (Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      const View * view = sfm_data_.views.at(itObs->first).get();
      if (!sfm_data_.IsPoseAndIntrinsicDefined(view))
        continue;
      const openMVG::cameras::IntrinsicBase * cam = sfm_data_.GetIntrinsics().at(view->id_intrinsic).get();
      const openMVG::geometry::Pose3 pose = sfm_data_.GetPoseOrDie(view);
      trianObj.add(
        cam->get_projective_equivalent(pose),
        cam->get_ud_pixel(itObs->second.x));
    }

    // Compute the 3D point
    const Vec3 X = trianObj.compute();
    if (trianObj.minDepth() > 0) // Keep the point only if it have a positive depth
    {
      vec_triangulated[iterCP->first] = X;
      vec_control_points[iterCP->first] = landmark.X;
      vec_triangulation_errors[iterCP->first] = trianObj.error()/(double)trianObj.size();
    }
    else
    {
        b_BA_Status = false;
//      std::cout << "Invalid triangulation" << std::endl;
//      return;
    }
  }

  // compute the similarity
  {
    // data conversion to appropriate container
    Mat x1(3, vec_control_points.size()),
        x2(3, vec_control_points.size());
    for (int i = 0; i < vec_control_points.size(); ++i)
    {
      Hash_Map<IndexT, Vec3>::iterator it_vcp = vec_control_points.begin();
      std::advance(it_vcp, i);
      Hash_Map<IndexT, Vec3>::iterator it_vt = vec_triangulated.begin();
      std::advance(it_vt, i);
      x1.col(i) = it_vt->second;//vec_triangulated[i];
      x2.col(i) = it_vcp->second;//vec_control_points[i];
    }

    std::cout
      << "Control points observation triangulations:\n"
      << x1 << std::endl << std::endl
      << "Control points coords:\n"
      << x2 << std::endl << std::endl;

    Vec3 t;
    Mat3 R;
    double S;
    if (openMVG::geometry::FindRTS(x1, x2, &S, &t, &R))
    {
      openMVG::geometry::Refine_RTS(x1,x2,&S,&t,&R);
      std::cout << "Found transform:\n"
        << " scale: " << S << "\n"
        << " rotation:\n" << R << "\n"
        << " translation: "<< t.transpose() << std::endl;


      //--
      // Apply the found transformation as a 3D Similarity transformation matrix // S * R * X + t
      //--

      const openMVG::geometry::Similarity3 sim(geometry::Pose3(R, -R.transpose() * t/S), S);
      openMVG::sfm::ApplySimilarity(sim, sfm_data_);

      Save(sfm_data_,
        stlplus::create_filespec(outPath, "structure_05_Similarity", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));

      // Display some statistics:
      std::stringstream os;
      for (Landmarks::const_iterator iterL = sfm_data_.control_points.begin();
        iterL != sfm_data_.control_points.end(); ++iterL)
      {
        const IndexT CPIndex = iterL->first;
        // If the control point has not been used, continue...
        if (vec_triangulation_errors.find(CPIndex) == vec_triangulation_errors.end())
          continue;

        os
          << "CP index: " << CPIndex << "\n"
          << "CP triangulation error: " << vec_triangulation_errors[CPIndex] << " pixel(s)\n"
          << "CP registration error: "
          << (sim(vec_triangulated[CPIndex]) - vec_control_points[CPIndex]).norm() << " user unit(s)"<< "\n\n";
      }
      std::cout << os.str();




//      QMessageBox msgBox;
//      msgBox.setText(QString::fromStdString(string_pattern_replace(os.str(), "\n", "<br>")));
//      msgBox.exec();
    }
    else
    {
        std::cout << "Registration failed. Please check your Control Points coordinates." << std::endl;
        b_BA_Status = false;
//      QMessageBox msgBox;
//      msgBox.setText("Registration failed. Please check your Control Points coordinates.");
//      msgBox.exec();
    }
  }

  //---
  // Bundle adjustment with GCP
  //---
  {
    using namespace openMVG::sfm;
    Bundle_Adjustment_Ceres::BA_Ceres_options options;
    Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
    Control_Point_Parameter control_point_opt(20.0, true);
    if (!bundle_adjustment_obj.Adjust(sfm_data_,
        Optimize_Options
        (
          cameras::Intrinsic_Parameter_Type::NONE, // Keep intrinsic constant
          Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
          Structure_Parameter_Type::ADJUST_ALL, // Adjust structure
          control_point_opt // Use GCP and weight more their observation residuals
          )
        )
      )
    {
        std::cout << "BA with GCP failed." << std::endl;
//      QMessageBox msgBox;
//      msgBox.setText("BA with GCP failed.");
//      msgBox.exec();
    }

    Save(sfm_data_,
      stlplus::create_filespec(outPath, "structure_06_GCP_1", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));

    // Remove outliers (max_angle, residual error)
    const size_t pointcount_initial = sfm_data_.structure.size();
    RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
    const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
    RemoveOutliers_AngleError(sfm_data_, 2.0);
    const size_t pointcount_angular_filter = sfm_data_.structure.size();
    std::cout << "Outlier removal (remaining #points):\n"
      << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
      << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
      << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

//    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(outPath, "structure_06_gcp_removed", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }


    if (!bundle_adjustment_obj.Adjust(sfm_data_,
        Optimize_Options
        (
          cameras::Intrinsic_Parameter_Type::NONE, // Keep intrinsic constant
          Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
          Structure_Parameter_Type::ADJUST_ALL, // Adjust structure
          control_point_opt // Use GCP and weight more their observation residuals
          )
        )
      )
    {
        std::cout << "BA with GCP failed." << std::endl;
//      QMessageBox msgBox;
//      msgBox.setText("BA with GCP failed.");
//      msgBox.exec();
    }


    Save(sfm_data_,
      stlplus::create_filespec(outPath, "structure_06_GCP", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));

  }


  return b_BA_Status;//b_BA_Status;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_init()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );

//  // Remove outliers (max_angle, residual error)
//  const size_t pointcount_initial1 = sfm_data_.structure.size();
//  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//  const size_t pointcount_pixelresidual_filter1 = sfm_data_.structure.size();
//  RemoveOutliers_AngleError(sfm_data_, 2.0);
//  const size_t pointcount_angular_filter1 = sfm_data_.structure.size();
//  std::cout << "Outlier removal (remaining #points):\n"
//    << "\t initial structure size #3DPoints: " << pointcount_initial1 << "\n"
//    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter1 << "\n"
//    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter1 << std::endl;


//  bool b_BA_Status// = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::NONE, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );


//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  if (b_BA_Status)
  {
    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }


  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  Optimize_Options ba_refine_options;
  if (ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
      ba_refine_options = Optimize_Options(
        ReconstructionEngine::intrinsic_refinement_options_,
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }else{
      ba_refine_options = Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }

  b_BA_Status = bundle_adjustment_obj.Adjust(sfm_data_, ba_refine_options);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  return b_BA_Status;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_init_BAStep()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );

//  // Remove outliers (max_angle, residual error)
//  const size_t pointcount_initial1 = sfm_data_.structure.size();
//  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//  const size_t pointcount_pixelresidual_filter1 = sfm_data_.structure.size();
//  RemoveOutliers_AngleError(sfm_data_, 2.0);
//  const size_t pointcount_angular_filter1 = sfm_data_.structure.size();
//  std::cout << "Outlier removal (remaining #points):\n"
//    << "\t initial structure size #3DPoints: " << pointcount_initial1 << "\n"
//    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter1 << "\n"
//    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter1 << std::endl;


//  bool b_BA_Status// = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::NONE, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );


//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  int initSteps = 30;
  for(int step = 0; step < initSteps; ++step)
  {
      bool b_BA_Status = bundle_adjustment_obj.Adjust_BAStep
        (
          sfm_data_,
          Optimize_Options(
            Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
            Extrinsic_Parameter_Type::ADJUST_ALL,
            Structure_Parameter_Type::ADJUST_ALL,
            Control_Point_Parameter(),
            this->b_use_motion_prior_)
        );
      std::stringstream ss;
      ss << step;
      std::string s;
      ss >> s;
      if (b_BA_Status && !sLogging_file_.empty())
      {
        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_"+s+"_step_init", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));
      }

  }

  return true;
//  if (b_BA_Status)
  {
    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }


  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  Optimize_Options ba_refine_options;
  if (ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
      ba_refine_options = Optimize_Options(
        ReconstructionEngine::intrinsic_refinement_options_,
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }else{
      ba_refine_options = Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }

  b_BA_Status = bundle_adjustment_obj.Adjust(sfm_data_, ba_refine_options);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  return b_BA_Status;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_step23(const std::map<int, bool>& fixedImgIdList, const std::map<int, bool>& fixedXIdList)
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;

  bool b_BA_Status = true;
  if (b_BA_Status)
  {
    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust_step23
      (
        sfm_data_,
        fixedImgIdList,
        fixedXIdList,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }


  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust_step23
      (
        sfm_data_,
        fixedImgIdList,
        fixedXIdList,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  Optimize_Options ba_refine_options;
  if (ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
      ba_refine_options = Optimize_Options(
        ReconstructionEngine::intrinsic_refinement_options_,
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }else{
      ba_refine_options = Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }

  b_BA_Status = bundle_adjustment_obj.Adjust_step23(sfm_data_, fixedImgIdList, fixedXIdList, ba_refine_options);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  return b_BA_Status;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_init_4_margeTest()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );

//  // Remove outliers (max_angle, residual error)
//  const size_t pointcount_initial1 = sfm_data_.structure.size();
//  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//  const size_t pointcount_pixelresidual_filter1 = sfm_data_.structure.size();
//  RemoveOutliers_AngleError(sfm_data_, 2.0);
//  const size_t pointcount_angular_filter1 = sfm_data_.structure.size();
//  std::cout << "Outlier removal (remaining #points):\n"
//    << "\t initial structure size #3DPoints: " << pointcount_initial1 << "\n"
//    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter1 << "\n"
//    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter1 << std::endl;


//  bool b_BA_Status// = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::NONE, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );


//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  if (b_BA_Status)
  {
    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }


  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  Optimize_Options ba_refine_options;
  if (ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
      ba_refine_options = Optimize_Options(
        ReconstructionEngine::intrinsic_refinement_options_,
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }else{
      ba_refine_options = Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }

  b_BA_Status = bundle_adjustment_obj.Adjust(sfm_data_, ba_refine_options);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  return b_BA_Status;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_init_water_c1()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );

//  // Remove outliers (max_angle, residual error)
//  const size_t pointcount_initial1 = sfm_data_.structure.size();
//  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//  const size_t pointcount_pixelresidual_filter1 = sfm_data_.structure.size();
//  RemoveOutliers_AngleError(sfm_data_, 2.0);
//  const size_t pointcount_angular_filter1 = sfm_data_.structure.size();
//  std::cout << "Outlier removal (remaining #points):\n"
//    << "\t initial structure size #3DPoints: " << pointcount_initial1 << "\n"
//    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter1 << "\n"
//    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter1 << std::endl;


//  bool b_BA_Status// = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::NONE, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );


//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  if (b_BA_Status)
  {
    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

//    b_BA_Status = bundle_adjustment_obj.Adjust_water_prepare
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_TRANSLATION,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
//    b_BA_Status = bundle_adjustment_obj.Adjust_water_prepare
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_ROTATION,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_R_Xi", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));

//    bundle_adjustment_obj.Adjust_water_c1_S(
//                sfm_data_,
//                Optimize_Options(
//                    Intrinsic_Parameter_Type::NONE,
//                    Extrinsic_Parameter_Type::ADJUST_ROTATION,
//                    Structure_Parameter_Type::ADJUST_ALL,
//                    Control_Point_Parameter(),
//                    this->b_use_motion_prior_
//                    )
//                );
    b_BA_Status = bundle_adjustment_obj.Adjust_water_c1_S//_water_prepare
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_RT_Xi", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
    std::cout << "return here" << std::endl;
//    return true;
//    {
//        // Remove outliers (max_angle, residual error)
//        const size_t pointcount_initial = sfm_data_.structure.size();
//        RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//        const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
//        RemoveOutliers_AngleError(sfm_data_, 2.0);
//        const size_t pointcount_angular_filter = sfm_data_.structure.size();
//        std::cout << "Outlier removal (remaining #points):\n"
//          << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
//          << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
//          << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;
//    }

    {
//        for (const auto & pose_it : sfm_data_.poses)
//        {
////          const IndexT indexPose = pose_it.first;

//          const Pose3 & pose = pose_it.second;
//          const Mat3 R = pose.rotation();
//          const Vec3 t = pose.translation();
//          std::cout << "id : " << pose_it.first << std::endl;
//          std::cout << "R : " << R << std::endl;
//          std::cout << "t : " << t << std::endl;
//          std::cout << "****-----------------****" << std::endl;
//        }

    }

    // - refine only Structure and Rotations & translations
//    b_BA_Status = bundle_adjustment_obj.Adjust_water_c1_S
//      (
//        sfm_data_,
//        Optimize_Options(
//          ReconstructionEngine::intrinsic_refinement_options_,  //Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    if (b_BA_Status && !sLogging_file_.empty())
//    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_020_refine_T_Xi_water", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }
//    // - refine only Structure and Rotations & translations
//    b_BA_Status = bundle_adjustment_obj.Adjust_water_c1_S
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    if (b_BA_Status && !sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_11_refine_R_Xi_water", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }
//    // - refine only Structure and Rotations & translations
//    b_BA_Status = bundle_adjustment_obj.Adjust_water_c1_S
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    if (b_BA_Status && !sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_12_refine_RT_Xi_water", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }
  }
  std::cout << "water C1" << std::endl;

//  // - refine only Structure and Rotations & translations
//  b_BA_Status = bundle_adjustment_obj.Adjust_water_c1
//    (
//      sfm_data_,
//      Optimize_Options(
//        ReconstructionEngine::intrinsic_refinement_options_,//Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_ALL,
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );
//  if (b_BA_Status && !sLogging_file_.empty())
//  {
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_RT_Xi", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
//  }

  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
//    b_BA_Status = bundle_adjustment_obj.Adjust_water_c1
//      (
//        sfm_data_,
//        Optimize_Options(
//          ReconstructionEngine::intrinsic_refinement_options_,
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    if (b_BA_Status && !sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

//  sfm_data_.poses.erase(15);
//  sfm_data_.views.erase(15);
//  for (auto & structure_landmark_it : sfm_data_.structure)
//  {
//    Observations & obs = structure_landmark_it.second.obs;

//    for (auto & obs_it : obs)
//    {
//      // Build the residual block corresponding to the track observation:
//        if(obs_it.first == 15)
//        {
//            obs.erase(obs_it.first);
//            break;
//        }
////      const View * view = sfm_data.views.at(obs_it.first).get();
//    }
//  }


  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_15_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }
//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  Optimize_Options ba_refine_options;
  if (ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
      ba_refine_options = Optimize_Options(
        ReconstructionEngine::intrinsic_refinement_options_,
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }else{
      ba_refine_options = Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }

  b_BA_Status = bundle_adjustment_obj.Adjust_water_c1_S(sfm_data_, ba_refine_options);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  return b_BA_Status;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_init_water2_c1()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_ALL, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );

//  b_BA_Status //= true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_ROTATION, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );

//  b_BA_Status// = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_ALL, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );


//  // Remove outliers (max_angle, residual error)
//  const size_t pointcount_initial1 = sfm_data_.structure.size();
//  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//  const size_t pointcount_pixelresidual_filter1 = sfm_data_.structure.size();
//  RemoveOutliers_AngleError(sfm_data_, 2.0);
//  const size_t pointcount_angular_filter1 = sfm_data_.structure.size();
//  std::cout << "Outlier removal (remaining #points):\n"
//    << "\t initial structure size #3DPoints: " << pointcount_initial1 << "\n"
//    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter1 << "\n"
//    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter1 << std::endl;




//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }
  {

//      std::ofstream out;
//      out.open("/home/guang/data/testWaterStep/imgs/m_out2/out_.txt");
//      std::cout <<"out ****-----------------****" << std::endl;
//      {
//          for (auto & structure_landmark_it : sfm_data_.structure)
//          {
//            const Observations & obs = structure_landmark_it.second.obs;

//            for (const auto & obs_it : obs)
//            {
//              // Build the residual block corresponding to the track observation:
//              const View * view = sfm_data_.views.at(obs_it.first).get();

//    //          // Each Residual block takes a point and a camera as input and outputs a 2
//    //          // dimensional residual. Internally, the cost function stores the observed
//    //          // image location and compares the reprojection against the observation.
//    //          ceres::CostFunction* cost_function =
//    //                  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_2::Create(obs_it.second.x);
//    //    //        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);


//              //#pragma omp critical
//              //      {

////              std::cout <<"out 1" << std::endl;
//              const Mat3 Ro = sfm_data_.poses[view->id_pose].rotation();
//              const Vec3 t = sfm_data_.poses[view->id_pose].translation();
//              double angleAxis[3];
//              ceres::RotationMatrixToAngleAxis((const double*)Ro.data(), angleAxis);

////              std::cout <<"out 2" << std::endl;
//                  const double * cam_R = angleAxis;
//                  const Vec3 cam_t = t;

//                  double pos_3dpoint[3];
////                  std::cout <<"out 2" << std::endl;
//                  pos_3dpoint[0] = structure_landmark_it.second.X(0);
//                  pos_3dpoint[1] = structure_landmark_it.second.X(1);
//                  pos_3dpoint[2] = structure_landmark_it.second.X(2);

//                  double pos_proj[3];
//                  double X[3];
//                  X[0] = pos_3dpoint[0];
//                  X[1] = pos_3dpoint[1];
//                  X[2] = 0.0;
//                  // Rotate the point according the camera rotation
//                  ceres::AngleAxisRotatePoint(cam_R, X, pos_proj);

////                  std::cout <<"out 3" << std::endl;
//                  // Apply the camera translation
//                  pos_proj[0] += cam_t(0);
//                  pos_proj[1] += cam_t(1);
//                  pos_proj[2] += cam_t(2);

//                  //prepare l
//                  double R[9];
//                  ceres::AngleAxisToRotationMatrix(cam_R, R);

////                  std::cout <<"out 4" << std::endl;
//                  std::vector<double> it_intrinsics = sfm_data_.intrinsics[view->id_intrinsic]->getParams();

//                  const double f = it_intrinsics[0]; //focal
//                  const double u = it_intrinsics[1]; //principal_point_x
//                  const double v = it_intrinsics[2]; //principal_point_y

////                  std::cout <<"out 5" << std::endl;
//                  // Transform the point from homogeneous to euclidean (undistorted point)
//                  const double x_i = u + f * (pos_proj[0] / pos_proj[2]);
//                  const double y_i = v + f * (pos_proj[1] / pos_proj[2]);

//                  const double P1 = f*R[0] + u*R[2];
//                  const double P2 = f*R[3] + u*R[5];
//                  const double P3 = f*R[6] + u*R[8];
//                  const double P4 = f*cam_t[0] + u*cam_t[2];
//                  const double P5 = f*R[1] + v*R[2];
//                  const double P6 = f*R[4] + v*R[5];
//                  const double P7 = f*R[7] + v*R[8];
//                  const double P8 = f*cam_t[1] + v*cam_t[2];
//                  const double P9 = R[2];
//                  const double P10 = R[5];
//                  const double P11 = R[8];
//                  const double P12 = cam_t[2];


////                  std::cout <<"out 6" << std::endl;

//                  const double xu = obs_it.second.x(0);
////                  std::cout <<"out 6 1" << std::endl;
//                  const double yu = obs_it.second.x(1);
////                  std::cout <<"out 6 2" << std::endl;

//                  const double l0 = P7 - P11*yu;
////                  std::cout <<"out 6 3" << std::endl;
//                  const double l1 = P11*xu - P3;
////                  std::cout <<"out 6 4" << std::endl;
//                  const double l2 = P3*yu - P7*xu;


////                  std::cout <<"out 7" << std::endl;

//                  double out_residuals = (x_i*l0 + y_i*l1 + l2)/sqrt(l0*l0 + l1*l1);
//              //    out_residuals[1] = T(0.0);
//              //    out_residuals[0] = abs(x_i*l0 + y_i*l1 + l2)/sqrt(l0*l0 + l1*l1);

////                  std::cout <<"out 8" << std::endl;
//                  if(out_residuals > 0.1)
//                  {
//                      out << "cam_R : " << cam_R[0] << " " << cam_R[1] << " " << cam_R[2] << std::endl;
//                      out << "R : " << R[0] << " " << R[1] << " " << R[2]
//                                          << R[3] << " " << R[4] << " " << R[5]
//                                          << R[6] << " " << R[7] << " " << R[8] << std::endl;
//                      out << "cam_t : " << cam_t(0) << " " << cam_t(1) << " " << cam_t(2) << std::endl;
//                      out << "f : " << f << " u : " << u << " v : " << v << std::endl;
//                      out << "xu : " << xu << " yu : " << yu << std::endl;
//                      out << "pos_3dpoint[0]: " << pos_3dpoint[0] << " pos_3dpoint[1]: " << pos_3dpoint[1] << " pos_3dpoint[2]: " << pos_3dpoint[2] << std::endl;
//                      out << "out_residuals[0]: "  << out_residuals << std::endl;
//                      out << "****-----------------****" << std::endl;
//                  }



//            }


//          }
//      }

//      out.close();
  }

  if (b_BA_Status)
  {
//    if (!sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

//    b_BA_Status = bundle_adjustment_obj.Adjust
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, //ReconstructionEngine::intrinsic_refinement_options_, //Intrinsic_Parameter_Type::NONE, //ReconstructionEngine::intrinsic_refinement_options_, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "test_01_refine_RT_Xi", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));

//    {
//        // Remove outliers (max_angle, residual error)
//        const size_t pointcount_initial = sfm_data_.structure.size();
//        RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//        const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
//        RemoveOutliers_AngleError(sfm_data_, 2.0);
//        const size_t pointcount_angular_filter = sfm_data_.structure.size();
//        std::cout << "Outlier removal (remaining #points):\n"
//          << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
//          << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
//          << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;
//        Save(sfm_data_,
//          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "test_02_refine_RT_Xi", "ply"),
//          ESfM_Data(EXTRINSICS | STRUCTURE));
//    }
//    b_BA_Status = bundle_adjustment_obj.Adjust
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, //ReconstructionEngine::intrinsic_refinement_options_, //Intrinsic_Parameter_Type::NONE, //ReconstructionEngine::intrinsic_refinement_options_, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "test_03_refine_RT_Xi", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));


    {

//        std::ofstream out;
//        out.open("/home/guang/data/testWaterStep/imgs/m_out2/out_0.txt");
//        std::cout <<"out ****-----------------****" << std::endl;
//        {
//            for (auto & structure_landmark_it : sfm_data_.structure)
//            {
//              const Observations & obs = structure_landmark_it.second.obs;

//              for (const auto & obs_it : obs)
//              {
//                // Build the residual block corresponding to the track observation:
//                const View * view = sfm_data_.views.at(obs_it.first).get();

//      //          // Each Residual block takes a point and a camera as input and outputs a 2
//      //          // dimensional residual. Internally, the cost function stores the observed
//      //          // image location and compares the reprojection against the observation.
//      //          ceres::CostFunction* cost_function =
//      //                  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_2::Create(obs_it.second.x);
//      //    //        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);


//                //#pragma omp critical
//                //      {

//  //              std::cout <<"out 1" << std::endl;
//                const Mat3 Ro = sfm_data_.poses[view->id_pose].rotation();
//                const Vec3 t = sfm_data_.poses[view->id_pose].translation();
//                double angleAxis[3];
//                ceres::RotationMatrixToAngleAxis((const double*)Ro.data(), angleAxis);

//  //              std::cout <<"out 2" << std::endl;
//                    const double * cam_R = angleAxis;
//                    const Vec3 cam_t = t;

//                    double pos_3dpoint[3];
//  //                  std::cout <<"out 2" << std::endl;
//                    pos_3dpoint[0] = structure_landmark_it.second.X(0);
//                    pos_3dpoint[1] = structure_landmark_it.second.X(1);
//                    pos_3dpoint[2] = structure_landmark_it.second.X(2);

//                    double pos_proj[3];
//                    double X[3];
//                    X[0] = pos_3dpoint[0];
//                    X[1] = pos_3dpoint[1];
//                    X[2] = 0.0;
//                    // Rotate the point according the camera rotation
//                    ceres::AngleAxisRotatePoint(cam_R, X, pos_proj);

//  //                  std::cout <<"out 3" << std::endl;
//                    // Apply the camera translation
//                    pos_proj[0] += cam_t(0);
//                    pos_proj[1] += cam_t(1);
//                    pos_proj[2] += cam_t(2);

//                    //prepare l
//                    double R[9];
//                    ceres::AngleAxisToRotationMatrix(cam_R, R);

//  //                  std::cout <<"out 4" << std::endl;
//                    std::vector<double> it_intrinsics = sfm_data_.intrinsics[view->id_intrinsic]->getParams();

//                    const double f = it_intrinsics[0]; //focal
//                    const double u = it_intrinsics[1]; //principal_point_x
//                    const double v = it_intrinsics[2]; //principal_point_y

//  //                  std::cout <<"out 5" << std::endl;
//                    // Transform the point from homogeneous to euclidean (undistorted point)
//                    const double x_i = u + f * (pos_proj[0] / pos_proj[2]);
//                    const double y_i = v + f * (pos_proj[1] / pos_proj[2]);

//                    const double P1 = f*R[0] + u*R[2];
//                    const double P2 = f*R[3] + u*R[5];
//                    const double P3 = f*R[6] + u*R[8];
//                    const double P4 = f*cam_t[0] + u*cam_t[2];
//                    const double P5 = f*R[1] + v*R[2];
//                    const double P6 = f*R[4] + v*R[5];
//                    const double P7 = f*R[7] + v*R[8];
//                    const double P8 = f*cam_t[1] + v*cam_t[2];
//                    const double P9 = R[2];
//                    const double P10 = R[5];
//                    const double P11 = R[8];
//                    const double P12 = cam_t[2];


//  //                  std::cout <<"out 6" << std::endl;

//                    const double xu = obs_it.second.x(0);
//  //                  std::cout <<"out 6 1" << std::endl;
//                    const double yu = obs_it.second.x(1);
//  //                  std::cout <<"out 6 2" << std::endl;

//                    const double l0 = P7 - P11*yu;
//  //                  std::cout <<"out 6 3" << std::endl;
//                    const double l1 = P11*xu - P3;
//  //                  std::cout <<"out 6 4" << std::endl;
//                    const double l2 = P3*yu - P7*xu;


//  //                  std::cout <<"out 7" << std::endl;

//                    double out_residuals = (x_i*l0 + y_i*l1 + l2)/sqrt(l0*l0 + l1*l1);
//                //    out_residuals[1] = T(0.0);
//                //    out_residuals[0] = abs(x_i*l0 + y_i*l1 + l2)/sqrt(l0*l0 + l1*l1);

//  //                  std::cout <<"out 8" << std::endl;
//                    if(out_residuals > 0.1)
//                    {
//                        out << "cam_R : " << cam_R[0] << " " << cam_R[1] << " " << cam_R[2] << std::endl;
//                        out << "R : " << R[0] << " " << R[1] << " " << R[2]
//                                            << R[3] << " " << R[4] << " " << R[5]
//                                            << R[6] << " " << R[7] << " " << R[8] << std::endl;
//                        out << "cam_t : " << cam_t(0) << " " << cam_t(1) << " " << cam_t(2) << std::endl;
//                        out << "f : " << f << " u : " << u << " v : " << v << std::endl;
//                        out << "xu : " << xu << " yu : " << yu << std::endl;
//                        out << "pos_3dpoint[0]: " << pos_3dpoint[0] << " pos_3dpoint[1]: " << pos_3dpoint[1] << " pos_3dpoint[2]: " << pos_3dpoint[2] << std::endl;
//                        out << "out_residuals[0]: "  << out_residuals << std::endl;
//                        out << "****-----------------****" << std::endl;
//                    }



//              }


//            }
//        }

//        out.close();
    }


    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust_water2_c1
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, //ReconstructionEngine::intrinsic_refinement_options_, //Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_1_refine_RT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    // - refine only Structure and Rotations & translations
//    b_BA_Status = bundle_adjustment_obj.Adjust_water2_c1
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_ROTATION,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    if (b_BA_Status && !sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_2_refine_RT_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

//    // - refine only Structure and Rotations & translations
//    b_BA_Status = bundle_adjustment_obj.Adjust_water2_c1
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    if (b_BA_Status && !sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_3_refine_RT_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }
  }
  std::cout << "water 2 C1" << std::endl;

  {

//      std::ofstream out;
//      out.open("/home/guang/data/testWaterStep/imgs/m_out2/out_1.txt");
//      std::cout <<"out ****-----------------****" << std::endl;
//      {
//          for (auto & structure_landmark_it : sfm_data_.structure)
//          {
//            const Observations & obs = structure_landmark_it.second.obs;

//            for (const auto & obs_it : obs)
//            {
//              // Build the residual block corresponding to the track observation:
//              const View * view = sfm_data_.views.at(obs_it.first).get();

//    //          // Each Residual block takes a point and a camera as input and outputs a 2
//    //          // dimensional residual. Internally, the cost function stores the observed
//    //          // image location and compares the reprojection against the observation.
//    //          ceres::CostFunction* cost_function =
//    //                  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_2::Create(obs_it.second.x);
//    //    //        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);


//              //#pragma omp critical
//              //      {

////              std::cout <<"out 1" << std::endl;
//              const Mat3 Ro = sfm_data_.poses[view->id_pose].rotation();
//              const Vec3 t = sfm_data_.poses[view->id_pose].translation();
//              double angleAxis[3];
//              ceres::RotationMatrixToAngleAxis((const double*)Ro.data(), angleAxis);

////              std::cout <<"out 2" << std::endl;
//                  const double * cam_R = angleAxis;
//                  const Vec3 cam_t = t;

//                  double pos_3dpoint[3];
////                  std::cout <<"out 2" << std::endl;
//                  pos_3dpoint[0] = structure_landmark_it.second.X(0);
//                  pos_3dpoint[1] = structure_landmark_it.second.X(1);
//                  pos_3dpoint[2] = structure_landmark_it.second.X(2);

//                  double pos_proj[3];
//                  double X[3];
//                  X[0] = pos_3dpoint[0];
//                  X[1] = pos_3dpoint[1];
//                  X[2] = 0.0;
//                  // Rotate the point according the camera rotation
//                  ceres::AngleAxisRotatePoint(cam_R, X, pos_proj);

////                  std::cout <<"out 3" << std::endl;
//                  // Apply the camera translation
//                  pos_proj[0] += cam_t(0);
//                  pos_proj[1] += cam_t(1);
//                  pos_proj[2] += cam_t(2);

//                  //prepare l
//                  double R[9];
//                  ceres::AngleAxisToRotationMatrix(cam_R, R);

////                  std::cout <<"out 4" << std::endl;
//                  std::vector<double> it_intrinsics = sfm_data_.intrinsics[view->id_intrinsic]->getParams();

//                  const double f = it_intrinsics[0]; //focal
//                  const double u = it_intrinsics[1]; //principal_point_x
//                  const double v = it_intrinsics[2]; //principal_point_y

////                  std::cout <<"out 5" << std::endl;
//                  // Transform the point from homogeneous to euclidean (undistorted point)
//                  const double x_i = u + f * (pos_proj[0] / pos_proj[2]);
//                  const double y_i = v + f * (pos_proj[1] / pos_proj[2]);

//                  const double P1 = f*R[0] + u*R[2];
//                  const double P2 = f*R[3] + u*R[5];
//                  const double P3 = f*R[6] + u*R[8];
//                  const double P4 = f*cam_t[0] + u*cam_t[2];
//                  const double P5 = f*R[1] + v*R[2];
//                  const double P6 = f*R[4] + v*R[5];
//                  const double P7 = f*R[7] + v*R[8];
//                  const double P8 = f*cam_t[1] + v*cam_t[2];
//                  const double P9 = R[2];
//                  const double P10 = R[5];
//                  const double P11 = R[8];
//                  const double P12 = cam_t[2];


////                  std::cout <<"out 6" << std::endl;

//                  const double xu = obs_it.second.x(0);
////                  std::cout <<"out 6 1" << std::endl;
//                  const double yu = obs_it.second.x(1);
////                  std::cout <<"out 6 2" << std::endl;

//                  const double l0 = P7 - P11*yu;
////                  std::cout <<"out 6 3" << std::endl;
//                  const double l1 = P11*xu - P3;
////                  std::cout <<"out 6 4" << std::endl;
//                  const double l2 = P3*yu - P7*xu;


////                  std::cout <<"out 7" << std::endl;

//                  double out_residuals = (x_i*l0 + y_i*l1 + l2)/sqrt(l0*l0 + l1*l1);
//              //    out_residuals[1] = T(0.0);
//              //    out_residuals[0] = abs(x_i*l0 + y_i*l1 + l2)/sqrt(l0*l0 + l1*l1);

////                  std::cout <<"out 8" << std::endl;
//                  if(out_residuals > 0.1)
//                  {
//                      out << "cam_R : " << cam_R[0] << " " << cam_R[1] << " " << cam_R[2] << std::endl;
//                      out << "R : " << R[0] << " " << R[1] << " " << R[2]
//                                          << R[3] << " " << R[4] << " " << R[5]
//                                          << R[6] << " " << R[7] << " " << R[8] << std::endl;
//                      out << "cam_t : " << cam_t(0) << " " << cam_t(1) << " " << cam_t(2) << std::endl;
//                      out << "f : " << f << " u : " << u << " v : " << v << std::endl;
//                      out << "xu : " << xu << " yu : " << yu << std::endl;
//                      out << "pos_3dpoint[0]: " << pos_3dpoint[0] << " pos_3dpoint[1]: " << pos_3dpoint[1] << " pos_3dpoint[2]: " << pos_3dpoint[2] << std::endl;
//                      out << "out_residuals[0]: "  << out_residuals << std::endl;
//                      out << "****-----------------****" << std::endl;
//                  }



//            }


//          }
//      }

//      out.close();
  }
  // - refine only Structure and Rotations & translations
//  b_BA_Status = bundle_adjustment_obj.Adjust_water2_c1
//    (
//      sfm_data_,
//      Optimize_Options(
//        ReconstructionEngine::intrinsic_refinement_options_,//Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_ALL,
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );
//  if (b_BA_Status && !sLogging_file_.empty())
//  {
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_RT_Xi", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
//  }

  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust_water2_c1
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    //if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  Optimize_Options ba_refine_options;
  if (ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
      ba_refine_options = Optimize_Options(
        ReconstructionEngine::intrinsic_refinement_options_,
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }else{
      ba_refine_options = Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }

  b_BA_Status = bundle_adjustment_obj.Adjust_water2_c1(sfm_data_, ba_refine_options);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  return b_BA_Status;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_init_water3_c1()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );

//  // Remove outliers (max_angle, residual error)
//  const size_t pointcount_initial1 = sfm_data_.structure.size();
//  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//  const size_t pointcount_pixelresidual_filter1 = sfm_data_.structure.size();
//  RemoveOutliers_AngleError(sfm_data_, 2.0);
//  const size_t pointcount_angular_filter1 = sfm_data_.structure.size();
//  std::cout << "Outlier removal (remaining #points):\n"
//    << "\t initial structure size #3DPoints: " << pointcount_initial1 << "\n"
//    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter1 << "\n"
//    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter1 << std::endl;


//  bool b_BA_Status// = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::NONE, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );


//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  if (b_BA_Status)
  {
//    if (!sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

//    b_BA_Status = bundle_adjustment_obj.Adjust
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
//    {
//        // Remove outliers (max_angle, residual error)
//        const size_t pointcount_initial = sfm_data_.structure.size();
//        RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//        const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
//        RemoveOutliers_AngleError(sfm_data_, 2.0);
//        const size_t pointcount_angular_filter = sfm_data_.structure.size();
//        std::cout << "Outlier removal (remaining #points):\n"
//          << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
//          << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
//          << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;
//        Save(sfm_data_,
//          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_000_refine_T_Xi", "ply"),
//          ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

    {
        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_0_refine_RT_Xi", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));

        SfM_Data_Structure_Computation_Blind structure_estimator(true);
        structure_estimator.triangulate(sfm_data_);

        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_1_refine_RT_Xi", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust_water3_c1
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );

    {
        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_0_refine_RT_Xi", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));

        SfM_Data_Structure_Computation_Blind structure_estimator(true);
        structure_estimator.triangulate(sfm_data_);

        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_1_refine_RT_Xi", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));
    }


    b_BA_Status = bundle_adjustment_obj.Adjust_water3_c1
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    {
        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_0_refine_RT_Xi", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));

        SfM_Data_Structure_Computation_Blind structure_estimator(true);
        structure_estimator.triangulate(sfm_data_);

        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_1_refine_RT_Xi", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));
    }
//    if (b_BA_Status && !sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

    b_BA_Status = bundle_adjustment_obj.Adjust_water3_c1
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    {
        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_0_refine_RT_Xi", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));

        SfM_Data_Structure_Computation_Blind structure_estimator(true);
        structure_estimator.triangulate(sfm_data_);

        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_1_refine_RT_Xi", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));
    }
//    if (b_BA_Status && !sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }
  }
//  std::cout << "water 3 C1 " << std::endl;

//  // - refine only Structure and Rotations & translations
//  b_BA_Status = bundle_adjustment_obj.Adjust_water3_c1
//    (
//      sfm_data_,
//      Optimize_Options(
//        ReconstructionEngine::intrinsic_refinement_options_,//Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_ALL,
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );
//  if (b_BA_Status && !sLogging_file_.empty())
//  {
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_RT_Xi", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
//  }

  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust_water3_c1
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  Optimize_Options ba_refine_options;
  if (ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
      ba_refine_options = Optimize_Options(
        ReconstructionEngine::intrinsic_refinement_options_,
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }else{
      ba_refine_options = Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }

  b_BA_Status = bundle_adjustment_obj.Adjust(sfm_data_, ba_refine_options);//_water3_c1(sfm_data_, ba_refine_options);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  return b_BA_Status;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_init_water_c4()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = true;
//          = bundle_adjustment_obj.Adjust


  if (b_BA_Status)
  {
//    if (!sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

//    b_BA_Status = bundle_adjustment_obj.Adjust
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_TRANSLATION,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
////    {
////        // Remove outliers (max_angle, residual error)
////        const size_t pointcount_initial = sfm_data_.structure.size();
////        RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
////        const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
////        RemoveOutliers_AngleError(sfm_data_, 2.0);
////        const size_t pointcount_angular_filter = sfm_data_.structure.size();
////        std::cout << "Outlier removal (remaining #points):\n"
////          << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
////          << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
////          << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;
////    }


//      b_BA_Status = bundle_adjustment_obj.Adjust//_water_c4_S
//        (
//          sfm_data_,
//          Optimize_Options(
//            ReconstructionEngine::intrinsic_refinement_options_,//Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//            Extrinsic_Parameter_Type::ADJUST_ALL,
//            Structure_Parameter_Type::ADJUST_ALL,
//            Control_Point_Parameter(),
//            this->b_use_motion_prior_)
//        );

//      b_BA_Status = bundle_adjustment_obj.Adjust//_water_c4_S
//        (
//          sfm_data_,
//          Optimize_Options(
//            ReconstructionEngine::intrinsic_refinement_options_, // Intrinsics are held as constant
//            Extrinsic_Parameter_Type::ADJUST_ALL,
//            Structure_Parameter_Type::ADJUST_ALL,
//            Control_Point_Parameter(),
//            this->b_use_motion_prior_)
//        );
      const size_t pointcount_initial = sfm_data_.structure.size();
      RemoveOutliers_PixelResidualError(sfm_data_, 8.0);//4.0);
      const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
      RemoveOutliers_AngleError(sfm_data_, 5.0);//2.0);
      const size_t pointcount_angular_filter = sfm_data_.structure.size();
      std::cout << "Outlier removal (remaining #points):\n"
        << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
        << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
        << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

//      // Check that poses & intrinsic cover some measures (after outlier removal)
//      const IndexT minPointPerPose = 12; // 6 min
//      const IndexT minTrackLength = 3; // 2 min
//      if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//      {
//        // TODO: must ensure that track graph is producing a single connected component

//        const size_t pointcount_cleaning = sfm_data_.structure.size();
//        std::cout << "Point_cloud cleaning:\n"
//          << "\t #3DPoints: " << pointcount_cleaning << "\n";
//      }

      b_BA_Status = bundle_adjustment_obj.Adjust//_water_c4_S
        (
          sfm_data_,
          Optimize_Options(
            Intrinsic_Parameter_Type::NONE,//ReconstructionEngine::intrinsic_refinement_options_,//Intrinsic_Parameter_Type::NONE,//ReconstructionEngine::intrinsic_refinement_options_, // Intrinsics are held as constant
            Extrinsic_Parameter_Type::ADJUST_ALL,
            Structure_Parameter_Type::ADJUST_ALL,
            Control_Point_Parameter(),
            this->b_use_motion_prior_)
        );

//    b_BA_Status = bundle_adjustment_obj.Adjust//_water_c4_S
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE,//ReconstructionEngine::intrinsic_refinement_options_,//Intrinsic_Parameter_Type::NONE,//ReconstructionEngine::intrinsic_refinement_options_, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_1_refine_RT_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));


      std::cout << "error 6" <<std::endl;
      {
//          SfM_Data newdata = sfm_data_;
////          newdata.intrinsics = sfmEngine.Get_SfM_Data().intrinsics;
//          Mat3 changeR;
//          Vec3 changet;
//          float s;
//          {
//              s = 1.0;
////              changeR << 0.999954402447, 0.00858933758, 0.004174199887,
////                         0.00193876868, -0.609908938408, 0.792469143867,
////                         0.009352667257, -0.792424976826, -0.60989767313;
////              changet << -0.004134478979, -0.811449706554, 1.615199804306;
////              changeR << 0.999785900116, 0.014877740294, 0.014379678294,
////                         0.0, -0.694968700409, 0.719039976597,
////                         0.020691115409, -0.718886077404, -0.694819927216;
////              changet << -0.000973686634, -1.160322427750, 0.586338579655;
////              changeR << 0.994793236256, -0.018639463931, -0.100195005536,
////                         0.000000000000, -0.983132660389, 0.182893991470,
////                         -0.101914025843, -0.181941702962, -0.978013694286;
////              changet << 0.092311173677, -0.344401419163, 1.838596105576;
//              //pot6
////              changeR << 0.999977469444, 0.004028763156, 0.005373591557,
////                         -0.00019429854, -0.782413601875, 0.622759103775,
////                         0.006713320035, -0.622746109962, -0.78239518404;
////              changet << -0.003047744278, -2.917632579803, 1.824462413788;
//              //pot7
////              changeR << 1.0, -0.000014351798, -0.000032206939,
////                         0.0, -0.913415312767, 0.407028824091,
////                         -0.000035259905, -0.407028824091, -0.913415312767;
////              changet << 0.000033166896, -0.462338685989, 1.979631304741;
//              //pot11
////              changeR << 1.0, 0.0, 0.0,
////                         0.0, 0.970297813416, -0.24191364646,
////                         0.0, 0.24191364646, 0.970297813416;
////              changet << 0.0, 0.132218018174, 0.055407147855;
//              //pot11
////              changeR << 0.999147057533, -0.009465177543, -0.040193554014,
////                         0.0, 0.973374664783, -0.229219928384,
////                         0.041292995214, 0.229024425149, 0.972544431686;
////              changet << 0.041424274445, 0.235334098339, 0.020509306341;


//              double angleAxis[3] = {0,45.0/180*3.14,0};
//              ceres::AngleAxisToRotationMatrix(angleAxis, (double*)changeR.data());
//              changet << 0.0, 0.0, 0.0;

////              double angleAxis[3];
////              ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
////              // angleAxis + translation
////              map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};



//              std::cout << "detR : " << changeR.determinant() << std::endl;

//              //
//              Mat3 changeZ;
//              changeZ << 1, 0, 0,
//                         0, 1, 0,
//                         0, 0, 1;
////              changeZ << -1, 0, 0,
////                         0, 1, 0,
////                         0, 0, -1;
//              //
//              changeR = changeZ * changeR;
//              //
//              changeR = changeR/s;
//              changet = changet/s;

//              //s = 0.98445;//0.30558845;//1.0;//3.1908;
//          }

//          for (const auto & pose_it : newdata.poses)
//          {
//            const IndexT indexPose = pose_it.first;

//            const Pose3 & pose = pose_it.second;
//            const Mat3 R = pose.rotation();
//            const Vec3 C = pose.center();

//            Mat3 Rnew = R * changeR.inverse();
////            Vec3 Cnew = changet + Rnew.inverse() * R * C;
//            Vec3 Cnew = changeR * s * C + changet;

//            newdata.poses[indexPose] = Pose3(Rnew, Cnew);

//          }

//          newdata.structure = newdata.structure;
//          //set all X
//          for(Landmarks::iterator itX = newdata.structure.begin();
//              itX != newdata.structure.end(); ++itX)
//          {
//              const Vec3 & X = itX->second.X;

//              Vec3 newX = changeR * s * X + changet;
////                  Vec3 newX = R_gc * X + t_gc;
////                  Vec3 newX = R0.inverse() * X;
//              //
//              itX->second.X = newX;
////              for(Observations::iterator itObs = itX->second.obs.begin();
////                  itObs != itX->second.obs.end(); ++itObs)
////              {
////                  const View * view = newdata.views.at(itObs->first).get();
////                  const openMVG::cameras::Pinhole_Intrinsic * cam = dynamic_cast<const openMVG::cameras::Pinhole_Intrinsic*>(newdata.intrinsics[view->id_intrinsic].get());
////                  Vec2 ud_x = cam->get_ud_pixel(itObs->second.x);
////                  itObs->second.x = ud_x;
////              }
//          }

//          //-- Export to disk computed scene (data & visualizable results)
//          std::cout << "...Export SfM_Data to disk." << std::endl;

//          Save(newdata,
//            stlplus::create_filespec(newdata.s_root_path, "sfm_data_", ".json"),
//            ESfM_Data(ALL));

//          Save(newdata,
//            stlplus::create_filespec(newdata.s_root_path, "sfm_data_", ".ply"),
//            ESfM_Data(ALL));

//          Save(newdata,
//            stlplus::create_filespec(newdata.s_root_path, "INTRINSICS", ".json"),
//            ESfM_Data(INTRINSICS));
//          std::cout << "error 7" <<std::endl;
//          sfm_data_ = newdata;

//          return true;

      }

      {
          SfM_Data newdata = sfm_data_;
          std::cout << "newdata " << newdata.poses.size() << std::endl;
          //set all X
          for(Landmarks::iterator itX = newdata.structure.begin();
              itX != newdata.structure.end(); ++itX)
          {
              for(Observations::iterator itObs = itX->second.obs.begin();
                  itObs != itX->second.obs.end(); ++itObs)
              {
                  const View * view = newdata.views.at(itObs->first).get();
                  const openMVG::cameras::Pinhole_Intrinsic * cam = dynamic_cast<const openMVG::cameras::Pinhole_Intrinsic*>(newdata.intrinsics[view->id_intrinsic].get());
                  Vec2 ud_x = cam->get_ud_pixel(itObs->second.x);
                  itObs->second.x = ud_x;
              }
          }
          sfm_data_ = newdata;
      }




//      return true;
//      b_BA_Status = bundle_adjustment_obj.Adjust//_water_c4_S
//        (
//          sfm_data_,
//          Optimize_Options(
//            Intrinsic_Parameter_Type::NONE,//ReconstructionEngine::intrinsic_refinement_options_, // Intrinsics are held as constant
//            Extrinsic_Parameter_Type::ADJUST_ALL,
//            Structure_Parameter_Type::ADJUST_ALL,
//            Control_Point_Parameter(),
//            this->b_use_motion_prior_)
//        );
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_2_refine_RT_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));

    // - refine only Structure and Rotations & translations
//    b_BA_Status = bundle_adjustment_obj.Adjust_water_c1_S
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::NONE,//Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          true)
////          this->b_use_motion_prior_)
//      );

    {
        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_0_refine_RT_Xi", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));
        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_0_refine_RT_Xi", "json"),
          ESfM_Data(EXTRINSICS));

        b_BA_Status = bundle_adjustment_obj.Adjust//_water_c1_pi
          (
            sfm_data_,
            Optimize_Options(
              Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
              Extrinsic_Parameter_Type::ADJUST_ALL,
              Structure_Parameter_Type::ADJUST_ALL,
              Control_Point_Parameter(),
              false)//true)//false)//
    //          this->b_use_motion_prior_)
          );
        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_1_refine_RT_Xi", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));
        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_1_refine_RT_Xi", "json"),
          ESfM_Data(EXTRINSICS));

//        b_BA_Status = bundle_adjustment_obj.Adjust_water_c4_S
//          (
//            sfm_data_,
//            Optimize_Options(
//              ReconstructionEngine::intrinsic_refinement_options_,//Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//              Extrinsic_Parameter_Type::ADJUST_ALL,
//              Structure_Parameter_Type::ADJUST_ALL,
//              Control_Point_Parameter(),
//              true)
//    //          this->b_use_motion_prior_)
//          );

//        Save(sfm_data_,
//          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_2_refine_RT_Xi", "ply"),
//          ESfM_Data(EXTRINSICS | STRUCTURE));
    }
    return true;

    b_BA_Status = bundle_adjustment_obj.Adjust_water_c1_S
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
//    {
//        Save(sfm_data_,
//          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_0_refine_RT_Xi", "ply"),
//          ESfM_Data(EXTRINSICS | STRUCTURE));

//        SfM_Data_Structure_Computation_Blind structure_estimator(true);
//        structure_estimator.triangulate(sfm_data_);

//        Save(sfm_data_,
//          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_1_refine_RT_Xi", "ply"),
//          ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

//    b_BA_Status = bundle_adjustment_obj.Adjust//_water5_c1
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::NONE,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    {
//        Save(sfm_data_,
//          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_0_refine_RT_Xi", "ply"),
//          ESfM_Data(EXTRINSICS | STRUCTURE));

//        SfM_Data_Structure_Computation_Blind structure_estimator(true);
//        structure_estimator.triangulate(sfm_data_);

//        Save(sfm_data_,
//          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_1_refine_RT_Xi", "ply"),
//          ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

//    if (b_BA_Status && !sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }
  }
  std::cout << "water 3 C1 " << std::endl;


//  // - refine only Structure and Rotations & translations
//  b_BA_Status = bundle_adjustment_obj.Adjust_water3_c1
//    (
//      sfm_data_,
//      Optimize_Options(
//        ReconstructionEngine::intrinsic_refinement_options_,//Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_ALL,
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );
//  if (b_BA_Status && !sLogging_file_.empty())
//  {
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_RT_Xi", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
//  }

  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust_water_c4_S
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  // Check that poses & intrinsic cover some measures (after outlier removal)
  const IndexT minPointPerPose = 12; // 6 min
  const IndexT minTrackLength = 3; // 2 min
  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
  {
    // TODO: must ensure that track graph is producing a single connected component

    const size_t pointcount_cleaning = sfm_data_.structure.size();
    std::cout << "Point_cloud cleaning:\n"
      << "\t #3DPoints: " << pointcount_cleaning << "\n";
  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  Optimize_Options ba_refine_options;
  if (ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
      ba_refine_options = Optimize_Options(
        ReconstructionEngine::intrinsic_refinement_options_,
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }else{
      ba_refine_options = Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }

  b_BA_Status = bundle_adjustment_obj.Adjust_water_c4_S(sfm_data_, ba_refine_options);//_water5_c1(sfm_data_, ba_refine_options);//
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  return b_BA_Status;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_init_water5_c1()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = true;
//          = bundle_adjustment_obj.Adjust


  if (b_BA_Status)
  {
//    if (!sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

//    b_BA_Status = bundle_adjustment_obj.Adjust
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
//    {
//        // Remove outliers (max_angle, residual error)
//        const size_t pointcount_initial = sfm_data_.structure.size();
//        RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//        const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
//        RemoveOutliers_AngleError(sfm_data_, 2.0);
//        const size_t pointcount_angular_filter = sfm_data_.structure.size();
//        std::cout << "Outlier removal (remaining #points):\n"
//          << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
//          << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
//          << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;
//    }

    {
        SfM_Data newdata = sfm_data_;
        newdata.intrinsics = sfm_data_.intrinsics;
        Mat3 changeR;
        Vec3 changet;
        float s;
        {

            s = 0.985698;
            changeR << 0.985211133957, 0.001476573758, 0.030937716365,
                       0.004451564513, 0.967604577541, -0.187941178679,
                       -0.030651364475, 0.187988087535, 0.967120110989;
            changet <<  -0.039469640702,  0.225849807262,  -0.146059274673;



            Mat3 changeZ;
            changeZ << -1, 0, 0,
                       0, 1, 0,
                       0, 0, -1;


            changeR = changeZ * changeR;

            changeR = changeR/s;
            changet = changet/s;

            //s = 0.98445;//0.30558845;//1.0;//3.1908;
        }

        for (const auto & pose_it : sfm_data_.poses)
        {
          const IndexT indexPose = pose_it.first;

          const Pose3 & pose = pose_it.second;
          const Mat3 R = pose.rotation();
          const Vec3 C = pose.center();

          Mat3 Rnew = R * changeR.inverse();
//            Vec3 Cnew = changet + Rnew.inverse() * R * C;
          Vec3 Cnew = changeR * s * C + changet;

          newdata.poses[indexPose] = Pose3(Rnew, Cnew);

        }

        newdata.structure = sfm_data_.structure;
        //set all X
        for(Landmarks::iterator itX = newdata.structure.begin();
            itX != newdata.structure.end(); ++itX)
        {
            const Vec3 & X = itX->second.X;

            Vec3 newX = changeR * s * X + changet;
//                  Vec3 newX = R_gc * X + t_gc;
//                  Vec3 newX = R0.inverse() * X;
            itX->second.X = newX;
        }

        sfm_data_ = newdata;

        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_RT_Xi", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));

    }
    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust_water_c1_S
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );

    {
        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_0_refine_RT_Xi", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));

//        SfM_Data_Structure_Computation_Blind structure_estimator(true);
//        structure_estimator.triangulate(sfm_data_);

        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_1_refine_RT_Xi", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));
    }
//    return false;

//    b_BA_Status = bundle_adjustment_obj.Adjust_water_c1_S
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    {
//        Save(sfm_data_,
//          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_0_refine_RT_Xi", "ply"),
//          ESfM_Data(EXTRINSICS | STRUCTURE));

//        SfM_Data_Structure_Computation_Blind structure_estimator(true);
//        structure_estimator.triangulate(sfm_data_);

//        Save(sfm_data_,
//          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_1_refine_RT_Xi", "ply"),
//          ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

//    b_BA_Status = bundle_adjustment_obj.Adjust//_water5_c1
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::NONE,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    {
//        Save(sfm_data_,
//          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_0_refine_RT_Xi", "ply"),
//          ESfM_Data(EXTRINSICS | STRUCTURE));

//        SfM_Data_Structure_Computation_Blind structure_estimator(true);
//        structure_estimator.triangulate(sfm_data_);

//        Save(sfm_data_,
//          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_1_refine_RT_Xi", "ply"),
//          ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

//    if (b_BA_Status && !sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }
  }
  std::cout << "water 3 C1 " << std::endl;


//  // - refine only Structure and Rotations & translations
//  b_BA_Status = bundle_adjustment_obj.Adjust_water3_c1
//    (
//      sfm_data_,
//      Optimize_Options(
//        ReconstructionEngine::intrinsic_refinement_options_,//Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_ALL,
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );
//  if (b_BA_Status && !sLogging_file_.empty())
//  {
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_RT_Xi", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
//  }

  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust_water_c1_S
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  Optimize_Options ba_refine_options;
  if (ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
      ba_refine_options = Optimize_Options(
        ReconstructionEngine::intrinsic_refinement_options_,
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }else{
      ba_refine_options = Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }

  b_BA_Status = bundle_adjustment_obj.Adjust_water_c1_S(sfm_data_, ba_refine_options);//_water5_c1(sfm_data_, ba_refine_options);//
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  return b_BA_Status;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_init_water5_c4()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = true;
//          = bundle_adjustment_obj.Adjust


  if (b_BA_Status)
  {
//    if (!sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

//    b_BA_Status = bundle_adjustment_obj.Adjust
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
//    {
//        // Remove outliers (max_angle, residual error)
//        const size_t pointcount_initial = sfm_data_.structure.size();
//        RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//        const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
//        RemoveOutliers_AngleError(sfm_data_, 2.0);
//        const size_t pointcount_angular_filter = sfm_data_.structure.size();
//        std::cout << "Outlier removal (remaining #points):\n"
//          << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
//          << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
//          << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;
//    }

    {
        SfM_Data newdata = sfm_data_;
        newdata.intrinsics = sfm_data_.intrinsics;
        Mat3 changeR;
        Vec3 changet;
        float s;
        {

            s = 0.985698;
            changeR << 0.985211133957, 0.001476573758, 0.030937716365,
                       0.004451564513, 0.967604577541, -0.187941178679,
                       -0.030651364475, 0.187988087535, 0.967120110989;
            changet <<  -0.039469640702,  0.225849807262,  -0.146059274673;



            Mat3 changeZ;
            changeZ << -1, 0, 0,
                       0, 1, 0,
                       0, 0, -1;


            changeR = changeZ * changeR;

            changeR = changeR/s;
            changet = changet/s;

            //s = 0.98445;//0.30558845;//1.0;//3.1908;
        }

        for (const auto & pose_it : sfm_data_.poses)
        {
          const IndexT indexPose = pose_it.first;

          const Pose3 & pose = pose_it.second;
          const Mat3 R = pose.rotation();
          const Vec3 C = pose.center();

          Mat3 Rnew = R * changeR.inverse();
//            Vec3 Cnew = changet + Rnew.inverse() * R * C;
          Vec3 Cnew = changeR * s * C + changet;

          newdata.poses[indexPose] = Pose3(Rnew, Cnew);

        }

        newdata.structure = sfm_data_.structure;
        //set all X
        for(Landmarks::iterator itX = newdata.structure.begin();
            itX != newdata.structure.end(); ++itX)
        {
            const Vec3 & X = itX->second.X;

            Vec3 newX = changeR * s * X + changet;
//                  Vec3 newX = R_gc * X + t_gc;
//                  Vec3 newX = R0.inverse() * X;
            itX->second.X = newX;
        }

        sfm_data_ = newdata;

        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_RT_Xi", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));

    }
    // - refine only Structure and Rotations & translations
      b_BA_Status = bundle_adjustment_obj.Adjust
        (
          sfm_data_,
          Optimize_Options(
            Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
            Extrinsic_Parameter_Type::ADJUST_ALL,
            Structure_Parameter_Type::ADJUST_ALL,
            Control_Point_Parameter(),
            this->b_use_motion_prior_)
        );

    b_BA_Status = bundle_adjustment_obj.Adjust_water_c4_S
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );

    {
        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_0_refine_RT_Xi", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));

//        SfM_Data_Structure_Computation_Blind structure_estimator(true);
//        structure_estimator.triangulate(sfm_data_);

        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_1_refine_RT_Xi", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));
    }
//    return false;

    b_BA_Status = bundle_adjustment_obj.Adjust_water_c4_S
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
//    {
//        Save(sfm_data_,
//          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_0_refine_RT_Xi", "ply"),
//          ESfM_Data(EXTRINSICS | STRUCTURE));

//        SfM_Data_Structure_Computation_Blind structure_estimator(true);
//        structure_estimator.triangulate(sfm_data_);

//        Save(sfm_data_,
//          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_1_refine_RT_Xi", "ply"),
//          ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

//    b_BA_Status = bundle_adjustment_obj.Adjust//_water5_c1
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::NONE,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    {
//        Save(sfm_data_,
//          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_0_refine_RT_Xi", "ply"),
//          ESfM_Data(EXTRINSICS | STRUCTURE));

//        SfM_Data_Structure_Computation_Blind structure_estimator(true);
//        structure_estimator.triangulate(sfm_data_);

//        Save(sfm_data_,
//          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_1_refine_RT_Xi", "ply"),
//          ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

//    if (b_BA_Status && !sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }
  }
  std::cout << "water 3 C1 " << std::endl;


//  // - refine only Structure and Rotations & translations
//  b_BA_Status = bundle_adjustment_obj.Adjust_water3_c1
//    (
//      sfm_data_,
//      Optimize_Options(
//        ReconstructionEngine::intrinsic_refinement_options_,//Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_ALL,
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );
//  if (b_BA_Status && !sLogging_file_.empty())
//  {
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_RT_Xi", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
//  }

  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust_water_c4_S
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  Optimize_Options ba_refine_options;
  if (ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
      ba_refine_options = Optimize_Options(
        ReconstructionEngine::intrinsic_refinement_options_,
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }else{
      ba_refine_options = Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }

  b_BA_Status = bundle_adjustment_obj.Adjust_water_c4_S(sfm_data_, ba_refine_options);//_water5_c1(sfm_data_, ba_refine_options);//
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  return b_BA_Status;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_init_water6_c1()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = true;
//          = bundle_adjustment_obj.Adjust


  if (b_BA_Status)
  {
//    if (!sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

//    b_BA_Status = bundle_adjustment_obj.Adjust
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
//    {
//        // Remove outliers (max_angle, residual error)
//        const size_t pointcount_initial = sfm_data_.structure.size();
//        RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//        const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
//        RemoveOutliers_AngleError(sfm_data_, 2.0);
//        const size_t pointcount_angular_filter = sfm_data_.structure.size();
//        std::cout << "Outlier removal (remaining #points):\n"
//          << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
//          << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
//          << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;
//    }

//      const Mat3 R210239 = sfm_data_.poses[210239].rotation();
//      Vec3 C210239 = sfm_data_.poses[210239].center();
//      C210239(0) += 10;
//      sfm_data_.poses[210239] = Pose3(R210239, C210239);

//       const Mat3 R210237 = sfm_data_.poses[210237].rotation();
//       Vec3 C210237 = sfm_data_.poses[210237].center();
//       C210237(0) += 10;
//       sfm_data_.poses[210237] = Pose3(R210237, C210237);

    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust_water6_c1
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );

    {
        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_0_refine_RT_Xi", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));

        SfM_Data_Structure_Computation_Blind structure_estimator(true);
        structure_estimator.triangulate(sfm_data_);

        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_1_refine_RT_Xi", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    //if (b_BA_Status && !sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_1_refine_RT_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }
//    return false;

//    b_BA_Status = bundle_adjustment_obj.Adjust_water6_c1
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    //if (b_BA_Status && !sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_2_refine_RT_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

//    b_BA_Status = bundle_adjustment_obj.Adjust_water6_c1
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    //if (b_BA_Status && !sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_3_refine_RT_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }
  }
  std::cout << "water 3 C1 " << std::endl;


  {

//      std::ofstream out;
//      out.open("/home/guang/data/testWaterStep/imgs/m_out6/out_1.txt");
//      std::cout <<"out ****-----------------****" << std::endl;
//      {
//          for (auto & structure_landmark_it : sfm_data_.structure)
//          {
//            const Observations & obs = structure_landmark_it.second.obs;

//            for (const auto & obs_it : obs)
//            {
//              // Build the residual block corresponding to the track observation:
//                if(obs_it.first != 210573)
//                {
//                    continue;

//                }
//              const View * view = sfm_data_.views.at(obs_it.first).get();

//    //          // Each Residual block takes a point and a camera as input and outputs a 2
//    //          // dimensional residual. Internally, the cost function stores the observed
//    //          // image location and compares the reprojection against the observation.
//    //          ceres::CostFunction* cost_function =
//    //                  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2_water_C1_2::Create(obs_it.second.x);
//    //    //        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);


//              //#pragma omp critical
//              //      {

////              std::cout <<"out 1" << std::endl;
//              const Mat3 Ro = sfm_data_.poses[view->id_pose].rotation();
//              const Vec3 t = sfm_data_.poses[view->id_pose].translation();
//              double angleAxis[3];
//              ceres::RotationMatrixToAngleAxis((const double*)Ro.data(), angleAxis);

////              std::cout <<"out 2" << std::endl;
//                  const double * cam_R = angleAxis;
//                  const Vec3 cam_t = t;

//                  double pos_3dpoint[3];
////                  std::cout <<"out 2" << std::endl;
//                  pos_3dpoint[0] = structure_landmark_it.second.X(0);
//                  pos_3dpoint[1] = structure_landmark_it.second.X(1);
//                  pos_3dpoint[2] = structure_landmark_it.second.X(2);

//                  double pos_proj[3];
//                  double X[3];
//                  X[0] = pos_3dpoint[0];
//                  X[1] = pos_3dpoint[1];
//                  X[2] = pos_3dpoint[2];//0,0
//                  // Rotate the point according the camera rotation
//                  ceres::AngleAxisRotatePoint(cam_R, X, pos_proj);

////                  std::cout <<"out 3" << std::endl;
//                  // Apply the camera translation
//                  pos_proj[0] += cam_t(0);
//                  pos_proj[1] += cam_t(1);
//                  pos_proj[2] += cam_t(2);

//                  //prepare l
//                  double R[9];
//                  ceres::AngleAxisToRotationMatrix(cam_R, R);

////                  std::cout <<"out 4" << std::endl;
//                  std::vector<double> it_intrinsics = sfm_data_.intrinsics[view->id_intrinsic]->getParams();

//                  const double f = it_intrinsics[0]; //focal
//                  const double u = it_intrinsics[1]; //principal_point_x
//                  const double v = it_intrinsics[2]; //principal_point_y

////                  std::cout <<"out 5" << std::endl;
//                  // Transform the point from homogeneous to euclidean (undistorted point)
//                  const double x_i = u + f * (pos_proj[0] / pos_proj[2]);
//                  const double y_i = v + f * (pos_proj[1] / pos_proj[2]);

//                  const double P1 = f*R[0] + u*R[2];
//                  const double P2 = f*R[3] + u*R[5];
//                  const double P3 = f*R[6] + u*R[8];
//                  const double P4 = f*cam_t[0] + u*cam_t[2];
//                  const double P5 = f*R[1] + v*R[2];
//                  const double P6 = f*R[4] + v*R[5];
//                  const double P7 = f*R[7] + v*R[8];
//                  const double P8 = f*cam_t[1] + v*cam_t[2];
//                  const double P9 = R[2];
//                  const double P10 = R[5];
//                  const double P11 = R[8];
//                  const double P12 = cam_t[2];

//                  const double xu = obs_it.second.x(0);
//                  const double yu = obs_it.second.x(1);

//                  const double nc_x = P3/P11;
//                  const double nc_y = P7/P11;

//                  const double A = yu - nc_y;
//                  const double B = nc_x - xu;
//                  const double C = xu*nc_y - yu*nc_x;

//                  const double D = abs(A*x_i + B*y_i + C) / sqrt(A*A+B*B);
//                  //
//                  const double l_ncxi_x = x_i - nc_x;
//                  const double l_ncxi_y = y_i - nc_y;
//                  const double normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
//              //    const T normal_ncxi = (l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
//                  //
//                  const double l_ncxu_x = xu - nc_x;
//                  const double l_ncxu_y = yu - nc_y;
//                  const double normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);

//                  double out_residuals[2];
//                  out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi);
//                  out_residuals[1] = (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);

//                  {
//                      out << "***----------------***" << std::endl;
//                      out << "R : " << cam_R[0] << " " << cam_R[1] << " " << cam_R[2] << std::endl;
//                      out << "t : " << cam_t[0] << " " << cam_t[1] << " " << cam_t[2] << std::endl;
//                      out << "f : " << f << " u : " << u << " v : " << v << std::endl;
//                      out << "x_i : " << x_i << " y_i : " << y_i << std::endl;
//                      out << "xu : " << xu << " yu : " << yu << std::endl;
//                      out << "nc_x : " << nc_x << " nc_y : " << nc_y << std::endl;
//                      out << "l_ncxi_x : " << l_ncxi_x << " l_ncxi_y : " << l_ncxi_y << std::endl;
//                      out << "l_ncxu_x : " << l_ncxu_x << " l_ncxu_y : " << l_ncxu_y << std::endl;
//                      out << "normal_ncxi : " << normal_ncxi << " normal_ncxu : " << normal_ncxu << std::endl;
//                      out << "out_residuals[0] : " << out_residuals[0] << " out_residuals[1] : " << out_residuals[1] << std::endl;
//                      out << "X : " << X[0] << " " << X[1] << " " << X[2] << std::endl;
//                      out << "      ---------" << std::endl;
//                  }

//                  {

//                      double pos_3dpoint[3];
//    //                  std::cout <<"out 2" << std::endl;
//                      pos_3dpoint[0] = structure_landmark_it.second.X(0);
//                      pos_3dpoint[1] = structure_landmark_it.second.X(1);
//                      pos_3dpoint[2] = structure_landmark_it.second.X(2);

//                      double pos_proj[3];
//                      double X[3];
//                      X[0] = pos_3dpoint[0];
//                      X[1] = pos_3dpoint[1];
//                      X[2] = 0.0;
//                      ceres::AngleAxisRotatePoint(cam_R, X, pos_proj);

//                  //    // Apply the camera translation
//                      pos_proj[0] += cam_t[0];
//                      pos_proj[1] += cam_t[1];
//                      pos_proj[2] += cam_t[2];

//                  //    //prepare l
//                      double R[9];
//                  //    ceres::MatrixAdapter<T, 3, 3>(R);
//                      ceres::AngleAxisToRotationMatrix(cam_R, R);

//                      std::vector<double> it_intrinsics = sfm_data_.intrinsics[view->id_intrinsic]->getParams();

//                      const double f = it_intrinsics[0]; //focal
//                      const double u = it_intrinsics[1]; //principal_point_x
//                      const double v = it_intrinsics[2]; //principal_point_y

//                      // Transform the point from homogeneous to euclidean (undistorted point)
//                      const double x_i = u + f * (pos_proj[0] / pos_proj[2]);
//                      const double y_i = v + f * (pos_proj[1] / pos_proj[2]);

//                      const double P1 = f*R[0] + u*R[2];
//                      const double P2 = f*R[3] + u*R[5];
//                      const double P3 = f*R[6] + u*R[8];
//                      const double P4 = f*cam_t[0] + u*cam_t[2];
//                      const double P5 = f*R[1] + v*R[2];
//                      const double P6 = f*R[4] + v*R[5];
//                      const double P7 = f*R[7] + v*R[8];
//                      const double P8 = f*cam_t[1] + v*cam_t[2];
//                      const double P9 = R[2];
//                      const double P10 = R[5];
//                      const double P11 = R[8];
//                      const double P12 = cam_t[2];

//                      const double xu = obs_it.second.x(0);
//                      const double yu = obs_it.second.x(1);

//                      const double nc_x = P3/P11;
//                      const double nc_y = P7/P11;
//                  //    const T nc_x = P3;
//                  //    const T nc_y = P7;
//                      const double A = yu - nc_y;
//                      const double B = nc_x - xu;
//                      const double C = xu*nc_y - yu*nc_x;

//                      const double D = abs(A*x_i + B*y_i + C) / sqrt(A*A+B*B);
//                      //
//                      const double l_ncxi_x = x_i - nc_x;
//                      const double l_ncxi_y = y_i - nc_y;
//                      const double normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
//                  //    const T normal_ncxi = (l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
//                      //
//                      const double l_ncxu_x = xu - nc_x;
//                      const double l_ncxu_y = yu - nc_y;
//                      const double normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);

//                      double out_residuals[2];
//                      out_residuals[0] = (l_ncxu_x/normal_ncxu - l_ncxi_x/normal_ncxi);
//                      out_residuals[1] = (l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi);

//                      out << "R : " << cam_R[0] << " " << cam_R[1] << " " << cam_R[2] << std::endl;
//                      out << "t : " << cam_t[0] << " " << cam_t[1] << " " << cam_t[2] << std::endl;
//                      out << "f : " << f << " u : " << u << " v : " << v << std::endl;
//                      out << "x_i : " << x_i << " y_i : " << y_i << std::endl;
//                      out << "xu : " << xu << " yu : " << yu << std::endl;
//                      out << "nc_x : " << nc_x << " nc_y : " << nc_y << std::endl;
//                      out << "l_ncxi_x : " << l_ncxi_x << " l_ncxi_y : " << l_ncxi_y << std::endl;
//                      out << "l_ncxu_x : " << l_ncxu_x << " l_ncxu_y : " << l_ncxu_y << std::endl;
//                      out << "normal_ncxi : " << normal_ncxi << " normal_ncxu : " << normal_ncxu << std::endl;
//                      out << "out_residuals[0] : " << out_residuals[0] << " out_residuals[1] : " << out_residuals[1] << std::endl;
//                      out << "X : " << X[0] << " " << X[1] << " " << X[2] << std::endl;
//                  }


////                  if(out_residuals > 0.1)
////                  {
////                      out << "cam_R : " << cam_R[0] << " " << cam_R[1] << " " << cam_R[2] << std::endl;
//////                      out << "R : " << R[0] << " " << R[1] << " " << R[2]
//////                                          << R[3] << " " << R[4] << " " << R[5]
//////                                          << R[6] << " " << R[7] << " " << R[8] << std::endl;
////                      out << "cam_t : " << cam_t(0) << " " << cam_t(1) << " " << cam_t(2) << std::endl;
////                      out << "f : " << f << " u : " << u << " v : " << v << std::endl;
////                      out << "xu : " << xu << " yu : " << yu << std::endl;
////                      out << "pos_3dpoint[0]: " << pos_3dpoint[0] << " pos_3dpoint[1]: " << pos_3dpoint[1] << " pos_3dpoint[2]: " << pos_3dpoint[2] << std::endl;
////                      out << "out_residuals[0]: "  << out_residuals << std::endl;
////                      out << "out_residuals[0]_new: "  << out_residuals_new << std::endl;
////                      out << "****-----------------****" << std::endl;
////                  }



//            }


//          }
//      }

//      out.close();
  }

//  // - refine only Structure and Rotations & translations
//  b_BA_Status = bundle_adjustment_obj.Adjust_water3_c1
//    (
//      sfm_data_,
//      Optimize_Options(
//        ReconstructionEngine::intrinsic_refinement_options_,//Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_ALL,
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );
//  if (b_BA_Status && !sLogging_file_.empty())
//  {
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_RT_Xi", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
//  }

  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust_water6_c1
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  Optimize_Options ba_refine_options;
  if (ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
      ba_refine_options = Optimize_Options(
        ReconstructionEngine::intrinsic_refinement_options_,
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }else{
      ba_refine_options = Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }

  b_BA_Status = bundle_adjustment_obj.Adjust_water6_c1(sfm_data_, ba_refine_options);//(sfm_data_, ba_refine_options);//
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  return b_BA_Status;
}



bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_init_water4_c1()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );

//  // Remove outliers (max_angle, residual error)
//  const size_t pointcount_initial1 = sfm_data_.structure.size();
//  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//  const size_t pointcount_pixelresidual_filter1 = sfm_data_.structure.size();
//  RemoveOutliers_AngleError(sfm_data_, 2.0);
//  const size_t pointcount_angular_filter1 = sfm_data_.structure.size();
//  std::cout << "Outlier removal (remaining #points):\n"
//    << "\t initial structure size #3DPoints: " << pointcount_initial1 << "\n"
//    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter1 << "\n"
//    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter1 << std::endl;


//  bool b_BA_Status// = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::NONE, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );


//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  if (b_BA_Status)
  {
//    if (!sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

//    b_BA_Status = bundle_adjustment_obj.Adjust
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
//    {
//        // Remove outliers (max_angle, residual error)
//        const size_t pointcount_initial = sfm_data_.structure.size();
//        RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//        const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
//        RemoveOutliers_AngleError(sfm_data_, 2.0);
//        const size_t pointcount_angular_filter = sfm_data_.structure.size();
//        std::cout << "Outlier removal (remaining #points):\n"
//          << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
//          << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
//          << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;
//        Save(sfm_data_,
//          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_remove", "ply"),
//          ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust_water4_c1
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }
  std::cout << "water 3 C1 " << std::endl;

//  // - refine only Structure and Rotations & translations
//  b_BA_Status = bundle_adjustment_obj.Adjust_water3_c1
//    (
//      sfm_data_,
//      Optimize_Options(
//        ReconstructionEngine::intrinsic_refinement_options_,//Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_ALL,
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );
//  if (b_BA_Status && !sLogging_file_.empty())
//  {
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_RT_Xi", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
//  }

  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust_water4_c1
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  Optimize_Options ba_refine_options;
  if (ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
      ba_refine_options = Optimize_Options(
        ReconstructionEngine::intrinsic_refinement_options_,
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }else{
      ba_refine_options = Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }

  b_BA_Status = bundle_adjustment_obj.Adjust_water4_c1(sfm_data_, ba_refine_options);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  return b_BA_Status;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_init_water4_c1_2()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );

//  // Remove outliers (max_angle, residual error)
//  const size_t pointcount_initial1 = sfm_data_.structure.size();
//  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//  const size_t pointcount_pixelresidual_filter1 = sfm_data_.structure.size();
//  RemoveOutliers_AngleError(sfm_data_, 2.0);
//  const size_t pointcount_angular_filter1 = sfm_data_.structure.size();
//  std::cout << "Outlier removal (remaining #points):\n"
//    << "\t initial structure size #3DPoints: " << pointcount_initial1 << "\n"
//    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter1 << "\n"
//    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter1 << std::endl;


//  bool b_BA_Status// = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::NONE, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );


//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  if (b_BA_Status)
  {
//    if (!sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

//    b_BA_Status = bundle_adjustment_obj.Adjust
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
//    {
//        // Remove outliers (max_angle, residual error)
//        const size_t pointcount_initial = sfm_data_.structure.size();
//        RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//        const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
//        RemoveOutliers_AngleError(sfm_data_, 2.0);
//        const size_t pointcount_angular_filter = sfm_data_.structure.size();
//        std::cout << "Outlier removal (remaining #points):\n"
//          << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
//          << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
//          << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;
//        Save(sfm_data_,
//          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_remove", "ply"),
//          ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust_water4_c1
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, //ReconstructionEngine::intrinsic_refinement_options_,//Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );

    {
        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_0_refine_RT_Xi", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));

        SfM_Data_Structure_Computation_Blind structure_estimator(true);
        structure_estimator.triangulate(sfm_data_);

        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_1_refine_RT_Xi", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));
    }
    //if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }
  std::cout << "water 3 C1 Structure_Parameter_Type::NONE" << std::endl;

//  // - refine only Structure and Rotations & translations
//  b_BA_Status = bundle_adjustment_obj.Adjust_water3_c1
//    (
//      sfm_data_,
//      Optimize_Options(
//        ReconstructionEngine::intrinsic_refinement_options_,//Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_ALL,
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );
//  if (b_BA_Status && !sLogging_file_.empty())
//  {
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_RT_Xi", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
//  }

  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust_water4_c1
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  // Check that poses & intrinsic cover some measures (after outlier removal)
  const IndexT minPointPerPose = 12; // 6 min
  const IndexT minTrackLength = 3; // 2 min
  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
  {
    // TODO: must ensure that track graph is producing a single connected component

    const size_t pointcount_cleaning = sfm_data_.structure.size();
    std::cout << "Point_cloud cleaning:\n"
      << "\t #3DPoints: " << pointcount_cleaning << "\n";
  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  Optimize_Options ba_refine_options;
  if (ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
      ba_refine_options = Optimize_Options(
        ReconstructionEngine::intrinsic_refinement_options_,
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }else{
      ba_refine_options = Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }

  b_BA_Status = bundle_adjustment_obj.Adjust_water4_c1(sfm_data_, ba_refine_options);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  return b_BA_Status;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_init_water4_c4()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );

//  // Remove outliers (max_angle, residual error)
//  const size_t pointcount_initial1 = sfm_data_.structure.size();
//  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//  const size_t pointcount_pixelresidual_filter1 = sfm_data_.structure.size();
//  RemoveOutliers_AngleError(sfm_data_, 2.0);
//  const size_t pointcount_angular_filter1 = sfm_data_.structure.size();
//  std::cout << "Outlier removal (remaining #points):\n"
//    << "\t initial structure size #3DPoints: " << pointcount_initial1 << "\n"
//    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter1 << "\n"
//    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter1 << std::endl;


//  bool b_BA_Status// = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::NONE, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );


//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  if (b_BA_Status)
  {
//    if (!sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

      std::cout << "water 4 c4 1" << std::endl;
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,//Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
    {
        std::cout << "water 4 c4 2" << std::endl;
        // Remove outliers (max_angle, residual error)
        const size_t pointcount_initial = sfm_data_.structure.size();
        RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
        const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
        RemoveOutliers_AngleError(sfm_data_, 2.0);
        const size_t pointcount_angular_filter = sfm_data_.structure.size();
        std::cout << "Outlier removal (remaining #points):\n"
          << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
          << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
          << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;
        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_remove", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    std::cout << "water 4 c4 3" << std::endl;
    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust_water4_c4
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,//Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }
  std::cout << "water 3 C1 " << std::endl;

//  // - refine only Structure and Rotations & translations
//  b_BA_Status = bundle_adjustment_obj.Adjust_water3_c1
//    (
//      sfm_data_,
//      Optimize_Options(
//        ReconstructionEngine::intrinsic_refinement_options_,//Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_ALL,
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );
//  if (b_BA_Status && !sLogging_file_.empty())
//  {
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_RT_Xi", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
//  }

  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
      std::cout << "water 4 c4 4" << std::endl;
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust_water4_c4
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  std::cout << "water 4 c4 5" << std::endl;
  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  std::cout << "water 4 c4 6" << std::endl;
  Optimize_Options ba_refine_options;
  if (ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
      ba_refine_options = Optimize_Options(
        ReconstructionEngine::intrinsic_refinement_options_,
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }else{
      ba_refine_options = Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }

  std::cout << "water 4 c4 8" << std::endl;
  b_BA_Status = bundle_adjustment_obj.Adjust_water4_c4(sfm_data_, ba_refine_options);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }
  std::cout << "water 4 c4 8" << std::endl;

  return b_BA_Status;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_init_6()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = true;

  std::vector<int> fileIdList;
  fileIdList.push_back(1);
  fileIdList.push_back(4);
  fileIdList.push_back(2);
  fileIdList.push_back(5);
  fileIdList.push_back(3);
  fileIdList.push_back(6);
  std::map<int, bool> controlImgMap;
  std::map<int, bool> controlStructureMap;
  SfM_Data lastData = sfm_data_;
  for(std::size_t t = 0; t < 6; ++t)
//  for(std::size_t t = 0; t < 3; ++t)
  {
      std::map<int, bool> computeImgMap;
      std::map<int, bool> computeImgStructureMap_temp;
//      int t = 0;
      /// step 0
      ///
      int fId = fileIdList[t];
      std::ifstream inF;
      std::stringstream ss;
      ss << fId;
      std::string s;
      ss >> s;
      std::string fPath = sfm_data_.s_root_path +"/image_seq" + s + ".res";
      std::cout << "open file " << fPath << std::endl;
      inF.open(fPath);
      if(inF)
      {
          int count = 0;
          int line = 0;
          inF >> count;
          while(line<count && !inF.eof())
          {
              int thisId;
              inF >> thisId;
              std::pair<int, bool> tempPair;
              tempPair.first = thisId;
              tempPair.second = true;
              computeImgMap.insert(tempPair);
              ++line;
          }

      }else{
          std::cout << "Open file " << fPath << " failed" << std::endl;
          return false;
      }
      inF.close();
      {
//          std::string s2;
//          if(t == 0)
//          {
//              s2 = "4";
//          }else if(t == 1)
//          {
//              s2 = "5";
//          }else if(t == 2)
//          {
//              s2 = "6";
//          }
//          std::string fPath = sfm_data_.s_root_path +"/image_seq"+s2+".res";
//          std::cout << "open file " << fPath << std::endl;
//          inF.open(fPath);
//          if(inF)
//          {
//              int count = 0;
//              int line = 0;
//              inF >> count;
//              while(line<count && !inF.eof())
//              {
//                  int thisId;
//                  inF >> thisId;
//                  std::pair<int, bool> tempPair;
//                  tempPair.first = thisId;
//                  tempPair.second = true;
//                  computeImgMap.insert(tempPair);
//                  ++line;
//              }

//          }else{
//              std::cout << "Open file " << fPath << " failed" << std::endl;
//              return false;
//          }
//          inF.close();
      }
      //find structure id that need to compute

      /// step 1
      ///
      std::cout << "r " << s << " step 1 " << std::endl;
      Landmarks newLandmarks_step1;
      b_BA_Status = bundle_adjustment_obj.Adjust_6
        (
          controlImgMap,
          controlStructureMap,
          computeImgMap,
          computeImgStructureMap_temp,
          newLandmarks_step1,
          lastData,
          Optimize_Options(
            Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant//ReconstructionEngine::intrinsic_refinement_options_,//
            Extrinsic_Parameter_Type::ADJUST_ALL,
            Structure_Parameter_Type::ADJUST_ALL,
            Control_Point_Parameter(),
            this->b_use_motion_prior_)
        );
      SfM_Data newData_step1;
      newData_step1 = lastData;
//      newData_step1.intrinsics = sfm_data_.intrinsics;
//      for(std::size_t idView = 0; idView < sfm_data_.views.size(); ++idView)
//      {
//          Views::const_iterator itV = sfm_data_.views.begin();
//          std::advance(itV, idView);
//          if(computeImgMap.find(itV->second.get()->id_pose) == computeImgMap.end()
//          && controlImgMap.find(itV->second.get()->id_pose) == controlImgMap.end())
//          {
//              continue;
//          }
//          newData_step1.poses[itV->second.get()->id_pose] = sfm_data_.poses[itV->second.get()->id_pose];
//          newData_step1.views[itV->first] = sfm_data_.views[itV->first];
//      }
      newData_step1.structure = newLandmarks_step1;
      if (!sLogging_file_.empty())
      {
        Save(lastData,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "iterator"+s+"_step1_"+"structure_00_refine_T_Xi", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));
        Save(newData_step1,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "iterator"+s+"_step1_"+"structure_00_refine_T_Xi_new_data", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));
      }

      /// step 2
      ///
      std::cout << "r " << s << " step 2 " << std::endl;
      Landmarks newLandmarks_step2;
      SfM_Data newData_step2;
      if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
          // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
          //if(t == 0 || t == 2 || t == 4)
          {
              b_BA_Status = bundle_adjustment_obj.Adjust_6
                (
                  controlImgMap,
                  controlStructureMap,
                  computeImgMap,
                  computeImgStructureMap_temp,
                  newLandmarks_step2,
                  newData_step1,
                  Optimize_Options(
                    ReconstructionEngine::intrinsic_refinement_options_,
                    Extrinsic_Parameter_Type::ADJUST_ALL,
                    Structure_Parameter_Type::ADJUST_ALL,
                    Control_Point_Parameter(),
                    this->b_use_motion_prior_)
                );
          }
          //else{
//              b_BA_Status = bundle_adjustment_obj.Adjust_6
//                (
//                  controlImgMap,
//                  controlStructureMap,
//                  computeImgMap,
//                  computeImgStructureMap_temp,
//                  newLandmarks_step2,
//                  newData_step1,
//                  Optimize_Options(
//                    Intrinsic_Parameter_Type::NONE,//ReconstructionEngine::intrinsic_refinement_options_,//
//                    Extrinsic_Parameter_Type::ADJUST_ALL,
//                    Structure_Parameter_Type::ADJUST_ALL,
//                    Control_Point_Parameter(),
//                    this->b_use_motion_prior_)
//                );
          //}

          newData_step2 = newData_step1;
//          newData_step2.intrinsics = newData_step1.intrinsics;
//          for(std::size_t idView = 0; idView < newData_step1.views.size(); ++idView)
//          {
//              Views::const_iterator itV = newData_step1.views.begin();
//              std::advance(itV, idView);
//              if(computeImgMap.find(itV->second.get()->id_pose) == computeImgMap.end()
//              && controlImgMap.find(itV->second.get()->id_pose) == controlImgMap.end())
//              {
//                  continue;
//              }
//              newData_step2.poses[itV->second.get()->id_pose] = newData_step1.poses[itV->second.get()->id_pose];
//              newData_step2.views[itV->first] = newData_step1.views[itV->first];
//          }
          newData_step2.structure = newLandmarks_step2;
          if (b_BA_Status && !sLogging_file_.empty())
          {
            Save(newData_step1,
              stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "iterator"+s+"_step2_"+"structure_02_refine_KRT_Xi", "ply"),
              ESfM_Data(EXTRINSICS | STRUCTURE));
            Save(newData_step2,
              stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "iterator"+s+"_step2_"+"structure_02_refine_KRT_Xi_new_data", "ply"),
              ESfM_Data(EXTRINSICS | STRUCTURE));
          }
      }else{
          newData_step2 = newData_step1;
      }

      /// step 3
      ///
      // Remove outliers (max_angle, residual error)
      std::cout << "r " << s << " step 3 " << std::endl;
      const size_t pointcount_initial = newData_step2.structure.size();
      RemoveOutliers_PixelResidualError(newData_step2, 4.0);
      const size_t pointcount_pixelresidual_filter = newData_step2.structure.size();
      RemoveOutliers_AngleError(newData_step2, 2.0);
      const size_t pointcount_angular_filter = newData_step2.structure.size();
      std::cout << "Outlier removal (remaining #points):\n"
        << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
        << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
        << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

      if (!sLogging_file_.empty())
      {
        Save(newData_step2,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "iterator"+s+"_step3_"+"structure_03_outlier_removed_new_data", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));
      }

      /// step 4
      ///
      // --
      // Final BA. We refine one more time,
      // since some outlier have been removed and so a better solution can be found.
      //--
      std::cout << "r " << s << " step 4 " << std::endl;
      Optimize_Options ba_refine_options;
      if (ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE){// && (t == 0 || t == 2 || t == 4)) {
          ba_refine_options = Optimize_Options(
            ReconstructionEngine::intrinsic_refinement_options_,
            Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
            Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
            Control_Point_Parameter(),
            this->b_use_motion_prior_);
      }else{
          ba_refine_options = Optimize_Options(
            Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
            Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
            Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
            Control_Point_Parameter(),
            this->b_use_motion_prior_);
      }
      Landmarks newLandmarks_step4;
      std::map<int, bool> computeImgStructureMap_res;
      b_BA_Status = bundle_adjustment_obj.Adjust_6(controlImgMap, controlStructureMap,computeImgMap, computeImgStructureMap_res, newLandmarks_step4, newData_step2, ba_refine_options);
      SfM_Data newData_step4;
      newData_step4 = newData_step2;
      newData_step4.structure = newLandmarks_step4;
      if (b_BA_Status && !sLogging_file_.empty())
      {
        Save(newData_step4,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "iterator"+s+"_step4_"+"structure_04_outlier_removed_new_data", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));
      }

      if(t == 5)
      {
          sfm_data_ = newData_step4;
          break;
      }
      /// step 5 prepare for next interate
      ///
      std::cout << "r " << s << " step 5 " << std::endl;
      lastData.intrinsics = newData_step4.intrinsics;
      #pragma omp parallel for
      for(std::size_t idV = 0; idV < newData_step4.views.size(); ++idV)
      {
          Views::const_iterator itV = newData_step4.views.begin();
          std::advance(itV, idV);
          //
          lastData.views[itV->first] = itV->second;//newData_step4.views[itV->first];
          lastData.poses[itV->second.get()->id_pose] = newData_step4.poses[itV->second.get()->id_pose];
      }
      #pragma omp parallel for
      for(std::size_t idS = 0; idS < newData_step4.structure.size(); ++idS)
      {
          Landmarks::const_iterator itS = newData_step4.structure.begin();
          std::advance(itS, idS);
          //
          lastData.structure[itS->first] = itS->second;
      }


      #pragma omp parallel for
      for(std::size_t idCI = 0; idCI < computeImgMap.size(); ++idCI)
      {
          std::map<int, bool>::const_iterator itCI = computeImgMap.begin();
          std::advance(itCI, idCI);
          if(controlImgMap.find(itCI->first) == controlImgMap.end())
          {
              #pragma omp critical
              {
                  controlImgMap.insert(*itCI);
              }
          }
      }
      std::cout << "computeImgMap : " << computeImgMap.size() << std::endl;
      std::cout << "controlImgMap : " << controlImgMap.size() << std::endl;
      #pragma omp parallel for
      for(std::size_t idCS = 0; idCS < computeImgStructureMap_res.size(); ++idCS)
      {
          std::map<int, bool>::const_iterator itCS = computeImgStructureMap_res.begin();
          std::advance(itCS, idCS);
          if(controlStructureMap.find(itCS->first) == controlStructureMap.end())
          {
              #pragma omp critical
              {
                  controlStructureMap.insert(*itCS);
              }
          }
      }
      std::cout << "computeImgStructureMap_res : " << computeImgStructureMap_res.size() << std::endl;
      std::cout << "controlStructureMap : " << controlStructureMap.size() << std::endl;
      std::cout << "r " << s << " finished " << std::endl;

      if(t == 0)
      {
          SfM_Data newSfMData;
          newSfMData = newData_step4;
//          std::stringstream bIdss;
//          bIdss << bId;
//          std::string bIds;
//          bIdss >> bIds;

          std::map<int, Vec3> c_gps;
          std::ifstream in_gps;
          std::string inGpsPath = newData_step4.s_root_path + "/n_gps_imu_321.res";
          std::cout << "inGpsPath : " << inGpsPath << std::endl;
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

          double length_gps = 0.0;
          double length_c = 0.0;
//          sfm::Poses::iterator itEnd = sfm_data.poses.end();
//          -- itEnd;
          int length_N = 0;
          for(std::size_t poseListId1 = 0; poseListId1 != computeImgMap.size(); ++poseListId1)
    //      for(sfm::Poses::iterator itPose = sfm_data.poses.begin();
    //          itPose != itEnd; ++itPose)
          {
              int pId1 = computeImgMap[poseListId1];
              const Vec3 c1 = newData_step4.poses[pId1].center();

              if (pId1/100000 != 2)
              {
                  continue;
              }
              const Vec3 cg1 = c_gps.find(pId1)->second;

    //          sfm::Poses::iterator itPose2 = itPose;
    //          ++ itPose2;
    //          for(;itPose2 != sfm_data.poses.end(); ++itPose2)
              for(std::size_t poseListId2 = poseListId1; poseListId2 != computeImgMap.size(); ++poseListId2)
              {
                  int pId2 = computeImgMap[poseListId2];
                  const Vec3 c2 = newData_step4.poses[pId2].center();

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

          //compute H -> Xg = H * Xc
          std::vector<Vec3> Xc, Xg;
          for(std::size_t poseListId = 0; poseListId != computeImgMap.size(); ++poseListId)
          {
              int pId = computeImgMap[poseListId];
              if (pId/100000 != 2)
              {
                  continue;
              }
              //prepare Xg
              Vec3 Xg_temp;
              Xg_temp = c_gps.find(pId)->second;
              Xg.push_back(Xg_temp);

              //prepare Xc
              Vec3 Xc_temp;
              Xc_temp = s * newData_step4.poses[pId].center();//itPose->second.center();
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
          //
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
          //
          Mat3 H_gc;
          H_gc << 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0;
          for(std::size_t i = 0; i < Xc.size(); ++i)
          {
              H_gc += Xg_bary[i] * Xc_bary[i].transpose();
          }
          //
          Eigen::JacobiSVD<Mat3> svd_gc(H_gc, Eigen::ComputeFullU | Eigen::ComputeFullV);
          Mat3 s_gc = svd_gc.matrixU();
          Mat3 d_gc = svd_gc.matrixV();
          //
          Mat3 R_gc;
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
          //
          Vec3 t_gc;
          t_gc << 0.0, 0.0, 0.0;
          for(std::size_t i = 0; i < Xc.size(); ++i)
          {
              Vec3 t_gc_temp;
              t_gc_temp = Xg[i] - R_gc * Xc[i];
              t_gc += t_gc_temp;
          }
          t_gc /= Xc.size();
          //
          std::cout << "s: " << s << std::endl;
          std::cout << "R: " << R_gc << std::endl;
          std::cout << "t: " << t_gc << std::endl;
    //      for (Poses::const_iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose)
          for(std::size_t poseListId = 0; poseListId != computeImgMap.size(); ++poseListId)
          {
              const Pose3 & pose = newData_step4.poses[computeImgMap[poseListId]];
              const Mat3 R = pose.rotation();
              const Vec3 t = pose.translation();
              const Vec3 C = pose.center();
              Mat3 newR = R * R_gc.inverse();
              Vec3 newC = R_gc * s * C + t_gc;
              newSfMData.poses[computeImgMap[poseListId]] = Pose3(newR, newC);
          }
          //set all X
          for(Landmarks::const_iterator itX = newData_step4.structure.begin();
              itX != newData_step4.structure.end(); ++itX)
          {
              const Vec3 & X = itX->second.X;
              Vec3 newX = R_gc * s* X + t_gc;
              newSfMData.structure[itX->first].X = newX;
          }
          Save(newSfMData,
            stlplus::create_filespec(newSfMData.s_root_path, "new_cloud_and_poses", ".ply"),
            ESfM_Data(ALL));
          Save(newSfMData,
               stlplus::create_filespec(newSfMData.s_root_path, "new_cloud_and_poses", ".bin"),
               ESfM_Data(ALL));

          std::cout << "save newSfMData finish" << std::endl;

          newData_step4 = newSfMData;
      }



  }







//  if (b_BA_Status)
//  {
//    if (!sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }

//    // - refine only Structure and Rotations & translations
//    b_BA_Status = bundle_adjustment_obj.Adjust
//      (
//        sfm_data_,
//        Optimize_Options(
//          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    if (b_BA_Status && !sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }
//  }


//  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
//    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
//    b_BA_Status = bundle_adjustment_obj.Adjust
//      (
//        sfm_data_,
//        Optimize_Options(
//          ReconstructionEngine::intrinsic_refinement_options_,
//          Extrinsic_Parameter_Type::ADJUST_ALL,
//          Structure_Parameter_Type::ADJUST_ALL,
//          Control_Point_Parameter(),
//          this->b_use_motion_prior_)
//      );
//    if (b_BA_Status && !sLogging_file_.empty())
//    {
//      Save(sfm_data_,
//        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
//    }
//  }

//  // Remove outliers (max_angle, residual error)
//  const size_t pointcount_initial = sfm_data_.structure.size();
//  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
//  RemoveOutliers_AngleError(sfm_data_, 2.0);
//  const size_t pointcount_angular_filter = sfm_data_.structure.size();
//  std::cout << "Outlier removal (remaining #points):\n"
//    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
//    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
//    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

//  if (!sLogging_file_.empty())
//  {
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
//  }

//  // --
//  // Final BA. We refine one more time,
//  // since some outlier have been removed and so a better solution can be found.
//  //--
//  Optimize_Options ba_refine_options;
//  if (ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
//      ba_refine_options = Optimize_Options(
//        ReconstructionEngine::intrinsic_refinement_options_,
//        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
//        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_);
//  }else{
//      ba_refine_options = Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
//        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_);
//  }

//  b_BA_Status = bundle_adjustment_obj.Adjust(sfm_data_, ba_refine_options);
//  if (b_BA_Status && !sLogging_file_.empty())
//  {
//    Save(sfm_data_,
//      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
//      ESfM_Data(EXTRINSICS | STRUCTURE));
//  }

  return b_BA_Status;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_change_step(std::map<int, Pose3> &fixIdList)
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );

//  // Remove outliers (max_angle, residual error)
//  const size_t pointcount_initial1 = sfm_data_.structure.size();
//  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//  const size_t pointcount_pixelresidual_filter1 = sfm_data_.structure.size();
//  RemoveOutliers_AngleError(sfm_data_, 2.0);
//  const size_t pointcount_angular_filter1 = sfm_data_.structure.size();
//  std::cout << "Outlier removal (remaining #points):\n"
//    << "\t initial structure size #3DPoints: " << pointcount_initial1 << "\n"
//    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter1 << "\n"
//    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter1 << std::endl;


//  bool b_BA_Status// = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::NONE, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );


//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  if (b_BA_Status)
  {
    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust_fix_some
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_),
        fixIdList
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }


  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust_fix_some
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_),
        fixIdList
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  Optimize_Options ba_refine_options;
  if (ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
      ba_refine_options = Optimize_Options(
        ReconstructionEngine::intrinsic_refinement_options_,
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }else{
      ba_refine_options = Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }

  b_BA_Status = bundle_adjustment_obj.Adjust_fix_some(sfm_data_, ba_refine_options, fixIdList);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  return b_BA_Status;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_for_large()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );

//  // Remove outliers (max_angle, residual error)
//  const size_t pointcount_initial1 = sfm_data_.structure.size();
//  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//  const size_t pointcount_pixelresidual_filter1 = sfm_data_.structure.size();
//  RemoveOutliers_AngleError(sfm_data_, 2.0);
//  const size_t pointcount_angular_filter1 = sfm_data_.structure.size();
//  std::cout << "Outlier removal (remaining #points):\n"
//    << "\t initial structure size #3DPoints: " << pointcount_initial1 << "\n"
//    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter1 << "\n"
//    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter1 << std::endl;


//  bool b_BA_Status// = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::NONE, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );


//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  if (b_BA_Status)
  {
    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust_down
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );

    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi_down", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    b_BA_Status = bundle_adjustment_obj.Adjust_front
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );

    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi_front", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    b_BA_Status = bundle_adjustment_obj.Adjust_back
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );

    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi_back", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );

    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi_down_front_back", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

  }


  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  Optimize_Options ba_refine_options;
  if (ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
      ba_refine_options = Optimize_Options(
        ReconstructionEngine::intrinsic_refinement_options_,
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }else{
      ba_refine_options = Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_);
  }


  b_BA_Status = bundle_adjustment_obj.Adjust(sfm_data_, ba_refine_options);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  return b_BA_Status;
}



bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_init_fix2()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );

//  // Remove outliers (max_angle, residual error)
//  const size_t pointcount_initial1 = sfm_data_.structure.size();
//  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
//  const size_t pointcount_pixelresidual_filter1 = sfm_data_.structure.size();
//  RemoveOutliers_AngleError(sfm_data_, 2.0);
//  const size_t pointcount_angular_filter1 = sfm_data_.structure.size();
//  std::cout << "Outlier removal (remaining #points):\n"
//    << "\t initial structure size #3DPoints: " << pointcount_initial1 << "\n"
//    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter1 << "\n"
//    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter1 << std::endl;


//  bool b_BA_Status// = true;
//          = bundle_adjustment_obj.Adjust
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::NONE, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );


//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  if (b_BA_Status)
  {
    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust_fix2
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }


  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust_fix2
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  const Optimize_Options ba_refine_options(
    ReconstructionEngine::intrinsic_refinement_options_,
    Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
    Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
    Control_Point_Parameter(),
    this->b_use_motion_prior_);

  b_BA_Status = bundle_adjustment_obj.Adjust_fix2(sfm_data_, ba_refine_options);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  return b_BA_Status;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_init_old()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = true;
//          = bundle_adjustment_obj.Adjust_old
//    (
//      sfm_data_,
//      Optimize_Options(
//        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::NONE, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        this->b_use_motion_prior_)
//    );

  if (b_BA_Status)
  {
    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust_old
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust_old
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

//  // Check that poses & intrinsic cover some measures (after outlier removal)
//  const IndexT minPointPerPose = 12; // 6 min
//  const IndexT minTrackLength = 3; // 2 min
//  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
//  {
//    // TODO: must ensure that track graph is producing a single connected component

//    const size_t pointcount_cleaning = sfm_data_.structure.size();
//    std::cout << "Point_cloud cleaning:\n"
//      << "\t #3DPoints: " << pointcount_cleaning << "\n";
//  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  const Optimize_Options ba_refine_options(
    ReconstructionEngine::intrinsic_refinement_options_,
    Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
    Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
    Control_Point_Parameter(),
    this->b_use_motion_prior_);

  b_BA_Status = bundle_adjustment_obj.Adjust_old(sfm_data_, ba_refine_options);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  return b_BA_Status;
}


// Adjust the scene for single cameras(down)
bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_single_fix_fuv()
{

    Bundle_Adjustment_Ceres bundle_adjustment_obj;
    bool b_BA_Status;
    {
                std::cout<< "initial intrinsics :"<<std::endl;
                    Hash_Map<IndexT, std::vector<double> > map_intrinsics;
                    map_intrinsics[0] = sfm_data_.intrinsics[0]->getParams();
                    std::cout<< "map_intrinsics " << map_intrinsics[0][0]<<" "<<map_intrinsics[0][1]<<" "<<
                             map_intrinsics[0][2]<<" "<<map_intrinsics[0][3]<<" "<<map_intrinsics[0][4]<<" "<<
                             map_intrinsics[0][5]<<" "<<map_intrinsics[0][6]<<" "<<map_intrinsics[0][7]<<std::endl;
    }
    b_BA_Status = bundle_adjustment_obj.Adjust_single_fix_fuv
      (
        sfm_data_,
        Optimize_Options(
        Intrinsic_Parameter_Type::NONE, //
        Extrinsic_Parameter_Type::ADJUST_TRANSLATION, //
        Structure_Parameter_Type::ADJUST_ALL)
      );

    if (b_BA_Status)
    {
      if (!sLogging_file_.empty())
      {
        Save(sfm_data_,
         stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_01", "ply"),
         ESfM_Data(EXTRINSICS | STRUCTURE));
      }

      b_BA_Status = bundle_adjustment_obj.Adjust_single_fix_fuv
        (
          sfm_data_,
          Optimize_Options(
          Intrinsic_Parameter_Type::NONE, //
          Extrinsic_Parameter_Type::ADJUST_ALL, //
          Structure_Parameter_Type::ADJUST_ALL)
        );

      if (!sLogging_file_.empty() && b_BA_Status)
      {
        Save(sfm_data_,
         stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_02", "ply"),
         ESfM_Data(EXTRINSICS | STRUCTURE));
      }

      b_BA_Status = bundle_adjustment_obj.Adjust_single_fix_fuv
        (
          sfm_data_,
          Optimize_Options(
          Intrinsic_Parameter_Type::ADJUST_DISTORTION,//ReconstructionEngine::intrinsic_refinement_options_,//Intrinsic_Parameter_Type::NONE, //
          Extrinsic_Parameter_Type::ADJUST_ALL, //
          Structure_Parameter_Type::ADJUST_ALL)
        );


      {
                  std::cout<< "result intrinsics 1 :"<<std::endl;
                      Hash_Map<IndexT, std::vector<double> > map_intrinsics;
                      map_intrinsics[0] = sfm_data_.intrinsics[0]->getParams();
                      std::cout<< "map_intrinsics " << map_intrinsics[0][0]<<" "<<map_intrinsics[0][1]<<" "<<
                               map_intrinsics[0][2]<<" "<<map_intrinsics[0][3]<<" "<<map_intrinsics[0][4]<<" "<<
                               map_intrinsics[0][5]<<" "<<map_intrinsics[0][6]<<" "<<map_intrinsics[0][7]<<std::endl;
      }

      if (!sLogging_file_.empty() && b_BA_Status)
      {
        Save(sfm_data_,
         stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_03", "ply"),
         ESfM_Data(EXTRINSICS | STRUCTURE));
      }

      // Remove outliers (max_angle, residual error)
      const size_t pointcount_initial = sfm_data_.structure.size();
      RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
      const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
      RemoveOutliers_AngleError(sfm_data_, 2.0);
      const size_t pointcount_angular_filter = sfm_data_.structure.size();
      std::cout << "Outlier removal (remaining #points):\n"
        << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
        << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
        << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

      if (!sLogging_file_.empty())
      {
        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_04", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));
      }


      //return b_BA_Status;

      b_BA_Status = bundle_adjustment_obj.Adjust_single_fix_fuv
        (
          sfm_data_,
          Optimize_Options(
            Intrinsic_Parameter_Type::ADJUST_DISTORTION,//ReconstructionEngine::intrinsic_refinement_options_,//Intrinsic_Parameter_Type::NONE, //ReconstructionEngine::intrinsic_refinement_options_,  //Intrinsic_Parameter_Type::NONE, //
            Extrinsic_Parameter_Type::ADJUST_ALL,
            Structure_Parameter_Type::ADJUST_ALL)
        );
          if (b_BA_Status && !sLogging_file_.empty())
          {
            Save(sfm_data_,
              stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_05", "ply"),
              ESfM_Data(EXTRINSICS | STRUCTURE));
          }

          {
                      std::cout<< "result intrinsics 2 :"<<std::endl;
                          Hash_Map<IndexT, std::vector<double> > map_intrinsics;
                          map_intrinsics[0] = sfm_data_.intrinsics[0]->getParams();
                          std::cout<< "map_intrinsics " << map_intrinsics[0][0]<<" "<<map_intrinsics[0][1]<<" "<<
                                   map_intrinsics[0][2]<<" "<<map_intrinsics[0][3]<<" "<<map_intrinsics[0][4]<<" "<<
                                   map_intrinsics[0][5]<<" "<<map_intrinsics[0][6]<<" "<<map_intrinsics[0][7]<<std::endl;
          }
    }

    return b_BA_Status;
}

// Adjust the scene for single cameras(down), change function, c1
bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_down_gps_change()
{

    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_00", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));

    std::vector<double> out_transformation_R_imu(3);
    std::vector<double> out_tramsformation_x_gps(3);

    Bundle_Adjustment_Ceres bundle_adjustment_obj;
    bool b_BA_Status;
    b_BA_Status = bundle_adjustment_obj.Adjust_down_gps_cail_first_change
      (
        sfm_data_,
        Optimize_Options(
        Intrinsic_Parameter_Type::NONE, //econstructionEngine::intrinsic_refinement_options_,//
        Extrinsic_Parameter_Type::ADJUST_ALL, //
        Structure_Parameter_Type::NONE),//ADJUST_ALL
        out_transformation_R_imu,
        out_tramsformation_x_gps
      );
    if (b_BA_Status)
    {
      if (!sLogging_file_.empty())
      {
        Save(sfm_data_,
         stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_01", "ply"),
         ESfM_Data(EXTRINSICS | STRUCTURE));
      }

      b_BA_Status = bundle_adjustment_obj.Adjust_down_gps_cail_second_change
        (
          sfm_data_,
          Optimize_Options(
          Intrinsic_Parameter_Type::NONE, //
          Extrinsic_Parameter_Type::ADJUST_ALL, //
          Structure_Parameter_Type::ADJUST_ALL),
          out_transformation_R_imu,
          out_tramsformation_x_gps
        );

      if (!sLogging_file_.empty() && b_BA_Status)
      {
        Save(sfm_data_,
         stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_02", "ply"),
         ESfM_Data(EXTRINSICS | STRUCTURE));
      }

      // Remove outliers (max_angle, residual error)
      const size_t pointcount_initial = sfm_data_.structure.size();
      RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
      const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
      RemoveOutliers_AngleError(sfm_data_, 2.0);
      const size_t pointcount_angular_filter = sfm_data_.structure.size();
      std::cout << "Outlier removal (remaining #points):\n"
        << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
        << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
        << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

      if (!sLogging_file_.empty())
      {
        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_03", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));
      }


      //return b_BA_Status;

      b_BA_Status = bundle_adjustment_obj.Adjust_down_gps_cail_second_change
        (
          sfm_data_,
          Optimize_Options(
            Intrinsic_Parameter_Type::NONE, //ReconstructionEngine::intrinsic_refinement_options_,  //Intrinsic_Parameter_Type::NONE, //
            Extrinsic_Parameter_Type::ADJUST_ALL,
            Structure_Parameter_Type::ADJUST_ALL),
            out_transformation_R_imu,
            out_tramsformation_x_gps
        );
          if (b_BA_Status && !sLogging_file_.empty())
          {
            Save(sfm_data_,
              stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_04", "ply"),
              ESfM_Data(EXTRINSICS | STRUCTURE));
          }
    }

    return b_BA_Status;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_threecamera_gps_change()
{

    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_00", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));

    std::vector<double> out_transformation_R_imu(3);
    std::vector<double> out_tramsformation_x_gps(3);
    std::vector<double> out_transformation_br(6);
    std::vector<double> out_transformation_fr(6);

    Bundle_Adjustment_Ceres bundle_adjustment_obj;
    bool b_BA_Status;
    b_BA_Status = bundle_adjustment_obj.Adjust_threecamera_gps_cail_first_change_c4
      (
        sfm_data_,
        Optimize_Options(
        Intrinsic_Parameter_Type::NONE, //econstructionEngine::intrinsic_refinement_options_,//
        Extrinsic_Parameter_Type::ADJUST_ALL, //
        Structure_Parameter_Type::NONE),//ADJUST_ALL
        out_transformation_R_imu,
        out_tramsformation_x_gps,
        out_transformation_br,
        out_transformation_fr
      );
    if (b_BA_Status)
    {
      if (!sLogging_file_.empty())
      {
        Save(sfm_data_,
         stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_01", "ply"),
         ESfM_Data(EXTRINSICS | STRUCTURE));
      }

      b_BA_Status = bundle_adjustment_obj.Adjust_threecamera_gps_cail_second_change_c4
        (
          sfm_data_,
          Optimize_Options(
          Intrinsic_Parameter_Type::NONE, //
          Extrinsic_Parameter_Type::ADJUST_ALL, //
          Structure_Parameter_Type::ADJUST_ALL),
          out_transformation_R_imu,
          out_tramsformation_x_gps,
          out_transformation_br,
          out_transformation_fr
        );

      if (!sLogging_file_.empty() && b_BA_Status)
      {
        Save(sfm_data_,
         stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_02", "ply"),
         ESfM_Data(EXTRINSICS | STRUCTURE));
      }

      // Remove outliers (max_angle, residual error)
      const size_t pointcount_initial = sfm_data_.structure.size();
      RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
      const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
      RemoveOutliers_AngleError(sfm_data_, 2.0);
      const size_t pointcount_angular_filter = sfm_data_.structure.size();
      std::cout << "Outlier removal (remaining #points):\n"
        << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
        << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
        << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

      if (!sLogging_file_.empty())
      {
        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_03", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));
      }


      //return b_BA_Status;

      b_BA_Status = bundle_adjustment_obj.Adjust_threecamera_gps_cail_second_change_c4
        (
          sfm_data_,
          Optimize_Options(
            Intrinsic_Parameter_Type::NONE, //ReconstructionEngine::intrinsic_refinement_options_,  //Intrinsic_Parameter_Type::NONE, //
            Extrinsic_Parameter_Type::ADJUST_ALL,
            Structure_Parameter_Type::ADJUST_ALL),
            out_transformation_R_imu,
            out_tramsformation_x_gps,
            out_transformation_br,
            out_transformation_fr
        );
          if (b_BA_Status && !sLogging_file_.empty())
          {
            Save(sfm_data_,
              stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_04", "ply"),
              ESfM_Data(EXTRINSICS | STRUCTURE));
          }
    }

    return b_BA_Status;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_threecamera_gps_cail_new
(
        PoseStatusMap &poseStatusMap,
        const std::vector<double>& transformation_R_imu,
        const std::vector<double>& tramsformation_x_gps,
        const std::vector<double>& transformation_br,
        const std::vector<double>& transformation_fr,
        Hash_Map<IndexT, std::vector<double> >& C_gps_Map,
        Hash_Map<IndexT, std::vector<double> >& R_imu_Map)
{

    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_00", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));

    std::vector<double> out_transformation_R_imu = transformation_R_imu;
    std::vector<double> out_tramsformation_x_gps = tramsformation_x_gps;
    std::vector<double> out_transformation_br = transformation_br;
    std::vector<double> out_transformation_fr = transformation_fr;

    Bundle_Adjustment_Ceres bundle_adjustment_obj;
    bool b_BA_Status;
    b_BA_Status = bundle_adjustment_obj.Adjust_threecamera_gps_cail_first_cail_new_c1
      (
        sfm_data_,
        Optimize_Options(
        Intrinsic_Parameter_Type::NONE, //econstructionEngine::intrinsic_refinement_options_,//
        Extrinsic_Parameter_Type::NONE, //
        Structure_Parameter_Type::ADJUST_ALL),//ADJUST_ALL
        poseStatusMap,
        out_transformation_R_imu,
        out_tramsformation_x_gps,
        out_transformation_br,
        out_transformation_fr,
        C_gps_Map,
        R_imu_Map
      );
    if (b_BA_Status)
    {
      if (!sLogging_file_.empty())
      {
        Save(sfm_data_,
         stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_01", "ply"),
         ESfM_Data(EXTRINSICS | STRUCTURE));
      }

      b_BA_Status = bundle_adjustment_obj.Adjust_threecamera_gps_cail_first_cail_new_c1
        (
          sfm_data_,
          Optimize_Options(
          Intrinsic_Parameter_Type::NONE, //econstructionEngine::intrinsic_refinement_options_,//
          Extrinsic_Parameter_Type::ADJUST_ALL, //
          Structure_Parameter_Type::ADJUST_ALL),//ADJUST_ALL
          poseStatusMap,
          out_transformation_R_imu,
          out_tramsformation_x_gps,
          out_transformation_br,
          out_transformation_fr,
          C_gps_Map,
          R_imu_Map
        );

      if (!sLogging_file_.empty() && b_BA_Status)
      {
        Save(sfm_data_,
         stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_02", "ply"),
         ESfM_Data(EXTRINSICS | STRUCTURE));
      }

      // Remove outliers (max_angle, residual error)
      const size_t pointcount_initial = sfm_data_.structure.size();
      RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
      const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
      RemoveOutliers_AngleError(sfm_data_, 2.0);
      const size_t pointcount_angular_filter = sfm_data_.structure.size();
      std::cout << "Outlier removal (remaining #points):\n"
        << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
        << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
        << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

      if (!sLogging_file_.empty())
      {
        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_03", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));
      }


      //return b_BA_Status;

      b_BA_Status = bundle_adjustment_obj.Adjust_threecamera_gps_cail_first_cail_new_c1
        (
          sfm_data_,
          Optimize_Options(
          Intrinsic_Parameter_Type::NONE, //econstructionEngine::intrinsic_refinement_options_,//
          Extrinsic_Parameter_Type::ADJUST_ALL, //
          Structure_Parameter_Type::ADJUST_ALL),//ADJUST_ALL
          poseStatusMap,
          out_transformation_R_imu,
          out_tramsformation_x_gps,
          out_transformation_br,
          out_transformation_fr,
          C_gps_Map,
          R_imu_Map
        );
          if (b_BA_Status && !sLogging_file_.empty())
          {
            Save(sfm_data_,
              stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_04", "ply"),
              ESfM_Data(EXTRINSICS | STRUCTURE));
          }
    }

    return b_BA_Status;
}


bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_threecamera_gps_change_c4()
{

    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_00", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));

    std::vector<double> out_transformation_R_imu(3);
    std::vector<double> out_tramsformation_x_gps(3);
    std::vector<double> out_transformation_br(6);
    std::vector<double> out_transformation_fr(6);

    Bundle_Adjustment_Ceres bundle_adjustment_obj;
    bool b_BA_Status;
    b_BA_Status = bundle_adjustment_obj.Adjust_threecamera_gps_cail_first_change_c4
      (
        sfm_data_,
        Optimize_Options(
        Intrinsic_Parameter_Type::NONE, //econstructionEngine::intrinsic_refinement_options_,//
        Extrinsic_Parameter_Type::ADJUST_ALL, //
        Structure_Parameter_Type::NONE),//ADJUST_ALL
        out_transformation_R_imu,
        out_tramsformation_x_gps,
        out_transformation_br,
        out_transformation_fr
      );
    if (b_BA_Status)
    {
      if (!sLogging_file_.empty())
      {
        Save(sfm_data_,
         stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_01", "ply"),
         ESfM_Data(EXTRINSICS | STRUCTURE));
      }

      b_BA_Status = bundle_adjustment_obj.Adjust_threecamera_gps_cail_second_change_c4
        (
          sfm_data_,
          Optimize_Options(
          Intrinsic_Parameter_Type::NONE, //
          Extrinsic_Parameter_Type::ADJUST_ALL, //
          Structure_Parameter_Type::ADJUST_ALL),
          out_transformation_R_imu,
          out_tramsformation_x_gps,
          out_transformation_br,
          out_transformation_fr
        );

      if (!sLogging_file_.empty() && b_BA_Status)
      {
        Save(sfm_data_,
         stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_02", "ply"),
         ESfM_Data(EXTRINSICS | STRUCTURE));
      }

      // Remove outliers (max_angle, residual error)
      const size_t pointcount_initial = sfm_data_.structure.size();
      RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
      const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
      RemoveOutliers_AngleError(sfm_data_, 2.0);
      const size_t pointcount_angular_filter = sfm_data_.structure.size();
      std::cout << "Outlier removal (remaining #points):\n"
        << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
        << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
        << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

      if (!sLogging_file_.empty())
      {
        Save(sfm_data_,
          stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_03", "ply"),
          ESfM_Data(EXTRINSICS | STRUCTURE));
      }


      //return b_BA_Status;

      b_BA_Status = bundle_adjustment_obj.Adjust_threecamera_gps_cail_second_change_c4
        (
          sfm_data_,
          Optimize_Options(
            Intrinsic_Parameter_Type::NONE, //ReconstructionEngine::intrinsic_refinement_options_,  //Intrinsic_Parameter_Type::NONE, //
            Extrinsic_Parameter_Type::ADJUST_ALL,
            Structure_Parameter_Type::ADJUST_ALL),
            out_transformation_R_imu,
            out_tramsformation_x_gps,
            out_transformation_br,
            out_transformation_fr
        );
          if (b_BA_Status && !sLogging_file_.empty())
          {
            Save(sfm_data_,
              stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_use_test_04", "ply"),
              ESfM_Data(EXTRINSICS | STRUCTURE));
          }
    }

    return b_BA_Status;
}


void GlobalSfMReconstructionEngine_RelativeMotions::Compute_Relative_Rotations
(
  rotation_averaging::RelativeRotations & vec_relatives_R
)
{
  //
  // Build the Relative pose graph from matches:
  //
  /// pairwise view relation between poseIds
  using PoseWiseMatches = std::map< Pair, Pair_Set >;

  // List shared correspondences (pairs) between poses
  PoseWiseMatches poseWiseMatches;
  for (const auto & iterMatches : matches_provider_->pairWise_matches_)
  {
    const Pair pair = iterMatches.first;
    const View * v1 = sfm_data_.GetViews().at(pair.first).get();
    const View * v2 = sfm_data_.GetViews().at(pair.second).get();
    poseWiseMatches[Pair(v1->id_pose, v2->id_pose)].insert(pair);
  }

  C_Progress_display my_progress_bar( poseWiseMatches.size(),
      std::cout, "\n- Relative pose computation -\n" );

#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif
  // Compute the relative pose from pairwise point matches:
  for (int i = 0; i < static_cast<int>(poseWiseMatches.size()); ++i)
  {
    ++my_progress_bar;
    {
      PoseWiseMatches::const_iterator iter (poseWiseMatches.begin());
      std::advance(iter, i);
      const auto & relative_pose_iterator(*iter);
      const Pair relative_pose_pair = relative_pose_iterator.first;
      const Pair_Set & match_pairs = relative_pose_iterator.second;

      // If a pair has the same ID, discard it
      if (relative_pose_pair.first == relative_pose_pair.second)
      {
        continue;
      }

      // Select common bearing vectors
      if (match_pairs.size() > 1)
      {
        std::cerr << "Compute relative pose between more than two view is not supported" << std::endl;
        continue;
      }

      const Pair pairIterator = *(match_pairs.begin());

      const IndexT I = pairIterator.first;
      const IndexT J = pairIterator.second;

      const View * view_I = sfm_data_.views[I].get();
      const View * view_J = sfm_data_.views[J].get();

      // Check that valid cameras are existing for the pair of view
      if (sfm_data_.GetIntrinsics().count(view_I->id_intrinsic) == 0 ||
        sfm_data_.GetIntrinsics().count(view_J->id_intrinsic) == 0)
        continue;


      const IntrinsicBase * cam_I = sfm_data_.GetIntrinsics().at(view_I->id_intrinsic).get();
      const IntrinsicBase * cam_J = sfm_data_.GetIntrinsics().at(view_J->id_intrinsic).get();

      // Setup corresponding bearing vector
      const matching::IndMatches & matches = matches_provider_->pairWise_matches_.at(pairIterator);
      size_t nbBearing = matches.size();
      Mat x1(2, nbBearing), x2(2, nbBearing);
      nbBearing = 0;
      for (const auto & match : matches)
      {
        x1.col(nbBearing) = ((*cam_I)(cam_I->get_ud_pixel(features_provider_->feats_per_view[I][match.i_].coords().cast<double>()))).hnormalized();
        x2.col(nbBearing++) = ((*cam_J)(cam_J->get_ud_pixel(features_provider_->feats_per_view[J][match.j_].coords().cast<double>()))).hnormalized();
      }

      RelativePose_Info relativePose_info;
      // Compute max authorized error as geometric mean of camera plane tolerated residual error
      relativePose_info.initial_residual_tolerance = std::pow(
        cam_I->imagePlane_toCameraPlaneError(2.5) *
        cam_J->imagePlane_toCameraPlaneError(2.5),
        1./2.);

      // Since we use normalized features, we will use unit image size and intrinsic matrix:
      const std::pair<size_t, size_t> imageSize(1., 1.);
      const Mat3 K  = Mat3::Identity();

      if (!robustRelativePose(K, K, x1, x2, relativePose_info, imageSize, imageSize, 256))
      {
        continue;
      }
      const bool bRefine_using_BA = true;
      if (bRefine_using_BA)
      {
        // Refine the defined scene
        SfM_Data tiny_scene;
        tiny_scene.views.insert(*sfm_data_.GetViews().find(view_I->id_view));
        tiny_scene.views.insert(*sfm_data_.GetViews().find(view_J->id_view));
        tiny_scene.intrinsics.insert(*sfm_data_.GetIntrinsics().find(view_I->id_intrinsic));
        tiny_scene.intrinsics.insert(*sfm_data_.GetIntrinsics().find(view_J->id_intrinsic));

        // Init poses
        const Pose3 & Pose_I = tiny_scene.poses[view_I->id_pose] = Pose3(Mat3::Identity(), Vec3::Zero());
        const Pose3 & Pose_J = tiny_scene.poses[view_J->id_pose] = relativePose_info.relativePose;

        // Init structure
        const Mat34 P1 = cam_I->get_projective_equivalent(Pose_I);
        const Mat34 P2 = cam_J->get_projective_equivalent(Pose_J);
        Landmarks & landmarks = tiny_scene.structure;
        for (Mat::Index k = 0; k < x1.cols(); ++k)
        {
          const Vec2 x1_ = features_provider_->feats_per_view[I][matches[k].i_].coords().cast<double>();
          const Vec2 x2_ = features_provider_->feats_per_view[J][matches[k].j_].coords().cast<double>();
          Vec3 X;
          TriangulateDLT(P1, x1_.homogeneous(), P2, x2_.homogeneous(), &X);
          Observations obs;
          obs[view_I->id_view] = Observation(x1_, matches[k].i_);
          obs[view_J->id_view] = Observation(x2_, matches[k].j_);
          landmarks[k].obs = obs;
          landmarks[k].X = X;
        }
        // - refine only Structure and Rotations & translations (keep intrinsic constant)
        Bundle_Adjustment_Ceres::BA_Ceres_options options(false, false);
        options.linear_solver_type_ = ceres::DENSE_SCHUR;
        Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
        const Optimize_Options ba_refine_options
          (Intrinsic_Parameter_Type::NONE, // -> Keep intrinsic constant
          Extrinsic_Parameter_Type::ADJUST_ALL, // adjust camera motion
          Structure_Parameter_Type::ADJUST_ALL);// adjust scene structure
        if (bundle_adjustment_obj.Adjust(tiny_scene, ba_refine_options))
        {
          // --> to debug: save relative pair geometry on disk
          // std::ostringstream os;
          // os << relative_pose_pair.first << "_" << relative_pose_pair.second << ".ply";
          // Save(tiny_scene, os.str(), ESfM_Data(STRUCTURE | EXTRINSICS));
          //
          const Mat3 R1 = tiny_scene.poses[view_I->id_pose].rotation();
          const Mat3 R2 = tiny_scene.poses[view_J->id_pose].rotation();
          const Vec3 t1 = tiny_scene.poses[view_I->id_pose].translation();
          const Vec3 t2 = tiny_scene.poses[view_J->id_pose].translation();
          // Compute relative motion and save it
          Mat3 Rrel;
          Vec3 trel;
          RelativeCameraMotion(R1, t1, R2, t2, &Rrel, &trel);
          // Update found relative pose
          relativePose_info.relativePose = Pose3(Rrel, -Rrel.transpose() * trel);
        }
      }
#ifdef OPENMVG_USE_OPENMP
      #pragma omp critical
#endif
      {
        // Add the relative rotation to the relative 'rotation' pose graph
        using namespace openMVG::rotation_averaging;
          vec_relatives_R.emplace_back(
            relative_pose_pair.first, relative_pose_pair.second,
            relativePose_info.relativePose.rotation(),
            1.f);
      }
    }
  } // for all relative pose

  // Log input graph to the HTML report
  if (!sLogging_file_.empty() && !sOut_directory_.empty())
  {
    // Log a relative view graph
    {
//      std::set<IndexT> set_ViewIds;
//      std::transform(sfm_data_.GetViews().begin(), sfm_data_.GetViews().end(),
//        std::inserter(set_ViewIds, set_ViewIds.begin()), stl::RetrieveKey());
//      graph::indexedGraph putativeGraph(set_ViewIds, getPairs(matches_provider_->pairWise_matches_));
//      graph::exportToGraphvizData(
//        stlplus::create_filespec(sOut_directory_, "global_relative_rotation_view_graph"),
//        putativeGraph);
    }

    // Log a relative pose graph
    {
      std::set<IndexT> set_pose_ids;
      Pair_Set relative_pose_pairs;
      for (const auto & relative_R : vec_relatives_R)
      {
        const Pair relative_pose_indices(relative_R.i, relative_R.j);
        relative_pose_pairs.insert(relative_pose_indices);
        set_pose_ids.insert(relative_R.i);
        set_pose_ids.insert(relative_R.j);
      }
//      const std::string sGraph_name = "global_relative_rotation_pose_graph";
//      graph::indexedGraph putativeGraph(set_pose_ids, relative_pose_pairs);
//      graph::exportToGraphvizData(
//        stlplus::create_filespec(sOut_directory_, sGraph_name),
//        putativeGraph);
      using namespace htmlDocument;
      std::ostringstream os;

      os << "<br>" << "global_relative_rotation_pose_graph" << "<br>"
         << "<img src=\""
         << stlplus::create_filespec(sOut_directory_, "global_relative_rotation_pose_graph", "svg")
         << "\" height=\"600\">\n";

      html_doc_stream_->pushInfo(os.str());
    }
  }
}

} // namespace sfm
} // namespace openMVG
