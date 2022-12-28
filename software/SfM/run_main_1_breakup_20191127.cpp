// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/cameras/Cameras_Common_command_line_helper.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"
#include "openMVG/sfm/pipelines/global/sfm_global_engine_relative_motions.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_report.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/features/feature.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <memory>
#include <string>
#include <list>

//
#include <ceres/ceres.h>
#include "ceres/problem.h"
#include "ceres/solver.h"
#include <ceres/rotation.h>
#include <ceres/types.h>



using namespace openMVG;
using namespace openMVG::sfm;



///// Define a collection of landmarks are indexed by their TrackId
//using Landmarks = Hash_Map<IndexT, Landmark>;


void load_matches(const std::__cxx11::string file_path, Landmarks &inputLM);


struct ResidualErrorFunctor_Pinhole_Intrinsic_WaterXY
{
  ResidualErrorFunctor_Pinhole_Intrinsic_WaterXY(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum : uint8_t {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2
  };

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y] )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_extrinsics;
    const T * cam_C = &cam_extrinsics[3];

    //prepare l
    T R[9];
    ceres::AngleAxisToRotationMatrix(cam_R, R);

    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];

    T x1_c[3];
    x1_c[0] = (m_pos_2dpoint[0] - principal_point_x) /focal;
    x1_c[1] = (m_pos_2dpoint[1] - principal_point_y) /focal;
    x1_c[2] = (T) 1.0;

    const T Xrep_x = R[0]*x1_c[0] + R[1]*x1_c[1] + R[2] + cam_C[0];
    const T Xrep_y = R[3]*x1_c[0] + R[4]*x1_c[1] + R[5] + cam_C[1];
    //const T Xrep_z = R[6]*x1_c0 + R[7]*x1_c1 + R[8] + cam_C[2];


    const T L_XC_0 = cam_C[1] - pos_3dpoint[1];
    const T L_XC_1 = pos_3dpoint[0] - cam_C[0];
    const T L_XC_2 = cam_C[0]*pos_3dpoint[1] - cam_C[1]*pos_3dpoint[0];

    //const T d_Xrep_L = ceres::abs(Xrep_x*L_XC_0 + Xrep_y*L_XC_1 + L_XC_2);// / ceres::sqrt((L_XC_0*L_XC_0 + L_XC_1*L_XC_1));
    const T d_Xrep_L = ceres::abs(Xrep_x*L_XC_0 + Xrep_y*L_XC_1 + L_XC_2) / ceres::sqrt((L_XC_0*L_XC_0 + L_XC_1*L_XC_1));

    out_residuals[0] = d_Xrep_L;

    return true;
  }

  static int num_residuals() { return 1; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {

    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_WaterXY, 1, 3, 5, 2>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_WaterXY(observation.data())));
    }

  }

  const double * m_pos_2dpoint; // The 2D observation
};


void random_select(const SfM_Data &sfm_data, const std::shared_ptr<Matches_Provider> &matches_provider,
                   Views &selected_views, Matches_Provider &selected_matches);
void get_newX(const Landmarks &matches_LM, SfM_Data &sfm_data_init, std::vector<bool> &used_X_id_list);
void calculate_X(const std::vector<std::pair<Pose3, Vec2> > &pose_x_list, Vec3 &X_res );
int get_newC(const Landmarks &matches_LM, SfM_Data &sfm_data_init);
void calculate_C(const std::vector<std::pair<Vec3, Vec2> > &X_x_list, Pose3 &camera_res );

void run_process_horizontal(const SfM_Data &sfm_data_input, const std::shared_ptr<Matches_Provider> &matches_input_provider, const std::shared_ptr<Features_Provider> &features_input_provider,
                            const std::string sOutDir, SfM_Data &sfm_data_output);
bool optimisic_operation(const SfM_Data &sfm_data_input, SfM_Data &sfm_data_output);

int main(int argc, char **argv)
{
    //bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust_water_xy_z_n()

    CmdLine cmd;

    std::string sSfM_Data_Filename;
    std::string sMatchesDir;
    std::string sOutDir = "";
//    int iRotationAveragingMethod = int (ROTATION_AVERAGING_L2);
//    int iTranslationAveragingMethod = int (TRANSLATION_AVERAGING_SOFTL1);
//    std::string sIntrinsic_refinement_options = "ADJUST_ALL";
    bool b_use_motion_priors = false;

    cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
    cmd.add( make_option('m', sMatchesDir, "matchdir") );
    cmd.add( make_option('o', sOutDir, "outdir") );
//    cmd.add( make_option('r', iRotationAveragingMethod, "rotationAveraging") );
//    cmd.add( make_option('t', iTranslationAveragingMethod, "translationAveraging") );
//    cmd.add( make_option('f', sIntrinsic_refinement_options, "refineIntrinsics") );
//    cmd.add( make_switch('P', "prior_usage") );

    try {
      if (argc == 1) throw std::string("Invalid parameter.");
      cmd.process(argc, argv);
    } catch (const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--input_file] path to a SfM_Data scene\n"
      << "[-m|--matchdir] path to the matches that corresponds to the provided SfM_Data scene\n"
      << "[-o|--outdir] path where the output data will be stored\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
    }

    /// ***
    /// Load all data
    /// ***

    // Load input SfM_Data scene
    SfM_Data sfm_data;
    if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS))) {
      std::cerr << std::endl
        << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
      return EXIT_FAILURE;
    }

//    return EXIT_SUCCESS;

    // Init the regions_type from the image describer file (used for image regions extraction)
    using namespace openMVG::features;
    const std::string sImage_describer = stlplus::create_filespec(sMatchesDir, "image_describer", "json");
    std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
    if (!regions_type)
    {
      std::cerr << "Invalid: "
        << sImage_describer << " regions type file." << std::endl;
      return EXIT_FAILURE;
    }

    // Features reading
    std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
    if (!feats_provider->load(sfm_data, sMatchesDir, regions_type)) {
      std::cerr << std::endl
        << "Invalid features." << std::endl;
      return EXIT_FAILURE;
    }
    // Matches reading
    std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
    if // Try to read the two matches file formats
    (
      !(matches_provider->load(sfm_data, stlplus::create_filespec(sMatchesDir, "matches.e.txt")) ||
        matches_provider->load(sfm_data, stlplus::create_filespec(sMatchesDir, "matches.e.bin")))
    )
    {
      std::cerr << std::endl
        << "Invalid matches file." << std::endl;
      return EXIT_FAILURE;
    }

//    if (sOutDir.empty())  {
//      std::cerr << "\nIt is an invalid output directory" << std::endl;
//      return EXIT_FAILURE;
//    }

//    if (!stlplus::folder_exists(sOutDir))
//    {
//      if (!stlplus::folder_create(sOutDir))
//      {
//        std::cerr << "\nCannot create the output directory" << std::endl;
//      }
//    }




    Landmarks landmarks = sfm_data.structure;
    //load_matches(matches_file_path, landmarks);


    SfM_Data sfm_data_after_horizontal;
    run_process_horizontal(sfm_data, matches_provider, feats_provider, sOutDir, sfm_data_after_horizontal);




    Poses poses;
    int camera_num;

}

/// ***
/// select 3 cameras randomly
/// ***
void random_select(const SfM_Data &sfm_data, const std::shared_ptr<Matches_Provider> &matches_provider,
                   Views &selected_views, Matches_Provider &selected_matches)
{
    int camera_count = sfm_data.views.size();

    std::srand(time(NULL));
    int cId1 = rand() % (camera_count);

    std::srand(time(NULL));
//    int cId0 = rand() % (cId1);
    int cId0 = 0;//rand() % (cId1);

    std::srand(time(NULL));
    int cId2 = (rand() % (camera_count-1-cId1) ) + cId1 + 1 ;

    //
//    selected_views[cId0] = sfm_data.views[cId0].second;
//    selected_views[cId1] = sfm_data.views[cId1];
//    selected_views[cId2] = sfm_data.views[cId2];

    //
    std::vector<int> setected_cId_list;
    setected_cId_list.push_back(cId0);
    setected_cId_list.push_back(cId1);
    setected_cId_list.push_back(cId2);

    //
    for (const auto & iterMatches : matches_provider->pairWise_matches_)
    {
      const Pair pair = iterMatches.first;
      const int c1 = pair.first;
      const int c2 = pair.second;

      std::vector<int>::iterator find_c1 = std::find(setected_cId_list.begin(), setected_cId_list.end(), c1);
      if(find_c1 != setected_cId_list.end())
      {
          std::vector<int>::iterator find_c2 = std::find(setected_cId_list.begin(), setected_cId_list.end(), c2);
          if(find_c2 != setected_cId_list.end())
          {
              selected_matches.pairWise_matches_.insert(iterMatches);
          }
      }

    }

}

/// ***
/// Get init sub-scene
/// ***
bool get_init_subscene(const SfM_Data &sfm_data, const std::shared_ptr<Matches_Provider> &matches_provider, const std::shared_ptr<Features_Provider> &features_provider,
                       const std::string sOutDir,
                       SfM_Data &sfm_data_res)
{
    sfm_data_res = sfm_data;

    //openMVG::system::Timer timer;
    GlobalSfMReconstructionEngine_RelativeMotions sfmEngine(
      sfm_data_res,
      sOutDir+"/temp",
      stlplus::create_filespec(sOutDir+"/temp", "Reconstruction_Report.html"));

    // Configure the features_provider & the matches_provider
    sfmEngine.SetFeaturesProvider(features_provider.get());
    sfmEngine.SetMatchesProvider(matches_provider.get());

    const cameras::Intrinsic_Parameter_Type intrinsic_refinement_options =
      cameras::StringTo_Intrinsic_Parameter_Type("NONE");

    // Configure reconstruction parameters
    sfmEngine.Set_Intrinsics_Refinement_Type(intrinsic_refinement_options);
    bool b_use_motion_priors = false;
    sfmEngine.Set_Use_Motion_Prior(b_use_motion_priors);

    // Configure motion averaging method
    sfmEngine.SetRotationAveragingMethod(
      ERotationAveragingMethod(ROTATION_AVERAGING_L2));
    sfmEngine.SetTranslationAveragingMethod(
      ETranslationAveragingMethod(TRANSLATION_AVERAGING_SOFTL1));

  //  if (sfmEngine.Process_threeCameras_imu_gps_init_new(rtFilePath))//Process_threeCameras_change_step(rtFilePath, refSfmDataIdList))//Process_threeCameras_imu_gps_init_new())//Process_test(sMatchesDir))//_down_imu_gps(sMatchesDir))
    if(sfmEngine.Process())
    {
        return true;
    }else{
        return false;
    }

}

/// ***
/// Get the whole structure
/// ***
void get_tracks(const SfM_Data &sfm_data, const std::shared_ptr<Matches_Provider> &matches_provider, const std::shared_ptr<Features_Provider> &features_provider,
                Landmarks &structure_whole)
{
    using namespace openMVG::tracks;
    TracksBuilder tracksBuilder;
//#if defined USE_ALL_VALID_MATCHES // not used by default
    matching::PairWiseMatches pose_supported_matches;

    //for (const std::pair< Pair, matching::IndMatches > & match_info :  matches_provider_->pairWise_matches_)
    for (const auto & match_info :  matches_provider->pairWise_matches_)
    {
      const View * vI = sfm_data.GetViews().at(match_info.first.first).get();
      const View * vJ = sfm_data.GetViews().at(match_info.first.second).get();
      //if (sfm_data.IsPoseAndIntrinsicDefined(vI) && sfm_data.IsPoseAndIntrinsicDefined(vJ))
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
    Landmarks & structure = structure_whole;
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
        const features::PointFeature & pt = features_provider->feats_per_view.at(imaIndex)[featIndex];
        obs[imaIndex] = Observation(pt.coords().cast<double>(), featIndex);
      }
    }

}

void load_matches(const std::string file_path, Landmarks &inputLM)
{
    std::fstream in;
    in.open(file_path.c_str());
    if(!in)
    {
        std::cout << "file " << file_path << " opened error!" << std::endl;

        }

    int lm_count;
    in >> lm_count;

    int this_lm = 0;
    while(this_lm < lm_count || !in.eof())
    {
        int x_count = 0;
        in >> x_count;

        Landmark thisLM;

        Observations obs_list;

        int this_x = 0;
        while(this_x < x_count)
        {
            double camera_id, feature_id, x_x, x_y;
            in >> camera_id >> feature_id >> x_x >> x_y;


            Observation thisObs;
            thisObs.x = Vec2(x_x, x_y);
            thisObs.id_feat = feature_id;

            obs_list[camera_id] = thisObs;

        }

        thisLM.obs = obs_list;

        //
        inputLM[this_lm] = thisLM;



        ++this_lm;

    }

    in.close();


}

/// ***
/// The calcuation process of solving Z
/// ***
void run_process_horizontal(const SfM_Data &sfm_data_input, const std::shared_ptr<Matches_Provider> &matches_input_provider, const std::shared_ptr<Features_Provider> &features_input_provider,
                            const std::string sOutDir, SfM_Data &sfm_data_output)
{
//    std::vector<int> used_camera_id_list;
//    int new_X_num, new_C_num;


    Landmarks structure_whole;
    get_tracks(sfm_data_input, matches_input_provider, features_input_provider, structure_whole);

    std::vector<bool> used_X_id_list(structure_whole.size(), false);
    int X_threshold = structure_whole.size()/4*3;

    SfM_Data sfm_data_prepared;
    int all_X_number = 0;
    while(all_X_number < X_threshold)//new_C_num != 0 && used_camera_id_list.size()  )//!=  ) // camera_count > ... && number_X > ...
    {

        //
        Views selected_views;
        Matches_Provider selected_matches;
        random_select(sfm_data_input, matches_input_provider, selected_views, selected_matches);

        std::shared_ptr<Matches_Provider> selected_matches_ptr(new Matches_Provider(selected_matches));// = new Matches_Provider(selected_matches);


        // get init subscene
        SfM_Data sfm_data_init;              

        if(!get_init_subscene(sfm_data_input, selected_matches_ptr, features_input_provider, sOutDir, sfm_data_init))
        {
            continue;
        }

        //
        int new_camera_number = 0;

        //
        do{
            // X
            get_newX(structure_whole, sfm_data_init, used_X_id_list);
            // C
            new_camera_number = get_newC(structure_whole, sfm_data_init);
            //

        }while(new_camera_number!=0);

        //
        all_X_number = 0;
        for(const auto &flagX : used_X_id_list)
        {
            if(flagX == true)
            {
                ++all_X_number;
            }

        }

        sfm_data_prepared = sfm_data_init;
//        sfm_data_prepared.s_root_path = sfm_data_init.s_root_path;
//        sfm_data_prepared.poses = sfm_data_init.poses;
//        sfm_data_prepared.structure = sfm_data_init.structure;



    }

    Save(sfm_data_output,
      stlplus::create_filespec(sOutDir, "cloud_and_poses_horizontal_before_opt", ".ply"),
      ESfM_Data(ALL));


    /// optimisic
    ///

    // clean structure
    Landmarks::iterator thisX = sfm_data_prepared.structure.begin();
    while(thisX != sfm_data_prepared.structure.end())
    {
        if(used_X_id_list[thisX->first] == false)
        {
            thisX = sfm_data_prepared.structure.erase(thisX);
        }else{
            ++thisX;
        }

    }

    // optimisic
    optimisic_operation(sfm_data_prepared, sfm_data_output);

    // save horizontal results
    Save(
        sfm_data_output,
        stlplus::create_filespec( sfm_data_output.s_root_path, "sfm_data_horizontal.json" ).c_str(),
        ESfM_Data(ALL));
    Save(sfm_data_output,
      stlplus::create_filespec(sOutDir, "cloud_and_poses_horizontal", ".ply"),
      ESfM_Data(ALL));


}

void get_init_data()
{
    const int C1_id = 0;
    const int X_num = 9;

    std::map<int, std::map<int, Vec2> > x_selected;

    //random_select


    //
    Pose3 C1;

    Mat3 H1;
    H1 << C1.rotation()(0,0), C1.rotation()(0,1), C1.translation()(0),
          C1.rotation()(1,0), C1.rotation()(1,1), C1.translation()(0),
          C1.rotation()(2,0), C1.rotation()(2,1), C1.translation()(0);

    Vec3 r3(C1.rotation()(0,2), C1.rotation()(1,2), C1.rotation()(2,2));
    Mat3 r3cross;
    r3cross << 0.0, -r3(2), r3(1),
               -r3(2), 0.0, -r3(0),
               r3(1), -r3(0), 0.0;

    Mat3 AH1;
    AH1 = H1.transpose()*r3cross;

    std::vector<Vec2> Xxy_list;

    const Vec2 this_x1 = x_selected[0][0];

    // in C1
    std::vector<Vec3> X0_X1_1_factor_list;
    for (std::size_t XId = 0; XId < X_num; ++XId)
    {
        Vec2 this_x = x_selected[0][XId];

        Vec3 X0_X1_1_factor = Vec3(this_x(0)*AH1(0,0) + this_x(1)*AH1(0,1) + AH1(0,2),
                                   this_x(0)*AH1(1,0) + this_x(1)*AH1(1,1) + AH1(1,2),
                                   this_x(0)*AH1(2,0) + this_x(1)*AH1(2,1) + AH1(2,2));

        X0_X1_1_factor_list.push_back(X0_X1_1_factor);

//        Vec2 this_X0_X1_factor( , );
//        Xxy_list.push_back();

    }








}

//void random_select(const int camera_num, const Landmarks &matches_LM,  std::vector<std::pair<std::pair<int,int>, Vec2> > &x_selected, //<X ID, Camera ID>
//                   std::vector<int> &XId_selected, std::vector<int> &CId_selected)
//{
//    //C-5, X-9
//    std::list<int> random_selected_list;



//    //

//    //
//}

void get_newX(const Landmarks &matches_LM, SfM_Data &sfm_data_init, std::vector<bool> &used_X_id_list)
{
    const Poses &poses = sfm_data_init.poses;

    std::vector<bool> record_solved_pose(sfm_data_init.views.size(), false);
    for(const auto & pose_it : poses)
    {
        const IndexT indexPose = pose_it.first;
        record_solved_pose[indexPose] = true;
    }


    Landmarks new_structure_forXxy = sfm_data_init.structure;

    // find and calculate new Xxy
    for(const auto & itLM : matches_LM)
    {
        const Observations obs_list = itLM.second.obs;

        int num_solved_c = 0;

        std::vector<std::pair<Pose3, Vec2> > pose_x_used;

        for (const auto & itObs : obs_list)
        {
            if(record_solved_pose[itObs.first] == true)
            {
                pose_x_used.push_back(std::pair<Pose3, Vec2>(poses.at(itObs.first), itObs.second.x));
                ++num_solved_c;

            }
        }

        if (num_solved_c > 2)
        {
            Vec3 X_res;
            calculate_X(pose_x_used, X_res);

            //std::map<int, Vec3> &newX
            new_structure_forXxy[itLM.first].X = X_res;
            used_X_id_list[itLM.first] = true;

        }

    }

    // update
    sfm_data_init.structure = new_structure_forXxy;

}

void calculate_X(const std::vector<std::pair<Pose3, Vec2> > &pose_x_list, Vec3 &X_res ) // X^T H^T [n_infy]_\times x = 0
{

    Eigen::MatrixXd A(pose_x_list.size(), 3);
    std::size_t listId = 0;
    for (const auto & itpx : pose_x_list)
    {
        Vec3 x_img = Vec3(itpx.second(0), itpx.second(1), 1.0);

        Pose3 thisPose = itpx.first;
        Mat3 thisR = thisPose.rotation();
        Vec3 thist = thisPose.translation();

        Vec3 n_infy = Vec3(thisR(0,2), thisR(1,2), thisR(2,2));
        Mat3 H;
        H << thisR(0,0), thisR(0,1), thist(0),
             thisR(1,0), thisR(1,1), thist(1),
             thisR(2,0), thisR(2,1), thist(2);


        Vec3 l_img = n_infy.cross(x_img);

        Vec3 HTl = H.transpose() * l_img;

        A(listId,0) = HTl(0);
        A(listId,1) = HTl(1);
        A(listId,2) = HTl(2);


        ++listId;

    }

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3> > eigen_sovler(A.transpose()*A);

    eigen_sovler.eigenvalues();
    eigen_sovler.eigenvectors();

    //
    X_res(0) = eigen_sovler.eigenvectors()(0,0);
    X_res(1) = eigen_sovler.eigenvectors()(1,0);
    X_res(2) = eigen_sovler.eigenvectors()(2,0);

}

int get_newC(const Landmarks &matches_LM, SfM_Data &sfm_data_init)//(const Landmarks matches_LM, const int camera_num, Poses &poses)
{
    const Poses &poses = sfm_data_init.poses;

    std::vector<bool> record_solved_pose(sfm_data_init.views.size(), false);
    for(const auto & pose_it : poses)
    {
        const IndexT indexPose = pose_it.first;
        record_solved_pose[indexPose] = true;
    }

    int new_pose_number = 0;

    //
    for (std::size_t pId = 0; pId < sfm_data_init.views.size(); ++pId)
    {
        if (record_solved_pose[pId] == false)
        {
            std::vector<std::pair<Vec3, Vec2> > X_x_used;
            int num_unsolved_c = 0;

            //
            for(const auto & itLM : matches_LM)
            {
                const Observations obs_list = itLM.second.obs;

                for (const auto & itObs : obs_list)
                {
                    if(itObs.first == pId)
                    {
                        X_x_used.push_back(std::pair<Vec3, Vec2>(itLM.second.X, itObs.second.x));
                        ++num_unsolved_c;

                    }
                }

            }

            //
            if (num_unsolved_c > 2)// ?? for save a camera
            {
                Pose3 C_res;
                calculate_C(X_x_used, C_res);

                sfm_data_init.poses[pId] = C_res;

                ++new_pose_number;

            }
        }
    }

    return new_pose_number;

}

void calculate_C(const std::vector<std::pair<Vec3, Vec2> > &X_x_list, Pose3 &camera_res )
{
    Eigen::MatrixXd AH(X_x_list.size(), 9);
    std::size_t listId = 0;
    for (const auto & itXx : X_x_list)
    {
        Vec3 x = Vec3(itXx.second(0), itXx.second(1), 1.0);

        Vec3 X = itXx.first;

        //
        AH(listId,0) = X(0)*x(0);
        AH(listId,1) = X(0)*x(1);
        AH(listId,2) = X(0);
        AH(listId,3) = X(1)*x(0);
        AH(listId,4) = X(1)*x(1);
        AH(listId,5) = X(1);
        AH(listId,6) = x(0);
        AH(listId,7) = x(1);
        AH(listId,8) = 1.0;

        ++listId;

    }

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 9, 9> > eigen_sovler(AH.transpose()*AH);

    eigen_sovler.eigenvalues();
    eigen_sovler.eigenvectors();

    Mat3 solved_res;
    solved_res << eigen_sovler.eigenvectors()(0,0), eigen_sovler.eigenvectors()(1,0), eigen_sovler.eigenvectors()(2,0),
                  eigen_sovler.eigenvectors()(3,0), eigen_sovler.eigenvectors()(4,0), eigen_sovler.eigenvectors()(5,0),
                  eigen_sovler.eigenvectors()(6,0), eigen_sovler.eigenvectors()(7,0), eigen_sovler.eigenvectors()(8,0);

    Vec3 r1(eigen_sovler.eigenvectors()(3,0), eigen_sovler.eigenvectors()(4,0), eigen_sovler.eigenvectors()(5,0));
    r1 = r1/(r1.norm());

    Vec3 r2(eigen_sovler.eigenvectors()(0,0), eigen_sovler.eigenvectors()(1,0), eigen_sovler.eigenvectors()(2,0));
    r2 = -r2/(r2.norm());

    Vec3 r3 = r1.cross(r2);

    Mat3 r3cross;
    r3cross << 0.0, -r3(2), r3(1),
               -r3(2), 0.0, -r3(0),
               r3(1), -r3(0), 0.0;

    Vec3 rAH3(eigen_sovler.eigenvectors()(6,0), eigen_sovler.eigenvectors()(7,0), eigen_sovler.eigenvectors()(8,0));
    Vec3 t_res = r3cross.inverse() * rAH3 / (rAH3.norm());

    Mat3 R_res;
    R_res << r1(0), r2(0), r3(0),
             r1(1), r2(1), r3(1),
             r1(2), r2(2), r3(2);
    camera_res = Pose3(R_res, t_res);



}

double calculate_residuals_forXY()
{
    double residuals = 0;

    //




    return residuals;

}

bool optimisic_operation(const SfM_Data &sfm_data_input, SfM_Data &sfm_data_output)
{    
    sfm_data_output = sfm_data_input;

    //----------
    // Add camera parameters
    // - intrinsics
    // - poses [R|t]

    // Create residuals for each observation in the bundle adjustment problem. The
    // parameters for cameras and points are added automatically.
    //----------

{
//    double pose_center_robust_fitting_error = 0.0;
//    openMVG::geometry::Similarity3 sim_to_center;
//    bool b_usable_prior = false;
//    if (options.use_motion_priors_opt && sfm_data.GetViews().size() > 3)
//    {
//      // - Compute a robust X-Y affine transformation & apply it
//      // - This early transformation enhance the conditionning (solution closer to the Prior coordinate system)
//      {
//        // Collect corresponding camera centers
//        std::vector<Vec3> X_SfM, X_GPS;
//        for (const auto & view_it : sfm_data.GetViews())
//        {
//          const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
//          if (prior != nullptr && prior->b_use_pose_center_ && sfm_data.IsPoseAndIntrinsicDefined(prior))
//          {
//            X_SfM.push_back( sfm_data.GetPoses().at(prior->id_pose).center() );
//            X_GPS.push_back( prior->pose_center_ );
//          }
//        }
//        openMVG::geometry::Similarity3 sim;

////        // Compute the registration:
////        if (X_GPS.size() > 3)
////        {
////          const Mat X_SfM_Mat = Eigen::Map<Mat>(X_SfM[0].data(),3, X_SfM.size());
////          const Mat X_GPS_Mat = Eigen::Map<Mat>(X_GPS[0].data(),3, X_GPS.size());
////          geometry::kernel::Similarity3_Kernel kernel(X_SfM_Mat, X_GPS_Mat);
////          const double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
////          if (lmeds_median != std::numeric_limits<double>::max())
////          {
////            b_usable_prior = true; // PRIOR can be used safely

////            // Compute the median residual error once the registration is applied
////            for (Vec3 & pos : X_SfM) // Transform SfM poses for residual computation
////            {
////              pos = sim(pos);
////            }
////            Vec residual = (Eigen::Map<Mat3X>(X_SfM[0].data(), 3, X_SfM.size()) - Eigen::Map<Mat3X>(X_GPS[0].data(), 3, X_GPS.size())).colwise().norm();
////            std::sort(residual.data(), residual.data() + residual.size());
////            pose_center_robust_fitting_error = residual(residual.size()/2);

////            // Apply the found transformation to the SfM Data Scene
////            openMVG::sfm::ApplySimilarity(sim, sfm_data);

////            // Move entire scene to center for better numerical stability
////            Vec3 pose_centroid = Vec3::Zero();
////            for (const auto & pose_it : sfm_data.poses)
////            {
////              pose_centroid += (pose_it.second.center() / (double)sfm_data.poses.size());
////            }
////            sim_to_center = openMVG::geometry::Similarity3(openMVG::sfm::Pose3(Mat3::Identity(), pose_centroid), 1.0);
////            openMVG::sfm::ApplySimilarity(sim_to_center, sfm_data, true);
////          }
////        }
//      }
//    }
}
    ceres::Problem problem;

    // Data wrapper for refinement:
    Hash_Map<IndexT, std::vector<double> > map_intrinsics;
    Hash_Map<IndexT, std::vector<double> > map_poses;

    // Setup Poses data & subparametrization
    for (const auto & pose_it : sfm_data_input.poses)
    {
      const IndexT indexPose = pose_it.first;

      const Pose3 & pose = pose_it.second;
      const Mat3 R = pose.rotation();
      const Vec3 t = pose.translation();

      double angleAxis[3];
      ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
      // angleAxis + translation
      map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

      double * parameter_block = &map_poses[indexPose][0];
      problem.AddParameterBlock(parameter_block, 6);

      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);

//      if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
//      {
//        // set the whole parameter block as constant for best performance
//        problem.SetParameterBlockConstant(parameter_block);
//      }
//      else  // Subset parametrization
//      {
//        std::vector<int> vec_constant_extrinsic;
//        // If we adjust only the translation, we must set ROTATION as constant
//        if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
//        {
//          // Subset rotation parametrization
//          vec_constant_extrinsic.push_back(0);
//          vec_constant_extrinsic.push_back(1);
//          vec_constant_extrinsic.push_back(2);
//        }
//        // If we adjust only the rotation, we must set TRANSLATION as constant
//        if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
//        {
//          // Subset translation parametrization
//          vec_constant_extrinsic.push_back(3);
//          vec_constant_extrinsic.push_back(4);
//          vec_constant_extrinsic.push_back(5);
//        }
//        if (!vec_constant_extrinsic.empty())
//        {
//          ceres::SubsetParameterization *subset_parameterization =
//            new ceres::SubsetParameterization(6, vec_constant_extrinsic);
//          problem.SetParameterization(parameter_block, subset_parameterization);
//        }
//      }
    }

    // Setup Intrinsics data & subparametrization
    for (const auto & intrinsic_it : sfm_data_input.intrinsics)
    {
      const IndexT indexCam = intrinsic_it.first;

      if (isValid(intrinsic_it.second->getType()))
      {
        map_intrinsics[indexCam] = intrinsic_it.second->getParams();
        if (!map_intrinsics[indexCam].empty())
        {
          double * parameter_block = &map_intrinsics[indexCam][0];
          problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());

          // set the whole parameter block as constant for best performance
          problem.SetParameterBlockConstant(parameter_block);

//          if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
//          {
//            // set the whole parameter block as constant for best performance
//            problem.SetParameterBlockConstant(parameter_block);
//          }
//          else
//          {
//            const std::vector<int> vec_constant_intrinsic =
//              intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
//            if (!vec_constant_intrinsic.empty())
//            {
//              ceres::SubsetParameterization *subset_parameterization =
//                new ceres::SubsetParameterization(
//                  map_intrinsics[indexCam].size(), vec_constant_intrinsic);
//              problem.SetParameterization(parameter_block, subset_parameterization);
//            }
//          }

        }
      }
      else
      {
        std::cerr << "Unsupported camera type." << std::endl;
      }
    }

    // Set a LossFunction to be less penalized by false measurements
    //  - set it to nullptr if you don't want use a lossFunction.
    ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(Square(4.0)); //nullptr;

    // For all visibility add reprojections errors:
    for (auto & structure_landmark_it : sfm_data_output.structure)
    {
      const Observations & obs = structure_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data_input.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        // !!!
        ceres::CostFunction* cost_function =
                ResidualErrorFunctor_Pinhole_Intrinsic_WaterXY::Create(obs_it.second.x);


        if (cost_function)
        {
          if (!map_intrinsics[view->id_intrinsic].empty())
          {
            problem.AddResidualBlock(cost_function,
              p_LossFunction,
              &map_intrinsics[view->id_intrinsic][0],
              &map_poses[view->id_pose][0],
              structure_landmark_it.second.X.data());
          }
          else
          {
            problem.AddResidualBlock(cost_function,
              p_LossFunction,
              &map_poses[view->id_pose][0],
              structure_landmark_it.second.X.data());
          }
        }
      }
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());

    }

    //
    // Configure a BA engine and run it
    //  Make Ceres automatically detect the bundle structure.
  //  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
  //  linear_solver_type_ = ceres::SPARSE_SCHUR;
  //  linear_solver_type_ = ceres::DENSE_SCHUR;
  //  preconditioner_type_ = ceres::JACOBI;

    ceres::Solver::Options ceres_config_options;
    ceres_config_options.max_num_iterations = 500;
    ceres_config_options.preconditioner_type = ceres::JACOBI;
      //static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
    ceres_config_options.linear_solver_type =  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
      //static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
    ceres_config_options.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
      //static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
    ceres_config_options.minimizer_progress_to_stdout = 500;//ceres_options_.bVerbose_;
    ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
    ceres_config_options.num_threads = 6;//ceres_options_.nb_threads_;
    ceres_config_options.num_linear_solver_threads = 6;//ceres_options_.nb_threads_;
    ceres_config_options.parameter_tolerance = 0.000000001;//ceres_options_.parameter_tolerance_;

  //  std::cout << "start ceres solver" << std::endl;
    // Solve BA
    ceres::Solver::Summary summary;
    ceres::Solve(ceres_config_options, &problem, &summary);
//    if (ceres_options_.bCeres_summary_)
//      std::cout << summary.FullReport() << std::endl;

    // If no error, get back refined parameters
    if (!summary.IsSolutionUsable())
    {
//      if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
      return false;
    }
    else // Solution is usable
    {
      //if (ceres_options_.bVerbose_)
      {
        // Display statistics about the minimization
        std::cout << std::endl
          << "Bundle Adjustment statistics (approximated RMSE):\n"
          << " #views: " << sfm_data_output.views.size() << "\n"
          << " #poses: " << sfm_data_output.poses.size() << "\n"
          << " #intrinsics: " << sfm_data_output.intrinsics.size() << "\n"
          << " #tracks: " << sfm_data_output.structure.size() << "\n"
          << " #residuals: " << summary.num_residuals << "\n"
          << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
          << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
          << " BriefReport : " << summary.BriefReport() << "\n"
          << " Time (s): " << summary.total_time_in_seconds << "\n"
          << std::endl;
      }

      // Update camera poses with refined data
      //if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
      {
        for (auto & pose_it : sfm_data_output.poses)
        {
          const IndexT indexPose = pose_it.first;

          Mat3 R_refined;
          ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
          Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
          // Update the pose
          Pose3 & pose = pose_it.second;
          pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
        }
      }

//      // Update camera intrinsics with refined data
//      if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
//      {
//        for (auto & intrinsic_it : sfm_data.intrinsics)
//        {
//          const IndexT indexCam = intrinsic_it.first;

//          const std::vector<double> & vec_params = map_intrinsics[indexCam];
//          intrinsic_it.second->updateFromParams(vec_params);
//        }
//      }

      return true;
    }

}

///// ***
///// The calcuation process of solving Z
///// ***
//void run_process_vertical(const SfM_Data &sfm_data_input, const std::shared_ptr<Matches_Provider> &matches_input_provider, const std::shared_ptr<Features_Provider> &features_input_provider,
//                            const std::string sOutDir, SfM_Data &sfm_data_output)
//{


////    Landmarks structure_whole;
////    get_tracks(sfm_data_input, matches_input_provider, features_input_provider, structure_whole);

//    std::vector<bool> used_X_id_list(structure_whole.size(), false);
//    int X_threshold = structure_whole.size()/4*3;

//    SfM_Data sfm_data_prepared;
//    int all_X_number = 0;
//    while(all_X_number > X_threshold)//new_C_num != 0 && used_camera_id_list.size()  )//!=  ) // camera_count > ... && number_X > ...
//    {

//        //
//        Views selected_views;
//        Matches_Provider selected_matches;
//        random_select(sfm_data_input, matches_input_provider, selected_views, selected_matches);

//        // get init subscene
//        SfM_Data sfm_data_init;

//        if(!get_init_subscene(sfm_data_input, matches_input_provider, features_input_provider, sOutDir, sfm_data_init))
//        {
//            continue;
//        }

//        //
//        int new_camera_number = 0;

//        //
//        do{
//            // X
//            get_newX(structure_whole, sfm_data_init, used_X_id_list);
//            // C
//            new_camera_number = get_newC(structure_whole, sfm_data_init);
//            //

//        }while(new_camera_number!=0);

//        //
//        all_X_number = 0;
//        for(const auto &flagX : used_X_id_list)
//        {
//            if(flagX == true)
//            {
//                ++all_X_number;
//            }

//        }

//        sfm_data_prepared = sfm_data_init;

//    }



//    /// optimisic
//    ///

//    // clean structure
//    Landmarks::iterator thisX = sfm_data_prepared.structure.begin();
//    while(thisX != sfm_data_prepared.structure.end())
//    {
//        if(used_X_id_list[thisX->first] == false)
//        {
//            thisX = sfm_data_prepared.structure.erase(thisX);
//        }else{
//            ++thisX;
//        }

//    }

//    // optimisic
//    optimisic_operation(sfm_data_prepared, sfm_data_output);

//    // save horizontal results
//    Save(
//        sfm_data_output,
//        stlplus::create_filespec( sfm_data_output.s_root_path, "sfm_data_horizontal.json" ).c_str(),
//        ESfM_Data(ALL));
//    Save(sfm_data_output,
//      stlplus::create_filespec(sOutDir, "cloud_and_poses_horizontal", ".ply"),
//      ESfM_Data(ALL));


//}

///// ***
///// select  cameras randomly and  points randomly
///// ***
//void random_select_for_Z(const SfM_Data &sfm_data, const std::shared_ptr<Matches_Provider> &matches_provider,
//                   Views &selected_views, Matches_Provider &selected_matches)
//{
//    int camera_count = sfm_data.views.size();

//    std::srand(time(NULL));
//    int cId1 = rand() % (camera_count);

//    std::srand(time(NULL));
//    int cId0 = rand() % (cId1);

//    std::srand(time(NULL));
//    int cId2 = (rand() % (camera_count-1-cId1) ) + cId1 + 1 ;

//    //
////    selected_views[cId0] = sfm_data.views[cId0].second;
////    selected_views[cId1] = sfm_data.views[cId1];
////    selected_views[cId2] = sfm_data.views[cId2];

//    //
//    std::vector<int> setected_cId_list;
//    setected_cId_list.push_back(cId0);
//    setected_cId_list.push_back(cId1);
//    setected_cId_list.push_back(cId2);

//    //
//    for (const auto & iterMatches : matches_provider->pairWise_matches_)
//    {
//      const Pair pair = iterMatches.first;
//      const int c1 = pair.first;
//      const int c2 = pair.second;

//      std::vector<int>::iterator find_c1 = std::find(setected_cId_list.begin(), setected_cId_list.end(), c1);
//      if(find_c1 != setected_cId_list.end())
//      {
//          std::vector<int>::iterator find_c2 = std::find(setected_cId_list.begin(), setected_cId_list.end(), c2);
//          if(find_c2 != setected_cId_list.end())
//          {
//              selected_matches.pairWise_matches_.insert(iterMatches);
//          }
//      }

//    }

//}
