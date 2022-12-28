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

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/cameras/Camera_undistort_image.hpp"

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


void load_matches(const std::string file_path, Landmarks &inputLM);


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
//  {
//    //--
//    // Apply external parameters (Pose)
//    //--
////#pragma omp critical
////      {
//    const T * cam_R = cam_extrinsics;
//    const T * cam_t = &cam_extrinsics[3];


//    T pos_proj[3];
//    T X[2];
//    X[0] = pos_3dpoint[0];
//    X[1] = pos_3dpoint[1];
//    //X[2] = T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);//pos_3dpoint[2];//T(0.0);
//    // Rotate the point according the camera rotation
////    ceres::AngleAxisRotatePoint(cam_R, X, pos_proj);

////    //prepare R
//    const T theta = sqrt(cam_R[0]*cam_R[0] + cam_R[1]*cam_R[1] + cam_R[2]*cam_R[2]);
//    const T cosTheta = cos(theta);
//    const T sinTheta = sin(theta);
//    const T r[3] = {cam_R[0]/theta, cam_R[1]/theta, cam_R[2]/theta};
//    const T one_cosTheta = T(1.0) - cosTheta;

//    T R[9];
//    R[0] = cosTheta + one_cosTheta*r[0]*r[0];
//    R[1] = one_cosTheta*r[0]*r[1] + sinTheta*r[2];
//    R[2] = one_cosTheta*r[0]*r[2] - sinTheta*r[1];
//    R[3] = one_cosTheta*r[0]*r[1] - sinTheta*r[2];
//    R[4] = cosTheta + one_cosTheta*r[1]*r[1];
//    R[5] = one_cosTheta*r[1]*r[2] + sinTheta*r[0];
//    R[6] = one_cosTheta*r[0]*r[2] + sinTheta*r[1];
//    R[7] = one_cosTheta*r[1]*r[2] - sinTheta*r[0];
//    R[8] = cosTheta + one_cosTheta*r[2]*r[2];

//    const T& f = cam_intrinsics[OFFSET_FOCAL_LENGTH]; //focal
//    const T& u = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X]; //principal_point_x
//    const T& v = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y]; //principal_point_y
//    //
//    const T P1 = f*R[0] + u*R[2];
//    const T P2 = f*R[3] + u*R[5];
//    const T P3 = f*R[6] + u*R[8];
//    const T P4 = f*cam_t[0] + u*cam_t[2];
//    const T P5 = f*R[1] + v*R[2];
//    const T P6 = f*R[4] + v*R[5];
//    const T P7 = f*R[7] + v*R[8];
//    const T P8 = f*cam_t[1] + v*cam_t[2];
//    const T P9 = R[2];
//    const T P10 = R[5];
//    const T P11 = R[8];
//    const T P12 = cam_t[2];
//    //
//    const T x_i = (P1*X[0] + P2*X[1] + P3*X[2] + P4) / (P9*X[0] + P10*X[1] + P11*X[2] + P12);
//    const T y_i = (P5*X[0] + P6*X[1] + P7*X[2] + P8) / (P9*X[0] + P10*X[1] + P11*X[2] + P12);
//    //
//    const T xu = T(m_pos_2dpoint[0]);
//    const T yu = T(m_pos_2dpoint[1]);
//    //
//    const T nc_x = P3/P11;
//    const T nc_y = P7/P11;
//    //
////    const T l_ncxi_x = x_i - nc_x;
////    const T l_ncxi_y = y_i - nc_y;
////    const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
////    const T normal_ncxi = (l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);
//    //
//    const T l_ncxu_x = xu - nc_x;
//    const T l_ncxu_y = yu - nc_y;
//    const T normal_ncxu = sqrt(l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);
////    const T normal_ncxu = (l_ncxu_x*l_ncxu_x + l_ncxu_y*l_ncxu_y);

//    const T A = yu - nc_y;
//    const T B = nc_x - xu;
//    const T C = xu*nc_y - yu*nc_x;

//    const T D = ceres::abs(A*x_i + B*y_i + C) / sqrt(A*A+B*B);
//    ///
//    {
//        out_residuals[0] = D * normal_ncxu / T(2.0);
//    }

//    {
////        const T l_ncxi_x = x_i - nc_x;
////        const T l_ncxi_y = y_i - nc_y;
////        const T normal_ncxi = sqrt(l_ncxi_x*l_ncxi_x + l_ncxi_y*l_ncxi_y);

////        out_residuals[0] = abs(l_ncxi_x/normal_ncxi * l_ncxu_y/normal_ncxu - l_ncxi_y/normal_ncxi * l_ncxu_x/normal_ncxu) / T(2.0);

//    }
//    return true;
//  }
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

struct ResidualErrorFunctor_Pinhole_Intrinsic_WaterZ
{
  ResidualErrorFunctor_Pinhole_Intrinsic_WaterZ(const double* const pos_2dpoint, const double* const _water_plane)
  :m_pos_2dpoint(pos_2dpoint), water_plane(_water_plane)
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

    //
    //
    const T n = (T) 1.33;
    const T * Za = (T*) water_plane;

    T nl1[3];
    nl1[0] = R[0]*x1_c[0] + R[1]*x1_c[1] + R[2];
    nl1[1] = R[3]*x1_c[0] + R[4]*x1_c[1] + R[5];
    nl1[2] = R[6]*x1_c[0] + R[7]*x1_c[1] + R[8];

    T tan_a1 = ceres::tan(ceres::acos( nl1[2] / ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1] + nl1[2]*nl1[2])));
    T a1 = ceres::acos((nl1[2]) / ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1] + nl1[2]*nl1[2]));
    if (a1 > ceres::acos(-1)/2)
    {
        a1 = ceres::acos(-1) - a1;
        tan_a1 = ceres::tan(a1);
    }
    const T tan_a2 = ceres::tan(ceres::asin(ceres::sin(a1)/n));

    const T Zp = pos_3dpoint[2];
    const T r_in = tan_a1 * (Za[0]-cam_C[2]);

    const T r_xy = ceres::sqrt((pos_3dpoint[0]-cam_C[0])*(pos_3dpoint[0]-cam_C[0]) + (pos_3dpoint[1]-cam_C[1])*(pos_3dpoint[1]-cam_C[1]) );
    const T Xq[3] = {(Za[0]-cam_C[2])/nl1[2]*nl1[0]+cam_C[0], (Za[0]-cam_C[2])/nl1[2]*nl1[1]+cam_C[1], Za[0]};
    const T l_ref[3] = {nl1[0]/ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1]), nl1[1]/ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1]), (T)1.0/tan_a2*(ceres::sin(a1)/n)*nl1[2]/ceres::abs(nl1[2])};

//    if (Zp <= Za[0])
    if (r_in >= r_xy)
    {
        const T X_subtract_C[3] = {pos_3dpoint[0]-cam_C[0], pos_3dpoint[1]-cam_C[1], pos_3dpoint[2]-cam_C[2]};
        const T crossXC[3] = {X_subtract_C[1]*nl1[2]-X_subtract_C[2]*nl1[1], X_subtract_C[2]*nl1[0]-X_subtract_C[0]*nl1[2], X_subtract_C[0]*nl1[1]-X_subtract_C[1]*nl1[0]};
        const T D_in = ceres::sqrt(crossXC[0]*crossXC[0] + crossXC[1]*crossXC[1] + crossXC[2]*crossXC[2]) / ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1] + nl1[2]*nl1[2]);

        //out_residuals[1] = r_xy-(r_in-tan_a1*(Za[0]-Zp));
        out_residuals[1] = r_in - r_xy-tan_a1*(Za[0]-Zp);

    }else{
        const T n_ref[3] = {nl1[0]/ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1]), nl1[1]/ceres::sqrt(nl1[0]*nl1[0] + nl1[1]*nl1[1]), (T)1.0/tan_a2*(ceres::sin(a1)/n)*nl1[2]/ceres::abs(nl1[2])};
        const T X_subtract_C[3] = {pos_3dpoint[0]-Xq[0], pos_3dpoint[1]-Xq[1], pos_3dpoint[2]-Xq[2]};
        const T crossXC[3] = {X_subtract_C[1]*n_ref[2]-X_subtract_C[2]*n_ref[1], X_subtract_C[2]*n_ref[0]-X_subtract_C[0]*n_ref[2], X_subtract_C[0]*n_ref[1]-X_subtract_C[1]*n_ref[0]};
        const T D_ref = ceres::sqrt(crossXC[0]*crossXC[0] + crossXC[1]*crossXC[1] + crossXC[2]*crossXC[2]) / ceres::sqrt(n_ref[0]*n_ref[0] + n_ref[1]*n_ref[1] + n_ref[2]*n_ref[2]);
        const T r_ref = tan_a2 * (Zp-Za[0]);

        out_residuals[1] = r_xy-r_in-r_ref;

    }

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double & water_plane,
    const double weight = 0.0
  )
  {
//    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_WaterZ, 2, 3, 6, 3>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_WaterZ(observation.data(), &water_plane)));
    }
//    else
//    {
//      return
//        (new ceres::AutoDiffCostFunction
//          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_WaterZ>, 2, 3, 6, 3>
//          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_WaterZ>
//            (new ResidualErrorFunctor_Pinhole_Intrinsic_WaterZ(observation.data(), &water_plane), weight)));
//    }
  }

  const double * m_pos_2dpoint; // The 2D observation
  const double * water_plane;
};


void random_select(const SfM_Data &sfm_data, const std::vector<bool> &available_poses, std::vector<bool> &record_solved_poses);
void get_newX(const std::vector<bool> &record_solved_pose, SfM_Data &sfm_data_init, std::vector<bool> &used_X_id_list);
bool calculate_X(const std::vector<bool> &record_solved_pose, SfM_Data &sfm_data_input, Landmarks &new_structure_forXxy);
bool calculate_XC(const std::vector<bool> &record_solved_pose, SfM_Data &sfm_data_input, Landmarks &new_structure_forXxy);
int get_newC(const std::vector<bool> &record_solved_Xxy, std::vector<bool> &record_solved_pose, SfM_Data &sfm_data_init);
void calculate_C(const std::vector<std::pair<Vec3, Vec2> > &X_x_list, Pose3 &camera_res );
bool calculate_Z(const SfM_Data &sfm_data_input, const std::string sOutDir,
                 double &water_plane, double &water_n, SfM_Data &sfm_data_output);

void run_process_horizontal(const SfM_Data &sfm_data_input, const std::shared_ptr<Matches_Provider> &matches_input_provider, const std::shared_ptr<Features_Provider> &features_input_provider,
                            const std::string sOutDir, SfM_Data &sfm_data_output);
bool optimisic_operation(const SfM_Data &sfm_data_input, SfM_Data &sfm_data_output);
bool get_init_guess(const SfM_Data &sfm_data, const std::shared_ptr<Matches_Provider> &matches_provider, const std::shared_ptr<Features_Provider> &features_provider,
                    const std::string sOutDir, SfM_Data &sfm_data_res);

// Z
double random_select_waterPlane(const SfM_Data &sfm_data);
bool run_process_vertical_test(const SfM_Data &sfm_data_input, const std::string sOutDir, SfM_Data &sfm_data_output);
void run_process_vertical(const SfM_Data &sfm_data_input, const std::string sOutDir, SfM_Data &sfm_data_output);
void remove2Z0(const double &water_plane, SfM_Data &sfm_data);
//
bool Adjust_water_z_fixC1(const double & water_plane, SfM_Data & sfm_data);

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

    std::ofstream out_match;
    out_match.open("/home/guang/jiao/water_data/run_pic_9_undis_/feat.txt");


    int t = 0;
    for(auto thisM : matches_provider->pairWise_matches_)
    {

        if(t==604)
        {
            Pair imgPair = thisM.first;
            std::cout << "image1: " << imgPair.first << " image2: " << imgPair.second << std::endl;
            std::vector<matching::IndMatch> features = thisM.second;
            for (auto thisFP : features)
            {
                out_match << feats_provider->feats_per_view[imgPair.first][thisFP.i_].x() << " "
                          << feats_provider->feats_per_view[imgPair.first][thisFP.i_].y() << " "
                          << feats_provider->feats_per_view[imgPair.second][thisFP.j_].x() << " "
                          << feats_provider->feats_per_view[imgPair.second][thisFP.j_].y() << std::endl;
            }
            break;

        }

        ++t;



    }

    return 0;



    Landmarks landmarks = sfm_data.structure;
    //load_matches(matches_file_path, landmarks);


    SfM_Data sfm_data_after_horizontal;
    run_process_horizontal(sfm_data, matches_provider, feats_provider, sOutDir, sfm_data_after_horizontal);
    std::cout << "horizontal process finished!" << std::endl;


    SfM_Data sfm_data_after_vertical;
//    run_process_vertical_test(sfm_data_after_horizontal, sOutDir, sfm_data_after_vertical);
    run_process_vertical(sfm_data_after_horizontal, sOutDir, sfm_data_after_vertical);
    std::cout << "vertical process finished!" << std::endl;


    Poses poses;
    int camera_num;


    ///
    //-- Export to disk computed scene (data & visualizable results)
    std::cout << "...Export SfM_Data to disk." << std::endl;
    Save(sfm_data_after_vertical,
      stlplus::create_filespec(sOutDir, "sfm_data", ".bin"),
      ESfM_Data(ALL));
}

/// ***
/// select 3 cameras randomly
/// ***
void random_select(const SfM_Data &sfm_data, const std::vector<bool> &available_poses,
                   std::vector<bool> &record_solved_poses)
{
    int camera_count = sfm_data.poses.size();

    //    std::srand(time(NULL));
    //    int cId0 = rand() % (cId1);
    int cId0 = 0;//rand() % (cId1);

    std::srand(time(NULL));
    int cId1 = rand() % (camera_count);

//    std::srand(time(NULL));
    int cId2 = (rand() % (camera_count-1-cId1) ) + cId1 + 1 ;

    //
//    selected_views[cId0] = sfm_data.views[cId0].second;
//    selected_views[cId1] = sfm_data.views[cId1];
//    selected_views[cId2] = sfm_data.views[cId2];

    //
//    std::vector<int> setected_cId_list;
//    setected_cId_list.push_back(cId0);
//    setected_cId_list.push_back(cId1);
//    setected_cId_list.push_back(cId2);

    record_solved_poses[0] = true;

    Poses::const_iterator itPose = sfm_data.poses.begin();
    std::advance(itPose, cId1);
    record_solved_poses[itPose->first] = true;

    itPose = sfm_data.poses.begin();
    std::advance(itPose, cId2);
    record_solved_poses[itPose->first] = true;

//    record_solved_poses[cId0] = true;
//    record_solved_poses[cId1] = true;
//    record_solved_poses[cId2] = true;

    //
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

}

/// ***
/// calculate sub-scene
/// ***
bool calculate_subscene(const SfM_Data &sfm_data, const std::shared_ptr<Matches_Provider> &matches_provider, const std::shared_ptr<Features_Provider> &features_provider,
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
/// Get initial guess
/// ***
bool get_init_guess(const SfM_Data &sfm_data, const std::shared_ptr<Matches_Provider> &matches_provider, const std::shared_ptr<Features_Provider> &features_provider,
                    const std::string sOutDir, SfM_Data &sfm_data_res)
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
        sfm_data_res = sfmEngine.Get_SfM_Data();
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


//    Landmarks structure_whole;
//    get_tracks(sfm_data_input, matches_input_provider, features_input_provider, structure_whole);

    // get initial gass
    SfM_Data sfm_data_init;
    get_init_guess(sfm_data_input, matches_input_provider, features_input_provider, sOutDir, sfm_data_init);

    {
////        // clean structure
//        Landmarks::iterator thisX = sfm_data_init.structure.begin();
//        while(thisX != sfm_data_init.structure.end())
//        {
//            if((thisX)->second.X(2) > 1.04)
//            {
//                thisX = sfm_data_init.structure.erase(thisX);
//            }
////            else if((thisX)->second.X(1) < 0.4){
////                thisX = sfm_data_init.structure.erase(thisX);
////            }
//                else{
//                ++thisX;
//            }

//        }
//        std::cout << "clean structure finished" << std::endl;
    }
    {
//        // move plane
//        Mat3 transR;


//        transR << 1.0, 0.0, 0.0,
//                  0.0, 0.994692564011, -0.102891989052,
//                  0.0, 0.102891989052, 0.994692564011;


//        for (auto & pose_it : sfm_data_init.poses)
//        {
//          Pose3 & pose = pose_it.second;
//          Mat3 R = pose.rotation();
//          Vec3 C = pose.center();

//          Mat3 newR = R*transR.transpose();
//          Vec3 newC = transR*C;

//          pose = Pose3(newR, newC);

//        }
//        //
//        for (auto & structure_landmark_it : sfm_data_init.structure)
//        {
//            structure_landmark_it.second.X = transR * structure_landmark_it.second.X;
//        }
//        Save(sfm_data_init,
//          stlplus::create_filespec(sOutDir, "structure_water_init_plane", "ply"),
//          ESfM_Data(EXTRINSICS | STRUCTURE));
    }


    std::vector<bool> available_pose(sfm_data_input.views.size(), false);
    for (const auto &pose : sfm_data_init.poses)
    {
        available_pose[pose.first] = true;
    }

    std::vector<bool> used_X_id_list(sfm_data_input.structure.size(), false);
    int X_threshold = sfm_data_init.structure.size()/4*3;

    SfM_Data sfm_data_prepared;
    Landmarks::const_iterator X_end = sfm_data_init.structure.end();
    --X_end;
    const int X_count = (X_end->first);
    std::vector<bool> record_solved_X(X_count, false);
    int all_X_number = 0;
    while(all_X_number < X_threshold)//new_C_num != 0 && used_camera_id_list.size()  )//!=  ) // camera_count > ... && number_X > ...
    {

        sfm_data_prepared = sfm_data_init;

        std::vector<bool> record_solved_pose(sfm_data_init.views.size(), false);
        record_solved_X = std::vector<bool>(X_count, false);


        //
//        Views selected_views;
//        Matches_Provider selected_matches;
        random_select(sfm_data_prepared, available_pose, record_solved_pose);

//        std::shared_ptr<Matches_Provider> selected_matches_ptr(new Matches_Provider(selected_matches));// = new Matches_Provider(selected_matches);


        // get init subscene


//        if(!calculate_subscene(sfm_data_input, selected_matches_ptr, features_input_provider, sOutDir, sfm_data_init))
//        {
//            continue;
//        }

        //
        int new_camera_number = 0;


        //
        do{
            // X
            get_newX(record_solved_pose, sfm_data_prepared, record_solved_X);

            // C
            new_camera_number = get_newC(record_solved_X, record_solved_pose, sfm_data_prepared);

            //judge


        }while(new_camera_number!=0);

        get_newX(record_solved_pose, sfm_data_prepared, record_solved_X);

        //
        all_X_number = 0;
        for(const auto flagX : record_solved_X)
        {
            if(flagX == true)
            {
                ++all_X_number;
            }

        }

//        sfm_data_prepared = sfm_data_init;
//        sfm_data_prepared.s_root_path = sfm_data_init.s_root_path;
//        sfm_data_prepared.poses = sfm_data_init.poses;
//        sfm_data_prepared.structure = sfm_data_init.structure;



    }

    Save(sfm_data_prepared,
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
    std::cout << "clean structure finished" << std::endl;

    // optimisic
    optimisic_operation(sfm_data_prepared, sfm_data_output);

    std::cout << "horizontal optimisic finished" << std::endl;

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

void get_newX(const std::vector<bool> &record_solved_pose, SfM_Data &sfm_data_init, std::vector<bool> &used_X_id_list)
{
//    const Poses &poses = sfm_data_init.poses;

//    std::vector<bool> record_solved_pose(sfm_data_init.views.size(), false);
//    for(const auto & pose_it : poses)
//    {
//        const IndexT indexPose = pose_it.first;
//        record_solved_pose[indexPose] = true;
//    }

    Landmarks new_structure_forXxy; // = sfm_data_init.structure;


    // find and calculate new Xxy
    for(const auto & itLM : sfm_data_init.structure)
    {
        const Observations obs_list = itLM.second.obs;

//        int num_solved_c = 0;
//        std::vector<std::pair<Pose3, Vec2> > pose_x_used;

        Observations new_obs_forXxy;
        for (const auto & itObs : obs_list)
        {
            if(record_solved_pose[itObs.first] == true)
            {
                new_obs_forXxy[itObs.first] = itObs.second;
//                pose_x_used.push_back(std::pair<Pose3, Vec2>(poses.at(itObs.first), itObs.second.x));
//                ++num_solved_c;

            }
        }

        if (new_obs_forXxy.size() > 2) //2
        {
            new_structure_forXxy[itLM.first].obs = new_obs_forXxy;
            new_structure_forXxy[itLM.first].X = sfm_data_init.structure[itLM.first].X;

            used_X_id_list[itLM.first] = true;

        }

    }

    if (new_structure_forXxy.size() == 0 )
        return;

    calculate_XC(record_solved_pose, sfm_data_init, new_structure_forXxy);

    // update
//    sfm_data_init.structure = new_structure_forXxy;

}

//bool calculate_X(const std::vector<bool> &record_solved_pose, SfM_Data &sfm_data_input, Landmarks &new_structure_forXxy) // X^T H^T [n_infy]_\times x = 0
//{
////    sfm_data_output = sfm_data_input;

//     //----------
//     // Add camera parameters
//     // - intrinsics
//     // - poses [R|t]

//     // Create residuals for each observation in the bundle adjustment problem. The
//     // parameters for cameras and points are added automatically.
//     //----------

//     ceres::Problem problem;

//     // Data wrapper for refinement:
//     Hash_Map<IndexT, std::vector<double> > map_intrinsics;
//     Hash_Map<IndexT, std::vector<double> > map_poses;

//     // Setup Poses data & subparametrization
//     for (const auto & pose_it : sfm_data_input.poses)
//     {
//       const IndexT indexPose = pose_it.first;

//       if(record_solved_pose[indexPose] == false)
//       {
//           continue;
//       }

//       const Pose3 & pose = pose_it.second;
//       const Mat3 R = pose.rotation();
//       const Vec3 t = pose.translation();
//       const Vec3 C = pose.center();

//       double angleAxis[3];
//       ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
//       // angleAxis + translation
////       map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};
//       map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(0), C(1), C(2)};

//       double * parameter_block = &map_poses[indexPose][0];
//       problem.AddParameterBlock(parameter_block, 6);

//       // set the whole parameter block as constant for best performance
//       problem.SetParameterBlockConstant(parameter_block);

//     }

//     // Setup Intrinsics data & subparametrization
//     for (const auto & intrinsic_it : sfm_data_input.intrinsics)
//     {
//       const IndexT indexCam = intrinsic_it.first;

//       if (isValid(intrinsic_it.second->getType()))
//       {
//         map_intrinsics[indexCam] = intrinsic_it.second->getParams();
//         if (!map_intrinsics[indexCam].empty())
//         {
//           double * parameter_block = &map_intrinsics[indexCam][0];
//           problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());

//           // set the whole parameter block as constant for best performance
//           problem.SetParameterBlockConstant(parameter_block);

//         }
//       }
//       else
//       {
//         std::cerr << "Unsupported camera type." << std::endl;
//       }
//     }

//     // Set a LossFunction to be less penalized by false measurements
//     //  - set it to nullptr if you don't want use a lossFunction.
//     ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(Square(4.0)); //nullptr;

//     // For all visibility add reprojections errors:
//     for (auto & structure_landmark_it : new_structure_forXxy)
//     {
//       const Observations & obs = structure_landmark_it.second.obs;

//       for (const auto & obs_it : obs)
//       {
//         // Build the residual block corresponding to the track observation:
//         const View * view = sfm_data_input.views.at(obs_it.first).get();

//         // Each Residual block takes a point and a camera as input and outputs a 2
//         // dimensional residual. Internally, the cost function stores the observed
//         // image location and compares the reprojection against the observation.
//         // !!!
//         ceres::CostFunction* cost_function =
//                 ResidualErrorFunctor_Pinhole_Intrinsic_WaterXY::Create(obs_it.second.x);


//         if (cost_function)
//         {
//           if (!map_intrinsics[view->id_intrinsic].empty())
//           {
//             problem.AddResidualBlock(cost_function,
//               p_LossFunction,
//               &map_intrinsics[view->id_intrinsic][0],
//               &map_poses[view->id_pose][0],
//               structure_landmark_it.second.X.data());
//           }
//           else
//           {
//             problem.AddResidualBlock(cost_function,
//               p_LossFunction,
//               &map_poses[view->id_pose][0],
//               structure_landmark_it.second.X.data());
//           }
//         }
//       }
//       problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());

//     }

//     ceres::Solver::Options ceres_config_options;
//     ceres_config_options.max_num_iterations = 500;
//     ceres_config_options.preconditioner_type = ceres::JACOBI;
//       //static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
//     ceres_config_options.linear_solver_type =  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
//       //static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
//     ceres_config_options.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
//       //static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
//     ceres_config_options.minimizer_progress_to_stdout = 500;//ceres_options_.bVerbose_;
//     ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
//     ceres_config_options.num_threads = 6;//ceres_options_.nb_threads_;
//     ceres_config_options.num_linear_solver_threads = 6;//ceres_options_.nb_threads_;
//     ceres_config_options.parameter_tolerance = 0.000000001;//ceres_options_.parameter_tolerance_;

//   //  std::cout << "start ceres solver" << std::endl;
//     // Solve BA
//     ceres::Solver::Summary summary;
//     ceres::Solve(ceres_config_options, &problem, &summary);
// //    if (ceres_options_.bCeres_summary_)
// //      std::cout << summary.FullReport() << std::endl;

//     // If no error, get back refined parameters
//     if (!summary.IsSolutionUsable())
//     {
// //      if (ceres_options_.bVerbose_)
//       std::cout << "Bundle Adjustment failed." << std::endl;
//       return false;
//     }
//     else // Solution is usable
//     {
//       //if (ceres_options_.bVerbose_)
//       {
//         // Display statistics about the minimization
//         std::cout << std::endl
//           << "Bundle Adjustment statistics (approximated RMSE):\n"
//           << " #views: " << sfm_data_input.views.size() << "\n"
//           << " #poses: " << sfm_data_input.poses.size() << "\n"
//           << " #intrinsics: " << sfm_data_input.intrinsics.size() << "\n"
//           << " #tracks: " << new_structure_forXxy.size() << "\n"
//           << " #residuals: " << summary.num_residuals << "\n"
//           << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
//           << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
//           << " BriefReport : " << summary.BriefReport() << "\n"
//           << " Time (s): " << summary.total_time_in_seconds << "\n"
//           << std::endl;
//       }

//       // Update camera poses with refined data
//       //if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
//       {
//         for (auto & pose_it : sfm_data_input.poses)
//         {
//           const IndexT indexPose = pose_it.first;

//           Mat3 R_refined;
//           ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
////           Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
//           Vec3 C_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
//           // Update the pose
//           Pose3 & pose = pose_it.second;
////           pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
//           pose = Pose3(R_refined, C_refined);
//         }

//         for (auto & structure_landmark_it : new_structure_forXxy)
//         {
//             sfm_data_input.structure[structure_landmark_it.first].X = new_structure_forXxy[structure_landmark_it.first].X;
//         }

//       }

//       return true;
//     }


//}

int get_newC(const std::vector<bool> &record_solved_Xxy, std::vector<bool> &record_solved_pose, SfM_Data &sfm_data_init)//(const Landmarks matches_LM, const int camera_num, Poses &poses)
{
//    const Poses &poses = sfm_data_init.poses;

//    std::vector<bool> record_solved_pose(sfm_data_init.views.size(), false);
//    for(const auto & pose_it : poses)
//    {
//        const IndexT indexPose = pose_it.first;
//        record_solved_pose[indexPose] = true;
//    }

    int new_pose_number = 0;

    Landmarks new_structure_forCxy;
    std::vector<int> record_poseX_times(sfm_data_init.views.size(), 0);


    //
    for(const auto & itLM : sfm_data_init.structure)
    {
        if (record_solved_Xxy[itLM.first] == false)
            continue;

        const Observations obs_list = itLM.second.obs;

        Observations new_obs_forCxy;
        for (const auto & itObs : obs_list)
        {
            if(record_solved_pose[itObs.first] == true)
            {
                new_obs_forCxy[itObs.first] = itObs.second;
            }else{
                ++record_poseX_times[itObs.first];
            }
        }

        if (new_obs_forCxy.size() > 2)
        {
            new_structure_forCxy[itLM.first].obs = new_obs_forCxy;
            new_structure_forCxy[itLM.first].X = sfm_data_init.structure[itLM.first].X;

        }

    }

    //
    for (std::size_t pId = 0; pId < record_poseX_times.size(); ++pId)
    {
        if (record_poseX_times[pId] > 2)
        {
            record_solved_pose[pId] = true;
            ++new_pose_number;
        }

    }

    //
    for(const auto & itLM : sfm_data_init.structure)
    {
        if (record_solved_Xxy[itLM.first] == false)
            continue;

        const Observations obs_list = itLM.second.obs;

        //Observations new_obs_forCxy;
        for (const auto & itObs : obs_list)
        {

            if(record_poseX_times[itObs.first] > 2)
            {
                new_structure_forCxy[itLM.first].obs[itObs.first] = itObs.second;

            }
        }

    }

    if(new_structure_forCxy.size() == 0)
        return 0;

    calculate_XC(record_solved_pose, sfm_data_init, new_structure_forCxy);

    return new_pose_number;

}

bool calculate_XC(const std::vector<bool> &record_solved_pose, SfM_Data &sfm_data_input, Landmarks &new_structure_forXxy) // X^T H^T [n_infy]_\times x = 0
{
//    sfm_data_output = sfm_data_input;

     //----------
     // Add camera parameters
     // - intrinsics
     // - poses [R|t]

     // Create residuals for each observation in the bundle adjustment problem. The
     // parameters for cameras and points are added automatically.
     //----------

     ceres::Problem problem;

     // Data wrapper for refinement:
     Hash_Map<IndexT, std::vector<double> > map_intrinsics;
     Hash_Map<IndexT, std::vector<double> > map_poses;

     // Setup Poses data & subparametrization
     for (const auto & pose_it : sfm_data_input.poses)
     {
       const IndexT indexPose = pose_it.first;

       if(record_solved_pose[indexPose] == false)
       {
           continue;
       }

       const Pose3 & pose = pose_it.second;
       const Mat3 R = pose.rotation();
       const Vec3 C = pose.center();

       double angleAxis[3];
       ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
       // angleAxis + translation
       //map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(0), C(1), C(2)};
       map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(0), C(1)};

       double * parameter_block = &map_poses[indexPose][0];
       problem.AddParameterBlock(parameter_block, 5);

     }

//     for (const auto & pose_it : sfm_data_input.poses)
//     {
//       const IndexT indexPose = pose_it.first;

//       if(record_solved_pose[indexPose] == false)
//       {
//           continue;
//       }

//       const Pose3 & pose = pose_it.second;
//       const Mat3 R = pose.rotation();
//       const Vec3 t = pose.translation();

//       double angleAxis[3];
//       ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
//       // angleAxis + translation
//       map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

//       double * parameter_block = &map_poses[indexPose][0];
//       problem.AddParameterBlock(parameter_block, 6);

//       // set the whole parameter block as constant for best performance
//       problem.SetParameterBlockConstant(parameter_block);

//     }

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
           problem.AddParameterBlock(parameter_block, 3);

           // set the whole parameter block as constant for best performance
           problem.SetParameterBlockConstant(parameter_block);

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
     for (auto & structure_landmark_it : new_structure_forXxy)
     {
       const Observations & obs = structure_landmark_it.second.obs;

       for (const auto & obs_it : obs)
       {
         // Build the residual block corresponding to the track observation:
         const View * view = sfm_data_input.views.at(obs_it.first).get();

         if(record_solved_pose[view->id_pose] == false)
         {
             std::cout << "error here" << std::endl;
             continue;
         }
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
           << " #views: " << sfm_data_input.views.size() << "\n"
           << " #poses: " << sfm_data_input.poses.size() << "\n"
           << " #intrinsics: " << sfm_data_input.intrinsics.size() << "\n"
           << " #tracks: " << new_structure_forXxy.size() << "\n"
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
           for (auto & pose_it : sfm_data_input.poses)
           {
             const IndexT indexPose = pose_it.first;

             if(record_solved_pose[indexPose] == false)
             {
                 continue;
             }

             Mat3 R_refined;
             ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
             Pose3 & pose = pose_it.second;
             Vec3 C_init = pose.center();
//             Vec3 t_init = pose.translation();
             Vec3 C_refined(map_poses[indexPose][3], map_poses[indexPose][4], C_init(2));
             //Vec3 C_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
             // Update the pose
     //        Pose3 & pose = pose_it.second;
             pose = Pose3(R_refined, C_refined);

             //std::cout << "diff pose : " << init_map_poses[indexPose][3] - map_poses[indexPose][3] << " " << init_map_poses[indexPose][4] - map_poses[indexPose][4] << std::endl;
           }
         }

//       {
//         for (auto & pose_it : sfm_data_input.poses)
//         {
//           const IndexT indexPose = pose_it.first;

//           Mat3 R_refined;
//           ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
//           Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
//           // Update the pose
//           Pose3 & pose = pose_it.second;
//           pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
//         }

//         for (auto & structure_landmark_it : new_structure_forXxy)
//         {
//             sfm_data_input.structure[structure_landmark_it.first].X = new_structure_forXxy[structure_landmark_it.first].X;
//         }

//       }

       return true;
     }


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
      const Vec3 C = pose.center();

      double angleAxis[3];
      ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
      // angleAxis + translation
//      map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};
      map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(0), C(1)};

      double * parameter_block = &map_poses[indexPose][0];
      problem.AddParameterBlock(parameter_block, 5);

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
          Pose3 & pose = pose_it.second;
          Vec3 C_init = pose.center();
//             Vec3 t_init = pose.translation();
          Vec3 C_refined(map_poses[indexPose][3], map_poses[indexPose][4], C_init(2));
          //Vec3 C_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
          // Update the pose
  //        Pose3 & pose = pose_it.second;
          pose = Pose3(R_refined, C_refined);
        }
//        {
//          const IndexT indexPose = pose_it.first;

//          Mat3 R_refined;
//          ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
//          Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
//          // Update the pose
//          Pose3 & pose = pose_it.second;
//          pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
//        }


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

bool judge_stopZ(const SfM_Data &sfm_data, const double &water_n, const double &water_plane)
{
    if(water_n < 1.3 || water_n > 1.4)
    {
        return false;
    }

    int underwater_point_count = 0;
    for (const auto thisX : sfm_data.structure)
    {
        if(thisX.second.X(2) > water_plane)
        {
            ++underwater_point_count;
        }
    }

    if(underwater_point_count<sfm_data.structure.size()/10)
    {
        return false;
    }

    return true;

}

/// ***
/// The calcuation process of solving Z
/// ***
void run_process_vertical(const SfM_Data &sfm_data_input, const std::string sOutDir, SfM_Data &sfm_data_output)
{

    double water_n, water_plane;
    bool cal_flag, judge_flag;

    SfM_Data sfm_data_cal;

    do{
        water_n = 1.33;
        water_plane = random_select_waterPlane(sfm_data_input);

        cal_flag = calculate_Z(sfm_data_input, sOutDir, water_plane, water_n, sfm_data_cal);

        judge_flag = judge_stopZ(sfm_data_cal, water_n, water_plane);

        Save(sfm_data_cal,
             stlplus::create_filespec(sOutDir, "cloud_and_poses_vertical_temp", "ply"),
             ESfM_Data(ALL));

    }while(!cal_flag || !judge_flag);

    sfm_data_output = sfm_data_cal;

    std::cout << "vertical optimisic finished" << std::endl;

    // remove water plane to z=0
    remove2Z0(water_plane, sfm_data_output);

    Save(sfm_data_output,
      stlplus::create_filespec(sOutDir, "cloud_and_poses_vertical", ".ply"),
      ESfM_Data(ALL));

    std::stringstream wp_ss;
    wp_ss << water_plane;
    std::string wp_s;
    wp_ss >> wp_s;
    std::ofstream record_wp;
    record_wp.open(sOutDir + "/wp\ =\ " + wp_s);
    record_wp.close();

    std::stringstream n_ss;
    n_ss << water_n;
    std::string n_s;
    n_ss >> n_s;
    std::ofstream record_n;
    record_n.open(sOutDir + "/n\ =\ " + n_s);
    record_n.close();

}

void remove2Z0(const double &water_plane, SfM_Data &sfm_data)
{
    // transform Z-axis
    Mat3 transR;
    transR << 0.0, 1.0, 0.0,
              1.0, 0.0, 0.0,
              0.0, 0.0, -1.0;
//      transR << 1.0, 0.0, 0.0,
//                0.0, 1.0, 0.0,
//                0.0, 0.0, 1.0;

    for (auto & pose_it : sfm_data.poses)
    {
      Pose3 & pose = pose_it.second;
      Mat3 R = pose.rotation();
      Vec3 C = pose.center();
      C(2) = C(2) - water_plane;

      Mat3 newR = R*transR.transpose();
      Vec3 newC = transR*C;
      //newC(2) = newC(2) + water_plane;

      pose = Pose3(newR, newC);

    }
    //
    for (auto & structure_landmark_it : sfm_data.structure)
    {
        structure_landmark_it.second.X(2) = structure_landmark_it.second.X(2) - water_plane;
        structure_landmark_it.second.X = transR * structure_landmark_it.second.X;
        //structure_landmark_it.second.X(2) = structure_landmark_it.second.X(2) + water_plane;
    }

}

bool run_process_vertical_test(const SfM_Data &sfm_data_input, const std::string sOutDir, SfM_Data &sfm_data_output)
{
    sfm_data_output = sfm_data_input;

    std::cout << "with pose id: " << std::endl;
    for (const auto & pose_it : sfm_data_output.poses)
    {
      const IndexT indexPose = pose_it.first;
      std::cout << indexPose << " ";
    }
    std::cout << std::endl;

    {
//        // if need trans2plane
//        Mat3 transR;

//        transR << 0.0, 1.0, 0.0,
//                  1.0, 0.0, 0.0,
//                  0.0, 0.0, -1.0; // synthetic

//        std::cout << "transR : " << transR << std::endl;

//        for (auto & pose_it : sfm_data_output.poses)
//        {
//          Pose3 & pose = pose_it.second;
//          Mat3 R = pose.rotation();
//          Vec3 C = pose.center();

//          Mat3 newR = R*transR.transpose();
//          Vec3 newC = transR*C;

//          pose = Pose3(newR, newC);

//        }
//        //
//        for (auto & structure_landmark_it : sfm_data_output.structure)
//        {
//            structure_landmark_it.second.X = transR * structure_landmark_it.second.X;
//        }
//        Save(sfm_data_output,
//          stlplus::create_filespec(sOutDir, "structure_water_init_plane", "ply"),
//          ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    Bundle_Adjustment_Ceres bundle_adjustment_obj;
    bool b_BA_Status;


    //double water_plane = 1.22;//
    //double water_plane = 1.04;//run_6
//    double water_plane = 0.99;//run_6
    double water_plane = 1.223;
//    double water_plane = 1.006;//run_6
    //    double water_plane = -0.4;//run_6

//        b_BA_Status = bundle_adjustment_obj.Adjust_water_min_z
//          (
//            sfm_data_output,
//            water_plane,
//            Optimize_Options(
//              openMVG::cameras::Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//              Extrinsic_Parameter_Type::NONE,//ADJUST_TRANSLATION, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//              Structure_Parameter_Type::ADJUST_ALL,
//              Control_Point_Parameter(),
//              false)
//          );

//        std::cout << "water_plane : " << water_plane << std::endl;

      b_BA_Status = Adjust_water_z_fixC1(water_plane, sfm_data_output);

      Save(sfm_data_output,
        stlplus::create_filespec(sOutDir, "structure_water_refine_Z", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));

      //
//      b_BA_Status = bundle_adjustment_obj.Adjust_water_min_z
//        (
//          sfm_data_output,
//          water_plane,
//          Optimize_Options(
//            openMVG::cameras::Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//            Extrinsic_Parameter_Type::NONE,//ADJUST_TRANSLATION, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//            Structure_Parameter_Type::ADJUST_ALL,
//            Control_Point_Parameter(),
//            false)
//        );
      Save(sfm_data_output,
        stlplus::create_filespec(sOutDir, "structure_water_refine_Z2", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));

    std::cout << "water_plane : " << water_plane << std::endl;

//    b_BA_Status = bundle_adjustment_obj.Adjust_water_z_fixC1
//    (
//      sfm_data_output,
//      water_plane,
//      Optimize_Options(
//        openMVG::cameras::Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::NONE,//ADJUST_TRANSLATION, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        false)
//    );
    Save(sfm_data_output,
    stlplus::create_filespec(sOutDir, "structure_water_refine_Z3", "ply"),
    ESfM_Data(EXTRINSICS | STRUCTURE));

    double init_ref_n = 1.33;

//    b_BA_Status = bundle_adjustment_obj.Adjust_water_z
//    (
//      sfm_data_output,
//      water_plane,
//      Optimize_Options(
//        openMVG::cameras::Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::NONE,//ADJUST_TRANSLATION, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        false)
//    );
    Save(sfm_data_output,
    stlplus::create_filespec(sOutDir, "structure_water_refine_Z4", "ply"),
    ESfM_Data(EXTRINSICS | STRUCTURE));

    std::cout << "water_plane : " << water_plane << std::endl;


    b_BA_Status = bundle_adjustment_obj.Adjust_water_z_n
    (
      sfm_data_output,
      water_plane,
      init_ref_n,
      Optimize_Options(
        openMVG::cameras::Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::NONE,//ADJUST_TRANSLATION, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
        Structure_Parameter_Type::ADJUST_ALL,
        Control_Point_Parameter(),
        false)
    );
    Save(sfm_data_output,
    stlplus::create_filespec(sOutDir, "structure_water_refine_Zn", "ply"),
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
//      transR << 1.0, 0.0, 0.0,
//                0.0, 1.0, 0.0,
//                0.0, 0.0, 1.0;

      for (auto & pose_it : sfm_data_output.poses)
      {
        Pose3 & pose = pose_it.second;
        Mat3 R = pose.rotation();
        Vec3 C = pose.center();
        C(2) = C(2) - water_plane;

        Mat3 newR = R*transR.transpose();
        Vec3 newC = transR*C;
        //newC(2) = newC(2) + water_plane;

        pose = Pose3(newR, newC);

      }
      //
      for (auto & structure_landmark_it : sfm_data_output.structure)
      {
          structure_landmark_it.second.X(2) = structure_landmark_it.second.X(2) - water_plane;
          structure_landmark_it.second.X = transR * structure_landmark_it.second.X;
          //structure_landmark_it.second.X(2) = structure_landmark_it.second.X(2) + water_plane;
      }
      Save(sfm_data_output,
        stlplus::create_filespec(sOutDir, "structure_vertical_finished_Z", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
    return true;

}


bool calculate_Z(const SfM_Data &sfm_data_input, const std::string sOutDir,
                 double &water_plane, double &water_n, SfM_Data &sfm_data_output)
{
    sfm_data_output = sfm_data_input;

    std::cout << "with pose id: " << std::endl;
    for (const auto & pose_it : sfm_data_output.poses)
    {
      const IndexT indexPose = pose_it.first;
      std::cout << indexPose << " ";
    }
    std::cout << std::endl;

    {
//        // if need trans2plane
//        Mat3 transR;

//        transR << 0.0, 1.0, 0.0,
//                  1.0, 0.0, 0.0,
//                  0.0, 0.0, -1.0; // synthetic

//        std::cout << "transR : " << transR << std::endl;

//        for (auto & pose_it : sfm_data_output.poses)
//        {
//          Pose3 & pose = pose_it.second;
//          Mat3 R = pose.rotation();
//          Vec3 C = pose.center();

//          Mat3 newR = R*transR.transpose();
//          Vec3 newC = transR*C;

//          pose = Pose3(newR, newC);

//        }
//        //
//        for (auto & structure_landmark_it : sfm_data_output.structure)
//        {
//            structure_landmark_it.second.X = transR * structure_landmark_it.second.X;
//        }
//        Save(sfm_data_output,
//          stlplus::create_filespec(sOutDir, "structure_water_init_plane", "ply"),
//          ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    Bundle_Adjustment_Ceres bundle_adjustment_obj;
    bool b_BA_Status;


//        b_BA_Status = bundle_adjustment_obj.Adjust_water_min_z
//          (
//            sfm_data_output,
//            water_plane,
//            Optimize_Options(
//              openMVG::cameras::Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//              Extrinsic_Parameter_Type::NONE,//ADJUST_TRANSLATION, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//              Structure_Parameter_Type::ADJUST_ALL,
//              Control_Point_Parameter(),
//              false)
//          );

//        std::cout << "water_plane : " << water_plane << std::endl;

      b_BA_Status = Adjust_water_z_fixC1(water_plane, sfm_data_output);

      Save(sfm_data_output,
        stlplus::create_filespec(sOutDir, "structure_water_refine_Z", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));

      //
//      b_BA_Status = bundle_adjustment_obj.Adjust_water_min_z
//        (
//          sfm_data_output,
//          water_plane,
//          Optimize_Options(
//            openMVG::cameras::Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//            Extrinsic_Parameter_Type::NONE,//ADJUST_TRANSLATION, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//            Structure_Parameter_Type::ADJUST_ALL,
//            Control_Point_Parameter(),
//            false)
//        );
      Save(sfm_data_output,
        stlplus::create_filespec(sOutDir, "structure_water_refine_Z2", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));

    std::cout << "water_plane : " << water_plane << std::endl;

//    b_BA_Status = bundle_adjustment_obj.Adjust_water_z_fixC1
//    (
//      sfm_data_output,
//      water_plane,
//      Optimize_Options(
//        openMVG::cameras::Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::NONE,//ADJUST_TRANSLATION, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        false)
//    );
    Save(sfm_data_output,
    stlplus::create_filespec(sOutDir, "structure_water_refine_Z3", "ply"),
    ESfM_Data(EXTRINSICS | STRUCTURE));

//    double init_ref_n = 1.33;

//    b_BA_Status = bundle_adjustment_obj.Adjust_water_z
//    (
//      sfm_data_output,
//      water_plane,
//      Optimize_Options(
//        openMVG::cameras::Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
//        Extrinsic_Parameter_Type::NONE,//ADJUST_TRANSLATION, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
//        Structure_Parameter_Type::ADJUST_ALL,
//        Control_Point_Parameter(),
//        false)
//    );
    Save(sfm_data_output,
    stlplus::create_filespec(sOutDir, "structure_water_refine_Z4", "ply"),
    ESfM_Data(EXTRINSICS | STRUCTURE));

    std::cout << "water_plane : " << water_plane << std::endl;


    b_BA_Status = bundle_adjustment_obj.Adjust_water_z_n
    (
      sfm_data_output,
      water_plane,
      water_n,
      Optimize_Options(
        openMVG::cameras::Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::NONE,//ADJUST_TRANSLATION, // Rotations are held as constant//Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
        Structure_Parameter_Type::ADJUST_ALL,
        Control_Point_Parameter(),
        false)
    );
    Save(sfm_data_output,
    stlplus::create_filespec(sOutDir, "structure_water_refine_Zn", "ply"),
    ESfM_Data(EXTRINSICS | STRUCTURE));

    std::cout << "water_plane : " << water_plane << std::endl;
    std::cout << "Adjust_water_z_n n = " << water_n  << std::endl;
    std::cout << "Adjust_water_z_n finished" << std::endl;

    {
//      // transform Z-axis
//      Mat3 transR;
//      transR << 0.0, 1.0, 0.0,
//                1.0, 0.0, 0.0,
//                0.0, 0.0, -1.0;
////      transR << 1.0, 0.0, 0.0,
////                0.0, 1.0, 0.0,
////                0.0, 0.0, 1.0;

//      for (auto & pose_it : sfm_data_output.poses)
//      {
//        Pose3 & pose = pose_it.second;
//        Mat3 R = pose.rotation();
//        Vec3 C = pose.center();
//        C(2) = C(2) - water_plane;

//        Mat3 newR = R*transR.transpose();
//        Vec3 newC = transR*C;
//        //newC(2) = newC(2) + water_plane;

//        pose = Pose3(newR, newC);

//      }
//      //
//      for (auto & structure_landmark_it : sfm_data_output.structure)
//      {
//          structure_landmark_it.second.X(2) = structure_landmark_it.second.X(2) - water_plane;
//          structure_landmark_it.second.X = transR * structure_landmark_it.second.X;
//          //structure_landmark_it.second.X(2) = structure_landmark_it.second.X(2) + water_plane;
//      }
//      Save(sfm_data_output,
//        stlplus::create_filespec(sOutDir, "structure_vertical_finished_Z", "ply"),
//        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
    return true;

}

/// ***
/// select 3 cameras randomly
/// ***
double random_select_waterPlane(const SfM_Data &sfm_data)
{
    const int quantify_count = 1000;

    Landmarks::const_iterator itFirst = sfm_data.structure.begin();
    double maxXz = itFirst->second.X(2);
    double minXz = itFirst->second.X(2);//(*itFirst).second.X(2);
    for(const auto itX : sfm_data.structure)
    {
        if (itX.second.X(2) < minXz)
        {
            minXz = itX.second.X(2);
        }else if (itX.second.X(2) > maxXz)
        {
            maxXz = itX.second.X(2);
        }

    }

    const double quantify_interval = (maxXz-minXz)/(double)quantify_count;

    std::srand(time(NULL));
    const int quantify_num = rand() % (quantify_count);

    const double waterPlane = minXz + quantify_num * quantify_interval;

    return waterPlane;

}

bool Adjust_water_z_fixC1(const double & water_plane, SfM_Data & sfm_data)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses, init_map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 C = pose.center();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(0), C(1), C(2)};
    //map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], C(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);

    if (indexPose == 0 || indexPose == 1)
        problem.SetParameterBlockConstant(parameter_block);

  }
  init_map_poses = map_poses;

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();
      if (!map_intrinsics[indexCam].empty())
      {
        double * parameter_block = &map_intrinsics[indexCam][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());

        problem.SetParameterBlockConstant(parameter_block);

//        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
//        {
//          // set the whole parameter block as constant for best performance
//          problem.SetParameterBlockConstant(parameter_block);
//        }
//        else
//        {
//          const std::vector<int> vec_constant_intrinsic =
//            intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
//          if (!vec_constant_intrinsic.empty())
//          {
//            ceres::SubsetParameterization *subset_parameterization =
//              new ceres::SubsetParameterization(
//                map_intrinsics[indexCam].size(), vec_constant_intrinsic);
//            problem.SetParameterization(parameter_block, subset_parameterization);
//          }
//        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to nullptr if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(16.0);


  // For all visibility add reprojections errors:
  Landmarks save_structure;//sfm_data.structure
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

//        if (structure_landmark_it.second.X(2) > 2.53322799) //72
//            continue;
//        if (structure_landmark_it.second.X(2) > 1.0416) //73
//            continue;

//        if (structure_landmark_it.second.X(2) > 1.063) //71_1
//            continue;

    if (obs.size() < 3)
        continue;

    save_structure.insert(structure_landmark_it);
  }
  sfm_data.structure = save_structure;
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        ResidualErrorFunctor_Pinhole_Intrinsic_WaterZ::Create(obs_it.second.x, water_plane);
      problem.AddResidualBlock(cost_function,
        p_LossFunction,
        &map_intrinsics[view->id_intrinsic][0],
        &map_poses[view->id_pose][0],
        structure_landmark_it.second.X.data());


    }

  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
//  sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
//  linear_solver_type_ = ceres::SPARSE_SCHUR;
//  linear_solver_type_ = ceres::DENSE_SCHUR;
//  preconditioner_type_ = ceres::JACOBI;



  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 1500;
  ceres_config_options.preconditioner_type = ceres::JACOBI;
//    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type = ceres::DENSE_SCHUR;//  ceres::DENSE_SCHUR;//ceres::DENSE_SCHUR;//ceres::SPARSE_SCHUR;
//    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
//    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = 500;//ceres_options_.bVerbose_;
  ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;//ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
  ceres_config_options.num_threads = 4;//ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = 4;//ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = 0.000000001;//ceres_options_.parameter_tolerance_/10000;

//  std::cout << "start ceres solver" << std::endl;
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);


  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
//    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
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
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Pose3 & pose = pose_it.second;
        Vec3 C_init = pose.center();
        //Vec3 C_refined(C_init(0), C_init(1), map_poses[indexPose][3]);
        Vec3 C_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
//        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, C_refined);

        //std::cout << "diff pose : " << init_map_poses[indexPose][3] - map_poses[indexPose][3] << " " << init_map_poses[indexPose][4] - map_poses[indexPose][4] << std::endl;
      }
    }

    // Update camera intrinsics with refined data
    //if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second->updateFromParams(vec_params);
      }
    }

    // Structure is already updated directly if needed (no data wrapping)


    return true;
  }
}

/// ***
///
///
//bool saveData_for_show(const SfM_Data &sfm_data_undiswaterdata, const SfM_Data &sfm_data_diswaterdata, const std::string sOutDir)
//{
//    {
////        SfM_Data newdata = sfmEngine.Get_SfM_Data();
//        std::ofstream out_e, out_k;
//        std::string e_path = sfm_data_undiswaterdata.s_root_path+"/RC_new_withWater.txt";//undisWater.txt";//
//        std::string k_path = sfm_data_undiswaterdata.s_root_path+"/K_new_withWater.txt";//undisWater.txt";//
//        out_e.open(e_path);
//        out_k.open(k_path);
//        if(!out_e)
//        {
//            std::cout << "open file " << e_path << std::endl;
//            return EXIT_FAILURE;
//        }

//        for(int imgId = 0; imgId < sfmEngine.Get_SfM_Data().poses.size(); ++imgId)
//        {
//            const Pose3 pose = sfmEngine.Get_SfM_Data().poses.at(imgId);
//            Mat3 R = pose.rotation();
//            Vec3 C = pose.center();

//            out_e << R(0,0) << " " << R(0,1) << " " << R(0,2) << " " << C(0) << std::endl
//                  << R(1,0) << " " << R(1,1) << " " << R(1,2) << " " << C(1) << std::endl
//                  << R(2,0) << " " << R(2,1) << " " << R(2,2) << " " << C(2) << std::endl;

//            std::vector<double> cam_intrinsics = sfm_data.intrinsics[sfm_data.views[imgId]->id_intrinsic]->getParams();
//            out_k << cam_intrinsics[0] << " " << cam_intrinsics[1] << " " << cam_intrinsics[2] << std::endl;
//        }
//        out_e.close();
//        out_k.close();

//        //save for analyse error
//        std::ofstream out_X, out_xm;
//        out_X.open(sOutDir+"/X_withWater.txt");//undisWater.txt");//withWater.txt");//
//        out_xm.open(sfmEngine.Get_SfM_Data().s_root_path+"/xm_withWater.txt");//undisWater.txt");//withWater.txt");
//        for(Landmarks::const_iterator itX = undiswaterdata.structure.begin();
//            itX != undiswaterdata.structure.end(); ++itX)
//        {
//            //for pot17
////            if(itX->second.X(0)<-0.706062 || itX->second.X(0) > 0.618697 || itX->second.X(1)<-0.691933 || itX->second.X(1)>0.29423 || itX->second.X(2) > 2.551696)
////            {

////                continue;
////            }

////            if(itX->second.X(2) > 1.114)
////            {
////                continue;
////            }

//            out_xm << itX->second.obs.size() << " ";
//            for(Observations::const_iterator itObs = itX->second.obs.begin();
//                itObs != itX->second.obs.end(); ++itObs)
//            {
//                out_xm << itObs->first << " " << itObs->second.id_feat << " " << itObs->second.x(0) << " " << itObs->second.x(1) << " ";

//            }
//            out_xm << itX->second.X(0) << " "
//                   << itX->second.X(1) << " "
//                   << itX->second.X(2) << std::endl;
//            out_X << itX->second.obs.begin()->first << " "
//                  << itX->second.obs.begin()->second.id_feat << " "
//                  << itX->second.X(0) << " "
//                  << itX->second.X(1) << " "
//                  << itX->second.X(2) << std::endl;

//        }
//        out_X.close();
//        out_xm.close();
//    }

//}
