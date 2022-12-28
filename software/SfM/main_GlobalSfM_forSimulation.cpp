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

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include "openMVG/multiview/triangulation_nview.hpp"

#include <cstdlib>
#include <memory>
#include <string>

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <ceres/problem.h>
#include <ceres/solver.h>
#include <ceres/types.h>
#include <ceres/cost_function.h>
#include <list>

#include <Eigen/Dense>

using namespace openMVG;
using namespace openMVG::sfm;

#define WATER_N 1.33

//struct WaterPlaneCostFunction
//{
//  const std::list<std::pair<Eigen::Vector3d, Eigen::Vector3d> > Llist;
////  const std::vector<double> CJ_Inf;

//  WaterPlaneCostFunction
//  (
//    const std::list<std::pair<Eigen::Vector3d, Eigen::Vector3d> > & _Llist
//  ): Llist(_Llist)
//  {
//  }

//  template <typename T> bool
//  operator()
//  (
//    const T* const Za,
//    T* residuals
//  )
////  template <typename T>
////  bool
////  operator()
////  (
////    const double* const Za,
////    double* residuals
////  )
//  const
//  {

//      std::list<std::pair<Eigen::Matrix<T,4,1>, Eigen::Matrix<T,4,1> > > AllLine;
//      for(std::list<std::pair<Eigen::Vector3d, Eigen::Vector3d> >::const_iterator itL = Llist.begin();
//          itL != Llist.end(); ++itL)
//      {
//          Eigen::Vector3d C = itL->first;
//          Eigen::Vector3d X = itL->second;
//          const Eigen::Vector3d l_il = X - C;
//          const T lanmda = ((*(Za)) - (T)C(2))/(T)l_il(2);

////          Eigen::Matrix<T,3,1>;
////          r << *Za, *Za, *Za;
//          const T c_il_w_0 = lanmda * (T)l_il(0) + (T)C(0);
//          const T c_il_w_1 = lanmda * (T)l_il(1) + (T)C(1);
//          const T c_il_w_2 = lanmda * (T)l_il(2) + (T)C(2);
//          T l_ilx2 = (T)l_il(0)*(T)l_il(0);
//          T l_ily2 = (T)l_il(1)*(T)l_il(1);
//          T sin_in = sqrt((l_ilx2 + l_ily2)/(l_ilx2 + l_ily2 + l_il(2)*l_il(2)));
//          T ref_l_z = -sqrt( (l_ilx2 + l_ily2)/(sin_in*sin_in / (T)(WATER_N*WATER_N)) - l_ilx2 - l_ily2);
//          // get reslut
////          Eigen::Vector4d thisLinePlane1, thisLinePlane2;
//          Eigen::Matrix<T,4,1> thisLinePlane1, thisLinePlane2;
//          thisLinePlane1 << (T)l_il(1), (T)-l_il(0),(T)0.0, (T)l_il(0)*(T)c_il_w_1-(T)l_il(1)*(T)c_il_w_0;
//          thisLinePlane2 << (T)0.0, (T)ref_l_z, (T)-l_il(1), (T)l_il(1)*(T)c_il_w_2-(T)ref_l_z*(T)c_il_w_1;
//          std::pair<Eigen::Matrix<T,4,1>, Eigen::Matrix<T,4,1>> thisLine;
////              std::cout << sqrt(thisLinePlane1(0)*thisLinePlane1(0) + thisLinePlane1(1)*thisLinePlane1(1) + thisLinePlane1(2)*thisLinePlane1(2)) << std::endl;
////              std::cout << sqrt(thisLinePlane2(0)*thisLinePlane2(0) + thisLinePlane2(1)*thisLinePlane2(1) + thisLinePlane2(2)*thisLinePlane2(2)) << std::endl;
//          thisLine.first = thisLinePlane1 / sqrt(thisLinePlane1(0)*thisLinePlane1(0) + thisLinePlane1(1)*thisLinePlane1(1) + thisLinePlane1(2)*thisLinePlane1(2));
//          thisLine.second = thisLinePlane2 / sqrt(thisLinePlane2(0)*thisLinePlane2(0) + thisLinePlane2(1)*thisLinePlane2(1) + thisLinePlane2(2)*thisLinePlane2(2));

//          AllLine.push_back(thisLine);

//      }
//      //
////      Eigen::MatrixXd A(AllLine.size()*2, 4);
//      const int ws = AllLine.size()*2;
//      Eigen::Matrix<T, Eigen::Dynamic, 4> A;
////      Eigen::Matrix<T, 2, 4> A2;
////      Eigen::SparseMatrix<T,AllLine.size()*2, 4> A;
//      A = Eigen::Matrix<T, Eigen::Dynamic, 4>::Zero(AllLine.size()*2, 4);
//      int iter = 0;

////      typename std::list<std::pair<Eigen::Matrix<T,4,1>, Eigen::Matrix<T,4,1> > >::iterator my;// my;// itAllL1;// = AllLine.begin();

//      for(typename std::list<std::pair<Eigen::Matrix<T,4,1>, Eigen::Matrix<T,4,1> > >::iterator itAllL = AllLine.begin();
//          itAllL != AllLine.end(); ++itAllL, ++iter)//std::size_t iter = 0; iter < allLines.size(); ++iter)
//      {

//          A(iter*2, 0) = itAllL->first(0);
//          A(iter*2, 1) = itAllL->first(1);
//          A(iter*2, 2) = itAllL->first(2);
//          A(iter*2, 3) = itAllL->first(3);
//          //
//          A(iter*2+1, 0) = itAllL->second(0);
//          A(iter*2+1, 1) = itAllL->second(1);
//          A(iter*2+1, 2) = itAllL->second(2);
//          A(iter*2+1, 3) = itAllL->second(3);

//      }
//      A/1;
////      // do SVD
////      Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, 4>> svd(A, Eigen::ComputeThinU|Eigen::ComputeThinV);
////      Eigen::Matrix<T, 4, Eigen::Dynamic> V = svd.matrixV();
////      int col = 3;
////      Eigen::Matrix<T, 3,1> resX = Eigen::Matrix<T, 3,1>(V(0,col)/V(3,col), V(1,col)/V(3,col), V(2,col)/V(3,col));

////      double meanL = 0;
////      for(std::list<std::pair<Eigen::Vector3d, Eigen::Vector3d> >::const_iterator itL = Llist.begin();
////          itL != Llist.end(); ++itL)
////      {
////          Eigen::Vector3d C = itL->first;
////          Eigen::Vector3d X = itL->second;
////          Eigen::Vector3d l = X-C;
////          double D = (l.cross3(resX)).norm();
////          double nL = l.norm();
////          double L = D / nL;



////          meanL += L/Llist.size();
////      }



////      residuals[0] = T(meanL);//*(resI-resJ);

//    return true;
//  }
//};

//struct WaterPlaneCostFunction
//{
//  const std::vector<double> CI_Inf;
//  const std::vector<double> CJ_Inf;

//  WaterPlaneCostFunction
//  (
//    const std::vector<double> & _CI_Inf,
//    const std::vector<double> & _CJ_Inf
//  ): CI_Inf(_CI_Inf), CJ_Inf(_CJ_Inf)
//  {
//  }

//  template <typename T> bool
//  operator()
//  (
//    const T* const Za,
//    T* residuals
//  )
//  const
//  {
//      const T rpI = T(CI_Inf[0]);
//      const T ZcI = T(CI_Inf[1]);
//      const T ZpI_ = T(CI_Inf[2]);
//      const T lI_xy = T(CI_Inf[3]);
//      const T lI_z = T(CI_Inf[4]);
//      const T lI = lI_xy / lI_z;

//      const T rpJ = T(CJ_Inf[0]);
//      const T ZcJ = T(CJ_Inf[1]);
//      const T ZpJ_ = T(CJ_Inf[2]);
//      const T lJ_xy = T(CJ_Inf[3]);
//      const T lJ_z = T(CJ_Inf[4]);
//      const T lJ = lJ_xy / lJ_z;

//      const T water_n = T(1.34);
//      const T resI = ((water_n*water_n-T(1.0))*(rpI-(*Za-ZcI)*lI)*(rpI-(*Za-ZcI)*lI) + water_n*water_n*ZpI_*ZpI_);
//      const T resJ = ((water_n*water_n-T(1.0))*(rpJ-(*Za-ZcJ)*lJ)*(rpJ-(*Za-ZcJ)*lJ) + water_n*water_n*ZpJ_*ZpJ_);

//      residuals[0] = abs(resI-resJ);//*(resI-resJ);

//    return true;
//  }
//};
///
///
///
///
///
///
///
struct WaterPlaneCostFunction2
{
  std::vector<double> CI_Inf;
  std::vector<double> CJ_Inf;
  Vec3 CI, CJ;

  WaterPlaneCostFunction2
  (
    const std::vector<double> & _CI_Inf,
    const std::vector<double> & _CJ_Inf,
    const Vec3 & _CI,
    const Vec3 & _CJ
  ): CI_Inf(_CI_Inf), CJ_Inf(_CJ_Inf), CI(_CI), CJ(_CJ)
  {
  }

  template <typename T> bool
  operator()
  (
      const T* const Za,
      const T* const P_xy,
      T* residuals
  )
  const
  {
      const T rpI = sqrt((P_xy[0]-CI(0))*(P_xy[0]-CI(0)) + (P_xy[1]-CI(1))*(P_xy[1]-CI(1)));
//      const T rpI = T(CI_Inf[0]);
      const T ZcI = T(CI_Inf[1]);
      const T ZpI_ = T(CI_Inf[2]);
      const T lI_xy = T(CI_Inf[3]);
      const T lI_z = T(CI_Inf[4]);
      const T lI = lI_xy / lI_z;
      //
      const T rpJ = sqrt((P_xy[0]-CJ(0))*(P_xy[0]-CJ(0)) + (P_xy[1]-CJ(1))*(P_xy[1]-CJ(1)));
//      const T rpJ = T(CJ_Inf[0]);
      const T ZcJ = T(CJ_Inf[1]);
      const T ZpJ_ = T(CJ_Inf[2]);
      const T lJ_xy = T(CJ_Inf[3]);
      const T lJ_z = T(CJ_Inf[4]);
      const T lJ = lJ_xy / lJ_z;

      const T water_n = T(1.34);
      const T resI = ((water_n*water_n-T(1.0))*(rpI-(*Za-ZcI)*lI)*(rpI-(*Za-ZcI)*lI) + water_n*water_n*ZpI_*ZpI_);
      const T resJ = ((water_n*water_n-T(1.0))*(rpJ-(*Za-ZcJ)*lJ)*(rpJ-(*Za-ZcJ)*lJ) + water_n*water_n*ZpJ_*ZpJ_);

//      residuals[0] = abs(resI-resJ);//*(resI-resJ);
//      residuals[0] = resI-resJ;//*(resI-resJ);
      residuals[0] = abs(resI-resJ);//*(resI-resJ);


    return true;
  }
};
//
//
//
//
//
//
//
//
struct WaterPlaneCostFunction3
{
  std::vector<double> CI_Inf;
  std::vector<double> CJ_Inf;
  Vec3 imgIFeatC, imgJFeatC;

  WaterPlaneCostFunction3
  (
    const std::vector<double> & _CI_Inf,
    const std::vector<double> & _CJ_Inf,
    const Vec3 & _imgIFeatC,
    const Vec3 & _imgJFeatC
  ): CI_Inf(_CI_Inf), CJ_Inf(_CJ_Inf), imgIFeatC(_imgIFeatC), imgJFeatC(_imgJFeatC)
  {
  }

  template <typename T> bool
  operator()
  (
      const T* const Za,
      const T* const P_xy,
      const T* const cam_extrinsics_I, // R_t
      const T* const cam_extrinsics_J, // R_t
      T* residuals
  )
  const
  {
      // CI PI
      const T * cam_R_I = &cam_extrinsics_I[0];
      const T * cam_t_I = &cam_extrinsics_I[3];
      const T cam_R_transpose_I[3] = {-cam_R_I[0], -cam_R_I[1], -cam_R_I[2]};
      T pose_center_I[3];
      ceres::AngleAxisRotatePoint(cam_R_transpose_I, cam_t_I, pose_center_I);
      pose_center_I[0] *= T(-1);
      pose_center_I[1] *= T(-1);
      pose_center_I[2] *= T(-1);
      const T temp_CL_I[3] = {T(imgIFeatC(0))-cam_t_I[0], T(imgIFeatC(1))-cam_t_I[1], T(imgIFeatC(2))-cam_t_I[2]};
      T PCI[3];
      ceres::AngleAxisRotatePoint(cam_R_transpose_I, temp_CL_I, PCI);
      // CJ PI
      const T * cam_R_J = &cam_extrinsics_J[0];
      const T * cam_t_J = &cam_extrinsics_J[3];
      const T cam_R_transpose_J[3] = {-cam_R_J[0], -cam_R_J[1], -cam_R_J[2]};
      T pose_center_J[3];
      ceres::AngleAxisRotatePoint(cam_R_transpose_J, cam_t_J, pose_center_J);
      pose_center_J[0] *= T(-1);
      pose_center_J[1] *= T(-1);
      pose_center_J[2] *= T(-1);
      const T temp_CL_J[3] = {T(imgJFeatC(0))-cam_t_J[0], T(imgJFeatC(1))-cam_t_J[1], T(imgJFeatC(2))-cam_t_J[2]};
      T PCJ[3];
      ceres::AngleAxisRotatePoint(cam_R_transpose_J, temp_CL_J, PCJ);
      //
      // others
      const T rpI = sqrt((PCI[0]-pose_center_I[0])*(PCI[0]-pose_center_I[0]) + (PCI[1]-pose_center_I[1])*(PCI[1]-pose_center_I[1]));
//      const T rpI = T(CI_Inf[0]);
      const T ZcI = T(CI_Inf[1]);
      const T ZpI_ = T(CI_Inf[2]);
      const T lI_xy = T(CI_Inf[3]);
      const T lI_z = T(CI_Inf[4]);
      const T lI = lI_xy / lI_z;
      //
      const T rpJ = sqrt((PCJ[0]-pose_center_I[0])*(PCJ[0]-pose_center_I[0]) + (PCJ[1]-pose_center_I[1])*(PCJ[1]-pose_center_I[1]));
//      const T rpJ = T(CJ_Inf[0]);
      const T ZcJ = T(CJ_Inf[1]);
      const T ZpJ_ = T(CJ_Inf[2]);
      const T lJ_xy = T(CJ_Inf[3]);
      const T lJ_z = T(CJ_Inf[4]);
      const T lJ = lJ_xy / lJ_z;

      const T water_n = T(1.34);
      const T resI = ((water_n*water_n-T(1.0))*(rpI-(*Za-ZcI)*lI)*(rpI-(*Za-ZcI)*lI) + water_n*water_n*ZpI_*ZpI_);
      const T resJ = ((water_n*water_n-T(1.0))*(rpJ-(*Za-ZcJ)*lJ)*(rpJ-(*Za-ZcJ)*lJ) + water_n*water_n*ZpJ_*ZpJ_);

//      residuals[0] = abs(resI-resJ);//*(resI-resJ);
//      residuals[0] = resI-resJ;//*(resI-resJ);
      residuals[0] = abs(resI-resJ);//*(resI-resJ);


    return true;
  }
};


struct PoseCenterConstraintCostFunction
{
  Vec3 weight_;
  Vec3 pose_center_constraint_;

  PoseCenterConstraintCostFunction
  (
    const Vec3 & center,
    const Vec3 & weight
  ): weight_(weight), pose_center_constraint_(center)
  {
  }

  template <typename T> bool
  operator()
  (
    const T* const cam_extrinsics, // R_t
    T* residuals
  )
  const
  {
    const T * cam_R = &cam_extrinsics[0];
    const T * cam_t = &cam_extrinsics[3];
    const T cam_R_transpose[3] = {-cam_R[0], -cam_R[1], -cam_R[2]};

    T pose_center[3];
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R_transpose, cam_t, pose_center);
    pose_center[0] *= T(-1);
    pose_center[1] *= T(-1);
    pose_center[2] *= T(-1);

    residuals[0] = T(weight_[0]) * (pose_center[0] - T(pose_center_constraint_[0]));
    residuals[1] = T(weight_[1]) * (pose_center[1] - T(pose_center_constraint_[1]));
    residuals[2] = T(weight_[2]) * (pose_center[2] - T(pose_center_constraint_[2]));

    return true;
  }
};


int main(int argc, char **argv)
{
  using namespace std;
  std::cout << std::endl
//    << "-----------------------------------------------------------\n"
//    << "Global Structure from Motion:\n"
//    << "-----------------------------------------------------------\n"
//    << "Open Source implementation of the paper:\n"
//    << "\"Global Fusion of Relative Motions for "
//    << "Robust, Accurate and Scalable Structure from Motion.\"\n"
//    << "Pierre Moulon, Pascal Monasse and Renaud Marlet. "
//    << " ICCV 2013." << std::endl
//    << "------------------------------------------------------------"
    << std::endl;


  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sMatchesDir;
  std::string sOutDir = "";
  int iRotationAveragingMethod = int (ROTATION_AVERAGING_L1);//(ROTATION_AVERAGING_L2);
  int iTranslationAveragingMethod = int (TRANSLATION_AVERAGING_SOFTL1);
  std::string sIntrinsic_refinement_options = "NONE";
  bool b_use_motion_priors = false;
  std::string rtFilePath = "";
  bool bNeedSixteenth = true;

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('m', sMatchesDir, "matchdir") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('r', iRotationAveragingMethod, "rotationAveraging") );
  cmd.add( make_option('t', iTranslationAveragingMethod, "translationAveraging") );
  cmd.add( make_option('f', sIntrinsic_refinement_options, "refineIntrinsics") );
  cmd.add( make_switch('P', "prior_usage") );
  cmd.add( make_option('n', rtFilePath, "file path of R t file, n_gps_imu_321_x.res") );
  cmd.add( make_option('e', bNeedSixteenth, "if need resize images into one-sixteenth") );


  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch (const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_file] path to a SfM_Data scene\n"
    << "[-m|--matchdir] path to the matches that corresponds to the provided SfM_Data scene\n"
    << "[-o|--outdir] path where the output data will be stored\n"
    << "\n[Optional]\n"
    << "[-r|--rotationAveraging]\n"
      << "\t 1 -> L1 minimization\n"
      << "\t 2 -> L2 minimization (default)\n"
    << "[-t|--translationAveraging]:\n"
      << "\t 1 -> L1 minimization\n"
      << "\t 2 -> L2 minimization of sum of squared Chordal distances\n"
      << "\t 3 -> SoftL1 minimization (default)\n"
    << "[-f|--refineIntrinsics] Intrinsic parameters refinement option\n"
      << "\t ADJUST_ALL -> refine all existing parameters (default) \n"
      << "\t NONE -> intrinsic parameters are held as constant\n"
      << "\t ADJUST_FOCAL_LENGTH -> refine only the focal length\n"
      << "\t ADJUST_PRINCIPAL_POINT -> refine only the principal point position\n"
      << "\t ADJUST_DISTORTION -> refine only the distortion coefficient(s) (if any)\n"
      << "\t -> NOTE: options can be combined thanks to '|'\n"
      << "\t ADJUST_FOCAL_LENGTH|ADJUST_PRINCIPAL_POINT\n"
      <<      "\t\t-> refine the focal length & the principal point position\n"
      << "\t ADJUST_FOCAL_LENGTH|ADJUST_DISTORTION\n"
      <<      "\t\t-> refine the focal length & the distortion coefficient(s) (if any)\n"
      << "\t ADJUST_PRINCIPAL_POINT|ADJUST_DISTORTION\n"
      <<      "\t\t-> refine the principal point position & the distortion coefficient(s) (if any)\n"
    << "[-P|--prior_usage] Enable usage of motion priors (i.e GPS positions)\n"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  if (iRotationAveragingMethod < ROTATION_AVERAGING_L1 ||
      iRotationAveragingMethod > ROTATION_AVERAGING_L2 )  {
    std::cerr << "\n Rotation averaging method is invalid" << std::endl;
    return EXIT_FAILURE;
  }

  const cameras::Intrinsic_Parameter_Type intrinsic_refinement_options =
    cameras::StringTo_Intrinsic_Parameter_Type(sIntrinsic_refinement_options);
  if (intrinsic_refinement_options == static_cast<cameras::Intrinsic_Parameter_Type>(0) )
  {
    std::cerr << "Invalid input for Bundle Adjusment Intrinsic parameter refinement option" << std::endl;
    return EXIT_FAILURE;
  }

  if (iTranslationAveragingMethod < TRANSLATION_AVERAGING_L1 ||
      iTranslationAveragingMethod > TRANSLATION_AVERAGING_SOFTL1 )  {
    std::cerr << "\n Translation averaging method is invalid" << std::endl;
    return EXIT_FAILURE;
  }

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }


//  Save(sfm_data,
//    stlplus::create_filespec(sOutDir, "cloud_and_poses_testInit", ".ply"),
//    ESfM_Data(ALL));

  // Init the regions_type from the image describer file (used for image regions extraction)
//  using namespace openMVG::features;
//  const std::string sImage_describer = stlplus::create_filespec(sMatchesDir, "image_describer", "json");
//  std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
//  if (!regions_type)
//  {
//    std::cerr << "Invalid: "
//      << sImage_describer << " regions type file." << std::endl;
//    return EXIT_FAILURE;
//  }

  // Features reading
  std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
//  if (!feats_provider->load(sfm_data, sMatchesDir, regions_type)) {
//    std::cerr << std::endl
//      << "Invalid features." << std::endl;
//    return EXIT_FAILURE;
//  }
  // set feats matches extrinsics
  std::ifstream in_f, in_m, in_e;
  std::string f_path = sfm_data.s_root_path+"/x.txt";
  std::string m_path = sfm_data.s_root_path+"/m.txt";
  std::string e_path = sfm_data.s_root_path+"/RC.txt";
  in_f.open(f_path);
  in_m.open(m_path);
  in_e.open(e_path);
  if(!in_f)
  {
      std::cout << "open file " << f_path << std::endl;
      return EXIT_FAILURE;
  }
  if(!in_m)
  {
      std::cout << "open file " << m_path << std::endl;
      return EXIT_FAILURE;
  }
  if(!in_e)
  {
      std::cout << "open file " << e_path << std::endl;
      return EXIT_FAILURE;
  }
//  //
  const int imgSize = 18;//18;//18;
  int imgId = 0;
  while(!in_f.eof() && imgId < imgSize)
  {
      openMVG::features::PointFeatures eachImg;
      int thisImgId, featSize;
      in_f >> thisImgId >> featSize;
      for(std::size_t fId = 0; fId < featSize; ++fId)
      {
          double x, y;
          in_f >> x >> y;
          eachImg.push_back(openMVG::features::PointFeature(x,y));//[fId] = features::PointFeatures(x, y);
      }
      feats_provider->feats_per_view[thisImgId-1] = eachImg;
      ++imgId;
  }
  in_f.close();
  imgId = 0;
  while(!in_e.eof() && imgId < imgSize)
  {
      double R[9], C[3];

      in_e >> R[0] >> R[1] >> R[2] >> C[0]
           >> R[3] >> R[4] >> R[5] >> C[1]
           >> R[6] >> R[7] >> R[8] >> C[2];
      Mat3 poseR;
      poseR << R[0], R[1], R[2],
               R[3], R[4], R[5],
               R[6], R[7], R[8];
      Vec3 poseC;
      poseC << C[0], C[1], C[2];
      Mat3 transR;
      transR << 0.0, 1.0, 0.0,
                1.0, 0.0, 0.0,
                0.0, 0.0, -1.0;
//      sfm_data.poses[imgId] = Pose3(poseR, poseC);
//      std::cout << "error 1" << std::endl;
////      sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(sfm_data.GetViews().at(imgId).get());
////      View& view_it = sfm_data.GetViews().begin();
//      View& view_it = sfm_data.GetViews();
//      const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
////      sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(sfm_data.GetViews().at(imgId));
//      std::cout << "error 2" << std::endl;
////      prior->SetPoseCenterPrior(poseC, Vec3(1.0,1.0,1.0));
//      prior->center_weight_ = Vec3(1.0,1.0,1.0);
//      prior->pose_center_ = poseC;
//      std::cout << "error 3" << std::endl;
////      /sfm_data.views[imgId].SetPoseCenterPrior(poseC,Vec3(1.0,1.0,1.0));
      ++imgId;
  }
  in_e.close();

//  for()
//  {

//  }

  //set match
  matching::PairWiseMatches pairWise_matches_;
  const int pairSize = 149;//153;//149;//
  int pairId = 0;
  while (!in_m.eof() && pairId < pairSize)
  {

      int imgIdI, imgIdJ, listSize;

      in_m >> imgIdI >> imgIdJ >> listSize;
      Pair imgPair = Pair(imgIdI-1, imgIdJ-1);
      matching::IndMatches indM;
      for (std::size_t fId = 0; fId < listSize; ++fId)
      {
          int fIdI, fIdJ;
          in_m >> fIdI >> fIdJ;
          indM.push_back(matching::IndMatch(fIdI-1, fIdJ-1));
      }
      pairWise_matches_.insert(std::pair<Pair, matching::IndMatches>(imgPair, indM));
      //
      ++pairId;
  }
  in_m.close();
  std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
  matches_provider->pairWise_matches_ = pairWise_matches_;

//  //set match
//  matching::PairWiseMatches pairWise_matches_;
//  for (std::size_t imgIdI = 0; imgIdI < imgSize-1; ++imgIdI)
//  {
//      for (std::size_t imgIdJ = imgIdI+1; imgIdJ < imgSize; ++imgIdJ)
//      {
//          Pair imgPair = Pair(imgIdI, imgIdJ);
//          matching::IndMatches indM;
//          //
//          for (std::size_t fId = 0; fId < featSize; ++fId)
//          {
//              indM.push_back(matching::IndMatch(fId, fId));
//          }
//          //
//          pairWise_matches_.insert(std::pair<Pair, matching::IndMatches>(imgPair, indM));
//      }
//  }
//  std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
//  matches_provider->pairWise_matches_ = pairWise_matches_;










  // Matches reading
//  std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
//  if // Try to read the two matches file formats
//  (
//    !(matches_provider->load(sfm_data, stlplus::create_filespec(sMatchesDir, "matches.e.txt")) ||
//      matches_provider->load(sfm_data, stlplus::create_filespec(sMatchesDir, "matches.e.bin")))
//  )
//  {
//    std::cerr << std::endl
//      << "Invalid matches file." << std::endl;
//    return EXIT_FAILURE;
//  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutDir))
  {
    if (!stlplus::folder_create(sOutDir))
    {
      std::cerr << "\nCannot create the output directory" << std::endl;
    }
  }

  ///!!!
//  if(bNeedSixteenth)
//  {
//      for (auto itView = sfm_data.views.begin(); itView != sfm_data.views.end(); ++ itView)
//      {
//        const IndexT imgId = itView->first;
//        #pragma omp parallel for schedule(dynamic)
//        for(size_t i = 0; i < feats_provider->feats_per_view[imgId].size(); ++i)
//        {
//            openMVG::features::PointFeatures::iterator itPF = feats_provider->feats_per_view[imgId].begin();
//            std::advance(itPF, i);
//            (*itPF).x() *= 4.0;
//            (*itPF).y() *= 4.0;
//        }
//      }
//  }

  //---------------------------------------
  // Global SfM reconstruction process
  //---------------------------------------

  openMVG::system::Timer timer;
  GlobalSfMReconstructionEngine_RelativeMotions sfmEngine(
    sfm_data,
    sOutDir,
    stlplus::create_filespec(sOutDir, "Reconstruction_Report.html"));

  // Configure the features_provider & the matches_provider
  sfmEngine.SetFeaturesProvider(feats_provider.get());
  sfmEngine.SetMatchesProvider(matches_provider.get());

  // Configure reconstruction parameters
  sfmEngine.Set_Intrinsics_Refinement_Type(intrinsic_refinement_options);
  b_use_motion_priors = false;//true;// cmd.used('P');
  sfmEngine.Set_Use_Motion_Prior(b_use_motion_priors);

  // Configure motion averaging method
  sfmEngine.SetRotationAveragingMethod(
    ERotationAveragingMethod(iRotationAveragingMethod));
  sfmEngine.SetTranslationAveragingMethod(
    ETranslationAveragingMethod(iTranslationAveragingMethod));

  SfM_Data newdata, undiswaterdata;
  if (sfmEngine.Process_forSimulation())//_GCP_GPS_water(rtFilePath)) //Process_threeCameras_imu_gps_init_new())//Process_test(sMatchesDir))//_down_imu_gps(sMatchesDir))
  {
      {
          {
      //        SfM_Data newdata = sfmEngine.Get_SfM_Data();
              SfM_Data undiswaterdata = sfmEngine.Get_SfM_Data();
              std::ofstream out_e, out_k;
              std::string e_path = sfmEngine.Get_SfM_Data().s_root_path+"/RC_new_withWater.txt";//undisWater.txt";//
              std::string k_path = sfmEngine.Get_SfM_Data().s_root_path+"/K_new_withWater.txt";//undisWater.txt";//
              out_e.open(e_path);
              out_k.open(k_path);
              if(!out_e)
              {
                  std::cout << "open file " << e_path << std::endl;
                  return EXIT_FAILURE;
              }

              for(int imgId = 0; imgId < sfmEngine.Get_SfM_Data().poses.size(); ++imgId)
              {
                  const Pose3 pose = sfmEngine.Get_SfM_Data().poses.at(imgId);
                  Mat3 R = pose.rotation();
                  Vec3 C = pose.center();

                  out_e << R(0,0) << " " << R(0,1) << " " << R(0,2) << " " << C(0) << std::endl
                        << R(1,0) << " " << R(1,1) << " " << R(1,2) << " " << C(1) << std::endl
                        << R(2,0) << " " << R(2,1) << " " << R(2,2) << " " << C(2) << std::endl;

                  std::vector<double> cam_intrinsics = sfm_data.intrinsics[sfm_data.views[imgId]->id_intrinsic]->getParams();
                  out_k << cam_intrinsics[0] << " " << cam_intrinsics[1] << " " << cam_intrinsics[2] << std::endl;
              }
              out_e.close();
              out_k.close();

              //save for analyse error
              std::ofstream out_X, out_xm;
              out_X.open(sOutDir+"/X_withWater.txt");//undisWater.txt");//withWater.txt");//
              out_xm.open(sOutDir+"/xm_withWater.txt");//undisWater.txt");//withWater.txt");
              for(Landmarks::const_iterator itX = undiswaterdata.structure.begin();
                  itX != undiswaterdata.structure.end(); ++itX)
              {
                  out_xm << itX->second.obs.size() << " ";
                  for(Observations::const_iterator itObs = itX->second.obs.begin();
                      itObs != itX->second.obs.end(); ++itObs)
                  {
                      out_xm << itObs->first << " " << itObs->second.id_feat << " " << itObs->second.x(0) << " " << itObs->second.x(1) << " ";

                  }
                  out_xm << itX->second.X(0) << " "
                         << itX->second.X(1) << " "
                         << itX->second.X(2) << std::endl;
                  out_X << itX->second.obs.begin()->first << " "
                        << itX->second.obs.begin()->second.id_feat << " "
                        << itX->second.X(0) << " "
                        << itX->second.X(1) << " "
                        << itX->second.X(2) << std::endl;

              }
              out_X.close();
              out_xm.close();
          }
      }

      {
//      newdata = sfmEngine.Get_SfM_Data();
//      undiswaterdata = sfmEngine.Get_SfM_Data();
//      std::ofstream out_e;
//      std::string e_path = sfmEngine.Get_SfM_Data().s_root_path+"/RC_new_withWater.txt";//undisWater.txt";//
//      out_e.open(e_path);
//      if(!out_e)
//      {
//          std::cout << "open file " << e_path << std::endl;
//          return EXIT_FAILURE;
//      }

//      imgId = 0;
//      for(std::size_t imgId = 0; imgId < sfmEngine.Get_SfM_Data().poses.size(); ++imgId)
//      {
//          const Pose3 pose = sfmEngine.Get_SfM_Data().poses.at(imgId);
//          Mat3 R = pose.rotation();
//          Vec3 C = pose.center();

//          out_e << R(0,0) << " " << R(0,1) << " " << R(0,2) << " " << C(0) << std::endl
//                << R(1,0) << " " << R(1,1) << " " << R(1,2) << " " << C(1) << std::endl
//                << R(2,0) << " " << R(2,1) << " " << R(2,2) << " " << C(2) << std::endl;
//      }
//      out_e.close();

      }
    std::cout << std::endl << " Total Ac-Global-Sfm took (s): " << timer.elapsed() << std::endl;


    std::cout << "...Generating SfM_Report.html" << std::endl;
    Generate_SfM_Report(sfmEngine.Get_SfM_Data(),
      stlplus::create_filespec(sOutDir, "SfMReconstruction_Report.html"));

    //-- Export to disk computed scene (data & visualizable results)
    std::cout << "...Export SfM_Data to disk." << std::endl;
    Save(sfmEngine.Get_SfM_Data(),
      stlplus::create_filespec(sOutDir, "sfm_data", ".bin"),
      ESfM_Data(ALL));

    Save(sfmEngine.Get_SfM_Data(),
      stlplus::create_filespec(sOutDir, "cloud_and_poses", ".ply"),
      ESfM_Data(ALL));

    //save for analyse error
//    std::ofstream out_X, out_xm;
//    out_X.open(sOutDir+"/X_withWater.txt");//undisWater.txt");//withWater.txt");//
//    out_xm.open(sOutDir+"/xm_withWater.txt");//undisWater.txt");//withWater.txt");
//    for(Landmarks::const_iterator itX = undiswaterdata.structure.begin();
//        itX != undiswaterdata.structure.end(); ++itX)
//    {
//        out_xm << itX->second.obs.size() << " ";
//        for(Observations::const_iterator itObs = itX->second.obs.begin();
//            itObs != itX->second.obs.end(); ++itObs)
//        {
//            out_xm << itObs->first << " " << itObs->second.id_feat << " " << itObs->second.x(0) << " " << itObs->second.x(1) << " ";

//        }
//        out_xm << itX->second.X(0) << " "
//               << itX->second.X(1) << " "
//               << itX->second.X(2) << std::endl;
//        out_X << itX->second.obs.begin()->first << " "
//              << itX->second.obs.begin()->second.id_feat << " "
//              << itX->second.X(0) << " "
//              << itX->second.X(1) << " "
//              << itX->second.X(2) << std::endl;

//    }
//    out_X.close();
//    out_xm.close();
//    return EXIT_SUCCESS;
  }

  std::map<std::pair<int, int>, Vec3> P_Map;
  {
      //preparing P
      for(Landmarks::const_iterator itX = newdata.structure.begin();
          itX != newdata.structure.end(); ++itX)
      {
          Vec3 thisP = itX->second.X;
          for(Observations::const_iterator itObs = itX->second.obs.begin();
              itObs != itX->second.obs.end(); ++itObs)
          {
              std::pair<int, int> imgFeatPair;
              imgFeatPair.first = itObs->first;
              imgFeatPair.second = itObs->second.id_feat;
              //
              P_Map.insert(std::pair<std::pair<int, int>, Vec3>(imgFeatPair, thisP));

          }
      }

  }

  // underwater scene
  {

      const double n = 1.33;
      const double N = (n*n - 1.0) / (n*n);
      const double waterPlane = 0.4;//0.25;
      ceres::Problem problem;
      double resZa = 0.25;
      for(Landmarks::const_iterator itX = newdata.structure.begin();
          itX != newdata.structure.end(); ++itX)
      {
          Vec3 thisP = itX->second.X;
          if(thisP(2)>waterPlane)
          {
              continue;
          }

//          for(Observations::const_iterator itObs = itX->second.obs.begin();
//              itObs != itX->second.obs.end(); ++itObs)
//          {
//              std::pair<int, int> imgFeatPair;
//              imgFeatPair.first = itObs->first;
//              imgFeatPair.second = itObs->second.id_feat;
//              //
//              P_Map.insert(std::pair<std::pair<int, int>, Vec3>(imgFeatPair, thisP));

//          }

          //prepare line
//          std::vector<std::pair<Eigen::Vector4d, Eigen::Vector4d> > allLines(0); //lanmda*l + p ; <l, p>
          std::list<std::pair<Eigen::Vector4d, Eigen::Vector4d> > allLines;
          std::list<std::pair<Eigen::Vector3d, Eigen::Vector3d> > Llist;
//          std::cout << allLines.size() << std::endl;
          for(Observations::const_iterator itObs = itX->second.obs.begin();
              itObs != itX->second.obs.end(); ++itObs)
          {
              const Pose3 pose = newdata.poses[itObs->first];
              const Vec3 C = pose.center();
              Mat3 K;
              K << 3047.0, 0.0, 2000.0,
                   0.0, 3047.0, 1500.0,
                   0.0, 0.0, 1.0;
              const Vec3 X = pose.rotation().inverse()*(K.inverse() * Vec3(itObs->second.x(0), itObs->second.x(1), 1.0) - pose.translation());
              const Vec3 l_il = X - C;
              const double lanmda = (waterPlane - C(2))/l_il(2);
              const Vec3 c_il_w = lanmda * l_il + C;
              double l_ilx2 = l_il(0)*l_il(0);
              double l_ily2 = l_il(1)*l_il(1);
              double sin_in = sqrt((l_ilx2 + l_ily2)/(l_ilx2 + l_ily2 + l_il(2)*l_il(2)));
              double ref_l_z = -sqrt( (l_ilx2 + l_ily2)/(sin_in*sin_in / (WATER_N*WATER_N)) - l_ilx2 - l_ily2);
//              double ref_l_z = waterPlane-sqrt( (l_ilx2 + l_ily2)/(sin_in*sin_in / (WATER_N*WATER_N)) - l_ilx2 - l_ily2);
              // get reslut
              Eigen::Vector4d thisLinePlane1, thisLinePlane2;
              thisLinePlane1 << l_il(1), -l_il(0), 0.0, l_il(0)*c_il_w(1)-l_il(1)*c_il_w(0);
              thisLinePlane2 << 0.0, ref_l_z, -l_il(1), l_il(1)*c_il_w(2)-ref_l_z*c_il_w(1);
              std::pair<Eigen::Vector4d, Eigen::Vector4d> thisLine;
//              std::cout << sqrt(thisLinePlane1(0)*thisLinePlane1(0) + thisLinePlane1(1)*thisLinePlane1(1) + thisLinePlane1(2)*thisLinePlane1(2)) << std::endl;
//              std::cout << sqrt(thisLinePlane2(0)*thisLinePlane2(0) + thisLinePlane2(1)*thisLinePlane2(1) + thisLinePlane2(2)*thisLinePlane2(2)) << std::endl;
              thisLine.first = thisLinePlane1 / sqrt(thisLinePlane1(0)*thisLinePlane1(0) + thisLinePlane1(1)*thisLinePlane1(1) + thisLinePlane1(2)*thisLinePlane1(2));
              thisLine.second = thisLinePlane2 / sqrt(thisLinePlane2(0)*thisLinePlane2(0) + thisLinePlane2(1)*thisLinePlane2(1) + thisLinePlane2(2)*thisLinePlane2(2));

              allLines.push_back(thisLine);

              //
              std::pair<Eigen::Vector3d, Eigen::Vector3d> thisL;
              thisL.first = C;
              thisL.second = X;
              Llist.push_back(thisL);
          }
          // get resX
          {
//              Eigen::MatrixX4d A(allLines.size()*2, 4);
              Eigen::MatrixXd A(allLines.size()*2, 4);
              A = Eigen::MatrixXd::Zero(allLines.size()*2, 4);
              int iter = 0;
              for(std::list<std::pair<Eigen::Vector4d, Eigen::Vector4d> >::const_iterator itAllL = allLines.begin();
                  itAllL != allLines.end(); ++itAllL, ++iter)//std::size_t iter = 0; iter < allLines.size(); ++iter)
              {

                  A(iter*2, 0) = itAllL->first(0);
                  A(iter*2, 1) = itAllL->first(1);
                  A(iter*2, 2) = itAllL->first(2);
                  A(iter*2, 3) = itAllL->first(3);
                  //
                  A(iter*2+1, 0) = itAllL->second(0);
                  A(iter*2+1, 1) = itAllL->second(1);
                  A(iter*2+1, 2) = itAllL->second(2);
                  A(iter*2+1, 3) = itAllL->second(3);

              }
              // do SVD
              Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU|Eigen::ComputeThinV);
              Eigen::Matrix4Xd V = svd.matrixV();
    //                std::cout << " svd.matrixV() : " << svd.matrixV() << std::endl;
    //                std::cout << " V : " << std::endl << V << std::endl;
              int col = 3;//allLines.size()*2-1;
    //                std::cout << " col : " << col << std::endl;
              Vec3 resX = Vec3(V(0,col)/V(3,col), V(1,col)/V(3,col), V(2,col)/V(3,col));

              undiswaterdata.structure[itX->first].X = resX;
          }

          //for water plane
          {
//              ceres::CostFunction * cost_function =
//                new ceres::AutoDiffCostFunction<WaterPlaneCostFunction, 1, 1>(
//                  new WaterPlaneCostFunction(Llist));
//              problem.AddResidualBlock(cost_function, nullptr, &resZa);
          }

      }


      Save(undiswaterdata,
        stlplus::create_filespec(sOutDir, "cloud_and_poses_undiswater", ".ply"),
        ESfM_Data(ALL));

      //save for analyse error
      std::ofstream out_X, out_xm;
      out_X.open(sOutDir+"/X_undisWater.txt");
      out_xm.open(sOutDir+"/xm_undisWater.txt");
      for(Landmarks::const_iterator itX = undiswaterdata.structure.begin();
          itX != undiswaterdata.structure.end(); ++itX)
      {
          out_xm << itX->second.obs.size() << " ";
          for(Observations::const_iterator itObs = itX->second.obs.begin();
              itObs != itX->second.obs.end(); ++itObs)
          {
              out_xm << itObs->first << " " << itObs->second.id_feat << " " << itObs->second.x(0) << " " << itObs->second.x(1) << " ";

          }
          out_xm << itX->second.X(0) << " "
                 << itX->second.X(1) << " "
                 << itX->second.X(2) << std::endl;
          out_X << itX->second.obs.begin()->first << " "
                << itX->second.obs.begin()->second.id_feat << " "
                << itX->second.X(0) << " "
                << itX->second.X(1) << " "
                << itX->second.X(2) << std::endl;

      }
      out_X.close();
      out_xm.close();



  }
  return EXIT_SUCCESS;





  {
      double n = 1.33;
      double N = (n*n - 1.0) / (n*n);

      matching::PairWiseMatches matchesPairs = matches_provider->pairWise_matches_;

      //std::map<std::pair<>, >
//      newdata.structure;
      for(std::size_t idMP = 0; idMP < matchesPairs.size(); ++idMP)
      {
          matching::PairWiseMatches::const_iterator itMP = matchesPairs.begin();
          std::advance(itMP, idMP);
          int imgI = itMP->first.first;
          int imgJ = itMP->first.second;
//          std::cout << "imgI : " << imgI << " imgJ : " << imgJ << std::endl;
          //
          /// get features
          std::vector<double> WPlanList;
          std::vector<double> WPlanList_;
          std::vector<double> WPlanList__;
          //
          std::vector<std::vector<double> > imgIInfoList, imgJInfoList;
          std::vector<Vec3> imgICList, imgJCList;
          std::vector<std::vector<double>> PxyList;
          std::vector<std::vector<double>> EIList, EJList;
          std::vector<Vec3> imgPCIList, imgPCJList;

          for(std::size_t idMF = 0; idMF < itMP->second.size(); ++idMF)
          {
              matching::IndMatches::const_iterator itMF = itMP->second.begin();
              std::advance(itMF, idMF);
              //
              int featI = itMF->i_;
              int featJ = itMF->j_;
              //
              openMVG::features::PointFeature fI = feats_provider->feats_per_view[imgI][featI];
              openMVG::features::PointFeature fJ = feats_provider->feats_per_view[imgJ][featJ];
              //
              // I
              std::shared_ptr<cameras::IntrinsicBase> itIntrinsicI = newdata.intrinsics[newdata.views[imgI]->id_intrinsic];
              std::vector<double> intrinsicI = itIntrinsicI->getParams();
              Pose3 poseI = newdata.poses[imgI];
              Mat3 RI = poseI.rotation();
              Vec3 tI = poseI.translation();
              Vec3 CI = poseI.center();
              Vec3 reprojI = RI.transpose()*(Vec3(fI.x()-intrinsicI[1], fI.x()-intrinsicI[2], intrinsicI[0]) / intrinsicI[0] - tI);
              Vec2 lI = Vec2(sqrt((reprojI(0)-CI(0))*(reprojI(0)-CI(0)) + (reprojI(1)-CI(1))*(reprojI(1)-CI(1))), reprojI(2)-CI(2));
              Vec3 CLI = reprojI-CI;
              // J
              std::shared_ptr<cameras::IntrinsicBase> itIntrinsicJ = newdata.intrinsics[newdata.views[imgJ]->id_intrinsic];
              std::vector<double> intrinsicJ = itIntrinsicJ->getParams();
              Pose3 poseJ = newdata.poses[imgJ];
              Mat3 RJ = poseJ.rotation();
              Vec3 tJ = poseJ.translation();
              Vec3 CJ = poseJ.center();
              Vec3 reprojJ = RJ.transpose()*(Vec3(fJ.x()-intrinsicJ[1], fJ.x()-intrinsicJ[2], intrinsicJ[0]) / intrinsicJ[0] - tJ);
              Vec2 lJ = Vec2(sqrt((reprojJ(0)-CJ(0))*(reprojJ(0)-CJ(0)) + (reprojJ(1)-CJ(1))*(reprojJ(1)-CJ(1))), reprojJ(2)-CJ(2));
              Vec3 CLJ = reprojJ-CJ;
              //
              ///get P_xy
              ///
//                Eigen::MatrixXd functionA;
//                functionA.setZero(4,3);
//                functionA << CLI(1), -CLI(0), 0.0,
//                             CLI(2), 0.0, -CLI(0),
//                             CLJ(1), -CLJ(0), 0.0,
//                             CLJ(2), 0.0, -CLJ(0);
//                Eigen::Vector4d function_b;
////                Eigen::VectorXd function_b;
//                function_b << CLI(0)*CI(1)-CLI(1)*CI(0),
//                              CLI(0)*CI(2)-CLI(2)*CI(0),
//                              CLJ(0)*CJ(1)-CLJ(1)*CJ(0),
//                              CLJ(0)*CJ(2)-CLJ(2)*CJ(0);
//                functionA.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(-function_b);
//                Eigen::Vector3d testP = functionA.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(-function_b);
//                if(testP(2) > 0.25)
//                {
//                    continue;
//                }
//                const Vec2 P_xy = Vec2(testP(0), testP(1));
              //
//                openMVG::Triangulation trian_obj;
//                trian_obj.add(
//                  itIntrinsicI->get_projective_equivalent(poseI),
//                  Vec2(fI.x(), fI.y()));
//                trian_obj.add(
//                  itIntrinsicJ->get_projective_equivalent(poseJ),
//                  Vec2(fJ.x(), fJ.y()));
//                //
//                Vec3 X = trian_obj.compute();
//                const Vec2 P_xy = Vec2(X(0), X(1));
//                if(X(2) > 0.25)
//                {
//                    continue;
//                }
              //
              std::map<std::pair<int, int>, Vec3>::const_iterator itPI = P_Map.find(std::pair<int, int>(imgI, featI));
              std::map<std::pair<int, int>, Vec3>::const_iterator itPJ = P_Map.find(std::pair<int, int>(imgJ, featJ));
              if(itPI == P_Map.end() || itPJ == P_Map.end())
              {
                  continue;
              }
              const Vec2 P_xy = Vec2(itPI->second(0), itPI->second(1));
              if(itPI->second(2) > 0.25)
              {
                  continue;
              }
              //
//              if(itPI->second(0) < -0.366938 || itPI->second(0) > 0.154272 || itPI->second(1) < -1.183584 || itPI->second(1) > -0.388861)
////                if(itP->second(2) < -1.5 || itP->second(2) > -1.2)
//              {
//                  continue;

//              }
//              if(P_xy(0) < -0.366938 || P_xy(0) > 0.154272 || P_xy(1) < -1.183584 || P_xy(1) > -0.388861)
////                if(itP->second(2) < -1.5 || itP->second(2) > -1.2)
//              {
//                  continue;

//              }

              lI(0) = sqrt((P_xy(0)-CI(0))*(P_xy(0)-CI(0)) + (P_xy(1)-CI(1))*(P_xy(1)-CI(1)));
              lJ(0) = sqrt((P_xy(0)-CJ(0))*(P_xy(0)-CJ(0)) + (P_xy(1)-CJ(1))*(P_xy(1)-CJ(1)));
              // I
              double rpI = sqrt((P_xy(0)-CI(0))*(P_xy(0)-CI(0)) + (P_xy(1)-CI(1))*(P_xy(1)-CI(1)));
              double ZpI_ = CI(2) + lI(1)*rpI/lI(0);
              // J
              double rpJ = sqrt((P_xy(0)-CJ(0))*(P_xy(0)-CJ(0)) + (P_xy(1)-CJ(1))*(P_xy(1)-CJ(1)));
              double ZpJ_ = CJ(2) + lJ(1)*rpJ/lJ(0);
              //
              {

                  std::vector<double> imgJInfo;
                  imgJInfo.push_back(rpJ);
                  imgJInfo.push_back(CJ(2));
                  imgJInfo.push_back(ZpJ_);
                  imgJInfo.push_back(lJ(0));
                  imgJInfo.push_back(lJ(1));
                  imgJInfoList.push_back(imgJInfo);
                  std::vector<double> imgIInfo;
                  imgIInfo.push_back(rpI);
                  imgIInfo.push_back(CI(2));
                  imgIInfo.push_back(ZpI_);
                  imgIInfo.push_back(lI(0));
                  imgIInfo.push_back(lI(1));
                  imgIInfoList.push_back(imgIInfo);

                  imgJCList.push_back(CJ);
                  imgICList.push_back(CI);

                  std::vector<double> tempPxy;
                  tempPxy.push_back(P_xy(0));
                  tempPxy.push_back(P_xy(1));
                  PxyList.push_back(tempPxy);

                  double angleAxisI[3];
                  ceres::RotationMatrixToAngleAxis((const double*)RI.data(), angleAxisI);
                  std::vector<double> EIdata = {angleAxisI[0], angleAxisI[1], angleAxisI[2], tI(0), tI(1), tI(2)};
                  EIList.push_back(EIdata);
                  double angleAxisJ[3];
                  ceres::RotationMatrixToAngleAxis((const double*)RJ.data(), angleAxisJ);
                  std::vector<double> EJdata = {angleAxisJ[0], angleAxisJ[1], angleAxisJ[2], tI(0), tI(1), tI(2)};
                  EJList.push_back(EJdata);

                  imgPCIList.push_back((Vec3(fI.x()-intrinsicI[1], fI.x()-intrinsicI[2], intrinsicI[0]) / intrinsicI[0] - tI));
                  imgPCJList.push_back((Vec3(fJ.x()-intrinsicJ[1], fJ.x()-intrinsicJ[2], intrinsicJ[0]) / intrinsicJ[0] - tJ));


              }
              //
              // compute a b c
              double LI = lI(0)/lI(1);
              double LJ = lJ(0)/lJ(1);
              double lIZCI = LI*CI(2);
              double lJZCJ = LJ*CJ(2);
              double a = N*(LI*LI - LJ*LJ);
              double b = 2.0*N*(LJ*(lJZCJ+rpJ) - LI*(lIZCI+rpI)) - 2.0*(ZpI_-ZpJ_);
              double c = N*(lIZCI*(lIZCI+2.0*rpI) - lJZCJ*(lJZCJ+2.0*rpJ) + (rpI+rpJ)*(rpI-rpJ)) + ZpI_*ZpI_ - ZpJ_*ZpJ_;
              //
              double b2_4ac = b*b - 4.0*a*c;
              if(b2_4ac < 0)
              {
                  std::cout << "<0 featI " << featI << " featJ " << featJ << std::endl;
              }else{
                  double res1 = (-b + sqrt(b2_4ac)) / (2.0*a);
                  double res2 = (-b - sqrt(b2_4ac)) / (2.0*a);
//                    if(res1 > testP(2) && res2 > testP(2))
//                    if(res1 < testP(2) && res2 < testP(2))
//                    if(res1 < X(2) && res2 < X(2))
//                    if(res1 < testP(2) && res2 < testP(2))
                  double maxZp = ZpI_ > ZpJ_ ? ZpI_ : ZpJ_;
                  double minZp = ZpI_ > ZpJ_ ? ZpJ_ : ZpI_;

                  double tempRes, tempRes_;
                  if ((res1<maxZp)&&(res1>minZp))
                  {
                      tempRes = res2;
                      tempRes_ = res1;
//                        WPlanList_.push_back(res1);
//                        WPlanList.push_back(res2);
                  }else if ((res2<maxZp)&&(res2>minZp))
                  {
                      tempRes = res1;
                      tempRes_ = res2;
//                        WPlanList_.push_back(res2);
//                        WPlanList.push_back(res1);
                  }else{
                      continue;
                  }
                  if(tempRes>itPI->second(2))
                  {
                      WPlanList_.push_back(tempRes);
                      WPlanList.push_back(tempRes_);

                  }
              }

          }
          //
          ///
          ceres::Problem problem;
          double resZa = 0.25;
          if(imgJInfoList.size() < 20)
          {
              continue;
          }
          std::cout << "imgI : " << imgI << " imgJ : " << imgJ << std::endl;
          std::cout << "size : " << imgJInfoList.size() << std::endl;
          for (std::size_t idInfo = 0; idInfo < imgJInfoList.size(); ++idInfo)
          {
            {
//              ceres::CostFunction * cost_function =
//                new ceres::AutoDiffCostFunction<WaterPlaneCostFunction, 1, 1>(
//                  new WaterPlaneCostFunction(imgIInfoList[idInfo], imgJInfoList[idInfo]));
//              problem.AddResidualBlock(cost_function, nullptr, &resZa);
            }

              {
                  ceres::CostFunction * cost_function =
                    new ceres::AutoDiffCostFunction<WaterPlaneCostFunction2, 1, 1, 2>(
                      new WaterPlaneCostFunction2(imgIInfoList[idInfo], imgJInfoList[idInfo], imgICList[idInfo], imgJCList[idInfo] ));
                  problem.AddResidualBlock(cost_function, nullptr, &resZa, &PxyList[idInfo][0]);
              }

              {
//                  ceres::CostFunction * cost_function =
//                    new ceres::AutoDiffCostFunction<WaterPlaneCostFunction3, 1, 1, 2, 6, 6>(
//                      new WaterPlaneCostFunction3(imgIInfoList[idInfo], imgJInfoList[idInfo], imgPCIList[idInfo], imgPCJList[idInfo] ));
//                  problem.AddResidualBlock(cost_function, nullptr, &resZa, &PxyList[idInfo][0], &EIList[idInfo][0], &EJList[idInfo][0]);

              }
          }
          ceres::Solver::Options ceres_config_options;
          ceres_config_options.max_num_iterations = 500;
          ceres_config_options.preconditioner_type = ceres::IDENTITY;
//            static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
          ceres_config_options.linear_solver_type = ceres::DENSE_SCHUR;
//            static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
          ceres_config_options.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
//            static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
//          ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
          ceres_config_options.logging_type = ceres::PER_MINIMIZER_ITERATION;//ceres::SILENT;
          ceres_config_options.num_threads = 60;//ceres_options_.nb_threads_;
//          ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
          ceres_config_options.parameter_tolerance = 0.000000001;//ceres_options_.parameter_tolerance_;

          ceres::Solver::Summary summary;
          ceres::Solve(ceres_config_options, &problem, &summary);

          // If no error, get back refined parameters
//          if (!summary.IsSolutionUsable())
//          {
//            std::cout << "Bundle Adjustment failed." << std::endl;
//            return false;
//          }
//          else // Solution is usable
          {
              std::cout << "resZa " << resZa << std::endl;

          }

      }

  }

  return EXIT_SUCCESS;
}
