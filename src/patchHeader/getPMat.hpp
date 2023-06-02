#pragma once
#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/cameras/Cameras_Common_command_line_helper.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_report.hpp"
#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

using namespace openMVG;
using namespace openMVG::sfm;

//2023-3-4 这个函数是chatGPT写的
Mat34 getProjectionMatrix(const SfM_Data& sfm_data, const IndexT id_view)
{
  // Get the view corresponding to the given id_view
  const View* view = sfm_data.GetViews().at(id_view).get();

  // Get the intrinsic and extrinsic camera parameters
  const auto intrinsic = sfm_data.GetIntrinsics().at(view->id_intrinsic);
  const geometry::Pose3 pose = sfm_data.GetPoseOrDie(view);

  // Compute the projection matrix P
  const Mat3& K = intrinsic->getIntrinsicMatrix();
  const Mat3 R = pose.rotation();
  const Vec3 t = pose.translation();
  Mat34 P;
  P.block<3,3>(0,0) = K * R;
  P.col(3) = K * t;
  return P;
}

//这也是chatGPT写的，这写法可真6
// 将Mat34类型数据按行优先存储在std::vector中
std::vector<float> mat34ToVector(const Mat34& mat) {
    std::vector<float> vec(12);
    vec[0] = mat(0,0);
    vec[1] = mat(0,1);
    vec[2] = mat(0,2);
    vec[3] = mat(0,3);
    vec[4] = mat(1,0);
    vec[5] = mat(1,1);
    vec[6] = mat(1,2);
    vec[7] = mat(1,3);
    vec[8] = mat(2,0);
    vec[9] = mat(2,1);
    vec[10] = mat(2,2);
    vec[11] = mat(2,3);
    return vec;
}
