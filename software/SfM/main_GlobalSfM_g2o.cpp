// g2o - General Graph Optimization
// Copyright (C) 2011 H. Strasdat
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <Eigen/StdVector>
#include <iostream>
#include <stdint.h>

#ifdef _MSC_VER
#include <unordered_set>
#else
#include <tr1/unordered_set>
#endif

#include "third_party/g2o/g2o/core/sparse_optimizer.h"
#include "third_party/g2o/g2o/core/block_solver.h"
#include "third_party/g2o/g2o/core/solver.h"
#include "third_party/g2o/g2o/core/robust_kernel_impl.h"
#include "third_party/g2o/g2o/core/optimization_algorithm_levenberg.h"
#include "third_party/g2o/g2o/solvers/cholmod/linear_solver_cholmod.h"
#include "third_party/g2o/g2o/solvers/dense/linear_solver_dense.h"
#include "third_party/g2o/g2o/types/sba/types_six_dof_expmap.h"
//#include "g2o/math_groups/se3quat.h"
#include "third_party/g2o/g2o/solvers/structure_only/structure_only_solver.h"

#include "third_party/g2o/g2o/core/optimization_algorithm_gauss_newton.h"
#include "third_party/g2o/g2o/solvers/csparse/linear_solver_csparse.h"
#include "third_party/g2o/g2o/solvers/eigen/linear_solver_eigen.h"
#include "third_party/g2o/g2o/solvers/pcg/linear_solver_pcg.h"

//#include "g2o/core/base_vertex.h"
//#include "g2o/core/base_binary_edge.h"
//#include "g2o/types/slam3d/se3_ops.h"
//#include <Eigen/Geometry>

#include<unistd.h>
#include<time.h>

using namespace Eigen;
using namespace std;

//class VertexSE3Expmap : public g2o::BaseVertex<6, g2o::SE3Quat>{
//public:
//  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

//  VertexSE3Expmap();

//  bool read(std::istream& is);

//  bool write(std::ostream& os) const;

//  virtual void setToOriginImpl() {
//    _estimate = SE3Quat();
//  }

//  virtual void oplusImpl(const double* update_)  {
//    Map<const Vector6d> update(update_);
//    setEstimate(SE3Quat::exp(update)*estimate());
//  }
//};

class MY_EdgeProjectXYZ2UV : public  g2o::BaseBinaryEdge<2, Vector2d, g2o::VertexSBAPointXYZ, g2o::VertexSE3Expmap>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  MY_EdgeProjectXYZ2UV();

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  void computeError()  {
    const g2o::VertexSBAPointXYZ* v0 = static_cast<const g2o::VertexSBAPointXYZ*>(_vertices[0]); //X
    const g2o::VertexSE3Expmap* v1 = static_cast<const g2o::VertexSE3Expmap*>(_vertices[1]); //Pose

//    std::pair<Vector2d, int> p2d(_measurement);
    Vector2d obs(_measurement(0), _measurement(1));
//    double _iId;//(_measurement(2));
    int iId = v1->id() / 10000 -1;

    const double * cam_intrinsics
        = (double* )(parameter(iId));

    Vector3d X_c = v1->estimate().map(v0->estimate());

    // Transform the point from homogeneous to euclidean (undistorted point)
    const double x_u = X_c(0) / X_c(2);
    const double y_u = X_c(1) / X_c(2);

    const double& focal = cam_intrinsics[0];
    const double& principal_point_x = cam_intrinsics[1];
    const double& principal_point_y = cam_intrinsics[2];
//    const double& k1 = cam_intrinsics[3];
//    const double& k2 = cam_intrinsics[4];
//    const double& k3 = cam_intrinsics[5];
//    const double& t1 = cam_intrinsics[6];
//    const double& t2 = cam_intrinsics[7];

    // Apply distortion (xd,yd) = disto(x_u,y_u)
//    const double r2 = x_u*x_u + y_u*y_u;
//    const double r4 = r2 * r2;
//    const double r6 = r4 * r2;
//    const double r_coeff = (1.0 + k1*r2 + k2*r4 + k3*r6);
//    const double t_x = t2 * (r2 + 2.0 * x_u*x_u) + 2.0 * t1 * x_u * y_u;
//    const double t_y = t1 * (r2 + 2.0 * y_u*y_u) + 2.0 * t2 * x_u * y_u;
//    const double x_d = x_u * r_coeff + t_x;
//    const double y_d = y_u * r_coeff + t_y;

    // Apply focal length and principal point to get the final image coordinates
    const double projected_x = principal_point_x + focal * x_u;
    const double projected_y = principal_point_y + focal * y_u;

    _error = obs - Vector2d(projected_x, projected_y);


  }

  virtual void linearizeOplus();

  g2o::CameraParameters * _cam;
};


class Sample {
public:
  static int uniform(int from, int to);
  static double uniform();
  static double gaussian(double sigma);
};

static double uniform_rand(double lowerBndr, double upperBndr){
  return lowerBndr + ((double) std::rand() / (RAND_MAX + 1.0)) * (upperBndr - lowerBndr);
}

static double gauss_rand(double mean, double sigma){
  double x, y, r2;
  do {
    x = -1.0 + 2.0 * uniform_rand(0.0, 1.0);
    y = -1.0 + 2.0 * uniform_rand(0.0, 1.0);
    r2 = x * x + y * y;
  } while (r2 > 1.0 || r2 == 0.0);
  return mean + sigma * y * std::sqrt(-2.0 * log(r2) / r2);
}

int Sample::uniform(int from, int to){
  return static_cast<int>(uniform_rand(from, to));
}

double Sample::uniform(){
  return uniform_rand(0., 1.);
}

double Sample::gaussian(double sigma){
  return gauss_rand(0., sigma);
}


struct ImgInfo{
//    std::vector<Vector2d> point2dList;
    int intrinsicId;
};
typedef std::map<int, ImgInfo> ImgList;

typedef std::map<int, Vector3d> XList;

struct LandmarkInfo{
    std::vector<std::pair<int, Vector2d> > imgfeatIdList; //<imgId, feat>
    Vector3d X;
};
typedef std::map<int, LandmarkInfo> Landmarks;

typedef std::map<int, std::pair<Eigen::Quaterniond, Vector3d> > Extrinsics; //R t



int main(int argc, const char* argv[]){

//  if (argc<2)
//  {
//    cout << endl;
//    cout << "Please type: " << endl;
//    cout << "ba_demo [PIXEL_NOISE] [OUTLIER RATIO] [ROBUST_KERNEL] [STRUCTURE_ONLY] [DENSE]" << endl;
//    cout << endl;
//    cout << "PIXEL_NOISE: noise in image space (E.g.: 1)" << endl;
//    cout << "OUTLIER_RATIO: probability of spuroius observation  (default: 0.0)" << endl;
//    cout << "ROBUST_KERNEL: use robust kernel (0 or 1; default: 0==false)" << endl;
//    cout << "STRUCTURE_ONLY: performe structure-only BA to get better point initializations (0 or 1; default: 0==false)" << endl;
//    cout << "DENSE: Use dense solver (0 or 1; default: 0==false)" << endl;
//    cout << endl;
//    cout << "Note, if OUTLIER_RATIO is above 0, ROBUST_KERNEL should be set to 1==true." << endl;
//    cout << endl;
//    exit(0);
//  }

//  double PIXEL_NOISE = atof(argv[1]);
//  double OUTLIER_RATIO = 0.0;

//  if (argc>2)  {
//    OUTLIER_RATIO = atof(argv[2]);
//  }

  bool ROBUST_KERNEL = false;
//  ROBUST_KERNEL = true;
//  if (argc>3){
//    ROBUST_KERNEL = atoi(argv[3]) != 0;
//  }
  bool STRUCTURE_ONLY = false;
//  STRUCTURE_ONLY = true;
//  if (argc>4){
//    STRUCTURE_ONLY = atoi(argv[4]) != 0;
//  }

  bool DENSE = false;
//  if (argc>5){
//    DENSE = atoi(argv[5]) != 0;
//  }

//  cout << "PIXEL_NOISE: " <<  PIXEL_NOISE << endl;
//  cout << "OUTLIER_RATIO: " << OUTLIER_RATIO<<  endl;
//  cout << "ROBUST_KERNEL: " << ROBUST_KERNEL << endl;
//  cout << "STRUCTURE_ONLY: " << STRUCTURE_ONLY<< endl;
//  cout << "DENSE: "<<  DENSE << endl;


  g2o::SparseOptimizer optimizer;
  optimizer.setVerbose(false);
  g2o::BlockSolver_6_3::LinearSolverType * linearSolver;
//  if (DENSE) {
//    linearSolver= new g2o::LinearSolverDense<g2o
//        ::BlockSolver_6_3::PoseMatrixType>();
//  } else {
//    linearSolver
//        = new g2o::LinearSolverCholmod<g2o
//        ::BlockSolver_6_3::PoseMatrixType>();
//  }


//  linearSolver = new g2o::LinearSolverCholmod<g2o::BlockSolver_6_3::PoseMatrixType>();
//  linearSolver = new g2o::LinearSolverCCS<g2o::BlockSolver_6_3::PoseMatrixType>();
//  linearSolver = new g2o::LinearSolverEigen<g2o::BlockSolver_6_3::PoseMatrixType>();

//  linearSolver = new g2o::LinearSolverDense<g2o::BlockSolver_6_3::PoseMatrixType>();

//  linearSolver = new g2o::LinearSolverCSparse<g2o::BlockSolver_6_3::PoseMatrixType>();
  linearSolver = new g2o::LinearSolverPCG<g2o::BlockSolver_6_3::PoseMatrixType>();

  g2o::BlockSolver_6_3 * solver_ptr
      = new g2o::BlockSolver_6_3(linearSolver);

//  cerr << "Using CSPARSE Levenberg" << endl;
  g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
//  g2o::OptimizationAlgorithmWithHessian* solver = new g2o::OptimizationAlgorithmWithHessian(solver_ptr);
  optimizer.setAlgorithm(solver);


  std::string root_path = "/home/guang/data/forG2OStep";
//  std::string root_path = "/home/guang/data/testWaterDatum/imgs";
//  std::string root_path = "/home/guang/data/forG2OStep2";

  std::string matches_path = root_path + "/matches/";

  std::cout << "root_path : " << root_path << std::endl;
//  std::cout << "sparse with Hessian " << root_path << std::endl;
  std::cout << "LinearSolverPCG2_ 133 -1 test pcg maxIter OptimizationAlgorithmLevenberg " << root_path << std::endl;
  ImgList imgList;
  /// load datum
  ///
  //
  /// load landmarks
  ///
  Landmarks landmarks;
  //
  std::ifstream in_landmark;
  std::string in_landmark_path = root_path + "/landmarks.res";
  in_landmark.open(in_landmark_path);
  if(!in_landmark)
  {
      std::cout << "open " << in_landmark_path << " file failed!" << std::endl;
      return EXIT_FAILURE;
  }
  int count_landmarks;
  in_landmark >> count_landmarks;
  int line_landmark = 0;
  while(line_landmark < count_landmarks && !in_landmark.eof())
  {
      int tempLandmarkId;
      double tempX_x, tempX_y, tempX_z;
      in_landmark >> tempLandmarkId >> tempX_x >> tempX_y >> tempX_z;
      Vector3d tempX(tempX_x, tempX_y, tempX_z);
      LandmarkInfo thisInfo;
      thisInfo.X = tempX;
      {
          int count_obs;
          in_landmark >> count_obs;
          int line_obs = 0;
          while(line_obs < count_obs)
          {
              int tempImgId;
              double tempFeat_x, tempFeat_y;
              in_landmark >> tempImgId >> tempFeat_x >> tempFeat_y;
              Vector2d tempFeat(tempFeat_x, tempFeat_y);
              thisInfo.imgfeatIdList.push_back(std::pair<int, Vector2d>(tempImgId, tempFeat));
              //
              ImgInfo thisimgInfo;
        //      thisimgInfo.point2dList = thisp2dList;
              thisimgInfo.intrinsicId = tempImgId / 100000 - 1;
              imgList[tempImgId] = thisimgInfo;
              //
              ++line_obs;
          }
      }
      landmarks[tempLandmarkId] = thisInfo;
      //
      ++line_landmark;
  }
  in_landmark.close();
  //
  /// load extrinsics
  ///
  //
  //extrinsics a
//  typedef std::map<int, std::pair<Eigen:: Quaterniond, Vector3d> > Extrinsics;
  std::cout << "loading extrinsics" << std::endl;
  Extrinsics extrinsics;
  std::ifstream in_extrinsics;
  std::string in_extrinsics_path = root_path + "/extrinsics.res";
  in_extrinsics.open(in_extrinsics_path);
  if(!in_extrinsics)
  {
      std::cout << "open " << in_extrinsics_path << " file failed!" << std::endl;
      return EXIT_FAILURE;
  }
  int count_extrinsics;
  in_extrinsics >> count_extrinsics;
  int line_extrinsic = 0;
  while(line_extrinsic<count_extrinsics && !in_extrinsics.eof())
  {
      int eId;
//      Eigen::Quaterniond q;
      double q0, q1, q2, q3;
      Vector3d t;
      in_extrinsics >> eId >> q0 >> q1 >> q2 >> q3 >> t(0) >> t(1) >> t(2);
      Eigen::Quaterniond q(q0 , q1, q2, q3);
      extrinsics[eId] = std::pair<Eigen::Quaterniond, Vector3d>(q,t);

      ++line_extrinsic;
  }
  //
  //
  /// prepare intrinsics
  ///
  std::cout << "preparing intrinsics" << std::endl;
//  double focal_length0 = 3047.0014970924713;
//  Vector2d principal_point0 (2000.0, 1500.0);
  double focal_length0 = 10285.958774;
  Vector2d principal_point0 (3657.823643, 2458.943385);
  double focal_length1 = 7282.094708;
  Vector2d principal_point1 (3655.785805, 2434.360780);
  double focal_length2 = 10303.589972;
  Vector2d principal_point2 (3657.823643, 2458.943385);
  //
  g2o::CameraParameters * cam_params0
      = new g2o::CameraParameters (focal_length0, principal_point0, 0.);
  cam_params0->setId(0);
  g2o::CameraParameters * cam_params1
      = new g2o::CameraParameters (focal_length1, principal_point1, 0.);
  cam_params1->setId(1);
  g2o::CameraParameters * cam_params2
      = new g2o::CameraParameters (focal_length2, principal_point2, 0.);
  cam_params2->setId(2);

  if (!optimizer.addParameter(cam_params0)) {
    assert(false);
  }
  if (!optimizer.addParameter(cam_params1)) {
    assert(false);
  }
  if (!optimizer.addParameter(cam_params2)) {
    assert(false);
  }

  /// prepare landmark
  ///
  //

  /// set vertex and edge
  ///
  /// set pose vertex
  ///
//  vector<g2o::SE3Quat,
//      aligned_allocator<g2o::SE3Quat> > true_poses;
//  int vertex_id = 0;
  std::cout << "preparing extrinsics" << std::endl;
//  std::map<int, int> PidMap;
  for (size_t eListId = 0; eListId < extrinsics.size(); ++eListId)
  {
      Extrinsics::const_iterator itE = extrinsics.begin();
      std::advance(itE, eListId);

      g2o::SE3Quat pose(itE->second.first, itE->second.second);
      g2o::VertexSE3Expmap *v_se3 = new g2o::VertexSE3Expmap();
      v_se3->setId(itE->first);
      v_se3->setEstimate(pose);
//      if(eListId < 2)
//      {
//          v_se3->setFixed(true);
//      }else{
//          v_se3->setFixed(false);
//      }
      optimizer.addVertex(v_se3);
//      PidMap[itE->first] = eListId;
//    true_poses.push_back(pose);
//      vertex_id++;
  }
  //
  /// set X vertex
  ///
  Extrinsics::const_iterator itLastE = extrinsics.end();
  --itLastE;
  const int point_idstart = itLastE->first+1;//extrinsics.size();
  {
//  std::map<int, int> XidMap;
//  for (std::size_t XListId = 0; XListId < X_List.size(); ++XListId)
//  {
//      XList::const_iterator itX = X_List.begin();
//      std::advance(itX, XListId);
//      //
//      g2o::VertexSBAPointXYZ * v_p = new g2o::VertexSBAPointXYZ();
//      v_p->setId(point_idstart + XListId);
//      v_p->setMarginalized(true);
//      v_p->setEstimate(itX->second + Vector3d(Sample::gaussian(1), Sample::gaussian(1), Sample::gaussian(1)));
//      optimizer.addVertex(v_p);
//      XidMap[itX->first] = point_idstart + XListId;
//  }
  //
  ///set landmark
  ///
//  double sum_diff2 = 0;
//  int point_num = 0;
//  long compateCount = 0;
//  #pragma omp parallel for
//  for(std::size_t lListId = 0; lListId < landmarks.size(); ++lListId)
//  {
//      #pragma omp parallel for
//      for(std::size_t pListId = 0; pListId < itL->second.imgfeatIdList.size(); ++pListId)
//      {
//          ++compateCount;
//      }

//  }

//  return 0;
  }

  std::cout << "preparing landmark" << std::endl;
  {
//  int num1 = 0;
  //#pragma omp parallel for
//  //  for(std::size_t lListId = 0; lListId < landmarks.size(); ++lListId)
//    for(std::size_t lListId = 0; lListId < 10000; ++lListId)
//    {

//  //      std::cout << "preparing landmark 1" << std::endl;
//        Landmarks::const_iterator itL = landmarks.begin();
//        std::advance(itL, lListId);
//  //      std::cout << "preparing landmark 2" << std::endl;
//        //
//        /// set X vertex
//        g2o::VertexSBAPointXYZ * v_p = new g2o::VertexSBAPointXYZ();
//        int thisXId = point_idstart + lListId;
//        v_p->setId(thisXId);
//        v_p->setMarginalized(true);
//        v_p->setEstimate(itL->second.X);// + Vector3d(Sample::gaussian(1), Sample::gaussian(1), Sample::gaussian(1)));
//        #pragma omp critical
//        {
//          optimizer.addVertex(v_p);
//  //        std::cout << "preparing landmark3" << std::endl;
//        }
//  //      XidMap[itX->first] = point_idstart + lListId;
//        //
//  //      int XId = itL->second.XId;
//        //#pragma omp parallel for
//        for(std::size_t pListId = 0; pListId < itL->second.imgfeatIdList.size(); ++pListId)
//        {
//  //          std::cout << "preparing landmark4" << std::endl;
//            std::vector<std::pair<int, Vector2d> >::const_iterator itp = itL->second.imgfeatIdList.begin();
//            std::advance(itp, pListId);
//            int imgId = itp->first;
//            Vector2d feat = itp->second;

//            g2o::EdgeProjectXYZ2UV * e
//                = new g2o::EdgeProjectXYZ2UV();
//  //          e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertices().find(thisXId)->second)); //X
//  //          e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertices().find(PidMap[imgId])->second)); //extrinsic
//            //#pragma omp critical
//            {
//  //          e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertices().find(thisXId)->second)); //X
//  //          e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertices().find(imgId)->second)); //extrinsic
//  //              optimizer.vertices()[thisXId];
//            e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertices()[thisXId])); //X
//            e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertices()[imgId])); //extrinsic
//            }
//  //          std::cout << "preparing landmark 5" << std::endl;
//            //
//  //          Vector3d measure;//(feat(0), feat(1), double(imgId/100000-1));
//  //          measure.first = feat;
//  //          measure.second = imgId / 100000 - 1;
//  //          Vector2d img_p = feat;//imgList[imgId].point2dList[featId];
//  //          std::cout << "preparing landmark 5 1" << std::endl;
//  //          e->setMeasurement(img_p);
//  //          e->setMeasurement(img_p);
//            e->setMeasurement(feat);
//  //          std::cout << "preparing landmark 5 2" << std::endl;
//            e->information() = Matrix2d::Identity();
//  //          std::cout << "preparing landmark 5 3" << std::endl;
//  //          if (ROBUST_KERNEL)
//  //          {
//  //            g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
//  //            e->setRobustKernel(rk);
//  //          }
//            e->setParameterId(0, 0);
//  //          std::cout << "preparing landmark 5 4" << std::endl;
//            #pragma omp critical
//            {
//  //              std::cout << "preparing landmark 5 5" << std::endl;
//              optimizer.addEdge(e);
//  //            std::cout << "preparing landmark 6" << std::endl;
//            }
//            ++num1;
//        }
//    }
//  for(std::size_t lListId = 0; lListId < landmarks.size(); ++lListId)
  }


  //test
  {
  std::clock_t startTime = std::clock();
  sleep(1);
  std::clock_t stopTime = std::clock();
  std::cout << "time " << (double) (stopTime-startTime) << std::endl;
  }

  int subSize = 10000;
  int times = landmarks.size()/subSize;
  int subCount;
  int lastCount;
  if(times == 0)
  {
      subCount = landmarks.size();
      lastCount = landmarks.size();

  }else{
      subCount = subSize;
      lastCount = landmarks.size() - (subSize*times);
  }
  std::cout << "subSize : " << subSize << ", times : " << times << ", subCount : " << subCount << ", lastCount : " << lastCount << " ,  landmarks : " << landmarks.size() << std::endl;
//  Landmarks::const_iterator itLThis = landmarks.begin();;
  Landmarks::const_iterator itLLast = landmarks.begin();;
  double spendTimeAll = 0;
  for(std::size_t timeId = 0; timeId < times+1; ++timeId)
  {
//      std::clock_t startTime = std::clock();
      std::size_t endCount;
      if(timeId == times)
      {
         endCount = lastCount;
      }else{
           endCount = subCount;
      }
      #pragma omp parallel for
      for(std::size_t lListId = 0; lListId < endCount; ++lListId)
      {
          Landmarks::const_iterator itLThis = itLLast;
          std::advance(itLThis, lListId);
          //
          g2o::VertexSBAPointXYZ * v_p = new g2o::VertexSBAPointXYZ();
          int thisXId = point_idstart + timeId*subCount+lListId;
          v_p->setId(thisXId);
          v_p->setMarginalized(true);
          v_p->setEstimate(itLThis->second.X);// + Vector3d(Sample::gaussian(1), Sample::gaussian(1), Sample::gaussian(1)));
          //
          #pragma omp critical
          {
              optimizer.addVertex(v_p);
          }
          #pragma omp parallel for
          for(std::size_t pListId = 0; pListId < itLThis->second.imgfeatIdList.size(); ++pListId)
          {
              std::vector<std::pair<int, Vector2d> >::const_iterator itp = itLThis->second.imgfeatIdList.begin();
              std::advance(itp, pListId);
              int imgId = itp->first;
              Vector2d feat = itp->second;
              //
              g2o::EdgeProjectXYZ2UV * e
                  = new g2o::EdgeProjectXYZ2UV();
              #pragma omp critical
              {
    //          e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertices().find(thisXId)->second)); //X
    //          e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertices().find(imgId)->second)); //extrinsic
    //              optimizer.vertices()[thisXId];
                e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertices()[thisXId])); //X
                e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertices()[imgId])); //extrinsic
              }
              e->setMeasurement(feat);
              e->information() = Matrix2d::Identity();
    //          if (ROBUST_KERNEL)
    //          {
    //            g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
    //            e->setRobustKernel(rk);
    //          }
              e->setParameterId(0, 0);
//              std::cout <<  "error 14 " << std::endl;
    //          std::cout << "preparing landmark 5 4" << std::endl;
              #pragma omp critical
              {
                optimizer.addEdge(e);
              }
          }
      }
      ///!!!
      std::advance(itLLast, endCount);


//      std::clock_t stopTime = std::clock();
//      double spendTime1 = (double)(stopTime - startTime)/10000000.0;
//      spendTimeAll += spendTime1;
//      std::cout << "time " << timeId << " : " << spendTime1 << " s, spendTimeAll : " << spendTimeAll << " s, " << spendTimeAll/(double)(timeId+1)*(double)times/3600.0 << std::endl;

  }


//  return 0;

//  std::clock_t stopTime1 = std::clock();
//  double spendTime1 = (double)(stopTime1 - startTime)/1000.0;
//  std::cout << "time2 : " << spendTime1 << " s" << std::endl;
//  std::cout << "num1 : " << num1 << std::endl;

//  std::cout << "Lid : " << lListId << " pSize : " << itL->second.imgfeatIdList.size() << " time : " << spendTime << " s" << std::endl;

//#include "g2o/core/hyper_graph.h"
//  g2o::SparseOptimizer optimizer2;
//  std::set<g2o::Edge*> e = optimizer.edges();


  //
  /// compute
  ///
  std::cout << "start compution" << std::endl;
  cout << endl;
  optimizer.initializeOptimization();
  std::cout << "start compution 1" << std::endl;
  optimizer.setVerbose(true);
  std::cout << "start compution 2" << std::endl;
  if (STRUCTURE_ONLY){
      std::cout << "start compution 3" << std::endl;
    g2o::StructureOnlySolver<3> structure_only_ba;
    cout << "Performing structure-only BA:"   << endl;
    g2o::OptimizableGraph::VertexContainer points;
    for (g2o::OptimizableGraph::VertexIDMap::const_iterator it = optimizer.vertices().begin(); it != optimizer.vertices().end(); ++it) {
      g2o::OptimizableGraph::Vertex* v = static_cast<g2o::OptimizableGraph::Vertex*>(it->second);
      if (v->dimension() == 3)
        points.push_back(v);
    }

    structure_only_ba.calc(points, 10);
  }
  //optimizer.save("test.g2o");
  cout << endl;
  cout << "Performing full BA:" << endl;
  optimizer.optimize(133);
//  optimizer._algorithm->solve(19);
//  double lambda = optimizer._algorithm->properties().at(1)->second;

  cout << endl;

//  optimizer.save((root_path+"/res.g2o").c_str());

  //
  /// compute finished
  /// get result
  ///
  std::cout << "get result" << std::endl;
  XList resXList;
  Extrinsics resPoseList;
  #pragma omp parallel for
  for(std::size_t lListId = 0; lListId < landmarks.size(); ++lListId)
  {
//      Landmarks::const_iterator itL = landmarks.begin();
//      std::advance(itL, lListId);
      //
      int thisXId = point_idstart + lListId;
      // get X
      g2o::HyperGraph::VertexIDMap::iterator itX = optimizer.vertices().find(thisXId);
      if (itX == optimizer.vertices().end())
      {
        cerr << "X vertex " << itX->first << " not in graph!" << endl;
        exit(-1);
      }
      g2o::VertexSBAPointXYZ * resX = dynamic_cast< g2o::VertexSBAPointXYZ * > (itX->second);
      #pragma omp critical
      resXList[lListId] = resX->estimate();//Vector3d();

//      for(std::size_t pListId = 0; pListId < itL->second.imgfeatIdList.size(); ++pListId)
//      {
//          std::vector<std::pair<int, Vector2d> >::const_iterator itp = itL->second.imgfeatIdList.begin();
//          std::advance(itp, pListId);
//          int imgId = itp->first;
////          int featId = itp->second;
//          // get pose
//          g2o::HyperGraph::VertexIDMap::iterator itPose = optimizer.vertices().find(PidMap[imgId]);
//          if (itPose == optimizer.vertices().end())
//          {
//            cerr << "Pose vertex " << itPose->first << " not in graph!" << endl;
//            exit(-1);
//          }
//          g2o::VertexSE3Expmap * resPose = dynamic_cast< g2o::VertexSE3Expmap * > (itPose->second);
//          //
//          // output
//          Eigen::Quaterniond res_q;
//          Vector3d res_t;
//          res_q = Eigen::Quaterniond(resPose->estimate()[0], resPose->estimate()[1], resPose->estimate()[2], resPose->estimate()[3]);
//          res_t = Vector3d(resPose->estimate()[4], resPose->estimate()[5], resPose->estimate()[6]);
//          #pragma omp critical
//          resPoseList[imgId] = std::pair<Eigen::Quaterniond, Vector3d>(res_q, res_t);

//      }
  }
  //
  /// save result
  ///
  {
  std::cout << "save result" << std::endl;
  std::ofstream out_res;
  std::string out_res_path = root_path + "/res_for_show_2.txt";
  out_res.open(out_res_path);
  if(!out_res)
  {
      std::cout << "create " << out_res_path << " file failed!" << std::endl;
      return EXIT_FAILURE;
  }
  int subSize = 10000;
  int times = resXList.size()/subSize;
  int subCount;
  int lastCount;
  if(times == 0)
  {
      subCount = resXList.size();
      lastCount = resXList.size();

  }else{
      subCount = subSize;
      lastCount = resXList.size() - (subSize*times);
  }
  std::cout << "subSize : " << subSize << ", times : " << times << ", subCount : " << subCount << ", lastCount : " << lastCount << " ,  landmarks : " << landmarks.size() << std::endl;
//  Landmarks::const_iterator itLThis = landmarks.begin();;
  XList::const_iterator itLLast = resXList.begin();;
//  double spendTimeAll = 0;
  for(std::size_t timeId = 0; timeId < times+1; ++timeId)
  {
      std::size_t endCount;
      if(timeId == times)
      {
         endCount = lastCount;
      }else{
           endCount = subCount;
      }
      #pragma omp parallel for
      for(std::size_t itXlistId = 0; itXlistId < endCount; ++itXlistId)
      {
          XList::const_iterator itX = itLLast;
          std::advance(itX, itXlistId);
          #pragma omp critical
          {
            out_res << itX->second(0) << " " << itX->second(1) << " " << itX->second(2) << std::endl;
          }
      }
      std::advance(itLLast, endCount);
  }

  out_res.close();
  }




//  //  cout << "Point error before optimisation (inliers only): " << sqrt(sum_diff2/point_num) << endl;
//  point_num = 0;
//  sum_diff2 = 0;
//  for (tr1::unordered_map<int,int>::iterator it=pointid_2_trueid.begin();
//       it!=pointid_2_trueid.end(); ++it){
//    g2o::HyperGraph::VertexIDMap::iterator v_it
//        = optimizer.vertices().find(it->first);
//    if (v_it==optimizer.vertices().end()){
//      cerr << "Vertex " << it->first << " not in graph!" << endl;
//      exit(-1);
//    }
//    g2o::VertexSBAPointXYZ * v_p
//        = dynamic_cast< g2o::VertexSBAPointXYZ * > (v_it->second);
//    if (v_p==0){
//      cerr << "Vertex " << it->first << "is not a PointXYZ!" << endl;
//      exit(-1);
//    }
//    Vector3d diff = v_p->estimate()-true_points[it->second];
//    if (inliers.find(it->first)==inliers.end())
//      continue;
//    sum_diff2 += diff.dot(diff);
//    ++point_num;
//  }
//  cout << "Point error after optimisation (inliers only): " << sqrt(sum_diff2/point_num) << endl;
//  cout << endl;
}
