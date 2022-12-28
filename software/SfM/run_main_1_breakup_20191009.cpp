// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//#include "openMVG/cameras/Camera_Common.hpp"
//#include "openMVG/cameras/Cameras_Common_command_line_helper.hpp"
//#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"
//#include "openMVG/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"
//#include "openMVG/sfm/pipelines/global/sfm_global_engine_relative_motions.hpp"
//#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
//#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
//#include "openMVG/sfm/sfm_data.hpp"
//#include "openMVG/sfm/sfm_data_io.hpp"
//#include "openMVG/sfm/sfm_report.hpp"
//#include "openMVG/system/timer.hpp"


//#include "third_party/cmdLine/cmdLine.h"
//#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <memory>
#include <string>

//
#include <ceres/ceres.h>
#include "ceres/problem.h"
#include "ceres/solver.h"
#include <ceres/rotation.h>
#include <ceres/types.h>

using Vec3 = Eigen::Vector3d;
using Vec2 = Eigen::Vector2d;
using Mat3 = Eigen::Matrix<double, 3, 3>;
using Mat3X = Eigen::Matrix<double, 3, Eigen::Dynamic>;
using IndexT = uint32_t;

template<typename K, typename V>
using Hash_Map = std::map<K, V, std::less<K>,
  Eigen::aligned_allocator<std::pair<const K, V> > >;

struct Observation
{
  //Observation():id_feat(UndefinedIndexT) {  }
  //Observation(const Vec2 & p, IndexT idFeat): x(p), id_feat(idFeat) {}

  Vec2 x;
  IndexT id_feat;

//  // Serialization
//  template <class Archive>
//  void save( Archive & ar) const;

//  // Serialization
//  template <class Archive>
//  void load( Archive & ar);
};

using Observations = Hash_Map<IndexT, Observation>;

struct Landmark
{
  Vec3 X;
  Observations obs;
};

class Pose3
{
  protected:

    /// Orientation matrix
    Mat3 rotation_;

    /// Center of rotation
    Vec3 center_;

  public:

    /**
    * @brief Constructor
    * @param r Rotation
    * @param c Center
    * @note Default (without args) defines an Identity pose.
    */
    Pose3
    (
      const Mat3& r = std::move(Mat3::Identity()),
      const Vec3& c = std::move(Vec3::Zero())
    )
    : rotation_( r ), center_( c ) {}

    /**
    * @brief Get Rotation matrix
    * @return Rotation matrix
    */
    const Mat3& rotation() const
    {
      return rotation_;
    }

    /**
    * @brief Get Rotation matrix
    * @return Rotation matrix
    */
    Mat3& rotation()
    {
      return rotation_;
    }

    /**
    * @brief Get center of rotation
    * @return center of rotation
    */
    const Vec3& center() const
    {
      return center_;
    }

    /**
    * @brief Get center of rotation
    * @return Center of rotation
    */
    Vec3& center()
    {
      return center_;
    }

    /**
    * @brief Get translation vector
    * @return translation vector
    * @note t = -RC
    */
    inline Vec3 translation() const
    {
      return -( rotation_ * center_ );
    }


    /**
    * @brief Apply pose
    * @param p Point
    * @return transformed point
    */
    inline Mat3X operator () ( const Mat3X& p ) const
    {
      return rotation_ * ( p.colwise() - center_ );
    }


    /**
    * @brief Composition of poses
    * @param P a Pose
    * @return Composition of current pose and parameter pose
    */
    Pose3 operator * ( const Pose3& P ) const
    {
      return Pose3( rotation_ * P.rotation_, P.center_ + P.rotation_.transpose() * center_ );
    }


    /**
    * @brief Get inverse of the pose
    * @return Inverse of the pose
    */
    Pose3 inverse() const
    {
      return Pose3( rotation_.transpose(),  -( rotation_ * center_ ) );
    }


    /**
    * @brief Return the depth (distance) of a point respect to the camera center
    * @param X Input point
    * @return Distance to center
    */
    double depth( const Vec3 &X ) const
    {
      return ( rotation_ * ( X - center_ ) )[2];
    }

//    /**
//    * Serialization out
//    * @param ar Archive
//    */
//    template <class Archive>
//    inline void save( Archive & ar ) const;

//    /**
//    * @brief Serialization in
//    * @param ar Archive
//    */
//    template <class Archive>
//    inline void load( Archive & ar );
};

using Poses = Hash_Map<IndexT, Pose3>;



/// Define a collection of landmarks are indexed by their TrackId
using Landmarks = Hash_Map<IndexT, Landmark>;


void load_matches(const std::__cxx11::string file_path, Landmarks &inputLM);
void calculate_X(const std::vector<std::pair<Pose3, Vec2> > &pose_x_list, Vec3 &X_res );
void calculate_C(const std::vector<std::pair<Vec3, Vec2> > &X_x_list, Pose3 &camera_res );

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


int main(int argc, char **argv)
{
    std::string matches_file_path = "";
    Landmarks landmarks;
    load_matches(matches_file_path, landmarks);


    Poses poses;
    int camera_num;

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

void run_process()
{
    std::vector<int> used_camera_id_list;
    std::vector<int> used_X_id_list;

    int new_X_num, new_C_num;

    while(new_C_num != 0 && used_camera_id_list.size() )//!=  )
    {


    }

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

        Vec2 this_X0_X1_factor( , );
        Xxy_list.push_back();

    }








}

void random_select(const Landmarks &matches_LM,  std::vector<std::pair<std::pair<int,int>, Vec2> > &x_selected, //<X ID, Camera ID>
                   std::vector<int> &XId_selected, std::vector<int> &CId_selected)
{
    //C-5, X-9





}

void get_newX(const Poses &poses, const int camera_num, const Landmarks &matches_LM, std::map<int, Vec3> &newX)
{
    std::vector<bool> record_solved_pose(camera_num, false);
    for(const auto & pose_it : poses)
    {
        const IndexT indexPose = pose_it.first;
        record_solved_pose[indexPose] = true;
    }



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
            newX[itLM.first] = X_res;

        }


    }

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

void get_newC(const Landmarks matches_LM, const int camera_num, Poses &poses)
{
    std::vector<bool> record_solved_pose(camera_num, false);
    for(const auto & pose_it : poses)
    {
        const IndexT indexPose = pose_it.first;
        record_solved_pose[indexPose] = true;
    }

    //
    for (std::size_t pId = 0; pId < camera_num; ++pId)
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

                //std::map<int, Vec3> &newX
                poses.at(pId) = C_res;

            }
        }
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
