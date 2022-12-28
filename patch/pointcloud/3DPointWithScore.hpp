#pragma once
#include<Eigen/Core>

//带分数的点云的数据类型
class PointWithScore
{
protected:
    using Vec3D=Eigen::Vector3d;

    //点的三维坐标
    Vec3D coord_;
    //点的分数
    double score_;

public:
    double* getCoord()
    {
        return coord_.data();
    }

    double& getScore()
    {
        return score_;
    }
};
