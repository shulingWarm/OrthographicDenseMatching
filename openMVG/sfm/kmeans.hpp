#ifndef _KMEANS_HPP_
#define _KMEANS_HPP_
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
//赵志豪后来自己添加的一个头文件，从github上抄来的一个聚类算法
//为了完成快速dom算法里面像素块贴图的过程中颜色块选取的问题采取的策略

//这个namespace里面的东西完全是从github上的某个聚类算法里面抄的
namespace kmeans {
using namespace std;
class Point
{
private:
    int pointId, clusterId;
    int dimensions;
    vector<double> values;

    vector<double> lineToVec(string &line)
    {
        vector<double> values;
        string tmp = "";

        for (int i = 0; i < (int)line.length(); i++)
        {
            if ((48 <= int(line[i]) && int(line[i])  <= 57) || line[i] == '.' || line[i] == '+' || line[i] == '-' || line[i] == 'e')
            {
                tmp += line[i];
            }
            else if (tmp.length() > 0)
            {

                values.push_back(stod(tmp));
                tmp = "";
            }
        }
        if (tmp.length() > 0)
        {
            values.push_back(stod(tmp));
            tmp = "";
        }

        return values;
    }

public:
    Point(int id, string line)
    {
        pointId = id;
        values = lineToVec(line);
        dimensions = values.size();
        clusterId = 0; // Initially not assigned to any cluster
    }

    //用一个其它位置的点来初始化一个新的点
    Point(int id,const double* const ptData,int ptDim)
    {
        //记录点的标号
        pointId=id;
        //点坐标数据开辟空间
        values.reserve(ptDim);
        //遍历点数据的每个维度
        for(int dimCount=0;dimCount<ptDim;++dimCount)
            values.push_back(ptData[dimCount]);
        //记录点的维度
        dimensions=ptDim;
        //初始化聚类数据为空
        clusterId=0;
    }

    int getDimensions() { return dimensions; }

    int getCluster() { return clusterId; }

    int getID() { return pointId; }

    void setCluster(int val) { clusterId = val; }

    double getVal(int pos) { return values[pos]; }
};

class Cluster
{
private:
    int clusterId;
    vector<double> centroid;
    vector<Point> points;

public:
    Cluster(int clusterId, Point centroid)
    {
        this->clusterId = clusterId;
        for (int i = 0; i < centroid.getDimensions(); i++)
        {
            this->centroid.push_back(centroid.getVal(i));
        }
        this->addPoint(centroid);
    }

    void addPoint(Point p)
    {
        p.setCluster(this->clusterId);
        points.push_back(p);
    }

    bool removePoint(int pointId)
    {
        int size = points.size();

        for (int i = 0; i < size; i++)
        {
            if (points[i].getID() == pointId)
            {
                points.erase(points.begin() + i);
                return true;
            }
        }
        return false;
    }

    //获取中心点的数据
    //把待修改的点保存到传入的点里面
    //传入的点的数据类型必须支持用[]做随机访问
    template<typename PointType>
    void getCenterPoint(PointType& dstPt) const
    {
        //遍历中心点的每个维度
        for(uint dimCount=0;dimCount<centroid.size();++dimCount)
        {
            //记录当前维度对应的坐标
            dstPt[dimCount]=centroid[dimCount];
        }
    }

    void removeAllPoints() { points.clear(); }

    int getId() { return clusterId; }

    Point getPoint(int pos) { return points[pos]; }

    int getSize() { return points.size(); }

    double getCentroidByPos(int pos) { return centroid[pos]; }

    void setCentroidByPos(int pos, double val) { this->centroid[pos] = val; }
};

class KMeans
{
private:
    int K, iters, dimensions, total_points;
    vector<Cluster> clusters;
    string output_dir;

    void clearClusters()
    {
        for (int i = 0; i < K; i++)
        {
            clusters[i].removeAllPoints();
        }
    }

    int getNearestClusterId(Point point)
    {
        double sum = 0.0, min_dist;
        int NearestClusterId;
        if(dimensions==1) {
            min_dist = abs(clusters[0].getCentroidByPos(0) - point.getVal(0));
        }
        else
        {
          for (int i = 0; i < dimensions; i++)
          {
             sum += pow(clusters[0].getCentroidByPos(i) - point.getVal(i), 2.0);
             // sum += abs(clusters[0].getCentroidByPos(i) - point.getVal(i));
          }
          min_dist = sqrt(sum);
        }
        NearestClusterId = clusters[0].getId();

        for (int i = 1; i < K; i++)
        {
            double dist;
            sum = 0.0;

            if(dimensions==1) {
                  dist = abs(clusters[i].getCentroidByPos(0) - point.getVal(0));
            }
            else {
              for (int j = 0; j < dimensions; j++)
              {
                  sum += pow(clusters[i].getCentroidByPos(j) - point.getVal(j), 2.0);
                  // sum += abs(clusters[i].getCentroidByPos(j) - point.getVal(j));
              }

              dist = sqrt(sum);
              // dist = sum;
            }
            if (dist < min_dist)
            {
                min_dist = dist;
                NearestClusterId = clusters[i].getId();
            }
        }

        return NearestClusterId;
    }

public:
    //输出目录是外来代码的旧接口，在这里不使用，留空即可
    KMeans(int K, int iterations, string output_dir="")
    {
        this->K = K;
        this->iters = iterations;
        this->output_dir = output_dir;
    }

    void run(vector<Point> &all_points)
    {
        total_points = all_points.size();
        dimensions = all_points[0].getDimensions();

        // Initializing Clusters
        vector<int> used_pointIds;

        for (int i = 1; i <= K; i++)
        {
            while (true)
            {
                int index = rand() % total_points;

                if (find(used_pointIds.begin(), used_pointIds.end(), index) ==
                    used_pointIds.end())
                {
                    used_pointIds.push_back(index);
                    all_points[index].setCluster(i);
                    Cluster cluster(i, all_points[index]);
                    clusters.push_back(cluster);
                    break;
                }
            }
        }

        int iter = 1;
        while (true)
        {
            bool done = true;

            // Add all points to their nearest cluster
            #pragma omp parallel for reduction(&&: done) num_threads(16)
            for (int i = 0; i < total_points; i++)
            {
                int currentClusterId = all_points[i].getCluster();
                int nearestClusterId = getNearestClusterId(all_points[i]);

                if (currentClusterId != nearestClusterId)
                {
                    all_points[i].setCluster(nearestClusterId);
                    done = false;
                }
            }

            // clear all existing clusters
            clearClusters();

            // reassign points to their new clusters
            for (int i = 0; i < total_points; i++)
            {
                // cluster index is ID-1
                clusters[all_points[i].getCluster() - 1].addPoint(all_points[i]);
            }

            // Recalculating the center of each cluster
            for (int i = 0; i < K; i++)
            {
                int ClusterSize = clusters[i].getSize();

                for (int j = 0; j < dimensions; j++)
                {
                    double sum = 0.0;
                    if (ClusterSize > 0)
                    {
                        #pragma omp parallel for reduction(+: sum) num_threads(16)
                        for (int p = 0; p < ClusterSize; p++)
                        {
                            sum += clusters[i].getPoint(p).getVal(j);
                        }
                        clusters[i].setCentroidByPos(j, sum / ClusterSize);
                    }
                }
            }

            if (done || iter >= iters)
            {
                break;
            }
            iter++;
        }

        //往文件里面记录这是旧的github代码的写法，这里不需要
        return;

        ofstream pointsFile;
        pointsFile.open(output_dir + "/" + to_string(K) + "-points.txt", ios::out);

        for (int i = 0; i < total_points; i++)
        {
            pointsFile << all_points[i].getCluster() << endl;
        }

        pointsFile.close();

        // Write cluster centers to file
        ofstream outfile;
        outfile.open(output_dir + "/" + to_string(K) + "-clusters.txt");
        if (outfile.is_open())
        {
            for (int i = 0; i < K; i++)
            {
                for (int j = 0; j < dimensions; j++)
                {
                    outfile << clusters[i].getCentroidByPos(j) << " "; // Output to file
                }
                outfile << endl;
            }
            outfile.close();
        }
        else
        {
            cout << "Error: Unable to write to clusters.txt";
        }
    }

    //获取聚类结果最多的类标号
    //比如说第1,2,3个类分别聚合到了3,9,4个数据，那么返回1，表示第2个
    uint getMaxClusterId()
    {
        //初始化聚类信息最多的类名
        uint maxCluId=0;
        //初始化聚类最多的个数
        uint maxAccount=0;
        //遍历每个聚类结果
        for(uint clusterId=0;clusterId<clusters.size();++clusterId)
        {
            //获取当前位置的聚类个数
            uint tempAccount=clusters[clusterId].getSize();
            //判断是否为更大的聚类个数
            if(tempAccount>maxAccount)
            {
                //记录更大的聚类个数
                maxAccount=tempAccount;
                //记录更大的聚类个数对应的类名
                maxCluId=clusterId;
            }
        }
        //返回最大聚类个数对应的类名
        return maxCluId;
    }

    //获取聚类最大结果中的点的个数，并且记录最大类的点的坐标
    //点坐标用模板的形式记录，为了方便使用多种函数
    //模板类中的点类型必须支持用方括号访问，维度也需要自己负责
    template<typename PointType>
    int getMaxCluster(PointType& dstPt)
    {
        //获取聚类数量最多的类对应的类标号
        uint maxCluId=getMaxClusterId();
        //把待记录的点传入到最多的类里面，记录它的中心点
        clusters[maxCluId].getCenterPoint<PointType>(dstPt);
        //返回最多数量的聚类结果到底聚合到了多少个数据
        return clusters[maxCluId].getSize();
    }
};
}
#endif
