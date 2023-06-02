#pragma once
#include"openMVG/sfm/sfm_data.hpp"
#include<vector>
#include<algorithm>

//对sfm里面所有的点根据某种规则做排序
class SfmStructureSorter
{
public:
    using SfmData=openMVG::sfm::SfM_Data;
    using Landmark=openMVG::sfm::Landmark;
    using Landmarks=openMVG::sfm::Landmarks;
    //点指针的列表
    using PointPtrList=std::vector<Landmark*>;
    //对点做排序的函数
    typedef bool(*CmpFunc)(const Landmark*,const Landmark*);

    //根据obs对点做排序的方法
    //obs越大越好
    static bool cmpByObsNum(const Landmark* point1,const Landmark* point2)
    {
        return point1->obs.size()>point2->obs.size();
    }

    //从点云里面生成点指针列表
    void makePointPtrList(Landmarks& pointsInHash,PointPtrList& dstList)
    {
        //开辟空间
        dstList.reserve(pointsInHash.size());
        //遍历记录
        for(auto& eachPoint : pointsInHash)
        {
            dstList.push_back(&(eachPoint.second));
        }
    }

    //对点做排序
    void sortPoint(SfmData& sfm,CmpFunc cmpPoint)
    {
        //用另一个容器来存储点
        Landmarks savedPoints;
        savedPoints.swap(sfm.structure);
        //遍历这个点，生成点指针
        PointPtrList pointList;
        makePointPtrList(savedPoints,pointList);
        //对点指针做排序
        std::sort(pointList.begin(),pointList.end(),cmpPoint);
        //按照排序的顺序重新把点放到sfm里面
        for(unsigned idPoint=0;idPoint<pointList.size();++idPoint)
        {
            //直接记录点信息
            sfm.structure[idPoint]=*pointList[idPoint];
        }
    }
};
