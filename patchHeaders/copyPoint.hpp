#ifndef _COPY_POINT_HPP_
#define _COPY_POINT_HPP_
#include"findMapNewKey.hpp"
#include"addPointToSfm.hpp"

//把一个点列表复制到SfMData数据里面
template<class SfM_Data,//sfm的数据
         class Landmark,//点的数据类型
         class Landmarks, //点集的数据类型
         class AddOnePointFunctor=AddPointToSfm<SfM_Data,Landmark,Landmarks>
>
class CopySfmDataPoint
{
public:

    //向sfm_data数据里面添加一个点
    AddOnePointFunctor addOnePoint;

    //把一堆点复制到另一个sfmData里面
    void operator()(SfM_Data& dstSfmData,
                    const Landmarks& fromPoints)
    {
        //遍历点集里面所有的点
        for(auto& eachPoint : fromPoints)
        {
            //把当前点添加到sfm_data里面
            addOnePoint(dstSfmData,eachPoint.second);
        }
    }
};

#endif
