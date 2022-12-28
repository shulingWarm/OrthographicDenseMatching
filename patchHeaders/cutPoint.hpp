#ifndef _CUT_POINT_HPP_
#define _CUT_POINT_HPP_
#include"range2D.hpp"
#include"mvgGetPoint.hpp"

//2022-10-17 输入一个点云和一个范围，删除不在这个范围里面的点
//哈希表必须有迭代器的功能，因为需要对点做删除
template<class Landmark,class Hash_Map,class T=double,
         class JudgeRangeFunctor=JudgeInRange2D<T>, //用于判断点是否在范围内的仿函数
         class GetPoint2DFunctor=GetMvgPoint<Landmark,T,Hash_Map::key_type> //用于从点里面获取二维点信息的
         >
class CutPoint
{
public:
    //用于判断点是否在范围内的仿函数
    JudgeRangeFunctor judgeInRange;

    //从一个point里面获取一个二维点
    //2022-10-18 这里也有可能是一个三维点
    GetPoint2DFunctor getPoint2D;

    void operator()(Hash_Map& points,const T* range) const
    {
        //根据范围对点作删除操作
        for(auto iter=points.begin();iter!=points.end();)
        {
            //从当前的数据里面获取点
            T tempPoint[3];
            getPoint2D(*iter,tempPoint,3);
            //判断点是否在范围内
            if(judgeInRange(range,tempPoint))
            {
                //进入下一个点
                ++iter;
            }
            else
            {
                //把点删除
                iter=points.erase(iter);
            }
        }
    }
};

#endif
