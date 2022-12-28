#ifndef _MVG_GET_POINT_HPP_
#define _MVG_GET_POINT_HPP_
#include<unordered_map>

//针对openmvg的形式获取点的三维坐标的方法
//因为针对的是mvg,所以默认是double
//哈希表的第1个键值必须是unsigned
template<class Landmark,class T=double,class Key=unsigned>
class GetMvgPoint
{
public:
    //从一个openmvg里面的点获取坐标序列
    //计算结果会存储在T*里面
    //这里面需要根据长度来获取对应数量的维度
    //只能是1,2,3里面的一个
    //注意这里传进来的dstPoint的pair<unsigned,Landmark>
    void operator()(const std::pair<Key,Landmark>& dstPoint,T* ans,unsigned length=3) const
    {
        //获取点的引用
        auto& point=dstPoint.second.X;
        //依次往里面写入点
        for(unsigned id=0;id<length;++id)
        {
            ans[id]=point[id];
        }
    }
};

#endif
