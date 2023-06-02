#ifndef _RANGE_2D_HPP_
#define _RANGE_2D_HPP_

//判断某个点是否在特定的范围内
template<class T=double>
class JudgeInRange2D
{
public:
    //范围的安排默认是x_min x_max y_min y_max
    bool operator()(const T* rangeData,const T* point) const
    {
        //返回点在范围内
        return point[0]>=rangeData[0] &&
                point[0]<rangeData[1] &&
                point[1]>=rangeData[2] &&
                point[1]<rangeData[3];
    }
};

#endif
