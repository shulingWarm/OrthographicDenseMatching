#pragma once

//基本的二值比较，这里是小于
template<class T=double>
class CompareLess
{
public:
    bool operator()(const T* data1,const T* data2) const
    {
        return data1[0]<data2[0];
    }
};
