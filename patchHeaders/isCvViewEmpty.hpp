#pragma once

//针对DOM程序里面判断图片是否载入过的方法
//传入的是view里面的迭代器
template<class ViewPair,
         bool invFlag=false //如果这一项为true，则需要把返回结果反转
>
class isViewImageLoad
{
public:
    //判断一个相机是否载入过
    bool operator()(ViewPair& dstPair)
    {
        bool tempFlag=!dstPair.second->cvImg_.empty();
        //如果需要反转就取反
        if(invFlag)
        {
            tempFlag=!tempFlag;
        }
        return tempFlag;
    }
};
