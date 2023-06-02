#pragma once

//2023-2-23 用网格单元初始化DOM时的相关配置

class GridConfig
{
public:
    //等于0的时候不使用插值
    static int& useInterpolateFlag()
    {
        static int flag_=0;
        return flag_;
    }

    //2022-2-22 是否固定所有的高度，这意味着所有的网格单元直接用原有的高程参与计算
    static int& useCurrentHeightFlag()
    {
        static int useFlag_=0;
        return useFlag_;
    }
};
