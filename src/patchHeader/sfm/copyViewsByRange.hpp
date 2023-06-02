#pragma once
#include"patchHeader/sfm/copySfmView.hpp"
#include"patchHeader/sfm/getInfoFromSfm.hpp"
#include<memory>

//通过范围决定是否复制的view
//基本上也就这个头文件会用到这个类，因此就暂时把它放这了
//以后要是想复用可以单独把它弄成头文件
class CopyViewByRange : public CopySfmView
{
public:
    const double* range_;
    std::shared_ptr<SfmInterface> sfmInterface_;
    //选取view时的范围扩展 如果卡着范围拿view,会导致这个区域的侧视相机不足
    double rangeExtend_=400;

    //构造函数
    CopyViewByRange(const double* range,
                    std::shared_ptr<SfmInterface> sfmInterface) :
        range_(range), sfmInterface_(sfmInterface){}

    virtual bool judgeNeedCopy(unsigned idView, SfmData &sfm) override
    {
        //需要确保传入的相机的prior是在范围内的
        auto priorCenter=sfmInterface_->getPriorCenterByView(idView,sfm);
        return priorCenter[0]>=range_[0]-rangeExtend_ &&
                priorCenter[0]<range_[1]+rangeExtend_ &&
                priorCenter[1]>=range_[2]-rangeExtend_ &&
                priorCenter[1]<range_[3]+rangeExtend_;
    }
};
