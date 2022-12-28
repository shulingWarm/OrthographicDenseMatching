#pragma once
#include"patch/initDOM/mvgInitDomByCloud.hpp"

//把点云里面的数据复制到DOM上的时候，如果分数超过了阈值
//把阈值设置为阈值减小一点的数值，不然会导致DOM无法着色
class DomInitByCloudWithScoreDown : virtual public AddCloudToInitedDom
{
protected:
    //生成DOM时的期望阈值
    const double expScore_=NCC_THRE;
    //超过阈值的情况下写的数值 贴图的时候如果遇到这个数值就表示需要强制贴图
    const double replaceScore_=PRIOR_NCC;

    void replaceCellData(DomCell& cell,Point3D& point) override
    {
        //调用父类的数值记录
        AddCloudToInitedDom::replaceCellData(cell,point);
        //把当前的cell的可用次数改成1次
        cell.referencedTime_=1;
    }

public:

};
