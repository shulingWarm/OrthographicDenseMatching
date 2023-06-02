#pragma once
#include<opencv2/opencv.hpp>
#include"openMVG/sfm/sfm_data.hpp"

//给DOM的数据类型添加mask,把mask的地方标注为不可用
class AddMaskToDOM
{
public:
    using DomGrid=openMVG::sfm::DomInfo;

    //mask需要是一个灰度图,不做检查
    //默认mask和grid的大小一样大，不做检查
    void addMask(DomGrid& grid,cv::Mat& mask)
    {
        //遍历mask
        for(unsigned idRow=0;idRow<mask.rows;++idRow)
        {
            //当前行的指针
            auto rowPtr=mask.ptr(grid.domHeight_-idRow-1);
            //遍历每一列
            for(unsigned idCol=0;idCol<mask.cols;++idCol)
            {
                //判断是否超过了阈值
                if(rowPtr[idCol]>10)
                {
                    grid.getUnit(idCol,idRow).setValidFlag(true);
                }
                else
                {
                    grid.getUnit(idCol,idRow).setValidFlag(false);
                }
            }
        }
    }

    //通过路径给DOM添加mask
    void addMask(DomGrid& grid,const std::string& maskPath)
    {
        //如果路径是空的就不处理
        if(maskPath.size()==0) return;
        //读取图片
        cv::Mat mask=cv::imread(maskPath,cv::IMREAD_GRAYSCALE);
        //确保图片具有相同的大小
        cv::Mat resizedMask;
        cv::resize(mask,resizedMask,cv::Size(grid.domWidth_,grid.domHeight_));
        //给网格添加mask
        addMask(grid,resizedMask);
    }
};
