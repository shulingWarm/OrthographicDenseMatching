#pragma once
#include"openMVG/sfm/sfm_data.hpp"

//给一个已有的DOM添加底图的功能
class AddInitBackground
{
protected:
    //dom网格的数据类型
    using DomGrid=openMVG::sfm::DomInfo;
    //图片的数据类型
    using Image=cv::Mat;

    //把一个图片盖到dom的结果上
    void addBackground(DomGrid& grid,Image& img)
    {
        //把图片恢复到和DOM一样的大小
        cv::resize(img,grid.domResult_,cv::Size(grid.domWidth_,grid.domHeight_));
    }

public:
    //把一个底图添加到DOM里面
    void addBackgroundToDom(DomGrid& grid,
                            const std::string& imgPath //要添加的背景图的路径
    )
    {
        //读取图片
        Image tempImg=cv::imread(imgPath);
        //如果读取不到图片就算了
        if(tempImg.empty())
        {
            if(imgPath.size()>0)
                std::cerr<<"the background "<<imgPath<<" cannot be read"<<std::endl;
            return;
        }
        //把图片添加到dom结果里面作为背景
        addBackground(grid,tempImg);
    }
};
