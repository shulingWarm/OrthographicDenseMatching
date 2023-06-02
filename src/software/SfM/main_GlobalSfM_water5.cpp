// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/cameras/Cameras_Common_command_line_helper.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_report.hpp"
#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include"openMVG/sfm/switch.hpp"

#include <cstdlib>
#include <memory>
#include <string>


#include"patchHeader/domGrid/domIO.hpp"
#include"patchHeader/process/domGenerate.hpp"
#include"patchHeader/grid/initDomByInterpolate.hpp"
#include"patchHeader/grid/gridConfig.hpp"
#include"patchHeader/process/domGenerateInterpolate.hpp"

using namespace openMVG;
using namespace openMVG::sfm;

//从xyz文件获取海平面的高度
double getDepthFromXyzFile(std::string filePath)
{
    //新建输入流
    std::ifstream fileHandle;
    //打开文件
    fileHandle.open(filePath,std::ios::in);
    //读取空格的位置
    char buffer[50];
    fileHandle.getline(buffer,50,' ');
    fileHandle.getline(buffer,50,' ');
    fileHandle.getline(buffer,50,' ');
    //关闭文件
    fileHandle.close();
    //把内容转换为double数值并返回
    return atof(buffer);
}


int main(int argc, char **argv)
{
  using namespace std;
  std::cout << std::endl
//    << "-----------------------------------------------------------\n"
//    << "Global Structure from Motion:\n"
//    << "-----------------------------------------------------------\n"
//    << "Open Source implementation of the paper:\n"
//    << "\"Global Fusion of Relative Motions for "
//    << "Robust, Accurate and Scalable Structure from Motion.\"\n"
//    << "Pierre Moulon, Pascal Monasse and Renaud Marlet. "
//    << " ICCV 2013." << std::endl
//    << "------------------------------------------------------------"
    << std::endl;


  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sMatchesDir;
  std::string sOutDir = "";
  std::string sIntrinsic_refinement_options = "ADJUST_ALL";
  bool b_use_motion_priors = false;
  std::string rtFilePath = "";
  std::string pointRange="";
  std::string sfmRootPath="";

  std::string xyzPath;//来自mfile的xyz文件的路径
  //快速DOM的分辨率是稀疏点云数量的几倍
  double resolutionForDom=10.0f;
  //固定的像素对应的长度
  double fixPixelLength=0.04;
  //快速DOM输出图片的最大分辨率，输出的图片最大就是这个像素，太大就打不开了
  long int maxResolutionForDom=8000*6000;

  cmd.add( make_option('i', sSfM_Data_Filename, "input-file") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  //固定的每个像素对应的长度，如果不写这个的话，那参考分辨率的比值
  cmd.add(make_option('L',fixPixelLength,"pixel-length"));
  //特殊的运行选项，这个选项只负责显示点云的整体范围，不做其它处理
  cmd.add(make_switch('t',"test-range"));
  //添加点云处理范围的选项，添加这个选项的时候只处理特定范围的点云
  cmd.add(make_option('r',pointRange,"range"));
  //sfm里面图片使用的根目录
  cmd.add(make_option('d',sfmRootPath,"image-dir"));
  //是否在最后的结果里面保存网格数据
  cmd.add(make_switch('g',"save dom grid"));
  //使用初始的稀疏点云重新生成一个网格
  int initByInterpolateFlag=0;
  cmd.add(make_option('m',initByInterpolateFlag,"initByInterpolateFlag"));

  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch (const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    <<"[-i|--input-file] Input scene file by OpenMVG format.\n"
    <<"[-o|--outdir]DOM result output dictionary\n"
    <<"[-L|--pixel-length]DOM pixel equivalent length. Unit: meter  Default 0.04\n"
    <<"[-t|--test-range]Print the range of pointcloud without generating DOM.\n"
    <<"[-r|--range]Used for indicating the range of generated DOM.\n"
    <<"\t For example, \"10 20 30 40\" means the range x:10m~20m y:30m~40m\n"
    <<"[-d|--image-dir] Image root path in sfm data file. Default: the sfm data exist root path.\n"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  //2023-2-22 记录是否需要用多分辨率的形式初始化 DOM
  GridConfig::useInterpolateFlag()=initByInterpolateFlag;

  DomGenerateInterpolate generator;
  generator.generateWithIndicate(sSfM_Data_Filename,sOutDir,
                        pointRange,fixPixelLength,sfmRootPath,"");
  return 0;

//  // Load input SfM_Data scene
//  SfM_Data sfm_data;
//  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
//    std::cerr << std::endl
//      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
//    return EXIT_FAILURE;
//  }

//  //写入sfm_data里面的分辨率相关信息
//    sfm_data.domInfo_.resolutionRate_=resolutionForDom;
//    sfm_data.domInfo_.maxResolution_=maxResolutionForDom;
//    //记录固定的像素值的长度
//    sfm_data.domInfo_.pixelLength_=fixPixelLength;
//    //记录点的范围信息
//    sfm_data.cloudRangeInfo_=pointRange;
//    //记录中间过程的输出目录
//    sfm_data.midProceeFolder_=sOutDir;

//    //是否仅仅测试一下range不做别的
//    bool onlyTestRange=cmd.used('t');

//    //图片的根目录
//    if(sfmRootPath.size()>0)
//    {
//        sfm_data.s_root_path=sfmRootPath;
//    }

//    //2022-10-21 如果只使用下视相机就把下视相机看不到的点删除
//#ifdef USE_ONLY_MID_CAM
//    DeleteViews<
//            Landmarks,
//            CanBeDownViewSeen<
//                    Landmarks::value_type,FindViewId<
//                            Observations,JudgeDownViewiD<>
//                    >,
//                    true //表示如果找不到下视相机就删除
//            >
//    > viewDeleter;
//    viewDeleter(sfm_data.structure);
//#endif

//  sfm_data.loadViewImgs();
//  try {
//      sfm_data.denseDomLikeMvs(onlyTestRange);
//  } catch (int errorFlag) {
//      //这里只负责显示拿到的错误信息
//      std::cout<<errorFlag<<std::endl;
//      throw errorFlag;
//  }
//  if(onlyTestRange) return 0;
//  //保存最后的DOM结果
//  sfm_data.saveDomResult(stlplus::create_filespec(sOutDir, "domResult", ".bmp"));
//  //释放图片内存
//  sfm_data.releaseImgs();
//  //根据dom图里面的坐标重构点云
//  sfm_data.getZAsCloud();

//  Save(sfm_data,
//    stlplus::create_filespec(sOutDir, "cloud_and_poses", ".ply"),
//    ESfM_Data(STRUCTURE));

//  //判断是否需要保存网格文件信息
//  if(cmd.used('g'))
//  {
//      DomIO ioTool;
//      ioTool.saveDom(sfm_data.domInfo_,
//                     stlplus::create_filespec(sOutDir, "grid", ".bin"));
//  }

//  return EXIT_FAILURE;
}
