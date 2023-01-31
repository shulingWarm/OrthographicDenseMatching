#include<iostream>


#include"third_party/cmdLine.h"
#include"patch/process/domGenerate.hpp"

int main(int argc, char **argv)
{
  using namespace std;
  std::cout << std::endl;


  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sOutDir = "";
  std::string pointRange="";
  std::string sfmRootPath="";

  //固定的像素对应的长度
  double fixPixelLength=0.04;

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
    return -1;
  }

  DomGenerateByManager generator;
    generator.domGenerate(sSfM_Data_Filename,sOutDir,
                          pointRange,fixPixelLength,sfmRootPath);
    return 0;

}
