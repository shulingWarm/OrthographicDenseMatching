#include<iostream>
#include"patch/process/domGenerateWithIndicate.hpp"
#include "third_party/cmdLine/cmdLine.h"

int main(int argc, char **argv)
{
    CmdLine cmd;
    std::string inputSfm;//输入的sfm_data.bin
    cmd.add( make_option('i', inputSfm, "input-file") );
    std::string outputDir;//结果保存的位置
    cmd.add( make_option('o', outputDir, "input-file") );
    std::string existGird;//已有的参考网格
    cmd.add( make_option('g', existGird, "existGird") );
    std::string rangeStr;//已有的参考网格 -r "0 150 -150 0"
    cmd.add( make_option('r', rangeStr, "rangeStr") );
    double pixelLength=0.04;//空间分辨率
    cmd.add( make_option('L', pixelLength, "pixelLength") );
    std::string imgRootPath;//图片的根目录，如果觉得sfm_data.bin里的root_path不对再用这个
    cmd.add( make_option('d', imgRootPath, "imgRootPath") );
    std::string backgroundPath;//生成DOM时的背景图
    cmd.add( make_option('b', backgroundPath, "backgroundPath") );
    std::string maskPath="";//生成DOM的时候，有些区域可能不需要生成，因此可以输入一个mask把它禁用掉
    cmd.add( make_option('m', maskPath, "maskPath") );

    try {
      if (argc == 1) throw std::string("Invalid parameter.");
      cmd.process(argc, argv);
    } catch (const std::string& s) {
      std::cerr << s << std::endl;
      return EXIT_FAILURE;
    }

    //用已有的网格生成DOM
    DomGenerateWithIndicate indicateGenerateTool;
    indicateGenerateTool.generateWithIndicate(
                inputSfm,outputDir,rangeStr,pixelLength,imgRootPath,existGird,backgroundPath,
                maskPath);

	return 0;
}
