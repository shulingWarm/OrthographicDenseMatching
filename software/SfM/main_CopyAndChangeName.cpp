
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/image/image_io.hpp"

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::image;
using namespace openMVG::sfm;

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress.hpp"

#include <stdlib.h>
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <pwd.h>

#include "openMVG/system/timer.hpp"

bool
get_src_files_name(std::string root_path, std::vector<std::string>& fileNameList)
{
#ifdef _WIN32
    _finddata_t file;
    long lf;
    std::string src = this->srcDirPath + "\\*.*";
    if ((lf = _findfirst(src.c_str(), &file)) == -1)
    {
        std::cout << this->srcDirPath << " not found" << std::endl;
        return false;
    }
    else{
        while (_findnext(lf, &file) == 0)
        {
            if (strcmp(file.name, ".") == 0 || strcmp(file.name, "..") == 0)
                continue;
            fileNameList.push_back(file.name);
        }
    }


    _findclose(lf);
#else// Linux
    DIR *dir;
    struct dirent *ptr;

    if ((dir=opendir(root_path.c_str())) == NULL)
    {
        std::cout << root_path << " not found" << std::endl;
        return false;
    }

    while ((ptr=readdir(dir)) != NULL)
    {
        if((ptr->d_name == ".") || (ptr->d_name == ".."))  //current / parent
            continue;
        else if(ptr->d_type == 8)  //file
            fileNameList.push_back(ptr->d_name);
        else if(ptr->d_type == 10)  //link file
            continue;
        else if(ptr->d_type == 4)  //dir
            continue;
    }
    closedir(dir);

#endif

    return true;

}


struct Sortie
{
    int id;
    std::string sCBDir = "";
    std::string sCDDir = "";
    std::string sCFDir = "";
    std::string poseFile = "";
    int oldIdMin, oldIdMax, newIdMin, uselessPose, placesCount;

};
typedef std::vector<Sortie> SortieList;

int main(int argc, char *argv[]) {

  CmdLine cmd;
  std::string sSfM_Data_Filename = "";
  std::string sConfigFile = "";
  std::string sOutDir = "";
  bool justChangeName = false;

  cmd.add( make_option('i', sSfM_Data_Filename, "sfmdata") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('f', sConfigFile, "config file") );
  cmd.add( make_option('n', justChangeName, "just use changing images name") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr
      << "Export undistorted images related to a sfm_data file.\n"
      << "Usage: " << argv[0] << '\n'
      << "[-i|--sfmdata] filename, the SfM_Data file to convert\n"
      << "[-o|--outdir] path\n"
      << "[-b|--camera back dir] path\n"
      << "[-d|--camera down dir] path\n"
      << "[-f|--camera front dir] path\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }


  //input config file and set data
  SortieList sortieList;
  std::ifstream inConfig;
  inConfig.open(sConfigFile);
  if(inConfig)
  {
      int count;
      inConfig >> count;
      int line = 0;
      while(line < count && !inConfig.eof())
      {
          Sortie s;
          inConfig >> s.id;
          inConfig >> s.sCBDir;
          inConfig >> s.sCDDir;
          inConfig >> s.sCFDir;
          inConfig >> s.oldIdMin, s.oldIdMax, s.newIdMin, s.uselessPose;

          sortieList.push_back(s);

          ++line;
      }


  }else{
      std::cout << "open config file error, please check it!" << std::endl;
      return EXIT_FAILURE;

  }
  inConfig.close();



  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create( sOutDir );

  SfM_Data sfm_data;

  if(!justChangeName)
  {
      //load intrinsics
      if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(INTRINSICS))) {
        std::cerr << std::endl
          << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
        return EXIT_FAILURE;
      }
  }


  //std::vector<std::string> fileNameList;
  //get_src_files_name(sfm_data.s_root_path, fileNameList);
  //copy images
  for(std::size_t id = 0; id < sortieList.size(); ++id)
  {
      int oldIdMin = sortieList[id].oldIdMin;
      int oldIdMax = sortieList[id].oldIdMax;
      int newIdMin = sortieList[id].newIdMin;



      //for each-three image
      for(int i = oldIdMin; i < oldIdMax; ++i)
      {
          char sOld[5], sNew[5];

          switch (sortieList[id].placesCount) {
          case 2:
              sprintf(sOld, "%02d", i);
              break;
          case 3:
              sprintf(sOld, "%03d", i);
              break;
          case 4:
              sprintf(sOld, "%04d", i);
              break;
          case 5:
              sprintf(sOld, "%05d", i);
              break;
          case 6:
              sprintf(sOld, "%06d", i);
              break;
          default:
              std::cout << "get places count error, please check it!" << std::endl;
              return EXIT_FAILURE;
              break;
          }
          sprintf(sNew, "%05d", (i+newIdMin-oldIdMin));

          std::string oldBName = sortieList[id].sCBDir + "/DSC" + std::string(sOld) + ".JPG";
          std::string oldDName = sortieList[id].sCDDir + "/DSC" + std::string(sOld) + ".JPG";
          std::string oldFName = sortieList[id].sCFDir + "/DSC" + std::string(sOld) + ".JPG";
          std::string newBName = sOutDir + "/1DSC" + std::string(sNew) + ".JPG";
          std::string newDName = sOutDir + "/2DSC" + std::string(sNew) + ".JPG";
          std::string newFName = sOutDir + "/3DSC" + std::string(sNew) + ".JPG";

          if(!justChangeName)
          {
              Intrinsics::const_iterator iterIntrinsicB, iterIntrinsicD, iterIntrinsicF;
              iterIntrinsicB = sfm_data.GetIntrinsics().find(0);
              iterIntrinsicD = sfm_data.GetIntrinsics().find(1);
              iterIntrinsicF = sfm_data.GetIntrinsics().find(2);
              const IntrinsicBase * camB, * camD, * camF;
              camB = iterIntrinsicB->second.get();
              camD = iterIntrinsicD->second.get();
              camF = iterIntrinsicF->second.get();

              Image<RGBColor> imageB, image_udB, imageD, image_udD, imageF, image_udF;

              //b
              if (ReadImage( oldBName.c_str(), &imageB))
              {
                openMVG::system::Timer timer1;
                UndistortImage(imageB, camB, image_udB, BLACK);
                std::cout << i << " undis cost : " << timer1.elapsed() << std::endl;

                openMVG::system::Timer timer2;
                WriteImage(newBName.c_str(), image_udB);
                std::cout << i << " write cost : " << timer2.elapsed() << std::endl;
              }else{
                  std::cout << "image B : " << oldBName << " error" << std::endl;
              }
              //d
              if (ReadImage( oldDName.c_str(), &imageD))
              {
                openMVG::system::Timer timer1;
                UndistortImage(imageD, camD, image_udD, BLACK);
                std::cout << i << " undis cost : " << timer1.elapsed() << std::endl;

                openMVG::system::Timer timer2;
                WriteImage(newDName.c_str(), image_udD);
                std::cout << i << " write cost : " << timer2.elapsed() << std::endl;
              }else{
                  std::cout << "image D : " << oldDName << " error" << std::endl;
              }
              //f
              if (ReadImage( oldFName.c_str(), &imageF))
              {
                openMVG::system::Timer timer1;
                UndistortImage(imageF, camF, image_udF, BLACK);
                std::cout << i << " undis cost : " << timer1.elapsed() << std::endl;

                openMVG::system::Timer timer2;
                WriteImage(newFName.c_str(), image_udF);
                std::cout << i << " write cost : " << timer2.elapsed() << std::endl;
              }else{
                  std::cout << "image F : " << oldFName << " error" << std::endl;
              }

          }else{
              stlplus::file_copy(oldBName, newBName);
              stlplus::file_copy(oldDName, newDName);
              stlplus::file_copy(oldFName, newFName);
          }

      }

  }


  //change pose file
  std::vector<std::string> allPoseInfo;//from  file
  for(std::size_t id = 0; id < sortieList.size(); ++id)
  {
      std::ifstream inPoseFile;
      inPoseFile.open(sortieList[id].poseFile);
      if(inPoseFile)
      {
          int count = sortieList[id].oldIdMax - sortieList[id].oldIdMin + 1 + sortieList[id].placesCount;
          int line = 0;
          while(line < count)
          {
              std::string input;
              std::getline(inPoseFile, input);
//              inPoseFile >> input;
              if(line >= sortieList[id].placesCount)
              {
                  allPoseInfo.push_back(input);
              }

              ++line;
          }

      }else{
          std::cout << "open file " << sortieList[id].poseFile << " error, please check it!" << std::endl;
          return EXIT_FAILURE;
      }
      inPoseFile.close();

  }

  std::ofstream outPoseFile;
  std::string newPoseFilePath = sOutDir + "/GPSDiff_Extraction_merge.txt";
  outPoseFile.open(newPoseFilePath);
  if(outPoseFile)
  {
      for (int i = 0; i < allPoseInfo.size(); ++i)
      {
          outPoseFile << allPoseInfo[i] << std::endl;
      }
  }else{
      std::cout << "create file GPSDiff_Extraction_merge.txt failed !" << std::endl;
  }
  outPoseFile.close();

  return( EXIT_SUCCESS );

}
