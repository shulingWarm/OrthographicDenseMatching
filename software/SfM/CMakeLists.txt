
###
# Intrinsic image analysis and SfM_Data container initialization
###
ADD_EXECUTABLE(UAVP_main_SfMInit_ImageListing main_SfMInit_ImageListing.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_SfMInit_ImageListing
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  easyexif
  )

ADD_EXECUTABLE(UAVP_main_SfMInit_ImageListing_single main_SfMInit_ImageListing_single.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_SfMInit_ImageListing_single
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  easyexif
  )

#convert a v0.6 lists.txt file to the new sfm_data.X format
ADD_EXECUTABLE(UAVP_main_ConvertList main_ConvertList.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_ConvertList
  openMVG_system
  openMVG_features
  openMVG_sfm
  )

# Installation rules
SET_PROPERTY(TARGET UAVP_main_SfMInit_ImageListing PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_SfMInit_ImageListing DESTINATION bin/)
SET_PROPERTY(TARGET UAVP_main_SfMInit_ImageListing_single PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_SfMInit_ImageListing_single DESTINATION bin/)
SET_PROPERTY(TARGET UAVP_main_ConvertList PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_ConvertList DESTINATION bin/)

###
# Add executable that computes:
# - openMVG_main_ComputeFeatures: features and descriptors
# - openMVG_main_ComputeMatches: putative matches + geometric filtered matches
###

ADD_EXECUTABLE(UAVP_main_ComputeFeatures main_ComputeFeatures.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_ComputeFeatures
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_multiview
  openMVG_sfm
  stlplus
  vlsift
  )

ADD_EXECUTABLE(UAVP_main_ListMatchingPairs main_ListMatchingPairs.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_ListMatchingPairs
  openMVG_system
  openMVG_features
  openMVG_multiview
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_ComputeMatches main_ComputeMatches.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_ComputeMatches
  openMVG_system
  openMVG_features
  openMVG_multiview
  openMVG_sfm
  openMVG_matching_image_collection
  stlplus
  )

###
# Add executable that computes:
# - openMVG_main_ComputeFeatures_gpu: features and descriptors
###
ADD_EXECUTABLE(UAVP_main_ComputeFeatures_gpu main_ComputeFeatures_gpu.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_ComputeFeatures_gpu
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_multiview
  openMVG_sfm
  stlplus
  vlsift
  siftgpu
#  ${SIFTGPU_LIBS}
  IL
#  ${OpenCV_LIBS}
#  ${GLEW_LIBRARIES}
#  ${GLUT_LIBRARIES}
#  ${OPENGL_LIBRARIES}
  )

# Installation rules
SET_PROPERTY(TARGET UAVP_main_ComputeFeatures_gpu PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_ComputeFeatures_gpu DESTINATION bin/)

# Installation rules
SET_PROPERTY(TARGET UAVP_main_ComputeFeatures PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_ComputeFeatures DESTINATION bin/)
SET_PROPERTY(TARGET UAVP_main_ListMatchingPairs PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_ListMatchingPairs DESTINATION bin/)
SET_PROPERTY(TARGET UAVP_main_ComputeMatches PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_ComputeMatches DESTINATION bin/)

###
# SfM Pipelines
# - Incremental/Sequential
# - Global
###
# - Convert SfM_Data format (from one to another)
# - Pair filtering (from frustum)
# - Compute structure from known camera poses
# - Compute structure color
###
ADD_EXECUTABLE(UAVP_main_IncrementalSfM main_IncrementalSfM.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_IncrementalSfM
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

find_package(OpenCV REQUIRED)
ADD_EXECUTABLE(UAVP_main_RunAllStep main_RunAllStep.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_RunAllStep
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_multiview
  openMVG_sfm
  openMVG_matching_image_collection
  stlplus
  vlsift
  ${OpenCV_LIBS}
  )

ADD_EXECUTABLE(UAVP_main_RunGCP main_RunGCP.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_RunGCP
    openMVG_system
    openMVG_image
    openMVG_features
    openMVG_sfm
    stlplus
  )

ADD_EXECUTABLE(UAVP_main_GlobalSfM main_GlobalSfM.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_GlobalSfM_single main_GlobalSfM_single.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM_single
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_GlobalSfM_4 main_GlobalSfM_4.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM_4
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_GlobalSfM_4_margeTest main_GlobalSfM_4_margeTest.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM_4_margeTest
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_GlobalSfM_changeBAStep_save main_GlobalSfM_changeBAStep_save.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM_changeBAStep_save
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_GlobalSfM_changeBAStep_load main_GlobalSfM_changeBAStep_load.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM_changeBAStep_load
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_GlobalSfM_testWaterPlane main_GlobalSfM_testWaterPlane.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM_testWaterPlane
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  ${CERES_LIBRARIES}
  )
target_include_directories(
  UAVP_main_GlobalSfM_testWaterPlane
  PRIVATE
  ${CERES_INCLUDE_DIRS})


ADD_EXECUTABLE(Run_main_1 run_main_1.cpp)
TARGET_LINK_LIBRARIES(Run_main_1
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  ${CERES_LIBRARIES}
  )
target_include_directories(
  Run_main_1
  PRIVATE
  ${CERES_INCLUDE_DIRS})

ADD_EXECUTABLE(main_for_show main_for_show.cpp)
TARGET_LINK_LIBRARIES(main_for_show
    openMVG_system
    openMVG_image
    openMVG_features
    openMVG_sfm
    stlplus
    ${CERES_LIBRARIES}
    )
target_include_directories(
    main_for_show
    PRIVATE
    ${CERES_INCLUDE_DIRS})


ADD_EXECUTABLE(UAVP_main_GlobalSfM_water main_GlobalSfM_water.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM_water
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_GlobalSfM_water2 main_GlobalSfM_water2.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM_water2
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_GlobalSfM_water3 main_GlobalSfM_water3.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM_water3
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_GlobalSfM_water4 main_GlobalSfM_water4.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM_water4
    openMVG_system
    openMVG_image
    openMVG_features
    openMVG_sfm
    stlplus
    ${CERES_LIBRARIES}
    )
target_include_directories(
    UAVP_main_GlobalSfM_water4
    PRIVATE
    ${CERES_INCLUDE_DIRS})

ADD_EXECUTABLE(DOMGenerate main_GlobalSfM_water5.cpp)
TARGET_LINK_LIBRARIES(DOMGenerate
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

#dom数据的保存和读取测试
ADD_EXECUTABLE(domIOTest domIOTest.cpp)
TARGET_LINK_LIBRARIES(domIOTest
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

#用已有的DOM指导新的DOM生成
ADD_EXECUTABLE(generateWithIndicate generateWithIndicate.cpp)
TARGET_LINK_LIBRARIES(generateWithIndicate
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

#把DOM转换成点云的测试代码
ADD_EXECUTABLE(convertDom2CloudTest convertDom2CloudTest.cpp)
TARGET_LINK_LIBRARIES(convertDom2CloudTest
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

#对一个已经初始化过的点云再叠加一些点云
ADD_EXECUTABLE(addCloudToDomTest addCloudToDomTest.cpp)
TARGET_LINK_LIBRARIES(addCloudToDomTest
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_GlobalSfM_water6 main_GlobalSfM_water6.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM_water6
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

#set g2o info
#LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/third_party/third_party/g2o/cmake_modules)
#find_package(G2O REQUIRED)
#if(G2O_FOUND)
#    include_directories(${G2O_INCLOUDE_DIR})
#    message("G2O lib id found:" ${G2O_INCLOUDE_DIR})
#endif(G2O_FOUND)
#LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)
#find_package(Eigen3 REQUIRED)
#find_package(CSparse REQUIRED)
#find_package(Cholmod REQUIRED)
#include_directories(${CSPARSE_INCLUDE_DIR})
#include_directories(${EIGEN3_INCLUDE_DIR})
#include_directories(${CHOLMOD_INCLUDE_DIR})



#SET(G2O_LIBS g2o_cli g2o_ext_freeglut_minimal g2o_simulator g2o_solver_slam2d_linear g2o_types_icp g2o_types_slam2d g2o_core g2o_interface g2o_solver_csparse g2o_solver_structure_only g2o_types_sba g2o_types_slam3d g2o_csparse_extension g2o_opengl_helper g2o_solver_dense g2o_stuff g2o_types_sclam2d g2o_parser g2o_solver_pcg g2o_types_data g2o_types_sim3 cxsparse)

#LINK_DIRECTORIES($(PROJECT_SOURCE_DIR)/third_party/g2o/lib)
#include_directories($(PROJECT_SOURCE_DIR)/third_party/g2o)
#include_directories($(PROJECT_SOURCE_DIR)/third_party)
#include_directories($(PROJECT_SOURCE_DIR)/third_party/g2o/g2o)
#include_directories($(PROJECT_SOURCE_DIR)/third_party/g2o/lib)
#ADD_EXECUTABLE(UAVP_main_GlobalSfM_g2o main_GlobalSfM_g2o.cpp)
#include_directories($(PROJECT_SOURCE_DIR)/third_party/g2o/)
#LINK_DIRECTORIES($(PROJECT_SOURCE_DIR)/third_party/g2o/lib)
#    "/home/guang/software/v201806172215/sparse/src/third_party/g2o/bin")
#FIND_LIBRARY(G2O_USE g2o_use "/home/guang/software/v201806172215/sparse/src/third_party/g2o")
#SET(G2O_FOLDER "/home/guang/software/v201806172215/sparse/src/third_party/g2o")
#link_libraries(g2o_core)
#TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM_g2o
#  openMVG_system
#  openMVG_image
#  openMVG_features
#  openMVG_sfm
#  stlplus
  #${PROJECT_SOURCE_DIR}/third_party/g2o/lib/libg2o.so
  #g2o
  #"/home/guang/software/v201806172215/sparse/src/third_party/g2o/lib/libg2o.a"
  #${G2O_LIBS}
  #${G2O_FOLDER}/lib/libg2o_core.so
  #${G2O_FOLDER}/lib/libg2o_core.a
  #core
  #solver_cholmod
  #solver_csparse
  #solver_pcg
  #solver_eigen
  #types_sba
  #"/home/guang/software/v201806172215/sparse/src/third_party/g2o/lib/libg2o_core.so"
  #"/home/guang/software/v201806172215/sparse/src/third_party/g2o/lib/libg2o_core.a"
  #g2o_core
  #libg2o_core.a
  #../../third_party/g2o/lib/libg2o_core.a
  #solver_cholmod
  #types_sba
#  )
#target_include_directories(UAVP_main_GlobalSfM_g2o
#    PRIVATE
#    ${G2O_INCLUDE_DIRS})
#SET_PROPERTY(TARGET UAVP_main_GlobalSfM_g2o PROPERTY FOLDER OpenMVG/software)
#INSTALL(TARGETS UAVP_main_GlobalSfM_g2o DESTINATION bin/)

ADD_EXECUTABLE(UAVP_main_GlobalSfM_test main_GlobalSfM_test.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM_test
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_GlobalSfM_forSimulation main_GlobalSfM_forSimulation.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM_forSimulation
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  ${CERES_LIBRARIES}
)
target_include_directories(
UAVP_main_GlobalSfM_forSimulation
PRIVATE
${CERES_INCLUDE_DIRS})

ADD_EXECUTABLE(UAVP_main_GlobalSfM_forSimulation_afterOptim main_GlobalSfM_forSimulation_afterOptim.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM_forSimulation_afterOptim
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  ${CERES_LIBRARIES}
)
target_include_directories(
UAVP_main_GlobalSfM_forSimulation_afterOptim
PRIVATE
${CERES_INCLUDE_DIRS})

ADD_EXECUTABLE(UAVP_main_GlobalSfM_AfterMerge main_GlobalSfM_AfterMerge.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM_AfterMerge
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_GlobalSfM_mergeForGrid main_GlobalSfM_mergeForGrid.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM_mergeForGrid
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_GlobalSfM_mergeWithGrid main_GlobalSfM_mergeWithGrid.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM_mergeWithGrid
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_GlobalSfM_sub main_GlobalSfM_sub.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM_sub
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_GlobalSfM_step22 main_GlobalSfM_step22.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM_step22
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_GlobalSfM_Second main_GlobalSfM_Second.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GlobalSfM_Second
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_GetFile main_GetFile.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GetFile
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  ${CERES_LIBRARIES}
  )
target_include_directories(
  UAVP_main_GetFile
  PRIVATE
  ${CERES_INCLUDE_DIRS})

ADD_EXECUTABLE(UAVP_main_PreparingGCPs main_PreparingGCPs.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_PreparingGCPs
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_GetSparse main_GetSparse.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_GetSparse
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_MergeSparse main_MergeSparse.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_MergeSparse
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )


ADD_EXECUTABLE(UAVP_main_ConvertSfM_DataFormat main_ConvertSfM_DataFormat.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_ConvertSfM_DataFormat
  openMVG_system
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_FrustumFiltering main_FrustumFiltering.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_FrustumFiltering
  openMVG_system
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_ComputeStructureFromKnownPoses main_ComputeStructureFromKnownPoses.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_ComputeStructureFromKnownPoses
  openMVG_system
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_ComputeSfM_DataColor main_ComputeSfM_DataColor.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_ComputeSfM_DataColor
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

ADD_EXECUTABLE(UAVP_main_PointsFiltering main_PointsFiltering.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_PointsFiltering
  openMVG_system
  openMVG_sfm
  stlplus
  )

# Installation rules
SET_PROPERTY(TARGET UAVP_main_IncrementalSfM PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_IncrementalSfM DESTINATION bin/)
SET_PROPERTY(TARGET UAVP_main_GlobalSfM PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_GlobalSfM DESTINATION bin/)
SET_PROPERTY(TARGET UAVP_main_GlobalSfM_single PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_GlobalSfM_single DESTINATION bin/)
SET_PROPERTY(TARGET UAVP_main_GlobalSfM_4 PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_GlobalSfM_4 DESTINATION bin/)
SET_PROPERTY(TARGET UAVP_main_GlobalSfM_water PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_GlobalSfM_water DESTINATION bin/)
SET_PROPERTY(TARGET UAVP_main_GetFile PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_GetFile DESTINATION bin/)
SET_PROPERTY(TARGET UAVP_main_GetSparse PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_GetSparse DESTINATION bin/)
SET_PROPERTY(TARGET UAVP_main_ConvertSfM_DataFormat PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_ConvertSfM_DataFormat DESTINATION bin/)
SET_PROPERTY(TARGET UAVP_main_FrustumFiltering PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_FrustumFiltering DESTINATION bin/)
SET_PROPERTY(TARGET UAVP_main_ComputeStructureFromKnownPoses PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_ComputeStructureFromKnownPoses DESTINATION bin/)
SET_PROPERTY(TARGET UAVP_main_ComputeSfM_DataColor PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_ComputeSfM_DataColor DESTINATION bin/)
SET_PROPERTY(TARGET UAVP_main_PointsFiltering PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_PointsFiltering DESTINATION bin/)

###
# SfM tools to visualize feature tracking data
###

# - View extracted Keypoints
#
ADD_EXECUTABLE(UAVP_main_exportKeypoints main_exportKeypoints.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_exportKeypoints
  openMVG_system
  openMVG_features
  openMVG_sfm
  stlplus
  )


# - View computed matches (putatives, geometrics) per image pair
#
ADD_EXECUTABLE(UAVP_main_exportMatches main_exportMatches.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_exportMatches
  openMVG_system
  openMVG_features
  openMVG_sfm
  stlplus
  )

# - View tracks per image pair
#
ADD_EXECUTABLE(UAVP_main_exportTracks main_exportTracks.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_exportTracks
  openMVG_system
  openMVG_features
  openMVG_sfm
  stlplus
  )

# - Export undistorted images related to a sfm_data file
#
ADD_EXECUTABLE(UAVP_main_ExportUndistortedImages main_ExportUndistortedImages.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_ExportUndistortedImages
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

# - Export undistorted images related to a sfm_data file
#
ADD_EXECUTABLE(UAVP_main_CopyAndChangeName main_CopyAndChangeName.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_CopyAndChangeName
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

# installation rules
SET_PROPERTY(TARGET UAVP_main_exportKeypoints PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_exportKeypoints DESTINATION bin/)
SET_PROPERTY(TARGET UAVP_main_exportMatches PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_exportMatches DESTINATION bin/)
SET_PROPERTY(TARGET UAVP_main_exportTracks PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_exportTracks DESTINATION bin/)
SET_PROPERTY(TARGET UAVP_main_ExportUndistortedImages PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_ExportUndistortedImages DESTINATION bin/)
SET_PROPERTY(TARGET UAVP_main_CopyAndChangeName PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_CopyAndChangeName DESTINATION bin/)

###
# SfM export to X
# - PMVS + bundler data format
# - OpenMVS
# - CMPMVS
# - Meshlab
# - MVE (File format v2)
###

# - Export a SfM openMVG scene to PMVS format
#
#ADD_EXECUTABLE(UAVP_main_openMVG2PMVS main_openMVG2PMVS.cpp)
#TARGET_LINK_LIBRARIES(UAVP_main_openMVG2PMVS
#  openMVG_system
#  openMVG_image
#  openMVG_features
#  openMVG_sfm
#  stlplus
#  )

# - Export a SfM openMVG scene to openMVS format
#
ADD_EXECUTABLE(UAVP_main_Transform main_openMVG2openMVS.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_Transform
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  )

find_package(OpenCV REQUIRED)
ADD_EXECUTABLE(UAVP_main_Transform_change main_openMVG2openMVS_change.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_Transform_change
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  ${OpenCV_LIBS}
  )

# - Export a SfM openMVG scene to CMPMVS format
#
#ADD_EXECUTABLE(UAVP_main_openMVG2CMPMVS main_openMVG2CMPMVS.cpp)
#TARGET_LINK_LIBRARIES(UAVP_main_openMVG2CMPMVS
#  openMVG_system
#  openMVG_image
#  openMVG_features
#  openMVG_sfm
#  stlplus
#  )


# - Export a SfM openMVG scene to MVE(v2) format
#
#ADD_EXECUTABLE(UAVP_main_openMVG2MVE2 main_openMVG2MVE2.cpp)
#TARGET_LINK_LIBRARIES(UAVP_main_openMVG2MVE2
#  openMVG_system
#  openMVG_image
#  openMVG_features
#  openMVG_sfm
#  stlplus
#  )

# - Export a SfM openMVG scene to meshlab scene with rasters
# -
#ADD_EXECUTABLE(UAVP_main_openMVG2MESHLAB main_openMVG2MESHLAB.cpp)
#TARGET_LINK_LIBRARIES(UAVP_main_openMVG2MESHLAB
#  openMVG_system
#  openMVG_image
#  openMVG_features
#  openMVG_sfm
#  stlplus
#  )

# - Export a SfM openMVG scene to mvs-texturing scene folder
# -
#ADD_EXECUTABLE(UAVP_main_openMVG2MVSTEXTURING main_openMVG2MVSTEXTURING.cpp)
#TARGET_LINK_LIBRARIES(UAVP_main_openMVG2MVSTEXTURING
#  openMVG_system
#  openMVG_image
#  openMVG_features
#  openMVG_sfm
#  stlplus
#  )

# - Export SfM openMVG camera scene as triangle meshes
# -
#ADD_EXECUTABLE(UAVP_main_ExportCameraFrustums main_ExportCameraFrustums.cpp)
#TARGET_LINK_LIBRARIES(UAVP_main_ExportCameraFrustums
#  openMVG_system
#  openMVG_features
#  openMVG_sfm
#  stlplus
#  )

# - Export SfM openMVG scene to NVM (N-View Match) file
#ADD_EXECUTABLE(UAVP_main_openMVG2NVM main_openMVG2NVM.cpp)
#TARGET_LINK_LIBRARIES( UAVP_main_openMVG2NVM
#  openMVG_system
#  openMVG_image
#  openMVG_features
#  openMVG_sfm
#  stlplus )

# installation rules
#SET_PROPERTY(TARGET UAVP_main_openMVG2PMVS PROPERTY FOLDER OpenMVG/software)
#INSTALL(TARGETS UAVP_main_openMVG2PMVS DESTINATION bin/)
SET_PROPERTY(TARGET UAVP_main_Transform PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_Transform DESTINATION bin/)
#SET_PROPERTY(TARGET UAVP_main_openMVG2CMPMVS PROPERTY FOLDER OpenMVG/software)
#INSTALL(TARGETS UAVP_main_openMVG2CMPMVS DESTINATION bin/)
#SET_PROPERTY(TARGET UAVP_main_openMVG2MESHLAB PROPERTY FOLDER OpenMVG/software)
#INSTALL(TARGETS UAVP_main_openMVG2MESHLAB DESTINATION bin/)
#SET_PROPERTY(TARGET UAVP_main_openMVG2MVSTEXTURING PROPERTY FOLDER OpenMVG/software)
#INSTALL(TARGETS UAVP_main_openMVG2MVSTEXTURING DESTINATION bin/)
#SET_PROPERTY(TARGET UAVP_main_openMVG2MVE2 PROPERTY FOLDER OpenMVG/software)
#INSTALL(TARGETS UAVP_main_openMVG2MVE2 DESTINATION bin/)
#SET_PROPERTY(TARGET UAVP_main_ExportCameraFrustums PROPERTY FOLDER OpenMVG/software)
#INSTALL(TARGETS UAVP_main_ExportCameraFrustums DESTINATION bin/)
#SET_PROPERTY(TARGET UAVP_main_openMVG2NVM PROPERTY FOLDER OpenMVG/software)
#INSTALL(TARGETS UAVP_main_openMVG2NVM DESTINATION bin/)

# -
# Use openCV Feature detector to describe the images
# -
IF(OpenMVG_USE_OPENCV)
  ADD_EXECUTABLE(UAVP_main_ComputeFeatures_OpenCV
    main_ComputeFeatures_OpenCV.cpp)
  TARGET_LINK_LIBRARIES(UAVP_main_ComputeFeatures_OpenCV
    openMVG_system
    openMVG_image
    openMVG_features
    openMVG_multiview
    openMVG_sfm
    stlplus
    vlsift
    ${OpenCV_LIBS})
  target_include_directories(UAVP_main_ComputeFeatures_OpenCV PRIVATE ${OpenCV_INCLUDE_DIRS})
  target_compile_definitions(UAVP_main_ComputeFeatures PRIVATE USE_OCVSIFT)
  
  SET_PROPERTY(TARGET UAVP_main_ComputeFeatures_OpenCV PROPERTY FOLDER OpenMVG/software)
  SET_PROPERTY(TARGET UAVP_main_ComputeFeatures_OpenCV PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
  INSTALL(TARGETS UAVP_main_ComputeFeatures_OpenCV DESTINATION bin/)
ENDIF(OpenMVG_USE_OPENCV)


###
# SfM Research tools
#
# - Quality comparison against a GT camera path (MultiView Evaluation dataset)
###
ADD_EXECUTABLE(UAVP_main_evalQuality main_evalQuality.cpp)
TARGET_LINK_LIBRARIES(UAVP_main_evalQuality
  openMVG_system
  openMVG_features
  openMVG_sfm
  stlplus
  )

#installation rules
SET_PROPERTY(TARGET UAVP_main_evalQuality PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS UAVP_main_evalQuality DESTINATION bin/)

##
# Export tutorial with valid path thanks to a configuration file
##
SET(OPENMVG_SOFTWARE_SFM_SRC_DIR ${CMAKE_CURRENT_LIST_DIR})
IF(MSVC OR APPLE)
  SET(OPENMVG_SOFTWARE_SFM_BUILD_DIR "${EXECUTABLE_OUTPUT_PATH}/${CMAKE_BUILD_TYPE}")
ELSE()
  SET(OPENMVG_SOFTWARE_SFM_BUILD_DIR ${EXECUTABLE_OUTPUT_PATH})
ENDIF()

SET(OPENMVG_TUTORIAL_IN_FILE "${CMAKE_CURRENT_LIST_DIR}/tutorial_demo.py.in")
CONFIGURE_FILE(${OPENMVG_TUTORIAL_IN_FILE}
    "${CMAKE_CURRENT_BINARY_DIR}/tutorial_demo.py" @ONLY)

SET(OPENMVG_TUTORIAL_IN_FILE "${CMAKE_CURRENT_LIST_DIR}/SfM_GlobalPipeline.py.in")
CONFIGURE_FILE(${OPENMVG_TUTORIAL_IN_FILE}
    "${CMAKE_CURRENT_BINARY_DIR}/SfM_GlobalPipeline.py" @ONLY)

SET(OPENMVG_TUTORIAL_IN_FILE "${CMAKE_CURRENT_LIST_DIR}/SfM_SequentialPipeline.py.in")
CONFIGURE_FILE(${OPENMVG_TUTORIAL_IN_FILE}
    "${CMAKE_CURRENT_BINARY_DIR}/SfM_SequentialPipeline.py" @ONLY)

add_subdirectory(clustering)
