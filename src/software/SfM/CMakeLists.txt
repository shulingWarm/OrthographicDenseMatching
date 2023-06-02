include_directories(${OpenCV_INCLUDE_DIRS})

ADD_EXECUTABLE(DOMGenerate main_GlobalSfM_water5.cpp)
TARGET_LINK_LIBRARIES(DOMGenerate
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  stlplus
  ${OpenCV_LIBS}
  )
