#cmake的最小版本
cmake_minimum_required(VERSION 3.2)
#工程名称
project(DOMGenerate)
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/third_party)
#链接opencv,openmvg改过之后需要使用
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
#编译openmvg的读取用的库
add_subdirectory(./third_party/openMVG)
#添加可执行文件
add_executable(DOMGenerate DOMGenerate.cpp)
target_link_libraries(DOMGenerate openMVGWorld ${OpenCV_LIBS})
#用于测试openmvg编译的代码
add_executable(loadCompile loadCompile.cpp)
target_link_libraries(loadCompile openMVGWorld)
