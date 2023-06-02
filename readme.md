# 介绍

本项目是这一论文的公开源码：[论文链接](https://www.mdpi.com/2072-4292/15/1/177)

## 功能说明

### 输入数据

1. OpenMVG框架下的稀疏点云(包括点云和相机内外参数)

2. 去畸变的图片

3. 所需的空间分辨率

### 输出数据

数字正射影像结果，本质上是一个RGB图像。

# 编译和构建

## 依赖的第三方库

需要满足OpenMVG框架的环境依赖

```shell
sudo apt-get install libpng-dev libjpeg-dev libtiff-dev libxxf86vm1 libxxf86vm-dev libxi-dev libxrandr-dev llvm clang
```

## 编译

可以使用cmake直接编译

```shell
mkdir build
cmake ..
make
```

## 运行

用于生成的正射影像的executable是domGenerate，可以根据代码中argument提示来运行。
