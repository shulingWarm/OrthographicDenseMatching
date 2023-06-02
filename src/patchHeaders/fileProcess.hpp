#pragma once
#include"absProcess.hpp"
#include<string>

//以文件作为输入和输出的抽象过程
template<class TInput,//运行的输入数据
         class TOutput,//输出时的数据
         class InputFuncotr,//输入过程的处理
         class ProcessFunctor,//从输入到输出的处理
         class OutputFunctor //用于保存的输出结果
>
class FileAbstractProcess : AbstractProcess<TInput,TOutput,InputFuncotr,ProcessFunctor,OutputFunctor>
{
public:
    //设置输入文件的路径
    void setInputPath(const std::string& filePath)
    {
        getInput.filePath_=filepath;
    }

    //设置输出文件的路径
    void setPutputPath(const std::string& filePath)
    {
        makeOutput.filePath_=filePath;
    }
};
