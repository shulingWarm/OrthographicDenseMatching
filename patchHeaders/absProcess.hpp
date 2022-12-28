#pragma once

//基本的抽象过程
template<class TInput,//运行的输入数据
         class TOutput,//输出时的数据
         class InputFuncotr,//输入过程的处理
         class ProcessFunctor,//从输入到输出的处理
         class OutputFunctor //用于保存的输出结果
>
class AbstractProcess
{
public:

    InputFuncotr getInput;
    ProcessFunctor process;
    OutputFunctor makeOutput;

    //完整的运行流程
    void operator()()
    {
        //新建一个流程文件对象
        TInput processObj;
        //新建一个输出流程的对象
        TOutput outputObj;
        //获取输入数据
        getInput(processObj);
        //从输入数据到输出数据
        process(processObj,outputObj);
        //处理输出结果
        makeOutput(outputObj);
    }
};
