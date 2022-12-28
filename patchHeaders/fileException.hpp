#pragma once
#include"throwString.hpp"

//抛出一个文件读取的异常
template<class ThrowFunctor=ThrowString>
class FilePathError
{
public:
    //基本的抛出异常的方式
    ThrowFunctor throwMethod;

    void operator()(const std::string& invalidPath)
    {
        //用基本的方式抛出异常
        throwMethod("invalid file: "+invalidPath);
    }
};
