#pragma once
#include<string>
#include<iostream>

//最通常的处理抛出字符串异常的内容
class ThrowString
{
public:
    void operator()(const std::string& str)
    {
        std::cerr<<str<<std::endl;
        throw -1;
    }
};
