#pragma once

//用一个已有的DOM结果初始化一个DOM
template<class DomGrid>
class InitDomByExist
{
public:
    //只是一个虚接口，用一个已有的DOM构造一个新的DOM
    virtual void initDom(DomGrid& srcGrid,DomGrid& dstGrid)=0;
};
