#include<iostream>
#include"patch/initDOM/transferDom2DomByCloud.hpp"
#include<vector>

using DomGrid=std::vector<double>;
using PointCloud=DomGrid;

//把DOM转换成点云的接口实现
class TransDom2CloudImple : virtual public TransDom2Cloud<DomGrid,PointCloud>
{
public:
    //把一个DOM转换成点云的形式
    virtual void transDom2Cloud(DomGrid& srcGrid,PointCloud& dstCloud) override
    {
        std::cout<<"1"<<std::endl;
    }
};

//用点云初始化一个DOM的接口实现
class InitDomByCloudImple : virtual public InitDomByCloud<DomGrid,PointCloud>
{
public:
    //用一个点云来初始化一个DOM
    virtual void initDomByCloud(DomGrid& grid,const PointCloud& refCloud) override
    {
        std::cout<<2<<std::endl;
    }
};

//通过点云转换DOM的实现
class TransferDomByInterCloud : public MakeDomByTransferCloud<DomGrid,PointCloud>,
        public TransDom2CloudImple,
        public InitDomByCloudImple
{

};


int main()
{
    //新建一个通过点云转换的实现
    TransferDomByInterCloud intrTransfer;
    std::vector<double> vec1,vec2;
    intrTransfer.initDom(vec1,vec2);
	return 0;
}
