#pragma once
#include"openMVG/sfm/sfm_data.hpp"
#include"patch/domGrid/cellInterface.hpp"

//sfm里面的dom
using SfmDomCell=openMVG::sfm::DomUnit;

//把sfm里面的cell套上cell接口的皮
class CellBySfm : public SfmDomCell,
        public CellInterface<double>
{
public:
    double& getCellHeight() override
    {
        return this->z_;
    }

    //获取分数
    double& getCellScore() override
    {
        return this->nccScore_;
    }
};
