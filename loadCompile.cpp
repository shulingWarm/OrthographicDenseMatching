#include<iostream>
#include"openMVG/sfm/sfm_data_io.hpp"
#include"openMVG/sfm/sfm_data.hpp"

int main()
{
    openMVG::sfm::SfM_Data sfm;
    openMVG::sfm::Load(sfm,"",openMVG::sfm::ESfM_Data::ALL);
	return 0;
}
