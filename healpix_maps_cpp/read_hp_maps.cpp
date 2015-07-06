#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <iostream>
#include <string>

int main(void)
{
/*
	// define particulars 
	Healpix_Map<double> map1;
	std::string const filename1("/share/data1/manera/4sree/Y1A1_SPT_frac_o.256_t.32768_EQU_z_V2.fits");
	
	// read the map
	read_Healpix_map_from_fits(filename1,map1);
	
	for(int i=0;i<100;++i)
	{
		std::cout<<i<<"\t"<<map1[i]<<std::endl;
	}
*/	
	Healpix_Map<double> map2;
	std::string const filename2("/share/data1/manera/4sree/Y1A1_SPT_maskfootprint_o.256_t.32768_EQU_all_V2new.fits");
	
	// read the map
	read_Healpix_map_from_fits(filename2,map2);
	
	for(int i=0;i<100;++i)
	{
		std::cout<<i<<"\t"<<map2[i]<<std::endl;
	}
	
	return 0;
}
