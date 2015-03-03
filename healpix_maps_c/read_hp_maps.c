#include <stdlib.h>
#include <stdio.h>
#include <chealpix.h>

int main(void)
{
	/* define the particulars */
	long const nside = 256;
	long const npix = nside2npix(nside);
	/*char* infile = "/share/data1/manera/4sree/Y1A1_SPT_frac_o.256_t.32768_EQU_z_V2.fits";*/
	char* infile = "/share/data1/manera/4sree/Y1A1_SPT_maskfootprint_o.256_t.32768_EQU_all_V2new.fits";
	char coordsys[10];
	char ordering[10];
	float* hpmap = (float *)malloc(npix*sizeof(float));
	
	/* read the map as float array */
	hpmap = read_healpix_map(infile,&nside,coordsys,ordering);
	
	for(int i=0;i<100;++i)
	{
		printf("%d\t%e\n",i,hpmap[i]);
	}
	
	/* free memory */
	free(hpmap);
	
	return 0;
}
