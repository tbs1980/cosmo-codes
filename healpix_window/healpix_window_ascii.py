import sys
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

# a dictionary of pixel window functions
hpix_dict={}
hpix_dict[2]="pixel_window_n0002.fits"
hpix_dict[4]="pixel_window_n0004.fits"
hpix_dict[8]="pixel_window_n0008.fits"
hpix_dict[16]="pixel_window_n0016.fits"
hpix_dict[32]="pixel_window_n0032.fits"
hpix_dict[64]="pixel_window_n0064.fits"
hpix_dict[128]="pixel_window_n0128.fits"
hpix_dict[256]="pixel_window_n0256.fits"
hpix_dict[512]="pixel_window_n0512.fits"
hpix_dict[1024]="pixel_window_n1024.fits"
hpix_dict[2048]="pixel_window_n2048.fits"
hpix_dict[4096]="pixel_window_n4096.fits"
hpix_dict[8192]="pixel_window_n8192.fits"


def make_window_func_file(hpix_data_dir,nside,num_bins=1):
    filename=hpix_data_dir+hpix_dict[nside]
    w=hp.read_cl(filename)
    l=np.arange(2*nside+1)
    wl=np.zeros(2*nside+1)
    wl=w[0][0:2*nside+1]
    outfile=str(hpix_dict[nside])
    outfile=outfile.replace("fits","dat")
    wind_func_profile = np.zeros([2*nside+1,num_bins+1])
    wind_func_profile[:,0] = l
    
    format_str = '%d'  
    for i in range(num_bins):
        wind_func_profile[:,i+1] = wl
        format_str += ',%0.6e'
    
    np.savetxt(outfile,wind_func_profile,delimiter=",",fmt=format_str)




if __name__ == "__main__":
    if len(sys.argv) == 4:
        hpix_data_dir=sys.argv[1]
        nside=int(sys.argv[2])
        num_bins = int(sys.argv[3])

        make_window_func_file(hpix_data_dir,nside,num_bins)
    else:
        print "usage: python ",sys.argv[0],"<healpix-data-dir> <nside> <num_bins>"
        print "example: python",sys.argv[0], "/arxiv/libraries/source/Healpix_3.11/data/ 256 5"
