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


def make_window_func_file(hpix_data_dir,nside):
    filename=hpix_data_dir+hpix_dict[nside]
    w=hp.read_cl(filename)
    l=np.arange(2*nside+1)
    wl=np.zeros(2*nside+1)
    wl=w[1][0:2*nside+1]# 1 implies polar
    outfile=str(hpix_dict[nside])
    outfile=outfile.replace("fits","dat")
    np.savetxt(outfile,np.asarray([l,wl,wl]).T,delimiter=",",fmt="%d,%0.12e,%0.12e") # we need two columns




if __name__ == "__main__":
    if len(sys.argv) == 3:
        hpix_data_dir=sys.argv[1]
        nside=int(sys.argv[2])

        make_window_func_file(hpix_data_dir,nside)
    else:
        print "usage: python ",sys.argv[0],"<healpix-data-dir> <nside>"
        print "example: python",sys.argv[0], "/arxiv/libraries/source/Healpix_3.11/data 256"
