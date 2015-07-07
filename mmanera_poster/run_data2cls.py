import numpy as np
import healpy as hp
import sys
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def ascii2map(asciiinpath,mask,z1,z2):
    """
    Create a helpix map from ascii file
    """
    # conversion factor from degrees to radians
    deg2rad = np.pi/180.
    n_side = hp.npix2nside(len(mask))

    data = np.zeros(len(mask))

    logger.info("Reading the ascii file "+asciiinpath)
    logger.info("Accumulating objects between "+str(z1)+" and "+str(z2))
    with open(asciiinpath) as f:
        i = 0
        for line in f:
            # get the entries from a line
            ents = line.split()
            ra_val = float(ents[0])
            dec_val = float(ents[1])
            z_val = float(ents[2])

            # is the object in the given z range
            if z_val > z1 and z_val < z2 :
                # convert to theta and phi
                theta = -deg2rad*dec_val + np.pi/2.
                phi = deg2rad*(ra_val - 180.)

                try:
                    pix = hp.ang2pix(n_side,theta,phi)
                    # is the pixel inside the mask?
                    if mask[pix] > 0:
                        data[pix] += 1
                except:
                    log_str = " Unexpected value of theta = " + str(theta) + " and "+str(phi)+ " for object " + str(i)
                    logger.info(log_str)
                    pass

            # increase the object count
            i = i + 1

    # return the data
    return data

def data2Cl(data,mask):
    """
    Compute the power spectrum of the data
    """

def data2Wtheta(data,mask):
    """
    Compute the correlation function of the data
    """

def run_data2cls(asciiinpath, z1, z2, maskpath, clpath, wthetapath):
    """
    A functon for computing w(theta) and C(l) from input ascii file
    """

    # create a healpix maps with nside identical to the input mask
    logger.info("Reading the mask "+maskpath)
    mask = hp.read_map(maskpath)

    # read the ascii file and populate the pixels
    data = ascii2map(asciiinpath,mask, z1, z2)

    cl = data2Cl(data,mask)

if __name__ == "__main__":
    if len(sys.argv) == 7:
        asciiinpath = sys.argv[1]
        z1 = float(sys.argv[2])
        z2 = float(sys.argv[3])
        maskpath = sys.argv[4]
        clpath = sys.argv[5]
        wthetapath = sys.argv[6]
        run_data2cls(asciiinpath, z1, z2, maskpath, clpath, wthetapath)
    else:
        print "usage python run_data2cls.py  asciiinpath z1 z2 maskpath clpath wthetapath"
