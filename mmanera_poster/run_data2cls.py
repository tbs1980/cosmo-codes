import numpy as np
import healpy as hp
import sys
import logging
from subprocess import call

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

spice_exe='spice'
spice_data='spice_data.fits'
spice_noise='spice_noise.fits'
spice_mask='spice_mask.fits'
spice_dl='spice_dl.dat'
spice_nl='spice_nl.dat'
spice_bl='spice_bl.dat'
spice_crr='spice_crr.dat'

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

def data2ClWt(data,maskpath, clpath, wthetapath):
    """
    Compute the power spectrum of the data
    """
    # get the correct nside
    logger.info("Reading the mask again for power spectrum  ")
    mask = hp.read_map(maskpath)
    nside = (hp.npix2nside(len(mask)))
    logger.info("Setting nside =  "+str(nside))

    apodize = True

    logger.info("Estimating apodiation sigma  ")
    #compute the power spectrum of the mask
    call([spice_exe,'-mapfile',maskpath,'-corfile',spice_crr])

    #load the crr file
    W_t_in = np.loadtxt(spice_crr,skiprows=1)
    theta = W_t_in[:,0]
    W_t = W_t_in[:,2]

    if apodize :
        logger.info("Computing the apodisation sigma")
        ap_sigma = np.max(theta[np.log10(np.abs(W_t))>-2.5])*180./np.pi

        if ap_sigma <= 0.:
            raise RuntimeError("ap_sigma should be >0")

        print "ap sigma = ",ap_sigma

    # remove the correlation file
    call(['rm',spice_crr])

    #compute the powe spectrum of the data
    logger.info("Computing the power spectrum")
    if apodize :
        call([spice_exe,'-mapfile',spice_data,'-maskfile',maskpath,'-clfile',
            clpath,'-corfile',wthetapath,'-nlmax',str(2*nside),'-verbosity','NO',
            '-thetamax', str(ap_sigma+0.1), '-apodizesigma', str(ap_sigma)])
    else:
        call([spice_exe,'-mapfile',spice_data,'-maskfile',spice_mask,'-clfile',
            clpath,'-corfile',wthetapath,'-nlmax',str(2*nside),'-verbosity','NO'])


    # read the power spectrum
    #D_l = np.loadtxt(spice_dl,skiprows=1)[:,1]

    # read the correlation function
    #W_t = np.loadtxt(spice_crr,skiprows=1)

    # delete files
    #call(['rm',spice_dl])
    #call(['rm',spice_crr])

    #return (D_l,W_t)


def run_data2cls(asciiinpath, z1, z2, maskpath, clpath, wthetapath):
    """
    A functon for computing w(theta) and C(l) from input ascii file
    """

    # create a healpix maps with nside identical to the input mask
    logger.info("Reading the mask "+maskpath)
    mask = hp.read_map(maskpath)

    # read the ascii file and populate the pixels
    data = ascii2map(asciiinpath,mask, z1, z2)

    # write map to disc for spice
    hp.write_map(spice_data,data)

    # compute the power spectrum and correlation fucntion
    data2ClWt(data,maskpath, clpath, wthetapath)

    # delete stuff that are not required
    call(['rm',spice_data])

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
