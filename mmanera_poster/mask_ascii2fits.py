import numpy as np
import healpy as hp
import sys
import logging


def ascii2fits(file_1,nside_1,file_2,nside_2):
    """
    A function for converting ascii map to fits with a mask.
    @param file_1 input ascii file with useful pixel numbers
    @param nside_1 nside of the input ascii file
    @param file_2 output file path
    @param nside_2 nside of the output map
    """

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    if(nside_2 > nside_1):
        raise RuntimeError("nside of the input should be greater than the output")

    logger.info("Reading data from "+file_1)
    pixs = np.loadtxt(file_1,skiprows=1)

    logger.info("Creating healpix mask with nside "+str(nside_2))
    mask = np.zeros(hp.nside2npix(nside_1))
    for pix in pixs:
        mask[int(pix)] = 1

    ud_mask = hp.ud_grade(mask,nside_2)

    logger.info("Writing output to "+file_2)
    hp.write_map(file_2,ud_mask)

if __name__ == "__main__":
    if len(sys.argv) == 5:
        file_1 = sys.argv[1]
        nside_1 = int(sys.argv[2])
        file_2 = sys.argv[3]
        nside_2 = int(sys.argv[4])
        ascii2fits(file_1,nside_1,file_2,nside_2)
    else:
        print "usage python mask_ascii2fits.py  filepath1 nside1 filepath2 nside_2"
