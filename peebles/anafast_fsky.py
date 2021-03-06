from subprocess import call
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import sys
import time


spice_data='spice_data.fits'
spice_data_alm='spice_data_alm.fits'
spice_noise='spice_noise.fits'
spice_noise_alm='spice_noise_alm.fits'
spice_mask='spice_mask.fits'
spice_dl='spice_dl.dat'
spice_nl='spice_nl.dat'
spice_bl='spice_bl.dat'
spice_crr='spice_crr.dat'

def get_mask_file(inv_noise_map):
    mask=np.zeros(np.shape(inv_noise_map))
    mask[inv_noise_map>0]=1

    return mask


def compute_pcl_estimate_with_weights(counts_file,weight_file,beam_file,invWeight=False):
    """
    Compute the anafast/fsky power spectrum from counts map and mask

    @param counts_file a healpix map with galaxy counts
    @param mask_file a healpix map with w
    """
    #read the data file
    d = hp.read_map(counts_file)

    #read the weights file
    weight = hp.read_map(weight_file)

    if d.shape != weight.shape :
        raise RuntimeError("data and weights have different dimensions")

    #force all the values < 0 in weights to zero
    weight[weight<0] = 0.

    #compute the useful area
    useful_pixels = len(weight[weight>0])

    fsky = useful_pixels/float(np.shape(weight)[0])

    print "fsky =",fsky

    #now apply the mask
    if invWeight == True:
    	weight[weight>0] = 1./weight[weight>0]
    
    d *= weight

    total_objects = np.sum(d[weight>0])
    print "total objects = ", total_objects


    N_bar = total_objects/float(useful_pixels) #only some pixels are useful
    print "N_bar = ",N_bar

    #make the data over-density
    d = (d-N_bar)/N_bar

    #make the values outside the mask to zero
    d[weight <= 0] = 0.

    shot_noise = 4*np.pi*fsky/total_objects
    print "shot noise = ", shot_noise

    nside=hp.npix2nside(np.shape(d)[0])

    #load the beam file
    B_l_in = np.loadtxt(beam_file,delimiter=",")
    B_l = B_l_in[:,1]

    #compute the powe spectrum of the data
    D_l = hp.anafast(d,lmax=2*nside)
    #apply the fksy correction
    D_l /=fsky

    #apply beam to the cls
    D_l /= B_l**2

    #compute the noise power spectrum
    N_l = np.ones(np.shape(D_l))*shot_noise

    #apply beam to the nls
    N_l /= B_l**2

    # subtract
    S_l = D_l - N_l


    return (D_l,N_l,S_l)

def compute_pcl_estimate(data_file,inv_noise_file,beam_file,num_samps,counts=False):
    #write the data file
    d = hp.read_map(data_file)

    msk = np.ones(np.shape(d))

    #create a mask file from inv_noise
    if inv_noise_file != None :
        inv_n = hp.read_map(inv_noise_file)
        msk = get_mask_file(inv_n)

    #compute the useful area
    useful_pixels = np.shape(msk)[0]

    if inv_noise_file != None :
        useful_pixels = len(msk[msk>0])

    fsky = 1.

    if inv_noise_file != None :
        fsky = useful_pixels/float(np.shape(msk)[0])

    print "fsky =",fsky

    total_objects = np.sum(d)

    if inv_noise_file != None :
        total_objects = np.sum(d[msk>0])

    print "total objects = ", total_objects

    shot_noise = 4*np.pi*fsky/total_objects

    N_bar = total_objects/float(len(msk)) #all pixels are usefule

    if inv_noise_file != None :
        N_bar = total_objects/float(len(msk[msk>0])) #only some pixels are useful

    # are the input maps counts rather than density?
    if counts == True:
        shot_noise*=(N_bar*N_bar)
        print "N_bar = ",N_bar

    print "shot noise = ", shot_noise

    if d.shape != d.shape :
        raise RuntimeError("data and noise have different dimensions")

    nside=hp.npix2nside(np.shape(d)[0])


    #load the beam file
    B_l_in = np.loadtxt(beam_file,delimiter=",")
    B_l = B_l_in[:,1]

    #compute the powe spectrum of the data
    D_l = hp.anafast(d,lmax=2*nside)
    #apply the fksy correction
    D_l /=fsky

    #apply beam to the cls
    D_l /= B_l**2

    #compute the noise power spectrum
    N_l = np.ones(np.shape(D_l))*shot_noise

    #apply beam to the nls
    N_l /= B_l**2

    # subtract
    S_l = D_l - N_l


    return (D_l,N_l,S_l)


def write_pcl(output_file,C_l,N_l,S_l):
  ell=np.arange(0,np.shape(S_l)[0])
  np.savetxt(output_file,np.asarray([ell,C_l,N_l,S_l]).T,delimiter=",")

if __name__ == "__main__":
  if len(sys.argv) == 6 :
    start_time = time.time()

    data_file = sys.argv[1]
    inv_noise_file =sys.argv[2]
    beam_file = sys.argv[3]
    output_file = sys.argv[4]
    num_samps = int(sys.argv[5])

    C_l,N_l,S_l = compute_pcl_estimate(data_file,inv_noise_file,beam_file,num_samps)
    write_pcl(output_file,C_l,N_l,S_l)

    print ""
    print (time.time() - start_time) / 60.0, 'minutes'
  else:
    print "usage: python ",sys.argv[0],"<data> <inv-noise-cov-mat> <beam-file> <output-cl-file> <n-samps>"
    print "example: python",sys.argv[0], "./data.fits ./invNoise.fits ./window_func_temp_ns128.bl ./base.pcl 100"
