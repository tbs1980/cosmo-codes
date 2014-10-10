from subprocess import call
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import sys
import time


spice_exe='spice'
spice_data_1='spice_data_1.fits'
spice_data_2='spice_data_2.fits'
spice_noise_1='spice_noise_1.fits'
spice_noise_2='spice_noise_2.fits'
spice_mask_1='spice_mask_1.fits'
spice_mask_2='spice_mask_2.fits'
spice_dl='spice_dl.dat'
spice_nl='spice_nl.dat'
spice_bl_1='spice_bl_1.dat'
spice_bl_2='spice_bl_2.dat'
spice_crr='spice_crr.dat'
def get_mask_file(inv_noise_map):
    mask=np.zeros(np.shape(inv_noise_map))
    mask[inv_noise_map>0]=1

    return mask

def compute_pcl_estimate(data_file_1,data_file_2,inv_noise_file_1,inv_noise_file_2,beam_file_1,beam_file_2,num_samps,apodize=False):
    #write the data file
    d_1 = hp.read_map(data_file_1)
    d_2 = hp.read_map(data_file_2)
    hp.write_map(spice_data_1,m=d_1)
    hp.write_map(spice_data_2,m=d_2)

    #create a mask file from inv_noise
    inv_n_1 = hp.read_map(inv_noise_file_1)
    inv_n_2 = hp.read_map(inv_noise_file_2)
    msk_1 = get_mask_file(inv_n_1)
    msk_2 = get_mask_file(inv_n_1)
    hp.write_map(spice_mask_1,m=msk_1)
    hp.write_map(spice_mask_2,m=msk_2)

    if d_1.shape != inv_n_1.shape :
        raise RuntimeError("data-1 and noise-1 have different dimensions")
    if d_2.shape != inv_n_2.shape :
        raise RuntimeError("data-2 and noise-2 have different dimensions")
    if d_1.shape != d_2.shape :
        raise RuntimeError("data-1 and data-2 have different dimensions")

    nside=hp.npix2nside(np.shape(d_1)[0])


    #compute the power spectrum of the mask
    call([spice_exe,'-mapfile',spice_mask_1,'-mapfile2',spice_mask_2,'-corfile',spice_crr])#do I need the beam file here?


    #load the crr file
    W_t_in = np.loadtxt(spice_crr,skiprows=1)
    theta = W_t_in[:,0]
    W_t = W_t_in[:,2]

    if apodize :
        ap_sigma = np.max(theta[np.log10(np.abs(W_t))>-3.5])*180./np.pi


        if ap_sigma <= 0.:
            raise RuntimeError("ap_sigma should be >0")

        print "ap sigma = ",ap_sigma


    #write the noise map-1
    n_1 = np.zeros(np.shape(inv_n_1))
    n_1[inv_n_1>0]  = 1./np.sqrt(inv_n_1[inv_n_1>0])

    #write the noise map-2
    n_2 = np.zeros(np.shape(inv_n_2))
    n_2[inv_n_2>0]  = 1./np.sqrt(inv_n_2[inv_n_2>0])

    #write the beam file-1
    B_l_in = np.loadtxt(beam_file_1,delimiter=",")
    np.savetxt(spice_bl_1,np.asarray([B_l_in[:,0],B_l_in[:,1]]).T,fmt='%d   %0.2f')
    #B_l = B_l_in[:,1]

    #write the beam file-2
    B_l_in = np.loadtxt(beam_file_2,delimiter=",")
    np.savetxt(spice_bl_2,np.asarray([B_l_in[:,0],B_l_in[:,1]]).T,fmt='%d   %0.2f')
    #B_l = B_l_in[:,1]

    #compute the powe spectrum of the data
    if apodize :
        call([spice_exe,'-mapfile',spice_data_1,'-mapfile2',spice_data_2,'-maskfile',spice_mask_1,'-maskfile2',spice_mask_2,'-clfile',
            spice_dl,'-corfile','NO','-nlmax',str(2*nside),'-verbosity','NO','-beam_file',spice_bl_1,'-beam_file2',spice_bl_1,
            '-thetamax', str(ap_sigma+0.1), '-apodizesigma', str(ap_sigma)])
    else:
        call([spice_exe,'-mapfile',spice_data_1,'-mapfile2',spice_data_2,'-maskfile',spice_mask_1,'-maskfile2',spice_mask_2,'-clfile',
            spice_dl,'-corfile','NO','-nlmax',str(2*nside),'-verbosity','NO','-beam_file',spice_bl_1,'-beam_file2',spice_bl_1])

    #delete data
    call(['rm',spice_data_1])
    call(['rm',spice_data_2])


    #read the power spectrum
    D_l = np.loadtxt(spice_dl,skiprows=1)[:,1]

    #apply beam to the cls
    #D_l /= B_l**2

    #compute the noise power spectrum using Monte Carlo
    N_l = np.zeros(np.shape(D_l))
    mu = np.zeros(np.shape(d_1))
    sig= np.ones(np.shape(d_1))
    for samp in range(num_samps):
        if samp % 100 == 0 :
            print "samples taken =",samp
        # draw a realisation from noise
        n_1_i = n_1*np.random.normal(mu,sig)
        n_2_i = n_2*np.random.normal(mu,sig)
        #write this to file
        hp.write_map(spice_noise_1,m=n_1_i)
        hp.write_map(spice_noise_2,m=n_2_i)

        # find the power spectrum of this realisation
        if apodize :
            call([spice_exe,'-mapfile',spice_noise_1,'-mapfile2',spice_noise_2,'-maskfile',spice_mask_1,'-maskfile2',spice_mask_2,'-clfile',
                spice_nl,'-corfile','NO','-nlmax',str(2*nside),'-verbosity','NO','-beam_file',spice_bl_1,'-beam_file2',spice_bl_1,
                '-thetamax', str(ap_sigma+0.1), '-apodizesigma', str(ap_sigma)])
        else:
            call([spice_exe,'-mapfile',spice_noise_1,'-mapfile2',spice_noise_2,'-maskfile',spice_mask_1,'-maskfile2',spice_mask_2,'-clfile',
                spice_nl,'-corfile','NO','-nlmax',str(2*nside),'-verbosity','NO','-beam_file',spice_bl_1,'-beam_file2',spice_bl_1])

        #read the power spectrum
        N_l_i = np.loadtxt(spice_nl,skiprows=1)[:,1]

        # accumulate
        N_l += N_l_i

        #delete the noise realisation
        call(['rm',spice_noise_1])
        call(['rm',spice_noise_2])

    N_l /= float(num_samps)

    #apply beam to the nls
    #N_l /= B_l**2

    # subtract
    S_l = D_l - N_l

    #delete the mask
    call(['rm',spice_mask_1])
    call(['rm',spice_mask_2])
    call(['rm',spice_bl_1])
    call(['rm',spice_bl_2])
    call(['rm',spice_crr])
    call(['rm',spice_dl])
    call(['rm',spice_nl])

    return (D_l,N_l,S_l)

def write_pcl(output_file,C_l,N_l,S_l):
  ell=np.arange(0,np.shape(S_l)[0])
  np.savetxt(output_file,np.asarray([ell,C_l,N_l,S_l]).T,delimiter=",")

if __name__ == "__main__":
  if len(sys.argv) == 10 :
    start_time = time.time()

    data_file_1 = sys.argv[1]
    data_file_2 = sys.argv[2]
    inv_noise_file_1 =sys.argv[3]
    inv_noise_file_2 =sys.argv[4]
    beam_file_1 = sys.argv[5]
    beam_file_2 = sys.argv[6]
    output_file = sys.argv[7]
    num_samps = int(sys.argv[8])
    apodize = ( sys.argv[9] == "True" )

    C_l,N_l,S_l = compute_pcl_estimate(data_file_1,data_file_2,inv_noise_file_1,inv_noise_file_2,beam_file_1,beam_file_2,num_samps,apodize)
    write_pcl(output_file,C_l,N_l,S_l)

    print ""
    print (time.time() - start_time) / 60.0, 'minutes'
  else:
    print "usage: python ",sys.argv[0],"<data-1> <data-2> <inv-noise-cov-mat-1> <inv-noise-cov-mat-2> <beam-file-1> <beam-file-2> <output-cl-file> <n-samps> <apodize?>"
    print "example: python",sys.argv[0], "./data_1.fits ./data_2.fits ./invNoise_1.fits ./invNoise_2.fits ./window_func_temp_ns128.bl ./window_func_temp_ns128.bl ./base_cross.pcl 100 False"
