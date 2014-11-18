from subprocess import call
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import sys
import time


map_to_alm='Map2Alm'
alm_to_cl='Alm2Cl'
spice_data='spice_data.fits'
spice_data_alm='spice_data_alm.fits'
spice_noise='spice_noise.fits'
spice_noise_alm='spice_noise_alm.fits'
spice_mask='spice_mask.fits'
spice_dl='spice_dl.dat'
spice_nl='spice_nl.dat'
spice_bl='spice_bl.dat'
spice_crr='spice_crr.dat'
spice_ilm_jlm='spice_ilm_jlm.dat'
def get_mask_file(inv_noise_map):
    mask=np.zeros(np.shape(inv_noise_map))
    mask[inv_noise_map>0]=1

    return mask

def compute_peebles_pcl_estimate(data_file,inv_noise_file,beam_file,num_samps):
    #write the data file
    d = hp.read_map(data_file)
    hp.write_map(spice_data,m=d)

    #create a mask file from inv_noise
    inv_n = hp.read_map(inv_noise_file)
    msk = get_mask_file(inv_n)
    hp.write_map(spice_mask,m=msk)

    if d.shape != inv_n.shape :
        raise RuntimeError("data and noise have different dimensions")

    nside=hp.npix2nside(np.shape(d)[0])


    #write the noise map
    n = np.zeros(np.shape(inv_n))
    n[inv_n>0]  = 1./np.sqrt(inv_n[inv_n>0])

    #write the beam file
    B_l_in = np.loadtxt(beam_file,delimiter=",")
    np.savetxt(spice_bl,np.asarray([B_l_in[:,0],B_l_in[:,1]]).T,fmt='%d   %0.2f')
    B_l = B_l_in[:,1]


    #compute the powe spectrum of the data
    call([map_to_alm,'-I',spice_data,'-O',spice_data_alm,'-L',str(2*nside),'-m',spice_mask])
    call([alm_to_cl,'-I',spice_data_alm,'-O',spice_dl,'-P','-m',spice_mask,'-G','-C',spice_ilm_jlm,'-M',spice_data,'-N',str(nside),'-L',str(2*nside+1)])

    call(['rm',spice_data_alm])
    call(['rm',spice_data])

    return

    #read the power spectrum
    D_l = np.loadtxt(spice_dl,skiprows=2)[:,1]

    #apply beam to the cls
    D_l /= B_l**2

    #compute the noise power spectrum using Monte Carlo
    N_l = np.zeros(np.shape(D_l))

    # subtract
    S_l = D_l - N_l

    #delete the mask
    call(['rm',spice_mask])
    call(['rm',spice_dl])
    call(['rm',spice_bl])

    return (D_l,N_l,S_l)

def compute_pcl_estimate(data_file,inv_noise_file,beam_file,num_samps):
    #write the data file
    d = hp.read_map(data_file)
    hp.write_map(spice_data,m=d)

    #create a mask file from inv_noise
    inv_n = hp.read_map(inv_noise_file)
    msk = get_mask_file(inv_n)
    hp.write_map(spice_mask,m=msk)

    if d.shape != inv_n.shape :
        raise RuntimeError("data and noise have different dimensions")

    nside=hp.npix2nside(np.shape(d)[0])


    #write the noise map
    n = np.zeros(np.shape(inv_n))
    n[inv_n>0]  = 1./np.sqrt(inv_n[inv_n>0])

    #write the beam file
    B_l_in = np.loadtxt(beam_file,delimiter=",")
    np.savetxt(spice_bl,np.asarray([B_l_in[:,0],B_l_in[:,1]]).T,fmt='%d   %0.2f')
    B_l = B_l_in[:,1]


    #compute the powe spectrum of the data
    call([map_to_alm,'-I',spice_data,'-O',spice_data_alm,'-L',str(2*nside),'-m',spice_mask])
    call([alm_to_cl,'-I',spice_data_alm,'-O',spice_dl,'-m',spice_mask,'-M',spice_data,'-N',str(nside),'-L',str(2*nside+1)])

    call(['rm',spice_data_alm])
    call(['rm',spice_data])

    #read the power spectrum
    D_l = np.loadtxt(spice_dl,skiprows=2)[:,1]

    #apply beam to the cls
    D_l /= B_l**2

    #compute the noise power spectrum using Monte Carlo
    N_l = np.zeros(np.shape(D_l))
    mu = np.zeros(np.shape(d))
    sig= np.ones(np.shape(d))
    for samp in range(num_samps):
        if samp % 100 == 0 :
            print "samples taken =",samp
        # draw a realisation from noise
        n_i = n*np.random.normal(mu,sig)
        #write this to file
        hp.write_map(spice_noise,m=n_i)

        # find the power spectrum of this realisation
        call([map_to_alm,'-I',spice_noise,'-O',spice_noise_alm,'-L',str(2*nside),'-m',spice_mask])
        call([alm_to_cl,'-I',spice_noise_alm,'-O',spice_nl,'-m',spice_mask,'-M',spice_noise,'-N',str(nside),'-L',str(2*nside+1)])

        #read the power spectrum
        N_l_i = np.loadtxt(spice_nl,skiprows=2)[:,1]

        # accumulate
        N_l += N_l_i

        #delete the noise realisation
        call(['rm',spice_noise])
        call(['rm',spice_noise_alm])

    N_l /= float(num_samps)

    #apply beam to the nls
    N_l /= B_l**2

    # subtract
    S_l = D_l - N_l

    #delete the mask
    call(['rm',spice_mask])
    call(['rm',spice_nl])
    call(['rm',spice_dl])
    call(['rm',spice_bl])

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
