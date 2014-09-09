import numpy as np
import healpy as hp
import sys
import time

def compute_pcl_estimate(data_file,inv_noise_file,beam_file):
  d = hp.read_map(data_file)
  inv_n = hp.read_map(inv_noise_file)
  n  = 1./np.sqrt(inv_n)

  nside = hp.npix2nside(np.shape(d)[0])

  # compute the total power spectrum
  C_l = hp.anafast(d,lmax=2*nside)

  num_samps = int(1000)

  N_l = np.zeros(np.shape(C_l))

  B_l_in = np.loadtxt(beam_file,delimiter=",")

  B_l = B_l_in[:,1]

  #apply beam to the cls
  C_l /= B_l**2

  # Monte Carlo
  mu = np.zeros(np.shape(d))
  for samp in range(num_samps):
    if samp % 100 == 0 :
      print "samples =",samp
    # draw a realisation from noise
    n_i = np.random.normal(mu,n)

    # find the power spectrum of this realisation
    N_l_i = hp.anafast(n_i,lmax=2*nside)

    # accumulate
    N_l += N_l_i

  N_l /= float(num_samps)

  #apply beam to the nls
  N_l /= B_l**2

  # subtract
  S_l = C_l - N_l

  return (C_l,N_l,S_l)

def write_pcl(output_file,C_l,N_l,S_l):
  ell=np.arange(0,np.shape(S_l)[0])
  np.savetxt(output_file,np.asarray([ell,C_l,N_l,S_l]).T,delimiter=",")


if __name__ == "__main__":
  if len(sys.argv) == 5 :
    start_time = time.time()

    data_file = sys.argv[1]
    inv_noise_file =sys.argv[2]
    beam_file = sys.argv[3]
    output_file = sys.argv[4]

    C_l,N_l,S_l = compute_pcl_estimate(data_file,inv_noise_file,beam_file)
    write_pcl(output_file,C_l,N_l,S_l)

    print ""
    print (time.time() - start_time) / 60.0, 'minutes'
  else:
    print "usage: python ",sys.argv[0],"<data> <inv-noise-cov-mat> <beam-file> <output-cl-file>"
    print "example: python",sys.argv[0], "./data.fits ./invNoise.fits ./window_func_temp_ns128.bl ./base.pcl "
