import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
import sys
sys.path.append("/arxiv/source_tree/cosmo-codes/spicy")

import spice_cl as scl


def write_pcl(output_file,S_l,S_l_error=None):
    ell=np.arange(0,np.shape(S_l)[0])
    if S_l_error==None :
        np.savetxt(output_file,np.asarray([ell,S_l]).T,delimiter=",")
    else:
        np.savetxt(output_file,np.asarray([ell,S_l,S_l_error]).T,delimiter=",")


############## files names and other inputs ####################################
data_dir="/arxiv/projects/LSS/DES_mocks_from_Marc/files_2014_09_22/data_test/"
data_file = data_dir+"filemap_base_counts.fits"
inv_noise_file = data_dir+"filemap_base_noise.fits"
beam_file = "/arxiv/projects/LSS/DES_mocks_from_Marc/window_funcs/window_func_temp_ns128.bl"
output_prefix = "/arxiv/projects/LSS/DES_mocks_from_Marc/files_2014_09_22/data_test/pcl/filemap"
num_samps=100

######################## compute the base power spectrum #######################
C_l_base,N_l_base,S_l_base = scl.compute_pcl_estimate(data_file,inv_noise_file,beam_file,num_samps)
#write the base power spectrum
output_file = output_prefix+"_base.pcl"
write_pcl(output_file,C_l_base)#base has no noise present


######################### compute the power spectrum of models #################


#model = "fix"
models = ["fix","std"]

for model in models:
    S_l_model = []
    for ind in range(101,126):
        data_file = data_dir+"filemap_"+model+"_counts_"+str(ind)+".fits"
        inv_noise_file = data_dir+"filemap_"+model+"_noise_"+str(ind)+".fits"
        print ""
        print "computing cls for ",data_file
        print ""
        C_l_i,N_l_i,S_l_i = scl.compute_pcl_estimate(data_file,inv_noise_file,beam_file,num_samps)
        S_l_model.append(S_l_i)


    #now find the mean and std dvn
    S_l_model = np.vstack(S_l_model)

    S_l_mean = np.zeros(np.shape(S_l_base))
    S_l_sigma = np.zeros(np.shape(S_l_base))

    for l in range(0,np.shape(S_l_base)[0]):
        S_l_mean[l] = np.mean(S_l_model[:,l])
        S_l_sigma[l] = sqrt(np.var(S_l_model[:,l]))

    #write to file
    output_file = output_prefix+"_"+model+".pcl"
    write_pcl(output_file,S_l_mean,S_l_sigma)

    ############################### polots #########################################
    ell=np.arange(0,np.shape(S_l_base)[0])
    plt.clf()
    plt.plot(ell,ell*(ell+1.)*C_l_base,label="base",color='r')# remember base has no noise present
    plt.errorbar(ell,ell*(ell+1.)*S_l_mean,yerr=ell*(ell+1.)*S_l_sigma,label=model,color='gray',fmt='+')
    plt.xlim(2,256)
    plt.xlabel(r"$\ell$")
    plt.ylabel(r"$\ell(\ell+1)C_{\ell}$")
    plt.grid()
    plt.legend(loc=0)
    output_file = output_prefix+"_"+model+"_plots.png"
    plt.savefig(output_file)
