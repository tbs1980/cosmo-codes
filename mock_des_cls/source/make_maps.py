# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import numpy as np
import healpy as hp
import sys
import time


# read the pixel values from Marc's file and create maps
# The format of the file is
# ipix, base err, a1 err, a2 err, a3 err, a4 err, b1 err, b2 err, b3 err, b4 err

NSIDE=128
NPIX=hp.nside2npix(NSIDE)

# make a dictionary for the columns corresponding to maps
headers = {}
headers["base"]=1
headers["a1"]=3
headers["a2"]=5
headers["a3"]=7
headers["a4"]=9
headers["b1"]=11
headers["b2"]=13
headers["b3"]=15
headers["b4"]=17

def read_data(file_name):
	print "Reading the file",file_name
	marc_data = np.loadtxt(file_name,skiprows=1)

	return marc_data

def make_maps(marc_data,map_id):
	"""
	Make maps by reading the file in the fomat
	ipix, base err, a1 err, a2 err, a3 err, a4 err, b1 err, b2 err, b3 err, b4 err

	@param map_id the id (string) of the map {base, a1, a2, a3, a4, b1, b2, b3, b4}
	"""


	print "getting data of ",map_id

	# what data do we need?
	col = headers[map_id]

	# get it
	data = marc_data[:,col]
	noise_sigma = marc_data[:,col+1]

	inverse_noise = np.zeros(NPIX)

	inverse_noise[noise_sigma>0] = 1./noise_sigma[noise_sigma>0]**2

	return (data,inverse_noise)

def write_data_and_noise(file_prefix,data,inverse_noise):
	"""
	Write the data and iverse noise matrix as a fits file

	@param file_prefix prefix to the files
	@param data data as a healpix map
	@param inverse_noise inverse noise as a healpix map
	"""

	data_file_name = file_prefix + "_data.fits"
	inv_noise_file_name = file_prefix + "_invNoise.fits"

	print "writing data to file",data_file_name

	hp.write_map(data_file_name,data)

	print "writing inv noise to file",inv_noise_file_name

	hp.write_map(inv_noise_file_name,inverse_noise)


#############################
if __name__ == "__main__":
	if len(sys.argv) == 2 :
		start_time = time.time()

		marc_file_name = sys.argv[1]
		output_file_prefix = marc_file_name.replace(".dat","_")

		marc_data=read_data(marc_file_name)
		print ""

		for map_id in ["base","a1","a2","a3","a4","b1","b2","b3","b4"]:
			data,inverse_noise=make_maps(marc_data,map_id)
			write_data_and_noise(output_file_prefix+map_id,data,inverse_noise)
			print ""

		print (time.time() - start_time) / 60.0, 'minutes'
	else:
		print "usage: python ",sys.argv[0],"<path-to-catalogue>"
		print "example: python",sys.argv[0], "/share/data1/sbalan/Misc/DES_mocks_from_Marc/desmock_test_filetable_0128.dat"
