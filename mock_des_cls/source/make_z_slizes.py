# we assume that the data is formatted as
# ra dec z_obs  Mr (at z=0.1)  g-r (at z=0.1),  g_pbs, r_obs, i_obs ,z_obs, Y_obs

# 122.339938 2.714809 1.49938 -20.0583 0.4302 23.9241 24.0125 24.0119 24.0060
# 122.340257 2.714721 1.49986 -16.7435 0.2219 27.3779 27.3258 27.2076 27.1580
# 122.344290 2.671008 1.49898 -17.3039 0.8330 26.0372 26.0523 25.9714 25.9543
# 122.329476 2.763764 1.49904 -18.8016 0.7857 23.8054 23.7478 23.6235 23.5723
# 122.350149 2.788570 1.49867 -19.5123 0.7961 22.0180 22.3088 22.5070 22.5401
# 122.349923 2.788709 1.49893 -17.3942 0.6991 26.0867 26.2018 26.2164 26.2206
# 122.350995 2.783779 1.49850 -16.9972 0.4215 27.0413 26.7375 26.3845 26.2715
# 122.361283 2.802072 1.49834 -20.8577 0.7985 22.8971 22.7835 22.6266 22.6059
# 122.361351 2.803728 1.49767 -19.9191 0.9377 23.0343 22.8713 22.6138 22.5634
# 122.360214 2.801781 1.49903 -19.4785 0.6393 23.6923 23.4421 23.1104 23.0319
# 122.333567 2.815541 1.49891 -18.8219 0.7533 24.3696 24.4226 24.3584 24.3550
# 122.333720 2.815412 1.49957 -17.6552 0.6007 25.6228 25.7226 25.7138 25.7161
# 122.332781 2.674860 1.49919 -19.4307 0.4451 23.7421 23.4469 23.0734 22.9824


import numpy as np
import healpy as hp
import sys
import time

degr2rad=np.pi/180.

def make_counts_maps(ra,dec,z,z_tics,mask):
    """
    Make counts maps from ra,dec and z data.

    @param ra Right assention of the object in degrees
    @param dec Declination of the object in degrees
    @param z Redshift of the object
    @param z_tics an array of points that define the z-slizes
    @param mask The mask to be applied. This should be a vecor or zero or one
           entries and it should have the same number of pxiels as the output map.

    Returns counts map as a numpy array of dimensions [npix,len(z_tics)-1]
    """

    print "\ncreating the counts map"

    # sanity checks
    if len(ra) != len(dec) : raise RuntimeError("ra and dec have different sizes")
    if len(ra) != len(z) : raise RuntimeError("ra and z have different sizes")
    if len(z_tics) <2 : raise RuntimeError("at least 2 tics required to define a bin")

    nside=0
    try:
        nside = hp.npix2nside(np.shape(mask)[0])
    except:
        raise RuntimeError("length of mask not a healpix npix")

    # define the maps
    nmaps = len(z_tics)-1
    npix = np.shape(mask)[0]
    counts = np.zeros([npix,nmaps])

    # go through each object and add them into the correct bins

    for i in range(len(ra)):
        # convert (ra,dec) to (theta,phi)
        theta=-degr2rad*dec[i] + np.pi/2.
        phi=degr2rad*ra[i]

        # find the pixel id
        pixid=hp.ang2pix(nside,theta,phi)

        # is the pixel in the mask?
        if mask[pixid]>0 :
            # find the bin to which the current object belong to
            binid=np.searchsorted(z_tics,z[i])-1

            # is the bind in the range?
            if binid >=0 and binid < len(z_tics):
                counts[pixid,binid] += 1

    # return the maps
    return counts

def get_rad_dec_z_from_cat(cat_file_name):
    """
    Extract rad,dec,z information from a catalogue file

    Parameters

    @param cat_file_name The name of the catalogue file

    Returns

    1) ra
    2) dec
    3) z

    """

    print "\ncreating the ra,de,z objects\n"

    ra=[]
    dec=[]
    z=[]

    # go through each line of the file and append the ra,dec and z
    with open(cat_file_name) as f:
        for line in f:
            strarr = np.array(line.split(), dtype=float)
            ra.append(strarr[0])
            dec.append(strarr[1])
            z.append(strarr[2])

    print "total number of objects = ",len(z)

    return (ra,dec,z)

def write_counts_maps(counts,z_tics,output_prefix):
    for i in range(np.shape(counts)[1]):
        file_name = output_prefix+"_z_"+str(z_tics[i])+"_"+str(z_tics[i+1])+".fits"
        print file_name
        hp.write_map(file_name,counts[:,i])

def main():
    if len(sys.argv) >= 5 :
        start_time = time.time()

        cat_file_name = sys.argv[1]
        nside = int(sys.argv[2])
        output_prefix = sys.argv[3]

        z_tics = []

        for i in range(4,len(sys.argv)):
            if float(sys.argv[i]) >= 0:
                z_tics.append(float(sys.argv[i]))
            else:
                raise RuntimeError("z_tics contain negative values")

        if list(reversed(sorted(z_tics)))==z_tics :
            raise RuntimeError("z_tics are not in ascending order")

        print "\ncataloge file :",cat_file_name
        print "nside :",nside
        print "output_prefix :",output_prefix
        print "z_tics :",z_tics,"\n"


        ra,dec,z = get_rad_dec_z_from_cat(cat_file_name)

        mask = np.ones(hp.nside2npix(nside))

        counts = make_counts_maps(ra,dec,z,z_tics,mask)

        write_counts_maps(counts,z_tics,output_prefix)

        print (time.time() - start_time) / 60.0, 'minutes'
    else:
        print "usage: python ",sys.argv[0],"<path-to-catalogue> <nside> <output-prefix> <redshift-tics>"
        print "example: python",sys.argv[0], "./marc_catalogue_p00.all ./des_mock 512 0. 2."


if __name__ == "__main__":
    main()
