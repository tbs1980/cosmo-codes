# we assume that the data is formatted as
# ra         dec      z_obs  Mr (at z=0.1)  g-r (at z=0.1),  g_pbs, r_obs, i_obs ,z_obs, Y_obs
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

# 0          1         2       3        4      5       6       7       8       9
# ra         dec       z       mag      g-r    g       r       i       z       Y
# 195.072734 42.050844 0.51554 -22.1868 0.9318 21.1540 19.7255 19.1130 18.7745 18.6667
# 195.083890 42.043812 0.51422 -20.1694 0.7102 22.5554 21.4422 21.0598 20.8313 20.7449
# 195.133624 42.081140 0.51284 -17.2598 0.4691 25.2233 24.3356 23.9847 23.7179 23.5936
# 195.087239 42.069510 0.51549 -16.5667 0.7799 26.2685 25.1710 24.7179 24.3803 24.2245
# 195.073770 42.058030 0.51505 -17.4363 0.4904 24.6235 23.8908 23.7240 23.6116 23.5612
# 195.052606 42.059209 0.51979 -16.5719 0.3555 24.9954 24.6090 24.5512 24.4978 24.4706
# 195.077689 42.055170 0.51931 -18.2074 0.5831 24.3963 23.3874 23.0332 22.7782 22.6612
# 195.050116 42.038413 0.51405 -16.1236 0.4008 25.4614 24.9657 24.9647 24.9799 24.9895
# 195.077707 42.062975 0.51655 -16.3672 0.1052 25.0689 24.7595 24.7410 24.7155 24.7008
# 195.071650 42.050773 0.51783 -17.6503 0.2886 24.2317 23.6746 23.5157 23.3889 23.3282



import numpy as np
import healpy as hp
import sys
import time

degr2rad=np.pi/180.

def get_bin_and_pix_ids(ra,dec,z,nside,z_tics):
    theta=-degr2rad*dec + np.pi/2.
    phi=degr2rad*ra

    pix_id=hp.ang2pix(nside,theta,phi)

    bin_id=np.searchsorted(z_tics,z)-1

    return (pix_id,bin_id)


def slice_on_the_fly(cat_file_name,nside,z_tics):
    """
    Extract rad,dec,z information from a catalogue file, apply cuts and then
    add the object to the appropriate z bin and return the counts.

    Parameters

    @param cat_file_name The name of the catalogue file
    @nside the resolution at which we should sample
    @z_tics an array of points that define the z-slizes

    """

    nmaps = len(z_tics)-1
    npix = hp.nside2npix(nside)
    counts = np.zeros([npix,nmaps])

    invalids=[]

    # go through each line of the file and append the ra,dec and z
    with open(cat_file_name) as f:
        for line in f:
            strarr = np.array(line.split(), dtype=float)
            ra = strarr[0]
            dec = strarr[1]
            redshift = strarr[2]
            g = strarr[5]
            r = strarr[6]
            i = strarr[7]
            z = strarr[8]

            # apply cuts
            # 18 < i  < 22.5 and 0 < g - r < 3 and 0 < r - i < 2 and 0 < i - z < 3

            if (i > 18. and i < 22.5
                and (g-r) > 0. and (g-r) < 3.
                and (r-i) > 0. and (r-i) < 2.
                and (i-z) > 0. and (i-z) < 3.):

                #print "\n accepted \n"
                #print i,(g-r),(r-i),(i-z)
                #print "-----------"

                # is the declination value correct
                if dec < -90 or dec > 90:
                    print "\ninvalid parameters for the object, dicarding this one"
                    print strarr
                    print ""
                    invalids.append(strarr)
                else:
                    pix_id,bin_id = get_bin_and_pix_ids(ra,dec,redshift,nside,z_tics)

                    # is the bind in the range?
                    if bin_id >=0 and bin_id < len(z_tics)-1:
                        counts[pix_id,bin_id] += 1
            #else:
                #print "object discared"
                #print i,(g-r),(r-i),(i-z)
                #print ""

    print "total invalis ",len(invalids)
    for invlid in invalids:
        print invalid


    return counts

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


        counts = slice_on_the_fly(cat_file_name,nside,z_tics)

        write_counts_maps(counts,z_tics,output_prefix)

        print (time.time() - start_time) / 60.0, 'minutes'
    else:
        print "usage: python ",sys.argv[0],"<path-to-catalogue> <nside> <output-prefix> <redshift-tics>"
        print "example: python",sys.argv[0], "./marc_catalogue_p00.all 512 ./des_mock 0. 2."


if __name__ == "__main__":
    main()
