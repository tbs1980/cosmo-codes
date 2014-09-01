#!/bin/bash
#PBS -N b4
#PBS -l nodes=1:ppn=4
#PBS -l mem=16gb
#PBS -l walltime=120:0:00
#PBS -j oe
#PBS -V

source /home/sbalan/login-scripts/binpaths.sh
source /home/sbalan/login-scripts/libpaths.sh

OMP_NUM_THREADS=4
OMP_DYNAMIC=FALSE
export OMP_DYNAMIC OMP_NUM_THREADS
ulimit -Ss unlimited

blackpearl_exe=/share/splinter/sbalan/source_tree/blackperl/build/Blackperl/Interface/Blackperl.exe
data=/share/data1/sbalan/Misc/DES_mocks_from_Marc/desmock_test_filetable_0128_b4_data.fits
inv_noise=/share/data1/sbalan/Misc/DES_mocks_from_Marc/desmock_test_filetable_0128_b4_invNoise.fits
output=/share/data1/sbalan/Misc/DES_mocks_from_Marc/cl_est/blackpearl/b4/b4
window_func=/share/data1/sbalan/Misc/DES_mocks_from_Marc/window_funcs/window_func_temp_ns128.bl

${blackpearl_exe} -m ${data} -n ${inv_noise} -d 1 -p 1 -a 1 -s 0 -f 0 -i 2 -x 256 -e 128 -g 0.175 -b 10000 -r 123 -o ${output} -l true -F ${window_func}

