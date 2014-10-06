#source_file=/arxiv/source_tree/cosmo-codes/mock_des_cls/source/simple_pcl.py
source_file=/arxiv/source_tree/cosmo-codes/spicy/spice_cl.py
data_file=/arxiv/projects/LSS/DES_mocks_from_Marc/files_2014_09_22/data_test/filemap_base_counts.fits
inv_noise_file=/arxiv/projects/LSS/DES_mocks_from_Marc/files_2014_09_22/data_test/filemap_base_noise.fits
beam_file=/arxiv/projects/LSS/DES_mocks_from_Marc/window_funcs/window_func_temp_ns128.bl
output_file=/arxiv/projects/LSS/DES_mocks_from_Marc/files_2014_09_22/data_test/pcl/filemap_base.pcl
num_samps=10
export OMP_NUM_THREADS=4

python ${source_file} ${data_file} ${inv_noise_file} ${beam_file} ${output_file} ${num_samps} False
