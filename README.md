# LLPS_codes

Zenodo DOI - https://doi.org/10.5281/zenodo.8147618

Instructions for running MATLAB codes-
Sumon Sahu, May 5th 2023


Make a .tiff stack from set of .nd2 or .czi images in FIJI.
Open average_size_dist_v2_092822.m to analyze the .tiff stack.
Provide path and stack name.
The diameter analysis data is saved in variable ‘diameter array’ and file ‘Stack_diameters.mat’
The number of droplets per image is saved in variable ‘num_dropt’ and file ‘Stack_num_drops.mat’


Make a .tiff stack from set of .nd2 or .czi images in FIJI.
Open radial_profile_code_with_comment_12142022.m to analyze the .tiff stack.
Provide path and stack name.
The intensity data is saved in variable ‘profile array’ and file ‘profile_array_Stack_#.mat’
The normalized radial data is saved in variable ‘norm_r_array’ and file ‘norm_r_array_Stack_#.mat’
Notes - For further analysis, save the profile array and norm_r_array matrices in an excel file.	

Open plot_file_for_radial_distribution.m to bin the radial intensity profiles.
Provide path, excel file name and row-columns of interest.
Output normalized bins along radial direction are stored in ‘r_norm’ and corresponding mean and std intensity values are stored in ‘Yn_array_mean_final’ and ‘Yn_array_std_final’. Plot these to see the intensity distribution.

The 'folder' contains sample files to run for test.
