#!/bin/bash

# List of available default parameters
# recon_matrix_x;
# recon_matrix_y;
# recon_matrix_z;
# FOV_x;
# FOV_y;
# FOV_z;
# acc_factor_PE1;
# acc_factor_PE2;
# reference_lines_PE1;
# reference_lines_PE2;
# input_data;
# reference_data;
# traj_data;

debug=false

if "$debug"; then
	set -x
fi

# Those two lines will be ignored by the BartGadget but allow you to execute
# the script outside of Gadgetron with little changes
# Possible inputs provided by the BartGadget are:
#   - $reference_data
#   - $input_data
#   - $traj_data

input_data="path/to/a/file"
reference_data="path/to/another/file"

bart cc -S $reference_data cc_mat
# Compress coils to 12 virtual channels
bart extract 4 0 11 cc_mat cc_mat_P

bart fmac -C -s 8 $reference_data cc_mat_P reference_data_cc
bart transpose 3 4 reference_data_cc cc_reference_data

bart fmac -C -s 8 $input_data cc_mat_P input_data_cc
bart transpose 3 4 input_data_cc cc_input_data

bart ecalib -c0.7 -k7 -r$reference_lines_PE1 -m4 -S cc_reference_data maps
bart pics -l1 -r0.1 -i150 cc_input_data maps ims

if "$debug";then
	set +x
fi
