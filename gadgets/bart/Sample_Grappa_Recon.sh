#!/bin/bash

debug=false

if "$debug"; then
	set -x
fi

bart cc -S reference_data cc_mat
# Compress coils to 12 virtual channels
bart extract 4 0 11 cc_mat cc_mat_P

bart fmac -C -s 8 reference_data cc_mat_P reference_data_cc
bart transpose 3 4 reference_data_cc cc_reference_data

bart fmac -C -s 8 input_data cc_mat_P input_data_cc
bart transpose 3 4 input_data_cc cc_input_data

bart ecalib -c0.65 -k6 -r92 -m2 -S cc_reference_data maps
bart rsense -l1 -r0.001 -i50 cc_input_data maps ims

if "$debug";then
	set +x
fi
