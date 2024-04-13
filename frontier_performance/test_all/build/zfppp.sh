
#! /bin/bash

jsrun -n1 -c1 -a1 -g1 -r1 ./ZFP_hdf5 -z -i ../../gs_dataset/SDRBENCH-CESM-ATM-26x1800x3600/GCLDLWP_1_26_1800_3600.f32 -c test.h5 -n 3 26 1800 3600 -r 11 -x serial -t s -m abs -e 0.1
rm test.h5
jsrun -n1 -c1 -a1 -g1 -r1 ./ZFP_hdf5 -z -i ../../gs_dataset/SDRBENCH-CESM-ATM-26x1800x3600/GCLDLWP_1_26_1800_3600.f32 -c test.h5 -n 3 26 1800 3600 -r 11 -x serial -t s -m abs -e 0.01
rm test.h5
jsrun -n1 -c1 -a1 -g1 -r1 ./ZFP_hdf5 -z -i ../../gs_dataset/SDRBENCH-CESM-ATM-26x1800x3600/GCLDLWP_1_26_1800_3600.f32 -c test.h5 -n 3 26 1800 3600 -r 11 -x serial -t s -m abs -e 1
rm test.h5
jsrun -n1 -c1 -a1 -g1 -r1 ./ZFP_hdf5 -z -i ../../gs_dataset/SDRBENCH-CESM-ATM-26x1800x3600/GCLDLWP_1_26_1800_3600.f32 -c test.h5 -n 3 26 1800 3600 -r 11 -x serial -t s -m abs -e 10
