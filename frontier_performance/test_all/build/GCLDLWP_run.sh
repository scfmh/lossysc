#!/bin/bash
#SBATCH -A CSC303
#SBATCH -J hdf5_single
#SBATCH -o %x-%j.out
#SBATCH -t 0:05:00
#SBATCH -p batch
#SBATCH -N 1

module load rocm/5.4.0
module load cmake
module load cray-mpich/8.1.17
module load cray-python/3.9.12.1
module load zstd

GCLDLWP=/lustre/orion/proj-shared/csc303/zq53/mgard_new/adios_more/gs_dataset/SDRBENCH-CESM-ATM-26x1800x3600/GCLDLWP_1_26_1800_3600.f32
gs=/lustre/orion/proj-shared/csc303/zq53/mgard_new/adios_more/gs_dataset/L512Du0.3Dv0.3F0.01K0.05steps600noise1e-07_4g.dat
sim=/lustre/orion/proj-shared/csc303/zq53/mgard_new/adios_more/hdf5_all/build/

SRUN='srun -n 1 -N 1 --ntasks-per-node=1'
GRUN='srun -N1 -n1 -c1 --gpus-per-node=1 --gpu-bind=closest'

$SRUN "$sim"ZFP_hdf5 -z -i $GCLDLWP -c "$sim"zfp_c1.h5 -n 3 26 1800 3600 -r 11 -x serial -t s -m abs -e 0.1 2>&1 | tee -a zfp_GCLDLWP_CPU_1.txt
#$GRUN "$sim"MGARD_hdf5 -z -i $GCLDLWP -c "$sim"mgard_g1.h5 -n 3 26 1800 3600 -t s -m abs -e 0.1 -s 0 -l 0 -d hip 2>&1 | tee -a mgard_GCLDLWP_GPU_1.txt
#$SRUN "$sim"SZ_hdf5_cpu -z -i $GCLDLWP -c "$sim"sz_c1.h5 -n 3 26 1800 3600 -t s -f "$sim"sz_config/sz_16.73116.config 2>&1 | tee -a sz_GCLDLWP_CPU_1.txt
#$SRUN "$sim"SZ_hdf5_cpu -z -i $GCLDLWP -c "$sim"sz_c1.h5 -n 3 26 1800 3600 -t s -f "$sim"sz_config/sz_1.673116.config 2>&1 | tee -a sz_GCLDLWP_CPU_1.txt
#$SRUN "$sim"SZ_hdf5_cpu -z -i $GCLDLWP -c "$sim"sz_c1.h5 -n 3 26 1800 3600 -t s -f "$sim"sz_config/sz_0.1673116.config 2>&1 | tee -a sz_GCLDLWP_CPU_1.txt
#$SRUN "$sim"SZ_hdf5_cpu -z -i $GCLDLWP -c "$sim"sz_c1.h5 -n 3 26 1800 3600 -t s -f "$sim"sz_config/sz_0.01673116.config 2>&1 | tee -a sz_GCLDLWP_CPU_1.txt
#$SRUN "$sim"MGARD_hdf5 -z -i $GCLDLWP -c "$sim"mgard_c1.h5 -n 3 26 1800 3600 -t s -m abs -e 0.1 -s 0 -l 0 -d serial 2>&1 | tee -a mgard_GCLDLWP_CPU_1.txt

#$SRUN "$sim"ZFP_hdf5_double -z -i $gs -c "$sim"zfp_c1.h5 -n 3 512 512 2048 -r 11 -x serial -t d -m abs -e 0.1 2>&1 | tee -a zfp_gs_CPU_1.txt
#$GRUN "$sim"MGARD_hdf5_double -z -i $gs -c "$sim"mgard_g1.h5 -n 3 512 512 2048 -t d -m abs -e 0.1 -s 0 -l 0 -d hip 2>&1 | tee -a mgard_gs_GPU_1.txt
#$SRUN "$sim"SZ_hdf5_cpu -z -i $gs -c "$sim"sz_c1.h5 -n 3 512 512 2048 -t d -f "$sim"sz_config/sz_0.001.config 2>&1 | tee -a sz_gs_CPU_1.txt
#$SRUN "$sim"SZ_hdf5_cpu -z -i $gs -c "$sim"sz_c1.h5 -n 3 512 512 2048 -t d -f "$sim"sz_config/sz_1e-6.config 2>&1 | tee -a sz_gs_CPU_1.txt
#$SRUN "$sim"SZ_hdf5_cpu -z -i $gs -c "$sim"sz_c1.h5 -n 3 512 512 2048 -t d -f "$sim"sz_config/sz_1e-9.config 2>&1 | tee -a sz_gs_CPU_1.txt
#$SRUN "$sim"SZ_hdf5_cpu -z -i $gs -c "$sim"sz_c1.h5 -n 3 512 512 2048 -t d -f "$sim"sz_config/sz_1e-12.config 2>&1 | tee -a sz_gs_CPU_1.txt
#$SRUN "$sim"MGARD_hdf5_double -z -i $gs -c "$sim"mgard_c1.h5 -n 3 512 512 2048 -t d -m abs -e 0.1 -s inf -l 0 -d serial 2>&1 | tee -a mgard_gs_CPU_1.txt
