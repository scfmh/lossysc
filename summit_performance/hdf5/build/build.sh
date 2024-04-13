#! /bin/bash
  
export CC=`which gcc`
export CXX=`which g++`
#export LD_LIBRARY_PATH="/home/forensics/qzl/cusz/model/cusz-latest/install/lib/cmake/CUSZ/../../../lib:/home/forensics/qzl/MGARD_gpu/MGARD/install-cuda-turing/lib:$LD_LIBRARY_PATH"
#echo $LD_LIBRARY_PATH
cmake .. -DUSE_CUDA=ON\
         -DCMAKE_PREFIX_PATH="/ccs/home/zq53/lossy/SZ/cusz-latest/install/lib64/cmake/CUSZ;/ccs/home/zq53/lossy/mgard/mgard_7_23/MGARD/install-cuda-summit;/ccs/home/zq53/lossy/mgard/mgard_7_23/MGARD/install-cuda-summit/lib/cmake/zstd;/ccs/home/zq53/lossy/zfp/zfp_clean/zfp/install/lib64/cmake/zfp;/gpfs/alpine/proj-shared/csc143/zq53/mgard_new/adios_more/SZ/install/share/SZ/cmake"\
	-DCMAKE_INCLUDE_PATH="/ccs/home/zq53/lossy/zfp/zfp_clean/zfp/install/include;/ccs/home/zq53/lossy/SZ/cusz-latest/install/include/cusz;/ccs/home/zq53/IO_lib/hdf5_1_14_2/hdf5/install/include;/gpfs/alpine/proj-shared/csc143/zq53/mgard_new/adios_more/SZ/install/include"\
	-DCMAKE_LIBRARY_PATH="/ccs/home/zq53/lossy/SZ/cusz-latest/install/lib64;/ccs/home/zq53/IO_lib/hdf5_1_14_2/hdf5/install/lib;/ccs/home/zq53/lossy/mgard/mgard_7_23/MGARD/install-cuda-summit/lib;/gpfs/alpine/proj-shared/csc143/zq53/mgard_new/adios_more/SZ/install/lib64"
