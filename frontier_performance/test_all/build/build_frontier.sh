#! /bin/bash
  
#export CC=`which gcc`
#export CXX=`which g++`
module load rocm/5.4.0
module load cray-mpich/8.1.17

export CC=amdclang
export CXX=amdclang++


cmake .. -DUSE_HIP=ON\
         -DCMAKE_PREFIX_PATH="/MGARD/install-path;/MGARD/install-path/lib/cmake/zstd;/ZFP-install-path/lib64/cmake/zfp;/SZ-install-path"\
         -DCMAKE_INCLUDE_PATH="/SZ-install-path/include;/ZFP-install-path/include"\
	     -DCMAKE_LIBRARY_PATH="/MGARD/install-path/lib/libzstd.so"
