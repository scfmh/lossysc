#! /bin/bash

export CC=`which gcc`    
export CXX=`which g++`

cmake .. -DCMAKE_PREFIX_PATH="/ADIOS2/install-path/lib64/cmake/adios2;/MGARD/install-path/;/MGARD/install-path/lib/cmake/zstd;/ZFP/install-path/lib64/cmake/zfp;/CUSZ/install-path/lib64/cmake/CUSZ;/CUDA-install/cuda/11.4.2/include/cub/cmake"\
         -DCMAKE_INCLUDE_PATH="/ZFP/install-path/include;/CUSZ/install-path/include/cusz;"\
         -DCMAKE_LIBRARY_PATH="/CUSZ/install-path/lib64"

