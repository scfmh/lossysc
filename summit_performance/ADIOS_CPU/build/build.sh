#! /bin/bash

cmake .. -DCMAKE_PREFIX_PATH="/MGARD/install-path/;/MGARD/install-path/lib/cmake/zstd;/ZFP/install-path/lib64/cmake/zfp;/SZ/install-path;/ADIOS2/install-path/lib64/cmake/adios2"\
         -DCMAKE_INCLUDE_PATH="/SZ/install-path/include;/ZFP/install-path/include"

