#! /bin/bash

export CC=`which gcc`    
export CXX=`which g++`

cmake .. -DCMAKE_PREFIX_PATH="/home/forensics/qzl/MGARD_gpu/MGARD/install-cuda-turing;/home/forensics/qzl/MGARD_gpu/MGARD/install-cuda-turing/lib/cmake/zstd"\

