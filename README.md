epository for Testing Lossy Compressors on Different HPC Systems

This repository is set up for testing various lossy compressors under different high-performance computing (HPC) systems. Before testing, the following dependency libraries and their versions must be downloaded and used:

- **ADIOS2 (v2.9.0)**: [https://github.com/ornladios/ADIOS2.git](https://github.com/ornladios/ADIOS2.git)
- **HDF5 (1.14.1)**: [https://github.com/HDFGroup/hdf5.git](https://github.com/HDFGroup/hdf5.git)
- **SZ (2.1.12.5)**: [https://github.com/szcompressor/SZ.git](https://github.com/szcompressor/SZ.git)
- **Cusz (0.3.1)**: [https://github.com/szcompressor/cuSZ.git](https://github.com/szcompressor/cuSZ.git)
- **ZFP (1.0.0)**: [https://github.com/LLNL/zfp.git](https://github.com/LLNL/zfp.git)
- **MGARD (1.5.0)**: [https://github.com/CODARcode/MGARD.git](https://github.com/CODARcode/MGARD.git)

## A. Performance on Summit

### Setup
1. Go to the folder `summit_performance` and select either ADIOS or HDF5.
2. Load the related libraries as follows:
    ```bash
    module load cmake/3.18.4
    module load spectrum-mpi/10.4.0.3-20210112
    module load python/3.8-anaconda3
    module load cuda/11.4
    module load gcc/9
    module load zstd
    ```
3. Navigate to the related folder, open the build, and edit `build.sh` based on your dependency path.
4. Run the build script:
    ```bash
    bash build.sh
    ```

### Testing Examples
For example, using a generic file name `gc=GCLDLWP_1_26_1800_3600.f32`:
- **MGARD GPU**:
    ```bash
    ./build/cpu-application-simulator_GCLDLWP -z 1 -i ${gc} -c output.bp -n 3 3600 1800 26 -t s -s 0 -l 0 -m abs -e 0.1 -d cuda
    ```
- **ZFP GPU**:
    ```bash
    ./build/SZ_cpu_simulation_cuda_iter_ICLDIWP -i ${gc} -c output.bp -n 3 3600 1800 26 -e 16.73116
    ```

## B. Performance on Frontier

### Setup
1. Load the related libraries as follows:
    ```bash
    module load rocm/5.1.0
    module load cmake
    module load cray-mpich/8.1.17
    module load cray-python/3.9.12.1
    module load zstd
    ```
2. Navigate to the related folder, open the build, and edit `build.sh` based on your dependency path.
3. Run the build script:
    ```bash
    bash build.sh
    ```

### Testing Examples
- **ZFP CPU**:
    ```bash
    ./build/ZFP_cpu_simulation_iter_exaalt -z 1 -i ${gc} -c output.bp -t s -n 3 3600 1800 26 -m abs -e 0.1
    ```
- **MGARD CPU**:
    ```bash
    ./build/cpu-application-simulator_exaalt -z 1 -i ${gc} -c output.bp -t s -n 3 3600 1800 26 -m abs -e 0.1 -l 0 -s 0 -d serial
    ```
- **SZ CPU**:
    ```bash
    ./build/SZ_cpu_simulation_ICLDIWP -z 1 -i ${gc} -c output.bp -n 3 3600 1800 26 -m abs -e 0.1 -f sz_16.73116.config
    ```

## C. Testing Parallel Performance

For Adios and HDF5, include the following preprocessor directives in your source files to enable MPI:
- **Adios**:
    ```c
    #define MYADIOS_USE_MPI 1
    #include <adios2.h>
    #if MYADIOS_USE_MPI
    #include <mpi.h>
    #endif
    ```
- **HDF5**:
    ```c
    #define HDF5_USE_MPI 1
    #include <hdf5.h>
    #if HDF5_USE_MPI
    #include <mpi.h>
    #include "H5FDmpio.h"
    #endif
    ```

## D. Testing GPU and CPU Power

### GPU Power
Use `nvprof` for profiling GPU applications for example:
```bash
nvprof --system-profiling on --print-gpu-trace --log-file mgard.log $mgard -z -i {gc} -c temp.mgard -t s -n 3 3600 1800 26 -m abs -e 16.7 -s inf -l 0 -h 0 -d cuda -v 3
```
### CPU Power
require root privileges, we can build the dependency based on their offical scripts.
in `power/SZ-2.1.12.5/sz/src/utility.c`, we add the function like `start_measure_cpu` and we can call them when we want to use them.

For example. we can use the command line 
```bash
sudo perf stat -e power/energy-cores/ -o sz_GCLDLWP.txt ./build/bin/sz -z -f -c sz.config16.7 -i GCLDLWP_1_26_1800_3600.f32 -3 3600 1800 26
```

## E. QoI testing
For MGARD: open the build and edit the build.sh details
 ```bash
 run bash build.sh
```
For SZ: build the SZ-qoi/SZ3 based on their offical scripts
Run: you can find the example 
   ```bash 
   $MGARD -z -i ${gc} -c aaa -t s -n 3 3600 1800 26 -m abs -e 16 -s -0.5 -d serial
     SZ: $SZ -f -i ${gc} -o test.dat.sz.out -3 26 1800 3600 -c sz.config_gcld1 -M ABS 133.6 -a
```

## F. Modeling

 run the compressor first and get the basic components,(modify the datapath and simulation path) like:
```bash
#! /bin/bash
Datapath=/home/forensics/qzl/gpu_model/EXAFEL
Sim=/home/forensics/qzl/gpu_model/MGARD/build-cuda-turing/mgard/bin
#16481.0321045
e_values=("1648.104" "164.8104" "16.48104" "1.648104")
#e_values=("1648.104")
for e in "${e_values[@]}"; do
nvprof $Sim/mgard-x -z -i ${Datapath}/SDRBENCH-EXAFEL-data-130x1480x1552.f32 -c aaa -m abs -e $e -t s -n 3 130 1480 1552 -s inf -l 0 -d cuda 2>&1 | tee EXAFEL_130_1480_1552_${e}_130_1480_1552.txt
nvprof $Sim/mgard-x -z -i ${Datapath}/SDRBENCH-EXAFEL-data-65x740x776.f32 -c aaa -m abs -e $e -t s -n 3 65 740 776 -s inf -l 0 -d cuda 2>&1 | tee EXAFEL_130_1480_1552_${e}_65_740_776.txt
nvprof $Sim/mgard-x -z -i ${Datapath}/SDRBENCH-EXAFEL-data-65x185x388.f32 -c aaa -m abs -e $e -t s -n 3 65 185 388 -s inf -l 0 -d cuda 2>&1 | tee EXAFEL_130_1480_1552_${e}_65_185_388.txt
nvprof $Sim/mgard-x -z -i ${Datapath}/SDRBENCH-EXAFEL-data-65x97x185.f32 -c aaa -m abs -e $e -t s -n 3 65 97 185 -s inf -l 0 -d cuda 2>&1 | tee EXAFEL_130_1480_1552_${e}_65_97_185.txt
nvprof $Sim/mgard-x -z -i ${Datapath}/SDRBENCH-EXAFEL-data-13x97x185.f32 -c aaa -m abs -e $e -t s -n 3 13 97 185 -s inf -l 0 -d cuda 2>&1 | tee EXAFEL_130_1480_1552_${e}_13_97_185.txt
done
```
and then 
```bash
run python mgard_gpu.py
```

