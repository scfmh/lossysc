#! /bin/bash
export LD_LIBRARY_PATH="/home/forensics/qzl/MGARD_gpu/MGARD/install-cuda-turing/lib":$LD_LIBRARY_PATH

SZ=/home/forensics/qzl/qoi/xinliang/SZ3/install/bin/sz
cfg=/home/forensics/qzl/qoi/xinliang/SZ3/test
Dataset=/home/forensics/qzl/GCLDLWP_all
MGARD=/home/forensics/qzl/qoi/build/qoi_linear

/home/forensics/qzl/qoi/build/qoi_linear -z -i ../../GCLDLWP_all/GCLDLWP_1_26_1800_3600.f32  -c aaa -t s -n 3 26 1800 3600 -m abs -e 0.00167 -s 2.5 -d cuda -q 2 -l 2

$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 8.35 -s -0.5 -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_gcld_a_0.5.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 0.0835 -s -0.5 -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_gcld_a_0.5.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 0.000835 -s -0.5 -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_gcld_a_0.5.txt

$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 8.35 -s 1 -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_gcld_a_1.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 0.0835 -s 1 -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_gcld_a_1.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 0.000835 -s 1 -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_gcld_a_1.txt

$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 8.35 -s 0 -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_gcld_a_0.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 0.0835 -s 0 -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_gcld_a_0.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 0.000835 -s 0 -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_gcld_a_0.txt

$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 8.35 -s inf -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_gcld_a_inf.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 0.0835 -s inf -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_gcld_a_inf.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 0.000835 -s inf -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_gcld_a_inf.txt

$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 8.35 -s -0.5 -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_gcld_a_0.5.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 0.0835 -s -0.5 -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_gcld_a_0.5.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 0.000835 -s -0.5 -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_gcld_a_0.5.txt

$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 8.35 -s 1 -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_gcld_a_1.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 0.0835 -s 1 -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_gcld_a_1.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 0.000835 -s 1 -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_gcld_a_1.txt

$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 8.35 -s 0 -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_gcld_a_0.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 0.0835 -s 0 -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_gcld_a_0.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 0.000835 -s 0 -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_gcld_a_0.txt

$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 8.35 -s inf -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_gcld_a_inf.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 0.0835 -s inf -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_gcld_a_inf.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 0.000835 -s inf -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_gcld_a_inf.txt
#########

$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 5.8605 -s -0.5 -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_0.5.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 0.058605 -s -0.5 -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_0.5.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 0.00058605 -s -0.5 -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_0.5.txt

$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 5.8605 -s 1 -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_1.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 0.058605 -s 1 -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_1.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 0.00058605 -s 1 -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_1.txt

$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 5.8605 -s 0 -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_0.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 0.058605 -s 0 -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_0.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 0.00058605 -s 0 -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_0.txt

$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 5.8605 -s inf -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_inf.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 0.058605 -s inf -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_inf.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 0.00058605 -s inf -d serial -q 2 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_inf.txt



$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 5.8605 -s -0.5 -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_0.5.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 0.058605 -s -0.5 -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_0.5.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 0.00058605 -s -0.5 -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_0.5.txt

$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 5.8605 -s 1 -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_1.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 0.058605 -s 1 -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_1.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 0.00058605 -s 1 -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_1.txt

$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 5.8605 -s 0 -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_0.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 0.058605 -s 0 -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_0.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 0.00058605 -s 0 -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_0.txt

$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 5.8605 -s inf -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_inf.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 0.058605 -s inf -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_inf.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 0.00058605 -s inf -d serial -q 2 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_inf.txt

