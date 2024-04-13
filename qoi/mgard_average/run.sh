#! /bin/bash
export LD_LIBRARY_PATH="/home/forensics/qzl/MGARD_gpu/MGARD/install-cuda-turing/lib":$LD_LIBRARY_PATH

SZ=/home/forensics/qzl/qoi/xinliang/SZ3/install/bin/sz
cfg=/home/forensics/qzl/qoi/xinliang/SZ3/test
Dataset=/home/forensics/qzl/GCLDLWP_all
MGARD=/home/forensics/qzl/qoi/build/qoi_average


$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 16.7 -s 0 -d serial -q 6 0 0 0 3600 1800 26 -l 2 2>&1 | tee -a mgard_3d_gcld_a_0.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 0.167 -s 0 -d serial -q 6 0 0 0 3600 1800 26 -l 2 2>&1 | tee -a mgard_3d_gcld_a_0.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 0.00167 -s 0 -d serial -q 6 0 0 0 3600 1800 26 -l 2 2>&1 | tee -a mgard_3d_gcld_a_0.txt

$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 16.7 -s 0 -d serial -q 4 0 0 46800 3600 -l 2 2>&1 | tee -a mgard_2d_gcld_a_0.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 0.167 -s 0 -d serial -q 4 0 0 46800 3600 -l 2 2>&1 | tee -a mgard_2d_gcld_a_0.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 0.00167 -s 0 -d serial -q 4 0 0 46800 3600 -l 2 2>&1 | tee -a mgard_2d_gcld_a_0.txt

#
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 16.7 -s -0.5 -d serial -q 6 0 0 0 3600 1800 26 -l 2 2>&1 | tee -a mgard_3d_gcld_a_5.txt

$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 0.167 -s -0.5 -d serial -q 6 0 0 0 3600 1800 26 -l 2 2>&1 | tee -a mgard_3d_gcld_a_5.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 0.00167 -s -0.5 -d serial -q 6 0 0 0 3600 1800 26 -l 2 2>&1 | tee -a mgard_3d_gcld_a_5.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 16.7 -s -0.5 -d serial -q 4 0 0 46800 3600 -l 2 2>&1 | tee -a mgard_2d_gcld_a_5.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 0.167 -s -0.5 -d serial -q 4 0 0 46800 3600 -l 2 2>&1 | tee -a mgard_2d_gcld_a_4.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 0.00167 -s -0.5 -d serial -q 4 0 0 46800 3600 -l 2 2>&1 | tee -a mgard_2d_gcld_a_5.txt
###
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 16.7 -s 1 -d serial -q 6 0 0 0 3600 1800 26 -l 2 2>&1 | tee -a mgard_3d_gcld_a_1.txt

$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 0.167 -s 1 -d serial -q 6 0 0 0 3600 1800 26 -l 2 2>&1 | tee -a mgard_3d_gcld_a_1.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 0.00167 -s 1 -d serial -q 6 0 0 0 3600 1800 26 -l 2 2>&1 | tee -a mgard_3d_gcld_a_1.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 16.7 -s 1 -d serial -q 4 0 0 46800 3600 -l 2 2>&1 | tee -a mgard_2d_gcld_a_1.txt

$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 0.167 -s 1 -d serial -q 4 0 0 46800 3600 -l 2 2>&1 | tee -a mgard_2d_gcld_a_1.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 0.00167 -s 1 -d serial -q 4 0 0 46800 3600 -l 2 2>&1 | tee -a mgard_2d_gcld_a_1.txt
#
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 16.7 -s inf -d serial -q 6 0 0 0 3600 1800 26 -l 2 2>&1 | tee -a mgard_3d_gcld_a_inf.txt

$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 0.167 -s inf -d serial -q 6 0 0 0 3600 1800 26 -l 2 2>&1 | tee -a mgard_3d_gcld_a_inf.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 3 26 1800 3600 -m abs -e 0.00167 -s inf -d serial -q 6 0 0 0 3600 1800 26 -l 2 2>&1 | tee -a mgard_3d_gcld_a_inf.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 16.7 -s inf -d serial -q 4 0 0 46800 3600 -l 2 2>&1 | tee -a mgard_2d_gcld_a_inf.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 0.167 -s inf -d serial -q 4 0 0 46800 3600 -l 2 2>&1 | tee -a mgard_2d_gcld_a_inf.txt
$MGARD -z -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -c aaa -t s -n 2 3600 46800 -m abs -e 0.00167 -s inf -d serial -q 4 0 0 46800 3600 -l 2 2>&1 | tee -a mgard_2d_gcld_a_inf.txt

#########
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 11.721 -s 0 -d serial -q 6 0 0 0 2338 979 109 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_0.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 0.11721 -s 0 -d serial -q 6 0 0 0 2338 979 109 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_0.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 0.0011721 -s 0 -d serial -q 6 0 0 0 2338 979 109 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_0.txt

$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 11.721 -s 0 -d serial -q 4 0 0 106711 2338 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_0.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 0.11721 -s 0 -d serial -q 4 0 0 106711 2338 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_0.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 0.0011721 -s 0 -d serial -q 4 0 0 106711 2338 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_0.txt



$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 11.721 -s -0.5 -d serial -q 6 0 0 0 2338 979 109 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_5.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 0.11721 -s -0.5 -d serial -q 6 0 0 0 2338 979 109 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_5.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 0.0011721 -s -0.5 -d serial -q 6 0 0 0 2338 979 109 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_5.txt

$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 11.721 -s -0.5 -d serial -q 4 0 0 106711 2338 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_5.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 0.11721 -s -0.5 -d serial -q 4 0 0 106711 2338 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_5.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 0.0011721 -s -0.5 -d serial -q 4 0 0 106711 2338 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_5.txt

$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 11.721 -s 1 -d serial -q 6 0 0 0 2338 979 109 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_1.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 0.11721 -s 1 -d serial -q 6 0 0 0 2338 979 109 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_1.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 0.0011721 -s 1 -d serial -q 6 0 0 0 2338 979 109 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_1.txt

$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 11.721 -s 1 -d serial -q 4 0 0 106711 2338 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_1.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 0.11721 -s 1 -d serial -q 4 0 0 106711 2338 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_1.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 0.0011721 -s 1 -d serial -q 4 0 0 106711 2338 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_1.txt

$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 11.721 -s inf -d serial -q 6 0 0 0 2338 979 109 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_inf.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 0.11721 -s inf -d serial -q 6 0 0 0 2338 979 109 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_inf.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 3 109 979 2338 -m abs -e 0.0011721 -s inf -d serial -q 6 0 0 0 2338 979 109 -l 2 2>&1 | tee -a mgard_3d_exaalt_a_inf.txt

$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 11.721 -s inf -d serial -q 4 0 0 106711 2338 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_inf.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 0.11721 -s inf -d serial -q 4 0 0 106711 2338 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_inf.txt
$MGARD -z -i "$Dataset"/dataset2-2338x106711.x.f32.dat -c aaa -t s -n 2 2338 106711 -m abs -e 0.0011721 -s inf -d serial -q 4 0 0 106711 2338 -l 2 2>&1 | tee -a mgard_2d_exaalt_a_inf.txt

