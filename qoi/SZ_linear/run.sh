#! /bin/bash
export LD_LIBRARY_PATH="/home/forensics/qzl/MGARD_gpu/MGARD/install-cuda-turing/lib":$LD_LIBRARY_PATH

SZ=/home/forensics/qzl/qoi/xinliang/SZ3/install/bin/sz
cfg=/home/forensics/qzl/qoi/xinliang/SZ3/linear
Dataset=/home/forensics/qzl/GCLDLWP_all
MGARD=/home/forensics/qzl/qoi/build/qoi_average

$SZ -f -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -o test.dat.sz.out -3 26 1800 3600 -c "$cfg"/sz.config_gcld1 -M ABS 133.6 -a 2>&1 | tee -a 3d_gcld_a.txt
$SZ -f -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -o test.dat.sz.out -3 26 1800 3600 -c "$cfg"/sz.config_gcld2 -M ABS 1.336 -a 2>&1 | tee -a 3d_gcld_a.txt
$SZ -f -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -o test.dat.sz.out -3 26 1800 3600 -c "$cfg"/sz.config_gcld3 -M ABS 0.01336 -a 2>&1 | tee -a 3d_gcld_a.txt

$SZ -f -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -o test.dat.sz.out -3 3600 1800 26 -c "$cfg"/sz.config_gcld1 -M ABS 133.6 -a 2>&1 | tee -a 3d_gcld_d.txt
$SZ -f -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -o test.dat.sz.out -3 3600 1800 26 -c "$cfg"/sz.config_gcld2 -M ABS 1.336 -a 2>&1 | tee -a 3d_gcld_d.txt
$SZ -f -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -o test.dat.sz.out -3 3600 1800 26 -c "$cfg"/sz.config_gcld3 -M ABS 0.01336 -a 2>&1 | tee -a 3d_gcld_d.txt

$SZ -f -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -o test.dat.sz.out -2 3600 46800 -c "$cfg"/sz.config_gcld1 -M ABS 133.6 -a 2>&1 | tee -a 2d_gcld_a.txt
$SZ -f -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -o test.dat.sz.out -2 3600 46800 -c "$cfg"/sz.config_gcld2 -M ABS 1.336 -a 2>&1 | tee -a 2d_gcld_a.txt
$SZ -f -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -o test.dat.sz.out -2 3600 46800 -c "$cfg"/sz.config_gcld3 -M ABS 0.01336 -a 2>&1 | tee -a 2d_gcld_a.txt

$SZ -f -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -o test.dat.sz.out -2 46800 3600 -c "$cfg"/sz.config_gcld1 -M ABS 133.6 -a 2>&1 | tee -a 2d_gcld_d.txt
$SZ -f -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -o test.dat.sz.out -2 46800 3600 -c "$cfg"/sz.config_gcld2 -M ABS 1.336 -a 2>&1 | tee -a 2d_gcld_d.txt
$SZ -f -i "$Dataset"/GCLDLWP_1_26_1800_3600.f32 -o test.dat.sz.out -2 46800 3600 -c "$cfg"/sz.config_gcld3 -M ABS 0.01336 -a 2>&1 | tee -a 2d_gcld_d.txt

$SZ -f -i "$Dataset"/dataset2-2338x106711.x.f32.dat -o test.dat.sz.out -3 2338 109 979 -c "$cfg"/sz.config_exaalt1 -M ABS 93.768 -a 2>&1 | tee -a 3d_exxalt_d.txt
$SZ -f -i "$Dataset"/dataset2-2338x106711.x.f32.dat -o test.dat.sz.out -3 2338 109 979 -c "$cfg"/sz.config_exaalt2 -M ABS 0.93768 -a 2>&1 | tee -a 3d_exxalt_d.txt
$SZ -f -i "$Dataset"/dataset2-2338x106711.x.f32.dat -o test.dat.sz.out -3 2338 109 979 -c "$cfg"/sz.config_exaalt3 -M ABS 0.0093768 -a 2>&1 | tee -a 3d_exxalt_d.txt

$SZ -f -i "$Dataset"/dataset2-2338x106711.x.f32.dat -o test.dat.sz.out -3 109 979 2338 -c "$cfg"/sz.config_exaalt1 -M ABS 93.768 -a 2>&1 | tee -a 3d_exxalt_a.txt
$SZ -f -i "$Dataset"/dataset2-2338x106711.x.f32.dat -o test.dat.sz.out -3 109 979 2338 -c "$cfg"/sz.config_exaalt2 -M ABS 0.93768 -a 2>&1 | tee -a 3d_exxalt_a.txt
$SZ -f -i "$Dataset"/dataset2-2338x106711.x.f32.dat -o test.dat.sz.out -3 109 979 2338 -c "$cfg"/sz.config_exaalt3 -M ABS 0.0093768 -a 2>&1 | tee -a 3d_exxalt_a.txt

$SZ -f -i "$Dataset"/dataset2-2338x106711.x.f32.dat -o test.dat.sz.out -2 2338 106711 -c "$cfg"/sz.config_exaalt1 -M ABS 93.768 -a 2>&1 | tee -a 2d_exxalt_a.txt
$SZ -f -i "$Dataset"/dataset2-2338x106711.x.f32.dat -o test.dat.sz.out -2 2338 106711 -c "$cfg"/sz.config_exaalt2 -M ABS 0.93768 -a 2>&1 | tee -a 2d_exxalt_a.txt
$SZ -f -i "$Dataset"/dataset2-2338x106711.x.f32.dat -o test.dat.sz.out -2 2338 106711 -c "$cfg"/sz.config_exaalt3 -M ABS 0.0093768 -a 2>&1 | tee -a 2d_exxalt_a.txt

$SZ -f -i "$Dataset"/dataset2-2338x106711.x.f32.dat -o test.dat.sz.out -2 106711 2338 -c "$cfg"/sz.config_exaalt1 -M ABS 93.768 -a 2>&1 | tee -a 2d_exxalt_d.txt
$SZ -f -i "$Dataset"/dataset2-2338x106711.x.f32.dat -o test.dat.sz.out -2 106711 2338 -c "$cfg"/sz.config_exaalt2 -M ABS 0.93768 -a 2>&1 | tee -a 2d_exxalt_d.txt
$SZ -f -i "$Dataset"/dataset2-2338x106711.x.f32.dat -o test.dat.sz.out -2 106711 2338 -c "$cfg"/sz.config_exaalt3 -M ABS 0.0093768 -a 2>&1 | tee -a 2d_exxalt_d.txt
##############


