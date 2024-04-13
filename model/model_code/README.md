A. General description:
general.py
gaussian_node.py is used for generation new huffman tree based on the sampe mean and var

B. SZ prediction:
extract_info.py get the basic description of the file.
parse_sz_information_nvprof.py gets the sample performance based on the nvprof
perf_sz_stat.py can get the extra information like hit ratio and quantiozation
1. for the SZ model, we need to meature the quantization distribution stat, so mean and variance will be copied to
   the host side, we will to #define ZHENLU_STATICS in the file src/hf/detail/hf_bookg.inl or not
2. if we define ZHENLU_STATICS, extra time will be added in the report, so we can not use it, but only for the stat 
   calculation, for consistency, we named it as following format
   stat format: Dataset_d1_d2_d3_err_stat_sam.txt
3. if we don't define it, the results are used for the sampled time and named as follows
   time format: Dataset_d1_d2_d3_err_time_sam.txt
   baselines are Dataset_d1_d2_d3_err_stat_1.txt and Dataset_d1_d2_d3_err_time_1.txt

C. MGARD prediction:
total_grid_calculate.py is used for the MGARD special grids with non scalability operations
parse_time_nvprof.py is used for the MGRAD unit time descriptions.
MGARD GPU model threads are related to the dimensions description, so the name should be slightly different with 
SZ, we also need to name it includes the sample demisions. 
1. for the MGARD model, we need to meature the quantization distribution stat, so mean and variance will be copied to
   the host side, we will to #define HF_STATICS in the file ./include/mgard-x/Lossless/ParallelHuffman/GetCodebook.hpp 
   or not (Quantization_static)
2. because we need to get the sample demisions, so we do not need to spificy each sample demisions and ignore the sample ratio
   e.g., sample dataset: GCLDLWP_26_1800_3600_err_13_900_1800.txt
         full dataset: GCLDLWP_26_1800_3600_err_26_1800_3600.txt