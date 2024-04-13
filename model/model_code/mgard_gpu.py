from gaussian_node import generate_newnode
import numpy as np
import math
import os
import csv
from parse_time_nvprof_com import timeparse_nvprof
from perf_statistic import informationparse
from extract_info import extract_info
from total_grid_calculate import total_time, LWPK_kernel, GPK_kernel, LPK_kernel1, LPK_kernel2, LPK_kernel3_IPK 
from total_grid_calculate import Quan_kernel, Hist_kernel, Encode_kernel, Deflate_kernel, Condense_kernel
title_b = ['Lwpk', 'GpkReo3D', 'Lpk1Reo3D', 'Lpk2Reo3D','Lpk3Reo3D', 'Ipk1Reo3D','Ipk2Reo3D', 'Ipk3Reo3D', 'Quantization','Histogram', 'FillArray', 'Sort',
                         'FirstNonzeroIndex','GenerateCL','GenerateCW', 'ReverseArray', 'ReorderByIndex', 'Encode', 'Deflate', 'Condense', 'Real_quan','Outlier_num','HtoD', 'memset', 'DtoH' ]

def basic_prediction(tracefile):
    result_all = {}
    newnode = 0
    name, d1, d2, d3, err_s, d1_s, d2_s, d3_s = extract_info(tracefile)
    block_info = informationparse(tracefile)
    time_data = timeparse_nvprof(tracefile)
    # timefile = statfile.replace('stat', 'time')
    quan_s = block_info['Real_quan']
    outlier_s = block_info['Outlier_point']
    mean_s = block_info['mean']
    var_s = block_info['var'] 
    hitratio = 1-outlier_s/(d1_s*d2_s*d3_s)
    sample_lengh = d1_s * d2_s * d3_s
    all_length = d1 * d2 * d3
    if d1_s * d2_s * d3_s == d1 * d2 * d3:
        result_all  = time_data
        result_all['Real_quan']= quan_s
        result_all['Outlier_num'] = outlier_s
        newnode = quan_s
        outlier = outlier_s
        sorted_result = {key: result_all[key] for key in title_b}
        result_all = sorted_result
        print("------", newnode, outlier)
        print (result_all)
    else:
        if mean_s == 1 and var_s == 0:
            newnode = int(quan_s * 1.2 )
        else:
            newnode = int(quan_s * 1.2 )
            #newnode = generate_newnode(block_info['mean'], block_info['var'], hitratio, all_length)
            if newnode > 8192:
                newnode = 8192
        # if hf_distribution== True:
        #     newnode = generate_newnode(stats['Mean'], stats['Variance'], hitratio, all_length)   
        # else:
        #     newnode = int(quan_s * 1.2 )
        print("quantization info", quan_s, newnode)
        outlier = all_length - all_length * hitratio
        result_all['Lwpk'] = time_data['Lwpk'] / LWPK_kernel(d1_s, d2_s, d3_s, block_info['Lwpk_tbz'], block_info['Lwpk_tby'], block_info['Lwpk_tbx']) * LWPK_kernel(d1, d2, d3, block_info['Lwpk_tbz'], block_info['Lwpk_tby'], block_info['Lwpk_tbx'])
        result_all['GpkReo3D'] = time_data['GpkReo3D'] / GPK_kernel(d1_s, d2_s, d3_s, block_info['Gpk_tbz'], block_info['Gpk_tby'], block_info['Gpk_tbx']) * GPK_kernel(d1, d2, d3, block_info['Gpk_tbz'], block_info['Gpk_tby'], block_info['Gpk_tbx'])
        result_all['Lpk1Reo3D'] = time_data['Lpk1Reo3D'] / LPK_kernel1(d1_s, d2_s, d3_s, block_info['Lpk1_tbz'], block_info['Lpk1_tby'], block_info['Lpk1_tbx']) * LPK_kernel1(d1, d2, d3, block_info['Lpk1_tbz'], block_info['Lpk1_tby'], block_info['Lpk1_tbx'])
        result_all['Lpk2Reo3D'] = time_data['Lpk2Reo3D'] / LPK_kernel2(d1_s, d2_s, d3_s, block_info['Lpk2_tbz'], block_info['Lpk2_tby'], block_info['Lpk2_tbx']) * LPK_kernel2(d1, d2, d3, block_info['Lpk2_tbz'], block_info['Lpk2_tby'], block_info['Lpk2_tbx'])
        total_grid_LPK3_s, total_grid_IPK1_s , total_grid_IPK2_s ,total_grid_IPK3_s = LPK_kernel3_IPK(d1_s, d2_s, d3_s, block_info['Lpk3_tbz'], block_info['Lpk3_tby'], block_info['Lpk3_tbx'], block_info['Ipk2_tbz'], block_info['Ipk2_tby'], block_info['Ipk2_tbx'])
        total_grid_LPK3, total_grid_IPK1 , total_grid_IPK2 ,total_grid_IPK3 = LPK_kernel3_IPK(d1, d2, d3, block_info['Lpk3_tbz'], block_info['Lpk3_tby'], block_info['Lpk3_tbx'], block_info['Ipk2_tbz'], block_info['Ipk2_tby'], block_info['Ipk2_tbx'])
        result_all['Lpk3Reo3D'] = time_data['Lpk3Reo3D'] / total_grid_LPK3_s * total_grid_LPK3
        result_all['Ipk1Reo3D'] = time_data['Ipk1Reo3D'] / total_grid_IPK1_s * total_grid_IPK1
        result_all['Ipk2Reo3D'] = time_data['Ipk2Reo3D'] / total_grid_IPK2_s * total_grid_IPK2
        result_all['Ipk3Reo3D'] = time_data['Ipk3Reo3D'] / total_grid_IPK3_s * total_grid_IPK3
        result_all['Quantization'] = time_data['Quantization'] / Quan_kernel(d1_s, d2_s, d3_s, block_info['Quan_tbz'], block_info['Quan_tby'], block_info['Quan_tbx']) * Quan_kernel(d1, d2, d3, block_info['Quan_tbz'], block_info['Quan_tby'], block_info['Quan_tbx']) 
        #result_all['Histogram'] = time_data['Histogram'] / Hist_kernel(d1_s, d2_s, d3_s) * Hist_kernel(d1, d2, d3)
        #result_all['Quantization'] = time_data['Quantization'] / sample_lengh * all_length
        result_all['Histogram'] = time_data['Histogram'] / sample_lengh * all_length
        result_all['Sort'] = time_data['Sort'] / (int(math.log2(quan_s))+1) * (int(math.log2(newnode)+1))
        result_all['GenerateCL'] =  time_data['GenerateCL'] / (quan_s * math.log2(quan_s)) * (newnode * math.log2(newnode))
        result_all['GenerateCW'] = time_data['GenerateCW'] / quan_s * newnode  # sample_lengh = d1_s * d2_s * d3_s
        result_all['FirstNonzeroIndex'] = time_data['FirstNonzeroIndex'] * newnode / quan_s
        result_all['FillArray'] = time_data['FillArray'] * newnode / quan_s
        result_all['ReverseArray'] = time_data['ReverseArray'] # if not accuracy, time_data['ReverseArray']  / quan_s * newnode instead
        result_all['ReorderByIndex'] = time_data['ReorderByIndex'] #  if not accuracy, time_data['ReorderByIndex']  / quan_s * newnode instead
        result_all['Encode'] = time_data['Encode'] / Encode_kernel(d1_s, d2_s, d3_s, block_info['Encode_tbx']) * Encode_kernel(d1, d2, d3, block_info['Encode_tbx'])
        #result_all['Deflate'] = time_data['Deflate'] * Deflate_kernel(all_length - outlier, block_info['PART_SIZE'], block_info['Deflate_tbx']) / Deflate_kernel(sample_lengh - outlier_s, block_info['PART_SIZE'], block_info['Deflate_tbx'])
        result_all['Deflate'] = time_data['Deflate'] / Deflate_kernel(quan_s, block_info['PART_SIZE'], block_info['Deflate_tbx']) * Deflate_kernel(newnode, block_info['PART_SIZE'], block_info['Deflate_tbx'])
        #result_all['Deflate'] = time_data['Deflate'] * Deflate_kernel(outlier, block_info['PART_SIZE'], block_info['Deflate_tbx']) / Deflate_kernel(outlier_s, block_info['PART_SIZE'], block_info['Deflate_tbx'])
        result_all['Condense'] = time_data['Condense']  / Condense_kernel(d1_s, d2_s, d3_s) * Condense_kernel(d1, d2, d3)
        result_all['Real_quan']= newnode
        result_all['Outlier_num'] = outlier
        result_all['HtoD'] = time_data['HtoD'] / sample_lengh * all_length
        result_all['memset'] = time_data['memset'] / sample_lengh * all_length
        result_all['DtoH'] = time_data['DtoH'] / sample_lengh / quan_s * all_length * newnode
        print(result_all)
    return name, d1_s, d2_s, d3_s, err_s, newnode, outlier, result_all

if __name__ == '__main__':
    result_dict = {}
    path = "/home/forensics/qzl/gpu_model/mgard_trace"
    dirs = os.listdir(path)
    title = ['Name','Lwpk', 'GpkReo3D', 'Lpk1Reo3D', 'Lpk2Reo3D','Lpk3Reo3D', 'Ipk1Reo3D','Ipk2Reo3D', 'Ipk3Reo3D', 'Quantization','Histogram', 'FillArray', 'Sort',
                 'FirstNonzeroIndex','GenerateCL','GenerateCW', 'ReverseArray', 'ReorderByIndex', 'Encode', 'Deflate', 'Condense', 'Real_quan','Outlier_num','HtoD', 'memset', 'DtoH' ]
    for file in dirs:
        if '.txt' in file:
            print(file)
            full_file = path + "/" + file   
            name_s, d1_s, d2_s, d3_s, err_s, newnode, outlier, result_all = basic_prediction(full_file)
            full_name = (str(name_s) + "_" + str(d1_s) + "_" + str(d2_s) + "_" + str(d3_s)).replace(path+"/","") 
            write_file = str(name_s).replace(path+"/","") +"_" + str(err_s) + ".csv"
            row_info = [str(full_name)]  + [result_all[key] for key in result_all]
            file_exists = False
            try:
                with open(write_file, 'r') as csvfile:
                    reader = csv.reader(csvfile)
                    for row in reader:
                        if row[0] == "Name":
                            file_exists = True
                            break
            except FileNotFoundError:
                pass
            with open(write_file, 'a', newline='') as csvfile:
                writer = csv.writer(csvfile)
    
                if not file_exists:
                    #title = ['Name'] + title
                    writer.writerow(title)
                
                writer.writerow(row_info)
