import pandas as pd
import math
import numpy as np

######## some special operations which have different number of threads will show in this file
######## others for the fixed threads show the scalability between the sample and full dataset
######## Fillarray, GetFirstNonzeroIndex, ReverseArray, ReorderByIndex are related to the huffman dict size
######## if size is fixed 8192, blocks are 9. tbx is related to the GPU. 
######## Condense tbx is 256 and will divide 80 for the basic

def total_time(d1, d2, d3):
    max_time = min(d1, d2, d3)
    #print(int(math.log2(max_time))+1)
    return int(math.log2(max_time))+1

# input dimensions d1, d2 ,d3 13 45 45 # 64 2 2 
# first and last one time and middle parts are double operations
def LWPK_kernel(d1, d2, d3, tz, ty, tx): 
    total_grid_LWPK = 0
    new_d1 = d1 
    new_d2 = d2 
    new_d3 = d3 
    total_grid_LWPK += math.ceil(new_d1 / tz) * math.ceil(new_d2 / ty)  * math.ceil(new_d3 /tx)
    for i in range(total_time(d1, d2, d3), 1, -1):
        new_d1 = int(new_d1 / 2) + 1
        new_d2 = int(new_d2 / 2) + 1
        new_d3 = int(new_d3 / 2) + 1
        total_grid_LWPK += math.ceil(new_d1 / tz) * math.ceil(new_d2 / ty) * math.ceil(new_d3 /tx) * 2
    new_d1 = int(new_d1 / 2) + 1
    new_d2 = int(new_d2 / 2) + 1
    new_d3 = int(new_d3 / 2) + 1
    total_grid_LWPK += math.ceil(new_d1 / tz)  * math.ceil(new_d2 / ty)  * math.ceil(new_d3 /tx)
    return total_grid_LWPK


def GPK_kernel(d1, d2, d3, tz, ty, tx): # input dimensions d1, d2 ,d3 13 45 45 # 64 2 2
    total_grid_GPK = 0
    new_d1 = d1 - 1
    new_d2 = d2 - 1
    new_d3 = d3 - 1
    total_grid_GPK += math.ceil(new_d1 / tz) * math.ceil(new_d2 / ty) * math.ceil(new_d3 /tx)
    for i in range(total_time(d1, d2, d3), 1, -1):
        new_d1 = math.ceil(new_d1 / 2)
        new_d2 = math.ceil(new_d2 / 2)
        new_d3 = math.ceil(new_d3 / 2)
        total_grid_GPK += math.ceil(new_d1 / tz) * math.ceil(new_d2 / ty) * math.ceil(new_d3 /tx)
    return total_grid_GPK

def LPK_kernel1(d1, d2, d3, tz, ty, tx): # input dimensions d1, d2 ,d3 13 45 45 # 4 4 16
    total_grid_LPK1 = 0
    new_d1 = d1
    new_d2 = d2
    new_d3 = int(d3 / 2) + 1
    total_grid_LPK1 += math.ceil(new_d1 / tz) * math.ceil(new_d2 / ty) * math.ceil(new_d3 /tx)
    for i in range(total_time(d1, d2, d3), 1, -1):
        new_d1 = int(new_d1 / 2) + 1
        new_d2 = int(new_d2 / 2) + 1
        new_d3 = int(new_d3 / 2) + 1
        total_grid_LPK1 += math.ceil(new_d1 / tz) * math.ceil(new_d2 / ty) * math.ceil(new_d3 /tx)
    return total_grid_LPK1

def LPK_kernel2(d1, d2, d3, tz, ty, tx): # input dimensions d1, d2 ,d3 13 45 45 # 4 4 16
    total_grid_LPK2 = 0
    new_d1 = d1
    new_d2 = int(d2 / 2) + 1
    new_d3 = int(d3 / 2) + 1
    total_grid_LPK2 += math.ceil(new_d1 / tz) * math.ceil(new_d2 / ty) * math.ceil(new_d3 /tx)
    for i in range(total_time(d1, d2, d3), 1, -1):
        new_d1 = int(new_d1 / 2) + 1
        new_d2 = int(new_d2 / 2) + 1
        new_d3 = int(new_d3 / 2) + 1
        total_grid_LPK2 += math.ceil(new_d1 / tz) * math.ceil(new_d2 / ty) * math.ceil(new_d3 /tx)
    return total_grid_LPK2
# input dimensions d1, d2 ,d3 13 45 45 # 4 4 8
#ipk1 1 4 32  ipk2 1 2 128 ipk3 1 2 128
def LPK_kernel3_IPK(d1, d2, d3, tz, ty, tx, tz2, ty2, tx2 ):
    total_grid_LPK3 = 0
    total_grid_IPK1 = 0
    total_grid_IPK2 = 0
    total_grid_IPK3 = 0
    new_d1 = int(d1 / 2) + 1
    new_d2 = int(d2 / 2) + 1
    new_d3 = int(d3 / 2) + 1
    total_grid_LPK3 += math.ceil(new_d1 / tz) * math.ceil(new_d2 / ty) * math.ceil(new_d3 /tx)
    total_grid_IPK1 += math.ceil(new_d1 / tz) * math.ceil(new_d2 / ty)
    total_grid_IPK2 += math.ceil(new_d1 / ty2) * math.ceil(new_d3 / tx2)
    total_grid_IPK3 += math.ceil(new_d2 / ty2) * math.ceil(new_d3 / tx2)
    for i in range(total_time(d1, d2, d3), 1, -1):
        new_d1 = int(new_d1 / 2) + 1
        new_d2 = int(new_d2 / 2) + 1
        new_d3 = int(new_d3 / 2) + 1
        total_grid_LPK3 += math.ceil(new_d1 / tz) * math.ceil(new_d2 / ty) * math.ceil(new_d3 /tx)
        total_grid_IPK1 += math.ceil(new_d1 / tz) * math.ceil(new_d2 / ty)
        total_grid_IPK2 += math.ceil(new_d1 / ty2) * math.ceil(new_d3 / tx2)
        total_grid_IPK3 += math.ceil(new_d2 / ty2) * math.ceil(new_d3 / tx2)
    return total_grid_LPK3, total_grid_IPK1 , total_grid_IPK2 ,total_grid_IPK3


def Quan_kernel(d1, d2, d3, tz, ty, tx ): # 4 , 4, 16
    total_grid_quan = 0
    total_grid_quan +=  math.ceil(d1 / tz) * math.ceil(d2 / ty) * math.ceil(d3 /tx)
    return total_grid_quan

def Hist_kernel(d1, d2, d3):
    if (d1 * d2 * d3) > 1024 * 80:
        total_grid_hist = 160
    else:
        total_grid_hist = 80
    return total_grid_hist

def Hist_kernel(d1, d2, d3):
    if (d1 * d2 * d3) > 1024 * 80:
        total_grid_hist = 160
    else:
        total_grid_hist = 80
    return total_grid_hist

def Encode_kernel(d1, d2, d3, tx ): # 4 , 4, 16
    total_grid_encode = 0
    total_grid_encode +=  math.ceil(d1 * d2 * d3 / tx)
    return total_grid_encode

#def Deflate_kernel(d1, d2, d3, PART_SIZE, tbx): ##PART_SIZE: 20480 tbx = 128
    #PART_SIZE = 20480
    #total_grid_deflate = int((int((d1 * d2 * d3 - 1) / PART_SIZE ) +1 ) / tbx ) + 1
    #return total_grid_deflate

def Deflate_kernel(real_q, PART_SIZE, tbx): ##PART_SIZE: 20480 tbx = 128
    PART_SIZE = 20480
    total_grid_deflate = int((int((real_q - 1) / PART_SIZE ) +1 ) / tbx ) + 1
    return total_grid_deflate

def Condense_kernel(d1, d2, d3):
    total_grid_condense = int(d1 * d2 * d3 / 80 / 256 ) + 1
    return total_grid_condense


# GPK_kernel(13, 45, 225, 2, 2, 64)
# LPK_kernel1(13, 45, 225, 4, 4, 16)
# LPK_kernel2(13, 45, 225, 4, 4, 16)
# LPK_kernel3_IPK(13, 45, 225, 4, 4, 8, 1,2,128)
# LPK_kernel3_IPK(13, 45, 45, 4, 4, 8, 1,2,128)
#
# GPK_kernel(26, 1800,3600, 2, 2, 64)
# LPK_kernel1(26, 1800,3600, 4, 4, 16)
# LPK_kernel2(26, 1800,3600, 4, 4, 16)
# LPK_kernel3_IPK(26, 1800,3600, 4, 4, 8, 1,2,128)
# Quan_kernel(26, 1800,3600, 4, 4, 16)

# Delate_kernel(26, 1800,3600,20480, 128)
# Delate_kernel(13, 45,45, 20480,128)
# Delate_kernel(13, 900,1800,20480, 128)
