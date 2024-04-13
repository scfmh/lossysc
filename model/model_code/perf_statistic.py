from general import remove_symbols
def informationparse(inputfile):
    result_all = {}
    Gpk_tbz, Gpk_tby, Gpk_tbx = 0, 0, 0
    Lpk1_tbz,  Lpk1_tby,  Lpk1_tbx = 0, 0, 0
    Lpk2_tbz,  Lpk2_tby,  Lpk2_tbx = 0, 0, 0
    Lpk3_tbz,  Lpk3_tby,  Lpk3_tbx = 0, 0, 0
    Ipk1_tbz,  Ipk1_tby,  Ipk1_tbx = 0, 0, 0
    Ipk2_tbz,  Ipk2_tby,  Ipk2_tbx = 0, 0, 0
    Ipk3_tbz,  Ipk3_tby,  Ipk3_tbx = 0, 0, 0
    Hist_tbz,  Hist_tby,  Hist_tbx =  0, 0, 0
    Quan_tbz,  Quan_tby,  Quan_tbx = 0, 0, 0
    mean, var = 0.0, 0.0
    Real_quan = 0
    Encode_tbx = 0
    Deflate_tbx = 0 
    Outlier_point = 0
    PART_SIZE = 0

    Gpk_done = False
    Lwpk_done = False
    Lpk1_done = False
    Lpk2_done = False
    Lpk3_done = False
    Ipk1_done = False
    Ipk2_done = False
    Ipk3_done = False
    Quan_done = False

    with open(inputfile, 'r') as rf:
        line = rf.readline()
        while line:
            if "Gpk block" in line and not Gpk_done:
                line = line.replace(",","")
                temp = list(line.split(" "))
                result_all['Gpk_tbz'], result_all['Gpk_tby'], result_all['Gpk_tbx'] = int(temp[3]),  int(temp[5]), int(temp[7])
                Gpk_done = True
            elif "Lwpk block" in line and not Lwpk_done:
                line = line.replace(",","")
                temp = list(line.split(" "))
                result_all['Lwpk_tbz'], result_all['Lwpk_tby'], result_all['Lwpk_tbx'] = int(temp[3]),  int(temp[5]), int(temp[7])
                Lwpk_done = True
            elif "Lpk1 block" in line and not Lpk1_done:
                line = line.replace(",","")
                temp = list(line.split(" "))
                result_all['Lpk1_tbz'], result_all['Lpk1_tby'], result_all['Lpk1_tbx'] = int(temp[3]),  int(temp[5]), int(temp[7])
                Lpk1_done = True
            elif "Lpk2 block" in line and not Lpk2_done:
                line = line.replace(",","")
                temp = list(line.split(" "))
                result_all['Lpk2_tbz'], result_all['Lpk2_tby'], result_all['Lpk2_tbx'] = int(temp[3]),  int(temp[5]), int(temp[7])
                Lpk2_done = True
            elif "Lpk3 block" in line and not Lpk3_done:
                line = line.replace(",","")
                temp = list(line.split(" "))
                result_all['Lpk3_tbz'], result_all['Lpk3_tby'], result_all['Lpk3_tbx'] = int(temp[3]),  int(temp[5]), int(temp[7])
                Lpk3_done = True
            elif "Ipk1 block" in line and not Ipk1_done:
                line = line.replace(",","")
                temp = list(line.split(" "))
                result_all['Ipk1_tbz'], result_all['Ipk1_tby'], result_all['Ipk1_tbx'] = int(temp[3]),  int(temp[5]), int(temp[7])
                Ipk1_done = True
            elif "Ipk2 block" in line and not Ipk2_done:
                line = line.replace(",","")
                temp = list(line.split(" "))
                result_all['Ipk2_tbz'], result_all['Ipk2_tby'], result_all['Ipk2_tbx'] = int(temp[3]),  int(temp[5]), int(temp[7])
                Ipk2_done = True
            elif "Ipk3 block" in line and not Ipk3_done:
                line = line.replace(",","")
                temp = list(line.split(" "))
                result_all['Ipk3_tbz'], result_all['Ipk3_tby'], result_all['Ipk3_tbx'] = int(temp[3]),  int(temp[5]), int(temp[7])
                Ipk3_done = True
            elif "LinearQuantizer block" in line:
                line = line.replace(",","")
                temp = list(line.split(" "))
                result_all['Quan_tbz'], result_all['Quan_tby'], result_all['Quan_tbx'] = int(temp[3]),  int(temp[5]), int(temp[7])
            elif "mean" in line and "var" in line:
                line = line.replace(",","")
                temp = list(line.split(" "))
                result_all['mean'], result_all['var'] = float(temp[1]),  float(temp[3])
            elif "Real quantization is" in line:
                temp = list(line.split(" "))
                result_all['Real_quan'] = int(temp[4])
            elif "EncodeFixedLen block" in line:
                temp = list(line.split(" "))
                result_all['Encode_tbx'] = int(temp[3])
            elif "Deflate block" in line:
                temp = list(line.split(" "))
                result_all['Deflate_tbx'] = int(temp[3])
            elif "Huffman block size" in line:
                line = line.replace("\n","")
                temp = list(line.split(" "))
                result_all['PART_SIZE'] = int(temp[3])
            elif "Outlier ratio:" in line:
                temp = list(line.split(" "))
                dividend_str, divisor_str = temp[2].split('/')
                result_all['Outlier_point'] = int(dividend_str)
            line = rf.readline()
    if "mean" not in result_all:
        result_all["mean"] = 1
        result_all['var'] = 0
        #print(result_all)
    return(result_all)

#informationparse("GCLDLWP_13_900_1800_3.txt")
