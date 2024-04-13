from general import remove_symbols
#desired_order = ['Lwpk', 'GpkReo3D', 'Lpk1Reo3D', 'Lpk2Reo3D','Lpk3Reo3D', 'Ipk1Reo3D','Ipk2Reo3D', 'Ipk3Reo3D', 'Quantization','Outlier','Histogram', 'FillArray', 'Sort',
                 #'FirstNonzeroIndex','GenerateCL','GenerateCW', 'ReverseArray', 'ReorderByIndex', 'Encode', 'Deflate', 'Condense']
desired_order = ['Lwpk', 'GpkReo3D', 'Lpk1Reo3D', 'Lpk2Reo3D','Lpk3Reo3D', 'Ipk1Reo3D','Ipk2Reo3D', 'Ipk3Reo3D', 'Quantization','Histogram', 'FillArray', 'Sort',
                 'FirstNonzeroIndex','GenerateCL','GenerateCW', 'ReverseArray', 'ReorderByIndex', 'Encode', 'Deflate', 'Condense','HtoD', 'memset', 'DtoH' ]
result_time = {}
def timeparse_nvprof(inputfile):
    Deflate_t = 0
    GenerateCL_t = 0
    Ipk1Reo3D_t = 0
    Ipk2Reo3D_t = 0
    Ipk3Reo3D_t = 0
    Sort_t = 0
    GenerateCW_t = 0
    Lpk1Reo3D_t = 0
    Lpk2Reo3D_t = 0
    Lpk3Reo3D_t = 0
    GpkReo3D_t = 0
    Histogram_t = 0
    Quantization_t = 0
    Quantization_status = 0
    Encode_t = 0
    Outlier_t = 0
    Lwpk_t = 0
    Lwpk_status = False
    Condense_t = 0
    ReverseArray_t = 0
    ReorderByIndex_t = 0
    GetFirstNonzeroIndex_t = 0
    FillArraySequence_t = 0
    HtoD_t = 0
    memset_t = 0
    DtoH_t = 0
    with open(inputfile, 'r') as rf:
        line = rf.readline()
        while line:
            if "DeflateFunctor" in line:
                line = remove_symbols(line)
                print(line)
                Deflate_t += float(list(line.split(" "))[2])
            elif "GenerateCLFunctor" in line:
                line = remove_symbols(line)
                GenerateCL_t += float(list(line.split(" "))[2])
            elif "Ipk1Reo3DFunctor" in line:
                line = remove_symbols(line)
                Ipk1Reo3D_t += float(list(line.split(" "))[2]) 
            elif "Ipk2Reo3DFunctor" in line:
                line = remove_symbols(line)
                Ipk2Reo3D_t += float(list(line.split(" "))[2]) 
            elif "Ipk3Reo3DFunctor" in line:
                line = remove_symbols(line)
                Ipk3Reo3D_t += float(list(line.split(" "))[2]) 
            elif "DeviceRadixSortOnesweepKernel" in line:
                line = remove_symbols(line)
                Sort_t = float(list(line.split(" "))[2])
            elif "GenerateCWFunctor" in line:
                line = remove_symbols(line)
                GenerateCW_t = float(list(line.split(" "))[2])
            elif "Lpk1Reo3DFunctor" in line:
                line = remove_symbols(line)
                Lpk1Reo3D_t += float(list(line.split(" "))[2]) 
            elif "Lpk2Reo3DFunctor" in line:
                line = remove_symbols(line)
                Lpk2Reo3D_t += float(list(line.split(" "))[2]) 
            elif "Lpk3Reo3DFunctor" in line:
                line = remove_symbols(line)
                Lpk3Reo3D_t += float(list(line.split(" "))[2]) 
            elif "GpkReo3DFunctor" in line:
                line = remove_symbols(line)
                GpkReo3D_t += float(list(line.split(" "))[2])
            elif "HistogramFunctor" in line:
                line = remove_symbols(line)
                Histogram_t += float(list(line.split(" "))[2])
            elif "LevelwiseLinearQuantizerNDFunctor" in line and Quantization_status == False:
                line = remove_symbols(line)
                Quantization_t += float(list(line.split(" "))[2])
                Quantization_status = True
            #elif "OutlierRestoreFunctor" in line:
                #line = remove_symbols(line)
                #Outlier_t += float(list(line.split(" "))[2])
            elif "EncodeFixedLenFunctor" in line:
                line = remove_symbols(line)
                Encode_t += float(list(line.split(" "))[2])
            elif "LwpkReoFunctor" in line:
                line = remove_symbols(line)
                Lwpk_t += float(list(line.split(" "))[2]) 
                # Lwpk_status = True
            elif "CondenseFunctor" in line:
                line = remove_symbols(line)
                Condense_t += float(list(line.split(" "))[2])
            elif "ReverseArrayFunctor" in line:
                line = remove_symbols(line)
                ReverseArray_t += float(list(line.split(" "))[2])
            elif "ReorderByIndexFunctor" in line:
                line = remove_symbols(line)
                ReorderByIndex_t += float(list(line.split(" "))[2])
            elif "GetFirstNonzeroIndexFunctor" in line:
                line = remove_symbols(line)
                GetFirstNonzeroIndex_t += float(list(line.split(" "))[2])
            elif "FillArraySequenceFunctor" in line:
                line = remove_symbols(line)
                FillArraySequence_t += float(list(line.split(" "))[2])
            elif "CUDA memcpy HtoD" in line:
                line = remove_symbols(line)
                HtoD_t += float(list(line.split(" "))[2])
            elif "CUDA memset" in line:
                line = remove_symbols(line)
                memset_t += float(list(line.split(" "))[2])
            elif "CUDA memcpy DtoH" in line:
                line = remove_symbols(line)
                DtoH_t += float(list(line.split(" "))[2])
            line = rf.readline()
            result_time['Deflate'] = Deflate_t
            result_time['GpkReo3D'] = GpkReo3D_t
            result_time['Lpk1Reo3D'] = Lpk1Reo3D_t
            result_time['Lpk2Reo3D'] = Lpk2Reo3D_t
            result_time['Lpk3Reo3D'] = Lpk3Reo3D_t
            result_time['Ipk1Reo3D'] = Ipk1Reo3D_t
            result_time['Ipk2Reo3D'] = Ipk2Reo3D_t
            result_time['Ipk3Reo3D'] = Ipk3Reo3D_t
            result_time['Quantization'] = Quantization_t
            result_time['Histogram'] = Histogram_t
            result_time['Sort'] = Sort_t
            result_time['GenerateCL'] = GenerateCL_t
            result_time['GenerateCW'] = GenerateCW_t
            result_time['Encode'] = Encode_t
            result_time['Lwpk'] = Lwpk_t
            result_time['FillArray'] = FillArraySequence_t
            result_time['FirstNonzeroIndex'] = GetFirstNonzeroIndex_t
            result_time['ReverseArray'] = ReverseArray_t
            result_time['Condense'] = Condense_t
            result_time['ReorderByIndex'] = ReorderByIndex_t
            result_time['HtoD'] = HtoD_t
            result_time['memset'] = memset_t
            result_time['DtoH'] = DtoH_t
            #result_time['Outlier'] = Outlier_t
            sorted_data = {key: result_time.get(key, 0) for key in desired_order}
        # print(GpkReo3D_t, Lpk1Reo3D_t, Lpk2Reo3D_t, Lpk3Reo3D_t, Ipk1Reo3D_t, Ipk2Reo3D_t, Ipk3Reo3D_t, Quantization_t, \
        #    Outlier_t, Histogram_t, Sort_t, GenerateCL_t, GenerateCW_t, Encode_t, Deflate_t)
        #print(sorted_data)
    return sorted_data
    # return sorted_data, GpkReo3D_t, Lpk1Reo3D_t, Lpk2Reo3D_t, Lpk3Reo3D_t, Ipk1Reo3D_t, Ipk2Reo3D_t, Ipk3Reo3D_t, Quantization_t, \
    #        Outlier_t, Histogram_t, Sort_t, GenerateCL_t, GenerateCW_t, Encode_t, Deflate_t


#timeparse_nvprof("GCLDLWP_13_900_1800_3.txt")
