#ifndef SZ3_SZ_HPP
#define SZ3_SZ_HPP


#include "SZ3/api/impl/SZImpl.hpp"
#include "SZ3/version.hpp"
#include <memory>

/**
 * API for compression
 * @tparam T source data type
 * @param conf compression configuration. Please update the config with 1). data dimension and shape and 2). desired settings.
 * @param data source data
 * @param outSize compressed data size in bytes
 * @return compressed data, remember to 'delete []' when the data is no longer needed.

The compression algorithms are:
ALGO_INTERP_LORENZO:
 The default algorithm in SZ3. It is the implementation of our ICDE'21 paper.
 The whole dataset will be compressed by interpolation or lorenzo predictor with auto-optimized settings.
ALGO_INTERP:
 The whole dataset will be compressed by interpolation predictor with default settings.
ALGO_LORENZO_REG:
 The whole dataset will be compressed by lorenzo and/or regression based predictors block by block with default settings.
 The four predictors ( 1st-order lorenzo, 2nd-order lorenzo, 1st-order regression, 2nd-order regression)
 can be enabled or disabled independently by conf settings (lorenzo, lorenzo2, regression, regression2).

Interpolation+lorenzo example:
SZ::Config conf(100, 200, 300); // 300 is the fastest dimension
conf.cmprAlgo = SZ::ALGO_INTERP_LORENZO;
conf.errorBoundMode = SZ::EB_ABS; // refer to def.hpp for all supported error bound mode
conf.absErrorBound = 1E-3; // absolute error bound 1e-3
char *compressedData = SZ_compress(conf, data, outSize);

Interpolation example:
SZ::Config conf(100, 200, 300); // 300 is the fastest dimension
conf.cmprAlgo = SZ::ALGO_INTERP;
conf.errorBoundMode = SZ::EB_REL; // refer to def.hpp for all supported error bound mode
conf.relErrorBound = 1E-3; // value-rang-based error bound 1e-3
char *compressedData = SZ_compress(conf, data, outSize);

Lorenzo/regression example :
SZ::Config conf(100, 200, 300); // 300 is the fastest dimension
conf.cmprAlgo = SZ::ALGO_LORENZO_REG;
conf.lorenzo = true; // only use 1st order lorenzo
conf.lorenzo2 = false;
conf.regression = false;
conf.regression2 = false;
conf.errorBoundMode = SZ::EB_ABS; // refer to def.hpp for all supported error bound mode
conf.absErrorBound = 1E-3; // absolute error bound 1e-3
char *compressedData = SZ_compress(conf, data, outSize);
 */
template<class T>
char *SZ_compress(SZ::Config &conf, T *data, size_t &outSize) {
    char *cmpData;
    if (conf.N == 1) {
        cmpData = SZ_compress_impl<T, 1>(conf, data, outSize);
    } else if (conf.N == 2) {
        cmpData = SZ_compress_impl<T, 2>(conf, data, outSize);
    } else if (conf.N == 3) {
        cmpData = SZ_compress_impl<T, 3>(conf, data, outSize);
    } else if (conf.N == 4) {
        cmpData = SZ_compress_impl<T, 4>(conf, data, outSize);
    } else {
        for (int i = 4; i < conf.N; i++) {
            conf.dims[3] *= conf.dims[i];
        }
        conf.dims.resize(4);
        conf.N = 4;
        cmpData = SZ_compress_impl<T, 4>(conf, data, outSize);
    }
    {
        //save config
        SZ::uchar *cmpDataPos = (SZ::uchar *) cmpData + outSize;
        conf.save(cmpDataPos);
        size_t newSize = (char *) cmpDataPos - cmpData;
        SZ::write(int(newSize - outSize), cmpDataPos);
        outSize = (char *) cmpDataPos - cmpData;
    }
    return cmpData;
}


/**
 * API for decompression
 * Similar with SZ_decompress(SZ::Config &conf, char *cmpData, size_t cmpSize)
 * The only difference is this one needs pre-allocated decData as input
 * @tparam T decompressed data type
 * @param conf compression configuration. Setting the correct config is NOT needed for decompression.
 * The correct config will be loaded from compressed data and returned.
 * @param cmpData compressed data
 * @param cmpSize compressed data size in bytes
 * @param decData pre-allocated memory space for decompressed data

 example:
 auto decompressedData = new float[100x200x300];
 SZ::Config conf;
 SZ_decompress(conf, char *cmpData, size_t cmpSize, decompressedData);

 */
template<class T>
void SZ_decompress(SZ::Config &conf, char *cmpData, size_t cmpSize, T *&decData) {
    {
        //load config
        int confSize;
        memcpy(&confSize, cmpData + (cmpSize - sizeof(int)), sizeof(int));
        SZ::uchar const *cmpDataPos = (SZ::uchar *) cmpData + (cmpSize - sizeof(int) - confSize);
        conf.load(cmpDataPos);
    }
    if (decData == nullptr) {
        decData = new T[conf.num];
    }
    if (conf.N == 1) {
        SZ_decompress_impl<T, 1>(conf, cmpData, cmpSize, decData);
    } else if (conf.N == 2) {
        SZ_decompress_impl<T, 2>(conf, cmpData, cmpSize, decData);
    } else if (conf.N == 3) {
        SZ_decompress_impl<T, 3>(conf, cmpData, cmpSize, decData);
    } else if (conf.N == 4) {
        SZ_decompress_impl<T, 4>(conf, cmpData, cmpSize, decData);
    } else {
        SZ_decompress_impl<T, 4>(conf, cmpData, cmpSize, decData);
    }
}

/**
 * API for decompression
 * @tparam T decompressed data type
 * @param conf compression configuration. Setting the correct config is NOT needed for decompression.
 * The correct config will be loaded from compressed data and returned.
 * @param cmpData compressed data
 * @param cmpSize compressed data size in bytes
 * @return decompressed data, remember to 'delete []' when the data is no longer needed.

 example:
 SZ::Config conf;
 float decompressedData = SZ_decompress(conf, char *cmpData, size_t cmpSize)
 */
template<class T>
T *SZ_decompress(SZ::Config &conf, char *cmpData, size_t cmpSize) {
    T *decData = nullptr;
    SZ_decompress<T>(conf, cmpData, cmpSize, decData);
    return decData;
}

#endif