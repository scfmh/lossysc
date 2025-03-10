#ifndef SZ3_IMPL_SZ_HPP
#define SZ3_IMPL_SZ_HPP

#include "SZ3/def.hpp"
#include "SZ3/api/impl/SZDispatcher.hpp"
#include "SZ3/api/impl/SZImplOMP.hpp"
#include <cmath>

template<class T, SZ::uint N>
char *SZ_compress_impl(SZ::Config &conf, T *data, size_t &outSize) {
#ifndef _OPENMP
    conf.openmp=false;
#endif
    if (conf.openmp) {
        return SZ_compress_OMP<T, N>(conf, data, outSize);
    } else {
        return SZ_compress_dispatcher<T, N>(conf, data, outSize);
    }
}


template<class T, SZ::uint N>
void SZ_decompress_impl(SZ::Config &conf, char *cmpData, size_t cmpSize, T *decData) {
#ifndef _OPENMP
    conf.openmp=false;
#endif
    if (conf.openmp) {
        SZ_decompress_OMP<T, N>(conf, cmpData, cmpSize, decData);
    } else {
        SZ_decompress_dispatcher<T, N>(conf, cmpData, cmpSize, decData);
    }
}

#endif