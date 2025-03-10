/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#ifndef MGARD_X_MEMORY_MANAGEMENT_HPP
#define MGARD_X_MEMORY_MANAGEMENT_HPP

#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <numeric>
#include <sstream>   // std::stringstream
#include <stdexcept> // std::runtime_error
#include <string>
#include <utility> // std::pair
#include <vector>
#include <type_traits>
// #include "MemoryManagement.h"

#define ANSI_RED "\x1b[31m"
#define ANSI_GREEN "\x1b[32m"
#define ANSI_RESET "\x1b[0m"

namespace mgard_x {

template<typename T1, typename T2>
auto add(T1 a, T2 b) {
    if constexpr (std::is_same_v<long long int, T1> || std::is_same_v<long long int, T2>) {
        using ResultType = long long int;
        std::cout << "sss" <<std::endl;
        return static_cast<ResultType>(a) + static_cast<ResultType>(b);
    } else if constexpr (std::is_same_v<int, T1> && std::is_same_v<int, T2>) {
        using ResultType = int;
        return static_cast<ResultType>(a) + static_cast<ResultType>(b);
    } else {
        using ResultType = decltype(a + b);
        return static_cast<ResultType>(a) + static_cast<ResultType>(b);
    }
}
template<typename T>
double calculateVariance(const T* values, const T* occurrences, int mean, int size) {
    double sumOfSquaredDifferences = 0;
    auto totalOccurrences = 0;

    for (int i = 0; i < size; ++i) {
        double difference = values[i] - mean;
        sumOfSquaredDifferences = add(sumOfSquaredDifferences ,difference * difference * occurrences[i]);
        totalOccurrences = add(totalOccurrences, occurrences[i]);
    }

    return sumOfSquaredDifferences / totalOccurrences;
}

template <typename SubArrayType>
int Quantization_mean(std::string name_key, std::string name_value, SubArrayType subArray_key, SubArrayType subArray_value, int first) {
  SIZE nrow = 1;
  SIZE ncol = 1;
  SIZE nfib = 1;

  DIM D = SubArrayType::NumDims;

  nfib = subArray_key.shape(D - 1);
  if (SubArrayType::NumDims >= 2)
    ncol = subArray_key.shape(D - 2);
  if (SubArrayType::NumDims >= 3)
    nrow = subArray_key.shape(D - 3);

  using T = typename SubArrayType::DataType;
  using DeviceType = typename SubArrayType::DevType;
  auto total =0;
  auto count = 0;
  int mean;
  //std::cout << "SubArray value: " << name_value << "(" << nrow << " * " << ncol << " * "
            //<< nfib << ") sizeof(T) = " << sizeof(T) << std::endl;
  T *v_key = new T[nrow * ncol * nfib];
  T *v_value = new T[nrow * ncol * nfib];
  int non_size = 8192 - first;
  T *real_key = new T[non_size];
  T *real_value = new T[non_size];
  DeviceRuntime<DeviceType>::SyncQueue(0);
  for (SIZE i = 0; i < nrow; i++) {
    MemoryManager<DeviceType>::CopyND(
        v_key + ncol * nfib * i, nfib,
        subArray_key.data() + subArray_key.lddv1() * subArray_key.lddv2() * i,
        subArray_key.lddv1(), nfib, ncol, 0);
   MemoryManager<DeviceType>::CopyND(
        v_value + ncol * nfib * i, nfib,
        subArray_value.data() + subArray_value.lddv1() * subArray_value.lddv2() * i,
        subArray_value.lddv1(), nfib, ncol, 0);
  }
  DeviceRuntime<DeviceType>::SyncQueue(0);
  for (int i = 0; i < nrow; i++) {
    //printf("[i = %d]\n", i);
    for (int j = 0; j < ncol; j++) {
      for (int k = 0; k < nfib; k++) {
        if (nfib * ncol * i + nfib * j + k >= first) {
              total = add(total, v_value[nfib * ncol * i + nfib * j + k] * v_key[nfib * ncol * i + nfib * j + k]);
              count = add(count, v_value[nfib * ncol * i + nfib * j + k]);
              real_key[nfib * ncol * i + nfib * j + k - first] = v_key[nfib * ncol * i + nfib * j + k];
              real_value[nfib * ncol * i + nfib * j + k - first] = v_key[nfib * ncol * i + nfib * j + k];
        } 
       }
      //std::cout << std::endl;
    }
    //std::cout << std::endl;
  }
  std::cout << "this is total " << total << "count " << count << std::endl;
  mean = int(total/count)+1;
  double std_value = calculateVariance(real_key, real_value, mean, non_size); 
  std::cout << "this is mean " << mean << "std " << std_value;
  std::cout << std::endl;
  delete[] v_key;
  delete[] v_value;
  return mean;
}


template <typename SubArrayType>
void PrintSubarray(std::string name, SubArrayType subArray) {
  // Handle<1, float> tmp_handle;

  SIZE nrow = 1;
  SIZE ncol = 1;
  SIZE nfib = 1;

  DIM D = SubArrayType::NumDims;

  nfib = subArray.shape(D - 1);
  if (SubArrayType::NumDims >= 2)
    ncol = subArray.shape(D - 2);
  if (SubArrayType::NumDims >= 3)
    nrow = subArray.shape(D - 3);

  using T = typename SubArrayType::DataType;
  using DeviceType = typename SubArrayType::DevType;

  std::cout << "SubArray: " << name << "(" << nrow << " * " << ncol << " * "
            << nfib << ") sizeof(T) = " << sizeof(T) << std::endl;
  // std::cout << name << "\n";

  T *v = new T[nrow * ncol * nfib];
  // cudaMemcpy3DAsyncHelper(tmp_handle, v, nfib * sizeof(T), nfib * sizeof(T),
  //                         ncol, subArray.data(), subArray.lddv1 * sizeof(T),
  //                         nfib * sizeof(T), subArray.lddv2, nfib * sizeof(T),
  //                         ncol, nrow, D2H, 0);
  DeviceRuntime<DeviceType>::SyncQueue(0);
  for (SIZE i = 0; i < nrow; i++) {
    MemoryManager<DeviceType>::CopyND(
        v + ncol * nfib * i, nfib,
        subArray.data() + subArray.lddv1() * subArray.lddv2() * i,
        subArray.lddv1(), nfib, ncol, 0);
  }
  DeviceRuntime<DeviceType>::SyncQueue(0);
  // tmp_handle.sync(0);

  for (int i = 0; i < nrow; i++) {
    printf("[i = %d]\n", i);
    for (int j = 0; j < ncol; j++) {
      for (int k = 0; k < nfib; k++) {
        // std::cout << "[ " << j << ", " << k <<" ]: ";
        if (std::is_same<T, std::uint8_t>::value) {
          std::cout << std::setw(8)
                    << (unsigned int)v[nfib * ncol * i + nfib * j + k] << " ";
        } else {
          std::cout << std::setw(8) << std::setprecision(6) << std::fixed
                    << v[nfib * ncol * i + nfib * j + k] << " ";
        }
        // std::cout << "\n";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  delete[] v;
}

template <typename SubArrayType>
void CompareSubarray(std::string name, SubArrayType subArray1,
                     SubArrayType subArray2) {
  // Handle<1, float> tmp_handle;

  DIM D = SubArrayType::NumDims;

  SIZE nrow = 1;
  SIZE ncol = 1;
  SIZE nfib = 1;

  nfib = subArray1.shape(D - 1);
  if (SubArrayType::NumDims >= 2)
    ncol = subArray1.shape(D - 2);
  if (SubArrayType::NumDims >= 3)
    nrow = subArray1.shape(D - 3);

  if (subArray1.shape(D - 1) != subArray2.shape(D - 1) ||
      subArray1.shape(D - 2) != subArray2.shape(D - 2) ||
      subArray1.shape(D - 3) != subArray2.shape(D - 3)) {
    std::cout << log::log_err << "CompareSubarray: shape mismatch!\n";
    exit(-1);
  }

  using T = typename SubArrayType::DataType;
  using DeviceType = typename SubArrayType::DevType;
  std::cout << "SubArray: " << name << "(" << nrow << " * " << ncol << " * "
            << nfib << ") sizeof(T) = " << sizeof(T) << std::endl;

  T *v1 = new T[nrow * ncol * nfib];
  T *v2 = new T[nrow * ncol * nfib];

  for (SIZE i = 0; i < nrow; i++) {
    MemoryManager<DeviceType>::CopyND(
        v1 + ncol * nfib * i, nfib,
        subArray1.data() + subArray1.lddv1() * subArray1.lddv2() * i,
        subArray1.lddv1(), nfib, ncol, 0);
    MemoryManager<DeviceType>::CopyND(
        v2 + ncol * nfib * i, nfib,
        subArray2.data() + subArray2.lddv1() * subArray2.lddv2() * i,
        subArray2.lddv1(), nfib, ncol, 0);
  }
  DeviceRuntime<DeviceType>::SyncQueue(0);

  T max_error = 0;
  bool pass = true;
  for (int i = 0; i < nrow; i++) {
    printf("[i = %d]\n", i);
    for (int j = 0; j < ncol; j++) {
      for (int k = 0; k < nfib; k++) {
        T a = v1[nfib * ncol * i + nfib * j + k];
        T b = v2[nfib * ncol * i + nfib * j + k];
        max_error = std::max(max_error, fabs(a - b));
        if (fabs(a - b) > 1e-5 * fabs(a)) {
          std::cout << ANSI_RED;
          pass = false;
        } else {
          std::cout << ANSI_GREEN;
        }
        if (std::is_same<T, std::uint8_t>::value) {
          std::cout << std::setw(8)
                    << (unsigned int)v2[nfib * ncol * i + nfib * j + k] << ", ";
        } else {
          std::cout << std::setw(8) << std::setprecision(6) << std::fixed
                    << v2[nfib * ncol * i + nfib * j + k] << ", ";
        }
        std::cout << ANSI_RESET;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  printf("Check: %d, max error: %f\n", pass, max_error);
  delete[] v1;
  delete[] v2;
}

template <typename SubArrayType1, typename SubArrayType2>
void CompareSubarray(std::string name, SubArrayType1 subArray1,
                     SubArrayType2 subArray2, bool print,
                     double error_thresold) {
  // Handle<1, float> tmp_handle;

  DIM D = SubArrayType1::NumDims;

  SIZE nrow = 1;
  SIZE ncol = 1;
  SIZE nfib = 1;

  nfib = subArray1.shape(D - 1);
  if (SubArrayType1::NumDims >= 2)
    ncol = subArray1.shape(D - 2);
  if (SubArrayType1::NumDims >= 3)
    nrow = subArray1.shape(D - 3);

  if (subArray1.shape(D - 1) != subArray2.shape[0] ||
      (SubArrayType1::NumDims >= 2 &&
       subArray1.shape(D - 2) != subArray2.shape[1]) ||
      (SubArrayType1::NumDims >= 3 &&
       subArray1.shape(D - 3) != subArray2.shape[2])) {
    std::cout << log::log_err << "CompareSubarray: shape mismatch!\n";
    exit(-1);
  }

  using T = typename SubArrayType1::DataType;
  using DeviceType = typename SubArrayType1::DevType;
  // std::cout << "SubArray: " << name << "(" << nrow << " * " << ncol << " * "
  // << nfib << ") sizeof(T) = "  <<sizeof(T) << std::endl;
  std::cout << name << "\n";

  T *v1 = new T[nrow * ncol * nfib];
  T *v2 = new T[nrow * ncol * nfib];
  for (SIZE i = 0; i < nrow; i++) {
    MemoryManager<DeviceType>::CopyND(
        v1 + ncol * nfib * i, nfib,
        subArray1.data() + subArray1.lddv1() * subArray1.lddv2() * i,
        subArray1.lddv1(), nfib, ncol, 0);
    MemoryManager<DeviceType>::CopyND(v2 + ncol * nfib * i, nfib,
                                      subArray2.data() +
                                          subArray2.lddv1 * subArray2.lddv2 * i,
                                      subArray2.lddv1, nfib, ncol, 0);
  }
  DeviceRuntime<DeviceType>::SyncQueue(0);

  T max_error = 0;
  bool pass = true;
  for (int i = 0; i < nrow; i++) {
    if (print)
      printf("[i = %d]\n", i);
    for (int j = 0; j < ncol; j++) {
      for (int k = 0; k < nfib; k++) {
        T a = v1[nfib * ncol * i + nfib * j + k];
        T b = v2[nfib * ncol * i + nfib * j + k];
        max_error = std::max(max_error, fabs(a - b));
        if (fabs(a - b) > error_thresold * fabs(a)) {
          if (print)
            std::cout << ANSI_RED;
          pass = false;
        } else {
          if (print)
            std::cout << ANSI_GREEN;
        }
        if (std::is_same<T, std::uint8_t>::value) {
          if (print)
            std::cout << std::setw(8)
                      << (unsigned int)v2[nfib * ncol * i + nfib * j + k]
                      << ", ";
        } else {
          if (print)
            std::cout << std::setw(8) << std::setprecision(6) << std::fixed
                      << v2[nfib * ncol * i + nfib * j + k] << ", ";
        }
        if (print)
          std::cout << ANSI_RESET;
      }
      if (print)
        std::cout << std::endl;
    }
    if (print)
      std::cout << std::endl;
  }
  if (print)
    std::cout << std::endl;

  printf("Check: %d, max error: %f\n", pass, max_error);
  delete[] v1;
  delete[] v2;
}

template <typename SubArrayType>
void CompareSubarray4D(SubArrayType subArray1, SubArrayType subArray2) {
  if (SubArrayType::NumDims != 4) {
    std::cout << log::log_err
              << "CompareSubarray4D expects 4D subarray type.\n";
    exit(-1);
  }

  DIM D = SubArrayType::NumDims;

  if (subArray1.shape(D - 4) != subArray2.shape(D - 4)) {
    std::cout << log::log_err << "CompareSubarray4D mismatch 4D size.\n";
    exit(-1);
  }

  using T = typename SubArrayType::DataType;
  SIZE idx[4] = {0, 0, 0, 0};
  for (SIZE i = 0; i < subArray1.shape(0); i++) {
    idx[3] = i;
    SubArrayType temp1 = subArray1;
    SubArrayType temp2 = subArray2;
    // Adding offset to the 4th dim. (slowest)
    temp1.offset_dim(0, i);
    temp2.offset_dim(0, i);
    // Make 3D slice on the other three dims
    CompareSubarray("4D = " + std::to_string(i), temp1.Slice3D(1, 2, 3),
                    temp2.Slice3D(1, 2, 3));
  }
}

template <typename SubArrayType>
void PrintSubarray4D(std::string name, SubArrayType subArray1) {
  if (SubArrayType::NumDims != 4) {
    std::cout << log::log_err << "PrintSubarray4D expects 4D subarray type.\n";
    exit(-1);
  }

  DIM D = SubArrayType::NumDims;

  std::cout << name << "\n";
  using T = typename SubArrayType::DataType;
  SIZE idx[4] = {0, 0, 0, 0};
  for (SIZE i = 0; i < subArray1.shape(D - 4); i++) {
    idx[3] = i;
    SubArrayType temp1 = subArray1;
    temp1.offset_dim(3, i);
    PrintSubarray("i = " + std::to_string(i), temp1.Slice3D(0, 1, 2));
  }
}

// print 3D CPU
template <typename T>
void verify_matrix(SIZE nrow, SIZE ncol, SIZE nfib, T *v, SIZE ldv1, SIZE ldv2,
                   std::string file_prefix, bool store, bool verify) {
  std::string filename = file_prefix + ".dat";
  if (store) {
    std::ofstream myfile;
    myfile.open(filename, std::ios::out | std::ios::binary);
    if (!myfile) {
      printf("Error: cannot write file\n");
      return;
    }
    myfile.write((char *)v, nrow * ncol * nfib * sizeof(T));
    myfile.close();
    if (!myfile.good()) {
      printf("Error occurred at write time!\n");
      return;
    }
  }
  if (verify) {
    std::fstream fin;
    fin.open(filename, std::ios::in | std::ios::binary);
    if (!fin) {
      printf("Error: cannot read file\n");
      return;
    }
    T *v2 = new T[nrow * ncol * nfib];
    fin.read((char *)v2, nrow * ncol * nfib * sizeof(T));
    fin.close();
    if (!fin.good()) {
      printf("Error occurred at reading time!\n");
      return;
    }

    bool mismatch = false;
    for (int i = 0; i < nrow; i++) {
      for (int j = 0; j < ncol; j++) {
        for (int k = 0; k < nfib; k++) {
          if (v[get_idx(ldv1, ldv2, i, j, k)] !=
              v2[get_idx(nfib, ncol, i, j, k)]) {
            std::cout << filename << ": ";
            printf("Mismatch[%d %d %d] %f - %f\n", i, j, k,
                   v[get_idx(ldv1, ldv2, i, j, k)],
                   v2[get_idx(nfib, ncol, i, j, k)]);
            mismatch = true;
          }
        }
      }
    }

    delete[] v2;
    if (mismatch)
      exit(-1);
  }
}

// print 3D GPU
template <typename T>
void verify_matrix_cuda(SIZE nrow, SIZE ncol, SIZE nfib, T *dv, SIZE lddv1,
                        SIZE lddv2, SIZE sizex, std::string file_prefix,
                        bool store, bool verify) {
  // std::cout << std::setw(10);
  // std::cout << std::setprecision(2) << std::fixed;
  if (store || verify) {
    // Handle<3, float> *tmp_handle = new Handle<3, float>();
    int queue_idx = 0;

    T *v = new T[nrow * ncol * nfib];
    // cudaMemcpy3DAsyncHelper(*tmp_handle, v, nfib * sizeof(T), nfib *
    // sizeof(T),
    //                         ncol, dv, lddv1 * sizeof(T), sizex * sizeof(T),
    //                         lddv2, nfib * sizeof(T), ncol, nrow, D2H,
    //                         queue_idx);
    // MemoryManager<CUDA>::CopyND(v, nfib, dv, lddv1,
    //                           nfib, ncol * nrow, 0);
    // DeviceRuntime<CUDA>::SyncQueue(0);
    // tmp_handle->sync(queue_idx);
    verify_matrix(nrow, ncol, nfib, v, nfib, ncol, file_prefix, store, verify);
    delete[] v;
    // delete tmp_handle;
  }
}

template <typename T, typename DeviceType>
void CompareSubArrays(SubArray<1, T, DeviceType> array1,
                      SubArray<1, T, DeviceType> array2) {
  SIZE n = array1.shape[0];
  using Mem = MemoryManager<DeviceType>;
  T *q1 = new T[n];
  T *q2 = new T[n];
  Mem::Copy1D(q1, array1.data(), n, 0);
  Mem::Copy1D(q2, array2.data(), n, 0);
  DeviceRuntime<DeviceType>::SyncQueue(0);
  bool pass = true;
  for (int i = 0; i < n; i++) {
    if (q1[i] != q2[i]) {
      pass = false;
      std::cout << "diff at " << i << "(" << q1[i] << ", " << q2[i] << ")\n";
    }
  }
  printf("pass: %d\n", pass);
  delete[] q1;
  delete[] q2;
}

template <DIM D, typename T, typename DeviceType>
void DumpSubArray(std::string name, SubArray<D, T, DeviceType> subArray) {
  SIZE nrow = 1;
  SIZE ncol = 1;
  SIZE nfib = 1;

  nfib = subArray.shape(D - 1);
  if (D >= 2)
    ncol = subArray.shape(D - 2);
  if (D >= 3)
    nrow = subArray.shape(D - 3);

  T *v = new T[nrow * ncol * nfib];
  DeviceRuntime<DeviceType>::SyncQueue(0);
  MemoryManager<DeviceType>::CopyND(v, nfib, subArray.data(),
                                    subArray.ld(D - 1), nfib, ncol * nrow, 0);
  DeviceRuntime<DeviceType>::SyncQueue(0);
  std::fstream myfile;
  myfile.open(name, std::ios::out | std::ios::binary);
  if (!myfile) {
    printf("Error: cannot open file\n");
    return;
  }
  myfile.write((char *)v, nrow * ncol * nfib * sizeof(T));
  myfile.close();
  if (!myfile.good()) {
    printf("Error occurred at write time!\n");
    return;
  }
  delete[] v;
}

template <DIM D, typename T, typename DeviceType>
void LoadSubArray(std::string name, SubArray<D, T, DeviceType> subArray) {

  SIZE nrow = 1;
  SIZE ncol = 1;
  SIZE nfib = 1;

  nfib = subArray.shape(D - 1);
  if (D >= 2)
    ncol = subArray.shape(D - 2);
  if (D >= 3)
    nrow = subArray.shape(D - 3);

  T *v = new T[nrow * ncol * nfib];

  std::fstream myfile;
  myfile.open(name, std::ios::in | std::ios::binary);
  if (!myfile) {
    printf("Error: cannot open file\n");
    return;
  }
  myfile.read((char *)v, nrow * ncol * nfib * sizeof(T));
  myfile.close();
  if (!myfile.good()) {
    printf("Error occurred at read time!\n");
    return;
  }
  MemoryManager<DeviceType>::CopyND(subArray.data(), subArray.ld(D - 1), v,
                                    nfib, nfib, ncol * nrow, 0);
  delete[] v;
}

} // namespace mgard_x

#endif
