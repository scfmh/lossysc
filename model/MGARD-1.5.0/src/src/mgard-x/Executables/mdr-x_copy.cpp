/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include <chrono>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <sys/stat.h>
#include <sys/types.h>

#include "compress_x.hpp"
#include "mdr_x.hpp"
#include "mgard-x/Utilities/ErrorCalculator.h"
// #include "compress_cuda.hpp"

using namespace std::chrono;

void print_usage_message(std::string error) {
  if (error.compare("") != 0) {
    std::cout << mgard_x::log::log_err << error << std::endl;
  }
  printf("Options\n\
\t -z: compress data\n\
\t\t -i <path to data file to be compressed>\n\
\t\t -c <path to compressed file>\n\
\t\t -t <s|d>: data type (s: single; d:double)\n\
\t\t -n <ndim>: total number of dimensions\n\
\t\t\t [dim1]: slowest dimention\n\
\t\t\t [dim2]: 2nd slowest dimention\n\
\t\t\t  ...\n\
\t\t\t [dimN]: fastest dimention\n\
\t\t -u <path to coordinate file>\n\
\t\t -m <abs|rel>: error bound mode (abs: abolute; rel: relative)\n\
\t\t -e <error>: error bound\n\
\t\t -s <smoothness>: smoothness parameter\n\
\t\t -l choose lossless compressor (0:Huffman 1:Huffman+LZ4 2:Huffman+Zstd)\n\
\t\t -d <auto|serial|openmp|cuda|hip|sycl>: device type\n\
\t\t -v enable verbose (show timing and statistics)\n\
\n\
\t -x: decompress data\n\
\t\t -c <path to compressed file>\n\
\t\t -o <path to decompressed file>\n\
\t\t -d <auto|serial|cuda|hip>: device type\n\
\t\t -v enable verbose (show timing and statistics)\n");
  exit(0);
}

bool has_arg(int argc, char *argv[], std::string option) {
  for (int i = 0; i < argc; i++) {
    if (option.compare(std::string(argv[i])) == 0) {
      return true;
    }
  }
  return false;
}

bool require_arg(int argc, char *argv[], std::string option) {
  for (int i = 0; i < argc; i++) {
    if (option.compare(std::string(argv[i])) == 0) {
      return true;
    }
  }
  print_usage_message("missing option: " + option + ".");
  return false;
}

std::string get_arg(int argc, char *argv[], std::string option) {
  if (require_arg(argc, argv, option)) {
    for (int i = 0; i < argc; i++) {
      if (option.compare(std::string(argv[i])) == 0) {
        return std::string(argv[i + 1]);
      }
    }
  }
  return std::string("");
}

int get_arg_int(int argc, char *argv[], std::string option) {
  if (require_arg(argc, argv, option)) {
    std::string arg;
    int i;
    for (i = 0; i < argc; i++) {
      if (option.compare(std::string(argv[i])) == 0) {
        arg = std::string(argv[i + 1]);
      }
    }
    try {
      int d = std::stoi(arg);
      return d;
    } catch (std::invalid_argument const &e) {
      print_usage_message("illegal argument for option " + option + ".");
      return 0;
    }
  }
  return 0;
}

std::vector<mgard_x::SIZE> get_arg_dims(int argc, char *argv[],
                                        std::string option) {
  std::vector<mgard_x::SIZE> shape;
  if (require_arg(argc, argv, option)) {
    std::string arg;
    int arg_idx = 0, i;
    for (i = 0; i < argc; i++) {
      if (option.compare(std::string(argv[i])) == 0) {
        arg = std::string(argv[i + 1]);
        arg_idx = i + 1;
      }
    }
    try {
      int d = std::stoi(arg);
      for (int i = 0; i < d; i++) {
        shape.push_back(std::stoi(argv[arg_idx + 1 + i]));
      }
      return shape;
    } catch (std::invalid_argument const &e) {
      print_usage_message("illegal argument for option " + option + ".");
      return shape;
    }
  }
  return shape;
}

double get_arg_double(int argc, char *argv[], std::string option) {
  if (require_arg(argc, argv, option)) {
    std::string arg;
    int i;
    for (i = 0; i < argc; i++) {
      if (option.compare(std::string(argv[i])) == 0) {
        arg = std::string(argv[i + 1]);
      }
    }
    try {
      double d = std::stod(arg);
      return d;
    } catch (std::invalid_argument const &e) {
      print_usage_message("illegal argument for option " + option + ".");
    }
  }
  return 0;
}

template <typename T> void min_max(size_t n, T *in_buff) {
  T min = std::numeric_limits<T>::infinity();
  T max = 0;
  for (size_t i = 0; i < n; i++) {
    if (min > in_buff[i]) {
      min = in_buff[i];
    }
    if (max < in_buff[i]) {
      max = in_buff[i];
    }
  }
  printf("Min: %f, Max: %f\n", min, max);
}

template <typename T> size_t readfile(std::string input_file, T *&in_buff) {
  //std::cout << mgard_x::log::log_info << "Loading file: " << input_file <<
  // "\n";

  FILE *pFile;
  pFile = fopen(input_file.c_str(), "rb");
  if (pFile == NULL) {
    std::cout << mgard_x::log::log_err << "file open error!\n";
    exit(1);
  }
  fseek(pFile, 0, SEEK_END);
  size_t lSize = ftell(pFile);
  rewind(pFile);
  in_buff = (T *)malloc(lSize);
  lSize = fread(in_buff, 1, lSize, pFile);
  fclose(pFile);
  //std::cout << "this is the readfile size " <<input_file << " "  << lSize << std::endl;
  // min_max(lSize/sizeof(T), in_buff);
  return lSize;
}

template <typename T>
void readfile(std::string input_file, std::vector<T> &in_buff) {
   //std::cout << mgard_x::log::log_info << "Loading file: " << input_file <<
   //"\n";

  FILE *pFile;
  pFile = fopen(input_file.c_str(), "rb");
  if (pFile == NULL) {
    std::cout << mgard_x::log::log_err << "file open error!\n";
    exit(1);
  }
  fseek(pFile, 0, SEEK_END);
  size_t lSize = ftell(pFile);
  rewind(pFile);
  in_buff.resize(lSize / sizeof(T));
  lSize = fread(in_buff.data(), 1, lSize, pFile);
  fclose(pFile);
}

template <typename T>
std::vector<T *> readcoords(const char *input_file, mgard_x::DIM D,
                            std::vector<mgard_x::SIZE> shape) {
  std::cout << mgard_x::log::log_info
            << "Loading coordinate file: " << input_file << "\n";
  FILE *pFile;
  pFile = fopen(input_file, "rb");
  if (pFile == NULL) {
    std::cout << mgard_x::log::log_err << "coordinate file open error!\n";
    exit(1);
  }
  fseek(pFile, 0, SEEK_END);
  size_t lSize = ftell(pFile);
  size_t expected_size = 0;
  for (mgard_x::DIM d = 0; d < D; d++) {
    expected_size += sizeof(T) * shape[d];
  }
  if (lSize < expected_size) {
    std::cout << mgard_x::log::log_err << "coordinate file read error!\n";
    exit(-1);
  }
  rewind(pFile);
  std::vector<T *> coords(D);
  for (mgard_x::DIM d = 0; d < D; d++) {
    coords[d] = (T *)malloc(shape[d]);
    lSize = fread(coords[d], sizeof(T), shape[d], pFile);
  }
  fclose(pFile);
  return coords;
}

template <typename T>
void writefile(std::string output_file, T *out_buff, size_t num_bytes) {
  FILE *file = fopen(output_file.c_str(), "w");
  fwrite(out_buff, 1, num_bytes, file);
  fclose(file);
}

template <typename T>
void print_statistics(double s, enum mgard_x::error_bound_type mode,
                      std::vector<mgard_x::SIZE> shape, T *original_data,
                      T *decompressed_data, T tol, bool normalize_coordinates) {
  mgard_x::SIZE n = 1;
  for (mgard_x::DIM d = 0; d < shape.size(); d++)
    n *= shape[d];
  std::cout.setf(std::ios::scientific);
  std::cout.precision(12);
  T actual_error = 0.0;
  /*std::cout << "original data" <<std::endl;
  for (int i = 0; i < n; i++)
  std::cout << original_data[i] << " ";
  std::cout << "decompressed data" <<std::endl;
  for (int i = 0; i < n; i++)
  std::cout << decompressed_data[i] << " "; */
  //std::cout<< "this is the size" << sizeof(original_data) << std::endl;
  std::cout << std::scientific;
  if (s == std::numeric_limits<T>::infinity()) {
    actual_error =
        mgard_x::L_inf_error(n, original_data, decompressed_data, mode);
    if (mode == mgard_x::error_bound_type::ABS) {
      std::cout << mgard_x::log::log_info
                << "Absoluate L_inf error: " << actual_error << " ("
                << (actual_error < tol ? "\e[32mSatisified\e[0m"
                                       : "\e[31mNot Satisified\e[0m")
                << ")"
                << "\n";
    } else if (mode == mgard_x::error_bound_type::REL) {
      std::cout << mgard_x::log::log_info
                << "Relative L_inf error: " << actual_error << " ("
                << (actual_error < tol ? "\e[32mSatisified\e[0m"
                                       : "\e[31mNot Satisified\e[0m")
                << ")"
                << "\n";
    }
  } else {
    actual_error = mgard_x::L_2_error(shape, original_data, decompressed_data,
                                      mode, normalize_coordinates);
    if (mode == mgard_x::error_bound_type::ABS) {
      std::cout << mgard_x::log::log_info
                << "Absoluate L_2 error: " << actual_error << " ("
                << (actual_error < tol ? "\e[32mSatisified\e[0m"
                                       : "\e[31mNot Satisified\e[0m")
                << ")"
                << "\n";
    } else if (mode == mgard_x::error_bound_type::REL) {
      std::cout << mgard_x::log::log_info
                << "Relative L_2 error: " << actual_error << " ("
                << (actual_error < tol ? "\e[32mSatisified\e[0m"
                                       : "\e[31mNot Satisified\e[0m")
                << ")"
                << "\n";
    }
  }

  std::cout << mgard_x::log::log_info
            << "MSE: " << mgard_x::MSE(n, original_data, decompressed_data)
            << "\n";
  std::cout << std::defaultfloat;
  std::cout << mgard_x::log::log_info
            << "PSNR: " << mgard_x::PSNR(n, original_data, decompressed_data)
            << "\n";

  if (actual_error > tol)
    exit(-1);
}

void create_dir(std::string name) {
  struct stat st = {0};
  if (stat(name.c_str(), &st) == -1) {
    mkdir(name.c_str(), 0700);
  }
}

void write_mdr(mgard_x::MDR::RefactoredMetadata &refactored_metadata,
               mgard_x::MDR::RefactoredData &refactored_data,
               std::string output) {
  create_dir(output);
  std::vector<mgard_x::Byte> serialized_metadata =
      refactored_metadata.Serialize();
  writefile(output + "/header", refactored_metadata.header.data(),
            refactored_metadata.header.size());
  writefile(output + "/metadata", serialized_metadata.data(),
            serialized_metadata.size());
  for (int subdomain_id = 0; subdomain_id < refactored_metadata.metadata.size();
       subdomain_id++) {
    for (int level_idx = 0;
         level_idx <
         refactored_metadata.metadata[subdomain_id].level_sizes.size();
         level_idx++) {
      for (int bitplane_idx = 0;
           bitplane_idx < refactored_metadata.metadata[subdomain_id]
                              .level_sizes[level_idx]
                              .size();
           bitplane_idx++) {
        std::string filename = "component_" + std::to_string(subdomain_id) +
                               "_" + std::to_string(level_idx) + "_" +
                               std::to_string(bitplane_idx);
        writefile(output + "/" + filename,
                  refactored_data.data[subdomain_id][level_idx][bitplane_idx],
                  refactored_metadata.metadata[subdomain_id]
                      .level_sizes[level_idx][bitplane_idx]);
      }
    }
  }
}

void read_mdr(mgard_x::MDR::RefactoredMetadata &refactored_metadata,
              mgard_x::MDR::RefactoredData &refactored_data,
              std::string input) {
  //std::cout<< input << std::endl;
  readfile(input + "/header", refactored_metadata.header);
  std::vector<mgard_x::Byte> serialized_metadata;
  readfile(input + "/metadata", serialized_metadata);
  refactored_metadata.Deserialize(serialized_metadata);
  refactored_data.Initialize(refactored_metadata.num_subdomains);
  for (int subdomain_id = 0; subdomain_id < refactored_metadata.metadata.size();
       subdomain_id++) {
    refactored_data.data[subdomain_id].resize(
        refactored_metadata.metadata[subdomain_id].level_sizes.size());
    for (int level_idx = 0;
         level_idx <
         refactored_metadata.metadata[subdomain_id].level_sizes.size();
         level_idx++) {
      refactored_data.data[subdomain_id][level_idx].resize(
          refactored_metadata.metadata[subdomain_id]
              .level_sizes[level_idx]
              .size());
      for (int bitplane_idx = 0;
           bitplane_idx < refactored_metadata.metadata[subdomain_id]
                              .level_sizes[level_idx]
                              .size();
           bitplane_idx++) {
        std::string filename = "component_" + std::to_string(subdomain_id) +
                               "_" + std::to_string(level_idx) + "_" +
                               std::to_string(bitplane_idx);
        mgard_x::SIZE level_size = readfile(
            input + "/" + filename,
            refactored_data.data[subdomain_id][level_idx][bitplane_idx]);
        if (level_size != refactored_metadata.metadata[subdomain_id]
                              .level_sizes[level_idx][bitplane_idx]) {
          std::cout << "mdr component size mismatch.";
          exit(-1);
        }
      }
    }
  }
}

int verbose_to_log_level(int verbose) {
  if (verbose == 0) {
    return mgard_x::log::ERR;
  } else if (verbose == 1) {
    return mgard_x::log::ERR | mgard_x::log::INFO;
  } else if (verbose == 2) {
    return mgard_x::log::ERR | mgard_x::log::TIME;
  } else if (verbose == 3) {
    return mgard_x::log::ERR | mgard_x::log::INFO | mgard_x::log::TIME;
  }
}

template <typename T>
int launch_refactor(mgard_x::DIM D, enum mgard_x::data_type dtype,
                    std::string input_file, std::string output_file,
                    std::vector<mgard_x::SIZE> shape, bool non_uniform,
                    const char *coords_file, int lossless,
                    int domain_decomposition,
                    enum mgard_x::device_type dev_type, int num_dev,
                    int verbose, bool prefetch,
                    mgard_x::SIZE max_memory_footprint) {

  mgard_x::Config config;
  config.log_level = verbose_to_log_level(verbose);
  config.decomposition = mgard_x::decomposition_type::MultiDim;
  if (domain_decomposition == 0) {
    config.domain_decomposition = mgard_x::domain_decomposition_type::MaxDim;
  } else {
    config.domain_decomposition = mgard_x::domain_decomposition_type::Block;
  }
  config.dev_type = dev_type;
  config.num_dev = num_dev;
  config.prefetch = prefetch;
  config.max_memory_footprint = max_memory_footprint;

  if (lossless == 0) {
    config.lossless = mgard_x::lossless_type::Huffman;
  } else if (lossless == 1) {
    config.lossless = mgard_x::lossless_type::Huffman_LZ4;
  } else if (lossless == 2) {
    config.lossless = mgard_x::lossless_type::Huffman_Zstd;
  } else if (lossless == 3) {
    config.lossless = mgard_x::lossless_type::CPU_Lossless;
  }

  size_t original_size = 1;
  for (mgard_x::DIM i = 0; i < D; i++)
    original_size *= shape[i];
  T *original_data;
  size_t in_size = 0;
  if (std::string(input_file).compare("random") == 0) {
    in_size = original_size * sizeof(T);
    original_data = new T[original_size];
    srand(7117);
    T c = 0;
    for (size_t i = 0; i < original_size; i++) {
      original_data[i] = rand() % 10 + 1;
    }
  } else {
    in_size = readfile(input_file, original_data);
    //std::cout << "this is correct input file;" << std::endl;
  }
  if (in_size != original_size * sizeof(T)) {
    std::cout << mgard_x::log::log_warn << "input file size mismatch "
              << in_size << " vs. " << original_size * sizeof(T) << "!\n";
  }

  mgard_x::MDR::RefactoredMetadata refactored_metadata;
  mgard_x::MDR::RefactoredData refactored_data;
  mgard_x::pin_memory(original_data, original_size * sizeof(T), config);
  std::vector<const mgard_x::Byte *> coords_byte;
  if (!non_uniform) {
    mgard_x::MDR::MDRefactor(D, dtype, shape, original_data,
                             refactored_metadata, refactored_data, config,
                             false);
  } else {
    std::vector<T *> coords;
    if (non_uniform) {
      coords = readcoords<T>(coords_file, D, shape);
    }
    for (auto &coord : coords) {
      coords_byte.push_back((const mgard_x::Byte *)coord);
    }
    mgard_x::MDR::MDRefactor(D, dtype, shape, original_data, coords_byte,
                             refactored_metadata, refactored_data, config,
                             false);
  }

  write_mdr(refactored_metadata, refactored_data, output_file);

  read_mdr(refactored_metadata, refactored_data, output_file);
  refactored_metadata.InitializeForReconstruction();
  mgard_x::MDR::MDRequest(refactored_metadata, config);
  //201667.424,
  //std::numeric_limits<double>::infinity(),
  //           mgard_x::error_bound_type::ABS, config);
  for (auto metadata : refactored_metadata.metadata) {
     metadata.PrintStatus();
   }

   mgard_x::MDR::ReconstructedData reconstructed_data;

   mgard_x::MDR::MDReconstruct(refactored_metadata,refactored_data,
                               reconstructed_data, config, false);
  print_statistics<double>(0, mgard_x::error_bound_type::ABS, shape, (double *)original_data,(double *)reconstructed_data.data[0], 0.0001,config.normalize_coordinates);
  mgard_x::unpin_memory(original_data, config);
  delete[](T *) original_data;

  return 0;
}

int launch_reconstruct(std::string input_file, std::string output_file,
                       std::string original_file, enum mgard_x::data_type dtype,
                       std::vector<mgard_x::SIZE> shape, double tol, double s,
                       enum mgard_x::error_bound_type mode,
                       bool adaptive_resolution,
                       enum mgard_x::device_type dev_type, int num_dev,
                       int verbose, bool prefetch) {

  mgard_x::Config config;
  config.log_level = verbose_to_log_level(verbose);
  config.dev_type = dev_type;
  config.num_dev;
  config.prefetch = prefetch;
  config.mdr_adaptive_resolution = adaptive_resolution;
  // config.collect_uncertainty = true;

  mgard_x::Byte *original_data;
  size_t in_size = 0;
  if (original_file.compare("none") != 0 || !config.mdr_adaptive_resolution) {
    if (original_file.compare("random") == 0) {
      size_t original_size = 1;
      for (mgard_x::DIM i = 0; i < shape.size(); i++)
        original_size *= shape[i];
      if (dtype == mgard_x::data_type::Float) {
        in_size = original_size * sizeof(float);
        original_data = (mgard_x::Byte *)new float[original_size];
        srand(7117);
        for (size_t i = 0; i < original_size; i++) {
          ((float *)original_data)[i] = rand() % 10 + 1;
        }
      } else if (dtype == mgard_x::data_type::Double) {
        in_size = original_size * sizeof(double);
        original_data = (mgard_x::Byte *)new double[original_size];
        srand(7117);
        for (size_t i = 0; i < original_size; i++) {
          ((double *)original_data)[i] = rand() % 10 + 1;
        }
      }

    } else {
      in_size = readfile(original_file, original_data);
    }
  }

  mgard_x::MDR::RefactoredMetadata refactored_metadata;
  mgard_x::MDR::RefactoredData refactored_data;
  read_mdr(refactored_metadata, refactored_data, input_file);
  //std::cout << "this is the total subdomains" << refactored_data.num_subdomains << std::endl;
  refactored_metadata.InitializeForReconstruction();


  //std::cout << "this is the error in the mdr-x.cpp, tol is " << tol << " s is " << s << std::endl; 
  for (auto &metadata : refactored_metadata.metadata) {
    metadata.requested_tol = tol;
    metadata.requested_s = s;
    //std::cout << "this is loop the error in the mdr-x.cpp, tol is " << metadata.requested_tol << " s is " << metadata.requested_s << std::endl;
  }
  mgard_x::MDR::MDRequest(refactored_metadata,config);
  for (auto &metadata : refactored_metadata.metadata) {
    metadata.PrintStatus();
    //metadata.PrintLevelSizes();
   
  }
  //refactored_metadata.DoneLoadingBitplans();

  mgard_x::MDR::ReconstructedData reconstructed_data;

  mgard_x::MDR::MDReconstruct(refactored_metadata, refactored_data,
                              reconstructed_data, config, false, original_data);
  int subdomain_id = 0;

  for (int subdomain_id = 0; subdomain_id < reconstructed_data.data.size();
       subdomain_id++) {
    //std::cout.setf(std::ios::scientific);
    //std::cout.precision(12);
    //std::cout << " look at reconstruct_data " << (double*)reconstructed_data.data[0][10] << std::endl;
    std::cout << "reconstructed_data " << subdomain_id << " : offset(";
    for (auto n : reconstructed_data.offset[subdomain_id]) {
      std::cout << n << " ";
    }
    std::cout << ") shape(";
    for (auto n : reconstructed_data.shape[subdomain_id]) {
      std::cout << n << " ";
    }
    std::cout << ")\n";
  }
  //std::cout << "ssss" << std::endl;
  
  writefile(output_file, reconstructed_data.data[0], shape[0]*shape[1]*shape[2]*8); 
  if (input_file.compare("none") != 0 && !config.mdr_adaptive_resolution) {
    if (dtype == mgard_x::data_type::Float) {
      print_statistics<float>(s, mode, shape, (float *)original_data,
                              (float *)reconstructed_data.data[0], tol,
                              config.normalize_coordinates);
    } else if (dtype == mgard_x::data_type::Double) {
      //std::cout<< reconstructed_data.shape[0][0] << " "<<reconstructed_data.shape[0][1] <<reconstructed_data.shape[0][2]<<" " << sizeof(reconstructed_data.data[0]) <<" " <<*(double *)reconstructed_data.data[0]<< std::endl; 
      print_statistics<double>(s, mode, shape, (double *)original_data,(double *)reconstructed_data.data[0], tol,config.normalize_coordinates);
      //print_statistics<double>(s, mode, shape, (double *)original_data,(double *)reconstructed_data.data[0], tol,true);
    }
  }
  /*std::cout << "ssss" << std::endl;
  for (int i = 0; i < shape[0]*shape[1]*shape[2]; i++)
  std ::cout <<reconstructed_data.data[0][i]<<std::endl;
  FILE * file = fopen(output_file.c_str(), "w");
  fwrite(original_data, 1,shape[0]*shape[1]*shape[2]*8, file);
  fclose(file);

  writefile(output_file, original_data, shape[0]*shape[1]*shape[2]*8);*/
  //writefile(output_file, reconstructed_data.data[0], shape[0]*shape[1]*shape[2]*8);
  //std::cout << "this is the length of the outputfile" << reconstructed_data.data[0].size()<<std::endl;
  return 0;
}

bool try_refactoring(int argc, char *argv[]) {
  if (!has_arg(argc, argv, "-z"))
    return false;
  std::cout << mgard_x::log::log_info << "mode: refactoring\n";
  std::string input_file = get_arg(argc, argv, "-i");
  std::string output_file = get_arg(argc, argv, "-c");

  std::cout << mgard_x::log::log_info << "original data: " << input_file
            << "\n";
  std::cout << mgard_x::log::log_info << "refactored data: " << output_file
            << "\n";

  enum mgard_x::data_type dtype;
  std::string dt = get_arg(argc, argv, "-t");
  if (dt.compare("s") == 0) {
    dtype = mgard_x::data_type::Float;
    std::cout << mgard_x::log::log_info << "data type: Single precision\n";
  } else if (dt.compare("d") == 0) {
    dtype = mgard_x::data_type::Double;
    std::cout << mgard_x::log::log_info << "data type: Double precision\n";
  } else
    print_usage_message("wrong data type.");

  mgard_x::DIM D = get_arg_int(argc, argv, "-n");
  std::vector<mgard_x::SIZE> shape = get_arg_dims(argc, argv, "-n");
  std::string shape_string = "shape (";
  for (mgard_x::DIM d = 0; d < shape.size(); d++)
    shape_string = shape_string + std::to_string(shape[d]) + " ";
  shape_string = shape_string + ")";

  bool non_uniform = false;
  std::string non_uniform_coords_file;
  if (has_arg(argc, argv, "-u")) {
    non_uniform = true;
    non_uniform_coords_file = get_arg(argc, argv, "-u");
    std::cout << mgard_x::log::log_info
              << "non-uniform coordinate file: " << non_uniform_coords_file
              << "\n";
  }

  int lossless_level = get_arg_int(argc, argv, "-l");
  if (lossless_level == 0) {
    std::cout << mgard_x::log::log_info << "lossless: Huffman\n";
  } else if (lossless_level == 1) {
    std::cout << mgard_x::log::log_info << "lossless: Huffman + LZ4\n";
  } else if (lossless_level == 2) {
    std::cout << mgard_x::log::log_info << "lossless: Huffman + Zstd\n";
  }

  enum mgard_x::device_type dev_type;
  std::string dev = get_arg(argc, argv, "-d");
  if (dev.compare("auto") == 0) {
    dev_type = mgard_x::device_type::AUTO;
    std::cout << mgard_x::log::log_info << "device type: AUTO\n";
  } else if (dev.compare("serial") == 0) {
    dev_type = mgard_x::device_type::SERIAL;
    std::cout << mgard_x::log::log_info << "device type: SERIAL\n";
  } else if (dev.compare("openmp") == 0) {
    dev_type = mgard_x::device_type::OPENMP;
    std::cout << mgard_x::log::log_info << "device type: OPENMP\n";
  } else if (dev.compare("cuda") == 0) {
    dev_type = mgard_x::device_type::CUDA;
    std::cout << mgard_x::log::log_info << "device type: CUDA\n";
  } else if (dev.compare("hip") == 0) {
    dev_type = mgard_x::device_type::HIP;
    std::cout << mgard_x::log::log_info << "device type: HIP\n";
  } else if (dev.compare("sycl") == 0) {
    dev_type = mgard_x::device_type::SYCL;
    std::cout << mgard_x::log::log_info << "device type: SYCL\n";
  } else {
    print_usage_message("wrong device type.");
  }

  int verbose = 0;
  if (has_arg(argc, argv, "-v")) {
    verbose = get_arg_int(argc, argv, "-v");
  }

  int num_dev = 1;
  if (has_arg(argc, argv, "-g")) {
    num_dev = get_arg_int(argc, argv, "-g");
  }

  int repeat = 1;
  if (has_arg(argc, argv, "-p")) {
    repeat = get_arg_int(argc, argv, "-p");
  }

  bool prefetch = true;
  if (has_arg(argc, argv, "-h")) {
    prefetch = get_arg_int(argc, argv, "-h") == 1 ? true : false;
  }

  mgard_x::SIZE max_memory_footprint =
      std::numeric_limits<mgard_x::SIZE>::max();
  if (has_arg(argc, argv, "-f")) {
    max_memory_footprint = (mgard_x::SIZE)get_arg_double(argc, argv, "-f");
  }

  int domain_decomposition = 0;
  if (has_arg(argc, argv, "-b")) {
    domain_decomposition = get_arg_int(argc, argv, "-b");
  }

  if (verbose)
    std::cout << mgard_x::log::log_info << "Verbose: enabled\n";
  for (int repeat_iter = 0; repeat_iter < repeat; repeat_iter++) {
    if (dtype == mgard_x::data_type::Double) {
      launch_refactor<double>(
          D, dtype, input_file.c_str(), output_file.c_str(), shape, non_uniform,
          non_uniform_coords_file.c_str(), lossless_level, domain_decomposition,
          dev_type, num_dev, verbose, prefetch, max_memory_footprint);
    } else if (dtype == mgard_x::data_type::Float) {
      launch_refactor<float>(
          D, dtype, input_file.c_str(), output_file.c_str(), shape, non_uniform,
          non_uniform_coords_file.c_str(), lossless_level, domain_decomposition,
          dev_type, num_dev, verbose, prefetch, max_memory_footprint);
    }
  }
  return true;
}

bool try_reconstruction(int argc, char *argv[]) {
  if (!has_arg(argc, argv, "-x"))
    return false;
  std::cout << mgard_x::log::log_info << "mode: reconstruction\n";
  std::string input_file = get_arg(argc, argv, "-c");
  std::string output_file = get_arg(argc, argv, "-o");
  std::cout << mgard_x::log::log_info << "refactored data: " << input_file
            << "\n";
  std::cout << mgard_x::log::log_info << "reconstructed data: " << output_file
            << "\n";
  std::string original_file = "none";
  enum mgard_x::data_type dtype;
  std::vector<mgard_x::SIZE> shape;
  if (has_arg(argc, argv, "-a")) {
    original_file = get_arg(argc, argv, "-a");
    std::string dt = get_arg(argc, argv, "-t");
    if (dt.compare("s") == 0) {
      dtype = mgard_x::data_type::Float;
      std::cout << mgard_x::log::log_info << "data type: Single precision\n";
    } else if (dt.compare("d") == 0) {
      dtype = mgard_x::data_type::Double;
      std::cout << mgard_x::log::log_info << "data type: Double precision\n";
    } else
      print_usage_message("wrong data type.");

    mgard_x::DIM D = get_arg_int(argc, argv, "-n");
    shape = get_arg_dims(argc, argv, "-n");
    std::string shape_string = "shape (";
    for (mgard_x::DIM d = 0; d < shape.size(); d++)
      shape_string = shape_string + std::to_string(shape[d]) + " ";
    shape_string = shape_string + ")";
    std::cout << mgard_x::log::log_info << "original data: " << original_file
              << "\n";
  }

  enum mgard_x::error_bound_type mode; // REL or ABS
  std::string em = get_arg(argc, argv, "-m");
  if (em.compare("rel") == 0) {
    mode = mgard_x::error_bound_type::REL;
    std::cout << mgard_x::log::log_info << "error bound mode: Relative\n";
  } else if (em.compare("abs") == 0) {
    mode = mgard_x::error_bound_type::ABS;
    std::cout << mgard_x::log::log_info << "error bound mode: Absolute\n";
  } else
    print_usage_message("wrong error bound mode.");

  double tol = get_arg_double(argc, argv, "-e");
  double s = get_arg_double(argc, argv, "-s");

  std::cout << std::scientific;
  std::cout << mgard_x::log::log_info << "error bound: " << tol << "\n";
  std::cout << std::defaultfloat;
  std::cout << mgard_x::log::log_info << "s: " << s << "\n";

  enum mgard_x::device_type dev_type;
  std::string dev = get_arg(argc, argv, "-d");
  if (dev.compare("auto") == 0) {
    dev_type = mgard_x::device_type::AUTO;
    std::cout << mgard_x::log::log_info << "device type: AUTO\n";
  } else if (dev.compare("serial") == 0) {
    dev_type = mgard_x::device_type::SERIAL;
    std::cout << mgard_x::log::log_info << "device type: SERIAL\n";
  } else if (dev.compare("openmp") == 0) {
    dev_type = mgard_x::device_type::OPENMP;
    std::cout << mgard_x::log::log_info << "device type: OPENMP\n";
  } else if (dev.compare("cuda") == 0) {
    dev_type = mgard_x::device_type::CUDA;
    std::cout << mgard_x::log::log_info << "device type: CUDA\n";
  } else if (dev.compare("hip") == 0) {
    dev_type = mgard_x::device_type::HIP;
    std::cout << mgard_x::log::log_info << "device type: HIP\n";
  } else if (dev.compare("sycl") == 0) {
    dev_type = mgard_x::device_type::HIP;
    std::cout << mgard_x::log::log_info << "device type: SYCL\n";
  } else {
    print_usage_message("wrong device type.");
  }

  int verbose = 0;
  if (has_arg(argc, argv, "-v")) {
    verbose = get_arg_int(argc, argv, "-v");
  }

  int num_dev = 1;
  if (has_arg(argc, argv, "-g")) {
    num_dev = get_arg_int(argc, argv, "-g");
  }

  int repeat = 1;
  if (has_arg(argc, argv, "-p")) {
    repeat = get_arg_int(argc, argv, "-p");
  }

  bool prefetch = true;
  if (has_arg(argc, argv, "-h")) {
    prefetch = get_arg_int(argc, argv, "-h") == 1 ? true : false;
  }

  bool adaptive_resolution = false;
  if (has_arg(argc, argv, "-q")) {
    adaptive_resolution = get_arg_int(argc, argv, "-q") == 1 ? true : false;
  }

  if (verbose)
    std::cout << mgard_x::log::log_info << "verbose: enabled.\n";
  for (int repeat_iter = 0; repeat_iter < repeat; repeat_iter++) {
    launch_reconstruct(input_file, output_file, original_file, dtype, shape,
                       tol, s, mode, adaptive_resolution, dev_type, num_dev,
                       verbose, prefetch);
  }
  return true;
}

int main(int argc, char *argv[]) {

  if (!try_refactoring(argc, argv) && !try_reconstruction(argc, argv)) {
    print_usage_message("");
  }
  return 0;
}
