
#include <chrono>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <cstring>
#include <limits>

#include "mgard/compress_x.hpp"
#include "mgard/mgard-x/Utilities/ErrorCalculator.h"

//#define MYADIOS_USE_MPI 1
#if MYADIOS_USE_MPI
#include <mpi.h>
#endif

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
\t\t -l choose lossless compressor (0:Huffman 1:Huffman+LZ4 3:Huffman+Zstd)\n\
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

template <typename T> void readfile(const char *input_file, T *in_buff, size_t read_size) {
  //std::cout << mgard_x::log::log_info << "Loading file: " << input_file << "\n";

  FILE *pFile;
  pFile = fopen(input_file, "rb");
  if (pFile == NULL) {
    std::cout << mgard_x::log::log_err << "file open error!\n";
    exit(1);
  }
  // fseek(pFile, 0, SEEK_END);
  // size_t lSize = ftell(pFile);
  // rewind(pFile);
  fread(in_buff, 1, read_size, pFile);
  fclose(pFile);
  // min_max(lSize/sizeof(T), in_buff);
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
void writefile(const char *output_file, size_t num_bytes, T *out_buff) {
  FILE *file = fopen(output_file, "w");
  fwrite(out_buff, 1, num_bytes, file);
  fclose(file);
}

// template <typename T>
void print_statistics(double s, enum mgard_x::error_bound_type mode,
                      std::vector<mgard_x::SIZE> shape, float *original_data,
                      float *decompressed_data, float tol, bool normalize_coordinates) {
  mgard_x::SIZE n = 1;
  for (mgard_x::DIM d = 0; d < shape.size(); d++)
    n *= shape[d];
  float actual_error = 0.0;
  std::cout << std::scientific;
  if (s == std::numeric_limits<float>::infinity()) {
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

  // if (actual_error > tol)
  //   exit(-1);
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
  } else {
    return mgard_x::log::ERR;
  }
}

// // template <typename T> 

struct Metrics {
    float relativeErrors;
    float rmse;
    float maxError;
    float avgDistance;
    float linearmaxError;
};

Metrics calculateMetrics(const float* original, const float* decompressed, size_t size, int k) {
    float relaError = 0.0;
    float rmseSum = 0.0;
    float maxError = 0.0;
    double avgOriginal = 0.0;
    double avgDecompressed = 0.0;
    float linearmaxError = 0.0;

    // Calculate dataRange
    float minOriginal = original[0];
    float maxOriginal = original[0];

    for (size_t i = 0; i < size; ++i) {
      minOriginal = std::min(minOriginal, original[i]);
      maxOriginal = std::max(maxOriginal, original[i]);
    }

    float dataRange = maxOriginal - minOriginal;

    for (size_t i = 0; i < size; ++i) {
        // relaError = std::max(std::abs(decompressed[i] - original[i]));


        float squaredError = std::pow(original[i] - decompressed[i], 2);
        rmseSum += squaredError;

        maxError = std::max(maxError, std::abs(original[i] - decompressed[i]));

        avgOriginal += original[i];
        avgDecompressed += decompressed[i];
    }

    float avgDistance = std::abs(avgOriginal / size - avgDecompressed / size);
    std::cout << "Max error " << maxError << " range " << dataRange << std::endl;

    Metrics metrics;
    metrics.relativeErrors = maxError/ dataRange;
    metrics.rmse = std::sqrt(rmseSum / size);
    metrics.maxError = maxError;
    metrics.avgDistance = avgDistance;
    metrics.linearmaxError = maxError * k;

    return metrics;
}


struct result {
	double write;
  double read;
  double comp;
  double decomp;
};

template <typename T>
result launch_compress(mgard_x::DIM D, enum mgard_x::data_type dtype,
                    const char *input_file, const char *output_file,
                    std::vector<mgard_x::SIZE> shape, bool non_uniform,
                    const char *coords_file, double tol, double s,
                    enum mgard_x::error_bound_type mode, int reorder,
                    int lossless, enum mgard_x::device_type dev_type,
                    int num_dev, int verbose, int num_sim_iterations,
                    bool use_compression, int accumulate_data, int compute_delay,
                    bool prefetch, mgard_x::SIZE max_memory_footprint, 
                    int k_parameter) {
   

  /*int comm_size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);*/
  int rank, comm_size;

  #if MYADIOS_USE_MPI 
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  #else
  rank = 0;
  comm_size = 1;
  #endif
  if (rank == 0) { 
  std::cout << ">>> use_compression: " << use_compression << "\n";
  std::cout << ">>> accumulate_data: " << accumulate_data << "\n";
  std::cout << ">>> prefetch: " << prefetch << "\n";
  std::cout << ">>> comm_size: " << comm_size << "\n";
  //std::cout << "this is one time launch_compress" << std::endl;
	}
  result res;

  mgard_x::Config config;
  config.log_level = verbose_to_log_level(verbose);
  config.decomposition = mgard_x::decomposition_type::MultiDim;
  config.dev_type = dev_type;
  config.num_dev = num_dev;
  config.reorder = reorder;
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
 
  float* original_data = new float[original_size];

  if (std::string(input_file).compare("random") == 0) {
    srand(7117);
    T c = 0;
    for (size_t i = 0; i < original_size*accumulate_data; i++) {
      original_data[i] = rand() % 10 + 1;
    }
  } else {
    readfile(input_file, original_data, original_size * sizeof(float));
  }
  original_size *= accumulate_data;
  shape[0] *= accumulate_data;
  if (rank == 0) {
    std::cout << mgard_x::log::log_info << "Data per rank: "
              << (double)original_size*sizeof(T)/1e9 << " GB\n";
  }
  
  //unsigned char * compressed_data;
  void *compressed_data = malloc(original_size * sizeof(float));
  size_t compressed_size = original_size * sizeof(float);
  //void *compressed_data = (void *)new T[original_size];
  //size_t compressed_size = original_size * sizeof(T);
  //mgard_x::pin_memory(original_data, original_size * sizeof(T), config);
  //mgard_x::pin_memory(compressed_data, original_size * sizeof(T), config);
  void *decompressed_data = malloc(original_size * sizeof(float));
  //mgard_x::pin_memory(decompressed_data, original_size * sizeof(T), config);
  // mgard_x::Hierarchy<D, dtype, dev_type> hierarchy(shape, coords, config);;
  // mgard_x::Compressor<D, dtype, dev_type> compressor(hierarchy, config);
  // mgard_x::Array<D, dtype, dev_type> InputArray(shape);
  // InputArray.load(original_data);
  mgard_x::Timer timer;
  // ********* Compression ********** //
    #if MYADIOS_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.start();
    mgard_x::compress(D, dtype, shape, tol, s, mode, original_data, compressed_data, compressed_size, config, true);
    std::cout << "----------" << std::endl;
    #if MYADIOS_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.end();
    double compression_time = timer.get();
    //std::cout << "rank " << rank << " compression time: " << compression_time << "\n";
    #if MYADIOS_USE_MPI
    double max_compression_time;
    MPI_Reduce(&compression_time, &max_compression_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
      res.comp = max_compression_time;
      std::cout << mgard_x::log::log_info << "Compression time: "
              << max_compression_time << "\n";
      std::cout << mgard_x::log::log_info << "Compression throughput: "
              << (double)original_size * sizeof(T) * comm_size/max_compression_time/1e9 << " GB/s.\n";
      std::cout << mgard_x::log::log_info << "Compression ratio: "
              << (double)original_size * sizeof(T) / (compressed_size) << "\n";
    }
    #else
    if (rank == 0) {
      res.comp = compression_time;
      std::cout << mgard_x::log::log_info << "Compression time: "
              << compression_time << "\n";
      std::cout << mgard_x::log::log_info << "Compression throughput: "
              << (double)original_size * sizeof(T) * comm_size/compression_time/1e9 << " GB/s.\n";
      std::cout << mgard_x::log::log_info << "Compression ratio: "
              << (double)original_size * sizeof(T) / (compressed_size) << "\n";
    }
    #endif
    timer.clear();
    #if MYADIOS_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.start();
    mgard_x::decompress(compressed_data, compressed_size, decompressed_data,
                        config, true);
    #if MYADIOS_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.end();
    double decompress_time = timer.get();
    #if MYADIOS_USE_MPI
    double max_decompress_time;
    //std::cout << "rank " << rank << " decompress time: " << decompress_time << "\n";
    MPI_Reduce(&decompress_time, &max_decompress_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
      res.decomp = max_decompress_time;
      std::cout << mgard_x::log::log_info << "Decompression time: "
              << max_decompress_time << "\n";
      std::cout << mgard_x::log::log_info << "Decompression throughput: "
              << (double)original_size * sizeof(T) * comm_size/max_decompress_time/1e9 << " GB/s.\n";
    }
    #else
    if (rank == 0) {
      res.decomp = decompress_time;
      std::cout << mgard_x::log::log_info << "Decompression time: "
              << decompress_time << "\n";
      std::cout << mgard_x::log::log_info << "Decompression throughput: "
              << (double)original_size * sizeof(T) * comm_size/decompress_time/1e9 << " GB/s.\n";
    }

    #endif
    timer.clear();
  #if MYADIOS_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif
  //mgard_x::unpin_memory(original_data, config);
  //mgard_x::unpin_memory(compressed_data, config);
  //mgard_x::unpin_memory(decompressed_data, config);
  /****************qoi a=verage***************/
  //float* float_decompressed_data = static_cast<float*>(decompressed_data);
  // for (int i = 0; i < original_size; i++)
  //  std::cout << " " << float_decompressed_data[i] << " ";
  //  for (int i = 0; i < original_size; i++)
  //  std::cout << " " << original_data[i] << " ";

  Metrics metrics = calculateMetrics(original_data, (float*)decompressed_data, original_size, k_parameter);
  std::cout << "Relative : " << metrics.relativeErrors << "\n";
  std::cout << "RMSE: " << metrics.rmse << "\n";
  std::cout << "Max Error: " << metrics.maxError << "\n";
  std::cout << "Avg Distance between Original and Decompressed: " << metrics.avgDistance << "\n";
  std::cout << "Qoi liner error" << metrics.linearmaxError << std::endl;
   if (metrics.linearmaxError < tol * k_parameter)
    std::cout << "********** Successful **********\n";
  else
    std::cout << "********** Fail with error preservation **********\n";
  print_statistics(s, mode, shape, original_data, (float *)decompressed_data,
                      tol, config.normalize_coordinates);
  
  delete[](double *) original_data;
  delete[](unsigned char *) compressed_data;
  //free(compressed_data);
  //free(decompressed_data);
  delete[](T *) decompressed_data;
  return res;
}

bool run(int argc, char *argv[]) {

  int comm_size, rank;
  #if MYADIOS_USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #else
  rank = 0;
  comm_size = 1;
  #endif

  //std::cout << "rank " << rank << " of " << comm_size << "\n";

  std::string input_file = get_arg(argc, argv, "-i");
  std::string output_file = get_arg(argc, argv, "-c");

  if (rank == 0) std::cout << mgard_x::log::log_info << "original data: " << input_file
            << "\n";
  if (rank == 0) std::cout << mgard_x::log::log_info << "compressed data: " << output_file
            << "\n";

  enum mgard_x::data_type dtype;
  std::string dt = get_arg(argc, argv, "-t");
  if (dt.compare("s") == 0) {
    dtype = mgard_x::data_type::Float;
    if (rank == 0) std::cout << mgard_x::log::log_info << "data type: Single precision\n";
  } else if (dt.compare("d") == 0) {
    dtype = mgard_x::data_type::Double;
    if (rank == 0) std::cout << mgard_x::log::log_info << "data type: Double precision\n";
  } else
    if (rank == 0) print_usage_message("wrong data type.");

  mgard_x::DIM D = get_arg_int(argc, argv, "-n");
  std::vector<mgard_x::SIZE> shape = get_arg_dims(argc, argv, "-n");
  if (rank == 0) {
    std::string shape_string = "shape (";
    for (mgard_x::DIM d = 0; d < shape.size(); d++)
      shape_string = shape_string + std::to_string(shape[d]) + " ";
    shape_string = shape_string + ")";
  }
  
  
  bool non_uniform = false;
  std::string non_uniform_coords_file;
  if (has_arg(argc, argv, "-u")) {
    non_uniform = true;
    non_uniform_coords_file = get_arg(argc, argv, "-u");
    if (rank == 0) std::cout << mgard_x::log::log_info
              << "non-uniform coordinate file: " << non_uniform_coords_file
              << "\n";
  }

  enum mgard_x::error_bound_type mode; // REL or ABS
  std::string em = get_arg(argc, argv, "-m");
  if (em.compare("rel") == 0) {
    mode = mgard_x::error_bound_type::REL;
    if (rank == 0) std::cout << mgard_x::log::log_info << "error bound mode: Relative\n";
  } else if (em.compare("abs") == 0) {
    mode = mgard_x::error_bound_type::ABS;
    if (rank == 0) std::cout << mgard_x::log::log_info << "error bound mode: Absolute\n";
  } else
    if (rank == 0) print_usage_message("wrong error bound mode.");

  double tol = get_arg_double(argc, argv, "-e");
  double s = get_arg_double(argc, argv, "-s");

  if (rank == 0) {
    std::cout << std::scientific;
    std::cout << mgard_x::log::log_info << "error bound: " << tol << "\n";
    std::cout << std::defaultfloat;
    std::cout << mgard_x::log::log_info << "s: " << s << "\n";
  }

  int reorder = 0;
  if (has_arg(argc, argv, "-r")) {
    reorder = get_arg_int(argc, argv, "-r");
  }

  int lossless_level = get_arg_int(argc, argv, "-l");
  if (lossless_level == 0) {
    if (rank == 0) std::cout << mgard_x::log::log_info << "lossless: Huffman\n";
  } else if (lossless_level == 1) {
    if (rank == 0) std::cout << mgard_x::log::log_info << "lossless: Huffman + LZ4\n";
  } else if (lossless_level == 2) {
    if (rank == 0) std::cout << mgard_x::log::log_info << "lossless: Huffman + Zstd\n";
  }

  enum mgard_x::device_type dev_type;
  std::string dev = get_arg(argc, argv, "-d");
  if (dev.compare("auto") == 0) {
    dev_type = mgard_x::device_type::AUTO;
    if (rank == 0) std::cout << mgard_x::log::log_info << "device type: AUTO\n";
  } else if (dev.compare("serial") == 0) {
    dev_type = mgard_x::device_type::SERIAL;
    if (rank == 0) std::cout << mgard_x::log::log_info << "device type: SERIAL\n";
  } else if (dev.compare("openmp") == 0) {
    dev_type = mgard_x::device_type::OPENMP;
    if (rank == 0) std::cout << mgard_x::log::log_info << "device type: OPENMP\n";
  } else if (dev.compare("cuda") == 0) {
    dev_type = mgard_x::device_type::CUDA;
    if (rank == 0) std::cout << mgard_x::log::log_info << "device type: CUDA\n";
  } else if (dev.compare("hip") == 0) {
    dev_type = mgard_x::device_type::HIP;
    if (rank == 0) std::cout << mgard_x::log::log_info << "device type: HIP\n";
  } else if (dev.compare("sycl") == 0) {
    dev_type = mgard_x::device_type::SYCL;
    if (rank == 0) std::cout << mgard_x::log::log_info << "device type: SYCL\n";
  } else {
    if (rank == 0) print_usage_message("wrong device type.");
  }
  
  int k_parameter = 1;
  if (has_arg(argc, argv, "-q")) {
    k_parameter = get_arg_int(argc, argv, "-q");
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
    if (rank == 0) std::cout << mgard_x::log::log_info << "iterations: " << repeat << "\n";
  }

  int max_accumulate_data = 1;
  if (has_arg(argc, argv, "-a")) {
    max_accumulate_data = get_arg_int(argc, argv, "-a");
    if (rank == 0) std::cout << mgard_x::log::log_info << "max_accumulate data: " << max_accumulate_data << "\n";
  }

  int compute_delay = 0;
  if (has_arg(argc, argv, "-k")) {
    compute_delay = get_arg_int(argc, argv, "-k");
    if (rank == 0) std::cout << mgard_x::log::log_info << "compute delay: " << compute_delay << "\n";
  }


  mgard_x::SIZE max_memory_footprint = std::numeric_limits<mgard_x::SIZE>::max();
  if (has_arg(argc, argv, "-f")) {
    max_memory_footprint = (mgard_x::SIZE)get_arg_double(argc, argv, "-f");
  }

  
	std::vector<double> write, read, comp_no_pref, decomp_no_pref, comp_w_pref, decomp_w_pref, write_comp, read_comp;
  bool prefetch, use_compression;
  //std::vector<double> ebs = {16.73116, 1.673116, 0.1673116, 0.01673116}; 
  //std::cout<<std::endl;
  // for (double tol : ebs) {
  
    result res_comp_w_prefetch;
    use_compression = true; 
    prefetch = false;
    if (dtype == mgard_x::data_type::Double) {
      res_comp_w_prefetch = launch_compress<double>(
        D, dtype, input_file.c_str(), output_file.c_str(), shape, non_uniform,
        non_uniform_coords_file.c_str(), tol, s, mode, reorder,
        lossless_level, dev_type, num_dev, verbose, repeat, use_compression,
        max_accumulate_data, compute_delay, prefetch, max_memory_footprint, k_parameter);
    } else if (dtype == mgard_x::data_type::Float) {
      res_comp_w_prefetch = launch_compress<float>(
        D, dtype, input_file.c_str(), output_file.c_str(), shape, non_uniform,
        non_uniform_coords_file.c_str(), tol, s, mode, reorder,
        lossless_level, dev_type, num_dev, verbose, repeat, use_compression,
        max_accumulate_data, compute_delay, prefetch, max_memory_footprint, k_parameter);
    }
 
    comp_w_pref.push_back(res_comp_w_prefetch.comp);
    decomp_w_pref.push_back(res_comp_w_prefetch.decomp);
    write_comp.push_back(res_comp_w_prefetch.write);
    read_comp.push_back(res_comp_w_prefetch.read);
  // }

  if (rank == 0) {
  for (double time : write) std::cout << time << ", ";
  std::cout << "\n";
  for (double time : read) std::cout << time << ", ";
  std::cout << "\n";
  for (double time : comp_no_pref) std::cout << time << ", ";
  std::cout << "\n";
  for (double time : decomp_no_pref) std::cout << time << ", ";
  std::cout << "\n";
  for (double time : comp_w_pref) std::cout << time << ", ";
  std::cout << "\n";
  for (double time : decomp_w_pref) std::cout << time << ", ";
  std::cout << "\n";
  for (double time : write_comp) std::cout << time << ", ";
  std::cout << "\n";
  for (double time : read_comp) std::cout << time << ", ";
  std::cout << "\n";
	}
  return true;
}


int main(int argc, char *argv[]) {
  int comm_size;
  int rank;
#if MYADIOS_USE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
#else
 rank = 0;
 comm_size = 1;
#endif

run(argc, argv);
  #if MYADIOS_USE_MPI
  MPI_Finalize();
  #endif
  //MPI_Finalize();
  return 0;

}
