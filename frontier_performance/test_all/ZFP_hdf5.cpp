#include <ios>       //std::ios_base::failure
#include <iostream>  //std::cout
#include <numeric>   //std::iota
#include <stdexcept> //std::invalid_argument std::exception
#include <vector>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>

#include "zfp.h"

#define HDF5_USE_MPI 1
#include <hdf5.h>

#if HDF5_USE_MPI
#include <mpi.h>
#include "H5FDmpio.h" 
#endif

#if SIZE_MAX == UCHAR_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG
   //#error "what is happening here?"
#endif

using namespace std::chrono;
enum device_type { SERIAL, CUDA};; // 0 is serial and 1 is cuda
class Timer {
public:
  void start() { err = clock_gettime(CLOCK_REALTIME, &start_time); }
  void end() {
    err = clock_gettime(CLOCK_REALTIME, &end_time);
    total_time +=
        (double)(end_time.tv_sec - start_time.tv_sec) +
        (double)(end_time.tv_nsec - start_time.tv_nsec) / (double)1000000000;
  }
  double get() {
    double time =
        (double)(end_time.tv_sec - start_time.tv_sec) +
        (double)(end_time.tv_nsec - start_time.tv_nsec) / (double)1000000000;
    clear();
    return time;
  }
  void clear() { total_time = 0; }
  void print(std::string s) {
    std::cout << s << " time: " << total_time << "s" << std::endl;
    clear();
  }

private:
  int err = 0;
  double total_time = 0;
  struct timespec start_time, end_time;
};


void print_usage_message(std::string error) {
  if (error.compare("") != 0) {
    std::cout  << error << std::endl;
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
\t\t -m <abs|rel>: error bound mode (abs: abolute; rel: relative)\n\
\t\t -r <rate>: fixed rate (CUDA used)\n\
\t\t -e <accuracy>: accuracy (SERIAL used)\n\
\t\t -x <device>: SERIAL|CUDA \n\
\n\
\t -x: decompress data\n\
\t\t -c <path to compressed file>\n\
\t\t -o <path to decompressed file>\n");
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

std::vector<size_t> get_arg_dims(int argc, char *argv[], std::string option) {
  std::vector<size_t> shape;
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
      size_t d = std::stoi(arg);
      for (size_t i = 0; i < d; i++) {
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

template <typename T> void readfile(const char *input_file, T *in_buff, size_t read_size) {
  //std::cout << mgard_x::log::log_info << "Loading file: " << input_file << "\n";

  FILE *pFile;
  pFile = fopen(input_file, "rb");
  if (pFile == NULL) {
    std::cout << "file open error!\n";
    exit(1);
  }
  // fseek(pFile, 0, SEEK_END);
  // size_t lSize = ftell(pFile);
  // rewind(pFile);
  fread(in_buff, 1, read_size, pFile);
  fclose(pFile);
  // min_max(lSize/sizeof(T), in_buff);
}

// void readfile(const char *input_file, float *in_buff, size_t read_size) {
//   //std::cout <<  "Loading file: " << input_file << "\n";
//   FILE *pFile;
//   pFile = fopen(input_file, "rb");
//   if (pFile == NULL) {
//     std::cout << "file open error!\n";
//     exit(1);
//   }
//   // fseek(pFile, 0, SEEK_END);
//   // size_t lSize = ftell(pFile);
//   // rewind(pFile);
//   fread(in_buff, 1, read_size, pFile);
//   fclose(pFile);
//   // min_max(lSize/sizeof(T), in_buff);
// }

void writefile(const char *output_file, size_t num_bytes, float *out_buff) {
  FILE *file = fopen(output_file, "w");
  fwrite(out_buff, 1, num_bytes, file);
  fclose(file);
  int file_descriptor = fileno(file);
   if (file_descriptor >= 0) {
     if (fsync(file_descriptor) == -1) {
            perror("fsync error");
     }
   }
}



struct result {
	double write;
  double read;
  double comp;
  double decomp;
};


template <typename T>
result  launch_compress(int D, const char *input_file, const char *output_file, std::vector<size_t> shape, 
double accuracy,int num_sim_iterations,bool use_compression, int compute_delay, zfp_type type, int fixed_rate, enum device_type dev_type){
  int rank, size;
  #if HDF5_USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Info info = MPI_INFO_NULL;
  #else
  rank = 0;
  size = 1;
  #endif
  result res;

  size_t original_size = 1;
  for (int i = 0; i < D; i++) 
    original_size *= shape[i];
  T* original_data = new T[original_size];
 /**********************read information *************/
  readfile(input_file, original_data, original_size * sizeof(T));
  T *d_original_data;
  unsigned char *d_compressed_data;

  /*this will list the general zfp information*/
  int status = 0;    /* return value: 0 = success */
  //  type = zfp_type_float;
  //zfp_type type;     /* array scalar type */
  zfp_field* field;  /* array meta data */
  zfp_stream* zfp;   /* compressed stream */
  //void* buffer;      /* storage for compressed stream */
  size_t bufsize;    /* byte size of compressed buffer */
  bitstream* stream; /* bit stream to write to or read from */
  size_t compressed_size;    /* byte size of compressed stream */

  //void *compressed_data = (void *)new T[original_size];
  unsigned char * compressed_data = (unsigned char *) new T[original_size];
  compressed_size = original_size * sizeof(T);
  T * decompressed_data = (T *) new T[original_size];
 
  if (rank == 0) 
  {
    std::cout << "Data per rank: " << (T)original_size*sizeof(T)/1e9 << " GB\n";
  }
  /*************************************************************************/
  /**/
  // H5Pset_all_coll_metadata_ops(plistId, false);  // Disable all metadata operations
  // H5Pset_coll_metadata_write(plistId, false);  // Disable collective metadata writes
  // H5Pset_sieve_buf_size(plistId, 0);  // Set the sieve buffer size to 0 to disable data buffering
  // H5Pset_cache(plistId, 0, 0, 0, 0.0);
  // H5Pset_buffer(plistId, 1, NULL, NULL);
  // H5Pset_all_coll_metadata_ops(plistId, false);  // Disable all metadata operations
  // H5Pset_coll_metadata_write(plistId, false);  // Disable collective metadata writes
  // H5Pset_sieve_buf_size(plistId, 0);  // Set the sieve buffer size to 0 to disable data buffering
  // #ifdef HDF5_USE_MPI
    // hid_t plistId = H5Pcreate(H5P_FILE_ACCESS);
    // // H5Pset_all_coll_metadata_ops(plistId, false);  // Disable all metadata operations
    // // H5Pset_coll_metadata_write(plistId, false);  // Disable collective metadata writes
    // // H5Pset_sieve_buf_size(plistId, 0);  // Set the sieve buffer size to 0 to disable data buffering
    // // H5Pset_cache(plistId, 0, 0, 0, 0.0);
    // // H5Pset_buffer(plistId, 1, NULL, NULL);
    // H5Pset_fapl_mpio(plistId, MPI_COMM_WORLD, info);//"SDS.h5"
    // /*******/
    // h5_file = H5Fcreate(output_file, H5F_ACC_TRUNC, H5P_DEFAULT, plistId);
    // #else
    // hid_t plistId = H5Pcreate(H5P_FILE_ACCESS);
    // // H5Pset_all_coll_metadata_ops(plistId, false);  // Disable all metadata operations
    // // H5Pset_coll_metadata_write(plistId, false);  // Disable collective metadata writes
    // // H5Pset_sieve_buf_size(plistId, 0);  // Set the sieve buffer size to 0 to disable data buffering
    // h5_file = H5Fcreate(output_file, H5F_ACC_TRUNC, H5P_DEFAULT, plistId);
    // #endif
    // H5Pset_buffer(plistId, 1, NULL, NULL);
  /**/

  /*#ifdef HDF5_USE_MPI
  hid_t plistId = H5Pcreate(H5P_FILE_ACCESS);
  size_t meta_block_size = 1 ; // 512 KB
  size_t sieve_buffer_size = 1 ; // 256 KB
  H5Pset_cache(plistId, 0, meta_block_size, sieve_buffer_size, 0);
  H5Pset_fapl_mpio(plistId, MPI_COMM_WORLD, info);//"SDS.h5"*/
  /*******/
  /*hid_t h5_file = H5Fcreate(output_file, H5F_ACC_TRUNC, H5P_DEFAULT, plistId);
  #else
  hid_t plistId = H5Pcreate(H5P_FILE_ACCESS);
  hid_t h5_file = H5Fcreate(output_file, H5F_ACC_TRUNC, H5P_DEFAULT, plistId);
  #endif*/

  /*#ifdef HDF5_USE_MPI
  H5Pclose(plistId);
  #endif*/
  size_t var_size;
  unsigned char *var_data;
  unsigned char *var_data1 = new unsigned char[sizeof(T) * original_size];
  
  hid_t h5_file;
  std::string out_fin;
  std::string prefix = "/lustre/orion/proj-shared/csc303/zq53/mgard_new/adios_more/hdf5_all/build/"; //this is used for the write file pat
  

// Get the current time before the operation
  Timer timer, timer_total;
  #ifdef HDF5_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif
  timer_total.start();
  for (int sim_iter = 0; sim_iter < num_sim_iterations; sim_iter++) 
  {
    sleep(compute_delay);
    if (use_compression) /***use compression or not ***********/
    {
      // ********* Compression ********** //
      #ifdef HDF5_USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
      #endif
      timer.start();
      if (dev_type == 0){
        if (D == 1) {
        field = zfp_field_1d(original_data, type, shape[0]);
      } else if (D == 2) {
        field = zfp_field_2d(original_data, type, shape[1], shape[0]);
      } else if (D == 3) {
        field = zfp_field_3d(original_data, type, shape[2], shape[1], shape[0]);
      } else if (D == 4) {
        field = zfp_field_4d(original_data, type, shape[3], shape[2], shape[1], shape[0]);
      } else {
        std::cout << "wrong D\n";
        exit(-1);
      }
      zfp = zfp_stream_open(NULL);
      zfp_stream_set_accuracy(zfp, accuracy);
      bufsize = zfp_stream_maximum_size(zfp, field);
      stream = stream_open(compressed_data, sizeof(T) * original_size);
      zfp_stream_set_bit_stream(zfp, stream);
      zfp_stream_rewind(zfp);
      compressed_size = zfp_compress(zfp, field);
      }
    #ifdef HDF5_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.end();
    double compression_time = timer.get();
    #ifdef HDF5_USE_MPI
    double max_compression_time;
    MPI_Reduce(&compression_time, &max_compression_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
      res.comp = max_compression_time;
      std::cout  << "Compression time: " << max_compression_time << "\n";
      std::cout  << "Compression throughput: " << (double)original_size * sizeof(T) * size/max_compression_time/1e9 << " GB/s.\n";
      std::cout  << "Compression ratio: "<< (double)original_size * sizeof(T) / (compressed_size) << "\n";
    }
    #else
      if (rank == 0)
      {
        res.comp = compression_time;
        std::cout  << "Compression time: " << compression_time << "\n";
        std::cout  << "Compression throughput: " << (double)original_size * sizeof(T) * size/compression_time/1e9 << " GB/s.\n";
        std::cout  << "Compression ratio: "<< (double)original_size * sizeof(T) / (compressed_size) << "\n";
      }
    #endif
    timer.clear();
    }
    // ********* Write********** //
    var_size = use_compression ? compressed_size : original_size * sizeof(T);
    var_data = use_compression ? (unsigned char *)compressed_data : (unsigned char *) original_data;
    //var_data1 = (unsigned char *) original_data;
    std::string var_name = "var_step" + std::to_string(sim_iter);
    /*******need? or not */
    size_t offset;
    size_t gloabl_size;
    #ifdef HDF5_USE_MPI
    MPI_Scan(&var_size, &offset, 1, my_MPI_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&var_size, &gloabl_size, 1, my_MPI_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
    offset -= compressed_size;
    hsize_t dimsf[1] = {gloabl_size};
    //hsize_t dimsf[1] = {1};
    hsize_t count[1] = {var_size};
    hsize_t start[1] = {offset};
    #else
     gloabl_size = var_size;
     hsize_t dimsf[1] = {gloabl_size};
    //hsize_t dimsf[1] = {1};
    hsize_t count[1] = {var_size};
    hsize_t start[1] = {0};
    #endif
    
    /***********/
    //timer.start();
    #ifdef HDF5_USE_MPI
    hid_t plistId = H5Pcreate(H5P_FILE_ACCESS);
    //H5Pset_fapl_mpio(plistId, MPI_COMM_WORLD, info);//"SDS.h5"
    /*size_t meta_block_size = 0 ; // 
    size_t sieve_buffer_size = 0 ; // 
    H5Pset_cache(plistId, 0, meta_block_size, sieve_buffer_size, 0);
    H5Pset_page_buffer_size(plistId, 0, 0, 0);
    H5Pset_fapl_core(plistId, 0, true);
    H5Pset_small_data_block_size(plistId, 0);
    H5Pset_meta_block_size(plistId, 0);*/
    H5Pset_fapl_mpio(plistId, MPI_COMM_WORLD, info);//"SDS.h5"
    /*******/
    std::string outputPath(output_file);
    std::string reducedPath = outputPath.substr(outputPath.rfind("/") + 1);
    if (dev_type == 0)
      out_fin = prefix + std::to_string(accuracy) + "_" + reducedPath;
    else if (dev_type == 1)
      out_fin = prefix + std::to_string(fixed_rate) + "_" + reducedPath; 
    hid_t h5_file = H5Fcreate(out_fin.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plistId);
    #else
    std::string outputPath(output_file);
    std::string reducedPath = outputPath.substr(outputPath.rfind("/") + 1);
    if (dev_type == 0)
      out_fin = prefix + std::to_string(accuracy) + "_" + reducedPath;
    else if (dev_type == 1)
      out_fin = prefix + std::to_string(fixed_rate) + "_" + reducedPath;
    hid_t plistId = H5Pcreate(H5P_FILE_ACCESS);
    /*size_t meta_block_size = 0 ; // 
    size_t sieve_buffer_size = 0 ; // 
    H5Pset_cache(plistId, 0, meta_block_size, sieve_buffer_size, 0);
    H5Pset_page_buffer_size(plistId, 0, 0, 0);
    H5Pset_fapl_core(plistId, 0, true);
    H5Pset_small_data_block_size(plistId, 0);
    H5Pset_meta_block_size(plistId, 0);*/
    hid_t h5_file = H5Fcreate(out_fin.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plistId);
    #endif
    H5Pclose(plistId); 
    hid_t gspace = H5Screate_simple(1, dimsf, NULL);  //(H5S_SCALAR)
    hid_t data_plist = H5Pcreate(H5P_DATASET_ACCESS);
    //H5Pset_chunk_cache(data_plist, 0, 0, 0); //H5Pset_chunk_cache(dapl_id, 12421, 16*1024*1024, H5D_CHUNK_CACHE_W0_DEFAULT);
    //H5Pset_buffer(data_plist, 100, NULL, NULL);
    hid_t create_plist = H5Pcreate(H5P_DATASET_CREATE);
   // H5Pset_alloc_time(create_plist, H5D_ALLOC_TIME_LATE);
   // H5Pset_fill_time(create_plist, H5D_FILL_TIME_ALLOC);
    //H5Pset_fill_time(create_plist, H5D_FILL_TIME_ALLOC);//H5D_FILL_TIME_NEVER
    #ifdef HDF5_USE_MPI
    hid_t dataset = H5Dcreate(h5_file, var_name.c_str(), H5T_NATIVE_UCHAR, gspace, H5P_DEFAULT, create_plist, data_plist);
    #else
    hid_t dataset = H5Dcreate(h5_file, var_name.c_str(), H5T_NATIVE_UCHAR, gspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    #endif
    std::cout << "ppppath " << out_fin.c_str() << " " << h5_file <<std::endl;
    /*****notatoion temp*/
    hid_t fspace = H5Dget_space(dataset);
    H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, count, NULL);
    hid_t mspace = H5Screate_simple(1, count, NULL);
    #ifdef HDF5_USE_MPI
    hid_t plistId1 = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plistId1, H5FD_MPIO_COLLECTIVE); // H5FD_MPIO_INDEPENDENT H5FD_MPIO_COLLECTIVE
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    //H5Fflush(h5_file, H5F_SCOPE_GLOBAL);
    timer.start();
    #ifdef HDF5_USE_MPI
    //H5Dwrite(dataset, H5T_NATIVE_UCHAR, H5S_ALL, fspace, plistId, var_data);

    H5Dwrite(dataset, H5T_NATIVE_UCHAR, mspace, fspace, plistId1, var_data);
    H5Fflush(h5_file, H5F_SCOPE_LOCAL); //H5F_SCOPE_LOCAL H5F_SCOPE_GLOBAL
    H5Pclose(plistId1);
    #else
    H5Dwrite(dataset, H5T_NATIVE_UCHAR, mspace, fspace, H5P_DEFAULT, var_data);
    H5Fflush(h5_file, H5F_SCOPE_LOCAL); //H5F_SCOPE_LOCAL H5F_SCOPE_GLOBAL
    #endif
    //std::cout << "- First buffer written" << std::endl;
    H5Dclose(dataset);//H5Sclose(plistId); 
    H5Sclose(gspace); H5Sclose(fspace); H5Sclose(mspace); H5Fclose(h5_file);
    //delete[](T *) original_data;
    #ifdef HDF5_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.end();
    double write_time = timer.get();
     #ifdef HDF5_USE_MPI
    double max_write_time;
    MPI_Reduce(&write_time, &max_write_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) 
    {
      res.write = max_write_time;
      std::cout << "Write time: "<< max_write_time << "\n";
      std::cout << "Write throughput: "<< (double)var_size* size/max_write_time/1e9 << " GB/s.\n";
    }
    #else
    if (rank == 0) 
    {
      res.write = write_time;
      std::cout << "Write time: "<< write_time << "\n";
      std::cout << "Write throughput: "<< (double)var_size* size/write_time/1e9 << " GB/s.\n";
    }
    #endif
    timer.clear();
  }
  //[](T *) original_data;
  #ifdef HDF5_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif
  for (int sim_iter = 0; sim_iter < num_sim_iterations; sim_iter++) 
  {
    sleep(compute_delay);
    //std::string var_name = "var_step" + std::to_string(0);
    // ********* Read ********** //
    std::string var_name = "var_step" + std::to_string(sim_iter);
    size_t offset;
    size_t gloabl_size;
    #ifdef HDF5_USE_MPI
    MPI_Scan(&var_size, &offset, 1, my_MPI_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&var_size, &gloabl_size, 1, my_MPI_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
    offset -= compressed_size;
    hsize_t dimsf[1] = {gloabl_size};
    hsize_t count[1] = {var_size};
    hsize_t start[1] = {offset};
    #else
    gloabl_size = var_size;
    hsize_t dimsf[1] = {gloabl_size};
    //hsize_t dimsf[1] = {1};
    hsize_t count[1] = {var_size};
    hsize_t start[1] = {0};
    #endif
    std::cout << "var_size "<<gloabl_size << "dimsf " << dimsf[0] << "count "<< count[0] << "start " << start[0]<< std::endl;
    timer.start();
    //hid_t plistId = H5Pcreate(H5P_FILE_ACCESS);
    //H5Pset_fapl_mpio(plistId, MPI_COMM_WORLD, info);
    hid_t file_id_r = H5Fopen(out_fin.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    //hid_t file_id_r = H5Fopen("16.731160_test.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    std::cout << "ppppath " << out_fin.c_str() << " " << file_id_r<<std::endl;
    hid_t dataset_r = H5Dopen(file_id_r, (var_name).c_str(), H5P_DEFAULT);
    std::cout << (var_name).c_str() << std::endl;
    // H5Dset_extent(dataset, dimsf);
    hid_t fspace_r = H5Dget_space(dataset_r);
    std::cout << "- Dataset get" << std::endl;
    H5Sselect_hyperslab(fspace_r, H5S_SELECT_SET, start, NULL, count, NULL);
    hid_t mspace_r = H5Screate_simple(1, count, NULL);
    #ifdef HDF5_USE_MPI
    hid_t plistId_tran = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plistId_tran, H5FD_MPIO_COLLECTIVE);

    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    // timer.start();
    #ifdef HDF5_USE_MPI
    // H5Dread(dataset, H5T_NATIVE_UCHAR, mspace, fspace, plistId, var_data);
    H5Drefresh(dataset_r);
    H5Dread(dataset_r, H5T_NATIVE_UCHAR, mspace_r, fspace_r, H5P_DEFAULT, var_data1);
    H5Pclose(plistId_tran);
    #else
    //H5Dread(dataset, H5T_NATIVE_UCHAR, H5S_ALL, fspace, H5P_DEFAULT, var_data);
    H5Drefresh(dataset_r);  
    std::cout << "strat success " <<std::endl;
    H5Dread(dataset_r, H5T_NATIVE_UCHAR, mspace_r, fspace_r, H5P_DEFAULT, var_data1);
    #endif
    H5Dclose(dataset_r);
    H5Sclose(fspace_r);
    H5Sclose(mspace_r);
    H5Fclose(file_id_r);
    //H5Fclose(h5_file);
    //delete[](T *) original_data;
    #ifdef HDF5_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    // H5Fclose(h5_file);
    timer.end();
    double read_time = timer.get();
    #ifdef HDF5_USE_MPI
    double max_read_time;
    MPI_Reduce(&read_time, &max_read_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
    res.read = max_read_time;
    std::cout <<  "Read time: " << max_read_time << "\n";
    std::cout <<  "Read throughput: "<< (double)var_size* size/max_read_time/1e9 << " GB/s.\n";
    }
    #else
    if (rank == 0) {
    res.read = read_time;
    std::cout <<  "Read time: " << read_time << "\n";
    std::cout <<  "Read throughput: "<< (double)var_size* size/read_time/1e9 << " GB/s.\n";
    }
    #endif
    timer.clear();
    if (use_compression) {
      // ********* Decompression ********** //
      #if HDF5_USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
      #endif
      timer.start();
      if (dev_type ==0){
        if (D == 1) {
        field = zfp_field_1d(decompressed_data, type, shape[0]);
        } else if (D == 2) {
        field = zfp_field_2d(decompressed_data, type, shape[1], shape[0]);
        } else if (D == 3) {
          field = zfp_field_3d(decompressed_data, type, shape[2], shape[1], shape[0]);
        } else if (D == 4) {
          field = zfp_field_4d(decompressed_data, type, shape[3], shape[2], shape[1], shape[0]);
        } else {
          std::cout << "wrong D\n";
          exit(-1);
        }
        // zfp = zfp_stream_open(NULL);
        zfp_stream_set_accuracy(zfp, accuracy);
        bufsize = zfp_stream_maximum_size(zfp, field);
        stream = stream_open(var_data, sizeof(T) * original_size);
        zfp_stream_set_bit_stream(zfp, stream);
        zfp_stream_rewind(zfp);
        size_t result = zfp_decompress(zfp, field);
      }

     //}
      #if HDF5_USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
      #endif
      timer.end();
      double decompress_time = timer.get();
      #ifdef HDF5_USE_MPI
      double max_decompress_time;
      MPI_Reduce(&decompress_time, &max_decompress_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
     
      if (rank == 0) 
      {
        res.decomp = max_decompress_time;
        std::cout <<  "Decompression time: "<< max_decompress_time << "\n";
        std::cout <<  "Decompression throughput: "<< (double)original_size * sizeof(T) * size/max_decompress_time/1e9 << " GB/s.\n";
      }
      #else
      if (rank == 0) 
      {
        res.decomp = decompress_time;
        std::cout <<  "Decompression time: "<< decompress_time << "\n";
        //std::cout <<  "Decompression1 time: "<< dcompression_time_1 / 1000 << "\n";
        std::cout <<  "Decompression throughput: "<< (double)original_size * sizeof(T) * size/decompress_time/1e9 << " GB/s.\n";
      }
      #endif
        
      timer.clear();
      zfp_field_free(field);
      zfp_stream_close(zfp);
      stream_close(stream);
      }
  }
  
  delete[](T *) original_data;
  delete[](unsigned char *) compressed_data;
  delete[] (T *)decompressed_data;
  //int c;
  //while ((c = getchar()) != '\n' && c != EOF) { }
  //return 0;
  return res;
}

bool run(int argc, char *argv[]) {
  int rank, size;
  int num_sim_iterations = 1;

  #if HDF5_USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  #else
  rank = 0;
  size = 1;
  #endif
  
  std::string input_file = get_arg(argc, argv, "-i");
  std::string output_file = get_arg(argc, argv, "-c");
  if (rank == 0) 
    std::cout <<  "original data: " << input_file << "\n";
  if (rank == 0) 
    std::cout <<  "compressed data: " << output_file << "\n";
  
  int DIM = get_arg_int(argc, argv, "-n");
  std::vector<size_t> shape = get_arg_dims(argc, argv, "-n");

  double accuracy = get_arg_double(argc, argv, "-e");
  int fixed_rate = get_arg_int(argc, argv, "-r");
  zfp_type type; 
  std::string dt = get_arg(argc, argv, "-t");
  if (dt.compare("s") == 0) {
      type = zfp_type_float;
      if (rank == 0) std::cout <<  "data type: Single precision\n";
    } else if (dt.compare("d") == 0) {
      type = zfp_type_double;
      if (rank == 0) std::cout << "data type: Double precision\n";
    } else
      if (rank == 0) print_usage_message("wrong data type.");

    int compute_delay = 0;
    if (has_arg(argc, argv, "-k")) {
      compute_delay = get_arg_int(argc, argv, "-k");
      if (rank == 0) std::cout << "compute delay: " << compute_delay << "\n";
    }

  
  device_type dev_type;
  std::string dev = get_arg(argc, argv, "-x");
  if (dev.compare("serial") == 0) {
    dev_type = SERIAL; // 
    if (rank == 0) std::cout << "device type: SERIAL\n";
  } else if (dev.compare("cuda") == 0) {
    dev_type = CUDA;
    if (rank == 0) std::cout <<  "device type: CUDA\n";
  } else {
    if (rank == 0) print_usage_message("wrong device type.");
  }

    
    std::vector<double> write, read, comp, decomp, write_comp, read_comp;
    bool use_compression;
    //std::vector<double> ebs = {16.73116, 1.673116, 0.1673116, 0.01673116};
    std::vector<double> ebs = {16.73116};
    std::vector<int> rate = {4, 7, 12, 15};
    result res_no_comp;
    for (int i=0; i<1;i++){
    use_compression = false;
    if (type == zfp_type_float)
    {
      res_no_comp = launch_compress<float>(
        DIM, input_file.c_str(), output_file.c_str(), shape, accuracy, 1, use_compression, compute_delay, type, fixed_rate, dev_type);}
    else if (type == zfp_type_double) 
    {
      res_no_comp = launch_compress<double>(
        DIM, input_file.c_str(), output_file.c_str(), shape, accuracy, 1, use_compression, compute_delay, type, fixed_rate, dev_type);}
    }
      if (dev_type == SERIAL){
      for (double accuracy : ebs) {
      if(rank == 0)
      std::cout << "this is the error "<<accuracy<<std::endl;
      result res_comp;
      use_compression = true; 
      if (type == zfp_type_float)
      {
      res_comp = launch_compress<float>(
          DIM, input_file.c_str(), output_file.c_str(), shape, accuracy, 1, use_compression, compute_delay, type, fixed_rate, dev_type); }
      else if (type == zfp_type_double)
      {
      res_comp = launch_compress<double>(
          DIM, input_file.c_str(), output_file.c_str(), shape, accuracy, 1, use_compression, compute_delay, type, fixed_rate, dev_type); }
      comp.push_back(res_comp.comp);
      decomp.push_back(res_comp.decomp);
      write_comp.push_back(res_comp.write);
      read_comp.push_back(res_comp.read);
    }} else if (dev_type == CUDA) 
     { for (int fixed_rate : rate) {
      if(rank == 0)
      std::cout << "this is the rate "<<fixed_rate<<std::endl;
      result res_comp;
      use_compression = true;
      if (type == zfp_type_float)
      {
      res_comp = launch_compress<float>(
          DIM, input_file.c_str(), output_file.c_str(), shape, accuracy, 1, use_compression, compute_delay, type, fixed_rate, dev_type); }
      else if (type == zfp_type_double)
      {
      res_comp = launch_compress<double>(
          DIM, input_file.c_str(), output_file.c_str(), shape, accuracy, 1, use_compression, compute_delay, type, fixed_rate, dev_type); }
      comp.push_back(res_comp.comp);
      decomp.push_back(res_comp.decomp);
      write_comp.push_back(res_comp.write);
      read_comp.push_back(res_comp.read);
    }
}

    return true;
}


int main(int argc, char *argv[]) {
  int size;
  int rank;

  #if HDF5_USE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  #else
  rank = 0;
  size = 1;
  #endif

  run(argc, argv);
  #if HDF5_USE_MPI
  MPI_Finalize();
  #endif
  return 0;
}
