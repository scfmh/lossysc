#include <ios>       //std::ios_base::failure
#include <iostream>  //std::cout
#include <numeric>   //std::iota
#include <stdexcept> //std::invalid_argument std::exception
#include <vector>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unistd.h>

#include "cusz.h"
#include "cuszapi.hh"
#include "utils/io.hh"
#include "cli/timerecord_viewer.hh"

//#define HDF5_USE_MPI 1
#include <hdf5.h>
#ifdef HDF5_USE_MPI
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

class Timer {
public:
  void start() { err = clock_gettime(CLOCK_REALTIME, &start_time); }
  void end() {
    err = clock_gettime(CLOCK_REALTIME, &end_time);
    total_time +=
        (double)(end_time.tv_sec - start_time.tv_sec) +
        (double)(end_time.tv_nsec - start_time.tv_nsec) / (double)1000000000;
        //std::cout << " time second " << start_time.tv_sec << " " << end_time.tv_sec <<" " << end_time.tv_sec - start_time.tv_sec << std::endl;
        //std::cout << " nano second " << end_time.tv_nsec << " "  << start_time.tv_nsec << " " << end_time.tv_nsec - start_time.tv_nsec << std::endl;
  }
  double get() {
    double time =
        (double)(end_time.tv_sec - start_time.tv_sec) +
        (double)( double(end_time.tv_nsec) - double(start_time.tv_nsec)) / (double)1000000000;
    //std::cout << " time second " << start_time.tv_sec << " " << end_time.tv_sec <<" " << end_time.tv_sec - start_time.tv_sec << std::endl;
       //std::cout << " nano second " << double (end_time.tv_nsec) << " "  <<double (start_time.tv_nsec) << " " << end_time.tv_nsec - start_time.tv_nsec << std::endl;
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
\t\t -e <error>: error tolerance\n\
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

void readfile(const char *input_file, float *in_buff, size_t read_size) {
  std::cout <<  "Loading file: " << input_file << "\n";
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

void writefile(const char *output_file, size_t num_bytes, float *out_buff) {
  FILE *file = fopen(output_file, "w");
  fwrite(out_buff, 1, num_bytes, file);
  fclose(file);
}



struct result {
	double write;
  double read;
  double comp;
  double decomp;
};


template <typename T>
result launch_compress(int DIM, const char *input_file, const char *output_file, std::vector<size_t> shape, 
 double tol,int num_sim_iterations,bool use_compression, int compute_delay){
 int rank, size;

 #ifdef HDF5_USE_MPI
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 MPI_Comm_size(MPI_COMM_WORLD, &size);
 MPI_Info info = MPI_INFO_NULL;
#else
 rank = 0;
 size = 1;
#endif

result res;

 auto original_size = 1; 
 for (int i = 0; i < DIM; i++)
   original_size *= shape[i];


 /**********************read information *************/
  T *d_original_data;
  unsigned char *d_compressed_data;
  if (rank == 0) 
 {
   std::cout << "Data per rank: " << (double)original_size*sizeof(float)/1e9 << " GB\n";
 }

/********************************************************************/
    cusz_header header;
    T * original_data = (T *) new T[original_size];
    io::read_binary_to_array(input_file, original_data, original_size);
    unsigned char * compressed_data = (unsigned char *) new T[original_size];
    size_t compressed_size = original_size * sizeof(T);
  
    size_t uncompressed_memlen = original_size * 1.03;
    size_t decompressed_memlen = uncompressed_memlen;

    /*cudaStream_t stream;
    cudaStreamCreate(&stream);

    cusz_framework* framework = new cusz_custom_framework{
        .pipeline     = Auto,
        .predictor    = cusz_custom_predictor{.type = LorenzoI},
        .quantization = cusz_custom_quantization{.radius = 512},
        .codec        = cusz_custom_codec{.type = Huffman}};

    cusz_compressor* comp       = cusz_create(framework, FP32);
    cusz_config*     config     = new cusz_config{.eb = tol, .mode = Abs};
    cusz_len uncomp_len;
      if (DIM == 1) {
      uncomp_len = cusz_len{shape[0], 1, 1, 1, 1.03};
    } else if (DIM == 2) {
      uncomp_len = cusz_len{shape[1], shape[0], 1, 1, 1.03};
    } else if (DIM == 3) {
      uncomp_len = cusz_len{shape[2], shape[1], shape[0], 1, 1.03};
    } else if (DIM == 4) {
      uncomp_len = cusz_len{shape[3], shape[2], shape[1], shape[0], 1.03};
    } else {
      std::cout << "wron demisions\n";
      exit(-1);
    }

    cusz_len   decomp_len = uncomp_len;*/

    cusz_len uncomp_len;
      if (DIM == 1) {
        uncomp_len = cusz_len{shape[0], 1, 1, 1, 1.03};
      } else if (DIM == 2) {
        uncomp_len = cusz_len{shape[1], shape[0], 1, 1, 1.03};
      } else if (DIM == 3) {
        uncomp_len = cusz_len{shape[2], shape[1], shape[0], 1, 1.03};
      } else if (DIM == 4) {
        uncomp_len = cusz_len{shape[3], shape[2], shape[1], shape[0], 1.03};
      } else {
        std::cout << "wron demisions\n";
        exit(-1);
      }

    cusz_len   decomp_len = uncomp_len;
    cusz::TimeRecord compress_timerecord;
    cusz::TimeRecord decompress_timerecord;

   

  /*************************************************************************/
    /*#ifdef HDF5_USE_MPI
    hid_t plistId = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plistId, MPI_COMM_WORLD, info);//"SDS.h5"
    hid_t h5_file = H5Fcreate(output_file, H5F_ACC_TRUNC, H5P_DEFAULT, plistId);
    #else
    hid_t h5_file = H5Fcreate(output_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    #endif
  
    #ifdef HDF5_USE_MPI
    H5Pclose(plistId);
    #endif*/

    size_t var_size;
    unsigned char *var_data;
    hid_t h5_file;
    std::string out_fin;
    std::string prefix = "/gpfs/alpine/proj-shared/csc143/zq53/mgard_new/adios_more/hdf5_all/build/"; //this is used for the write file path 
    Timer timer, timer_total;
    cudaEvent_t timer1, end;
    cudaEventCreate(&timer1);
    cudaEventCreate(&end);
    #ifdef HDF5_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    //timer_total.start();
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
      cudaEventRecord(timer1, 0);
      cudaStream_t stream;
      cudaStreamCreate(&stream);

      cusz_framework* framework = new cusz_custom_framework{
        .pipeline     = Auto,
        .predictor    = cusz_custom_predictor{.type = LorenzoI},
        .quantization = cusz_custom_quantization{.radius = 512},
        .codec        = cusz_custom_codec{.type = Huffman}};

      cusz_compressor* comp       = cusz_create(framework, FP32);
      cusz_config*     config     = new cusz_config{.eb = tol, .mode = Abs};

      cudaMalloc(&d_original_data, sizeof(T) * original_size);
      cudaMemcpy(d_original_data, original_data, sizeof(T) * original_size, cudaMemcpyHostToDevice);
      cusz_compress(comp, config,  d_original_data, uncomp_len, &d_compressed_data, &compressed_size, &header,(void*)&compress_timerecord, stream);
      //cusz_compress(comp, config, d_uncompressed, uncomp_len, &exposed_compressed, &compressed_len, &header,(void*)&compress_timerecord, stream);
      cudaStreamSynchronize(stream);
      cudaMemcpy(compressed_data, d_compressed_data, compressed_size, cudaMemcpyDeviceToHost);
      cudaFree(d_original_data);
      cudaFree(d_compressed_data);
      cusz_release(comp);
      cudaStreamDestroy(stream);
      cudaDeviceSynchronize();
      #ifdef HDF5_USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
      #endif
      timer.end();
      cudaEventRecord(end, 0);
      cudaEventSynchronize(end);
      double compression_time = timer.get();
      float compression_time_1;
      cudaEventElapsedTime(&compression_time_1, timer1, end);
      if (rank==0)
      {
      //std::cout << "rank " << rank << " compression time: " << compression_time << "\n";
      cusz::TimeRecordViewer::view_compression(&compress_timerecord, original_size * sizeof(T), compressed_size);
      }
      #ifdef HDF5_USE_MPI
      float max_compression_time_1;
      MPI_Reduce(&compression_time_1, &max_compression_time_1, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
      double max_compression_time;
      MPI_Reduce(&compression_time, &max_compression_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if (rank == 0)
      {
        res.comp = max_compression_time;
        std::cout  << "Compression time: " << max_compression_time << "\n";
        std::cout  << "Compression1 time: " << max_compression_time_1 / 1000<< "\n";
        std::cout  << "Compression throughput: " << (double)original_size * sizeof(T) * size/max_compression_time/1e9 << " GB/s.\n";
        std::cout  << "Compression ratio: "<< (double)original_size * sizeof(T) / (compressed_size) << "\n";
      }
      #else
        if (rank == 0)
        {
          res.comp = compression_time;
          std::cout  << "Compression time: " << compression_time << "\n";
          std::cout  << "Compression1 time: " << compression_time_1 / 1000<< "\n";
          std::cout  << "Compression throughput: " << (double)original_size * sizeof(T) * size/compression_time/1e9 << " GB/s.\n";
          std::cout  << "Compression ratio: "<< (double)original_size * sizeof(T) / (compressed_size) << "\n";
        }
      #endif
        timer.clear();
        cudaEventDestroy(timer1);
        cudaEventDestroy(end);
      } 


    // ********* Write********** //
    var_size = use_compression ? compressed_size : original_size * sizeof(T);
    var_data = use_compression ? (unsigned char *)compressed_data : (unsigned char *) original_data;
    std::string var_name = "var_step" + std::to_string(sim_iter);
    size_t offset;
    size_t gloabl_size;
    MPI_Scan(&var_size, &offset, 1, my_MPI_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&var_size, &gloabl_size, 1, my_MPI_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
    offset -= compressed_size;
    hsize_t dimsf[1] = {gloabl_size};
    hsize_t count[1] = {var_size};
    hsize_t start[1] = {offset};
    /***********/
    timer.start();
    #ifdef HDF5_USE_MPI
    hid_t plistId = H5Pcreate(H5P_FILE_ACCESS);
    /*size_t meta_block_size = 0 ; // 
    size_t sieve_buffer_size = 0 ; // 
    H5Pset_cache(plistId, 0, meta_block_size, sieve_buffer_size, 0);
    H5Pset_page_buffer_size(plistId, 0, 0, 0);
    H5Pset_fapl_core(plistId, 0, true);
    H5Pset_small_data_block_size(plistId, 0);
    H5Pset_meta_block_size(plistId, 0);*/
    H5Pset_fapl_mpio(plistId, MPI_COMM_WORLD, info);//"SDS.h5"
    std::string outputPath(output_file); 
    std::string reducedPath = outputPath.substr(outputPath.rfind("/") + 1);
    //out_fin =  prefix + std::to_string(tol) + "_" + output_file;
    out_fin =  prefix + std::to_string(tol) + "_" + reducedPath;
    std::cout << out_fin << std::endl;
    h5_file = H5Fcreate(out_fin.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plistId);
    #else
    hid_t plistId = H5Pcreate(H5P_FILE_ACCESS);
    std::string outputPath(output_file); 
    std::string reducedPath = outputPath.substr(outputPath.rfind("/") + 1);
    //out_fin =  prefix + std::to_string(tol) + "_" + output_file;
    out_fin =  prefix + std::to_string(tol) + "_" + reducedPath;
    std::cout << out_fin << std::endl;
    h5_file = H5Fcreate(out_fin.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plistId);
    #endif
 
    #ifdef HDF5_USE_MPI
    H5Pclose(plistId);
    #endif

    
    hid_t gspace = H5Screate_simple(1, dimsf, NULL);
    hid_t data_plist = H5Pcreate(H5P_DATASET_ACCESS);
    hid_t create_plist = H5Pcreate(H5P_DATASET_CREATE);
    // if not single, notation the next two
    //H5Pset_alloc_time(create_plist, H5D_ALLOC_TIME_LATE);
    //H5Pset_fill_time(create_plist, H5D_FILL_TIME_ALLOC);//H5D_FILL_TIME_NEVER
    hid_t dataset = H5Dcreate(h5_file, var_name.c_str(), H5T_NATIVE_UCHAR, gspace, H5P_DEFAULT, create_plist, data_plist);
    //hid_t dataset = H5Dcreate(h5_file, var_name.c_str(), H5T_NATIVE_UCHAR, gspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    /*****notatoion temp*/
    hid_t fspace = H5Dget_space(dataset);
    H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, count, NULL);
    //std::cout << "- First hyperslab selected" << std::endl;
    hid_t mspace = H5Screate_simple(1, count, NULL);
    #ifdef HDF5_USE_MPI
    hid_t plistId1 = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plistId1, H5FD_MPIO_COLLECTIVE);

    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    // timer.start();
    #ifdef HDF5_USE_MPI
    H5Dwrite(dataset, H5T_NATIVE_UCHAR, mspace, fspace, plistId1, var_data);
    H5Fflush(h5_file, H5F_SCOPE_LOCAL);
    H5Pclose(plistId1);
    #else
    H5Dwrite(dataset, H5T_NATIVE_UCHAR, mspace, fspace, H5P_DEFAULT, var_data);
    #endif
    //std::cout << "- First buffer written" << std::endl;
    H5Dclose(dataset); 
    H5Sclose(gspace); H5Sclose(fspace); H5Sclose(mspace); H5Fclose(h5_file);
    #ifdef HDF5_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.end();

    double write_time = timer.get();
    #ifdef HDF5_USE_MPI
    double max_write_time;
    MPI_Reduce(&write_time, &max_write_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
     res.write = max_write_time;
      std::cout << "Write time: "<< max_write_time << "\n";
      std::cout << "Write throughput: "<< (double)var_size* size/max_write_time/1e9 << " GB/s.\n";
    }
    #else
    if (rank == 0) {
     res.write = write_time;
      std::cout << "Write time: "<< write_time << "\n";
      std::cout << "Write throughput: "<< (double)var_size* size/write_time/1e9 << " GB/s.\n";
    }
    #endif
    timer.clear(); 

    }
    #ifdef HDF5_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    for (int sim_iter = 0; sim_iter < num_sim_iterations; sim_iter++) 
    {
      sleep(compute_delay);

      std::string var_name = "var_step" + std::to_string(sim_iter);
    size_t offset;
    size_t gloabl_size;
    MPI_Scan(&var_size, &offset, 1, my_MPI_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&var_size, &gloabl_size, 1, my_MPI_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
    offset -= compressed_size;
    hsize_t dimsf[1] = {gloabl_size};
    hsize_t count[1] = {var_size};
    hsize_t start[1] = {offset};
    timer.start();
    hid_t file_id = H5Fopen(out_fin.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t dataset_id = H5Dopen2(file_id, (var_name).c_str(), H5P_DEFAULT);
    //std::cout << "- Dataset opened" << std::endl;
    // H5Dset_extent(dataset, dimsf);
    hid_t fspace_id = H5Dget_space(dataset_id);
    //std::cout << "- Dataset get" << std::endl;
    H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, start, NULL, count, NULL);
    hid_t mspace_id = H5Screate_simple(1, count, NULL);
    #ifdef HDF5_USE_MPI
    hid_t plistId_tran = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plistId_tran, H5FD_MPIO_COLLECTIVE);

    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    // timer.start();
    #ifdef HDF5_USE_MPI
    // H5Dread(dataset, H5T_NATIVE_UCHAR, mspace, fspace, plistId, var_data);
    H5Drefresh(dataset_id);
    H5Dread(dataset_id, H5T_NATIVE_UCHAR, mspace_id, fspace_id, H5P_DEFAULT, var_data);
    H5Pclose(plistId_tran);
    #else
    //H5Dread(dataset, H5T_NATIVE_UCHAR, H5S_ALL, fspace, H5P_DEFAULT, var_data);
      H5Dread(dataset_id, H5T_NATIVE_UCHAR, mspace_id, fspace_id, H5P_DEFAULT, var_data);
    #endif
    H5Dclose(dataset_id);
    H5Sclose(fspace_id);
    H5Sclose(mspace_id);
    H5Fclose(file_id);
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
      cudaEvent_t timer1, end;
      cudaEventCreate(&timer1);
      cudaEventCreate(&end);
      if (use_compression) {
        // ********* Decompression ********** //
        cusz_framework *d_framework = cusz_default_framework();
        cusz_compressor *d_comp = cusz_create(d_framework, FP32);
        cusz_config*     d_config     = new cusz_config{.eb = tol, .mode = Abs};
        T * decompressed_data = (T *) new T[original_size];
        T *d_decompressed_data;
        cudaMalloc((void**)&d_decompressed_data, sizeof(T) * original_size);
        unsigned char * d_compressed_data1;
        cudaMalloc(&d_compressed_data1, compressed_size);
        cudaMemcpy(d_compressed_data1, compressed_data, compressed_size, cudaMemcpyHostToDevice);
        cudaStream_t d_stream;
        cudaStreamCreate(&d_stream);
        cudaDeviceSynchronize();

        #ifdef HDF5_USE_MPI
	      MPI_Barrier(MPI_COMM_WORLD);
        #endif
        timer.start();
        cudaEventRecord(timer1, 0);   
        cusz_decompress(d_comp, &header, d_compressed_data1, compressed_size, d_decompressed_data, 
        decomp_len,(void*)&decompress_timerecord, d_stream);

        cudaMemcpy(decompressed_data, d_decompressed_data, sizeof(T) * original_size,
             cudaMemcpyDeviceToHost);
        cudaFree(d_decompressed_data);
        cudaFree(d_compressed_data1);
        cusz_release(d_comp);
        cudaStreamDestroy(d_stream);
        #ifdef HDF5_USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
        #endif
        timer.end();
        cudaEventRecord(end, 0);
        cudaEventSynchronize(end);
        float dcompression_time_1;
        cudaEventElapsedTime(&dcompression_time_1, timer1, end);
        double decompress_time = timer.get();
        float max_dcompression_time_1;
         #ifdef HDF5_USE_MPI
        MPI_Reduce(&dcompression_time_1, &max_dcompression_time_1, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
        double max_decompress_time;
        MPI_Reduce(&decompress_time, &max_decompress_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank==0)
        cusz::TimeRecordViewer::view_decompression(&decompress_timerecord, original_size * sizeof(T)); 
        if (rank == 0) 
        {
          res.decomp = max_decompress_time;
          std::cout <<  "Decompression time: "<< max_decompress_time << "\n";
          std::cout <<  "Decompression1 time: "<< max_dcompression_time_1 / 1000 << "\n";
          std::cout <<  "Decompression throughput: "<< (double)original_size * sizeof(float) * size/max_decompress_time/1e9 << " GB/s.\n";
        }
        #else
        if (rank==0)
        cusz::TimeRecordViewer::view_decompression(&decompress_timerecord, original_size * sizeof(T)); 
        if (rank == 0) 
        {
          res.decomp = decompress_time;
          std::cout <<  "Decompression time: "<< decompress_time << "\n";
          //std::cout <<  "Decompression1 time: "<< dcompression_time_1 / 1000 << "\n";
          std::cout <<  "Decompression throughput: "<< (double)original_size * sizeof(float) * size/decompress_time/1e9 << " GB/s.\n";
        }
        #endif
          
        timer.clear();
        delete[](T *) decompressed_data;
      }
    }
     #ifdef HDF5_USE_MPI
     MPI_Barrier(MPI_COMM_WORLD);
     #endif

      delete[](T*) original_data;
      delete[](unsigned char *) compressed_data;
      //return 0;
      return res;
      }

bool run(int argc, char *argv[]) {
 // int rank, size;
  int num_sim_iterations = 1;

  int size, rank;
   #ifdef HDF5_USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #else
  size =1;
  rank = 0;
  #endif

 
  std::string input_file = get_arg(argc, argv, "-i");
  std::string output_file = get_arg(argc, argv, "-c");
  if (rank == 0) 
    std::cout <<  "original data: " << input_file << "\n";
  if (rank == 0) 
    std::cout <<  "compressed data: " << output_file << "\n";
  
  int DIM = get_arg_int(argc, argv, "-n");
  std::vector<size_t> shape = get_arg_dims(argc, argv, "-n");
  if (rank == 0) {
      std::string shape_string = "shape (";
      for (int d = 0; d < shape.size(); d++)
        shape_string = shape_string + std::to_string(shape[d]) + " ";
      shape_string = shape_string + ")";
    }
    double tol = get_arg_double(argc, argv, "-e");

    int compute_delay = 0;
    if (has_arg(argc, argv, "-k")) {
      compute_delay = get_arg_int(argc, argv, "-k");
      if (rank == 0) std::cout <<  "compute delay: " << compute_delay << "\n";
    }
 
  
	std::vector<double> write, read, comp, decomp, write_comp, read_comp;
  bool use_compression;
  result res_comp_wo;
  // std::vector<double> ebs = {1e14};
  std::vector<double> ebs = {16.73116, 1.673116, 0.01673116, 0.1673116};
  use_compression = false; 
  res_comp_wo = launch_compress<float>(
        DIM, input_file.c_str(), output_file.c_str(), shape, tol, 1, use_compression, compute_delay); 
  for (double tol : ebs) {
  result res_comp;
  if (rank = 0)
    std::cout << " The abs error is " << tol << std::endl;
  use_compression = true; 
  res_comp = launch_compress<float>(
      DIM, input_file.c_str(), output_file.c_str(), shape, tol, 1, use_compression, compute_delay); 
  }

  return true;
}


int main(int argc, char *argv[]) {
  int size;
  int rank;
  #ifdef HDF5_USE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #else
  rank = 0;
  size = 1;
  #endif
            
  run(argc, argv);
   #ifdef HDF5_USE_MPI
  MPI_Finalize();
  #endif
  return 0;
}
