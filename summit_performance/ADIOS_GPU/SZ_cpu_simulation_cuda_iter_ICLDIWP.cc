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

#include "cusz.h"
#include "cuszapi.hh"
#include "utils/io.hh"
#include "cli/timerecord_viewer.hh"

#define MYADIOS_USE_MPI 1
#include <adios2.h>
#if MYADIOS_USE_MPI
#include <mpi.h>
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

#if MYADIOS_USE_MPI
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
 rank = 0;
 size = 1;
#endif
 
/** ADIOS class factory of IO class objects */
#if MYADIOS_USE_MPI
 adios2::ADIOS adios(MPI_COMM_WORLD);
#else
 adios2::ADIOS adios;
#endif
 result res;
 adios2::IO adiosio = adios.DeclareIO("SimulationData");
 adiosio.SetEngine("BP4");

 adios2::Engine writer = adiosio.Open(output_file, adios2::Mode::Write);

 auto len = 1;
 for (int i = 0; i < DIM; i++)
   len *= shape[i];
  //std::cout << "origianl size is " << len << std::endl;
 //float* original_data = new float[original_size];
 /**********************read information *************/
 //readfile(input_file, original_data, original_size * sizeof(double));

  if (rank == 0) 
 {
   std::cout << "Data per rank: " << (double)len*sizeof(float)/1e9 << " GB\n";
 }

/********************************************************************/
    cusz_header header;
    uint8_t*    exposed_compressed;
    uint8_t* compressed;
    size_t      compressed_len;
  

    T *d_uncompressed, *h_uncompressed;
    T *d_decompressed, *h_decompressed;

    /* cuSZ requires a 3% overhead on device (not required on host). */
    size_t uncompressed_memlen = len * 1.03;
    size_t decompressed_memlen = uncompressed_memlen;

    // clang-format off
    cudaMalloc(     &d_uncompressed, sizeof(T) * uncompressed_memlen );
    cudaMallocHost( &h_uncompressed, sizeof(T) * len );
    cudaMalloc(     &d_decompressed, sizeof(T) * decompressed_memlen );
    cudaMallocHost( &h_decompressed, sizeof(T) * len );
    // clang-format on

    /* User handles loading from filesystem & transferring to device. */
    io::read_binary_to_array(input_file, h_uncompressed, len);
    cudaMemcpy(d_uncompressed, h_uncompressed, sizeof(T) * len, cudaMemcpyHostToDevice);

    /* a casual peek */
    //printf("peeking uncompressed data, 20 elements\n");
    //peek_devdata(d_uncompressed, 20);

    cudaStream_t stream;
    cudaStreamCreate(&stream);

    // using default
    //cusz_framework* framework = cusz_default_framework();
    // alternatively
    cusz_framework* framework = new cusz_custom_framework{
        .pipeline     = Auto,
        .predictor    = cusz_custom_predictor{.type = LorenzoI},
        .quantization = cusz_custom_quantization{.radius = 512},
        .codec        = cusz_custom_codec{.type = Huffman}};

    cusz_compressor* comp       = cusz_create(framework, FP32);
    cusz_config*     config     = new cusz_config{.eb = tol, .mode = Abs};
    //cusz_len         uncomp_len = cusz_len{3600, 1800, 26, 1, 1.03};
    //cusz_len         uncomp_len = cusz_len{3600, 1800, 26, 1, 1.03}; //SDRBENCH-EXAFEL-data-130x1480x1552.f32
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
    //cusz_len         uncomp_len = cusz_len{2338, 106711, 1, 1, 1.03}; //SDRBENCH-exaalt-helium/dataset2-2338x106711.x.f32.dat 109 979 2338
    cusz_len         decomp_len = uncomp_len;

    cusz::TimeRecord compress_timerecord;
    cusz::TimeRecord decompress_timerecord;
  
  /*************************************************************************/

 Timer timer, timer_total;
 cudaEvent_t timer1, end;
 cudaEventCreate(&timer1);
 cudaEventCreate(&end);
 MPI_Barrier(MPI_COMM_WORLD);
 //timer_total.start();
 for (int sim_iter = 0; sim_iter < num_sim_iterations; sim_iter++) 
 {
 // writer.BeginStep();
  
  sleep(compute_delay);
  if (use_compression) /***use compression or not ***********/
  {
   
  // ********* Compression ********** //
   #if MYADIOS_USE_MPI
   MPI_Barrier(MPI_COMM_WORLD);
   #endif
   timer.start();
   cudaEventRecord(timer1, 0);
    cusz_compress(comp, config, d_uncompressed, uncomp_len, &exposed_compressed, &compressed_len, &header,(void*)&compress_timerecord, stream);
   #if MYADIOS_USE_MPI
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
   cusz::TimeRecordViewer::view_compression(&compress_timerecord, len * sizeof(T), compressed_len);
   }
   #if MYADIOS_USE_MPI
   float max_compression_time_1;
   MPI_Reduce(&compression_time_1, &max_compression_time_1, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
   double max_compression_time;
   MPI_Reduce(&compression_time, &max_compression_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
   if (rank == 0)
   {
    res.comp = max_compression_time;
    std::cout  << "Compression time: " << max_compression_time << "\n";
    std::cout  << "Compression1 time: " << max_compression_time_1 / 1000<< "\n";
    std::cout  << "Compression throughput: " << (double)len * sizeof(T) * size/max_compression_time/1e9 << " GB/s.\n";
    std::cout  << "Compression ratio: "<< (double)len * sizeof(T) / (compressed_len) << "\n";
   }
   #else
   if (rank == 0)
   {
    res.comp = compression_time;
    std::cout  << "Compression time: " << compression_time << "\n";
    std::cout  << "Compression1 time: " << compression_time_1 / 1000<< "\n";
    std::cout  << "Compression throughput: " << (double)len * sizeof(T) * size/compression_time/1e9 << " GB/s.\n";
    std::cout  << "Compression ratio: "<< (double)len * sizeof(T) / (compressed_len) << "\n";
   }
   #endif
    timer.clear();
    cudaEventDestroy(timer1);
    cudaEventDestroy(end);
  } 


 // ********* Write********** //
 std::string var_name = "var_step" + std::to_string(sim_iter);
 size_t var_size = use_compression ? compressed_len : len * sizeof(T);
 unsigned char *var_data = use_compression ? (unsigned char *)exposed_compressed : (unsigned char *) d_uncompressed;
 adios2::Variable<unsigned char> simulation_var;
 simulation_var = adiosio.DefineVariable<unsigned char>(var_name,{var_size*size}, {var_size*rank}, {var_size});
 #if MYADIOS_USE_MPI
 MPI_Barrier(MPI_COMM_WORLD);
 #endif
 timer.start();
 writer.BeginStep();
 writer.Put<unsigned char>(simulation_var, var_data, adios2::Mode::Sync);
 writer.EndStep();
 writer.Close();
 //writer.Put<(unsigned char *)>(simulation_var, var_data, adios2::Mode::Sync);
 #if MYADIOS_USE_MPI
 MPI_Barrier(MPI_COMM_WORLD);
 #endif
 timer.end();
 double write_time = timer.get();
 #if MYADIOS_USE_MPI
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
 //writer.EndStep();

 }
 //writer.Close();
#if MYADIOS_USE_MPI
 MPI_Barrier(MPI_COMM_WORLD);
 #endif
 adios2::Engine reader = adiosio.Open(output_file, adios2::Mode::Read);
 for (int sim_iter = 0; sim_iter < num_sim_iterations; sim_iter++) 
 {
  sleep(compute_delay);

  //reader.BeginStep(); 
  std::string var_name = "var_step" + std::to_string(sim_iter);
  // ********* Read ********** //
  std::vector<unsigned char> var_data_vec;
  adios2::Variable<unsigned char> simulation_var;
  simulation_var = adiosio.InquireVariable<unsigned char>(var_name);
  adios2::Dims var_shape = simulation_var.Shape();
  size_t var_size = var_shape[0] / size;
  adios2::Box<adios2::Dims> sel({var_size*rank}, {var_size});
  simulation_var.SetSelection(sel);
  #if MYADIOS_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif
  timer.start();
  reader.BeginStep();
  reader.Get<unsigned char>(simulation_var, var_data_vec, adios2::Mode::Sync);
  reader.EndStep();
  reader.Close();
  #if MYADIOS_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif
  timer.end();
  double read_time = timer.get();
  #if MYADIOS_USE_MPI
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
    //compressed_size = var_shape[0] / size;
    //memcpy(compressed_data, var_data_vec.data(), compressed_size);
    timer.start();
    cudaMalloc(&compressed, compressed_len);
    cudaMemcpy(compressed, exposed_compressed, compressed_len, cudaMemcpyDeviceToDevice);
    #if MYADIOS_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    //timer.start();
    cudaEventRecord(timer1, 0);   
    cusz_decompress(
            comp, &header, exposed_compressed, compressed_len, d_decompressed, decomp_len,
            (void*)&decompress_timerecord, stream);
    #if MYADIOS_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.end();
    cudaEventRecord(end, 0);
    cudaEventSynchronize(end);
    float dcompression_time_1;
    cudaEventElapsedTime(&dcompression_time_1, timer1, end);
    double decompress_time = timer.get();
    #if MYADIOS_USE_MPI
    float max_dcompression_time_1;
    MPI_Reduce(&dcompression_time_1, &max_dcompression_time_1, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
    double max_decompress_time;
    MPI_Reduce(&decompress_time, &max_decompress_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank==0)
    cusz::TimeRecordViewer::view_decompression(&decompress_timerecord, len * sizeof(T)); 
    if (rank == 0) 
    {
      res.decomp = max_decompress_time;
      std::cout <<  "Decompression time: "<< max_decompress_time << "\n";
      std::cout <<  "Decompression1 time: "<< max_dcompression_time_1 / 1000 << "\n";
      std::cout <<  "Decompression throughput: "<< (double)len * sizeof(float) * size/max_decompress_time/1e9 << " GB/s.\n";
    }
    #else
    if (rank==0)
    cusz::TimeRecordViewer::view_decompression(&decompress_timerecord, len * sizeof(T));
    if (rank == 0)
    {
      res.decomp = decompress_time;
      std::cout <<  "Decompression time: "<< decompress_time << "\n";
      std::cout <<  "Decompression1 time: "<< dcompression_time_1 / 1000 << "\n";
      std::cout <<  "Decompression throughput: "<< (double)len * sizeof(float) * size/decompress_time/1e9 << " GB/s.\n";
    }
    #endif
    timer.clear();
  }
  //reader.EndStep(); 
 }
 //reader.Close();
 #if MYADIOS_USE_MPI
 MPI_Barrier(MPI_COMM_WORLD);
 #endif
 //delete[](double *) original_data;
 cusz_release(comp);
 cudaFree(compressed);
 cudaFreeHost(h_uncompressed);
 cudaFreeHost(h_decompressed);
 cudaFree(d_uncompressed);
 cudaFree(d_decompressed);
 // delete compressor;
 cudaStreamDestroy(stream);
 //delete[](unsigned char *) buffer;
 //delete[](double *) decompressed_data;
 //return 0;
 return res;
}

bool run(int argc, char *argv[]) {
 // int rank, size;
 int num_sim_iterations = 1;

/*#if MYADIOS_USE_MPI
 //int provided;
 // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
 //MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
 rank = 0;
 size = 1;
#endif*/

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

 
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
 
/*std::string dt = get_arg(argc, argv, "-t");
 if (dt.compare("s") == 0) {
    type = zfp_type_float;
    if (rank == 0) std::cout <<  "data type: Single precision\n";
  } else if (dt.compare("d") == 0) {
    type = zfp_type_double;
    if (rank == 0) std::cout << "data type: Double precision\n";
  } else
    if (rank == 0) print_usage_message("wrong data type.");*/

  
	std::vector<double> write, read, comp, decomp, write_comp, read_comp;
  bool use_compression;
  std::vector<double> ebs = {1e14};
  //for ( i = 0; i < 6; i++) {
  //for (double tol : ebs) { 
    //int i = 5;
  // if(accumulate_data[i] > max_accumulate_data) break; 
  //double result launch_compress(int DIM, const char *input_file, onst char *output_file, const char *cfgFile,std::vector<int> shape, 
 //double accuracy,int num_sim_iterations,bool use_compression, int compute_delay){
  /*result res_no_comp;
  use_compression = false;
    res_no_comp = launch_compress<float>(DIM, input_file.c_str(), output_file.c_str(), shape, tol, 1, use_compression, compute_delay);

  for (int i = 0; i < 4; i++) {
    write.push_back(res_no_comp.write);
    read.push_back(res_no_comp.read);     
  }*/
  //for (double tol : ebs) {
    result res_comp;
    if (rank = 0)
      std::cout << " The abs error is " << tol << std::endl;
    use_compression = true; 
    res_comp = launch_compress<float>(
        DIM, input_file.c_str(), output_file.c_str(), shape, tol, 1, use_compression, compute_delay); 
    /*comp.push_back(res_comp.comp);
    decomp.push_back(res_comp.decomp);
    write_comp.push_back(res_comp.write);
    read_comp.push_back(res_comp.read);*/
  //}

 /* if (rank == 0) {
  for (double time : write) std::cout << time << ", ";
  std::cout << "\n";
  for (double time : read) std::cout << time << ", ";
  std::cout << "\n";
  for (double time : comp) std::cout << time << ", ";
  std::cout << "\n";
  for (double time : decomp) std::cout << time << ", ";
  std::cout << "\n";
  for (double time : write_comp) std::cout << time << ", ";
  std::cout << "\n";
  for (double time : read_comp) std::cout << time << ", ";
  std::cout << "\n";
	}*/
  return true;
}


int main(int argc, char *argv[]) {
/*int size;
int rank;
#if MYADIOS_USE_MPI
 MPI_Init(&argc, &argv);
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
 rank = 0;
 size = 1;
#endif

//printf("rank %d of %d\n", rank, size);

  run(argc, argv);
  #if MYADIOS_USE_MPI
  MPI_Finalize();
  #endif*/
  //MPI_Finalize();
  MPI_Init(&argc, &argv);
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            
  //printf("rank %d of %d\n", rank, size);
  run(argc, argv);
  MPI_Finalize();
  return 0;
}
