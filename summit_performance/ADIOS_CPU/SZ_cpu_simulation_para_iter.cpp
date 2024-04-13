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

#include "sz.h"


#define MYADIOS_USE_MPI 1
#include <adios2.h>
#if MYADIOS_USE_MPI
#include <mpi.h>
#endif

using namespace std::chrono;

class SZ_Timer {
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




void Usage()
{
    std::cout << "\n";
    std::cout << "USAGE:\n";
    std::cout << "./helloBPSZ Nx sz_accuracy\n";
    std::cout << "\t Nx: size of float and double arrays to be compressed\n";
    std::cout << "\t sz_accuracy: absolute accuracy e.g. 0.1, 0.001, to skip "
                 "compression: -1\n\n";
}

void print_usage_message(std::string error) {
  if (error.compare("") != 0) {
    std::cout  << error << std::endl;
  }
  printf("Options\n\
\t -z: compress data\n\
\t\t -i <path to data file to be compressed>\n\
\t\t -c <path to compressed file>\n\
\t\t -n <ndim>: total number of dimensions\n\
\t\t\t [dim1]: slowest dimention\n\
\t\t\t [dim2]: 2nd slowest dimention\n\
\t\t\t  ...\n\
\t\t\t [dimN]: fastest dimention\n\
\t\t -m <abs|rel>: error bound mode (abs: abolute; rel: relative)\n\
\t\t -e <error>: error bound\n\
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

std::vector<int> get_arg_dims(int argc, char *argv[], std::string option) {
  std::vector<int> shape;
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

void readfile(const char *input_file, double *in_buff, size_t read_size) {
  //std::cout <<  "Loading file: " << input_file << "\n";
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

void writefile(const char *output_file, size_t num_bytes, double *out_buff) {
  FILE *file = fopen(output_file, "w");
  fwrite(out_buff, 1, num_bytes, file);
  fclose(file);
}

/****#define SZ_SCES 0  //successful
#define SZ_NSCS -1 //Not successful
#define SZ_FERR -2 //Failed to open input file******/
struct result {
	double write;
  double read;
  double comp;
  double decomp;
};



result launch_compress(int DIM, const char *input_file, const char *output_file, std::vector<int> shape, 
 double accuracy,int num_sim_iterations,bool use_compression, int compute_delay){
 int rank, size;

#if MYADIOS_USE_MPI
//  int provided;
//  // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
//  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
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

if (rank == 0) { 
  std::cout << ">>> use_compression: " << use_compression << "\n";
	}

 result res;
 adios2::IO io = adios.DeclareIO("SimulationData");
 io.SetEngine("BP4");
 //int status = 0;
 //std::cout <<"This is SZ configuration " << cfgFile << std::endl;; 
 //status = SZ_Init(cfgFile.c_str());

 sz_params sz;
 memset(&sz, 0, sizeof(sz_params));
 sz.sol_ID = SZ;
 sz.sampleDistance = 100;
 sz.quantization_intervals = 0;
 sz.max_quant_intervals = 65536;
 sz.predThreshold = 0.99;
 sz.szMode = SZ_BEST_COMPRESSION;
 sz.losslessCompressor = ZSTD_COMPRESSOR;
 sz.gzipMode = 1;
 sz.errorBoundMode = ABS;
 sz.absErrBound = 0.1;
 sz.relBoundRatio = 1E-5;
 sz.psnr = 80.0;
 sz.pw_relBoundRatio = 1E-5;
 sz.segment_size = static_cast<int>(std::pow(5., static_cast<double>(DIM)));
 sz.pwr_type = SZ_PWR_MIN_TYPE;
 SZ_Init_Params(&sz);

 adios2::Engine writer = io.Open(output_file, adios2::Mode::Write);

 size_t original_size = 1;
 for (int i = 0; i < DIM; i++)
   original_size *= shape[i];
 double* original_data = new double[original_size];
 /**********************read information *************/
 readfile(input_file, original_data, original_size * sizeof(double));
 //std::cout << "the original value is "<< original_data[11] << std::endl;
 //std::cout << "this is the original size "<< (sizeof(original_data)) <<std::endl; 
 if (rank == 0) 
 {
   std::cout << "Data per rank: " << (double)original_size*sizeof(double)/1e9 << " GB\n";
 }
 //void *compressed_data = (void *)new (double)[original_size];
 unsigned char * compressed_data;
 size_t compressed_size = original_size * sizeof(double);
 void *decompressed_data = malloc(original_size * sizeof(double));
 SZ_Timer timer, timer_total;
 #if MYADIOS_USE_MPI
 MPI_Barrier(MPI_COMM_WORLD);
 #endif
 timer_total.start();
 for (int sim_iter = 0; sim_iter < num_sim_iterations; sim_iter++) 
 {
  //writer.BeginStep();
  
  sleep(compute_delay);
  if (use_compression) /***use compression or not ***********/
  {
   
  // ********* Compression ********** //
   #if MYADIOS_USE_MPI
   MPI_Barrier(MPI_COMM_WORLD);
   #endif
   timer.start();
   /********* SZ details needs to be re-wroten ***************/
   /****SZ_FLOAT 0 SZ_DOUBLE 1 ,r5 is the 5th demision *******/
   //unsigned char *compressed_data = SZ_compress(1, original_data, &compressed_size, 0, 0, shape[2], shape[1], shape[0]); 
   compressed_data = SZ_compress(SZ_DOUBLE, original_data, &compressed_size, 0, 0, shape[2], shape[1], shape[0]);
   SZ_Finalize();
   #if MYADIOS_USE_MPI
   MPI_Barrier(MPI_COMM_WORLD);
   #endif
   timer.end();
   double compression_time = timer.get();
   #if MYADIOS_USE_MPI
   double max_compression_time;
   MPI_Reduce(&compression_time, &max_compression_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
   if (rank == 0)
   {
    res.comp = max_compression_time;
    std::cout  << "Compression time: " << max_compression_time << "\n";
    std::cout  << "Compression throughput: " << (double)original_size * sizeof(double) * size/max_compression_time/1e9 << " GB/s.\n";
    std::cout  << "Compression ratio: "<< (double)original_size * sizeof(double) / (compressed_size) << "\n";
   }
   #else
   if (rank == 0)
   {
    res.comp = compression_time;
    std::cout  << "Compression time: " << compression_time << "\n";
    std::cout  << "Compression throughput: " << (double)original_size * sizeof(double) * size/compression_time/1e9 << " GB/s.\n";
    std::cout  << "Compression ratio: "<< (double)original_size * sizeof(double) / (compressed_size) << "\n";
   }

   #endif
    timer.clear();
  } 
 
 //SZ_Finalize(); 
 // ********* Write********** //
 std::string var_name = "var_step" + std::to_string(sim_iter);
 //size_t var_size = compressed_size;
 size_t var_size = use_compression ? compressed_size : original_size * sizeof(double);
 unsigned char *var_data = use_compression ? (unsigned char *)compressed_data : (unsigned char *) original_data;
 //unsigned char *var_data = (unsigned char *)compressed_data;
 adios2::Variable<unsigned char> simulation_var;
 simulation_var = io.DefineVariable<unsigned char>(var_name,{var_size*size}, {var_size*rank}, {var_size});
 #if MYADIOS_USE_MPI
 MPI_Barrier(MPI_COMM_WORLD);
 #endif
 timer.start();
 writer.BeginStep();
 writer.Put<unsigned char>(simulation_var, var_data, adios2::Mode::Sync);
 writer.EndStep();
 writer.Close();
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
 //std::cout << "compressed_size is "<<compressed_size <<std::endl;
 //std::cout << "this is compress data "<< compressed_data << std::endl;
 //std::cout << "this is compress data "<< sizeof(*compressed_data) / sizeof(unsigned char) << std::endl; 
 MPI_Barrier(MPI_COMM_WORLD);
 timer_total.end();
 double total_time = timer_total.get();
 #if MYADIOS_USE_MPI
 double max_total_time;
 MPI_Reduce(&total_time, &max_total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
 if (rank == 0) 
 {
  std::cout << "Total time: " << max_total_time << "\n";
 }

 MPI_Barrier(MPI_COMM_WORLD);
 #else
 if (rank == 0)
 {
  std::cout << "Total time: " << total_time << "\n";
 }
 #endif
 timer_total.start();
 adios2::Engine reader = io.Open(output_file, adios2::Mode::Read);
 for (int sim_iter = 0; sim_iter < num_sim_iterations; sim_iter++) 
 {
  sleep(compute_delay);

  //reader.BeginStep(); 
  std::string var_name = "var_step" + std::to_string(sim_iter);
  // ********* Read ********** //
  std::vector<unsigned char> var_data_vec;
  adios2::Variable<unsigned char> simulation_var;
  simulation_var = io.InquireVariable<unsigned char>(var_name);
  adios2::Dims var_shape = simulation_var.Shape();
  size_t var_size = var_shape[0] / size;
  adios2::Box<adios2::Dims> sel({var_size*rank}, {var_size});
  simulation_var.SetSelection(sel);
  #if MYADIOS_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif
  timer.start();
  reader.BeginStep();
  //reader.Get<unsigned char *>(simulation_var, var_data_vec, adios2::Mode::Sync);
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
  //std::cout << "this is the data vect " << var_data_vec.size()<< std::endl;
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
    #if MYADIOS_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.start();
    SZ_decompress(1, compressed_data,compressed_size, 0, 0, shape[2], shape[1], shape[0]); 
    SZ_Finalize();
    #if MYADIOS_USE_MPI    
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.end();
    double decompress_time = timer.get();
    #if MYADIOS_USE_MPI
    double max_decompress_time;
    MPI_Reduce(&decompress_time, &max_decompress_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) 
    {
      res.decomp = max_decompress_time;
      std::cout <<  "Decompression time: "<< max_decompress_time << "\n";
      std::cout <<  "Decompression throughput: "<< (double)original_size * sizeof(double) * size/max_decompress_time/1e9 << " GB/s.\n";
    }
    #else
    if (rank == 0)
    {
      res.decomp = decompress_time;
      std::cout <<  "Decompression time: "<< decompress_time << "\n";
      std::cout <<  "Decompression throughput: "<< (double)original_size * sizeof(double) * size/decompress_time/1e9 << " GB/s.\n";
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
 timer_total.end();
 total_time = timer_total.get();
 timer_total.clear();
 #if MYADIOS_USE_MPI
 max_total_time;
 MPI_Reduce(&total_time, &max_total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
 if (rank == 0) 
 {
  std::cout << "Total time: " << max_total_time << "\n"; 
 }
 #else
 if (rank == 0)
 {
  std::cout << "Total time: " << total_time << "\n";
 }

 #endif

 delete[](double *) original_data;
 //delete[](unsigned char *) compressed_data;
 delete[](double *) decompressed_data;
 //return 0;
 return res;
}


bool run(int argc, char *argv[]) {
  int rank, size;
 int num_sim_iterations = 1;

#if MYADIOS_USE_MPI
 //int provided;
 // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
 //MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
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
 std::vector<int> shape = get_arg_dims(argc, argv, "-n");

  double accuracy = get_arg_double(argc, argv, "-e");

  int compute_delay = 0;
  if (has_arg(argc, argv, "-k")) {
    compute_delay = get_arg_int(argc, argv, "-k");
    if (rank == 0) std::cout << "compute delay: " << compute_delay << "\n";
  }


  
	std::vector<double> write, read, comp, decomp, write_comp, read_comp;
  bool use_compression;
  std::vector<double> ebs = {0.1, 1e-4, 1e-7, 1e-10};
  //for ( i = 0; i < 6; i++) {
  //for (double tol : ebs) { 
    //int i = 5;
  // if(accumulate_data[i] > max_accumulate_data) break; 
  //double result launch_compress(int DIM, const char *input_file, onst char *output_file, const char *cfgFile,std::vector<int> shape, 
 //double accuracy,int num_sim_iterations,bool use_compression, int compute_delay){
  result res_no_comp;
  use_compression = false;
    res_no_comp = launch_compress(DIM, input_file.c_str(), output_file.c_str(), shape, accuracy, 1, use_compression, compute_delay);

  /*for (int i = 0; i < 4; i++) {
    write.push_back(res_no_comp.write);
    read.push_back(res_no_comp.read);     
  }
  for (double accuracy : ebs) {
    std::cout << "this is the error "<<accuracy<<std::endl;
    result res_comp;
    use_compression = true; 
    res_comp = launch_compress(
        DIM, input_file.c_str(), output_file.c_str(), shape, accuracy, 1, use_compression, compute_delay); 
    comp.push_back(res_comp.comp);
    decomp.push_back(res_comp.decomp);
    write_comp.push_back(res_comp.write);
    read_comp.push_back(res_comp.read);
  }

  if (rank == 0) {
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
int size;
int rank;
#if MYADIOS_USE_MPI
 MPI_Init(&argc, &argv);
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
 rank = 0;
 size = 1;
#endif
  run(argc, argv);
  #if MYADIOS_USE_MPI
  MPI_Finalize();
  #endif
  //MPI_Finalize();
  return 0;
}


