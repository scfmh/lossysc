cmake_minimum_required(VERSION 3.19)

project(CompressionSimulationStudy CXX C)
find_package(mgard REQUIRED)
find_package(ADIOS2 REQUIRED)
#add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)



add_executable(cpu-application-simulator cpu-application-simulator.cpp)
target_link_libraries(cpu-application-simulator mgard::mgard adios2::adios2 MPI::MPI_CXX)
add_executable(cpu-application-simulator_GCLDLWP cpu-application-simulator_GCLDLWP.cpp)
target_link_libraries(cpu-application-simulator_GCLDLWP mgard::mgard adios2::adios2 MPI::MPI_CXX)

find_package(CUDAToolkit REQUIRED)


find_package(ZFP REQUIRED)
set(ZFP_INCLUDE_DIRS ${ZFP_DIR}/../../../include)
set(ZFP_LIBRARIS ${ZFP_DIR}/../../../lib64)
set (CMAKE_SHARED_LINKER_FLAGS "-Wl,--copy-dt-needed-entries")

add_executable(ZFP_cpu_simulation_iter_single ZFP_cpu_simulation_iter_single.cpp)
target_link_libraries(ZFP_cpu_simulation_iter_single adios2::adios2 CUDA::cudart MPI::MPI_CXX ${ZFP_LIBRARIES})
add_executable(ZFP_cpu_simulation_iter_double_general ZFP_cpu_simulation_iter_double_general.cpp)
target_link_libraries(ZFP_cpu_simulation_iter_double_general adios2::adios2 CUDA::cudart MPI::MPI_CXX ${ZFP_LIBRARIES})

add_executable(ZFP_raw ZFP_raw.cpp)
target_link_libraries(ZFP_raw adios2::adios2 CUDA::cudart MPI::MPI_CXX ${ZFP_LIBRARIES})

find_package(CUDAToolkit REQUIRED)


find_package(CUSZ REQUIRED)
set(CUSZ_INCLUDE_DIRS ${CUSZ_DIR}/../../../include/cusz)
set(CUSZ_LIBRARIS ${CUSZ_DIR}/../../../lib64/libcusz.so ${CUSZ_DIR}/../../../lib64/libparszcomp.so ${CUSZ_DIR}/../../../lib64/libparszargp.so)


add_executable(SZ_cpu_simulation_cuda_iter_ICLDIWP  SZ_cpu_simulation_cuda_iter_ICLDIWP.cc)
target_include_directories(SZ_cpu_simulation_cuda_iter_ICLDIWP PRIVATE ${CUSZ_INCLUDE_DIRS})
target_link_libraries(SZ_cpu_simulation_cuda_iter_ICLDIWP  CUDA::cudart adios2::adios2 MPI::MPI_CXX ${CUSZ_LIBRARIS})

