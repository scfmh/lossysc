cmake_minimum_required(VERSION 3.19)

project(CompressionSimulationStudy CXX C)
find_package(mgard REQUIRED)
find_package(MPI REQUIRED)
#add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)

#add_executable(hello_jobstep hello_jobstep.cpp)
#target_link_libraries(hello_jobstep mgard::mgard adios2::adios2 MPI::MPI_CXX)


add_executable(qoi_average qoi_average.cpp)
target_link_libraries(qoi_average mgard::mgard MPI::MPI_CXX)

add_executable(qoi_linear qoi_linear.cpp)
target_link_libraries(qoi_linear mgard::mgard MPI::MPI_CXX)
