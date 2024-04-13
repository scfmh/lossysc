cmake_minimum_required(VERSION 3.19)

project(CompressionSimulationStudy CXX C)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
#set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -w -fopenmp -target x86_64-pc-linux-gnu -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx90a")

if(NOT SZ_FOUND)
  if((NOT SZ_ROOT) AND (NOT (ENV{SZ_ROOT} STREQUAL "")))
    set(SZ_ROOT "$ENV{SZ_ROOT}")
  endif()
  if(SZ_ROOT)
    set(SZ_INCLUDE_OPTS HINTS ${SZ_ROOT}/include NO_DEFAULT_PATHS)
    set(SZ_LIBRARY_OPTS
      HINTS ${SZ_ROOT}/lib ${SZ_ROOT}/lib64
      NO_DEFAULT_PATHS
    )
  endif()

  find_path(SZ_INCLUDE_DIR sz.h PATH_SUFFIXES sz ${SZ_INCLUDE_OPTS})
  find_library(SZ_LIBRARY NAMES SZ ${SZ_LIBRARY_OPTS})
  find_library(ZLIB_LIBRARY NAMES z zlib ${SZ_LIBRARY_OPTS})
  find_library(ZSTD_LIBRARY NAMES zstd ${SZ_LIBRARY_OPTS})

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(SZ
    FOUND_VAR SZ_FOUND
    REQUIRED_VARS SZ_LIBRARY ZLIB_LIBRARY ZSTD_LIBRARY SZ_INCLUDE_DIR
  )
  if(SZ_FOUND)
    set(SZ_INCLUDE_DIRS ${SZ_INCLUDE_DIR})
    set(SZ_LIBRARIES ${SZ_LIBRARY} ${ZLIB_LIBRARY} ${ZSTD_LIBRARY})
    if(SZ_FOUND AND NOT TARGET SZ::SZ)
      add_library(SZ::SZ UNKNOWN IMPORTED)
      set_target_properties(SZ::SZ PROPERTIES
        IMPORTED_LOCATION             "${SZ_LIBRARY}"
        INTERFACE_LINK_LIBRARIES      "${SZ_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${SZ_INCLUDE_DIR}"
      )
    endif()
  endif()
endif()

#find_package(mgard REQUIRED)
find_package(ADIOS2 REQUIRED)
#add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)

add_executable(SZ_cpu_simulation SZ_cpu_simulation.cpp)
target_link_libraries(SZ_cpu_simulation SZ::SZ adios2::adios2 MPI::MPI_CXX)


add_executable(SZ_cpu_simulation_ICLDIWP SZ_cpu_simulation_ICLDIWP.cpp)
target_link_libraries(SZ_cpu_simulation_ICLDIWP SZ::SZ adios2::adios2 MPI::MPI_CXX)

add_executable(SZ_cpu_simulation_para_iter_single SZ_cpu_simulation_para_iter_single.cpp)
target_link_libraries(SZ_cpu_simulation_para_iter_single SZ::SZ adios2::adios2 MPI::MPI_CXX)

add_executable(SZ_cpu_simulation_para_iter SZ_cpu_simulation_para_iter.cpp)
target_link_libraries(SZ_cpu_simulation_para_iter SZ::SZ adios2::adios2 MPI::MPI_CXX)


find_package(ZFP REQUIRED)
set(ZFP_INCLUDE_DIRS ${ZFP_DIR}/../../../include)
set(ZFP_LIBRARIS ${ZFP_DIR}/../../../lib64)
add_executable(ZFP_cpu_simulation_iter ZFP_cpu_simulation_iter.cpp)
target_link_libraries(ZFP_cpu_simulation_iter adios2::adios2 MPI::MPI_CXX ${ZFP_LIBRARIES})
add_executable(ZFP_cpu_simulation_iter_exaalt ZFP_cpu_simulation_iter_exaalt.cpp)
target_link_libraries(ZFP_cpu_simulation_iter_exaalt adios2::adios2 MPI::MPI_CXX ${ZFP_LIBRARIES})

find_package(mgard REQUIRED)
add_executable(cpu-application-simulator cpu-application-simulator.cpp)
target_link_libraries(cpu-application-simulator mgard::mgard adios2::adios2 MPI::MPI_CXX)
add_executable(cpu-application-simulator_exaalt cpu-application-simulator_exaalt.cpp)
target_link_libraries(cpu-application-simulator_exaalt mgard::mgard adios2::adios2 MPI::MPI_CXX)
