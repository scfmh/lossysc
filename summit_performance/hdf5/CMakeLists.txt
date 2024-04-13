cmake_minimum_required(VERSION 3.19)

project(hdf5_compressor CXX C)

#set (CMAKE_CXX_STANDARD 17)
#set (CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(HDF5 REQUIRED)
find_package(MPI REQUIRED)

option(USE_CUDA "Enable CUDA support" OFF)
if(USE_CUDA)
  find_package(CUDA REQUIRED)
endif()

add_executable(write_1 write_1.cpp)
target_link_libraries(write_1 MPI::MPI_CXX CUDA::cudart ${ZFP_LIBRARIES} ${HDF5_LIBRARIES})
find_package(ZFP REQUIRED)
set(ZFP_INCLUDE_DIRS ${ZFP_DIR}/../../../include)
set(ZFP_LIBRARIS ${ZFP_DIR}/../../../lib64)
set (CMAKE_SHARED_LINKER_FLAGS "-Wl,--copy-dt-needed-entries")
if(USE_CUDA)
  find_package(CUDA REQUIRED)
  add_executable(ZFP_hdf5 ZFP_hdf5.cpp)
  target_link_libraries(ZFP_hdf5 MPI::MPI_CXX CUDA::cudart ${ZFP_LIBRARIES} ${HDF5_LIBRARIES})
  target_include_directories(ZFP_hdf5 PRIVATE ${HDF5_INCLUDE_DIRS})
  
  add_executable(ZFP_hdf5_double ZFP_hdf5_double.cpp)
  target_link_libraries(ZFP_hdf5_double MPI::MPI_CXX CUDA::cudart ${ZFP_LIBRARIES} ${HDF5_LIBRARIES})
  target_include_directories(ZFP_hdf5_double PRIVATE ${HDF5_INCLUDE_DIRS})
else()
  add_executable(ZFP_hdf5 ZFP_hdf5.cpp)
  target_link_libraries(ZFP_hdf5 MPI::MPI_CXX ${ZFP_LIBRARIES} ${HDF5_LIBRARIES})
  target_include_directories(ZFP_hdf5 PRIVATE ${HDF5_INCLUDE_DIRS})

  add_executable(ZFP_hdf5_double ZFP_hdf5_double.cpp)
  target_link_libraries(ZFP_hdf5_double MPI::MPI_CXX ${ZFP_LIBRARIES} ${HDF5_LIBRARIES})
  target_include_directories(ZFP_hdf5_double PRIVATE ${HDF5_INCLUDE_DIRS})
endif()

find_package(mgard REQUIRED)
add_executable(MGARD_hdf5 MGARD_hdf5.cpp)
target_link_libraries(MGARD_hdf5 mgard::mgard MPI::MPI_CXX ${HDF5_LIBRARIES})
target_include_directories(MGARD_hdf5 PRIVATE ${HDF5_INCLUDE_DIRS})

add_executable(MGARD_hdf5_double MGARD_hdf5_double.cpp)
target_link_libraries(MGARD_hdf5_double mgard::mgard MPI::MPI_CXX ${HDF5_LIBRARIES})
target_include_directories(MGARD_hdf5_double PRIVATE ${HDF5_INCLUDE_DIRS})

add_executable(MGARD_hdf5_special MGARD_hdf5_special.cpp)
target_link_libraries(MGARD_hdf5_special mgard::mgard MPI::MPI_CXX ${HDF5_LIBRARIES})
target_include_directories(MGARD_hdf5_special PRIVATE ${HDF5_INCLUDE_DIRS})


if(NOT SZ_FOUND)
	if((NOT SZ_ROOT) AND (NOT (ENV{SZ_ROOT} STREQUAL "")))
		set(SZ_ROOT "$ENV{SZ_ROOT}")
	endif()
	if(SZ_ROOT)
		set(SZ_INCLUDE_OPTS HINTS ${SZ_ROOT}/include NO_DEFAULT_PATHS)
		set(SZ_LIBRARY_OPTSHINTS ${SZ_ROOT}/lib ${SZ_ROOT}/lib64 NO_DEFAULT_PATHS)
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

add_executable(SZ_hdf5_cpu SZ_hdf5_cpu.cc)
target_link_libraries(SZ_hdf5_cpu ${HDF5_LIBRARIES} SZ::SZ MPI::MPI_CXX)
target_include_directories(SZ_hdf5_cpu PRIVATE ${HDF5_INCLUDE_DIRS})

#SZ test case 
find_package(CUSZ REQUIRED)
#list(APPEND CMAKE_MODULE_PATH "/home/forensics/qzl/cusz/model/cusz-latest/install/lib/cmake/CUSZ")
message(STATUS "CUSZ_INCLUDE_DIRS: ${CUSZ_INCLUDE_DIRS}")
message(STATUS "CUSZ_LIBRARIES: ${CUSZ_LIBRARIES}")
set(CUSZ_INCLUDE_DIRS ${CUSZ_DIR}/../../../include/cusz)
#set(CUSZ_LIBRARIES ${CUSZ_DIR}/../../../lib/libcusz.so ${CUSZ_DIR}/../../../lib/libparszcomp.so ${CUSZ_DIR}/../../../lib/libparszargp.so)
set(CUSZ_LIBRARIES ${CUSZ_DIR}/../../../lib64/libcusz.so ${CUSZ_DIR}/../../../lib64/libparszcomp.so)
message(STATUS "CUSZ_INCLUDE_DIRS: ${CUSZ_INCLUDE_DIRS}")
message(STATUS "CUSZ_LIBRARIES: ${CUSZ_LIBRARIES}")

# Define SZ executable
add_executable(SZ_hdf5_gpu SZ_hdf5_gpu.cc)

# Link your executable with the necessary libraries
target_link_libraries(SZ_hdf5_gpu
            ${HDF5_LIBRARIES}
            ${CUSZ_LIBRARIES}
            mgard::mgard
            CUDA::cudart
            MPI::MPI_CXX)
# Include necessary directories
target_include_directories(SZ_hdf5_gpu PRIVATE ${HDF5_INCLUDE_DIRS}  # Include HDF5 include directories 
                                             ${CUSZ_INCLUDE_DIRS} )
