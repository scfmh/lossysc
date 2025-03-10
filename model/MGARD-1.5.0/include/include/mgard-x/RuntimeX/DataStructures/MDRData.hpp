/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#ifndef MDR_X_MDR_DATA_HPP
#define MDR_X_MDR_DATA_HPP

#include "MDRMetadata.hpp"

namespace mgard_x {
namespace MDR {

template <typename DeviceType> class MDRData {
public:
  MDRData() {}
  // Data
  std::vector<std::vector<Array<1, Byte, DeviceType>>> compressed_bitplanes;

  MDRData(SIZE num_levels, SIZE num_bitplanes) {
    compressed_bitplanes.resize(num_levels);
    for (int level_idx = 0; level_idx < num_levels; level_idx++) {
      compressed_bitplanes[level_idx].resize(num_bitplanes);
    }
  }

  void Resize(SIZE num_levels, SIZE num_bitplanes) {
    if (compressed_bitplanes.size() != num_levels) {
      compressed_bitplanes.resize(num_levels);
    }
    for (int level_idx = 0; level_idx < num_levels; level_idx++) {
      if (compressed_bitplanes[level_idx].size() != num_bitplanes) {
        compressed_bitplanes[level_idx].resize(num_bitplanes);
      }
    }
  }

  void Resize(MDRMetadata &mdr_metadata) {
    compressed_bitplanes.resize(mdr_metadata.num_levels);
    for (int level_idx = 0; level_idx < mdr_metadata.num_levels; level_idx++) {
      compressed_bitplanes[level_idx].resize(mdr_metadata.num_bitplanes);
      for (int bitplane_idx = 0; bitplane_idx < mdr_metadata.num_bitplanes;
           bitplane_idx++) {
        compressed_bitplanes[level_idx][bitplane_idx].resize(
            {mdr_metadata.level_sizes[level_idx][bitplane_idx]});
      }
    }
  }

  void CopyFromRefactoredData(MDRMetadata &mdr_metadata,
                              std::vector<std::vector<Byte *>> &refactored_data,
                              int queue_idx) {
    Resize(mdr_metadata);
    for (int level_idx = 0; level_idx < mdr_metadata.num_levels; level_idx++) {
      for (int bitplane_idx =
               mdr_metadata.loaded_level_num_bitplanes[level_idx];
           bitplane_idx < mdr_metadata.requested_level_num_bitplanes[level_idx];
           bitplane_idx++) {
        //std::cout << "this is the level and bitplane index: " << level_idx << "and " << bitplane_idx << " ";
        MemoryManager<DeviceType>::Copy1D(
            compressed_bitplanes[level_idx][bitplane_idx].data(),
            refactored_data[level_idx][bitplane_idx],
            mdr_metadata.level_sizes[level_idx][bitplane_idx], queue_idx);
            //printf("\\x%.2x", compressed_bitplanes[level_idx][bitplane_idx].data());
      }
    } //compressed_bitplanes[level_idx][bitplane_idx].data()
    //std::cout << "this is the refactor_data in CopyFromRefactoredData. " << *(double *)refactored_data[1][0] << std::endl;
    mdr_metadata.DoneLoadingBitplans();
  }

  void CopyToRefactoredData(MDRMetadata &mdr_metadata,
                            std::vector<std::vector<Byte *>> &refactored_data,
                            int queue_idx) {
    refactored_data.resize(mdr_metadata.num_levels);
    for (int level_idx = 0; level_idx < mdr_metadata.num_levels; level_idx++) {
      refactored_data[level_idx].resize(mdr_metadata.num_bitplanes);
      for (int bitplane_idx = 0; bitplane_idx < mdr_metadata.num_bitplanes;
           bitplane_idx++) {
        MemoryManager<DeviceType>::MallocHost(
            refactored_data[level_idx][bitplane_idx],
            mdr_metadata.level_sizes[level_idx][bitplane_idx], queue_idx);
        MemoryManager<DeviceType>::Copy1D(
            refactored_data[level_idx][bitplane_idx],
            compressed_bitplanes[level_idx][bitplane_idx].data(),
            mdr_metadata.level_sizes[level_idx][bitplane_idx], queue_idx);
      }
    }
  }

  void VerifyLoadedBitplans(MDRMetadata &mdr_metadata) {
    // TODO: load
    std::cout << "this is verify load value" << std::endl;
    for (int level_idx = 0; level_idx < mdr_metadata.num_levels; level_idx++) {
      for (int bitplane_idx =
               mdr_metadata.prev_used_level_num_bitplanes[level_idx];
           bitplane_idx < mdr_metadata.loaded_level_num_bitplanes[level_idx];
           bitplane_idx++) {
        if (!compressed_bitplanes[level_idx][bitplane_idx]
                 .hasDeviceAllocation()) {
          log::err("Bitplane verification failed. level_idx(" +
                   std::to_string(level_idx) + ") bitplane_idx(" +
                   std::to_string(bitplane_idx) + ") No allocation found.\n");
          exit(-1);
        }
        if (compressed_bitplanes[level_idx][bitplane_idx].shape(0) !=
            mdr_metadata.level_sizes[level_idx][bitplane_idx]) {
          log::err("Bitplane verification failed. level_idx(" +
                   std::to_string(level_idx) + ") bitplane_idx(" +
                   std::to_string(bitplane_idx) + ") Size mismatch " +
                   std::to_string(
                       compressed_bitplanes[level_idx][bitplane_idx].shape(0)) +
                   " vs. " +
                   std::to_string(
                       mdr_metadata.level_sizes[level_idx][bitplane_idx]));
          exit(-1);
        }
      }
    }
  }
};

} // namespace MDR
} // namespace mgard_x

#endif
