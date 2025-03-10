#ifndef _MDR_REFACTOR_INTERFACE_HPP
#define _MDR_REFACTOR_INTERFACE_HPP

#include "../../RuntimeX/RuntimeX.h"
// #include "../DataStructures/MDRData.hpp"

namespace MDR {
namespace concepts {

// refactor: a general interface for scnetific data refactor
template <typename T> class RefactorInterface {
public:
  virtual ~RefactorInterface() = default;

  virtual void refactor(T const *data_, const std::vector<uint32_t> &dims,
                        uint8_t target_level, uint8_t num_bitplanes) = 0;

  virtual void write_metadata() const = 0;

  virtual void print() const = 0;
};
} // namespace concepts
} // namespace MDR

namespace mgard_x {
namespace MDR {
namespace concepts {

// refactor: a general interface for scnetific data refactor
template <DIM D, typename T_data, typename DeviceType> class RefactorInterface {
public:
  virtual ~RefactorInterface() = default;

  virtual void refactor(Array<D, T_data, DeviceType> &data_array,
                        uint8_t num_bitplanes, MDRMetadata &mdr_metadata,
                        MDRData<DeviceType> &mdr_data, int queue_idx) = 0;

  virtual void write_metadata(MDRMetadata &mdr_metadata) = 0;

  virtual void print() const = 0;
};
} // namespace concepts
} // namespace MDR
} // namespace mgard_x

#endif
