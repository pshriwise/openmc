#ifndef OPENMC_BANK_H
#define OPENMC_BANK_H

#include <cstdint>
#include <forward_list>
#include <mutex>

#include "openmc/openmp_interface.h"
#include "openmc/particle.h"
#include "openmc/position.h"
#include "openmc/shared_array.h"
#include "openmc/vector.h"

namespace openmc {

class SharedSecondaryBank;

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {

extern vector<SourceSite> source_bank;

extern SharedSecondaryBank shared_secondary_bank;

extern SharedArray<SourceSite> surf_source_bank;

extern SharedArray<SourceSite> fission_bank;

extern vector<int64_t> progeny_per_particle;

} // namespace simulation

//==============================================================================
// Classes
//==============================================================================

class SharedSecondaryBank {
public:
  using value_type = SourceSite;
  using const_iterator = std::forward_list<value_type>::const_iterator;

  void push_back(value_type secondary)
  {
    mutex_.lock();
    vec_.push_back(secondary);
    mutex_.unlock();
  }

  bool pop_back(value_type& x)
  {
    mutex_.lock();
    if (vec_.empty()) {
      mutex_.unlock();
      return false;
    }
    x = vec_.back();
    vec_.pop_back();
    mutex_.unlock();
    return true;
  }

  bool empty()
  {
    bool empty;
    mutex_.lock();
    empty = vec_.empty();
    mutex_.unlock();
    return empty;
  }

  void reserve(size_t size)
  {
    mutex_.lock();
    vec_.reserve(size);
    mutex_.unlock();
  }

  size_t size() {
    return vec_.size();
  }

private:
  std::vector<value_type> vec_;
  OpenMPMutex mutex_;
};

//==============================================================================
// Non-member functions
//==============================================================================

void sort_fission_bank();

void free_memory_bank();

void init_fission_bank(int64_t max);

void init_secondary_bank(int64_t max);

} // namespace openmc

#endif // OPENMC_BANK_H
