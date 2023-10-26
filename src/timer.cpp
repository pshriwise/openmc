#ifdef _OPENMP
#include <omp.h>
#endif
#include "openmc/timer.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {

Timer time_active;
Timer time_bank;
Timer time_bank_sample;
Timer time_bank_sendrecv;
Timer time_finalize;
Timer time_inactive;
Timer time_initialize;
Timer time_read_xs;
Timer time_statepoint;
Timer time_tallies;
Timer time_total;
Timer time_transport;
Timer time_event_init;
Timer time_event_calculate_xs;
Timer time_event_advance_particle;
Timer time_event_surface_crossing;
Timer time_event_collision;
Timer time_event_death;
std::vector<Timer> time_find_cell;
std::vector<Timer> time_distance_to_boundary;
std::vector<Timer> time_cross_surface;

} // namespace simulation

//==============================================================================
// Timer implementation
//==============================================================================

void Timer::start()
{
  running_ = true;
  start_ = clock::now();
}

void Timer::stop()
{
  elapsed_ = elapsed();
  running_ = false;
}

void Timer::reset()
{
  running_ = false;
  elapsed_ = 0.0;
}

double Timer::elapsed()
{
  if (running_) {
    std::chrono::duration<double> diff = clock::now() - start_;
    return elapsed_ + diff.count();
  } else {
    return elapsed_;
  }
}

//==============================================================================
// Non-member functions
//==============================================================================

void reset_timers()
{
#ifdef _OPENMP
  int n_threads = omp_get_max_threads();
#else
  int n_threads = 1;
#endif

  simulation::time_find_cell.resize(n_threads);
  for (auto& t : simulation::time_find_cell) t.reset();
  simulation::time_distance_to_boundary.resize(n_threads);
  for (auto& t : simulation::time_distance_to_boundary) t.reset();
  simulation::time_cross_surface.resize(n_threads);
  for (auto& t : simulation::time_cross_surface) t.reset();

  simulation::time_active.reset();
  simulation::time_bank.reset();
  simulation::time_bank_sample.reset();
  simulation::time_bank_sendrecv.reset();
  simulation::time_finalize.reset();
  simulation::time_inactive.reset();
  simulation::time_initialize.reset();
  simulation::time_read_xs.reset();
  simulation::time_statepoint.reset();
  simulation::time_tallies.reset();
  simulation::time_total.reset();
  simulation::time_transport.reset();
  simulation::time_event_init.reset();
  simulation::time_event_calculate_xs.reset();
  simulation::time_event_advance_particle.reset();
  simulation::time_event_surface_crossing.reset();
  simulation::time_event_collision.reset();
  simulation::time_event_death.reset();
}

} // namespace openmc
