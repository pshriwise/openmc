#include "openmc/finalize.h"

#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/cmfd_solver.h"
#include "openmc/constants.h"
#include "openmc/cross_sections.h"
#include "openmc/dagmc.h"
#include "openmc/eigenvalue.h"
#include "openmc/event.h"
#include "openmc/geometry.h"
#include "openmc/geometry_aux.h"
#include "openmc/material.h"
#include "openmc/mesh.h"
#include "openmc/message_passing.h"
#include "openmc/nuclide.h"
#include "openmc/photon.h"
#include "openmc/plot.h"
#include "openmc/random_lcg.h"
#include "openmc/random_ray/random_ray_simulation.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/source.h"
#include "openmc/surface.h"
#include "openmc/tallies/tally.h"
#include "openmc/thermal.h"
#include "openmc/timer.h"
#include "openmc/volume_calc.h"
#include "openmc/weight_windows.h"
#include "openmc/simulation_manager.h"

#include "xtensor/xview.hpp"

namespace openmc {

void free_memory()
{
  free_memory_geometry();
  free_memory_surfaces();
  free_memory_material();
  free_memory_volume();
  free_memory_simulation();
  free_memory_photon();
  free_memory_settings();
  free_memory_thermal();
  library_clear();
  nuclides_clear();
  free_memory_source();
  free_memory_mesh();
  free_memory_tally();
  free_memory_bank();
  free_memory_plot();
  free_memory_weight_windows();
  if (mpi::master) {
    free_memory_cmfd();
  }
  if (global_simulation.event_based()) {
    free_event_queues();
  }
}

} // namespace openmc

using namespace openmc;

int openmc_finalize()
{
  if (simulation::initialized)
    openmc_simulation_finalize();

  // Clear results
  openmc_reset();

  // Reset timers
  reset_timers();

  // Reset global variables
  global_simulation.set_assume_separate(false);
  global_simulation.set_check_overlaps(false);
  global_simulation.set_confidence_intervals(false);
  global_simulation.set_create_fission_neutrons(true);
  global_simulation.set_create_delayed_neutrons(true);
  global_simulation.set_electron_treatment(ElectronTreatment::LED);
  global_simulation.set_delayed_photon_scaling(true);
  global_simulation.set_energy_cutoff(std::array<double, 2>{0.0, 1000.0});
  global_simulation.set_time_cutoff(std::array<double, 2>{INFTY, INFTY});
  global_simulation.set_entropy_on(false);
  global_simulation.set_event_based(false);
  global_simulation.set_gen_per_batch(1);
  global_simulation.set_legendre_to_tabular(true);
  global_simulation.set_legendre_to_tabular_points(-1);
  global_simulation.set_material_cell_offsets(true);
  global_simulation.set_max_lost_particles(10);
  global_simulation.set_max_order(0);
  global_simulation.set_max_particles_in_flight(100000);
  global_simulation.set_max_particle_events(1'000'000);
  global_simulation.set_max_history_splits(10'000'000);
  global_simulation.set_max_tracks(1000);
  global_simulation.set_max_write_lost_particles(-1);
  global_simulation.set_n_log_bins(8000);
  global_simulation.set_n_inactive(0);
  global_simulation.set_n_particles(-1);
  global_simulation.set_output_summary(true);
  global_simulation.set_output_tallies(true);
  global_simulation.set_particle_restart_run(false);
  global_simulation.set_path_cross_sections("");
  global_simulation.set_path_input("");
  global_simulation.set_path_output("");
  global_simulation.set_path_particle_restart("");
  global_simulation.set_path_sourcepoint("");
  global_simulation.set_path_statepoint("");
  global_simulation.set_photon_transport(false);
  global_simulation.set_reduce_tallies(true);
  global_simulation.set_rel_max_lost_particles(1.0e-6);
  global_simulation.set_res_scat_on(false);
  global_simulation.set_res_scat_method(ResScatMethod::rvs);
  global_simulation.set_res_scat_energy_min(0.01);
  global_simulation.set_res_scat_energy_max(1000.0);
  global_simulation.set_restart_run(false);
  global_simulation.set_run_CE(true);
  global_simulation.set_run_mode(RunMode::UNSET);
  global_simulation.set_source_latest(false);
  global_simulation.set_source_rejection_fraction(0.05);
  global_simulation.set_source_separate(false);
  global_simulation.set_source_write(true);
  global_simulation.set_ssw_cell_id(C_NONE);
  global_simulation.set_ssw_cell_type(SSWCellType::None);
  global_simulation.set_ssw_max_particles(0);
  global_simulation.set_ssw_max_files(1);
  global_simulation.set_survival_biasing(false);
  global_simulation.set_temperature_default(293.6);
  global_simulation.set_temperature_method(TemperatureMethod::NEAREST);
  global_simulation.set_temperature_multipole(false);
  global_simulation.set_temperature_range(std::array<double, 2>{0.0, 0.0});
  global_simulation.set_temperature_tolerance(10.0);
  global_simulation.set_trigger_on(false);
  global_simulation.set_trigger_predict(false);
  global_simulation.set_trigger_batch_interval(1);
  global_simulation.set_uniform_source_sampling(false);
  global_simulation.set_ufs_on(false);
  global_simulation.set_urr_ptables_on(true);
  global_simulation.set_verbosity(7);
  global_simulation.set_weight_cutoff(0.25);
  global_simulation.set_weight_survive(1.0);
  global_simulation.set_weight_windows_file("");
  global_simulation.set_weight_windows_on(false);
  global_simulation.set_write_all_tracks(false);
  global_simulation.set_write_initial_source(false);
  global_simulation.set_cmfd_run(false);
  global_simulation.set_source_write_surf_id({});

  simulation::keff = 1.0;
  simulation::need_depletion_rx = false;
  simulation::ssw_current_file = 1;
  simulation::total_gen = 0;

  simulation::entropy_mesh = nullptr;
  simulation::ufs_mesh = nullptr;

  data::energy_max = {INFTY, INFTY};
  data::energy_min = {0.0, 0.0};
  data::temperature_min = 0.0;
  data::temperature_max = INFTY;
  model::root_universe = -1;
  model::plotter_seed = 1;
  openmc::openmc_set_seed(DEFAULT_SEED);
  openmc::openmc_set_stride(DEFAULT_STRIDE);

  // Deallocate arrays
  free_memory();

#ifdef LIBMESH
  global_simulation.set_libmesh_init(nullptr);
#endif

  // Free all MPI types
#ifdef OPENMC_MPI
  if (mpi::source_site != MPI_DATATYPE_NULL) {
    MPI_Type_free(&mpi::source_site);
  }
#endif

  openmc_reset_random_ray();

  return 0;
}

int openmc_reset()
{

  model::universe_cell_counts.clear();
  model::universe_level_counts.clear();

  for (auto& t : model::tallies) {
    t->reset();
  }

  // Reset global tallies
  simulation::n_realizations = 0;
  xt::view(simulation::global_tallies, xt::all()) = 0.0;

  simulation::k_col_abs = 0.0;
  simulation::k_col_tra = 0.0;
  simulation::k_abs_tra = 0.0;
  simulation::k_sum = {0.0, 0.0};
  simulation::satisfy_triggers = false;

  global_simulation.set_cmfd_run(false);

  simulation::n_lost_particles = 0;

  return 0;
}

int openmc_reset_timers()
{
  reset_timers();
  return 0;
}

int openmc_hard_reset()
{
  // Reset all tallies and timers
  openmc_reset();
  reset_timers();

  // Reset total generations and keff guess
  simulation::keff = 1.0;
  simulation::total_gen = 0;

  // Reset the random number generator state
  openmc::openmc_set_seed(DEFAULT_SEED);
  openmc::openmc_set_stride(DEFAULT_STRIDE);
  return 0;
}
