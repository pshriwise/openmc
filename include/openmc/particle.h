#ifndef OPENMC_PARTICLE_H
#define OPENMC_PARTICLE_H

//! \file particle.h
//! \brief Particle type

#include <cstdint>
#include <sstream>
#include <string>

#include "openmc/constants.h"
#include "openmc/memory.h" // for unique_ptr
#include "openmc/particle_data.h"
#include "openmc/position.h"
#include "openmc/random_lcg.h"
#include "openmc/tallies/filter_match.h"
#include "openmc/vector.h"

namespace openmc {

// Forward declare the Surface class for use in Particle::cross_vacuum_bc, etc.
class Surface;

/*
 * The Particle class encompasses data and methods for transporting particles
 * through their lifecycle. Its base class defines particle data layout in
 * memory. A more detailed description of the rationale behind this approach
 * can be found in particle_data.h.
 */

class Particle : public ParticleData {
public:
  //==========================================================================
  // Constructors

  Particle() = default;

  //! create a secondary particle
  //
  //! stores the current phase space attributes of the particle in the
  //! secondary bank and increments the number of sites in the secondary bank.
  //! \param wgt Weight of the secondary particle
  //! \param u Direction of the secondary particle
  //! \param E Energy of the secondary particle in [eV]
  //! \param type Particle type
  void create_secondary(double wgt, Direction u, double E, ParticleType type);

  //! initialize from a source site
  //
  //! initializes a particle from data stored in a source site. The source
  //! site may have been produced from an external source, from fission, or
  //! simply as a secondary particle.
  //! \param src Source site data
  void from_source(const SourceSite* src);

  // Coarse-grained particle events
  void event_calculate_xs();
  void event_advance();
  void event_delta_advance();
  void event_cross_surface();
  void event_collide();
  void event_revive_from_secondary();
  void event_death();

  //! Cross a surface and handle boundary conditions
  void cross_surface();

  //! Cross a vacuum boundary condition.
  //
  //! \param surf The surface (with the vacuum boundary condition) that the
  //!   particle struck.
  void cross_vacuum_bc(const Surface& surf);

  //! Cross a reflective boundary condition.
  //
  //! \param surf The surface (with the reflective boundary condition) that the
  //!   particle struck.
  //! \param new_u The direction of the particle after reflection.
  void cross_reflective_bc(const Surface& surf, Direction new_u);

  //! Cross a periodic boundary condition.
  //
  //! \param surf The surface (with the periodic boundary condition) that the
  //!   particle struck.
  //! \param new_r The position of the particle after translation/rotation.
  //! \param new_u The direction of the particle after translation/rotation.
  //! \param new_surface The signed index of the surface that the particle will
  //!   reside on after translation/rotation.
  void cross_periodic_bc(
    const Surface& surf, Position new_r, Direction new_u, int new_surface);

  //! mark a particle as lost and create a particle restart file
  //! \param message A warning message to display
  void mark_as_lost(const char* message);

  void mark_as_lost(const std::string& message)
  {
    mark_as_lost(message.c_str());
  }

  void mark_as_lost(const std::stringstream& message)
  {
    mark_as_lost(message.str());
  }

  //! create a particle restart HDF5 file
  void write_restart() const;
<<<<<<< HEAD
=======

  //! Gets the pointer to the particle's current PRN seed
  uint64_t* current_seed() {return seeds_ + stream_;}
  const uint64_t* current_seed() const {return seeds_ + stream_;}

  //==========================================================================
  // Data members

  // Cross section caches
  std::vector<NuclideMicroXS> neutron_xs_; //!< Microscopic neutron cross sections
  std::vector<ElementMicroXS> photon_xs_; //!< Microscopic photon cross sections
  MacroXS macro_xs_; //!< Macroscopic cross sections

  int64_t id_;  //!< Unique ID
  Type type_ {Type::neutron};   //!< Particle type (n, p, e, etc.)

  int n_coord_ {1};              //!< number of current coordinate levels
  int cell_instance_;            //!< offset for distributed properties
  std::vector<LocalCoord> coord_; //!< coordinates for all levels

  // Particle coordinates before crossing a surface
  int n_coord_last_ {1};      //!< number of current coordinates
  std::vector<int> cell_last_;  //!< coordinates for all levels

  // Energy data
  double E_;       //!< post-collision energy in eV
  double E_last_;  //!< pre-collision energy in eV
  int g_ {0};      //!< post-collision energy group (MG only)
  int g_last_;     //!< pre-collision energy group (MG only)

  // Other physical data
  double wgt_ {1.0};     //!< particle weight
  double mu_;      //!< angle of scatter
  bool alive_ {true};     //!< is particle alive?

  // Other physical data
  Position r_last_current_; //!< coordinates of the last collision or
                            //!< reflective/periodic surface crossing for
                            //!< current tallies
  Position r_last_;   //!< previous coordinates
  Direction u_last_;  //!< previous direction coordinates
  double wgt_last_ {1.0};   //!< pre-collision particle weight
  double wgt_absorb_ {0.0}; //!< weight absorbed for survival biasing

  // What event took place
  bool fission_ {false}; //!< did particle cause implicit fission
  TallyEvent event_;          //!< scatter, absorption
  int event_nuclide_;  //!< index in nuclides array
  int event_mt_;       //!< reaction MT
  int delayed_group_ {0};  //!< delayed group

  // Post-collision physical data
  int n_bank_ {0};        //!< number of fission sites banked
  int n_bank_second_ {0}; //!< number of secondary particles banked
  double wgt_bank_ {0.0}; //!< weight of fission sites banked
  int n_delayed_bank_[MAX_DELAYED_GROUPS];  //!< number of delayed fission
                                            //!< sites banked

  // Indices for various arrays
  int surface_ {0};         //!< index for surface particle is on
  int cell_born_ {-1};      //!< index for cell particle was born in
  int material_ {-1};       //!< index for current material
  int material_last_ {-1};  //!< index for last material

  // Boundary information
  BoundaryInfo boundary_;

  // Temperature of current cell
  double sqrtkT_ {-1.0};      //!< sqrt(k_Boltzmann * temperature) in eV
  double sqrtkT_last_ {0.0};  //!< last temperature

  // Statistical data
  int n_collision_ {0};  //!< number of collisions

  // Track output
  bool write_track_ {false};

  // Current PRNG state
  uint64_t seeds_[N_STREAMS]; // current seeds
  int      stream_;           // current RNG stream

  // Secondary particle bank
  std::vector<Particle::Bank> secondary_bank_;

  int64_t current_work_; // current work index

  std::vector<double> flux_derivs_;  // for derivatives for this particle

  std::vector<FilterMatch> filter_matches_; // tally filter matches

  std::vector<std::vector<Position>> tracks_; // tracks for outputting to file

  std::vector<NuBank> nu_bank_; // bank of most recently fissioned particles

  // Global tally accumulators
  double keff_tally_absorption_ {0.0};
  double keff_tally_collision_ {0.0};
  double keff_tally_tracklength_ {0.0};
  double keff_tally_leakage_ {0.0};

  bool trace_ {false};     //!< flag to show debug information

  double collision_distance_; // distance to particle's next closest collision

  int n_event_ {0}; // number of events executed in this particle's history

  // DagMC state variables
  #ifdef DAGMC
  moab::DagMC::RayHistory history_;
  Direction last_dir_;
  #endif

  int64_t n_progeny_ {0}; // Number of progeny produced by this particle
  bool delta_tracking_{false}; // !< Flag to indicate whether or not delta tracking is active
>>>>>>> Addging delta tracking flag to particle class.
};

//============================================================================
//! Functions
//============================================================================

std::string particle_type_to_str(ParticleType type);

ParticleType str_to_particle_type(std::string str);

} // namespace openmc

#endif // OPENMC_PARTICLE_H
