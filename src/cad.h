
#ifdef CAD

#ifndef CAD_H
#define CAD_H

#include "DagMC.hpp"
#include "thread_manager.hpp"
#include "openmc/cell.h"
#include "openmc/surface.h"

extern DagThreadManager* DTM;

extern "C" void load_cad_geometry_c();
extern "C" void free_memory_cad_c();

#endif // CAD_H

#endif
