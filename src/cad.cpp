
#ifdef CAD

#include "cad.h"

moab::DagMC* DAGMC;


void create_bounding_surface(moab::DagMC *DAGMC) {
  struct BBOX {
    double lower[3] = { openmc::INFTY,  openmc::INFTY,  openmc::INFTY};
    double upper[3] = {-openmc::INFTY, -openmc::INFTY, -openmc::INFTY};

    bool valid() {
      return ( lower[0] <= upper[0] &&
	       lower[1] <= upper[1] &&
	       lower[2] <= upper[2] );
    }
    
    void update(double x, double y, double z) {
      lower[0] = x < lower[0] ? x : lower[0];
      lower[1] = y < lower[1] ? y : lower[1];
      lower[2] = z < lower[2] ? z : lower[2];
      
      upper[0] = x > upper[0] ? x : upper[0];
      upper[1] = y > upper[1] ? y : upper[1];
      upper[2] = z > upper[2] ? z : upper[2];      
    }

    void update(double xyz[3]) {
      update(xyz[0], xyz[1], xyz[2]);
    }
  };

  int num_vols = DAGMC->num_entities(3);
  double vmin[3], vmax[3];
  BBOX box;
  moab::ErrorCode rval;
  for(int i = 0; i < num_vols; i++) {
    moab::EntityHandle vol = DAGMC->entity_by_index(3, i+1);

    rval = DAGMC->getobb(vol, vmin, vmax);
    MB_CHK_ERR_CONT(rval);

    box.update(vmin);
    box.update(vmax);
    
  }

  //create box elements in MOAB


  //start with vertices
  double vert_coords[3];
  std::vector<moab::EntityHandle> box_verts;
  moab::EntityHandle new_vert;
  
  moab::Interface* MBI = DAGMC->moab_instance();

  // walk low-z face, starting with lower-left corner
  vert_coords[0] = box.lower[0];
  vert_coords[1] = box.lower[1];
  vert_coords[2] = box.lower[2]; 
  rval = MBI->create_vertex(box.lower, new_vert);
  MB_CHK_ERR_CONT(rval);
  box_verts.push_back(new_vert);
  
  // update x
  vert_coords[0] = box.upper[0];
  rval = MBI->create_vertex(box.lower, new_vert);
  MB_CHK_ERR_CONT(rval);
  box_verts.push_back(new_vert);
  
  // update y
  vert_coords[1] = box.upper[1];
  rval = MBI->create_vertex(box.lower, new_vert);
  MB_CHK_ERR_CONT(rval);
  box_verts.push_back(new_vert);

  // revert x
  vert_coords[0] = box.lower[0];
  rval = MBI->create_vertex(box.lower, new_vert);
  MB_CHK_ERR_CONT(rval);
  box_verts.push_back(new_vert);

  // now walk the upper-z face, starting with upper-left corner
  vert_coords[0] = box.lower[0];
  vert_coords[1] = box.lower[1];
  vert_coords[2] = box.upper[2];
  rval = MBI->create_vertex(box.lower, new_vert);
  MB_CHK_ERR_CONT(rval);
  box_verts.push_back(new_vert);

  // update x
  vert_coords[0] = box.upper[0];
  rval = MBI->create_vertex(box.lower, new_vert);
  MB_CHK_ERR_CONT(rval);
  box_verts.push_back(new_vert);
  
  // update y
  vert_coords[1] = box.upper[1];
  rval = MBI->create_vertex(box.lower, new_vert);
  MB_CHK_ERR_CONT(rval);
  box_verts.push_back(new_vert);

  // revert x
  vert_coords[0] = box.lower[0];
  rval = MBI->create_vertex(box.lower, new_vert);
  MB_CHK_ERR_CONT(rval);
  box_verts.push_back(new_vert);

  // now we have 8 vertices to create triangles with
  moab::Range new_tris;
  moab::EntityHandle new_tri;
  moab::EntityHandle tri_conn[3];

  // lower-z face triangles
  tri_conn[0] = box_verts[1]; tri_conn[1] = box_verts[0]; tri_conn[2] = box_verts[2];
  rval = MBI->create_element(moab::MBTRI, tri_conn, 3, new_tri);
  MB_CHK_ERR_CONT(rval);
  new_tris.insert(new_tri);
  
  tri_conn[0] = box_verts[3]; tri_conn[1] = box_verts[2]; tri_conn[2] = box_verts[0];
  rval = MBI->create_element(moab::MBTRI, tri_conn, 3, new_tri);
  MB_CHK_ERR_CONT(rval);
  new_tris.insert(new_tri);

  // upper-z face triangles
  tri_conn[0] = box_verts[5]; tri_conn[1] = box_verts[6]; tri_conn[2] = box_verts[4];
  rval = MBI->create_element(moab::MBTRI, tri_conn, 3, new_tri);
  MB_CHK_ERR_CONT(rval);
  new_tris.insert(new_tri);

  tri_conn[0] = box_verts[7]; tri_conn[1] = box_verts[4]; tri_conn[2] = box_verts[6];
  rval = MBI->create_element(moab::MBTRI, tri_conn, 3, new_tri);
  MB_CHK_ERR_CONT(rval);
  new_tris.insert(new_tri);

  // lower-x face triangles
  tri_conn[0] = box_verts[0]; tri_conn[1] = box_verts[1]; tri_conn[2] = box_verts[4];
  rval = MBI->create_element(moab::MBTRI, tri_conn, 3, new_tri);
  MB_CHK_ERR_CONT(rval);
  new_tris.insert(new_tri);

  tri_conn[0] = box_verts[3]; tri_conn[1] = box_verts[4]; tri_conn[2] = box_verts[1];
  rval = MBI->create_element(moab::MBTRI, tri_conn, 3, new_tri);
  MB_CHK_ERR_CONT(rval);
  new_tris.insert(new_tri);

  // upper-x face triangles
  tri_conn[0] = box_verts[2]; tri_conn[1] = box_verts[3]; tri_conn[2] = box_verts[6];
  rval = MBI->create_element(moab::MBTRI, tri_conn, 3, new_tri);
  MB_CHK_ERR_CONT(rval);
  new_tris.insert(new_tri);

  tri_conn[0] = box_verts[7]; tri_conn[1] = box_verts[6]; tri_conn[2] = box_verts[3];
  rval = MBI->create_element(moab::MBTRI, tri_conn, 3, new_tri);
  MB_CHK_ERR_CONT(rval);
  new_tris.insert(new_tri);

  // lower-y face triangles
  tri_conn[0] = box_verts[0]; tri_conn[1] = box_verts[4]; tri_conn[2] = box_verts[3];
  rval = MBI->create_element(moab::MBTRI, tri_conn, 3, new_tri);
  MB_CHK_ERR_CONT(rval);
  new_tris.insert(new_tri);

  tri_conn[0] = box_verts[7]; tri_conn[1] = box_verts[3]; tri_conn[2] = box_verts[4];
  rval = MBI->create_element(moab::MBTRI, tri_conn, 3, new_tri);
  MB_CHK_ERR_CONT(rval);
  new_tris.insert(new_tri);

  // upper-y face triangles
  tri_conn[0] = box_verts[1]; tri_conn[1] = box_verts[2]; tri_conn[2] = box_verts[3];
  rval = MBI->create_element(moab::MBTRI, tri_conn, 3, new_tri);
  MB_CHK_ERR_CONT(rval);
  new_tris.insert(new_tri);

  tri_conn[0] = box_verts[6]; tri_conn[1] = box_verts[3]; tri_conn[2] = box_verts[2];
  rval = MBI->create_element(moab::MBTRI, tri_conn, 3, new_tri);
  MB_CHK_ERR_CONT(rval);
  new_tris.insert(new_tri);

  //create a new surface meshset
  moab::EntityHandle bounding_surf;
  rval = MBI->create_meshset(0, bounding_surf);
  MB_CHK_ERR_CONT(rval);

  moab::Tag geom_tag = DAGMC->geom_tag();
  moab::Tag sense_tag = DAGMC->

}

void load_cad_geometry_c()
{
  if(!DAGMC) {
    DAGMC = new moab::DagMC();
  }

  int32_t cad_univ_id = 0; // universe is always 0 for CAD
  
  moab::ErrorCode rval = DAGMC->load_file("dagmc.h5m");
  MB_CHK_ERR_CONT(rval);

  rval = DAGMC->init_OBBTree();
  MB_CHK_ERR_CONT(rval);

  std::vector< std::string > prop_keywords;
  prop_keywords.push_back("mat");

  std::map<std::string, std::string> ph;
  DAGMC->parse_properties(prop_keywords, ph, ":");
  MB_CHK_ERR_CONT(rval);
  
  // initialize cell objects
  openmc::n_cells = DAGMC->num_entities(3);
  for(int i = 0; i < openmc::n_cells; i++)
    {
      moab::EntityHandle vol_handle = DAGMC->entity_by_index(3, i+1);
      
      // set cell ids using global IDs
      openmc::CADCell* c = new openmc::CADCell();
      c->id = DAGMC->id_by_index(3, i+1);
      c->dagmc_ptr = DAGMC;
      c->universe = cad_univ_id; // set to zero for now

      c->fill = openmc::C_NONE;
      openmc::global_cells.push_back(c);
      openmc::cell_map[c->id] = c->id;

      // Populate the Universe vector and dict
      auto it = openmc::universe_map.find(cad_univ_id);
      if (it == openmc::universe_map.end()) {
	openmc::global_universes.push_back(new openmc::Universe());
	openmc::global_universes.back()-> id = cad_univ_id;
	openmc::global_universes.back()->cells.push_back(i);
	openmc::universe_map[cad_univ_id] = openmc::global_universes.size() - 1;
      }
      else {
	openmc::global_universes[it->second]->cells.push_back(i);
      }

      
      if(DAGMC->is_implicit_complement(vol_handle)) {
	// assuming implicit complement is void for now
        c->material.push_back(openmc::C_MATERIAL_VOID);	
	continue;
      }
      
      if(DAGMC->has_prop(vol_handle, "mat")){
	std::string mat_value;
	rval = DAGMC->prop_value(vol_handle, "mat", mat_value);
	MB_CHK_ERR_CONT(rval);	
	int mat_id = std::stoi(mat_value);
	c->material.push_back(mat_id);
      }
      else {
	std::cout << "Warning: volume without material found!" << std::endl;
      }
      
    }

  // initialize surface objects
  openmc::n_surfaces = DAGMC->num_entities(2);
  openmc::surfaces_c = new openmc::Surface*[openmc::n_surfaces];
  
  for(int i = 0; i < openmc::n_surfaces; i++)
    {
      // set cell ids using global IDs
      openmc::CADSurface* s = new openmc::CADSurface();
      s->id = DAGMC->id_by_index(2, i+1);
      s->dagmc_ptr = DAGMC;
      s->bc = openmc::BC_TRANSMIT;
      openmc::surfaces_c[i] = s;
      openmc::surface_map[s->id] = s->id;
    }

  return;
}

void free_memory_cad_c()
{
  delete DAGMC;
}

#endif
