import glob
from itertools import product
import os

import openmc
import openmc.lib
import numpy as np

import pytest
from tests.testing_harness import PyAPITestHarness
from tests.regression_tests import config

TETS_PER_VOXEL = 12

class UnstructuredMeshSourceTest(PyAPITestHarness):
    def __init__(self, statepoint_name, model, inputs_true, schemes):

        super().__init__(statepoint_name, model, inputs_true)
        self.schemes = schemes

    def _run_openmc(self):
        kwargs = {'openmc_exec' : 'openmc',
                  'event_based' : config['event'],
                  'tracks' : True}

        if config['mpi']:
            kwargs['mpi_args'] = [config['mpi'], '-n', config['mpi_np']]

        openmc.run(**kwargs)

    def _compare_results(self):
        # loop over the tracks and get data
        tracks = openmc.Tracks(filepath='tracks.h5')
        tracks_born = np.empty((len(tracks), 1))
        decimals = 0

        instances = np.zeros(1000)

        for i in range(0, len(tracks)):
            tracks_born[i] = tracks[i].particle_tracks[0].states['cell_id'][0]
            instances[int(tracks_born[i])-1] = instances[int(tracks_born[i])-1]+1

        if self.schemes == "file":
            assert(instances[0] > 0 and instances[1] > 0)
            assert(instances[0] > instances[1])

            for i in range(0, len(instances)):
                if i != 0 and i != 1:
                    assert(instances[i] == 0)

        else:
            assert(np.average(instances) == 10)
            assert(np.std(instances) < np.average(instances))
            assert(np.amax(instances) < 30)


    def _cleanup(self):
        super()._cleanup()
        output = glob.glob('tally*.vtk')
        output += glob.glob('tally*.e')
        for f in output:
            if os.path.exists(f):
                os.remove(f)




param_values = (['libmesh', 'moab'], # mesh libraries
                ['volume', 'file']) # Element weighting schemes

test_cases = []
# for i, (lib, holes) in enumerate(product(*param_values)):
for i, (lib, schemes) in enumerate(product(*param_values)):
    test_cases.append({'library' : lib,
                       'schemes' : schemes,
                       'inputs_true' : 'inputs_true{}.dat'.format(i)})


@pytest.mark.parametrize("test_opts", test_cases)
def test_unstructured_mesh(test_opts):

    openmc.reset_auto_ids()

    # skip the test if the library is not enabled
    if test_opts['library'] == 'moab' and not openmc.lib._dagmc_enabled():
        pytest.skip("DAGMC (and MOAB) mesh not enabled in this build.")

    if test_opts['library'] == 'libmesh' and not openmc.lib._libmesh_enabled():
        pytest.skip("LibMesh is not enabled in this build.")

    ### Materials ###
    materials = openmc.Materials()

    water_mat = openmc.Material(material_id=3, name="water")
    water_mat.add_nuclide("H1", 2.0)
    water_mat.add_nuclide("O16", 1.0)
    water_mat.set_density("atom/b-cm", 0.07416)
    materials.append(water_mat)

    materials.export_to_xml()

    ### Geometry ###
    dimen = 10
    size_hex = 20.0/dimen

    ### Geometry ###
    cell = np.empty((dimen, dimen, dimen), dtype=object)
    surfaces = np.empty((dimen+1, 3), dtype=object)

    geometry = openmc.Geometry()
    universe = openmc.Universe(universe_id=1, name="Contains all hexes")

    for i in range(0,dimen+1):
        surfaces[i][0] = openmc.XPlane(-10.0+i*size_hex, name="X plane at "+str(-10.0+i*size_hex))
        surfaces[i][1] = openmc.YPlane(-10.0+i*size_hex, name="Y plane at "+str(-10.0+i*size_hex))
        surfaces[i][2] = openmc.ZPlane(-10.0+i*size_hex, name="Z plane at "+str(-10.0+i*size_hex))

        surfaces[i][0].boundary_type = 'vacuum'
        surfaces[i][1].boundary_type = 'vacuum'
        surfaces[i][2].boundary_type = 'vacuum'

    for k in range(0,dimen):
        for j in range(0,dimen):
            for i in range(0,dimen):
                cell[i][j][k] = openmc.Cell(name=("x " + str(i) +" y " + str(j) +" z " + str(k)))
                cell[i][j][k].region = +surfaces[i][0] & -surfaces[i+1][0] & \
                                    +surfaces[j][1] & -surfaces[j+1][1] & \
                                    +surfaces[k][2] & -surfaces[k+1][2]
                cell[i][j][k].fill = water_mat
                universe.add_cell(cell[i][j][k])

    geometry = openmc.Geometry(universe)

    mesh_filename = "test_mesh_tets.e"

    uscd_mesh = openmc.UnstructuredMesh(mesh_filename, test_opts['library'])

    ### Tallies ###

    # create tallies
    tallies = openmc.Tallies()

    tally1 = openmc.Tally(1)
    tally1.scores = ['scatter', 'total', 'absorption']
    # Export tallies
    tallies = openmc.Tallies([tally1])
    tallies.export_to_xml()

    ### Settings ###
    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.particles = 10000
    settings.batches = 2

    # settings.create_fission_neutrons = False
    settings.tracks = [(1, 1, 1)]
    settings.max_tracks = 10000

    # source setup
    if test_opts['schemes'] == 'volume':
        space = openmc.stats.MeshIndependent(elem_weight_scheme=test_opts['schemes'], mesh=uscd_mesh)
    elif test_opts['schemes'] == 'file':
        array = np.zeros(12000)
        for i in range(0, 12):
            array[i] = 10
            array[i+12] = 2
        space = openmc.stats.MeshIndependent(elem_weight_scheme=test_opts['schemes'], weights_from_file=array, mesh=uscd_mesh)

    energy = openmc.stats.Discrete(x=[15.e+06], p=[1.0])
    source = openmc.Source(space=space, energy=energy)
    settings.source = source

    model = openmc.model.Model(geometry=geometry,
                               materials=materials,
                               tallies=tallies,
                               settings=settings)

    harness = UnstructuredMeshSourceTest('statepoint.2.h5',
                                   model,
                                   test_opts['inputs_true'],
                                   test_opts['schemes'])
    harness.main()
