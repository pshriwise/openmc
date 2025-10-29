from pathlib import Path

import openmc
import openmc.lib

import h5py
import numpy as np
import pytest

from tests.testing_harness import PyAPITestHarness, config

pytestmark = pytest.mark.skipif(
    not openmc.lib._xdg_enabled(),
    reason="XDG geometry is not enabled.")

@pytest.fixture
def model():
    openmc.reset_auto_ids()

    model = openmc.Model()

    # settings
    model.settings.batches = 5
    model.settings.inactive = 0
    model.settings.particles = 100

    source_box = openmc.stats.Box([-4, -4, -4],
                                  [ 4,  4,  4])
    source = openmc.IndependentSource(space=source_box)

    model.settings.source = source

    # geometry
    xdg_mesh = openmc.XDGMesh(Path("pincell.h5m"), library="moab")
    xdg_univ = openmc.XDGUniverse(xdg_mesh)
    xdg_univ.type = 'surface_mesh'

    model.geometry = openmc.Geometry(xdg_univ)

    # tally
    tally = openmc.Tally()
    tally.scores = ['total']
    tally.filters = [openmc.CellFilter(1)]
    model.tallies = [tally]

    # materials
    u235 = openmc.Material(name="fuel")
    u235.add_nuclide('U235', 1.0, 'ao')
    u235.set_density('g/cc', 11)
    u235.id = 40

    water = openmc.Material(name="water")
    water.add_nuclide('H1', 2.0, 'ao')
    water.add_nuclide('O16', 1.0, 'ao')
    water.set_density('g/cc', 1.0)
    water.add_s_alpha_beta('c_H_in_H2O')
    water.id = 41

    mats = openmc.Materials([u235, water])
    model.materials = mats

    return model


def test_missing_material_name(model):
    # remove the first material, which is identified by name in the DAGMC file
    model.materials = model.materials[1:]
    with pytest.raises(RuntimeError) as exec_info:
        model.run()
    exp_error_msg = "Material with name/ID 'fuel' not found for volume (cell) 1"
    assert exp_error_msg in str(exec_info.value)


@pytest.mark.parametrize("library,filename", [('moab', 'pincell.h5m'), ('libmesh', 'pincell-implicit.exo')])
def test_pincell(model, library, filename):
    for u in model.geometry.get_all_universes().values():
        if isinstance(u, openmc.DAGMCUniverse):
            u.library = library
            u.filename = filename

    harness = PyAPITestHarness('statepoint.5.h5', model, inputs_true=f'inputs_{library}.dat')
    harness.main()
