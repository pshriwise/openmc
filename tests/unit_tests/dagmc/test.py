import numpy as np
import pytest

import openmc
import openmc.lib

from tests import cdtemp

pytestmark = pytest.mark.skipif(
    not openmc.lib.feature_enabled('dagmc'),
    reason="DAGMC CAD geometry is not enabled.")


@pytest.fixture(scope="module", autouse=True)
def dagmc_model(dagmc_models):
    model = dagmc_models.legacy_pincell
    model.settings.temperature = {'tolerance': 50.0}
    model.settings.verbosity = 1
    model.materials[0].temperature = 320
    dagmc_universe = model.geometry.root_universe

    # check number of surfaces and volumes for this pincell model there should
    # be 5 volumes: two fuel regions, water, graveyard, implicit complement (the
    # implicit complement cell is created automatically at runtime)
    # and 21 surfaces: 3 cylinders (9 surfaces) and a bounding cubic shell
    # (12 surfaces)
    assert dagmc_universe.n_cells == 5
    assert dagmc_universe.n_surfaces == 21

    with cdtemp():
        model.export_to_xml()
        openmc.lib.init()
        yield

    openmc.lib.finalize()


@pytest.mark.parametrize("cell_id,exp_temp", ((1, 320.0),   # assigned by material
                                              (2, 300.0),   # assigned in dagmc file
                                              (3, 293.6)))  # assigned by default
def test_dagmc_temperatures(cell_id, exp_temp):
    cell = openmc.lib.cells[cell_id]
    assert np.isclose(cell.get_temperature(), exp_temp)
