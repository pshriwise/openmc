from collections.abc import Mapping
import os
import xml.etree.ElementTree as ET

import numpy as np
import pytest

import openmc
import openmc.lib


def test_geometry_derivative_results_and_reset(mpi_intracomm):
    cross_sections = os.environ.get('OPENMC_CROSS_SECTIONS')
    if cross_sections is None:
        pytest.skip('OPENMC_CROSS_SECTIONS is not set')
    if not ET.parse(cross_sections).getroot().findall('library'):
        pytest.skip('OPENMC_CROSS_SECTIONS does not list any libraries')

    openmc.reset_auto_ids()

    material = openmc.Material()
    material.add_nuclide('U235', 1.0)
    material.set_density('g/cm3', 10.0)

    surface = openmc.Sphere(r=1.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=material, region=-surface)

    model = openmc.Model()
    model.geometry.root_universe = openmc.Universe(cells=[cell])
    model.materials.append(material)

    tally = openmc.Tally(tally_id=1)
    tally.estimator = 'collision'
    tally.filters = [openmc.CellFilter([cell])]
    tally.scores = ['flux']
    model.tallies.append(tally)

    geom_deriv = openmc.GeometricDerivative(tally, cell, id=1)
    model.geometry_derivatives = [geom_deriv]

    model.settings.batches = 3
    model.settings.inactive = 0
    model.settings.particles = 10
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((0.0, 0.0, 0.0))
    )
    model.settings.output = {'summary': False}

    with openmc.lib.TemporarySession(
        model, intracomm=mpi_intracomm, output=False
    ):
        geom_derivatives = openmc.lib.geometry_derivatives
        assert isinstance(geom_derivatives, Mapping)
        assert len(geom_derivatives) == 1

        deriv = geom_derivatives[1]
        assert isinstance(deriv, openmc.lib.GeometryDerivative)
        assert deriv.id == 1
        assert deriv.tally_id == tally.id
        assert deriv.cell_id == cell.id
        assert deriv.tally is openmc.lib.tallies[tally.id]
        assert deriv.cell is openmc.lib.cells[cell.id]

        openmc.lib.run(output=False)

        assert deriv.surface_ids == [surface.id]
        assert deriv.geom_parameters.shape == (1, 4)
        assert deriv.results.shape == (1, 1, 2)
        assert deriv.mean.shape == (1, 1)

        deriv.reset()
        assert np.all(deriv.results == 0.0)

        openmc.lib.run(output=False)
        openmc.lib.reset()
        assert np.all(deriv.results == 0.0)
