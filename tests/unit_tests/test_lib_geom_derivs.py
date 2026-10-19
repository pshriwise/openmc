from collections.abc import Mapping
from copy import deepcopy
import os
import xml.etree.ElementTree as ET

import numpy as np
import pytest

import openmc
import openmc.lib


@pytest.fixture(autouse=True)
def require_cross_sections():
    cross_sections = os.environ.get('OPENMC_CROSS_SECTIONS')
    if cross_sections is None:
        pytest.skip('OPENMC_CROSS_SECTIONS is not set')
    if not ET.parse(cross_sections).getroot().findall('library'):
        pytest.skip('OPENMC_CROSS_SECTIONS does not list any libraries')


@pytest.mark.parametrize('parameters', [('r',), ('x0', 'y0', 'z0', 'r')])
def test_geometry_derivative_results_and_reset(mpi_intracomm, parameters):
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

    geom_deriv = openmc.GeometricDerivative(tally, {surface: parameters}, id=1)
    model.geometry_derivatives = [geom_deriv]

    model.settings.batches = 3
    model.settings.inactive = 0
    model.settings.particles = 100
    model.settings.run_mode = 'fixed source'
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
        assert deriv.surface_ids == [surface.id]
        assert deriv.tally is openmc.lib.tallies[tally.id]

        openmc.lib.run(output=False)

        assert deriv.surface_parameters == {
            surface.id: list(geom_deriv.surface_parameters[surface])}
        assert deriv.results.shape == (len(parameters), 1, 10)
        assert deriv.mean.shape == (len(parameters), 1)
        assert np.all(np.isfinite(deriv.results))
        assert np.any(deriv.mean != 0.0)

        deriv.reset()
        assert np.all(deriv.results == 0.0)

        openmc.lib.run(output=False)
        openmc.lib.reset()
        assert np.all(deriv.results == 0.0)


@pytest.mark.parametrize('event_based', [False, True])
@pytest.mark.parametrize('reverse_surfaces', [False, True])
def test_multiple_surfaces(sphere_model, mpi_intracomm, event_based,
                          reverse_surfaces):
    """Combined derivatives must retain each surface's independent results."""
    openmc.reset_auto_ids()
    model = deepcopy(sphere_model)
    material = model.materials[0]
    material.set_density('g/cm3', 10.0)
    # Deliberately choose IDs that differ from the internal surface indices.
    inner = openmc.Sphere(surface_id=17, r=1.0)
    outer = openmc.Sphere(surface_id=42, r=2.0, boundary_type='vacuum')
    cells = [openmc.Cell(fill=material, region=-inner),
             openmc.Cell(fill=material, region=+inner & -outer)]
    model.geometry = openmc.Geometry(cells)
    model.settings.batches = 3
    model.settings.event_based = event_based
    model.settings.output = {'summary': False}
    tally = openmc.Tally(
        filters=[openmc.CellFilter(cells),
                 openmc.EnergyFilter([0.0, 1.0e6, 20.0e6])],
        scores=['flux', 'total'], estimator='collision')
    model.tallies = [tally]

    parameters = [(outer, ['r', 'x0']), (inner, ['z0', 'r'])]
    if reverse_surfaces:
        parameters.reverse()
    combined = openmc.GeometricDerivative(tally, dict(parameters))
    separate = [openmc.GeometricDerivative(tally, {surface: params})
                for surface, params in parameters]
    model.geometry_derivatives = [combined, *separate]

    with openmc.lib.TemporarySession(
        model, intracomm=mpi_intracomm, output=False
    ):
        derivative = openmc.lib.geometry_derivatives[combined.id]
        assert derivative.surface_ids == [s.id for s, _ in parameters]
        assert list(derivative.surface_parameters.items()) == [
            (s.id, params) for s, params in parameters]

        openmc.lib.run(output=False)

        assert derivative.results.shape == (4, 8, 10)
        assert np.all(np.isfinite(derivative.results))
        for i, ((_, params), single) in enumerate(
                zip(parameters, separate)):
            reference = openmc.lib.geometry_derivatives[single.id]
            np.testing.assert_allclose(
                derivative.results[2*i:2*i + 2], reference.results)
            assert np.any(reference.mean[params.index('r')] != 0.0)

        derivative.reset()
        assert np.all(derivative.results == 0.0)
        assert np.any(reference.mean != 0.0)
        openmc.lib.reset()
        assert all(np.all(d.results == 0.0)
                   for d in openmc.lib.geometry_derivatives.values())
