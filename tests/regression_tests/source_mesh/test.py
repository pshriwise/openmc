
import numpy as np
import pytest

import openmc
import openmc.model

def test_source_mesh():

    # build a simple void model
    min, max = -10, 10
    box = openmc.model.RectangularParallelepiped(min, max, min, max, min, max, boundary_type='vacuum')

    geometry = openmc.Geometry([openmc.Cell(region=-box)])

    settings = openmc.Settings()
    settings.particles = 100
    settings.batches = 10
    settings.run_mode = 'fixed source'

    # define a 2 x 2 x 2 mesh
    mesh = openmc.RegularMesh.from_domain(-box, (2, 2, 2))

    strengths = np.ones((2, 2, 2))

    energy = openmc.stats.Discrete([1.e6], [1.0])

    angle = openmc.stats.Monodirectional((-1, -1, -1))

    sources = [openmc.Source(energy=energy, angle=angle, strength=0.0) for _ in range(strengths.size)]

    # create only one non-zero strength source in voxel (0, 0, 0)
    sources[0].strength = 1.0

    mesh_source = openmc.MeshSource(mesh, sources=sources)

    settings.source = mesh_source

    model = openmc.Model(geometry=geometry, settings=settings)

    # tally the flux on the mesh
    mesh_filter = openmc.MeshFilter(mesh)
    tally = openmc.Tally()
    tally.filters = [mesh_filter]
    tally.scores = ['flux']

    model.tallies = openmc.Tallies([tally])

    model.run()

