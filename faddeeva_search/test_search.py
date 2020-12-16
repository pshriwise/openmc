from pathlib import Path
import sys

import openmc
import pytest

current_dir = Path(__file__).parent

datasets = ('endfb_71', 'endfb_80')

def test_inf_media():
    # setup and infinite media problem
    model = openmc.model.Model()

    # create a material
    mat = openmc.Material()
    mat.add_nuclide('O16', 1.0, 'ao')

    cell = openmc.Cell(fill=mat)
    cell.region = openmc.rectangular_prism(1E5, 1E5, boundary_type='reflective')

    model.materials = [mat]
    model.materials.cross_sections = str(current_dir / "endfb71_hdf5" / "cross_sections.xml")

    model.geometry = openmc.Geometry([cell])

    model.settings.run_mode = 'fixed source'
    model.settings.particles = 10
    model.settings.batches = 10
    model.settings.temperature['multipole'] = True
    model.settings.resonance_scattering['enable'] = True


    model.export_to_xml()
    openmc.run()