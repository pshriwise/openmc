from pathlib import Path
import sys

import openmc
import pytest

current_dir = Path(__file__).parent

datasets = ('endfb71_hdf5', 'endfb80_hdf5')

# create a library from the cross sections XML (just ENDF/B 7.1 for now)
xs_xml = current_dir / datasets[0] / 'cross_sections.xml'
data_lib = openmc.data.DataLibrary.from_xml(xs_xml)


nuclides = []
for lib in data_lib.libraries:
    if lib['type'] in ('wmp', 'neutron'):
        nuclides.append(lib['materials'][0])


@pytest.mark.parametrize('nuclide', nuclides)
def test_inf_media(nuclide):
    # setup and infinite media problem
    model = openmc.model.Model()

    # create a material
    mat = openmc.Material()
    mat.add_nuclide(nuclide, 1.0, 'ao')

    cell = openmc.Cell(fill=mat)
    cell.region = openmc.rectangular_prism(1E5, 1E5, boundary_type='reflective')

    model.materials = [mat]
    xs_dir = str(current_dir / "endfb71_hdf5" / "cross_sections.xml")
    model.materials.cross_sections = xs_dir

    model.geometry = openmc.Geometry([cell])

    model.settings.run_mode = 'fixed source'
    model.settings.particles = 1000
    model.settings.batches = 10
    model.settings.temperature['multipole'] = True
    model.settings.temperature['default'] = 294.0
    model.settings.resonance_scattering['enable'] = True


    src_loc = openmc.stats.Point()
    src_e = openmc.stats.Discrete([10E6], [1.0])
    model.settings.source = openmc.Source(space=src_loc, energy=src_e)

    model.export_to_xml()
    openmc.run()

if __name__ == '__main__':
    test_inf_media('H1')