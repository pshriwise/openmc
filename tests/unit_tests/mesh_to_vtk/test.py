import difflib
import filecmp
import numpy as np
from pathlib import Path

import openmc
import pytest

from tests.regression_tests import config

pytest.importorskip('vtk')

def full_path(f):
    return Path(__file__).parent.absolute() / f

def diff_file(file1, file2):
    with open(file1) as fh:
        f1_text = fh.readlines()
    with open(file2) as fh:
        f2_text = fh.readlines()
    diff_lines = difflib.unified_diff(f1_text, f2_text)
    return ''.join(diff_lines)

# test meshes
reg_mesh = openmc.RegularMesh()
reg_mesh.lower_left = (0, 0, 0)
reg_mesh.upper_right = (20, 50, 50)
reg_mesh.dimension = (10, 20, 30)

rect_mesh = openmc.RectilinearMesh()
rect_mesh.x_grid = np.linspace(0, 10, 5)
rect_mesh.y_grid = np.logspace(np.log10(5), np.log10(20), 10)
rect_mesh.z_grid = np.linspace(1, 100, 20)

cyl_mesh = openmc.CylindricalMesh()
cyl_mesh.r_grid = np.linspace(0, 5, 8)
cyl_mesh.phi_grid = np.linspace(0, 2 * np.pi, 4)
cyl_mesh.z_grid = np.linspace(0, 2, 6)

sphere_mesh = openmc.SphericalMesh()
sphere_mesh.r_grid = np.linspace(0, 5, 3)
sphere_mesh.theta_grid = np.linspace(0, 0.5 * np.pi, 4)
sphere_mesh.phi_grid = np.linspace(0, 2*np.pi, 8)
sphere_mesh.write_vtk_mesh("spherical-curvilinear.vtk")
sphere_mesh.write_vtk_mesh("spherical-linear.vtk", curvilinear=False)

def mesh_data(mesh_dims):
    data = 100 * np.arange(np.prod(mesh_dims))
    return data.reshape(*mesh_dims)

# update mesh files if needed
if config['update']:
    reg_mesh.write_vtk_mesh('regular.vtk')
    rect_mesh.write_vtk_mesh('rectilinear.vtk')
    cyl_mesh.write_vtk_mesh('cylindrical-linear.vtk', curvilinear=False)
    cyl_mesh.write_vtk_mesh('cylindrical-curvilinear.vtk')
    sphere_mesh.write_vtk_mesh('spherical-linear.vtk', curvilinear=False)
    sphere_mesh.write_vtk_mesh('spherical-curvilinear.vtk')

    data = {'ascending_data': mesh_data(cyl_mesh.dimension)}
    cyl_mesh.write_vtk_mesh('cyl-data.vtk', data=data)


test_data = ((reg_mesh, False, 'regular'),
             (rect_mesh, False, 'rectilinear'),
             (cyl_mesh, False, 'cylindrical-linear'),
             (cyl_mesh, True, 'cylindrical-curvilinear'),
             (sphere_mesh, False, 'spherical-linear'),
             (sphere_mesh, True, 'spherical-curvilinear'))

@pytest.mark.parametrize('mesh_params',
                         test_data,
                         ids=lambda params: params[2])
def test_mesh_write_vtk(mesh_params, run_in_tmpdir):
    mesh, curvilinear, filename = mesh_params

    test_data = full_path(filename + ".vtk")
    filename = filename + "-actual.vtk"
    # write the mesh file and compare to the expected version
    mesh.write_vtk_mesh(filename, curvilinear=curvilinear)

    try:
        assert filecmp.cmp(test_data, filename)
    except AssertionError as e:
        diff = diff_file(test_data, filename)
        raise AssertionError(diff) from e

# check data writing
def test_mesh_write_vtk_data(run_in_tmpdir):

    data = {'ascending_data': mesh_data(cyl_mesh.dimension)}
    filename_expected = full_path('cyl-data.vtk')
    filename_actual = 'cyl-data-actual.vtk'
    cyl_mesh.write_vtk_mesh(filename_actual, data=data)

    try:
        assert filecmp.cmp(filename_actual, filename_expected)
    except AssertionError as e:
        diff = diff_file(filename_expected, filename_actual)
        raise AssertionError(diff) from e



