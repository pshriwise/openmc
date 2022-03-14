import filecmp
import numpy as np


import openmc
import pytest

pytest.importorskip('vtk')

reg_mesh = openmc.RegularMesh()
reg_mesh.lower_left = (0, 0, 0)
reg_mesh.upper_right = (20, 50, 50)
reg_mesh.dimension = (10, 20, 30)

reg_mesh.write_vtk_mesh("regular.vtk")

rect_mesh = openmc.RectilinearMesh()
rect_mesh.x_grid = np.linspace(0, 10, 5)
rect_mesh.y_grid = np.logspace(np.log10(5), np.log10(20), 10)
rect_mesh.z_grid = np.linspace(1, 100, 20)
rect_mesh.write_vtk_mesh("rectilinear.vtk")

cyl_mesh = openmc.CylindricalMesh()
cyl_mesh.r_grid = np.linspace(0, 5, 8)
cyl_mesh.phi_grid = np.linspace(0, 2 * np.pi, 4)
cyl_mesh.z_grid = np.linspace(0, 2, 6)
cyl_mesh.write_vtk_mesh("cylindrical-curvilinear.vtk")
cyl_mesh.write_vtk_mesh("cylindrical-linear.vtk", curvilinear=False)

sphere_mesh = openmc.SphericalMesh()
sphere_mesh.r_grid = np.linspace(0, 5, 3)
sphere_mesh.theta_grid = np.linspace(0, 0.5 * np.pi, 4)
sphere_mesh.phi_grid = np.linspace(0, 2*np.pi, 8)
sphere_mesh.write_vtk_mesh("spherical-curvilinear.vtk")
sphere_mesh.write_vtk_mesh("spherical-linear.vtk", curvilinear=False)

rand = 100 * np.arange(np.prod(cyl_mesh.dimension)).reshape(*cyl_mesh.dimension)

data = {'random_data': rand}
cyl_mesh.write_vtk_mesh('cyl.vtk', data=data)

def id_fn(params):
    return params[1]

test_data = ((cyl_mesh, "cylindrical-linear"),)

@pytest.mark.parametrize('mesh_params', test_data, ids=id_fn)
def test_mesh_write_vtk(mesh_params):
    mesh, filename = mesh_params

    test_data = filename + ".vtk"
    filename = filename + "-actual.vtk"
    # write the mesh file and compare to the expected version
    mesh.write_vtk_mesh(filename, curvilinear=False)

    assert filecmp.cmp(test_data, filename)




