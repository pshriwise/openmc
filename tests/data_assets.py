"""Canonical locations of shared test data.

This module records the locations of resources used in regression and unit tests,
such as DAGMC geometries, unstructured meshes, and weight-window data. Tests
can access these paths through the session-scoped ``dagmc_files``, ``umesh_files``,
and ``ww_files`` fixtures provided here.

Paths are absolute so they remain valid when tests change working directories
or model operations relocate the XML. ``PyAPITestHarness._normalize_data_paths``
rewrites paths under ``TESTS_DIR`` as paths relative to that directory before
comparing inputs, keeping ``inputs_true.dat`` reference files portable.

Tests that require a local copy of these files, such as the external C++ DAGMC driver that
loads ``dagmc.h5m`` directly, handle copying and cleanup in their own fixtures.
"""

from pathlib import Path

TESTS_DIR = Path(__file__).parent

_DAGMC_DIR = TESTS_DIR / 'regression_tests' / 'dagmc'
_UMESH_DIR = TESTS_DIR / 'regression_tests' / 'unstructured_mesh'
_UNIT_DAGMC_DIR = TESTS_DIR / 'unit_tests' / 'dagmc'

# DAGMC geometries
DAGMC_LEGACY_PINCELL = _DAGMC_DIR / 'legacy' / 'dagmc.h5m'
DAGMC_REFL_PINCELL = _DAGMC_DIR / 'refl' / 'dagmc.h5m'
DAGMC_UNIVERSES = _DAGMC_DIR / 'universes' / 'dagmc.h5m'
DAGMC_UWUW_METADATA = _DAGMC_DIR / 'uwuw' / 'dagmc.h5m'
DAGMC_BROKEN_MODEL = _UNIT_DAGMC_DIR / 'broken_model.h5m'
DAGMC_SPHERE_RAD_5 = _UNIT_DAGMC_DIR / 'dagmc_sphere_r5.h5m'
DAGMC_TETS_NO_GRAVEYARD = _UNIT_DAGMC_DIR / 'dagmc_tetrahedral_no_graveyard.h5m'
DAGMC_NESTED_SHELLS = (TESTS_DIR / 'unit_tests' / 'weightwindows' / 'dagmc' /
                       'nested_shell_geometry.h5m')

# Unstructured meshes. Both the MOAB and libMesh readers accept the '.exo'
# extension directly, so the '.e' aliases these meshes used to carry are not
# needed.
UMESH_TETS_H5M = (TESTS_DIR / 'regression_tests' / 'external_moab' /
                   'test_mesh_tets.h5m')
UMESH_TETS_VTK = _UMESH_DIR / 'test_mesh_dagmc_tets.vtk'
UMESH_TETS_EXO = _UMESH_DIR / 'test_mesh_tets.exo'
UMESH_TETS_W_HOLES_EXO = _UMESH_DIR / 'test_mesh_tets_w_holes.exo'
UMESH_HEXES_EXO = _UMESH_DIR / 'test_mesh_hexes.exo'

# Weight window bounds shared between unit and regression tests
WW_N = TESTS_DIR / 'regression_tests' / 'weightwindows' / 'ww_n.txt'
WW_P = TESTS_DIR / 'regression_tests' / 'weightwindows' / 'ww_p.txt'

# Groupings behind the 'dagmc_files' / 'umesh_files' / 'ww_files' fixtures.
# The attribute name a test uses is the key here.
DAGMC_FILES = {
    'legacy': DAGMC_LEGACY_PINCELL,
    'refl': DAGMC_REFL_PINCELL,
    'universes': DAGMC_UNIVERSES,
    'uwuw': DAGMC_UWUW_METADATA,
    'broken': DAGMC_BROKEN_MODEL,
    'sphere_r5': DAGMC_SPHERE_RAD_5,
    'tets_no_graveyard': DAGMC_TETS_NO_GRAVEYARD,
    'nested_shells': DAGMC_NESTED_SHELLS,
}

UMESH_FILES = {
    'tets_moab': UMESH_TETS_H5M,
    'dagmc_tets': UMESH_TETS_VTK,
    'tets': UMESH_TETS_EXO,
    'tets_w_holes': UMESH_TETS_W_HOLES_EXO,
    'hexes': UMESH_HEXES_EXO,
}

WW_FILES = {
    'neutron': WW_N,
    'photon': WW_P,
}
