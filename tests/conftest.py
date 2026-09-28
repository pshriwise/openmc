import os
import hashlib

import pytest
import openmc
import openmc.lib

from tests import data_assets
from tests.regression_tests import config as regression_config

# MD5 hash of the official NNDC HDF5 cross_sections.xml file.
# Generated via: md5sum /path/to/nndc_hdf5/cross_sections.xml
_NNDC_XS_MD5 = "2d00773012eda670bc9f95d96a31c989"

# Collected during pytest_configure, displayed at start and end of session
_environment_warnings = []


def _check_build_environment():
    """Check STRICT_FP and cross section data, collecting any warnings."""
    if not openmc.lib.feature_enabled('strict_fp'):
        _environment_warnings.append(
            "OpenMC was NOT built with -DOPENMC_ENABLE_STRICT_FP=on. "
            "Regression test results may not match reference values due to "
            "compiler floating-point optimizations. Rebuild with "
            "-DOPENMC_ENABLE_STRICT_FP=on for reproducible results."
        )

    xs_path = os.environ.get("OPENMC_CROSS_SECTIONS")
    if not xs_path:
        _environment_warnings.append(
            "OPENMC_CROSS_SECTIONS environment variable is not set. "
            "Regression tests require the NNDC HDF5 cross section data."
        )
    elif not os.path.isfile(xs_path):
        _environment_warnings.append(
            f"OPENMC_CROSS_SECTIONS ({xs_path}) is not a valid file path. "
            "Regression tests require the NNDC HDF5 cross section data."
        )
    else:
        with open(xs_path, "rb") as f:
            md5 = hashlib.md5(f.read()).hexdigest()
        if md5 != _NNDC_XS_MD5:
            _environment_warnings.append(
                f"OPENMC_CROSS_SECTIONS ({xs_path}) does not match the "
                "official NNDC HDF5 dataset. Regression tests expect the "
                "NNDC data; results may differ with other cross section "
                "libraries."
            )


def pytest_addoption(parser):
    parser.addoption('--exe')
    parser.addoption('--mpi', action='store_true')
    parser.addoption('--mpiexec')
    parser.addoption('--mpi-np')
    parser.addoption('--update', action='store_true')
    parser.addoption('--build-inputs', action='store_true')
    parser.addoption('--event', action='store_true')


def pytest_configure(config):
    opts = ['exe', 'mpi', 'mpiexec', 'mpi_np', 'update', 'build_inputs', 'event']
    for opt in opts:
        if config.getoption(opt) is not None:
            regression_config[opt] = config.getoption(opt)

    _check_build_environment()


def _print_environment_warnings(terminalreporter):
    """Print environment warnings as a visible section."""
    if _environment_warnings:
        terminalreporter.section("OpenMC Environment Warnings")
        for msg in _environment_warnings:
            terminalreporter.line(f"WARNING: {msg}", yellow=True)
        terminalreporter.line("")


def pytest_sessionstart(session):
    """Print environment warnings at the start of the test session."""
    _print_environment_warnings(session.config.pluginmanager.get_plugin(
        "terminalreporter"))


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    """Reprint environment warnings at the end so they aren't missed."""
    _print_environment_warnings(terminalreporter)


@pytest.fixture
def run_in_tmpdir(tmpdir):
    orig = tmpdir.chdir()
    try:
        yield
    finally:
        orig.chdir()

@pytest.fixture(scope="module")
def endf_data():
    return os.environ['OPENMC_ENDF_DATA']


class FileGroup:
    """Paths to a group of test data files, reached by attribute.

    Paths are absolute, so a test can use one regardless of the directory it
    runs in and regardless of whether OpenMC later relocates the model XML
    (``Model.convert_to_multigroup`` and ``Model.plot()`` both move it into a
    temporary directory). ``PyAPITestHarness._get_inputs`` rewrites these
    absolute paths to repo-relative ones before comparing against
    ``inputs_true.dat``, so reference files stay machine independent.
    """

    def __init__(self, files):
        self._files = files

    def __getattr__(self, name):
        if name.startswith('_'):
            raise AttributeError(name)
        try:
            return self._files[name]
        except KeyError:
            raise AttributeError(
                f"no test data file named '{name}'; available: "
                f"{', '.join(sorted(self._files))}") from None

    def __dir__(self):
        return [*super().__dir__(), *self._files]


def _file_group_fixture(group_name, files):
    """Build the '<group>_files' fixture.

    Session-scoped so that module- and session-scoped fixtures can use it; it
    only hands back constants, with nothing to set up or tear down.
    """
    @pytest.fixture(name=f'{group_name}_files', scope='session')
    def _fixture():
        return FileGroup(files)

    return _fixture


dagmc_files = _file_group_fixture('dagmc', data_assets.DAGMC_FILES)
umesh_files = _file_group_fixture('umesh', data_assets.UMESH_FILES)
ww_files = _file_group_fixture('ww', data_assets.WW_FILES)


class DAGMCModels:
    """DAGMC model builders exposed as attributes, each returning a fresh model."""

    def __init__(self, files):
        self._files = files

    @property
    def legacy_pincell(self):
        """A fuel-and-water pincell with a box source and a cell total tally."""
        openmc.reset_auto_ids()
        model = openmc.Model()
        model.settings.batches = 5
        model.settings.inactive = 0
        model.settings.particles = 100
        model.settings.source = openmc.IndependentSource(
            space=openmc.stats.Box([-4, -4, -4], [4, 4, 4]))

        model.geometry = openmc.Geometry(
            openmc.DAGMCUniverse(self._files.legacy))

        tally = openmc.Tally()
        tally.scores = ['total']
        tally.filters = [openmc.CellFilter(1)]
        model.tallies = [tally]

        fuel = openmc.Material(name='no-void fuel', material_id=40)
        fuel.add_nuclide('U235', 1.0, 'ao')
        fuel.set_density('g/cc', 11)

        water = openmc.Material(name='water', material_id=41)
        water.add_nuclide('H1', 2.0, 'ao')
        water.add_nuclide('O16', 1.0, 'ao')
        water.set_density('g/cc', 1.0)
        water.add_s_alpha_beta('c_H_in_H2O')

        model.materials = openmc.Materials([fuel, water])
        return model


@pytest.fixture(scope='session')
def dagmc_models(dagmc_files):
    """Shared DAGMC model builders for fixtures of any scope.

    Each attribute access creates an independent ``openmc.Model`` that callers
    can modify. Callers manage ID resets and library initialization/finalization.
    """
    return DAGMCModels(dagmc_files)


@pytest.fixture(scope='session', autouse=True)
def resolve_paths():
    with openmc.config.patch('resolve_paths', False):
        yield


@pytest.fixture(scope='session', autouse=True)
def disable_depletion_multiprocessing_under_mpi():
    """Fork-based depletion multiprocessing may deadlock if MPI is active."""
    if not regression_config['mpi']:
        yield
        return

    from openmc.deplete import pool

    original_setting = pool.USE_MULTIPROCESSING
    pool.USE_MULTIPROCESSING = False
    try:
        yield
    finally:
        pool.USE_MULTIPROCESSING = original_setting
