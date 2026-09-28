from pathlib import Path
import os
import shutil
import subprocess
import textwrap

import openmc
import openmc.lib
import pytest

from tests.regression_tests import config
from tests.testing_harness import PyAPITestHarness

pytestmark = pytest.mark.skipif(
    not openmc.lib.feature_enabled('dagmc'),
    reason="DAGMC is not enabled.")

# Test that an external DAGMC instance can be passed in through the C API

@pytest.fixture
def cpp_driver(request):
    """Compile the external source"""

    # Get build directory and write CMakeLists.txt file
    openmc_dir = Path(str(request.config.rootdir)) / 'build'
    with open('CMakeLists.txt', 'w') as f:
        f.write(textwrap.dedent("""
            cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
            project(openmc_cpp_driver CXX)
            add_executable(main main.cpp)
            find_package(OpenMC REQUIRED HINTS {})
            target_link_libraries(main OpenMC::libopenmc)
            target_compile_features(main PUBLIC cxx_std_14)
            set(CMAKE_CXX_FLAGS "-pedantic-errors")
            add_compile_definitions(OPENMC_DAGMC_ENABLED=1)
            """.format(openmc_dir)))

    # Create temporary build directory and change to there
    local_builddir = Path('build')
    local_builddir.mkdir(exist_ok=True)
    os.chdir(local_builddir)

    mpi_arg = "On" if config['mpi'] else "Off"

    try:
        # Run cmake/make to build the shared libary
        subprocess.run(['cmake', os.path.pardir, f'-DOPENMC_USE_MPI={mpi_arg}'], check=True)
        subprocess.run(['make'], check=True)
        os.chdir(os.path.pardir)

        yield "./build/main"

    finally:
        # Remove local build directory when test is complete
        shutil.rmtree('build')
        os.remove('CMakeLists.txt')

@pytest.fixture
def dagmc_h5m_in_cwd(dagmc_files):
    """Place the geometry alongside the test.

    Unlike every other DAGMC test, this one drives OpenMC from main.cpp, which
    calls ``load_file("dagmc.h5m")`` directly rather than going through the
    model XML. The C++ driver therefore needs the file in its working
    directory regardless of the path recorded in the XML.
    """
    dest = Path(__file__).parent / dagmc_files.legacy.name
    shutil.copy(dagmc_files.legacy, dest)
    try:
        yield dest
    finally:
        dest.unlink(missing_ok=True)


@pytest.fixture
def model(dagmc_models):
    model = dagmc_models.legacy_pincell
    model.settings.temperature['default'] = 293
    return model

class ExternalDAGMCTest(PyAPITestHarness):
    def __init__(self, executable, statepoint_name, model):
        super().__init__(statepoint_name, model)
        self.executable = executable

    def _run_openmc(self):
        """
        Just test if results generated with the external C++ API are
        self-consistent with the internal python API.
        We generate the "truth" results with the python API but
        the main test produces the results file by running the
        executable compiled from main.cpp. This future-proofs
        the test - we only care that the two routes are equivalent.
        """
        if config['update']:
            # Generate the results file with internal python API
            openmc.run(openmc_exec=config['exe'], event_based=config['event'])
        elif config['mpi']:
            mpi_args = [config['mpiexec'], '-n', config['mpi_np']]
            # Run main cpp executable with MPI
            openmc.run(openmc_exec=self.executable,
                       mpi_args=mpi_args,
                       event_based=config['event'])
        else:
            # Run main cpp executable
            openmc.run(openmc_exec=self.executable,
                       event_based=config['event'])

def test_external_dagmc(cpp_driver, model):
    harness = ExternalDAGMCTest(cpp_driver, 'statepoint.5.h5', model)
    harness.main()
