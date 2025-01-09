from contextlib import contextmanager
import os
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

import openmc
from .checkvalue import PathLike


@contextmanager
def change_directory(working_dir: PathLike | None = None, *, tmpdir: bool = False):
    """Context manager for executing in a provided working directory

    Parameters
    ----------
    working_dir : path-like
        Directory to switch to.
    tmpdir : bool
        Whether to use a temporary directory instead of a specific working directory

    """
    orig_dir = Path.cwd()

    # Set up temporary directory if requested
    if tmpdir:
        tmp = TemporaryDirectory()
        working_dir = tmp.name
    elif working_dir is None:
        raise ValueError('Must pass working_dir argument or specify tmpdir=True.')

    working_dir = Path(working_dir)
    working_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(working_dir)
    try:
        yield
    finally:
        os.chdir(orig_dir)
        if tmpdir:
            tmp.cleanup()


def input_path(filename: PathLike) -> Path:
    """Return a path object for an input file based on global configuration

    Parameters
    ----------
    filename : PathLike
        Path to input file

    Returns
    -------
    pathlib.Path
        Path object

    """
    if openmc.config['resolve_paths']:
        return Path(filename).resolve()
    else:
        return Path(filename)



def rotation_matrix(phi, theta, psi, degrees=True):
    """Generate a rotation matrix from rotations around the x, y, and z axes

    Parameters
    ----------
    phi : float
        Rotation around the x-axis.
    theta : float
        Rotation around the y-axis.
    psi : float
        Rotation around the z-axis.
    degrees : bool
        Whether input vaalues are in degrees or radians.


    Returns
    -------
    numpy.ndarray
        A three by three NumPy array containing the rotation matrix.

    """

    if degrees:
        phi, theta, psi = np.array((phi, theta, psi)) * (-np.pi/180.)

    c3, s3 = np.cos(phi), np.sin(phi)
    c2, s2 = np.cos(theta), np.sin(theta)
    c1, s1 = np.cos(psi), np.sin(psi)
    return np.array([[c1*c2, c1*s2*s3 - c3*s1, s1*s3 + c1*c3*s2],
                    [c2*s1, c1*c3 + s1*s2*s3, c3*s1*s2 - c1*s3],
                    [-s2, c2*s3, c2*c3]])
