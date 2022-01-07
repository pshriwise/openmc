from contextlib import contextmanager
from copy import copy
import os
from shutil import copy
import shutil
import tempfile


@contextmanager
def cdtemp(files=None):
    """Context manager to change to/return from a tmpdir.

    Parameters
    ----------

    files : Iterable of file locations (str or Path-like)
        Set of files to copy into the temporary directory
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        cwd = os.getcwd()
        for file in files:
            copy(file, tmpdir, follow_symlinks=True)
        try:
            os.chdir(tmpdir)
            yield
        finally:
            os.chdir(cwd)
