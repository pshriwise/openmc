from collections.abc import Mapping
from ctypes import c_double, c_int, c_int32, c_size_t, POINTER
from weakref import WeakValueDictionary

from numpy.ctypeslib import as_array

from openmc.exceptions import AllocationError, InvalidIDError
from . import _dll
from .cell import cells
from .core import _FortranObjectWithID
from .error import _error_handler
from .tally import tallies


__all__ = ['GeometryDerivative', 'geometry_derivatives']


_dll.openmc_get_geometry_derivative_index.argtypes = [
    c_int32, POINTER(c_int32)]
_dll.openmc_get_geometry_derivative_index.restype = c_int
_dll.openmc_get_geometry_derivative_index.errcheck = _error_handler
_dll.openmc_geometry_derivatives_size.restype = c_size_t
_dll.openmc_geometry_derivative_get_id.argtypes = [
    c_int32, POINTER(c_int32)]
_dll.openmc_geometry_derivative_get_id.restype = c_int
_dll.openmc_geometry_derivative_get_id.errcheck = _error_handler
_dll.openmc_geometry_derivative_get_tally_id.argtypes = [
    c_int32, POINTER(c_int32)]
_dll.openmc_geometry_derivative_get_tally_id.restype = c_int
_dll.openmc_geometry_derivative_get_tally_id.errcheck = _error_handler
_dll.openmc_geometry_derivative_get_cell_id.argtypes = [
    c_int32, POINTER(c_int32)]
_dll.openmc_geometry_derivative_get_cell_id.restype = c_int
_dll.openmc_geometry_derivative_get_cell_id.errcheck = _error_handler
_dll.openmc_geometry_derivative_get_surface_ids.argtypes = [
    c_int32, POINTER(POINTER(c_int32)), POINTER(c_size_t)]
_dll.openmc_geometry_derivative_get_surface_ids.restype = c_int
_dll.openmc_geometry_derivative_get_surface_ids.errcheck = _error_handler
_dll.openmc_geometry_derivative_parameters.argtypes = [
    c_int32, POINTER(POINTER(c_double)), POINTER(c_size_t*2)]
_dll.openmc_geometry_derivative_parameters.restype = c_int
_dll.openmc_geometry_derivative_parameters.errcheck = _error_handler
_dll.openmc_geometry_derivative_results.argtypes = [
    c_int32, POINTER(POINTER(c_double)), POINTER(c_size_t*3)]
_dll.openmc_geometry_derivative_results.restype = c_int
_dll.openmc_geometry_derivative_results.errcheck = _error_handler
_dll.openmc_geometry_derivative_reset.argtypes = [c_int32]
_dll.openmc_geometry_derivative_reset.restype = c_int
_dll.openmc_geometry_derivative_reset.errcheck = _error_handler


class GeometryDerivative(_FortranObjectWithID):
    """Geometry derivative stored internally.

    This class exposes a geometry derivative that is stored internally in the
    OpenMC library. To obtain a view of a geometry derivative with a given ID,
    use the :data:`openmc.lib.geometry_derivatives` mapping.

    Parameters
    ----------
    uid : int or None
        Unique ID of the geometry derivative.
    new : bool
        Geometry derivatives cannot be allocated through the C API. This
        argument must remain False when `index` is None.
    index : int or None
        Index in the geometry derivatives array.

    Attributes
    ----------
    cell : openmc.lib.Cell
        Cell to which the geometry derivative is applied.
    cell_id : int
        ID of the cell to which the geometry derivative is applied.
    geom_parameters : numpy.ndarray
        Geometry parameter accumulator data.
    id : int
        ID of the geometry derivative.
    mean : numpy.ndarray
        Batch-averaged geometry derivative results.
    results : numpy.ndarray
        Raw geometry derivative results. The last axis stores the current
        batch accumulator and accumulated sum.
    surface_ids : list of int
        IDs of surfaces used by the geometry derivative.
    tally : openmc.lib.Tally
        Tally to which the geometry derivative is applied.
    tally_id : int
        ID of the tally to which the geometry derivative is applied.

    """
    __instances = WeakValueDictionary()

    def __new__(cls, uid=None, new=False, index=None):
        mapping = geometry_derivatives
        if index is None:
            if new:
                raise NotImplementedError(
                    'Geometry derivative allocation is not implemented.')
            index = mapping[uid]._index

        if index not in cls.__instances:
            instance = super().__new__(cls)
            instance._index = index
            cls.__instances[index] = instance

        return cls.__instances[index]

    @property
    def id(self):
        geom_deriv_id = c_int32()
        _dll.openmc_geometry_derivative_get_id(self._index, geom_deriv_id)
        return geom_deriv_id.value

    @property
    def tally_id(self):
        tally_id = c_int32()
        _dll.openmc_geometry_derivative_get_tally_id(self._index, tally_id)
        return tally_id.value

    @property
    def tally(self):
        return tallies[self.tally_id]

    @property
    def cell_id(self):
        cell_id = c_int32()
        _dll.openmc_geometry_derivative_get_cell_id(self._index, cell_id)
        return cell_id.value

    @property
    def cell(self):
        return cells[self.cell_id]

    @property
    def surface_ids(self):
        surface_ids = POINTER(c_int32)()
        n = c_size_t()
        _dll.openmc_geometry_derivative_get_surface_ids(
            self._index, surface_ids, n)
        return [surface_ids[i] for i in range(n.value)]

    @property
    def geom_parameters(self):
        data = POINTER(c_double)()
        shape = (c_size_t*2)()
        _dll.openmc_geometry_derivative_parameters(self._index, data, shape)
        return as_array(data, tuple(shape))

    @property
    def results(self):
        data = POINTER(c_double)()
        shape = (c_size_t*3)()
        _dll.openmc_geometry_derivative_results(self._index, data, shape)
        return as_array(data, tuple(shape))

    @property
    def mean(self):
        n = self.tally.num_realizations
        sum_ = self.results[:, :, 1]
        if n > 0:
            return sum_ / n
        else:
            return sum_.copy()

    def reset(self):
        """Reset geometry derivative results."""
        _dll.openmc_geometry_derivative_reset(self._index)


class _GeometryDerivativeMapping(Mapping):
    def __getitem__(self, key):
        index = c_int32()
        try:
            _dll.openmc_get_geometry_derivative_index(key, index)
        except (AllocationError, InvalidIDError) as e:
            raise KeyError(str(e))
        return GeometryDerivative(index=index.value)

    def __iter__(self):
        for i in range(len(self)):
            yield GeometryDerivative(index=i).id

    def __len__(self):
        return _dll.openmc_geometry_derivatives_size()

    def __repr__(self):
        return repr(dict(self))

    def __delitem__(self, key):
        raise NotImplementedError(
            "GeometryDerivative object remove not implemented")


geometry_derivatives = _GeometryDerivativeMapping()
