from ctypes import (c_bool, c_char, c_int, c_size_t, c_int32, c_uint8,
                    c_double, Structure, POINTER)

from . import _dll
from .error import _error_handler

import numpy as np


class _Position(Structure):
    """Definition of an xyz location in space with underlying c-types

    C-type Attributes
    -----------------
    x : c_double
        Position's x value (default: 0.0)
    y : c_double
        Position's y value (default: 0.0)
    z : c_double
        Position's z value (default: 0.0)
    """
    _fields_ = [('x', c_double),
                ('y', c_double),
                ('z', c_double)]

    def __getitem__(self, idx):
        if idx == 0:
            return self.x
        elif idx == 1:
            return self.y
        elif idx == 2:
            return self.z
        else:
            raise IndexError("{} index is invalid for _Position".format(idx))

    def __setitem__(self, idx, val):
        if idx == 0:
            self.x = val
        elif idx == 1:
            self.y = val
        elif idx == 2:
            self.z = val
        else:
            raise IndexError("{} index is invalid for _Position".format(idx))

    def __repr__(self):
        return "({}, {}, {})".format(self.x, self.y, self.z)


class _RGBColor(Structure):
    """Structure defining an RGB color

    C-type Attributes
    -----------------
    red : c_unit8
        Color's red value
    green : c_uint8
        Color's green value
    blue : c_uint8
        Color's blue value
    """
    _fields_ = [('red', c_uint8),
                ('green', c_uint8),
                ('blue', c_uint8)]

    def __repr__(self):
        return "RGB Color: {}, {}, {}".format(self.red, self.green. self.blue)

    @property
    def red(self):
        return self.red

    @red.setter
    def red(self, val):
        self.red = val

    @property
    def green(self):
        return self.green

    @green.setter
    def green(self, val):
        self.green = val

    @property
    def blue(self):
        return self.blue

    @property
    def blue(self, val):
        self.blue = val


class _PlotBase(Structure):
    """A structure defining a 2-D geometry slice with underlying c-types

    C-Type Attributes
    -----------------
    origin_ : openmc.lib.plot._Position
        A position defining the origin of the plot.
    width_ : openmc.lib.plot._Position
        The width of the plot along the x, y, and z axes, respectively
    basis_ : c_int
        The axes basis of the plot view.
    pixels_ : c_size_t[3]
        The resolution of the plot in the horizontal and vertical dimensions
    level_ : c_int
        The universe level for the plot view

    Attributes
    ----------
    origin : tuple or list of float
        Origin (center) of the plot
    width : float
        The horizontal dimension of the plot in geometry units (cm)
    height : float
        The vertical dimension of the plot in geometry units (cm)
    basis : string
        One of {'xy', 'xz', 'yz'} indicating the horizontal and vertical
        axes of the plot.
    h_res : int
        The horizontal resolution of the plot in pixels
    v_res : int
        The vertical resolution of the plot in pixels
    level : int
        The universe level for the plot (default: -1 -> all universes shown)
    """
    _fields_ = [('origin_', _Position),
                ('width_', _Position),
                ('basis_', c_int),
                ('pixels_', 3*c_size_t),
                ('color_overlaps_', c_bool),
                ('level_', c_int)]

    def __init__(self):
        self.level_ = -1
        self.color_overlaps_ = False

    @property
    def origin(self):
        return self.origin_

    @property
    def width(self):
        return self.width_.x

    @property
    def height(self):
        return self.width_.y

    @property
    def basis(self):
        if self.basis_ == 1:
            return 'xy'
        elif self.basis_ == 2:
            return 'xz'
        elif self.basis_ == 3:
            return 'yz'

        raise ValueError("Plot basis {} is invalid".format(self.basis_))

    @basis.setter
    def basis(self, basis):
        if isinstance(basis, str):
            valid_bases = ('xy', 'xz', 'yz')
            basis = basis.lower()
            if basis not in valid_bases:
                raise ValueError("{} is not a valid plot basis.".format(basis))

            if basis == 'xy':
                self.basis_ = 1
            elif basis == 'xz':
                self.basis_ = 2
            elif basis == 'yz':
                self.basis_ = 3
            return

        if isinstance(basis, int):
            valid_bases = (1, 2, 3)
            if basis not in valid_bases:
                raise ValueError("{} is not a valid plot basis.".format(basis))
            self.basis_ = basis
            return

        raise ValueError("{} of type {} is an"
                         " invalid plot basis".format(basis, type(basis)))

    @property
    def h_res(self):
        return self.pixels_[0]

    @property
    def v_res(self):
        return self.pixels_[1]

    @property
    def level(self):
        return int(self.level_)

    @property
    def color_overlaps(self):
        return self.color_overlaps_

    @color_overlaps.setter
    def color_overlaps(self, color_overlaps):
        self.color_overlaps_ = color_overlaps

    @origin.setter
    def origin(self, origin):
        self.origin_.x = origin[0]
        self.origin_.y = origin[1]
        self.origin_.z = origin[2]

    @width.setter
    def width(self, width):
        self.width_.x = width

    @height.setter
    def height(self, height):
        self.width_.y = height

    @h_res.setter
    def h_res(self, h_res):
        self.pixels_[0] = h_res

    @v_res.setter
    def v_res(self, v_res):
        self.pixels_[1] = v_res

    @level.setter
    def level(self, level):
        self.level_ = level

    @property
    def color_overlaps(self):
        return self.color_overlaps_

    @color_overlaps.setter
    def color_overlaps(self, val):
        self.color_overlaps_ = val

    def __repr__(self):
        out_str = ["-----",
                   "Plot:",
                   "-----",
                   "Origin: {}".format(self.origin),
                   "Width: {}".format(self.width),
                   "Height: {}".format(self.height),
                   "Basis: {}".format(self.basis),
                   "HRes: {}".format(self.h_res),
                   "VRes: {}".format(self.v_res),
                   "Color Overlaps: {}".format(self.color_overlaps),
                   "Level: {}".format(self.level)]
        return '\n'.join(out_str)


class Plot(_PlotBase):
    """A structure defining a plot image to be generated by OpenMC

    C-Type Attributes
    -----------------
    id_ : c_int
        The plot ID
    type_ : c_int
        The type of plot. Either 'slice' (1) or 'voxel' (2)
    color_by_ : c_int
        How to color the plot. Either 0 (cell) or 1 (material)
    meshlines_width_ : c_int
        Thickness of the meshlines drawn on the plot.
    meshlines_color_ : openmc.lib.plot._RGBColor
        Color of the meshlines to draw on the plot
    not_found_ : openmc.lib.plot._RGBColor
        Color of regions outside the model
    overlap_color_ : openmc.lib.plot._RGBColor
        Color of overlapping regions in the geometry
    colors_ : Pointer to an array of openmc.lib.plot._RGBColor
        The colors used for each cell or material in the plot
    path_plot_ : Pointer to a character array
        Path to the file to be generated by this plot

    Attributes
    ----------
    id : int
        The plot ID
    type : str
        The type of plot. Either 'slice' or 'voxel'.
    color_by : str
        How to color the plot. Either 'cell' or 'material'.
    meshlines_width : int
        Thickness of the meshlines drawn on the plot.
    meshlines_color : length 3 iterable of int
        Color of the meshlines to draw on the plot
    not_found : length 3 iterable of ints
        Color of regions outside the model
    overlap_color : length 3 iterable of ints
        Color of overlapping regions in the geometry
    colors : Dict-like
        Mapping from cell or material to colors in the plot
    path_plot : str
        Path to the file to be generated by this plot
    """
    _fields_ = [('id_', c_int),
                ('type_', c_int),
                ('color_by_', c_int),
                ('meshlines_width_', c_int),
                ('meshlines_color_', _RGBColor),
                ('not_found_', _RGBColor),
                ('overlap_color_', _RGBColor),
                ('colors_', POINTER(_RGBColor)),
                ('path_plot_', POINTER(c_char))]

    def __init__(self, id=None):
        super().__init__()

        if id:
            self.id = id

    @property
    def id(self):
        return self.id_

    @property
    def type(self):
        return self.type_

    @property
    def color_by(self):
        return self.color_by_

    @property
    def meshlines_witdh(self):
        return self.meshlines_width_

    @property
    def index_meshlines_width_(self):
        return self.meshliens_width_

    @property
    def meshlines_color(self):
        return self.meshlines_color_

    @property
    def not_found_color(self):
        return self.not_found_

    @property
    def overlap_color(self):
        return self.overlap_color_

    @property
    def colors(self):
        return self.colors_

    @property
    def plot_path(self):
        return self.plot_path_

_dll.openmc_id_map.argtypes = [POINTER(_PlotBase), POINTER(c_int32)]
_dll.openmc_id_map.restype = c_int
_dll.openmc_id_map.errcheck = _error_handler


def id_map(plotbase):
    """
    Generate a 2-D map of cell and material IDs. Used for in-memory image
    generation.

    Parameters
    ----------
    plot : openmc.lib.plot._PlotBase
        Object describing the slice of the model to be generated

    Returns
    -------
    id_map : numpy.ndarray
        A NumPy array with shape (vertical pixels, horizontal pixels, 3) of
        OpenMC property ids with dtype int32

    """
    img_data = np.zeros((plotbase.v_res, plotbase.h_res, 3),
                        dtype=np.dtype('int32'))
    _dll.openmc_id_map(plotbase, img_data.ctypes.data_as(POINTER(c_int32)))
    return img_data


_dll.openmc_property_map.argtypes = [POINTER(_PlotBase), POINTER(c_double)]
_dll.openmc_property_map.restype = c_int
_dll.openmc_property_map.errcheck = _error_handler


def property_map(plotbase):
    """
    Generate a 2-D map of cell temperatures and material densities. Used for
    in-memory image generation.

    Parameters
    ----------
    plot : openmc.lib.plot._PlotBase
        Object describing the slice of the model to be generated

    Returns
    -------
    property_map : numpy.ndarray
        A NumPy array with shape (vertical pixels, horizontal pixels, 2) of
        OpenMC property ids with dtype float

    """
    prop_data = np.zeros((plotbase.v_res, plotbase.h_res, 2))
    _dll.openmc_property_map(plotbase, prop_data.ctypes.data_as(POINTER(c_double)))
    return prop_data

_dll.openmc_create_image.argtypes = [POINTER(Plot)]
_dll.openmc_create_image.restype = c_int
_dll.openmc_create_image.errcheck = _error_handler

def create_image(plot, filename):
    """
    Write the plot to file

    Parameters
    ----------
    plot : openmc.lib.plot.Plot
        Plot to be generated and written to file
    filename : str or Path
        Name of the image to be generated
    """

    old_path = plot.plot_path
    plot.plot_path = filename

    _dll.openmc_create_image(POINTER(plot))



