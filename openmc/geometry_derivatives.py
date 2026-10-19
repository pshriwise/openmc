
from __future__ import annotations
from collections.abc import Iterable, Mapping
import lxml.etree as ET

import openmc
import openmc.checkvalue as cv

from .mixin import IDManagerMixin


class GeometricDerivative(IDManagerMixin):
    """Class used to define the geometric derivative for a tally

    Parameters
    ----------
    id : int
        Unique ID for the geometric derivative
    tally : openmc.Tally
        Tally to which the geometric derivative is applied
    surface_parameters : dict
        Mapping of surface objects/IDs to iterables of surface parameters to differentiate with respect to.
    """
    next_id = 1
    used_ids = set()

    def __init__(self, tally : openmc.Tally, surface_parameters: dict, id: int = None):
        self.id = id
        self.tally = tally
        self.surface_parameters = surface_parameters

    @property
    def surface_parameters(self):
        return self._surface_parameters

    @surface_parameters.setter
    def surface_parameters(self, value):
        cv.check_type('Surface parameters', value, Mapping)
        value = dict(value)

        for key, val in value.items():
            cv.check_type('Surface', key, (openmc.Surface, int))

            if val == ['all']:
                if isinstance(key, openmc.Surface):
                    val = value[key] = key.differentiable_params()
                else:
                    raise ValueError(f"Surface parameters for surface ID {key} cannot be set to ['all'] without a Surface object.")
            cv.check_type('Surface parameters', val, Iterable)
            cv.check_iterable_type('Surface parameters', val, str)

            for param in val:
                if isinstance(key, openmc.Surface):
                    if param not in key.differentiable_params():
                        raise ValueError(f"Invalid surface parameter: {param}")

        self._surface_parameters = value

    def n_bins(self):
        """Return the number of bins in the geometric derivative

        Returns
        -------
        int
            Number of bins in the geometric derivative
        """
        # the number of bins in this geometric deriative
        # is the number of tally bins multiplied by the number
        # of differentiable surface parameters

        n_surface_parameters = sum(len(params) for params in self.surface_parameters.values())
        return self.tally.num_filter_bins * n_surface_parameters

    def to_xml_element(self):
        """Create a 'derivative' element to be written to an XML file.

        Returns
        -------
        lxml.etree._Element
            XML element

        """
        element = ET.Element("geometry_derivative")
        element.set("id", str(self.id))
        element.set("tally", str(self.tally.id))

        # create a subelement for each set of surface parameters (one per surface)
        for surface, params in self.surface_parameters.items():
            surface_elem = ET.SubElement(element, "surface")
            surface_elem.set("id", str(surface.id if isinstance(surface, openmc.Surface) else surface))
            surface_elem.set("parameters", " ".join(params))

        return element

    @classmethod
    def from_xml_element(cls, elem, tallies=None, surfaces=None):
        """Generate a geometric derivative from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element
        tallies : dict or None
            A dictionary with tally IDs as keys and tally instances as values that
            have already been read from XML. Pre-existing tallies are used
            and new tallies are added to when creating geometric derivative objects.
        surfaces : dict or None
            A dictionary with surface IDs as keys and surface instances as values that
            have already been read from XML. Pre-existing surfaces are used
            and new surfaces are added to when creating geometric derivative objects.

        Returns
        -------
        openmc.GeometricDerivative
            Geometric derivative object

        """
        id = int(elem.get("id"))
        tally = tallies[int(elem.get("tally"))]
        out = cls(tally, {}, id=id)
        surface_parameters = {}
        for surface_elem in elem.findall("surface"):
            surface_id = int(surface_elem.get("id"))
            surface = surfaces[surface_id]
            params = surface_elem.get("parameters").split()
            surface_parameters[surface] = params

        out.surface_parameters = surface_parameters

        return out


class GeometricDerivatives(cv.CheckedList):
    """Collection of GeometricDerivatives used for an OpenMC simulation.

    A staging container for serialization of GeometricDerivatives. It can be
    thought of as a normal Python list where each member is a
    :class:`GeometricDerivative`. It behaves like a list as the following
    example demonstrates:

    >>> t1 = openmc.Tally()
    >>> s1 = openmc.Sphere()
    >>> gd1 = openmc.GeometricDerivative(t1, {s1: ['r']})
    >>> gd2 = openmc.GeometricDerivative(t1, {s1: ['x0', 'y0', 'z0']})
    >>> gds = openmc.GeometricDerivatives([gd1])
    >>> gds.append(gd2)

    Parameters
    ----------
    geometric_derivatives : Iterable of openmc.GeometricDerivative
        GeometricDerivatives to add to the collection

    """

    def __init__(self, geometric_derivatives=None):
        super().__init__(GeometricDerivative, 'geometric derivatives collection')
        if geometric_derivatives is not None:
            self += geometric_derivatives

    def to_xml_element(self):
        """Creates a 'derivatives' element to be written to an XML file.
        """
        element = ET.Element("geometry_derivatives")
        for gd in self:
            element.append(gd.to_xml_element())
        return element

    @classmethod
    def from_xml_element(cls, elem, tallies=None, surfaces=None):
        """Generate geometric derivatives from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element
        tallies : dict or None
            A dictionary with tally IDs as keys and tally instances as values that
            have already been read from XML. Pre-existing tallies are used
            and new tallies are added to when creating geometric derivative objects.
        surfaces : dict or None
            A dictionary with surface IDs as keys and surface instances as values that
            have already been read from XML. Pre-existing surfaces are used
            and new surfaces are added to when creating geometric derivative objects.

        Returns
        -------
        openmc.GeometricDerivatives
            Collection of geometric derivative objects

        """
        gds = cls()
        for gd_elem in elem.findall("geometry_derivative"):
            gds.append(GeometricDerivative.from_xml_element(gd_elem, tallies, surfaces))
        return gds
