from collections.abc import Iterable
from numbers import Real, Integral

from xml.etree import ElementTree as ET
import numpy as np
from numpy.core.function_base import linspace

from openmc.filter import _PARTICLES
from openmc.mesh import MeshBase, RectilinearMesh
import openmc.checkvalue as cv

from .._xml import clean_indentation, get_text
from ..mixin import IDManagerMixin


class WeightWindowSettings(IDManagerMixin):
    """ A class to handle the creation of a set of specific weight window
    paramaters - a variance reduction class may have several of these

    Parameters
    ----------
    id : int
       Unique identifier for the weight window settings. If not
       specified an identifier will automatically be assigned.
    particle_types : str (default is 'netron')
        Particle type the weight windows will apply to
    energy_bins : Iterable of Real
        A list of values for which each successive pair constitutes a range of
        energies in [eV] for a single bin
    lower_ww_bounds : Iterable of Real
        A list of values for which each value is the lower bound of a weight
        window
    upper_bound_ratio : float
        Ratio of the lower to upper weight window bounds
    upper_ww_bounds : Iterable of Real
        A list of values for which each value is the upper bound of a weight
        window
    survival_ratio : float (default is 3)
        Ratio of the lower weight window bound to the survival weight for
        rouletting
    max_split : int (default is 10)
        Maximum allowable number of particles when splitting
    weight_cutoff : float (default is 1E-38)
        Threshold below which particles will be terminated

    Attributes
    ----------
    id : int
       Unique identifier for the weight window settings.
    particle_type : str
        Particle type the weight windows apply to
    energy_bins : Iterable of Real
        A list of values for which each successive pair constitutes a range of
        energies in [eV] for a single bin
    lower_ww_bounds : Iterable of Real
        A list of values for which each value is the lower bound of a weight
        window
    upper_ww_bounds : Iterable of Real
        A list of values for which each value is the upper bound of a weight
        window
    survival_ratio : float
        Ratio of the lower weight window bound to the survival weight for
        rouletting
    max_split : int
        Maximum allowable number of particles when splitting
    weight_cutoff : float
        Threshold below which particles will be terminated
    """
    next_id = 1
    used_ids = set()

    def __init__(self,
                id=None,
                particle_type=None,
                energy_bins=None,
                lower_ww_bounds=None,
                upper_bound_ratio=None,
                upper_ww_bounds=None,
                survival_ratio=3,
                max_split=10,
                weight_cutoff=1.e-38):

        self.id = id

        if particle_type:
            self.particle_type = particle_type
        else:
            self.particle_type = 'neutron'

        self.energy_bins = energy_bins
        self.lower_ww_bounds = lower_ww_bounds

        cv.check_length('Lower window bounds',
        self.lower_ww_bounds,
        len(self.energy_bins))

        if upper_ww_bounds is not None and upper_bound_ratio:
            msg = ("Exactly one of uppwer_ww_bounds and "
                   "upper_bound_ratio must be present.")
            raise ValueError(msg)

        if upper_ww_bounds is None and upper_bound_ratio is None:
            msg = ("Exactly one of uppwer_ww_bounds and "
                   "upper_bound_ratio must be present.")
            raise ValueError(msg)

        if upper_bound_ratio:
            self.upper_ww_bounds = \
                [lb * upper_bound_ratio for lb in self.lower_ww_bounds]

        if upper_ww_bounds is not None:
            self.upper_ww_bounds = upper_ww_bounds

        if len(self.lower_ww_bounds) != len(self.upper_ww_bounds):
            msg = ('Size of the lower and upper weight '
                    'window bounds do not match')
            raise ValueError(msg)

        self.survival_ratio = survival_ratio
        self.max_split = max_split
        self.weight_cutoff = weight_cutoff

    def __repr__(self):
        string = type(self).__name__ + '\n'
        string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
        string += '{0: <16}{1}{2}\n'.format('\tParticle Type', '=\t', self._particle_type)
        string += '{0: <16}{1}{2}\n'.format('\tEnergy Bins', '=\t', self._energy_bins)
        string += '{0: <16}{1}{2}\n'.format('\tLower WW Bounds', '=\t', self._lower_ww_bounds)
        string += '{0: <16}{1}{2}\n'.format('\tUpper WW Bounds', '=\t', self._upper_ww_bounds)
        string += '{0: <16}{1}{2}\n'.format('\tSurvival Ratio', '=\t', self._survival_ratio)
        string += '{0: <16}{1}{2}\n'.format('\tMax Split', '=\t', self._max_split)
        string += '{0: <16}{1}{2}\n'.format('\tWeight Cutoff', '=\t', self._weight_cutoff)
        return string

    @property
    def particle_type(self):
        return self._particle_type

    @particle_type.setter
    def particle_type(self, pt):
        cv.check_value('Particle type', pt, list(_PARTICLES))
        self._particle_type = pt

    @property
    def energy_bins(self):
        return self._energy_bins

    @energy_bins.setter
    def energy_bins(self, bins):
        cv.check_type('Energy bins', bins, Iterable)
        cv.check_iterable_type('Energy value', bins, Real)
        self._energy_bins = np.array(bins)

    @property
    def lower_ww_bounds(self):
        return self._lower_ww_bounds

    @lower_ww_bounds.setter
    def lower_ww_bounds(self, bounds):
        cv.check_type('Lower WW bounds', bounds, Iterable)
        cv.check_iterable_type('Weight window bound', bounds, Real)
        self._lower_ww_bounds = np.array(bounds)

    @property
    def upper_ww_bounds(self):
        return self._upper_ww_bounds

    @upper_ww_bounds.setter
    def upper_ww_bounds(self, bounds):
        cv.check_type('Upper WW bounds', bounds, Iterable)
        cv.check_iterable_type('Weight window bound', bounds, Real)
        self._upper_ww_bounds = np.array(bounds)

    @property
    def survival_ratio(self):
        return self._survival_ratio

    @survival_ratio.setter
    def survival_ratio(self, val):
        cv.check_type('Survival ratio', val, Real)
        if cv.check_greater_than('Survival ratio', val, 1.0, True):
            raise ValueError('Survival ratio cannot be less than one.')
        self._survival_ratio = val

    @property
    def max_split(self):
        return self._max_split

    @max_split.setter
    def max_split(self, val):
        cv.check_type('Max split', val, Integral)
        self._max_split = val

    @property
    def weight_cutoff(self):
        return self._weight_cutoff

    @weight_cutoff.setter
    def weight_cutoff(self, cutoff):
        cv.check_type('Weight cutoff', cutoff, Real)
        if cv.check_greater_than('Weight cutoff', cutoff, 0.0, True):
            raise ValueError('Weight cutoff must be greater than zero.')
        self._weight_cutoff = cutoff

    def scale_bounds(self, factor):
        """Scale the weight window bounds by some factor
        """
        # scale lower bounds
        for i, val in enumerate(self.lower_ww_bounds):
            self.lower_ww_bounds[i] = val * factor
        # scale upper bounds
        for i, val in enumerate(self.upper_ww_bounds):
            self.upper_ww_bounds[i] = val * factor

    def to_xml_element(self):
        """Return an XML representation of the weight windows

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing the weight window information
        """
        element = ET.Element('settings')

        element.set('id', str(self._id))

        subelement = ET.SubElement(element, 'particle_type')
        subelement.text = self.particle_type
        clean_indentation(subelement, level=2)

        subelement = ET.SubElement(element, 'energy_bins')
        subelement.text = ' '.join(str(e) for e in self.energy_bins)
        clean_indentation(subelement, level=2)

        subelement = ET.SubElement(element, 'lower_ww_bounds')
        subelement.text = ' '.join(str(b) for b in self.lower_ww_bounds)
        clean_indentation(subelement, level=2)

        subelement = ET.SubElement(element, 'upper_ww_bounds')
        subelement.text = ' '.join(str(b) for b in self.upper_ww_bounds)
        clean_indentation(subelement, level=2)

        subelement = ET.SubElement(element, 'survival_ratio')
        subelement.text = str(self.survival_ratio)
        clean_indentation(subelement, level=2)

        subelement = ET.SubElement(element, 'max_split')
        subelement.text = str(self.max_split)
        clean_indentation(subelement, level=2)

        subelement = ET.SubElement(element, 'weight_cutoff')
        subelement.text = str(self.weight_cutoff)
        clean_indentation(subelement, level=2)

        return element

    @classmethod
    def from_xml_element(cls, elem):
        """Generate weight window settings from an XML element

        Parameters
        ----------
        elem : xml.etree.ElementTree.Element
            XML element

        Returns
        -------
        openmc.WeightWindowSettings
            Weight window settings object
        """

        id = int(get_text(elem, 'id'))

        particle_type = get_text(elem, 'particle_type')
        ebins = \
            np.array([float(b) for b in get_text(elem, 'energy_bins').split()])

        lower_ww_bounds = \
            np.array([float(l) for l in get_text(elem, 'lower_ww_bounds').split()])

        upper_ww_bounds = \
            np.array([float(u) for u in get_text(elem, 'upper_ww_bounds').split()])

        survival_ratio = float(get_text(elem, 'survival_ratio'))

        max_split = int(get_text(elem, 'max_split'))

        weight_cutoff = float(get_text(elem, 'weight_cutoff'))

        ww_settings = cls(id=id,
        particle_type=particle_type,
        energy_bins=ebins,
        lower_ww_bounds=lower_ww_bounds,
        upper_ww_bounds=upper_ww_bounds,
        survival_ratio=survival_ratio,
        max_split=max_split,
        weight_cutoff=weight_cutoff)

        return ww_settings

    @classmethod
    def from_hdf5(cls, group):
        """Create weight window settings from HDF5 group

        Parameters
        ----------
        group : h5py.Group
            Group in HDF5 file

        Returns
        -------
        openmc.WeightwindowSettings
            A weight window settings object
        """

        id = int(group.name.split('/')[-1].lstrip('parameters'))

        ptype = group['particle_type'][()].decode()
        ebins = group['energy_bins'][()]
        lower_ww_bounds = group['lower_ww_bounds'][()]
        upper_ww_bounds = group['upper_ww_bounds'][()]
        survival_ratio = group['survival_ratio'][()]
        max_split = group['max_split'][()]
        weight_cutoff = group['weight_cutoff'][()]

        wws = cls(id=id,
        particle_type=ptype,
        energy_bins=ebins,
        lower_ww_bounds=lower_ww_bounds,
        upper_ww_bounds=upper_ww_bounds,
        survival_ratio=survival_ratio,
        max_split=max_split,
        weight_cutoff=weight_cutoff)

        return wws


class WeightWindowDomain(IDManagerMixin):
    """A class specifying a weight window domain as the combination
       of a mesh and weight window settings

    Parameter
    ---------
    mesh : openmc.MeshBase
        Mesh for the weight windows
    settings : openmc.WeightWindowSettings
        Settings for the weight window domains
    id : int
       Unique identifier for the weight window domain. If not
       specified an identifier will automatically be assigned.

    Attributes
    ----------
    id : int
        Unique identifier for the weight window domain.
    mesh : openmc.MeshBase
        Mesh for the weight windows
    settings : openmc.WeightWindowSettings
        Settings for the weight window domains

    """
    next_id = 1
    used_ids = set()

    def __init__(self, id=None, mesh=None, settings=None):

        self.id = id

        if not mesh:
            msg = ('Mesh object is required to '
            'create a weight window domain.')
            raise ValueError(msg)
        self.mesh = mesh

        if not settings:
            msg = ('Weight window settings are required to '
            'create a weight window domain.')
            raise ValueError(msg)
        self.settings = settings

    def __repr__(self):
        string = type(self).__name__ + '\n'
        string += '{0: <16}{1}{2}\n\n'.format('\tID', '=\t', self._id)
        string += '{0: <16}{1}{2}\n'.format('\tSettings:', '\n\t', self.settings)
        string += '{0: <16}{1}{2}\n'.format('\tMesh:', '\n\t', self.mesh)
        return string

    @property
    def mesh(self):
        return self._mesh

    @mesh.setter
    def mesh(self, mesh):
        cv.check_type('Weight window mesh', mesh, MeshBase)
        self._mesh = mesh

    @property
    def settings(self):
        return self._settings

    @settings.setter
    def settings(self, settings):
        cv.check_type('Weight window settings', settings, WeightWindowSettings)
        self._settings = settings

    def to_xml_element(self):
        element = ET.Element("domain")
        element.set('id', str(self.id))

        mesh_element = ET.SubElement(element, 'mesh')
        mesh_element.text = str(self.mesh.id)

        settings_element = ET.SubElement(element, 'settings')
        settings_element.text = str(self.settings.id)

        return element

    @staticmethod
    def wwinp(filename):
        """
        Returns the next value in the wwinp file.

        filename : str or pathlib.Path
            Location of the wwinp file
        """
        fh = open(filename, 'r')

        # read the first line of the file and
        # keep only the first four entries
        while(True):
            line = next(fh)
            if line and not line.startswith('c'):
                break

        values = line.strip().split()[:4]
        for value in values:
            yield value

        # the remainder of the file can be read as
        # sequential values
        while(True):
            line = next(fh)
            # skip empty or commented lines
            if not line or line.startswith('c'):
                continue
            values = line.strip().split()
            for value in values:
                yield value

    @classmethod
    def from_wwinp(cls, fh):
        """Reads a wwinp file into WeightWindowDomain's

        Parameters
        ----------
        fh : TextIOWrapper
            File handle of the MCNP weight window file, opened to the start of the block.

        Returns
        -------
        dict
            A mapping of the parameters read from the block
        """

        # create generator for getting the next parameter from the file
        wwinp = WeightWindowDomain.wwinp(fh)

        # first parameter, if, of wwinp file is unused
        next(wwinp)

        # check time parameter
        if int(float(next(wwinp))) > 1:
            raise ValueError('Time-dependent weight windows are not yet supported.')

        # number of particles
        ni = int(float(next(wwinp)))

        # read the mesh type
        nr = int(float(next(wwinp)))

        if nr != 10:
            # TODO: read the first entry by default and display a warning
            raise ValueError('Cylindrical meshes are not currently supported')

        # read the number of energy groups for each particle
        nes = [int(next(wwinp)) for _ in range(ni)]

        if len(nes) == 1:
            particles = ['neutron']
        elif len(nes) == 2:
            particles = ['neutron', 'photon']
        else:
            msg = ('More than two particle types are present. '
                   'Only neutron and photon weight windows will be read.')
            raise Warning(msg)

        # read number of fine mesh elements in each coarse element
        nfx = int(float(next(wwinp)))
        nfy = int(float(next(wwinp)))
        nfz = int(float(next(wwinp)))

        # read the mesh origin
        llc = tuple(float(next(wwinp)) for _ in range(3))

        # read the number of coarse mesh elements
        ncx = int(float(next(wwinp)))
        ncy = int(float(next(wwinp)))
        ncz = int(float(next(wwinp)))

        # skip the nwg value defining the geometry type, we already know this
        next(wwinp)


        # read the coordinates for each dimension into a rectilinear mesh
        mesh = RectilinearMesh()
        mesh.x_grid = cls._read_mesh_coords(wwinp, ncx)
        mesh.y_grid = cls._read_mesh_coords(wwinp, ncy)
        mesh.z_grid = cls._read_mesh_coords(wwinp, ncz)

        dims = ('x', 'y', 'z')

        # check consistency of mesh coordinates
        for dim, header_val, mesh_val in zip(dims, llc, mesh.llc):
            if header_val != mesh_val:
                msg = ('The {} corner of the mesh ({}) does not match '
                       'the value read in block 1 of the wwinp file ({})')
                raise ValueError(msg.format(dim, mesh_val, header_val))

        mesh_dims = mesh.dimension
        for dim, header_val, mesh_val in zip(dims, (nfx, nfy, nfz), mesh_dims):
            if header_val != mesh_val:
                msg = ('Total number of mesh elements read in the {} '
                       'direction ({}) is inconsistent with the '
                       'number read in block 1 of the wwinp file ({})')
                raise ValueError(msg.format(dim, mesh_val, header_val))

        # total number of fine mesh elements
        nft = nfx * nfy * nfz
        # read energy bins and weight window values for each particle
        ww_settings = []
        for particle, ne in zip(particles, nes):
            # read energy
            e_groups = np.asarray(list(float(next(wwinp)) for _ in range(ne)))

            # adjust energy from MeV to eV
            e_groups *= 1E6

            # create an array for weight window lower bounds
            ww_lb = np.zeros((ne, nft))
            for e in range(ne):
                ww_lb[e, :] = [float(next(wwinp)) for _ in range(nft)]

            settings = WeightWindowSettings(id=None,
                                            particle_type=particle,
                                            energy_bins=e_groups,
                                            # TODO: check order
                                            lower_ww_bounds=ww_lb.flatten(),
                                            upper_bound_ratio=5.0)
            ww_settings.append(settings)

        # tie weight window settings and mesh together in WeightWindowDomains
        return [cls(id=None, mesh=mesh, settings=s) for s in ww_settings]

    @staticmethod
    def _read_mesh_coords(wwinp, n_coarse_bins):
        """Read a set of mesh coordinates from a wwinp file

        Parameters
        ----------
        wwinp : IOTextWrapper
            wwinp file handle that will read the first value of the coordinates next.
        n_coarse_bins : int
            The number of coarse mesh bins to be read.
        """

        coords = [float(next(wwinp))]

        for _ in range(n_coarse_bins):
            # TODO: These are setup to read according to the MCNP5 format
            sx = int(float(next(wwinp))) # number of fine mesh elements in between
            px = float(next(wwinp))  # value of next coordinate
            qx = next(wwinp)  # this value is unused

            # append the fine mesh coordinates for this coarse element
            coords += list(np.linspace(coords[-1], px, sx + 1))[1:]

        return np.asarray(coords)

