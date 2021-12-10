from collections.abc import Iterable
from pathlib import Path

from xml.etree import ElementTree as ET

from openmc.mesh import MeshBase
import openmc.checkvalue as cv

from .._xml import clean_indentation, get_text


class VarianceReduction():
    """ A class to handle various Variance Reduction operations
    """
    def __init__(self):
        self._weight_window_domains = []

    def __repr__(self):
        string = type(self).__name__ +  '\n\n'
        string += '{0: <16}\n'.format('Weight Window Domains:')
        for wwd in self._weight_window_domains:
            string += '\t{}\n'.format(wwd)
        return string

    @property
    def weight_window_domains(self):
        return self._weight_window_domains

    @weight_window_domains.setter
    def weight_window_domains(self, domains):
        cv.check_type('Weight window domains', domains, Iterable)
        cv.check_iterable_type('Weight window domains',
                               domains,
                               WeightWindowDomain)
        self._weight_window_domains = domains

    def export_to_xml(self, path='variance_reduction.xml'):
        """Exports the variance reduction parameters to an XML file

        Parameters
        ----------
        path : str
            Path to file to write. Defaults to 'variance_reduction.xml'.
        """
        # Check if path is a directory
        p = Path(path)
        if p.is_dir():
            p /= 'variance_reduction.xml'

        with open(str(p), 'w', encoding='utf-8',
                  errors='xmlcharrefreplace') as fh:

            # Write the header and the opening tag for the root element.
            fh.write("<?xml version='1.0' encoding='utf-8'?>\n")

            # Create XML representation
            root_element = ET.Element("variance_reduction")

            if self.weight_window_domains:
                ww_element = ET.SubElement(root_element, 'weight_windows')

                for domain in self.weight_window_domains:
                    domain_element = domain.to_xml_element()
                    clean_indentation(domain_element, level=1)
                    ww_element.append(domain_element)

                    settings_element = domain.settings.to_xml_element()
                    clean_indentation(settings_element, level=1)
                    ww_element.append(settings_element)

                    mesh_element = domain.mesh.to_xml_element()
                    clean_indentation(mesh_element, level=1)
                    root_element.append(mesh_element)

                clean_indentation(ww_element)

            clean_indentation(root_element)

            # Write the XML Tree
            tree = ET.ElementTree(root_element)
            tree.write(str(p), xml_declaration=True, encoding='utf-8')

    @classmethod
    def from_xml(cls, path='variance_reduction.xml'):
        """Generate variance reduction parameters from XML file

        Parameters
        ----------
        path : str, optional
            Path to the variance reduction XML file

        Returns
        -------
        openmc.VarianceReduction
            VarianceReduction object

        """
        tree = ET.parse(path)
        root = tree.getroot()

        vr = cls()

        # read any meshes in the file first
        meshes = {}
        for mesh_elem in root.findall('mesh'):
            try:
                mesh = MeshBase.from_xml(mesh_elem)
            except AttributeError as e:
                msg = ('Can only read Regular or Unstructured meshes from XML')
                raise e(msg)

            meshes[mesh.id] = mesh

        # get the weight window node
        ww_elem = root.find('weight_windows')

        if ww_elem:
            ww_settings = {}
            for settings_elem in ww_elem.findall('settings'):
                settings = \
                        WeightWindowSettings.from_xml_element(settings_elem)
                ww_settings[settings.id] = settings

            # read weight window domains
            for domain_elem in ww_elem.findall('domain'):
                id = int(get_text(domain_elem, 'id'))
                mesh_id = int(get_text(domain_elem, 'mesh'))
                settings_id = int(get_text(domain_elem, 'settings'))

                domain = WeightWindowDomain(id=id,
                mesh=meshes[mesh_id],
                settings=ww_settings[settings_id])

                vr.weight_window_domains.append(domain)

        return vr

    @classmethod
    def from_hdf5(cls, group, meshes):
        """Generate a variance reduction class from HDF5 group

        Paramters
        ----------
        group : h5py.Group
            HDF5 group containing variance reduction information
        meshes: Mapping
            Set of meshes available in the model

        Returns
        -------
        openmc.VarianceReduction
            A variance reduction object
        """
        vr = cls()

        ww_settings_group = group['weight_window_settings']
        ww_settings = {}
        for settings_group in ww_settings_group.values():
            settings = WeightWindowSettings.from_hdf5(settings_group)
            ww_settings[settings.id] = settings

        ww_domains_group = group['weight_window_domains']
        for domains_group in ww_domains_group.values():
            id = int(domains_group.name.split('/')[-1].lstrip('domain'))
            mesh_id = domains_group['mesh'][()]
            settings_id = domains_group['settings'][()]

            domain = WeightWindowDomain(id,
                                        meshes[mesh_id],
                                        ww_settings[settings_id])
            vr.weight_window_domains.append(domain)

        return vr


