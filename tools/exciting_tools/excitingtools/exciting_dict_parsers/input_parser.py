""" Parsers for input.xml.

TODO(Fabian): Issues 117 & 121:
As more sub-elements are implemented in the input files, also add parsers here
"""
import copy
from typing import Tuple
from xml.etree import ElementTree

from excitingtools.parser_utils.parser_decorators import xml_root
from excitingtools.parser_utils.parser_utils import find_element, convert_string_dict


@xml_root
def parse_input_xml(root):
    """ Parse an input.xml file into dictionary. """
    assert root.tag == "input"
    return parse_element_xml(root)


def get_root_from_tag(root: ElementTree.Element, tag: str = None) -> Tuple[ElementTree.Element, str]:
    """ Get the root from a tag.
    :param tag: tag of interest
    :param root: xml root containing the tag (or having the specified tag as tag)
    :returns: the tag and the found root, if tag was None returns the tag of the given root
    """
    if tag is None:
        return root, root.tag

    root = find_element(root, tag)
    if root is None:
        raise ValueError(f"Your specified input has no tag {tag}.")

    return root, root.tag


@xml_root
def parse_element_xml(root, tag: str = None) -> dict:
    """ Parse a xml element into dictionary. Can be input.xml root or a subelement of it.
    Put the attributes simply in dict and add recursively the subtrees and nested dicts.
    Note: Parses all attributes into strings, would be nice to have as their actual data type but this
     requires some more overhead. Maybe look at 'ast.literal_eval()'.
    :param tag: the tag to parse
    :param root: the xml root containing the tag
    :returns: the parsed dictionary """
    root, tag = get_root_from_tag(root, tag)

    if tag in special_tags_to_parse_map.keys():
        return special_tags_to_parse_map[tag](root)

    element_dict = convert_string_dict(copy.deepcopy(root.attrib))

    subelements = list(root)
    for subelement in subelements:
        element_dict[subelement.tag] = parse_element_xml(subelement)

    return element_dict


@xml_root
def parse_structure(root) -> dict:
    """ Parse exciting input.xml structure element into python dictionary.
    :param root: Input for the parser.
    :returns: Dictionary containing the structure input element attributes and subelements. Looks like:
        {'atoms': List of atoms with atom positions in fractional coordinates,
         'lattice': List of 3 lattice vectors, 'species_path': species_path as string,
         'crystal_properties': dictionary with the crystal_properties,
         'species_properties': dictionary with the species_properties,
         all additional keys are structure attributes}
    """
    structure = find_element(root, 'structure')
    structure_properties = structure.attrib
    species_path = structure_properties.pop('speciespath')
    crystal = structure.find('crystal')
    crystal_properties = crystal.attrib
    lattice = []
    for base_vect in crystal.findall('basevect'):
        lattice.append([float(x) for x in base_vect.text.split()])

    atoms = []
    species_properties = {}
    for species in structure.findall('species'):
        species_attributes = species.attrib
        species_file = species_attributes.pop('speciesfile')
        species_symbol = species_file[:-4]
        species_properties[species_symbol] = species_attributes
        for atom in species:
            atom_attributes = atom.attrib
            coord = [float(x) for x in atom_attributes.pop('coord').split()]
            atom_dict = {'species': species_symbol, 'position': coord}
            atom_dict.update(atom_attributes)
            atoms.append(atom_dict)

    return {
        'atoms': atoms,
        'lattice': lattice,
        'species_path': species_path,
        'crystal_properties': crystal_properties,
        'species_properties': species_properties,
        **structure_properties
    }


# special tag to parse function map or lambda if one-liner
# necessary for tags which doesn't contain simply xml attributes and subtrees
special_tags_to_parse_map = {"title": lambda root: root.text,
                             "structure": parse_structure,
                             "qpointset": lambda root: [[float(x) for x in qpoint.text.split()] for qpoint in root],
                             "plan": lambda root: [doonly.attrib['task'] for doonly in root]}
