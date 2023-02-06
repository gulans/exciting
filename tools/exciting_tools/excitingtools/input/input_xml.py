"""Generate an exciting XML input tree.
"""
import copy
from pathlib import Path
from typing import Union
from xml.etree import ElementTree

from excitingtools.input.base_class import ExcitingInput
from excitingtools.input.structure import ExcitingStructure
from excitingtools.input.ground_state import ExcitingGroundStateInput
from excitingtools.input.xs import ExcitingXSInput
from excitingtools.input.xml_utils import xml_tree_to_pretty_str, prettify_tag_attributes
from excitingtools.exciting_obj_parsers.input_xml import parse_input_xml
from excitingtools.utils.jobflow_utils import special_serialization_attrs


# mapping from xml tag to python constructor for deserialisation
element_constructor = {'groundstate': ExcitingGroundStateInput,
                       'structure': ExcitingStructure,
                       'xs': ExcitingXSInput}


class ExcitingInputXML(ExcitingInput):
    _default_filename = "input.xml"

    def __init__(self,
                 structure: ExcitingStructure,
                 groundstate: ExcitingGroundStateInput,
                 **kwargs):
        """
        Container for a complete input xml file.

        :param structure: Structure containing lattice vectors and atomic positions.
        :param groundstate: exciting ground state input object.
        Optional arguments in kwargs:
        :param title: optional title for input file.
        :param xs: exciting xs input object.
        """
        self.structure = structure
        self.groundstate = groundstate
        self.title = kwargs.pop('title', "Exciting Calculation")
        for key, value in kwargs.items():
            if not isinstance(value, ExcitingInput):
                raise ValueError(f"Only ExcitingInput is allowed, key {key} has a value with type {type(value)}.")
            setattr(self, key, value)

    def as_dict(self) -> dict:
        """Returns a dictionary representing the current object for later recreation.
        The special serialise attributes are required for recognition by monty and jobflow.
        """
        serialise_attrs = special_serialization_attrs(self)
        inputs = vars(self)
        title = inputs.pop("title")
        input_dict = {key: value.as_dict() for key, value in inputs.items()}
        return {**serialise_attrs, "title": title, **input_dict}

    @classmethod
    def from_dict(cls, d: dict):
        my_dict = copy.deepcopy(d)
        # Remove key value pairs needed for workflow programs
        serialise_keys = special_serialization_attrs(cls)
        for key in serialise_keys:
            my_dict.pop(key, None)

        input_attributes = {"title": my_dict.pop("title")}
        # all other keys are related to ExcitingInput objects
        for element_name, field in my_dict.items():
            input_attributes[element_name] = element_constructor[element_name].from_dict(field)

        return cls(**input_attributes)

    @classmethod
    def from_xml(cls, input_file):
        return cls(**parse_input_xml(input_file))

    def initialise_input_xml(self) -> ElementTree.Element:
        """Initialise input.xml element tree for exciting.

        Information on the schema reference ignored, but could be reintroduced for validation purposes.
        root.set(
           '{http://www.w3.org/2001/XMLSchema-instance}noNamespaceSchemaLocation',
           'http://xml.exciting-code.org/excitinginput.xsd')

        :return: Element tree root.
        """
        root = ElementTree.Element('input')
        ElementTree.SubElement(root, 'title').text = self.title
        return root

    def to_xml(self) -> ElementTree.Element:
        """Compose XML ElementTrees from exciting input classes to create an input XML tree.

        :return: Input XML tree, with sub-elements inserted.
        """
        root = self.initialise_input_xml()

        elements = copy.deepcopy(vars(self))
        elements.pop("title")
        for element in elements.values():
            root.append(element.to_xml())

        return root

    def to_xml_str(self) -> str:
        """Compose XML ElementTrees from exciting input classes to create an input xml string.

        :return: Input XML tree as a string, with pretty formatting.
        """
        xml_tree = self.to_xml()
        tags_to_prettify = ["\t<structure", "\t\t<crystal", "\t\t<species", "\t\t\t<atom", "\t<groundstate", "\t<xs",
                            "\t\t<screening", "\t\t<BSE", "\t\t<energywindow"]
        input_xml_str = prettify_tag_attributes(xml_tree_to_pretty_str(xml_tree), tags_to_prettify)
        return input_xml_str

    def write(self, filename: Union[str, Path] = _default_filename):
        """Writes the xml string to file.

        :param filename: name of the file.
        """
        with open(filename, "w") as fid:
            fid.write(self.to_xml_str())
