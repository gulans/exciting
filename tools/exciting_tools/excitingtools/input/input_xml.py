"""Generate an exciting XML input tree.
"""
from pathlib import Path
from typing import Union
from xml.etree import ElementTree

from excitingtools.input.base_class import ExcitingXMLInput, AbstractExcitingInput
from excitingtools.input.xml_utils import xml_tree_to_pretty_str, prettify_tag_attributes


class ExcitingTitleInput(AbstractExcitingInput):
    """ Holds only the title but for consistency reasons as class. """
    name = "title"

    def __init__(self, title: str):
        self.title = title

    def to_xml(self, **kwargs) -> ElementTree:
        """ Puts title to xml, only the text is title. """
        title_tree = ElementTree.Element(self.name)
        title_tree.text = self.title
        return title_tree


class ExcitingInputXML(ExcitingXMLInput):
    """
    Container for a complete input xml file.
    """
    name = "input"
    _default_filename = "input.xml"

    def set_title(self, title: str):
        """ Set a new title. """
        self.__dict__["title"].title = title

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
