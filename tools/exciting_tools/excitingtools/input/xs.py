"""Module for class of exciting xs (excited states).
https://exciting.wikidot.com/ref:xs
"""
from typing import Optional, List, Union
from xml.etree import ElementTree

import numpy as np

from excitingtools.input.base_class import AbstractExcitingInput, ExcitingXMLInput
from excitingtools.utils.dict_utils import check_valid_keys
from excitingtools.utils.utils import list_to_str
from excitingtools.utils.valid_attributes import valid_plan_entries


class ExcitingXSBSEInput(ExcitingXMLInput):
    """
    Class for exciting BSE Input
    """
    name = "BSE"


class ExcitingXSScreeningInput(ExcitingXMLInput):
    """
    Class for exciting Screening Input
    """
    name = "screening"


class ExcitingXSEnergywindowInput(ExcitingXMLInput):
    """
    Class for exciting Energywindow Input
    """
    name = "energywindow"


class ExcitingXSQpointsetInput(AbstractExcitingInput):
    """
    Class for exciting Qpointset Input
    """
    name = "qpointset"

    def __init__(self, qpointset: Optional[Union[np.ndarray, List[List[float]]]] = np.array([0.0, 0.0, 0.0])):
        """
        Qpointset should be passed either as numpy array or as a list of lists, so either
        np.array([[0., 0., 0.], [0.0, 0.0, 0.01], ...])
        or
        [[0., 0., 0.], [0.0, 0.0, 0.01], ...]
        """
        self.qpointset = qpointset

    def to_xml(self) -> ElementTree.Element:
        """ Special implementation of to_xml for the qpointset element. """
        qpointset = ElementTree.Element(self.name)
        for qpoint in self.qpointset:
            ElementTree.SubElement(qpointset, 'qpoint').text = list_to_str(qpoint)

        return qpointset


class ExcitingXSPlanInput(AbstractExcitingInput):
    """
    Class for exciting Plan Input
    """
    name = "plan"

    def __init__(self, plan: List[str]):
        """
        Plan doonly elements are passed as a List of strings in the order exciting shall execute them:
            ['bse', 'xseigval', ...]
        """
        check_valid_keys(plan, valid_plan_entries, self.name)
        self.plan = plan

    def to_xml(self) -> ElementTree.Element:
        """ Special implementation of to_xml for the plan element. """
        plan = ElementTree.Element(self.name)
        for task in self.plan:
            ElementTree.SubElement(plan, 'doonly', task=task)

        return plan


class ExcitingXSInput(ExcitingXMLInput):
    """ Class allowing to write attributes to XML."""

    # TODO(Fabian): Add all the other subelements, see http://exciting.wikidot.com/ref:xs
    # Issue 121: https://git.physik.hu-berlin.de/sol/exciting/-/issues/121
    name = "xs"
