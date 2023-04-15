"""Module for class of exciting xs (excited states).
http://exciting.wikidot.com/ref:xs
"""
from typing import Optional, List, Union
from xml.etree import ElementTree

import numpy as np

from excitingtools.input.base_class import AbstractExcitingInput, ExcitingXMLInput
from excitingtools.utils.dict_utils import check_valid_keys
from excitingtools.utils.utils import list_to_str


class ExcitingXSBSEInput(ExcitingXMLInput):
    """
    Class for exciting BSE Input
    """
    name = "BSE"
    _valid_attributes = {'aresbse', 'blocks', 'bsedirsing', 'bsetype', 'checkposdef', 'chibar0', 'chibar0comp',
                         'chibarq', 'coupling', 'cuttype', 'distribute', 'econv', 'eecs', 'efind', 'fbzq',
                         'iqmtrange', 'lmaxdielt', 'measure', 'nexc', 'ngridksub', 'nleblaik', 'nosym', 'nstlbse',
                         'nstlxas', 'outputlevel', 'reducek', 'rgkmax', 'sciavbd', 'sciavqbd', 'sciavqhd',
                         'sciavqwg', 'sciavtype', 'scrherm', 'vkloff', 'writehamhdf5', 'writepotential', 'xas',
                         'xasatom', 'xasedge', 'xasspecies', 'xes'}


class ExcitingXSScreeningInput(ExcitingXMLInput):
    """
    Class for exciting Screening Input
    """
    name = "screening"
    _valid_attributes = {'do', 'intraband', 'nempty', 'ngridk', 'nosym', 'reducek', 'rgkmax', 'screentype',
                         'tr', 'vkloff'}


class ExcitingXSEnergywindowInput(ExcitingXMLInput):
    """
    Class for exciting Energywindow Input
    """
    name = "energywindow"
    _valid_attributes = {'intv', 'points'}


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
    _valid_plan_elements = {'xsgeneigvec', 'tetcalccw', 'writepmatxs', 'writeemat', 'df', 'df2', 'idf', 'scrgeneigvec',
                            'scrtetcalccw', 'scrwritepmat', 'screen', 'scrcoulint', 'exccoulint', 'bse', 'bsegenspec',
                            'writebevec', 'writekpathweights', 'bsesurvey', 'kernxc_bse', 'writebandgapgrid',
                            'write_wfplot', 'write_screen', 'writepmat', 'dielectric', 'writepmatasc', 'pmatxs2orig',
                            'writeoverlapxs', 'writeematasc', 'writepwmat', 'emattest', 'x0toasc', 'x0tobin',
                            'fxc_alda_check', 'kernxc_bse3', 'testxs', 'xsestimate', 'testmain', 'excitonWavefunction',
                            'portstate(1)', 'portstate(2)', 'portstate(-1)', 'portstate(-2)'}

    def __init__(self, plan: List[str]):
        """
        Plan doonly elements are passed as a List of strings in the order exciting shall execute them:
            ['bse', 'xseigval', ...]
        """
        check_valid_keys(plan, self._valid_plan_elements, self.name)
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
    _valid_attributes = {'bfieldc', 'broad', 'dbglev', 'dfoffdiag', 'dogroundstate', 'emattype', 'emaxdf',
                         'epsdfde', 'fastpmat', 'gqmax', 'gqmaxtype', 'lmaxapwwf', 'lmaxemat', 'maxscl', 'nempty',
                         'ngridk', 'ngridq', 'nosym', 'pwmat', 'reducek', 'reduceq', 'rgkmax', 'scissor', 'skipgnd',
                         'swidth', 'tappinfo', 'tevout', 'vkloff', 'writexsgrids', 'xstype'}
    _valid_subtrees = {"screening", "BSE", "qpointset", "energywindow", "plan"}
    _mandatory_keys = {"xstype"}
