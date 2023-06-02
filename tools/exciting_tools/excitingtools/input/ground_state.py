"""Module for class of exciting ground state.

The ground state classes could be generated automatically.
"""
from excitingtools.input.base_class import ExcitingXMLInput


class ExcitingGSSpinInput(ExcitingXMLInput):
    """
    Class for exciting spin input.
    """
    name = "spin"


class ExcitingGSSolverInput(ExcitingXMLInput):
    """
    Class for exciting solver input.
    """
    name = "solver"


class ExcitingGSDFTD2parametersInput(ExcitingXMLInput):
    """
    Class for exciting DFTD2parameters input.
    """
    name = "DFTD2parameters"


class ExcitingGSTSvdWparametersInput(ExcitingXMLInput):
    """
    Class for exciting TSvdWparameters input.
    """
    name = "TSvdWparameters"


class ExcitingGSDfthalfInput(ExcitingXMLInput):
    """
    Class for exciting dfthalf input.
    """
    name = "dfthalf"


class ExcitingGSHybridInput(ExcitingXMLInput):
    """
    Class for exciting Hybrid input.
    """
    name = "Hybrid"


class ExcitingGSSiriusInput(ExcitingXMLInput):
    """
    Class for exciting sirius input.
    """
    name = "sirius"


class ExcitingGSOEPInput(ExcitingXMLInput):
    """
    Class for exciting OEP input.
    """
    name = "OEP"


class ExcitingGSOutputInput(ExcitingXMLInput):
    """
    Class for exciting output input.
    """
    name = "output"


class ExcitingGSLibxcInput(ExcitingXMLInput):
    """
    Class for exciting libxc input.
    """
    name = "libxc"


class ExcitingGroundStateInput(ExcitingXMLInput):
    """
    Class for exciting groundstate input.
    """
    # Reference: http://exciting.wikidot.com/ref:groundstate
    name = 'groundstate'
