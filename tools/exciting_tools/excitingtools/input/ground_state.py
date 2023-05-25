"""Module for class of exciting ground state.

Ideally the input keywords (class attributes) should be parsed from the schema BUT
because excitingtools will also be available as a standalone package, one would need
to have a copy of the schema XML in excitingtools, which is kept synchronised with
the <EXCITINGROOT>/xml/.
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


class ExcitingGroundStateInput(ExcitingXMLInput):
    """
    Class for exciting groundstate input.
    """
    # Reference: http://exciting.wikidot.com/ref:groundstate
    name = 'groundstate'
