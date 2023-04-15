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
    _valid_attributes = {"bfieldc", "fixspin", "momfix", "nosv", "realspace", "reducebf", "spinorb",
                         "spinsprl", "taufsm", "vqlss"}


class ExcitingGSSolverInput(ExcitingXMLInput):
    """
    Class for exciting solver input.
    """
    name = "solver"
    _valid_attributes = {"ArpackImproveInverse", "ArpackLinSolve", "ArpackShift", "ArpackUserDefinedShift",
                         "DecompPrec", "epsarpack", "evaltol", "packedmatrixstorage", "type"}


class ExcitingGroundStateInput(ExcitingXMLInput):
    """
    Class for exciting groundstate input.
    """

    # Reference: http://exciting.wikidot.com/ref:groundstate
    name = 'groundstate'
    _valid_attributes = {'CoreRelativity', 'ExplicitKineticEnergy', 'PrelimLinSteps', 'ValenceRelativity', 'autokpt',
                         'beta0', 'betadec', 'betainc', 'cfdamp', 'chgexs', 'deband', 'dipolecorrection',
                         'dipoleposition', 'dlinengyfermi', 'do', 'energyref', 'epsband', 'epschg', 'epsengy',
                         'epsforcescf', 'epsocc', 'epspot', 'fermilinengy', 'findlinentype', 'fracinr', 'frozencore',
                         'gmaxvr', 'isgkmax', 'ldapu', 'lmaxapw', 'lmaxinr', 'lmaxmat', 'lmaxvr', 'lorecommendation',
                         'lradstep', 'maxscl', 'mixer', 'mixerswitch', 'modifiedsv', 'msecStoredSteps', 'nempty',
                         'ngridk', 'niterconvcheck', 'nktot', 'nosource', 'nosym', 'nprad', 'npsden', 'nwrite',
                         'outputlevel', 'ptnucl', 'radialgridtype', 'radkpt', 'reducek', 'rgkmax', 'scfconv', 'stype',
                         'swidth', 'symmorph', 'tevecsv', 'tfibs', 'tforce', 'tpartcharges', 'useDensityMatrix',
                         'vdWcorrection', 'vkloff', 'xctype'}
    _valid_subtrees = {"spin", "solver"}
