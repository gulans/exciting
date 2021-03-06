"""
API for selecting a parser, given a file name and parsing the data
"""
import os
import sys
from xml.etree.ElementTree import ParseError

from . import groundStateParser
from . import propertiesParser
from . import BSEParser
from . import GWParser
from . import RT_TDDFTParser
from .ErroneousFileError import ErroneousFileError

from excitingtools.dict_utils import container_converter

# Map file name to parser function
_file_to_parser = {
    'INFO.OUT': groundStateParser.parse_info_out,
    'info.xml': groundStateParser.parse_info_xml,
    'atoms.xml': groundStateParser.parse_atoms,
    'evalcore.xml': groundStateParser.parse_evalcore,
    'eigval.xml': groundStateParser.parse_eigval,
    'geometry.xml': groundStateParser.parse_geometry,
    'RHO3D.xml': propertiesParser.parse_plot_3D,
    'VCL3D.xml': propertiesParser.parse_plot_3D,
    'VXC3D.xml': propertiesParser.parse_plot_3D,
    'WF3D.xml': propertiesParser.parse_plot_3D,
    'ELF3D.xml': propertiesParser.parse_plot_3D,
    'EF3D.xml': propertiesParser.parse_plot_3D,
    'LSJ.xml': propertiesParser.parse_LSJ,
    'EFG.xml': propertiesParser.parse_EFG,
    'mossbauer.xml': propertiesParser.parse_mossbauer,
    'expiqr.xml': propertiesParser.parse_expiqr,
    'effmass.xml': propertiesParser.parse_effmass,
    'bandstructure.xml': propertiesParser.parse_bandstructure,
    'dos.xml': propertiesParser.parse_dos,
    'KERR.OUT': propertiesParser.parse_kerr,
    'EPSILON_11.OUT': propertiesParser.parse_epsilon,
    'EPSILON_12.OUT': propertiesParser.parse_epsilon,
    'EPSILON_33.OUT': propertiesParser.parse_epsilon,
    'CHI_111.OUT': propertiesParser.parse_chi,
    'ELNES.OUT': propertiesParser.parse_elnes,
    'SEEBECK_11.OUT': propertiesParser.parse_seebeck,
    'ELECTCOND_11.OUT': propertiesParser.parse_seebeck,
    'THERMALCOND_11.OUT': propertiesParser.parse_seebeck,
    'Z_11.OUT': propertiesParser.parse_seebeck,
    'ldos.out': propertiesParser.parse_ldos,
    'band_edges.out': propertiesParser.parse_band_edges,
    'spintext.xml': propertiesParser.parse_spintext,
    'POLARIZATION.OUT': propertiesParser.parse_polarization,
    'TDOS_WANNIER.OUT': propertiesParser.parse_tdos_wannier,
    'WANNIER_INFO.OUT': propertiesParser.parse_wannier_info,
    'coreoverlap.xml': propertiesParser.parse_core_overlap,  
    'EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT': BSEParser.parse_EPSILON_NAR,
    'EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT': BSEParser.parse_EPSILON_NAR,
    'EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT': BSEParser.parse_EPSILON_NAR,
    'EPSILON_NAR_FXCMB1_OC11_QMT001.OUT': BSEParser.parse_EPSILON_NAR,
    'EPSILON_NAR_FXCMB1_OC22_QMT001.OUT': BSEParser.parse_EPSILON_NAR,
    'EPSILON_NAR_FXCMB1_OC33_QMT001.OUT': BSEParser.parse_EPSILON_NAR,
    'EPSILON_NAR_NLF_FXCMB1_OC11_QMT001.OUT': BSEParser.parse_EPSILON_NAR,
    'EPSILON_NAR_NLF_FXCMB1_OC22_QMT001.OUT': BSEParser.parse_EPSILON_NAR,
    'EPSILON_NAR_NLF_FXCMB1_OC33_QMT001.OUT': BSEParser.parse_EPSILON_NAR,
    'EPSILON_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT': BSEParser.parse_EPSILON_NAR,
    'EPSILON_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT': BSEParser.parse_EPSILON_NAR,
    'EPSILON_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT': BSEParser.parse_EPSILON_NAR,
    'LOSS_NAR_FXCMB1_OC11_QMT001.OUT': BSEParser.parse_LOSS_NAR,
    'LOSS_NAR_FXCMB1_OC22_QMT001.OUT': BSEParser.parse_LOSS_NAR,
    'LOSS_NAR_FXCMB1_OC33_QMT001.OUT': BSEParser.parse_LOSS_NAR,
    'LOSS_NAR_NLF_FXCMB1_OC11_QMT001.OUT': BSEParser.parse_LOSS_NAR,
    'LOSS_NAR_NLF_FXCMB1_OC22_QMT001.OUT': BSEParser.parse_LOSS_NAR,
    'LOSS_NAR_NLF_FXCMB1_OC33_QMT001.OUT': BSEParser.parse_LOSS_NAR,
    'LOSS_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT': BSEParser.parse_LOSS_NAR,
    'LOSS_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT': BSEParser.parse_LOSS_NAR,
    'LOSS_BSE-singlet-TDA-BAR_SCR-full_OC33.out': BSEParser.parse_LOSS_NAR,
    'EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT': BSEParser.parse_EXCITON_NAR_BSE,
    'EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT': BSEParser.parse_EXCITON_NAR_BSE,
    'EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT': BSEParser.parse_EXCITON_NAR_BSE,
    'EXCITON_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT': BSEParser.parse_EXCITON_NAR_BSE,
    'EXCITON_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT': BSEParser.parse_EXCITON_NAR_BSE,
    'EXCITON_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT': BSEParser.parse_EXCITON_NAR_BSE,
    'EFERMI_GW.OUT': GWParser.parse_efermi_gw,
    'EVALQP.DAT': GWParser.parse_evalqp,
    'VXCNN.DAT': GWParser.parse_vxcnn,
    'EPS00_GW.OUT': GWParser.parse_eps00_gw,
    'JIND.OUT': RT_TDDFTParser.parse_jind,
    'NEXC.OUT': RT_TDDFTParser.parse_nexc,
    'ETOT_RTTDDFT.OUT': RT_TDDFTParser.parse_etot,
    'EIGVAL_': RT_TDDFTParser.parse_eigval_screenshots,
    'PROJ_': RT_TDDFTParser.parse_proj_screenshots
}


def parser_chooser(full_file_name: str):
    """
    Selects parser according to the name of the input file
    then returns the result of the parser.

    Could probably implement with decorators.

    param: str, full_file_name: file name prepended by full path
    return: parsed data
    """
    file_name = os.path.split(full_file_name)[1]

    if ('EIGVAL_' in file_name) or ('PROJ_' in file_name):
        file_name = file_name.split('_')[0] + '_'

    files_with_parsers = [name for name in _file_to_parser.keys()]

    if not (file_name in files_with_parsers):
        sys.exit("File does not have a parser:" + file_name)

    if not os.path.exists(full_file_name):
        raise OSError('File path not valid:' + full_file_name)

    parser = _file_to_parser[file_name]
    data = parser(full_file_name)
    return container_converter(data)
