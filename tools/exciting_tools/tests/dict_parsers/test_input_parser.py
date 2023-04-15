"""
Test for the input.xml file parser
"""
import pytest

from excitingtools.exciting_dict_parsers.input_parser import parse_element_xml, parse_structure, parse_input_xml

reference_input_str = """<?xml version="1.0" encoding="UTF-8"?>
<input>
  
  <title>Lithium Fluoride BSE</title>
  
  <structure speciespath="." autormt="false" epslat="1.0d-6">
    <crystal scale="1.0" stretch="1.0">
      <basevect>3.80402 3.80402 0.00000</basevect>
      <basevect>3.80402 0.00000 3.80402</basevect>
      <basevect>0.00000 3.80402 3.80402</basevect>
    </crystal>
    <species speciesfile="Li.xml" rmt="1.5">
      <atom coord="0.0000  0.0000  0.0000" bfcmt="0.0 0.0 0.0"/>
    </species>
    <species speciesfile="F.xml">
      <atom coord="0.5000  0.5000  0.5000" lockxyz="false true false"/>
    </species>
  </structure>
  
  <groundstate xctype="GGA_PBE" ngridk="4 4 4" epsengy="1d-7" outputlevel="high">
  </groundstate>

  <xs xstype="BSE" 
   ngridq="3 3 3"
   vkloff="0.05 0.15 0.25" 
   nempty="1"
   broad="0.0073499"
   nosym="true">

   <energywindow intv="0.0 1.0" 
    points="50" />

   <screening screentype="full"
    nempty="115" />

   <BSE bsetype="singlet"
    nstlbse="1 5 1 2" 
    aresbse="false"/>

   <qpointset>
      <qpoint>0.0 0.0 0.0</qpoint>
   </qpointset>
   
   <plan>
      <doonly task="screen" />
      <doonly task="bse" />
   </plan>
  </xs>
  
</input>
"""


def test_parse_title():
    assert parse_element_xml(reference_input_str, tag="title") == "Lithium Fluoride BSE"


def test_parse_groundstate():
    ground_state = parse_element_xml(reference_input_str, tag="groundstate")
    assert ground_state == {
        'xctype': 'GGA_PBE', 'ngridk': [4, 4, 4],
        'epsengy': 1e-7, 'outputlevel': 'high'
    }


def test_parse_groundstate_from_gs_root():
    ground_state = parse_element_xml('<groundstate xctype="GGA_PBE" ngridk="4 4 4" epsengy="1d-7" outputlevel="high"/>',
                                     tag="groundstate")
    assert ground_state == {
        'xctype': 'GGA_PBE', 'ngridk': [4, 4, 4],
        'epsengy': 1e-7, 'outputlevel': 'high'
    }


def test_parse_structure():
    structure = parse_structure(reference_input_str)
    structure_ref = {
        'atoms': [{'species': 'Li', 'position': [0.0, 0.0, 0.0],
                   'bfcmt': '0.0 0.0 0.0'},
                  {'species': 'F', 'position': [0.5, 0.5, 0.5],
                   'lockxyz': 'false true false'}],
        'lattice': [[3.80402, 3.80402, 0.0],
                    [3.80402, 0.0, 3.80402],
                    [0.0, 3.80402, 3.80402]],
        'species_path': '.',
        'crystal_properties': {'scale': '1.0', 'stretch': '1.0'},
        'species_properties': {'Li': {'rmt': '1.5'}, 'F': {}},
        'autormt': 'false',
        'epslat': '1.0d-6',
    }
    assert structure_ref == structure


def test_parse_xs():
    xs = parse_element_xml(reference_input_str, tag="xs")
    xs_ref = {
        'xstype': 'BSE',
        'ngridq': [3, 3, 3],
        'vkloff': [0.05, 0.15, 0.25],
        'nempty': 1,
        'broad': 0.0073499,
        'nosym': True,
        'energywindow': {'intv': [0.0, 1.0], 'points': 50},
        'screening': {'screentype': 'full', 'nempty': 115},
        'BSE': {'bsetype': 'singlet', 'nstlbse': [1, 5, 1, 2], 'aresbse': False},
        'qpointset': [[0.0, 0.0, 0.0]],
        'plan': ['screen', 'bse']
    }
    assert xs_ref == xs
    assert isinstance(xs["ngridq"][0], int)


def test_parse_input_xml():
    parsed_data = parse_element_xml(reference_input_str)
    assert set(parsed_data.keys()) == {'title', 'groundstate', 'structure', 'xs'}


def test_parse_input_xml_directly():
    parsed_data = parse_input_xml(reference_input_str)
    assert set(parsed_data.keys()) == {'title', 'groundstate', 'structure', 'xs'}


def test_parse_missing_tag():
    with pytest.raises(ValueError, match="Your specified input has no tag missing_tag"):
        parse_element_xml(reference_input_str, "missing_tag")
