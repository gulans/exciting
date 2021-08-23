"""
Module containing parsers for exciting properties
"""
import xml.etree.cElementTree as ET
from xml.etree.ElementTree import ParseError
import numpy as np
import os


def parse_plot_3d(name: str) -> dict:
    """
    Parser for RHO3D.xml, VCL3D.xml, VXC3D.xml, WF3D.xml, ELF3D.xml, EF3D.xmlit

    :param str name: File name
    :return dict output: Parsed data
    """
    root = ET.parse(name)
    plot_3d = {"title": root.find("title").text}
    grid = root.find("grid").attrib
    axis = []
    for ax in root.find("grid").findall("axis"):
        axis.append(ax.attrib)
    grid["axis"] = {}
    for item in axis:
        name = item["name"]
        grid["axis"][name] = item
    grid["value"] = root.find("grid").find("value").attrib
    plot_3d["grid"] = grid

    func_node = root.find("function")
    function = root.find("function").attrib
    row0 = []
    for node in func_node.findall("row"):
        row = node.attrib
        row1 = []
        for nod in node:
            row2 = nod.attrib
            row2["data"] = nod.text.split()
            row1.append(row2)
        row["row"] = {}
        for item in row1:
            name = item["index"]
            row["row"][name] = item
        row0.append(row)
    function["row"] = {}
    for item in row0:
        name = item["index"]
        function["row"][name] = item
    plot_3d["function"] = function

    return plot_3d


def parse_lsj(name: str) -> dict:
    """
    Parser for LSJ.xml

    :param str name: File name
    :return dict output: Parsed data
    """
    root = ET.parse(name).getroot()
    LSJ = {}
    species = []
    for node in root.findall("species"):
        spec = node.attrib
        atom = []
        for nod in node.findall("atom"):
            at = nod.attrib
            at["L"] = nod.find("L").text.split()
            at["S"] = nod.find("S").text.split()
            at["J"] = nod.find("J").text.split()
            atom.append(at)
        spec["atom"] = {}
        for item in atom:
            name = item["n"]
            spec["atom"][name] = item
        species.append(spec)
    LSJ["species"] = {}
    for item in species:
        name = item['n']
        LSJ['species'][name] = item

    return LSJ


def parse_efg(name: str) -> dict:
    """
    Parser for EFG.xml

    Returns a dictionary of the form:
      data = {'species1': {'chemicalSymbol': chemicalSymbol,
                          'atom1': { 'trace': trace,
                                     'efg: efg,
                                     'eigenvalues': eigenvalues
                                   },
                          'atom2': {},
                          }
              'species2':...
            }

    :param str name: File name
    :return dict output: Parsed data
    """
    # If one were to pass an string containing XML
    # tree = ET.ElementTree(ET.fromstring(xmlstring))
    # root = tree.getroot()

    root = ET.parse(name).getroot()
    data = {}

    for species in root.findall('species'):
        species_key = species.tag + str(species.attrib['n'])
        data[species_key] = {'chemicalSymbol': species.attrib['chemicalSymbol']}

        for atom in species.findall('atom'):
            atom_key = atom.tag + atom.attrib['n']

            for efg_tensor in atom.findall("EFG-tensor"):
                efg = np.empty(shape=(3,3))

                for i, line in enumerate(efg_tensor.findall("line")):
                    efg[i, :] = [float(x) for x in line.text.split()]

                for eigenvalues in efg_tensor.findall("eigenvalues"):
                    eig = [float(e) for e in eigenvalues.text.split()]

            data[species_key][atom_key] = {'trace': float(efg_tensor.attrib['trace']),
                                           'efg': efg,
                                           'eigenvalues': eig
                                           }
    return data


def parse_mossbauer(name: str) -> dict:
    """
    Parser for mossbauer.xml

    :param str name: File name
    :return dict output: Parsed data
    """
    root = ET.parse(name).getroot()
    mossbauer = {}
    species = []
    for node in root.findall("species"):
        spec = node.attrib
        atom = []
        for nod in node:
            atom.append(nod.attrib)
        spec["atom"] = {}
        for item in atom:
            name = item["n"]
            spec["atom"][name] = item
        species.append(spec)
    mossbauer["species"] = {}
    for item in species:
        name = item["n"]
        mossbauer["species"][name] = item

    return mossbauer


def parse_expiqr(name: str) -> dict:
    """
    Parser for expiqr.xml

    :param str name: File name
    :return dict output: Parsed data
    """
    root = ET.parse(name).getroot()
    expiqr = {}
    expiqr["q-vector"] = root.find("q-vector").attrib
    kgrid = {}
    for k in root.find("k-grid"):
        kgrid = k.attrib
        states = []
        for st in k.findall("state"):
            state = st.attrib
            states.append(state)
            statesj = []
            for s in st.findall("state"):
                statej = s.attrib
                statesj.append(statej)
            state["state"] = {}
            for item in statesj:
                name = item["j"]
                state["state"][name] = item
        kgrid["state"] = {}
        for item in states:
            name = item["i"]
            kgrid["state"][name] = item
    expiqr["k-grid"] = kgrid
    return expiqr


def parse_effmass(name: str) -> dict:
    """
    Parser for effmass.xml

    :param str name: File name
    :return dict output: Parsed data
    """
    root = ET.parse(name).getroot()
    effmass = {}
    effmass["k-point"] = root.find("k-point").attrib
    state = []
    for node in root.findall("state"):
        st = node.attrib
        evdk = node.find("evdk").attrib
        matrix1 = []
        for line in node.find("evdk").findall("line"):
            matrix1.append(line.text.split())
        evdk["matrix"] = matrix1
        emt = node.find("emt").attrib
        matrix2 = []
        for line in node.find("emt").findall("line"):
            matrix2.append(line.text.split())
        emt["matrix"] = matrix2
        st["evdk"] = evdk
        st["emt"] = emt
        state.append(st)
    effmass["state"] = {}
    for item in state:
        name = item["n"]
        effmass["state"][name] = item

    return effmass


def parse_bandstructure(name: str) -> dict:
    """
    Parser for bandstructure.xml

    :param str name: File name
    :return dict output: Parsed data
    """
    root = ET.parse(name).getroot()
    bandstructure = {}
    bandstructure["title"] = root.find("title").text
    bands = []
    for node in root.findall("band"):
        band = node.attrib
        points = []
        for nod in node.findall("point"):
            point = nod.attrib
            points.append(point)
        bands.append(band)
        band["point"] = {}
        i = 1
        for item in points:
            name = str(i)
            band["point"][name] = item
            i = i+1
    bandstructure["band"] = {}
    j = 1
    for item in bands:
        name = str(j)
        bandstructure["band"][name] = item
        j = j+1
    return bandstructure


def parse_dos(name: str) -> dict:
    """
    Parser for dos.xml

    :param str name: File name
    :return dict output: Parsed data
    """
    root = ET.parse(name).getroot()
    dos = {}
    dos["title"] = root.find("title").text
    totaldos = root.find("totaldos").attrib
    dos["totaldos"] = totaldos
    diagram  = root.find("totaldos").find("diagram").attrib
    dos["totaldos"]["diagram"] = diagram
    points = []
    for node in root.find("totaldos").find("diagram").findall("point"):
        point = node.attrib
        points.append(point)
    dos["totaldos"]["diagram"]["point"] = {}
    i = 1
    for item in points:
        name = str(i)
        dos["totaldos"]["diagram"]["point"][name] = item
        i = i+1
    return dos


def parse_kerr(name: str) -> dict:
    """
    Parser for KERR.OUT

    :param str name: File name
    :return dict output: Parsed data
    """
    try:
        data = np.genfromtxt(name, skip_header=1)
    except:
        raise ParseError

    out = {'energy': data[:, 0],
           're':     data[:, 1],
           'im':     data[:, 2]}

    return out


def parse_epsilon(name: str) -> dict:
    """
    Parser for EPSILON_ij.OUT

    :param str name: File name
    :return dict output: Parsed data
    """
    try:
        data = np.genfromtxt(name, skip_header=1)
    except:
        raise ParseError
    out = {'energy': data[:, 0], 're': data[:, 1], 'im': data[:, 2]}
    return out


def parse_chi(name: str) -> dict:
    """
    Parser for CHI_111.OUT

    :param str name: File name
    :return dict output: Parsed data
    """
    try:
        data = np.genfromtxt(name, skip_header=1)
    except:
        raise ParseError
    out = {'energy':  data[:, 0],
           're':      data[:, 1],
           'im':      data[:, 2],
           'modulus': data[:, 3]}
    return out


def parse_elnes(name: str) -> dict:
    """
    Parser for ELNES.OUT

    :param str name: File name
    :return dict output: Parsed data
    """
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {'energy': data[:, 0], 'elnes': data[:, 1]}
    return out


def parse_seebeck(name: str) -> dict:
    """
    Parser for SEEBECK_11.OUT

    :param str name: File name
    :return dict output: Parsed data
    """
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {'temperature': data[:, 0],
           'mu': data[:, 1],
           're': data[:, 2],
           'im': data[:, 3]}

    return out


def parse_ldos(name: str) -> dict:
    """
    Parser for ldos.out

    :param str name: File name
    :return dict output: Parsed data
    """
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {'energy': data[:, 0], 'ldos': data[:, 1]}

    return out


def parse_band_edges(name: str) -> dict:
    """
    Parser for band_edges.out

    Keys
    * c_axis corresponds to the linear grid along the magnitude
      of the c vector of the unit cell.
    * VBM = Valence band maximum
    * CBm = Conduction band minimum

    :param str name: File name
    :return dict output: Parsed data
    """
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {'c_axis': data[:, 0], 'VBM': data[:, 1], 'CBm': data[:, 2]}

    return out


def parse_spintext(name: str) -> list:
    """
    Parse spintext.xml

    TODO(Bene) Issue 87 Refactor to return a dict

    Each element of the list contains a dict with keys:
     ['ist', 'k-point', 'spin', 'energy']

    :param str name: Path to the spintext.xml that will be parsed
    :return dict spintext: List that holds the parsed spintexture.xml
    """
    #parse file
    file_name = 'spintext.xml'
    if name.split('/')[-1] != file_name:
        name = os.path.join(name, file_name)

    tree_spin = ET.parse(name)
    root_spin = tree_spin.getroot()

    spintext = []
    for band in root_spin.findall("band"):
        k_point = []
        spin = []
        energy = []

        for val in band.findall("k-point"):
            k_point.append([float(k) for k in val.attrib["vec"].split()])
            spin.append([float(s) for s in val.attrib["spin"].split()])
            energy.append(float(val.attrib["energy"]))

        spintext.append({"ist": int(band.attrib["ist"]),
                         "k-point":k_point,
                         "spin": spin,
                         "energy": energy}
                        )
    return spintext


def parse_polarization(name: str) -> dict:
    """
    Parser for POLARIZATION.OUT

    :param str name: File name
    :return dict output: Parsed data
    """
    file = open(name)
    lines = []
    for i in range(len(open(name).readlines())):
        line = next(file)
        if '#' not in line:
            lines.append(line.split())
    polarization = {'total': lines[0], 'electronic': lines[1], 'ionic': lines[2]}
    return polarization


def parse_tdos_wannier(name: str) -> dict:
    """
    Parser for TDOS_WANNIER.OUT

    :param str name: File name
    :return dict output: Parsed data
    """
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {'energy': data[:, 0], 'dos': data[:, 1]}

    return out


def parse_wannier_info(name: str) -> dict:
    """
    Parser for WANNIER_INFO.OUT

    :param str name: File name
    :return dict output: Parsed data
    """
    file = open(name)

    # Extract data
    lines = []
    data = []
    total = []
    start = False
    for i, line in enumerate(file.readlines()):
        if "* Wannier functions" in line:
            start = True
        if start:
            lines.append(line)
    for i in range(len(lines)):
        if lines[i].strip().startswith( '1'):
            for j in range (4):
                data.append(lines[i+j].split())
        if lines[i].strip().startswith( '5'):
            for j in range (4):
                data.append(lines[i+j].split())
        if lines[i].strip().startswith( 'total'):
            total.append(lines[i].split())
    file.close()

    # Package data into dictionary
    n_wannier = len(data)
    localisation_center = np.empty(shape=(n_wannier,3))
    wannier = {'n_wannier': n_wannier, 'Omega': [], 'Omega_I': [], 'Omega_D': [], 'Omega_OD': []}

    for i, item in enumerate(data):
        localisation_center[i, :] = [float(x) for x in item[1:4]]
        wannier['Omega'].append(float(item[4]))
        wannier['Omega_I'].append(float(item[5]))
        wannier['Omega_D'].append(float(item[6]))
        wannier['Omega_OD'].append(float(item[7]))

    wannier['localisation_center'] = localisation_center

    totals = {'Omega': [], 'Omega_I': [], 'Omega_D': [], 'Omega_OD': []}
    for j, item in enumerate(total):
        totals['Omega'].append(float(item[1]))
        totals['Omega_I'].append(float(item[2]))
        totals['Omega_D'].append(float(item[3]))
        totals['Omega_OD'].append(float(item[4]))

    wannier['total'] = totals

    return wannier


def parse_core_overlap(name: str) -> dict:
    """
    Parser for coreoverlap.xml

    Parsed dictionary has the structure:

        output = {'nkpt':  nkpt
                  'nstfv': nstfv
                  'ncg':   ncg
                  'kpoints': [{'index': index, 'pairs': pairs},
                              {'index': index, 'pairs': pairs},
                              ...]
                  }

    where output['kpoints'][ik]['pairs'] =
      [{'ist1': '1', 'ist2': '1', 'de': 12.97849772, 'overlap': 3.35753859e-07},
       {'ist1': '1', 'ist2': '2', 'de': 12.97849772, 'overlap': 3.35753859e-07},
       ...
       n_pairs]

    :param str name: File name
    :return dict output: Parsed data
    """
    try:
        tree = ET.parse(name)
    except:
        raise ParseError

    root = tree.getroot()
    core_overlap = {'nkpt': int(root.attrib['nkpt']),
                    'nstfv': int(root.attrib['nstfv']),
                    'ncg': int(root.attrib['ncg'])}

    k_points = []
    for k_point in root:
        kpt = {"index": int(k_point.attrib['index'])}
        pairs = []
        for pair_xml in k_point:
            pair = pair_xml.attrib
            pair['ist1'] = int(pair['ist1'])
            pair['ist2'] = int(pair['ist2'])
            pair["de"] = float(pair["de"])
            pair["overlap"] = float(pair["overlap"].split()[0])**2 + float(pair["overlap"].split()[1])**2
            pairs.append(pair)
        kpt["pairs"] = pairs
        k_points.append(kpt)
    core_overlap["kpoints"] = k_points

    return core_overlap
