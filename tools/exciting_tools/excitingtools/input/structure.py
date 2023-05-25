"""Structure class, mirroring that of exciting's structure XML sub-tree.
https://exciting.wikidot.com/ref:structure
"""
from pathlib import Path
from typing import Optional, Union, List, Dict
from xml.etree import ElementTree

from excitingtools.input.base_class import ExcitingXMLInput
from excitingtools.structure.lattice import check_lattice, check_lattice_vector_norms
from excitingtools.utils import valid_attributes
from excitingtools.utils.dict_utils import check_valid_keys
from excitingtools.utils.utils import list_to_str

# Set of all elements
all_species = {'Ni', 'La', 'K', 'Xe', 'Ag', 'Bk', 'Co', 'Md', 'Lu', 'Ar',
               'Bi', 'Cm', 'H', 'Yb', 'Zn', 'Te', 'I', 'Cl', 'As', 'Mg',
               'No', 'Ta', 'N', 'Ac', 'Y', 'At', 'Tb', 'Tc', 'Au', 'O',
               'Lr', 'In', 'Ge', 'Re', 'Pm', 'Gd', 'Kr', 'Po', 'Sc', 'Rf',
               'Sb', 'Rb', 'Ru', 'Dy', 'Ho', 'Ra', 'Se', 'Sr', 'Fr', 'Ga',
               'Fe', 'Es', 'Si', 'Pr', 'Pd', 'Er', 'Rn', 'Ir', 'He', 'Eu',
               'Pt', 'Pu', 'Sn', 'Pb', 'Hf', 'Fm', 'Rh', 'Sm', 'Pa', 'Hg',
               'Os', 'B', 'U', 'Zr', 'Cf', 'C', 'Na', 'Li', 'Mo', 'Cs',
               'Al', 'V', 'Cd', 'Tm', 'Tl', 'Ba', 'Ce', 'W', 'Am', 'Cr',
               'Nb', 'Mn', 'S', 'Ca', 'Be', 'Br', 'Th', 'Ti', 'Np', 'Ne',
               'P', 'Cu', 'F', 'Nd'}


class ExcitingStructureCrystalInput(ExcitingXMLInput):
    """
    Class for exciting structure crystal input.
    """
    name = "crystal"


class ExcitingStructureSpeciesInput(ExcitingXMLInput):
    """
    Class for exciting structure species input.
    """
    name = "species"


class ExcitingStructure(ExcitingXMLInput):
    """ Class allowing exciting XML structure to be written from python data.

    TODO(Fabian/Alex) 117. Implement all remaining attributes:
     All elements are species-specific. They should be passed like:
     species_properties = {'S': {'LDAplusU':{'J': J, 'U': U, 'l': l}} }
     Element: LDAplusU: J, U, l
     Element: dfthalfparam: ampl, cut, exponent
     Element: shell: ionization, number
    """
    name = "structure"
    # Path type
    path_type = Union[str, Path]

    # Mandatory attribute "coord" taken out because it's specified inside the atoms
    _valid_atom_attributes = set(valid_attributes.atom_valid_attributes) - {"coord"}

    def __init__(self,
                 atoms,
                 lattice: Optional[list] = None,
                 species_path: Optional[path_type] = './',
                 crystal_properties: Optional[Union[dict, ExcitingStructureCrystalInput]] = None,
                 species_properties: Optional[Dict[str, Union[dict, ExcitingStructureSpeciesInput]]] = None,
                 **kwargs):
        """ Initialise instance of ExcitingStructure.

        TODO(Alex) Issue 117. Create our own class with a subset of methods common to ASE' Atom()
          Then we can have a single API for this init. If ASE is used, xAtom() is just a wrapper of
          Atom(), else we have some light methods.
        TODO(Alex/Fabian) Issue 117.
          structure_attributes and crystal_attributes could equally be kwargs.
          Consider changing or extending before the first major version.

        :param atoms: Atoms object of type ase.atoms.Atoms or of the form List[dict], for example:
         atoms = [{'species': 'X', 'position': [x, y, z]}, ...].
         Each dict can also optionally contain the _valid_atom_attributes:
         {'species': 'X', 'position': [x, y, z],
           'bfcmt': [bx, by, bz], 'lockxyz': [lx, ly, lz], 'mommtfix': [mx, my, mz]}.
        If atoms are defined with ASE, optional atomic_properties cannot be specified.
        Eventually, the list of atoms will be replaced with our custom class, which will extend ase.Atoms()
        with the additional, optional attributes.

        :param lattice [a, b, c], where a, b and c are lattice vectors with 3 components.
         For example, a = [ax, ay, az]. Only required if one does not pass an ase Atoms object.
        :param species_path: Optional path to the location of species file/s.
        :param crystal_properties: Optional crystal properties. See _valid_crystal_attributes
        :param species_properties: Optional species properties, defined as:
        {'species1': {'rmt': rmt_value}, 'species2': {'rmt': rmt_value}}
        :param kwargs: Optional structure properties. Passed as kwargs. See _valid_attributes
        """
        if isinstance(species_path, Path):
            species_path = species_path.as_posix()
        super().__init__(speciespath=species_path, **kwargs)
        self.structure_attributes = {"speciespath": species_path, **kwargs}

        if isinstance(atoms, list) and lattice is None:
            raise ValueError("If atoms is a list, lattice must be passed as a separate argument.")

        # Simple container for atoms, as documented in __init__.
        if isinstance(atoms, list):
            check_lattice(lattice)
            check_lattice_vector_norms(lattice)
            self.lattice = lattice
            self.species = [atom['species'].capitalize() for atom in atoms]
            self.positions = [atom['position'] for atom in atoms]
            self.atom_properties = self._init_atom_properties(atoms)
        else:
            self.lattice, self.species, self.positions = self._init_lattice_species_positions_from_ase_atoms(atoms)
            self.atom_properties = [{}] * len(self.species)

        self.unique_species = sorted(set(self.species))

        # Catch symbols that are not valid elements
        check_valid_keys(self.unique_species, all_species, name='Species input')

        # Optional properties
        self.crystal_properties = self._initialise_subelement_attribute(ExcitingStructureCrystalInput,
                                                                        crystal_properties or {})
        self.species_properties = dict(self._init_species_properties(species_properties))

    @staticmethod
    def _init_lattice_species_positions_from_ase_atoms(atoms) -> tuple:
        """ Initialise lattice, species and positions from an ASE Atoms Object.

        Duck typing for atoms, such that ASE is not a hard dependency.

        :param atoms: ASE Atoms object.
        :return  Lattice, species, positions: Lattice, species and positions
        """
        try:
            cell = atoms.get_cell()
            # Convert to consistent form, [a, b, c], where a = [ax, ay, az]
            lattice = [list(cell[i, :]) for i in range(0, 3)]
            species = [x.capitalize() for x in atoms.get_chemical_symbols()]
            positions = atoms.get_positions()
            return lattice, species, positions
        except AttributeError:
            message = "atoms must either be an ase.atoms.Atoms object or List[dict], of the form" \
                      "[{'species': 'X', 'position': [x, y, z]}, ...]."
            raise AttributeError(message)

    def _init_atom_properties(self, atoms: List[dict]) -> List[dict]:
        """ Initialise atom_properties.

        For atoms that contain optional atomic properties, store them as
        dicts in a list of len(n_atoms). Atoms with none of these properties
        will be defined as empty dicts.

        For each element of atoms, one must have  {'species': 'X', 'position': [x, y, z]}  and
        may have the additional attributes: {'bfcmt': [bx, by, bz], 'lockxyz': [lx, ly, lz], 'mommtfix': [mx, my, mz]}.
        Extract the optional attributes and return in `atom_properties`, with string values.

        :param atoms: List container.
        :return atom_properties: List of atom properties. Each element is a dict.
        and the dict value has been converted to string - ready for XML usage.
        """
        atom_properties: List[dict] = []

        for atom in atoms:
            optional_property_keys = set(atom.keys()) & self._valid_atom_attributes
            optional_atom = {key: atom[key] for key in optional_property_keys}
            optional_properties = {}
            for key, value in optional_atom.items():
                optional_properties[key] = self._attributes_to_input_str[type(value)](value)
            atom_properties.append(optional_properties)

        return atom_properties

    def _init_species_properties(self, species_properties: dict) -> Dict[str, ExcitingXMLInput]:
        """ Initialise species_properties.

        For species without properties, return empty_properties: {'S': {}, 'Al': {}}.

        :param species_properties: Species properties
        :return Dict of ExitingXMLInput-species_properties.
        """
        if species_properties is None:
            species_properties = {}

        for species in self.unique_species:
            props = species_properties.get(species) or {}
            props["speciesfile"] = species + '.xml'
            yield species, self._initialise_subelement_attribute(ExcitingStructureSpeciesInput, props)

    def _group_atoms_by_species(self) -> dict:
        """Get the atomic indices for atoms of each species.

        For example, for:
          species = ['Cd', 'S', 'Cd]
        return:
          indices = {'Cd': [1, 3], 'S' : [2]}

        :return dict indices: Indices of atoms in species and positions
        """
        indices = {}
        for x in self.unique_species:
            indices[x] = [i for i, element in enumerate(self.species) if element == x]
        return indices

    def _xml_atomic_subtree(self, x: str, species, atomic_indices):
        """ Add the required atomic positions and any optional attributes, per species.

        :param x: Species
        :param species: Empty SubElement for species x, which gets filled
        """
        for iatom in atomic_indices[x]:
            coord_str = list_to_str(self.positions[iatom])
            ElementTree.SubElement(species, "atom", coord=coord_str, **self.atom_properties[iatom]).text = ' '

    def to_xml(self) -> ElementTree.Element:
        """Convert structure attributes to XML ElementTree
        Makes use of the to_xml() function of the ExitingXMLInput class to convert values to string.

        Expect to return an XML structure which looks like:
          <structure speciespath="./">

           <crystal scale="1.00" scale="1.00" >
             <basevect>1.0 1.0 0.0</basevect>
             <basevect>1.0 0.0 1.0</basevect>
             <basevect>0.0 1.0 1.0</basevect>
           </crystal>

           <species speciesfile="Al.xml">
             <atom coord="0.0  0.0  0.0"> </atom>
           </species>

          </structure>

        :return ET structure: Element tree containing structure attributes.
        """
        structure_attributes = {
            key: self._attributes_to_input_str[type(value)](value) for key, value in self.structure_attributes.items()
        }
        structure = ElementTree.Element(self.name, **structure_attributes)
        structure.text = ' '

        # Lattice vectors
        crystal = self.crystal_properties.to_xml()
        structure.append(crystal)
        for vector in self.lattice:
            ElementTree.SubElement(crystal, "basevect").text = list_to_str(vector)

        # Species tags
        atomic_indices = self._group_atoms_by_species()
        for x in self.unique_species:
            species = self.species_properties[x].to_xml()
            structure.append(species)
            self._xml_atomic_subtree(x, species, atomic_indices)

        return structure
