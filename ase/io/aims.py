import os
import time
import warnings

import numpy as np
from pathlib import Path

from ase import Atoms, Atom
from ase.calculators.calculator import kpts2mp
from ase.calculators.singlepoint import SinglePointDFTCalculator
from ase.constraints import FixAtoms, FixCartesian
from ase.data import atomic_numbers
from ase.units import Ang, fs
from ase.utils import reader, writer, lazyproperty


v_unit = Ang / (1000.0 * fs)


# Read aims geometry files
@reader
def read_aims(fd, apply_constraints=True):
    """Import FHI-aims geometry type files.

    Reads unitcell, atom positions and constraints from
    a geometry.in file.

    If geometric constraint (symmetry parameters) are in the file
    include that information in atoms.info["symmetry_block"]
    """

    lines = fd.readlines()
    return parse_geometry_lines(lines, apply_constraints=apply_constraints)


def parse_geometry_lines(lines, apply_constraints=True):

    from ase import Atoms
    from ase.constraints import (
        FixAtoms,
        FixCartesian,
        FixScaledParametricRelations,
        FixCartesianParametricRelations,
    )

    atoms = Atoms()

    positions = []
    cell = []
    symbols = []
    velocities = []
    magmoms = []
    symmetry_block = []
    charges = []
    fix = []
    fix_cart = []
    xyz = np.array([0, 0, 0])
    i = -1
    n_periodic = -1
    periodic = np.array([False, False, False])
    cart_positions, scaled_positions = False, False
    for line in lines:
        inp = line.split()
        if inp == []:
            continue
        if inp[0] in ["atom", "atom_frac"]:

            if inp[0] == "atom":
                cart_positions = True
            else:
                scaled_positions = True

            if xyz.all():
                fix.append(i)
            elif xyz.any():
                fix_cart.append(FixCartesian(i, xyz))
            floatvect = float(inp[1]), float(inp[2]), float(inp[3])
            positions.append(floatvect)
            symbols.append(inp[4])
            magmoms.append(0.0)
            charges.append(0.0)
            xyz = np.array([0, 0, 0])
            i += 1

        elif inp[0] == "lattice_vector":
            floatvect = float(inp[1]), float(inp[2]), float(inp[3])
            cell.append(floatvect)
            n_periodic = n_periodic + 1
            periodic[n_periodic] = True

        elif inp[0] == "initial_moment":
            magmoms[-1] = float(inp[1])

        elif inp[0] == "initial_charge":
            charges[-1] = float(inp[1])

        elif inp[0] == "constrain_relaxation":
            if inp[1] == ".true.":
                fix.append(i)
            elif inp[1] == "x":
                xyz[0] = 1
            elif inp[1] == "y":
                xyz[1] = 1
            elif inp[1] == "z":
                xyz[2] = 1

        elif inp[0] == "velocity":
            floatvect = [v_unit * float(l) for l in inp[1:4]]
            velocities.append(floatvect)

        elif inp[0] in [
            "symmetry_n_params",
            "symmetry_params",
            "symmetry_lv",
            "symmetry_frac",
        ]:
            symmetry_block.append(" ".join(inp))

    if xyz.all():
        fix.append(i)
    elif xyz.any():
        fix_cart.append(FixCartesian(i, xyz))

    if cart_positions and scaled_positions:
        raise Exception(
            "Can't specify atom positions with mixture of "
            "Cartesian and fractional coordinates"
        )
    elif scaled_positions and periodic.any():
        atoms = Atoms(symbols, scaled_positions=positions, cell=cell, pbc=periodic)
    else:
        atoms = Atoms(symbols, positions)

    if len(velocities) > 0:
        if len(velocities) != len(positions):
            raise Exception("Number of positions and velocities have to coincide.")
        atoms.set_velocities(velocities)

    fix_params = []

    if len(symmetry_block) > 5:
        params = symmetry_block[1].split()[1:]

        lattice_expressions = []
        lattice_params = []

        atomic_expressions = []
        atomic_params = []

        n_lat_param = int(symmetry_block[0].split(" ")[2])

        lattice_params = params[:n_lat_param]
        atomic_params = params[n_lat_param:]

        for ll, line in enumerate(symmetry_block[2:]):
            expression = " ".join(line.split(" ")[1:])
            if ll < 3:
                lattice_expressions += expression.split(",")
            else:
                atomic_expressions += expression.split(",")

        fix_params.append(
            FixCartesianParametricRelations.from_expressions(
                list(range(3)),
                lattice_params,
                lattice_expressions,
                use_cell=True,
            )
        )

        fix_params.append(
            FixScaledParametricRelations.from_expressions(
                list(range(len(atoms))), atomic_params, atomic_expressions
            )
        )

    if any(magmoms):
        atoms.set_initial_magnetic_moments(magmoms)
    if any(charges):
        atoms.set_initial_charges(charges)

    if periodic.any():
        atoms.set_cell(cell)
        atoms.set_pbc(periodic)
    if len(fix):
        atoms.set_constraint([FixAtoms(indices=fix)] + fix_cart + fix_params)
    else:
        atoms.set_constraint(fix_cart + fix_params)

    if fix_params and apply_constraints:
        atoms.set_positions(atoms.get_positions())
    return atoms


# Write aims geometry files
@writer
def write_aims(
    fd,
    atoms,
    scaled=False,
    geo_constrain=False,
    velocities=False,
    ghosts=None,
    info_str=None,
    wrap=False,
):
    """Method to write FHI-aims geometry files.

    Writes the atoms positions and constraints (only FixAtoms is
    supported at the moment).

    Args:
        fd: file object
            File to output structure to
        atoms: ase.atoms.Atoms
            structure to output to the file
        scaled: bool
            If True use fractional coordinates instead of Cartesian coordinates
        symmetry_block: list of str
            List of geometric constraints as defined in:
            https://arxiv.org/abs/1908.01610
        velocities: bool
            If True add the atomic velocity vectors to the file
        ghosts: list of Atoms
            A list of ghost atoms for the system
        info_str: str
            A string to be added to the header of the file
        wrap: bool
            Wrap atom positions to cell before writing
    """

    from ase.constraints import FixAtoms, FixCartesian

    if scaled and not np.all(atoms.pbc):
        warnings.warn("Setting scaled to False for molecular system.")
        scaled = False

    if geo_constrain:
        if not scaled and np.all(atoms.pbc):
            warnings.warn(
                "Setting scaled to True because a symmetry_block is detected."
            )
            scaled = True
        elif not np.all(atoms.pbc):
            warnings.warn(
                "Parameteric constraints can only be used in periodic systems."
            )
            geo_constrain = False

    fd.write("#=======================================================\n")
    if hasattr(fd, "name"):
        fd.write("# FHI-aims file: " + fd.name + "\n")
    fd.write("# Created using the Atomic Simulation Environment (ASE)\n")
    fd.write("# " + time.asctime() + "\n")

    # If writing additional information is requested via info_str:
    if info_str is not None:
        fd.write("\n# Additional information:\n")
        if isinstance(info_str, list):
            fd.write("\n".join(["#  {}".format(s) for s in info_str]))
        else:
            fd.write("# {}".format(info_str))
        fd.write("\n")

    fd.write("#=======================================================\n")

    i = 0
    if atoms.get_pbc().any():
        for n, vector in enumerate(atoms.get_cell()):
            fd.write("lattice_vector ")
            for i in range(3):
                fd.write("%16.16f " % vector[i])
            fd.write("\n")
    fix_cart = np.zeros([len(atoms), 3])

    if atoms.constraints:
        for constr in atoms.constraints:
            if isinstance(constr, FixAtoms):
                fix_cart[constr.index] = [1, 1, 1]
            elif isinstance(constr, FixCartesian):
                fix_cart[constr.index] = -constr.mask.astype(int) + 1

    if ghosts is None:
        ghosts = np.zeros(len(atoms))
    else:
        assert len(ghosts) == len(atoms)

    wrap = wrap and not geo_constrain
    scaled_positions = atoms.get_scaled_positions(wrap=wrap)

    for i, atom in enumerate(atoms):
        if ghosts[i] == 1:
            atomstring = "empty "
        elif scaled:
            atomstring = "atom_frac "
        else:
            atomstring = "atom "
        fd.write(atomstring)
        if scaled:
            for pos in scaled_positions[i]:
                fd.write("%16.16f " % pos)
        else:
            for pos in atom.position:
                fd.write("%16.16f " % pos)
        fd.write(atom.symbol)
        fd.write("\n")
        # (1) all coords are constrained:
        if fix_cart[i].all():
            fd.write("    constrain_relaxation .true.\n")
        # (2) some coords are constrained:
        elif fix_cart[i].any():
            xyz = fix_cart[i]
            for n in range(3):
                if xyz[n]:
                    fd.write("    constrain_relaxation %s\n" % "xyz"[n])
        if atom.charge:
            fd.write("    initial_charge %16.6f\n" % atom.charge)
        if atom.magmom:
            fd.write("    initial_moment %16.6f\n" % atom.magmom)

        # Write velocities if this is wanted
        if velocities and atoms.get_velocities() is not None:
            fd.write(
                "    velocity {:.16f} {:.16f} {:.16f}\n".format(
                    *atoms.get_velocities()[i] / v_unit
                )
            )

    if geo_constrain:
        for line in get_sym_block(atoms):
            fd.write(line)


def get_sym_block(atoms):
    """Get symmetry block for Parametric constraints in atoms.constraints"""
    from ase.constraints import (
        FixScaledParametricRelations,
        FixCartesianParametricRelations,
    )

    # Initialize param/expressions lists
    atomic_sym_params = []
    lv_sym_params = []
    atomic_param_constr = np.zeros((len(atoms),), dtype="<U100")
    lv_param_constr = np.zeros((3,), dtype="<U100")

    # Populate param/expressions list
    for constr in atoms.constraints:
        if isinstance(constr, FixScaledParametricRelations):
            atomic_sym_params += constr.params

            if np.any(atomic_param_constr[constr.indices] != ""):
                warnings.warn(
                    "multiple parametric constraints defined for the same "
                    "atom, using the last one defined"
                )

            atomic_param_constr[constr.indices] = [
                ", ".join(expression) for expression in constr.expressions
            ]
        elif isinstance(constr, FixCartesianParametricRelations):
            lv_sym_params += constr.params

            if np.any(lv_param_constr[constr.indices] != ""):
                warnings.warn(
                    "multiple parametric constraints defined for the same "
                    "lattice vector, using the last one defined"
                )

            lv_param_constr[constr.indices] = [
                ", ".join(expression) for expression in constr.expressions
            ]

    if np.all(atomic_param_constr == "") and np.all(lv_param_constr == ""):
        return []

    # Check Constraint Parameters
    if len(atomic_sym_params) != len(np.unique(atomic_sym_params)):
        warnings.warn(
            "Some parameters were used across constraints, they will be "
            "combined in the aims calculations"
        )
        atomic_sym_params = np.unique(atomic_sym_params)

    if len(lv_sym_params) != len(np.unique(lv_sym_params)):
        warnings.warn(
            "Some parameters were used across constraints, they will be "
            "combined in the aims calculations"
        )
        lv_sym_params = np.unique(lv_sym_params)

    if np.any(atomic_param_constr == ""):
        raise IOError(
            "FHI-aims input files require all atoms have defined parametric "
            "constraints"
        )

    cell_inds = np.where(lv_param_constr == "")[0]
    for ind in cell_inds:
        lv_param_constr[ind] = "{:.16f}, {:.16f}, {:.16f}".format(*atoms.cell[ind])

    n_atomic_params = len(atomic_sym_params)
    n_lv_params = len(lv_sym_params)
    n_total_params = n_atomic_params + n_lv_params

    sym_block = []
    if n_total_params > 0:
        sym_block.append("#" + "=" * 55 + "\n")
        sym_block.append("# Parametric constraints\n")
        sym_block.append("#" + "=" * 55 + "\n")
        sym_block.append(
            "symmetry_n_params {:d} {:d} {:d}\n".format(
                n_total_params, n_lv_params, n_atomic_params
            )
        )
        sym_block.append(
            "symmetry_params %s\n" % " ".join(lv_sym_params + atomic_sym_params)
        )

        for constr in lv_param_constr:
            sym_block.append("symmetry_lv {:s}\n".format(constr))

        for constr in atomic_param_constr:
            sym_block.append("symmetry_frac {:s}\n".format(constr))
    return sym_block


# Write aims control.in files
@writer
def write_control(fd, atoms, parameters, debug=False):
    parameters = dict(parameters)
    lim = "#" + "=" * 79

    if parameters["xc"] == "LDA":
        parameters["xc"] = "pw-lda"

    cubes = parameters.pop("cubes", None)

    fd.write(lim + "\n")
    for line in [
        "FHI-aims file",
        "Created using the Atomic Simulation Environment (ASE)",
        time.asctime(),
    ]:
        fd.write("# " + line + "\n")
    if debug:
        fd.write("# \n# List of parameters used to initialize the calculator:")
        for p, v in parameters.items():
            s = "#     {}:{}\n".format(p, v)
            fd.write(s)
    fd.write(lim + "\n")

    assert not ("kpts" in parameters and "k_grid" in parameters)
    assert not ("smearing" in parameters and "occupation_type" in parameters)

    for key, value in parameters.items():
        if key == "kpts":
            mp = kpts2mp(atoms, parameters["kpts"])
            fd.write("%-35s%d %d %d\n" % (("k_grid",) + tuple(mp)))
            dk = 0.5 - 0.5 / np.array(mp)
            fd.write("%-35s%f %f %f\n" % (("k_offset",) + tuple(dk)))
        elif key == "species_dir":
            continue
        elif key == "plus_u":
            continue
        elif key == "smearing":
            name = parameters["smearing"][0].lower()
            if name == "fermi-dirac":
                name = "fermi"
            width = parameters["smearing"][1]
            fd.write("%-35s%s %f" % ("occupation_type", name, width))
            if name == "methfessel-paxton":
                order = parameters["smearing"][2]
                fd.write(" %d" % order)
            fd.write("\n")
        elif key == "output":
            for output_type in value:
                fd.write("%-35s%s\n" % (key, output_type))
        elif key == "vdw_correction_hirshfeld" and value:
            fd.write("%-35s\n" % key)
        elif isinstance(value, bool):
            fd.write("%-35s.%s.\n" % (key, str(value).lower()))
        elif isinstance(value, (tuple, list)):
            fd.write("%-35s%s\n" % (key, " ".join(str(x) for x in value)))
        elif isinstance(value, str):
            fd.write("%-35s%s\n" % (key, value))
        else:
            fd.write("%-35s%r\n" % (key, value))

    if cubes:
        cubes.write(fd)

    fd.write(lim + "\n\n")
    write_species(fd, atoms, parameters)


def translate_tier(tier):
    if tier.lower() == "first":
        return 1
    elif tier.lower() == "second":
        return 2
    elif tier.lower() == "third":
        return 3
    elif tier.lower() == "fourth":
        return 4
    else:
        return -1


def write_species(fd, atoms, parameters):
    parameters = dict(parameters)
    species_path = parameters.get("species_dir")
    if species_path is None:
        species_path = os.environ.get("AIMS_SPECIES_DIR")
    if species_path is None:
        raise RuntimeError(
            "Missing species directory!  Use species_dir "
            + "parameter or set $AIMS_SPECIES_DIR environment variable."
        )

    species_path = Path(species_path)

    species = set(atoms.symbols)

    tier = parameters.pop("tier", None)

    if tier is not None:
        if isinstance(tier, int):
            tierlist = np.ones(len(species), "int") * tier
        elif isinstance(tier, list):
            assert len(tier) == len(species)
            tierlist = tier

    for i, symbol in enumerate(species):
        path = species_path / ("%02i_%s_default" % (atomic_numbers[symbol], symbol))
        reached_tiers = False
        with open(path) as species_fd:
            for line in species_fd:
                if tier is not None:
                    if "First tier" in line:
                        reached_tiers = True
                        targettier = tierlist[i]
                        foundtarget = False
                        do_uncomment = True
                    if reached_tiers:
                        line, foundtarget, do_uncomment = format_tiers(
                            line, targettier, foundtarget, do_uncomment
                        )
                fd.write(line)

        if tier is not None and not foundtarget:
            raise RuntimeError(
                "Basis tier %i not found for element %s" % (targettier, symbol)
            )

        if parameters.get("plus_u") is not None:
            if symbol in parameters.plus_u:
                fd.write("plus_u %s \n" % parameters.plus_u[symbol])


def format_tiers(line, targettier, foundtarget, do_uncomment):
    if "meV" in line:
        assert line[0] == "#"
        if "tier" in line and "Further" not in line:
            tier = line.split(" tier")[0]
            tier = tier.split('"')[-1]
            current_tier = translate_tier(tier)
            if current_tier == targettier:
                foundtarget = True
            elif current_tier > targettier:
                do_uncomment = False
        else:
            do_uncomment = False
        outputline = line
    elif do_uncomment and line[0] == "#":
        outputline = line[1:]
    elif not do_uncomment and line[0] != "#":
        outputline = "#" + line
    else:
        outputline = line
    return outputline, foundtarget, do_uncomment


# Read aims.out files
scalar_property_to_line_key = {
    "free_energy": ["| Electronic free energy"],
    "number_of_iterations": ["| Number of self-consistency cycles"],
    "magnetic_moment": ["N_up - N_down"],
    "n_atoms": ["| Number of atoms"],
    "n_bands": [
        "Number of Kohn-Sham states",
        "Reducing total number of  Kohn-Sham states",
        "Reducing total number of Kohn-Sham states",
    ],
    "n_electrons": ["The structure contains"],
    "n_kpts": ["| Number of k-points"],
    "n_spins": ["| Number of spin channels"],
    "electronic_temp": ["Occupation type:"],
    "fermi_energy": ["| Chemical potential (Fermi level)"],
}


class AimsOutChunk:
    """Base calss for AimsOutChunks"""

    def __init__(self, lines):
        """Constructor

        Parameters
        ----------
        lines: list of str
            The lines for the chunk
        """
        self.lines = lines

    def reverse_search_for(self, keys, line_start=0):
        """Find the last time one of the keys appears in self.lines

        Parameters
        ----------
        keys: list of str
            The key strings to search for in self.lines
        line_start: int
            The lowest index to search for in self.lines

        Returns
        -------
        int
            The last time one of the keys appears in self.lines
        """
        for ll, line in enumerate(self.lines[line_start:][::-1]):
            if any([key in line for key in keys]):
                return len(self.lines) - ll - 1

        return len(self.lines)

    def search_for_all(self, key, line_start=0, line_end=-1):
        """Find the all times the key appears in self.lines

        Parameters
        ----------
        key: str
            The key string to search for in self.lines
        line_start: int
            The first line to start the search from
        line_end: int
            The last line to end the search at

        Returns
        -------
        list of ints
            All times the key appears in the lines
        """
        line_index = []
        for ll, line in enumerate(self.lines[line_start:line_end]):
            if key in line:
                line_index.append(ll + line_start)
        return line_index

    def parse_scalar(self, property):
        """Parse a scalar property from the chunk

        Parameters
        ----------
        property: str
            The property key to parse

        Returns
        -------
        float
            The scalar value of the property
        """
        line_start = self.reverse_search_for(scalar_property_to_line_key[property])

        if line_start >= len(self.lines):
            return None

        line = self.lines[line_start]
        return float(line.split(":")[-1].strip().split()[0])


class AimsOutHeaderChunk(AimsOutChunk):
    """The header of the aims.out file containint general information"""

    def __init__(self, lines):
        """Constructor

        Parameters
        ----------
        lines: list of str
            The lines inside the aims.out header
        """
        super().__init__(lines)
        self._k_points = None
        self._k_point_weights = None

    def _parse_constraints(self):
        """Parse the constraints from the aims.out file

        Constraints for the lattice vectors are not supported.
        """

        line_inds = self.search_for_all("Found relaxation constraint for atom")
        if len(line_inds) == 0:
            return None

        fix = []
        fix_cart = []
        for ll in line_inds:
            line = self.lines[ll]
            xyz = [0, 0, 0]
            ind = int(line.split()[5][:-1]) - 1
            if "All coordinates fixed" in line:
                if ind not in fix:
                    fix.append(ind)
            if "coordinate fixed" in line:
                coord = line.split()[6]
                if coord == "x":
                    xyz[0] = 1
                elif coord == "y":
                    xyz[1] = 1
                elif coord == "z":
                    xyz[2] = 1
                keep = True
                for n, c in enumerate(fix_cart):
                    if ind == c.index:
                        keep = False
                        break
                if keep:
                    fix_cart.append(FixCartesian(ind, xyz))
                else:
                    fix_cart[n].mask[xyz.index(1)] = 0
        if len(fix) > 0:
            fix_cart.append(FixAtoms(indices=fix))

        return fix_cart

    def _parse_initial_cell(self):
        """Parse the initial cell from the aims.out file"""
        line_start = self.reverse_search_for(["| Unit cell:"])
        if line_start < len(self.lines):
            cell = [
                [float(inp) for inp in line.split()[-3:]]
                for line in self.lines[line_start + 1:line_start + 4]
            ]
        else:
            cell = None
        return cell

    def _parse_initial_atoms(self):
        """Create an atoms object for the initial geometry.in structure from the aims.out file"""
        line_start = self.reverse_search_for(["Atomic structure:"]) + 2
        if line_start >= len(self.lines):
            raise IOError("No structure information is inside the chunk.")

        cell = self.initial_cell
        atoms = Atoms()
        for line in self.lines[line_start:line_start + self.n_atoms]:
            inp = line.split()
            atoms.append(Atom(inp[3], (float(inp[4]), float(inp[5]), float(inp[6]))))
        assert len(atoms) == self.n_atoms

        if cell:
            atoms.set_cell(cell)
            atoms.set_pbc([True, True, True])
        atoms.set_constraint(self.constraints)

        return atoms

    def _parse_md(self):
        """Determine if the calculation is a molecular dynamics calculation or not"""
        return len(self.lines) > self.reverse_search_for(
            ["Complete information for previous time-step:"]
        )

    def _parse_relaxation(self):
        """Determine if the calculation is a geometry optimization or not"""
        return len(self.lines) > self.reverse_search_for(["Geometry relaxation:"])

    def _parse_k_points(self):
        """Get the list of k-points used in the calculation"""
        n_kpts = self.parse_scalar("n_kpts")
        if n_kpts is None:
            return
        n_kpts = int(n_kpts)

        line_start = self.reverse_search_for(["| K-points in task"]) + 1
        line_end = self.reverse_search_for(["| k-point:"]) + 1
        if line_end - line_start != n_kpts:
            return

        self._k_points = np.zeros((n_kpts, 3))
        self._k_point_weights = np.zeros((n_kpts))
        for kk, line in enumerate(self.lines[line_start:line_end]):
            self._k_points[kk] = [float(inp) for inp in line.split()[4:7]]
            self._k_point_weights[kk] = float(line.split()[-1])

    @lazyproperty
    def initial_atoms(self):
        """The initial structure defined in the geoemtry.in file"""
        return self._parse_initial_atoms()

    @lazyproperty
    def initial_cell(self):
        """The initial lattice vectors defined in the geoemtry.in file"""
        return self._parse_initial_cell()

    @lazyproperty
    def constraints(self):
        """The relaxation constraints for the calculation"""
        return self._parse_constraints()

    @lazyproperty
    def is_md(self):
        """True if the calculation is a molecular dynamics calculation"""
        return self._parse_md()

    @lazyproperty
    def is_relaxation(self):
        """True if the calculation is a geometry optimization"""
        return self._parse_relaxation()

    @lazyproperty
    def n_atoms(self):
        """The number of atoms for the material"""
        n_atoms = self.parse_scalar("n_atoms")
        if n_atoms is None:
            raise IOError("No information about the number of atoms in the header.")
        return int(n_atoms)

    @lazyproperty
    def n_bands(self):
        """The number of Kohn-Sham states for the chunk"""
        line_start = self.reverse_search_for(scalar_property_to_line_key["n_bands"])

        if line_start >= len(self.lines):
            raise IOError(
                "No information about the number of Kohn-Sham states in the header."
            )

        line = self.lines[line_start]
        if "| Number of Kohn-Sham states" in line:
            return int(line.split(":")[-1].strip().split()[0])

        return int(line.split()[-1].strip()[:-1])

    @lazyproperty
    def n_electrons(self):
        """The number of electrons for the chunk"""
        line_start = self.reverse_search_for(scalar_property_to_line_key["n_electrons"])

        if line_start >= len(self.lines):
            raise IOError("No information about the number of electrons in the header.")

        line = self.lines[line_start]
        return int(float(line.split()[-2]))

    @lazyproperty
    def n_k_points(self):
        """The number of k_ppoints for the calculation"""
        n_kpts = self.parse_scalar("n_kpts")
        if n_kpts is None:
            return None

        return int(n_kpts)

    @lazyproperty
    def n_spins(self):
        """The number of spin channels for the chunk"""
        n_spins = self.parse_scalar("n_spins")
        if n_spins is None:
            raise IOError(
                "No information about the number of spin channels in the header."
            )
        return int(n_spins)

    @lazyproperty
    def electronic_temperature(self):
        """The electronic temperature for the chunk"""
        line_start = self.reverse_search_for(
            scalar_property_to_line_key["electronic_temp"]
        )
        if line_start >= len(self.lines):
            return 0.10

        line = self.lines[line_start]
        return float(line.split("=")[-1].strip().split()[0])

    @property
    def k_points(self):
        """All k-points listed in the calculation"""
        if self._k_points is None:
            self._parse_k_points()
        return self._k_points

    @property
    def k_point_weights(self):
        """The k-point weights for the calculation"""
        if self._k_point_weights is None:
            self._parse_k_points()
        return self._k_point_weights

    @lazyproperty
    def header_summary(self):
        """Dictionary summarizing the information inside the header"""
        return {
            "initial_atoms": self.initial_atoms,
            "initial_cell": self.initial_cell,
            "constraints": self.constraints,
            "is_relaxation": self.is_relaxation,
            "is_md": self.is_md,
            "n_atoms": self.n_atoms,
            "n_bands": self.n_bands,
            "n_electrons": self.n_electrons,
            "n_spins": self.n_spins,
            "electronic_temperature": self.electronic_temperature,
            "n_k_points": self.n_k_points,
            "k_points": self.k_points,
            "k_point_weights": self.k_point_weights,
        }


class AimsOutCalcChunk(AimsOutChunk):
    """A part of the aims.out file correponding to a single structure"""

    def __init__(self, lines, header):
        """Constructor

        Parameters
        ----------
        lines: list of str
            The lines used for the structure
        header: dict
            A summary of the relevant information from the aims.out header
        """
        super().__init__(lines)
        self._header = header.header_summary
        self._atoms = None
        self._hirshfeld_charges = None
        self._hirshfeld_atomic_dipoles = None
        self._hirshfeld_volumes = None
        self._hirshfeld_dipole = None
        self._eigenvalues = None
        self._occupancies = None

    def _parse_atoms(self):
        """Create an atoms object for the subsequent structures calculated in the aims.out file"""
        start_keys = [
            "Atomic structure (and velocities) as used in the preceding time step",
            "Updated atomic structure",
            "Atomic structure that was used in the preceding time step of the wrapper",
        ]
        line_start = self.reverse_search_for(start_keys)
        if line_start >= len(self.lines):
            self._atoms = self.initial_atoms
            return self._atoms
        line_end = self.reverse_search_for(
            ['Writing the current geometry to file "geometry.in.next_step"'], line_start
        )
        cell = []
        velocities = []
        atoms = Atoms()
        for line in self.lines[line_start + 1:line_end]:
            if "lattice_vector   " in line:
                cell.append([float(inp) for inp in line.split()[1:]])
            elif "atom   " in line:
                line_split = line.split()
                atoms.append(
                    Atom(line_split[4], tuple([float(inp) for inp in line_split[1:4]]))
                )
            elif "velocity   " in line:
                velocities.append([float(inp) for inp in line.split()[1:]])

        assert len(atoms) == self.n_atoms
        assert (len(velocities) == self.n_atoms) or (len(velocities) == 0)
        if len(cell) == 3:
            atoms.set_cell(np.array(cell))
            atoms.set_pbc([True, True, True])
        elif len(cell) != 0:
            raise IOError("Parsed geometry has incorrect number of lattice vectors.")

        if len(velocities) > 0:
            atoms.set_velocities(np.array(velocities))
        atoms.set_constraint(self.constraints)

        self._atoms = atoms
        return atoms

    def _parse_forces(self):
        """Parse the forces from the aims.out file"""
        line_start = self.reverse_search_for(["Total atomic forces"]) + 1
        if line_start >= len(self.lines):
            return

        return np.array(
            [
                [float(inp) for inp in line.split()[-3:]]
                for line in self.lines[line_start:line_start + self.n_atoms]
            ]
        )

    def _parse_stresses(self):
        """Parse the stresses from the aims.out file"""
        line_start = 3 + self.reverse_search_for(
            ["Per atom stress (eV) used for heat flux calculation"]
        )
        if line_start >= len(self.lines):
            return
        stresses = []
        for line in self.lines[line_start:line_start + self.n_atoms]:
            xx, yy, zz, xy, xz, yz = [float(d) for d in line.split()[2:8]]
            stresses.append([xx, yy, zz, yz, xz, xy])

        return np.array(stresses)

    def _parse_stress(self):
        """Parse the stress from the aims.out file"""
        from ase.stress import full_3x3_to_voigt_6_stress

        line_start = self.reverse_search_for(
            ["Analytical stress tensor - Symmetrized"]
        )  # Offest to relevant lines
        if line_start >= len(self.lines):
            return

        stress = [
            [float(inp) for inp in line.split()[2:5]]
            for line in self.lines[line_start + 5:line_start + 8]
        ]
        return full_3x3_to_voigt_6_stress(stress)

    def _parse_metallic(self):
        """Checks the outputfile to see if the chunk corresponds to a metallic system"""
        line_start = self.reverse_search_for(
            [
                "material is metallic within the approximate finite broadening function (occupation_type)"
            ]
        )
        return line_start != len(self.lines)

    def _parse_energy(self):
        """Parse the energy from the aims.out file"""
        if self._atoms is None:
            self._parse_atoms()

        if np.all(self._atoms.pbc) and self.is_metallic:
            line_ind = self.reverse_search_for(["Total energy corrected"])
        else:
            line_ind = self.reverse_search_for(["Total energy uncorrected"])
        if line_ind == len(self.lines):
            raise IOError("No energy is associated with the structure.")

        return float(self.lines[line_ind].split()[5])

    def _parse_dipole(self):
        """Parse the electric dipole moment from the aims.out file file."""
        line_start = self.reverse_search_for(["Total dipole moment [eAng]"])
        if line_start >= len(self.lines):
            return

        line = self.lines[line_start]
        return np.array([float(inp) for inp in line.split()[6:9]])

    def _parse_hirshfeld(self):
        """Parse the Hirshfled charges volumes, and dipole moments from the ouput"""
        if self._atoms is None:
            self._parse_atoms()

        line_start = self.reverse_search_for(
            ["Performing Hirshfeld analysis of fragment charges and moments."]
        )
        if line_start >= len(self.lines):
            return None

        line_inds = self.search_for_all("Hirshfeld charge", line_start, -1)
        self._hirshfeld_charges = np.array(
            [float(self.lines[ind].split(":")[1]) for ind in line_inds]
        )

        line_inds = self.search_for_all("Hirshfeld volume", line_start, -1)
        self._hirshfeld_volumes = np.array(
            [float(self.lines[ind].split(":")[1]) for ind in line_inds]
        )

        line_inds = self.search_for_all("Hirshfeld dipole vector", line_start, -1)
        self._hirshfeld_atomic_dipoles = np.array(
            [
                [float(inp) for inp in self.lines[ind].split(":")[1].split()]
                for ind in line_inds
            ]
        )

        if not np.any(self._atoms.pbc):
            self._hirshfeld_dipole = np.sum(
                self._hirshfeld_charges.reshape((-1, 1)) * self._atoms.get_positions(),
                axis=1,
            )

    def _parse_eigenvalues(self):
        """Parse the eigenvalues and occupancies of the system. If eigenvalue for a particular k-point is not present in the output file then set it to np.nan"""
        if self._atoms is None:
            self._parse_atoms()

        # line_start = self.reverse_search_for(
        #     [
        #         "Convergence:    q app. |  density  | eigen (eV) | Etot (eV)",
        #         "Begin self-consistency iteration #",
        #     ]
        # )
        line_start = self.reverse_search_for(
            ["Writing Kohn-Sham eigenvalues."]
        )
        if line_start >= len(self.lines):
            return
        line_end = min(
            self.reverse_search_for(["Self-consistency cycle converged."], line_start),
            self.reverse_search_for(
                [
                    "What follows are estimated values for band gap, HOMO, LUMO, etc.",
                    "Current spin moment of the entire structure :",
                ],
                line_start,
            ),
        )
        n_kpts = self.n_k_points if np.all(self._atoms.pbc) else 1
        if n_kpts is None:
            return

        self._eigenvalues = np.full((n_kpts, self.n_bands, self.n_spins), np.nan)
        self._occupancies = np.full((n_kpts, self.n_bands, self.n_spins), np.nan)

        occupation_block_start = self.search_for_all(
            "State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]",
            line_start,
            line_end,
        )
        kpt_def = self.search_for_all("K-point: ", line_start, line_end)

        if self.n_k_points:
            kpt_inds = [int(self.lines[ll].split()[1]) - 1 for ll in kpt_def]
        else:
            kpt_inds = [0]

        assert len(kpt_inds) == len(occupation_block_start)
        spins = [0] * len(occupation_block_start)

        if self.n_spins == 2:
            spin_def = self.search_for_all("Spin-", line_start, line_end)
            assert len(spin_def) == len(occupation_block_start)

            spins = [int("Spin-down eigenvalues:" in self.lines[ll]) for ll in spin_def]

        for occ_start, kpt_ind, spin in zip(occupation_block_start, kpt_inds, spins):
            for ll, line in enumerate(
                self.lines[occ_start + 1:occ_start + self.n_bands + 1]
            ):
                line = line.replace("**************", "         10000")
                line = line.replace("***************", "          10000")
                line = line.replace("****************", "           10000")
                self._eigenvalues[kpt_ind, ll, spin] = float(line.split()[3])
                self._occupancies[kpt_ind, ll, spin] = float(line.split()[1])

    @property
    def atoms(self):
        """Convert AimsOutChunk to Atoms object and add all non-standard outputs to atoms.info"""
        if self._atoms is None:
            self._atoms = self._parse_atoms()

        self._atoms.calc = SinglePointDFTCalculator(
            self._atoms,
            energy=self.energy,
            free_energy=self.free_energy,
            forces=self.forces,
            stress=self.stress,
            stresses=self.stresses,
            magmom=self.magmom,
            dipole=self.dipole,
        )

        self._atoms.info["fermi_energy"] = self.E_f
        self._atoms.info["n_iter"] = self.n_iter
        self._atoms.info["hirshfeld_charges"] = self.hirshfeld_charges
        self._atoms.info["hirshfeld_dipole"] = self.hirshfeld_dipole
        self._atoms.info["hirshfeld_volumes"] = self.hirshfeld_volumes
        self._atoms.info["hirshfeld_atomic_dipoles"] = self.hirshfeld_atomic_dipoles
        self._atoms.info["eigenvalues"] = self.eigenvalues
        self._atoms.info["occupancies"] = self.occupancies

        return self._atoms

    @property
    def results(self):
        """Convert an AimsOutChunk to a Results Dictionary"""
        results = {
            "energy": self.energy,
            "free_energy": self.free_energy,
            "forces": self.forces,
            "stress": self.stress,
            "stresses": self.stresses,
            "magmom": self.magmom,
            "dipole": self.dipole,
            "fermi_energy": self.E_f,
            "n_iter": self.n_iter,
            "hirshfeld_charges": self.hirshfeld_charges,
            "hirshfeld_dipole": self.hirshfeld_dipole,
            "hirshfeld_volumes": self.hirshfeld_volumes,
            "hirshfeld_atomic_dipoles": self.hirshfeld_atomic_dipoles,
            "eigenvalues": self.eigenvalues,
            "occupancies": self.occupancies,
        }

        del_keys = []
        for key, val in results.items():
            if val is None:
                del_keys.append(key)

        for key in del_keys:
            del results[key]

        return results

    # Properties from the aims.out header
    @lazyproperty
    def initial_atoms(self):
        """The initial structure defined in the geoemtry.in file"""
        return self._header["initial_atoms"]

    @lazyproperty
    def initial_cell(self):
        """The initial lattice vectors defined in the geoemtry.in file"""
        return self._header["initial_cell"]

    @lazyproperty
    def constraints(self):
        """The relaxation constraints for the calculation"""
        return self._header["constraints"]

    @lazyproperty
    def n_atoms(self):
        """The number of atoms for the material"""
        return self._header["n_atoms"]

    @lazyproperty
    def n_bands(self):
        """The number of Kohn-Sham states for the chunk"""
        return self._header["n_bands"]

    @lazyproperty
    def n_electrons(self):
        """The number of electrons for the chunk"""
        return self._header["n_electrons"]

    @lazyproperty
    def n_spins(self):
        """The number of spin channels for the chunk"""
        return self._header["n_spins"]

    @lazyproperty
    def electronic_temperature(self):
        """The electronic temperature for the chunk"""
        return self._header["electronic_temperature"]

    @lazyproperty
    def n_k_points(self):
        """The number of electrons for the chunk"""
        return self._header["n_k_points"]

    @lazyproperty
    def k_points(self):
        """The number of spin channels for the chunk"""
        return self._header["k_points"]

    @lazyproperty
    def k_point_weights(self):
        """k_point_weights electronic temperature for the chunk"""
        return self._header["k_point_weights"]

    @lazyproperty
    def forces(self):
        """The forces for the chunk"""
        return self._parse_forces()

    @lazyproperty
    def stresses(self):
        """The atomic stresses for the chunk"""
        return self._parse_stresses()

    @lazyproperty
    def stress(self):
        """The stress for the chunk"""
        return self._parse_stress()

    # Properties for the specfic structure
    @lazyproperty
    def energy(self):
        """The energy for the chunk"""
        return self._parse_energy()

    @lazyproperty
    def free_energy(self):
        """The free energy for the chunk"""
        return self.parse_scalar("free_energy")

    @lazyproperty
    def n_iter(self):
        """The number of SCF iterations needed to converge the SCF cycle for the chunk"""
        return self.parse_scalar("number_of_iterations")

    @lazyproperty
    def magmom(self):
        """The magnetic moment for the chunk"""
        return self.parse_scalar("magnetic_moment")

    @lazyproperty
    def E_f(self):
        """The Fermi energy for the chunk"""
        return self.parse_scalar("fermi_energy")

    @lazyproperty
    def dipole(self):
        """The electronic dipole moment for the chunk"""
        return self._parse_dipole()

    @lazyproperty
    def is_metallic(self):
        """True if the chunk is for a metallic material"""
        return self._parse_metallic()

    @lazyproperty
    def converged(self):
        """True if the chunk is a fully converged final structure"""
        return (len(self.lines) > 0) and ("Have a nice day." in self.lines[-5:])

    @property
    def hirshfeld_charges(self):
        """The Hirshfeld charges for the chunk"""
        if self._hirshfeld_charges is None:
            self._parse_hirshfeld()
        return self._hirshfeld_charges

    @property
    def hirshfeld_atomic_dipoles(self):
        """The Hirshfeld atomic dipole moments for the chunk"""
        if self._hirshfeld_atomic_dipoles is None:
            self._parse_hirshfeld()
        return self._hirshfeld_atomic_dipoles

    @property
    def hirshfeld_volumes(self):
        """The Hirshfeld volume for the chunk"""
        if self._hirshfeld_volumes is None:
            self._parse_hirshfeld()
        return self._hirshfeld_volumes

    @property
    def hirshfeld_dipole(self):
        """The Hirshfeld systematic dipole moment for the chunk"""
        if self._atoms is None:
            self._parse_atoms()

        if not np.any(self._atoms.pbc) and self._hirshfeld_dipole is None:
            self._parse_hirshfeld()
        return self._hirshfeld_dipole

    @property
    def eigenvalues(self):
        """All outputted eigenvalues for the system"""
        if self._eigenvalues is None:
            self._parse_eigenvalues()
        return self._eigenvalues

    @property
    def occupancies(self):
        """All outputted occupancies for the system"""
        if self._occupancies is None:
            self._parse_eigenvalues()
        return self._occupancies


def get_header_chunk(fd):
    """Returns the header information from the aims.out file"""
    header = []
    line = ""

    # Stop the header once the first SCF cycle begins
    while (
        "Convergence:    q app. |  density  | eigen (eV) | Etot (eV)" not in line
        and "Begin self-consistency iteration #" not in line
    ):
        try:
            line = next(fd).strip()  # Raises StopIteration on empty file
        except StopIteration:
            return
        header.append(line)
    return AimsOutHeaderChunk(header)


def get_aims_out_chunks(fd, header_chunk):
    """Yield unprocessed chunks (header, lines) for each AimsOutChunk image."""
    try:
        line = next(fd).strip()  # Raises StopIteration on empty file
    except StopIteration:
        return

    # If the calculation is relaxation the updated structural information occurs before the re-initialization
    if header_chunk.is_relaxation:
        chunk_end_line = (
            "Geometry optimization: Attempting to predict improved coordinates."
        )
    else:
        chunk_end_line = "Begin self-consistency loop: Re-initialization"

    # If SCF is not converged then do not treat the next chunk_end_line as a new chunk until after the SCF is re-initialized
    ignore_chunk_end_line = False
    while True:
        try:
            line = next(fd).strip()  # Raises StopIteration on empty file
        except StopIteration:
            break

        lines = []
        while chunk_end_line not in line or ignore_chunk_end_line:
            lines.append(line)
            # If SCF cycle not converged, don't end chunk on next Re-initialization
            if (
                "Self-consistency cycle not yet converged - restarting mixer to attempt better convergence."
                in line
            ):
                ignore_chunk_end_line = True
            elif "Begin self-consistency loop: Re-initialization" in line:
                ignore_chunk_end_line = False

            try:
                line = next(fd).strip()
            except StopIteration:
                break

        yield AimsOutCalcChunk(lines, header_chunk)


@reader
def read_aims_output(fd, index=-1):
    """Import FHI-aims output files with all data available, i.e.
    relaxations, MD information, force information etc etc etc."""
    header_chunk = get_header_chunk(fd)
    chunks = list(get_aims_out_chunks(fd, header_chunk))
    if not chunks[-1].converged:
        warnings.warn(
            "This FHI-aims calculation did not fully converge, check the results and output file carefully."
        )

    # Relaxations have an additional fotter chunk due to how it is split
    if header_chunk.is_relaxation:
        images = [chunk.atoms for chunk in chunks[:-1]]
    else:
        images = [chunk.atoms for chunk in chunks]
    return images[index]


@reader
def read_aims_results(fd, index=-1):
    """Import FHI-aims output files and summarize all relevant information into a dictionary"""
    header_chunk = get_header_chunk(fd)
    chunks = list(get_aims_out_chunks(fd, header_chunk))
    if not chunks[-1].converged:
        warnings.warn(
            "This FHI-aims calculation did not fully converge, check the results and output file carefully."
        )

    # Relaxations have an additional fotter chunk due to how it is split
    if header_chunk.is_relaxation and (index == -1):
        return chunks[-2].results

    return chunks[index].results
