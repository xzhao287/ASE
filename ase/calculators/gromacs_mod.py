"""This module defines an ASE interface to GROMACS.

http://www.gromacs.org/
It is VERY SLOW compared to standard Gromacs
(due to slow formatted io required here).

Mainly intended to be the MM part in the ase QM/MM

Markus.Kaukonen@iki.fi

To be done:
1) change the documentation for the new file-io-calculator (test works now)
2) change gromacs program names
-now:     hard coded
-future:  set as dictionary in params_runs

"""

import os, datetime
import subprocess
from glob import glob

import numpy as np
from ase import Atoms
from ase import units
from ase.calculators.calculator import (CalculatorSetupError, FileIOCalculator,
                                        all_changes)
from ase.calculators.calculator import Calculator
from ase.io.gromos import read_gromos, write_gromos
from ase.data import chemical_symbols
from ase import units
import pdb

def parse_gromacs_version(output):
    import re
    match = re.search(r'GROMACS version\:\s*(\S+)', output, re.M)
    return match.group(1)


def get_gromacs_version(executable):
    output = subprocess.check_output([executable, '--version'],
                                     encoding='utf-8')
    return parse_gromacs_version(output)


def do_clean(name='#*'):
    """ remove files matching wildcards """
    myfiles = glob(name)
    for myfile in myfiles:
        try:
            os.remove(myfile)
        except OSError:
            pass

def read_gro(fileobj=None):
    """Read gromos geometry files (.gro).
    Reads:
    atom positions,
    and simulation cell (if present)
    tries to set atom types
    """
    # add a sanity check see if the gro file exists
    linenum = 0
    symbols = []
    coords = []
    velocities = []
    atomnum = 0
    if fileobj == None:
        pass # add sanity check
    with open(fileobj, 'r') as gro:
        for line in gro:
            linenum += 1
            if linenum == 2:
                atomnum = int(line)
                print(f'The number of atoms in gro is {atomnum}')
            elif linenum > 2 and len(line)>40:
                # need to check the units later
                x = 10*float(line[20:28])
                y = 10*float(line[28:36])
                z = 10*float(line[36:44])
                symbol = str(line[8:16]).strip()[0]
                coords.append([x, y, z])
                velocities.append([0, 0, 0])
                symbols.append(symbol)
            elif linenum > atomnum + 2:
                pass # use this case to set cell
    gro.close()
    gmx_system = Atoms(symbols=symbols,
                    positions=coords,
                    velocities=velocities)
    return gmx_system

def write_g96(g96_name, gro_name, atoms):
    gro = open(gro_name, 'r')
    g96 = open(g96_name, 'w')
    header = "TITLE\n"+\
                "t=   0.00000 step= 0\n"+\
                "END\n"+\
                "TIMESTEP\n"+\
                "              0       0.000000\n"+\
                "END\n"+\
                "POSITION\n" # to mod to write time & step
    g96.write(header)
    linenum = 0
    coords = atoms.get_positions() / 10.0
    for line in gro:
        linenum += 1
        atomnum = linenum-3
        if linenum > 2 and len(line)>40:
            newline = f'{line[:20]}'
            coord_str  = f'  {coords[atomnum][0]:15.9f}'
            coord_str += f'{coords[atomnum][1]:15.9f}'
            coord_str += f'{coords[atomnum][2]:15.9f}\n'
            newline += coord_str
            g96.write(newline)
    g96.write("END\n"+"BOX\n")
    cell = atoms.get_cell()
    # print(np.array(cell[0][0]))
    g96.write(str(f'{cell[0][0]/10:15.9f}{cell[1][1]/10:15.9f}{cell[2][2]/10:15.9f}'))
    g96.write("\nEND\n")


class Gromacs(FileIOCalculator):

    """Class for doing GROMACS calculations.
    Before running a gromacs calculation you must prepare the input files
    separately (pdb2gmx and grompp for instance.)

    Input parameters for gromacs runs (the .mdp file)
    are given in self.params and can be set when initializing the calculator
    or by method set_own.
    for example::

        CALC_MM_RELAX = Gromacs()
        CALC_MM_RELAX.set_own_params('integrator', 'steep',
                                     'use steepest descent')

    Run command line arguments for gromacs related programs:
    pdb2gmx, grompp, mdrun, energy, traj.  These can be given as::

        CALC_MM_RELAX = Gromacs()
        CALC_MM_RELAX.set_own_params_runs('force_field', 'oplsaa')
    """

    implemented_properties = ['energy', 'forces']
    discard_results_on_any_change = True

    default_parameters = dict(
        define='-DFLEXIBLE',
        integrator='cg',
        nsteps='10000',
        nstfout='10',
        nstlog='10',
        nstenergy='10',
        nstlist='10',
        ns_type='grid',
        pbc='xyz',
        rlist='1.15',
        coulombtype='PME-Switch',
        rcoulomb='0.8',
        vdwtype='shift',
        rvdw='0.8',
        rvdw_switch='0.75',
        DispCorr='Ener'
        )

    def __init__(self, restart=None,
                 ignore_bad_restart_file=FileIOCalculator._deprecated,
                 label='gromacs', atoms=None,
                 do_qmmm=False, clean=True,
                 water_model='tip3p', force_field='oplsaa', command=None,
                 ASEs_gmxfolder=None,
                 cname_g96=None,
                 cname_gro=None,
                 pname=None,
                 fname=None,
                 nname=None,
                 **kwargs):
        """Construct GROMACS-calculator object.

        Parameters
        ==========
        label: str
            Prefix to use for filenames (label.in, label.txt, ...).
            Default is 'gromacs'.

        do_qmmm : bool
            Is gromacs used as mm calculator for a qm/mm calculation

        clean :     bool
            Remove gromacs backup files
            and old gormacs.* files

        water_model: str
            Water model to be used in gromacs runs (see gromacs manual)

        force_field: str
            Force field to be used in gromacs runs

        command : str
            Gromacs executable; if None (default), choose available one from
            ('gmx', 'gmx_d', 'gmx_mpi', 'gmx_mpi_d')
        """

        self.do_qmmm = do_qmmm
        self.water_model = water_model
        self.force_field = force_field
        self.clean = clean
        self.params_doc = {}
        # add comments for gromacs input file
        self.params_doc['define'] = \
            'flexible/ rigid water'
        self.params_doc['integrator'] = \
            'md: molecular dynamics(Leapfrog), \n' + \
            '; md-vv: molecular dynamics(Velocity Verlet), \n' + \
            '; steep: steepest descent minimization, \n' + \
            '; cg: conjugate cradient minimization \n'

        self.positions = None
        self.atoms = None

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, command=command,
                                  **kwargs)
        self.set(**kwargs)
        # default values for runtime parameters
        # can be changed by self.set_own_params_runs('key', 'value')
        if ASEs_gmxfolder is None:
            ASEs_gmxfolder = './'
        if cname_g96 is None:
            cname_g96 = 'conf'
        if cname_gro is None:
            cname_gro = 'conf'
        if fname is None:
            fname = 'params'
        if pname is None:
            pname = 'topol'
        if nname is None:
            nname = 'index'

        self.ASEs_gmxfolder = ASEs_gmxfolder
        self.g96name = ASEs_gmxfolder + cname_g96+'.g96'
        self.groname = ASEs_gmxfolder + cname_gro+'.gro'
        self.fname = ASEs_gmxfolder + fname+'.mdp'
        self.pname = ASEs_gmxfolder + pname+'.top'
        self.nname = ASEs_gmxfolder + nname+'.ndx'
        self.sname = ASEs_gmxfolder + label
        # later add sanity check if needed files are present, for QM/MM calcs index files are also needed

        self.params_runs = {}
        self.params_runs['index_filename'] = 'index.ndx'
        self.params_runs['init_structure'] = self.label + '.pdb'
        self.params_runs['water'] = self.water_model
        self.params_runs['force_field'] = self.force_field

        # these below are required by qm/mm
        self.topology_filename = self.label + '.top'

        # clean up gromacs backups
        if self.clean:
            do_clean('gromacs.???')

        # write input files for gromacs program energy
        self.write_energy_files()

        if self.do_qmmm:
            self.parameters['integrator'] = 'md'
            self.parameters['nsteps'] = '0'

    def write_energy_files(self):
        """write input files for gromacs force and energy calculations
        for gromacs program energy"""
        filename = self.ASEs_gmxfolder+'inputGenergy.txt'
        with open(filename, 'w') as output:
            output.write('Potential  \n')
            output.write('   \n')
            output.write('   \n')

        filename = self.ASEs_gmxfolder+'inputGtraj.txt'
        with open(filename, 'w') as output:
            output.write('System  \n')
            output.write('   \n')
            output.write('   \n')

    def _get_name(self):
        return 'Gromacs'

    def _execute_gromacs(self, command):
        """ execute gmx command
        Parameters
        ----------
        command : str
        """
        if self.command:
            subprocess.check_call(self.command + ' ' + command, shell=True)
        else:
            raise CalculatorSetupError('Missing gromacs executable')

    def gen_tpr(self, f_name, c_name, p_name, n_name, s_name):
        """ Generates input file for a gromacs mdrun
        based on structure file and topology file
        resulting file is self.label + '.tpr
        """

        # generate gromacs run input file (gromacs.tpr)
        # print(s_name+'*')
        try:
            for f in glob(s_name+'*'):
                os.remove(f)
        except OSError:
            pass

        subcmd = 'grompp'
        command = ' '.join([
            subcmd,
            '-f', f_name,
            '-c', c_name,
            '-p', p_name,
            '-n', n_name,
            '-o', s_name + '.tpr',
            '-maxwarn', '10',
            self.params_runs.get('extra_grompp_parameters', ''),
            f'> {s_name}.{subcmd}out 2>&1'])
        self._execute_gromacs(command)

    def run_deffnm(self, s_name):
        """ runs a gromacs-mdrun with the
        current atom-configuration """

        # clean up gromacs backups
        # if self.clean:
        #     do_clean('#*')
        log_name = s_name+'.mdout'
        log = open(log_name, 'w')
        log.write('Starting Time: '+str(datetime.datetime.now())+'\n')
        subcmd = 'mdrun -deffnm'
        command = [subcmd]
        command += [
            s_name,
            f'>> {log_name} 2>&1']
        command = ' '.join(command)
        self._execute_gromacs(command)
        log.close()

    def run_with_grompp(self, f_name, c_name, p_name, n_name, s_name):
        self.gen_tpr(f_name,
                     c_name,
                     p_name,
                     n_name,
                     s_name)
        self.run_deffnm(s_name)

    def calculate(self, atoms=None, properties=['energy', 'forces'],
                  system_changes=all_changes):
        #pdb.set_trace()
        Calculator.calculate(self, atoms, properties, system_changes)
        """ runs a gromacs-mdrun and
        gets energy and forces
        rest below is to make gromacs calculator
        compactible with ase-Calculator class

        atoms: Atoms object
            Contains positions, unit-cell, ...
        properties: list of str
            List of what needs to be calculated.  Can be any combination
            of 'energy', 'forces'
        system_changes: list of str
            List of what has changed since last calculation.  Can be
            any combination of these five: 'positions', 'numbers', 'cell',
            'pbc', 'initial_charges' and 'initial_magmoms'.

        """
        print("calculator was called")
        g96_name = self.g96name
        gro_name = self.groname
        # print(f'{g96_name=}', f'{gro_name=}')
        # print(f'{self.fname}')
        write_g96(g96_name=g96_name, gro_name=gro_name, atoms=atoms)

        self.run_with_grompp(self.fname,
                             self.g96name,
                             self.pname,
                             self.nname,
                             self.sname)

        subcmd = 'energy'
        command = ' '.join([
            subcmd,
            '-f', self.ASEs_gmxfolder+self.label + '.edr',
            '-o', self.ASEs_gmxfolder+self.label + '_energy.xvg',
            f'< {self.ASEs_gmxfolder}'+'inputGenergy.txt',
            '-dp '
            f'> {self.ASEs_gmxfolder+self.label}.edrout 2>&1'])
        self._execute_gromacs(command)
        with open(self.ASEs_gmxfolder+self.label + '_energy.xvg') as fd:
            lastline = fd.readlines()[-1]
            energy = float(lastline.split()[1])
            print(' energy * units.kJ / units.mol')
            print(energy, units.kJ, units.mol)
            print(units.nm, units.fs)
        # We go for ASE units !
        self.results['energy'] = energy * units.kJ / units.mol
        # energies are about 100 times bigger in Gromacs units
        # when compared to ase units

        subcmd = 'traj'
        command = ' '.join([
            subcmd,
            '-f', self.ASEs_gmxfolder+self.label + '.trr',
            '-s', self.ASEs_gmxfolder+self.label + '.tpr',
            '-of', self.ASEs_gmxfolder+self.label + '_force.xvg',
            '-fp'
            f'< {self.ASEs_gmxfolder}'+'inputGtraj.txt',
            f'> {self.ASEs_gmxfolder+self.label}_force.xvgout 2>&1'])
        self._execute_gromacs(command)
        with open(self.ASEs_gmxfolder+self.label + '_force.xvg') as fd:
            lastline = fd.readlines()[-1]
            forces = np.array([float(f) for f in lastline.split()[1:]])
        # We go for ASE units !gromacsForce.xvg
        tmp_forces = forces / units.nm * units.kJ / units.mol
        tmp_forces = np.reshape(tmp_forces, (-1, 3))
        self.results['forces'] = tmp_forces



