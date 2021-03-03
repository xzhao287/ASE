import pytest
import numpy as np

from ase.build import bulk, molecule

import sys

from ase.optimize import BFGS
from ase.calculators.espresso import Espresso
from ase.calculators.socketio import SocketIOCalculator
from ase.constraints import ExpCellFilter
from ase.units import Ry


abinit_boilerplate = dict(
    ionmov=28,
    expert_user=1,
    optcell=2,
    tolmxf=1e-300,
    ntime=100_000,
    ecutsm=0.5,
)


commands = dict(
    espresso='{exe} < PREFIX.pwi --ipi {unixsocket}:UNIX > PREFIX.pwo',
    abinit='{exe} PREFIX.in --ipi {unixsocket}:UNIX > PREFIX.log',
)


calc = pytest.mark.calculator
@calc('espresso', ecutwfc=200 / Ry)
@calc('abinit', ecut=200, **abinit_boilerplate)
def test_socketio_espresso(factory):
    atoms = bulk('Si')

    name = factory.name
    exe = factory.factory.executable  # XXX ugly
    unixsocket = f'ase_test_socketio_{name}'

    espresso = factory.calc(
        kpts=[2, 2, 2],
    )
    template = commands[name]
    command = template.format(exe=exe, unixsocket=unixsocket)
    espresso.command = command

    atoms.rattle(stdev=.2, seed=42)

    opt = BFGS(ExpCellFilter(atoms))

    with pytest.warns(UserWarning, match='Subprocess exited'):
        with SocketIOCalculator(espresso, unixsocket=unixsocket) as calc:
            atoms.calc = calc
            for _ in opt.irun(fmax=0.05):
                e = atoms.get_potential_energy()
                fmax = max(np.linalg.norm(atoms.get_forces(), axis=0))
                print(e, fmax)
