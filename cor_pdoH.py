from ase.io import read, write
from ase.optimize import BFGS
from gpaw import GPAW, LCAO, MixerDif
from ase.constraints import FixAtoms

# Load last frame of the trajectory (the most recent relaxed structure)
slab = read('Pd_OH.traj@-1')

# Re-apply constraints (must redo because traj doesn’t store them fully)
c = FixAtoms(indices=[atom.index for atom in slab if atom.tag > 2])
slab.set_constraint(c)

# Re-attach calculator
calc = GPAW(
    mode=LCAO(),
    basis={'Pd': 'sz(dzp)', 'O': 'szp(dzp)', 'H': 'szp(dzp)'},
    xc='PBE',
    kpts={'density': 2.0, 'gamma': True},
    spinpol=True,
    occupations={'name': 'fermi-dirac', 'width': 0.10},
    txt='Pd_OH_restart.txt',
    symmetry={'point_group': False},
    mixer=MixerDif(beta=0.05, nmaxold=10, weight=100.0)
)
slab.calc = calc

# Resume optimization
dyn = BFGS(slab, trajectory='Pd_OH_restart.traj', logfile='Pd_OH_restart.log')
dyn.run(fmax=0.1)   # continue until forces meet threshold

# Save again
calc.write('Pd_OH_restart.gpw', mode='all')
write('Pd_OH_restart.traj', slab)

