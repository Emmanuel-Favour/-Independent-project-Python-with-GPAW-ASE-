from ase.io import read, write
from ase.build import add_adsorbate
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from gpaw import GPAW, LCAO, MixerDif
from ase import Atoms

# Read in the optimized Pd slab from the previous calculation
slab = read('pd_slab_opt.gpw')

# Place OH further from surface to reduce initial strong repulsion
cell_vectors = slab.get_cell()
x_center = cell_vectors[0, 0] / 2.0 + cell_vectors[1, 0] / 2.0
y_center = cell_vectors[0, 1] / 2.0 + cell_vectors[1, 1] / 2.0
fcc_site_coords = (x_center, y_center)

# Place OH higher (3.0 Å instead of 2.0)
add_adsorbate(slab, Atoms('OH', [(0, 0, 0.0), (0, 0, -1.0)]),
              height=3.0, position=fcc_site_coords)

# Fix bottom layers
c = FixAtoms(indices=[atom.index for atom in slab if atom.tag > 2])
slab.set_constraint(c)

# GPAW calculator with simpler basis + looser settings for speed/stability
calc = GPAW(
    mode=LCAO(),
   basis={'Pd': 'sz(dzp)', 'O': 'szp(dzp)', 'H': 'szp(dzp)'},  # bigger basis = smoother forces
    xc='PBE',
    kpts={'density': 2.0, 'gamma': True},         # slightly denser k-mesh
    spinpol=True,
    occupations={'name': 'fermi-dirac', 'width': 0.10},  # more smearing helps convergence
    txt='Pd_OH.txt',
    symmetry={'point_group': False},
    restart='pd_slab_opt.gpw',
    mixer=MixerDif(beta=0.05, nmaxold=10, weight=100.0)  # gentler, more stable mixing
)
slab.calc = calc

# Relax structure with looser convergence
dyn = BFGS(slab, trajectory='Pd_OH.traj', logfile='Pd_OH.log')
dyn.run(fmax=0.1)  # relaxed criterion

# Save results
calc.write('Pd_OH.gpw', mode='all')
write('Pd_OH.traj', slab)
print("Relaxation finished. Pd+OH structure saved.")