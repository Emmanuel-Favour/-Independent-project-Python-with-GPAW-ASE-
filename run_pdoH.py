from ase.io import read, write
from ase.build import molecule, add_adsorbate
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from gpaw import GPAW, PW, FermiDirac
from gpaw import GPAW, LCAO, MixerDif
from ase import Atoms

# Load your optimized Pd slab
slab = read('pd_slab_opt.gpw')

# --- This block directly fixes the TypeError ---
# It replaces the 'position="fcc"' with manual coordinates.
# For a 3x3 supercell of an fcc(111) surface, the high-symmetry fcc site
# can be found at the center of the supercell's xy-plane.
cell_vectors = slab.get_cell()
x_center = cell_vectors[0, 0] / 2.0 + cell_vectors[1, 0] / 2.0
y_center = cell_vectors[0, 1] / 2.0 + cell_vectors[1, 1] / 2.0
fcc_site_coords = (x_center, y_center)
# --- End of the TypeError fix ---

# Place OH with O pointing towards the surface
# The `position` argument now takes the manually specified (x, y) coordinates.
add_adsorbate(slab, Atoms('OH', [(0, 0, 0.0), (0, 0, -1.0)]), height=2.0, position=fcc_site_coords)

# Fix the bottom layers, as in the slab optimization
# The tags were set by the fcc111 builder and persisted in the gpw file.
c = FixAtoms(indices=[atom.index for atom in slab if atom.tag > 2])
slab.set_constraint(c)
# Set up the GPAW calculator, matching the efficient LCAO settings used for the slab.
# Use the slab calculation as a restart to speed up convergence.
calc = GPAW(
    mode=LCAO(),
    basis={'Pd': 'sz(dzp)', 'O': 'szp(dzp)', 'H': 'szp(dzp)'}, # Stick to these if possible
    # basis={'Pd': 'sz', 'O': 'sz', 'H': 'sz'}, # Use for fastest, least accurate runs
    xc='PBE',
    kpts={'density': 1.5, 'gamma': True}, # Reduced k-point density for speed
    spinpol=True,
    occupations={'name': 'fermi-dirac', 'width': 0.05}, # Slightly smaller smearing
    txt='Pd_OH.txt',
    symmetry={'point_group': False},
    restart='pd_slab_opt.gpw',
    mixer=MixerDif(beta=0.1, nmaxold=5, weight=50.0) # Robust mixer for stability
)
slab.calc = calc

# Relax structure
dyn = BFGS(slab, trajectory='Pd_OH.traj', logfile='Pd_OH.log')
dyn.run(fmax=0.01)

# Save results
calc.write('Pd_OH.gpw', mode='all')
write('Pd_OH.traj', slab)
print("Relaxation finished. Pd+OH structure saved.")

