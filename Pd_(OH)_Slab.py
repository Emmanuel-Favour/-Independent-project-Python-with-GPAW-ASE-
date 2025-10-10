from ase.build import add_adsorbate
from gpaw import GPAW, PW
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
import numpy as np
from ase import Atoms
# Corrected import statement: import GPAW directly from gpaw
from ase.io import read, write
from ase.build import bulk
from gpaw import GPAW, LCAO

# Read in the optimized Pd slab from the previous calculation
# This will load the optimized geometry and converged state from the LCAO run.
slab = read('projects/gpaw_calculations/pd_slab_opt.gpw') 
import os
print("Current working directory:", os.getcwd())
print("Files in this directory:", os.listdir())
import os
print(os.getcwd())
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
    basis={'Pd': 'sz', 'O': 'szp', 'H': 'szp'},
    xc='LDA',
    kpts={'density': 2.5, 'gamma': True},
    spinpol=True, # Spin polarization is needed due to the OH radical
    occupations={'name': 'fermi-dirac', 'width': 0.1},
    txt='Pd_OH_ads_opt_lcao.log',
    symmetry={'point_group': False},
    restart='projects/gpaw_calculations/pd_slab_opt.gpw'
)
slab.calc = calc

# Run the optimization
dyn = QuasiNewton(slab, trajectory='Pd_OH_ads_opt_lcao.traj', logfile='Pd_OH_ads_opt_lcao.log')
dyn.run(fmax=0.2, steps=50)

# Save the final optimized structure and energy
slab.write('Pd_OH_ads_opt_lcao.gpw')
E_test = slab.get_potential_energy()
print("Single-point Pd+OH energy:", E_test)
# E_slab_oh = slab.get_potential_energy()
# print(f"Optimized energy of the Pd+OH system: {E_slab_oh:.4f} eV")

# # Save energy
# with open('energies.txt', 'a') as f:
#     f.write(f"E_slab_oh (LCAO): {E_slab_oh:.4f}\n")

# slab.write('Pd_OH_ads_opt.traj')
# slab.write('Pd_OH_ads_opt.xyz')

# # Full GPAW restart file
# calc.write('Pd_OH_ads_opt.gpw', mode='all')

# # Energy
# E_slab_oh = slab.get_potential_energy()
# print(f"Optimized energy of the Pd+OH system: {E_slab_oh:.4f} eV")

# with open('energies.txt', 'a') as f:
#     f.write(f"E_slab_oh: {E_slab_oh:.4f}\n")