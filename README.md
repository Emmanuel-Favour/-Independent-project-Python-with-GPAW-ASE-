# Density Functional Theory Study of OH Adsorption on Pd Surface (GPAW/ASE)

This project investigates the adsorption behavior of hydroxyl (OH) radicals on a palladium (Pd) surface using Density Functional Theory (DFT) as implemented in GPAW with the Atomic Simulation Environment (ASE).
The study aims to compute the binding energy of OH on Pd(111) and understand the surface–adsorbate interaction at the atomic level, relevant to heterogeneous catalysis and fuel cell reactions.

### Project Overview
The adsorption of OH species on transition metal surfaces plays a critical role in several catalytic processes, including oxygen reduction, methanol oxidation, and hydrogen evolution reactions.
This project combines computational modeling and electronic structure theory to:
*Determine the equilibrium lattice constant of bulk Pd.
* Construct and relax a Pd(111) surface slab.
* Optimize the OH radical in gas phase.
* Compute the adsorbed Pd–OH system energy.
* Derive the OH binding energy using total energy differences.

### Computational Methodology
All calculations were carried out using GPAW (LCAO mode) within the PBE generalized gradient approximation (GGA) for exchange–correlation.
The ASE Python library was used for model construction, geometry optimization, and post-processing.
* Pd Slab Construction
    * The equilibrium lattice constant of bulk Pd was first optimized using an Equation of State (EOS) calculation, yielding
    a₀ = 3.8404 Å, which closely matches experimental and reported DFT values for fcc Pd.
    * Using this lattice constant, a Pd(111) surface was cleaved from the bulk structure.
    * Supercell size: 3 × 3 × 4
    * Vacuum spacing: 10 Å to eliminate slab–slab interactions.
    * The bottom three atomic layers were fixed during geometry optimization using FixAtoms to simulate bulk-like behavior, while the top layer was fully relaxed.
 * Exchange–Correlation and Basis
   
| Parameter                       | Setting                                                |
| ------------------------------- | ------------------------------------------------------ |
| Exchange–Correlation Functional | PBE                                                    |
| Basis                           | LCAO with `'sz(dzp)'` for Pd and `'szp(dzp)'` for O, H |
| Spin Polarization               | Enabled                                                |
| Smearing                        | Fermi–Dirac (width = 0.10 eV)                          |
| K-point Sampling                | Density = 2.0 Å⁻¹, Gamma-centered                      |
| Mixer                           | `MixerDif(beta=0.05, nmaxold=10, weight=100.0)`        |

* Optimization Parameters
  
| Parameter                       | Setting                                                |
| ------------------------------- | ------------------------------------------------------ |
| Exchange–Correlation Functional | PBE                                                    |
| Basis                           | LCAO with `'sz(dzp)'` for Pd and `'szp(dzp)'` for O, H |
| Spin Polarization               | Enabled                                                |
| Smearing                        | Fermi–Dirac (width = 0.10 eV)                          |
| K-point Sampling                | Density = 2.0 Å⁻¹, Gamma-centered                      |
| Mixer                           | `MixerDif(beta=0.05, nmaxold=10, weight=100.0)`        |

* Energy Calculations
  
| System         | Total Energy (eV) |
| -------------- | ----------------- |
| Pd slab        | -95.7328          |
| OH (gas phase) | -7.2882           |
| Pd + OH        | -104.0449         |

The binding energy of OH on Pd(111) was evaluated using:

Δ𝐸 binding = 𝐸 slab+OH − (𝐸slab + 𝐸OH)
Δ𝐸 binding = −104.0449 − (−95.7328 + −7.2882) = −1.0239eV

This negative value indicates stable and exothermic adsorption, confirming that OH strongly binds to the Pd(111) surface.


<img width="1309" height="786" alt="image" src="https://github.com/user-attachments/assets/b94ee7e2-e2d4-4e3a-83c3-c975e2c0227b" />

*The image visualizes the optimized adsorption configuration of OH on the Pd(111) surface, confirming that the molecule preferentially binds at the atop site with an adsorption energy of –1.02 eV. This interaction alters the local electron distribution of surface Pd atoms, which is a key factor in determining Pd’s catalytic behavior in reactions such as oxygen reduction and hydrogen evolution.*


### Key Observations
* The Pd(111) surface remained structurally stable during relaxation.
* The optimized OH orientation favored the fcc hollow site, consistent with literature.
* Forces converged smoothly using BFGS, confirming an optimized geometry.
* The computed binding energy (–1.02 eV) aligns with DFT-reported trends for OH adsorption on transition metal surfaces.

### Files and Outputs
| File                | Description                                                                                                                |
| ------------------- | -------------------------------------------------------------------------------------------------------------------------- |
| `Pd_Lattice2.ipynb` | Jupyter notebook containing code for Pd lattice constant optimization, Pd slab construction, and OH gas-phase optimization |
| `pd_slab_opt.gpw`   | Optimized Pd slab checkpoint file (contains wavefunctions and charge density)                                              |
| `pd_slab_opt.txt`   | GPAW text output for Pd slab optimization                                                                                  |
| `pd_slab_opt.traj`  | ASE trajectory file of Pd slab relaxation                                                                                  |
| `pd_slab_opt.log`   | Log file recording optimization steps and forces for Pd slab                                                               |
| `Pd_OH.gpw`         | Final Pd–OH system checkpoint file after adsorption optimization                                                           |
| `Pd_OH.txt`         | GPAW text output for Pd–OH adsorption                                                                                      |
| `Pd_OH.traj`        | ASE trajectory file for Pd–OH adsorption relaxation                                                                        |
| `Pd_OH.log`         | Optimization log showing convergence behavior for Pd–OH system                                                             |
| `OH_gas_opt.gpw`    | Optimized OH gas-phase structure checkpoint                                                                                |
| `OH_gas_opt.txt`    | GPAW text output for gas-phase OH optimization                                                                             |
| `OH_gas_opt.traj`   | Trajectory file for OH gas-phase relaxation                                                                                |
| `OH_gas_opt.log`    | Log file of forces and energies during OH optimization                                                                     |
| `rerun_pdoH.py`     | Python script used for rerunning or continuing Pd–OH adsorption calculations                                               |


### Conclusion
This study demonstrates how first-principles DFT calculations can accurately predict adsorption energetics on metallic surfaces.
By systematically optimizing the lattice, surface, and adsorbate, and reusing converged wavefunctions, the workflow achieves physically consistent and computationally efficient results.

Binding Energy (ΔE₍binding₎): –1.0239 eV
This insight contributes to understanding catalytic activity trends on Pd surfaces — foundational for energy conversion and electrocatalytic applications.

### Tools & Environment
| Tool                  | Version / Notes  |
| --------------------- | ---------------- |
| GPAW                  | 24.x (LCAO mode) |
| ASE                   | 3.23+            |
| Python                | 3.10+            |
| Parallelization       | 4 CPU cores    |
| Vacuum                | 10 Å             |
| Convergence Threshold | fmax = 0.1 eV/Å  |

### Citation / Reference
J. Enkovaara et al., J. Phys.: Condens. Matter 22, 253202 (2010) — GPAW
A. Hjorth Larsen et al., J. Phys.: Condens. Matter 29, 273002 (2017) — ASE
