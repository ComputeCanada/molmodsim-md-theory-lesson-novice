---
title: "Electrostatic Interactions"
teaching: 15
exercises: 0
questions:
- "How electrostatic interactions are calculated in periodic systems?"
objectives:
- "Learn what parameters control the accuracy of electrostatic calculations"
keypoints:
- "Calculation of electrostatic potentials is the most time consuming part of any MD simulation"
- "Long-range part of electrostatic interacions is calculated by approximating Coulomb potentials on a grid" 
- "Denser grid increases accuracy, but significantly slows down simulation"
---

Due to the long-range behavior of Coulomb interactions, the task of computing Coulomb potentials is often the most time consuming part of any MD simulation.  Therefore, fast and efficient algorithms are required to accelerate these calculations.

Particle Mesh Ewald (PME) method is the most widely used method using the Ewald decomposition technique, the potential is divided into two parts; a real space part (computed in real space) and a Fourier space part (computed in Fourier space). The real space sum is short-ranged and as for the LJ potentials, it can be truncated when sufficiently decayed. The Fourier space part, on the other hand, is long-ranged but smooth and periodic. Therefore, its Fourier spectrum decays rapidly in Fourier space. This sum can then be accelerated by utilizing fast Fourier transforms (FFT).  The method requires periodic boundary conditions and charge neutrality of the molecular system in order to accurately calculate the total Coulombic interaction.

Parameters controlling PME calculations are listed in the Table below.

| Variable \ MD package | GROMACS                  | NAMD                      | AMBER                  |
|-----------------------|--------------------------|---------------------------|------------------------|
| Fourier grid spacing  | fourierspacing (1.2)     | PMEGridSpacing  (1.5)     |                        |
| Grid Dimension X      | fourier-nx               | PMEGridSizeX              |  nfft1                 |
| Grid Dimension Y      | fourier-ny               | PMEGridSizeY              |  nfft2                 |
| Grid Dimension Z      | fourier-nz               | PMEGridSizeZ              |  nfft3                 |
| Direct space tolerance| ewald-rtol ($$10^{-5}$$) | PMETolerance ($$10^{-5}$$)| dsum_tol ($$10^{-6}$$) |
| Interpolaton order    | pme-order (4)            | PMEInterpOrder (4)        | order (4)              |

**Grid dimension** values give the size of the charge grid (upon which the reciprocal sums are interpolated) in each dimension. Higher values lead to higher accuracy (when the direct space tolerance is also lowered) but considerably slow the calculation. Generally it has been found that reasonable results are obtained when grid dimension values are approximately equal to periodic box size, leading to a grid spacing of 1.0 Å. Significant performance enhancement in the calculation of the fast Fourier transform is obtained by having each of the values be a product of powers of 2, 3, and/or 5.  If the values are not given, programs will chose values to meet these criteria.

**Interpolation order** is the order of the B-spline interpolation. The higher the order, the better the accuracy (unless the charge grid is too coarse). The minimum order is 3. An order of 4 (the default) implies a cubic spline approximation which is a good standard value.  