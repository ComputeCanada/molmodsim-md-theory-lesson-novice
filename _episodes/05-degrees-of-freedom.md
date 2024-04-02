---
title: "Degrees of Freedom"
teaching: 10
exercises: 0
questions:
- "How is kinetic energy contained and distributed in a dynamic molecular system"
- "Why constraints are used in MD simulations, and how they can affect dynamics"
objectives:
- ""
keypoints:
- "Degrees of Freedom in Rigid bodies."
- "Constraints decrease the number of degrees of freedom"
- "Imposing constraints can affect simulation outcome"
---
### Molecular degrees of freedom 
- The number of unique ways a molecule may move (increase its kinetic energy).
- Describe how kinetic energy is contained and distributed in a molecule. 
- Translational, rotational, and vibrational components. 

#### Equipartition theorem
- The energy is shared equally among the degrees of freedom (equipartition theorem). 
- Each degree of freedom has an average energy of $$\frac{1}{2}k_BT$$ and contributes $$\frac{1}{2}k_B$$ to the system's heat capacity.

#### Degrees of freedom and thermodynamics properties
- When kinetic energy is applied to simulation systems containing molecules with different degrees of freedom, the temperature increase will vary.
- If energy is spread over many places (degrees of freedom), the temperature will change less.
- The number of degrees of freedom is an essential quantity for estimating various thermodynamic variables for a simulation system (such as heat capacity, entropy, and temperature).

## Translational degrees of freedom
- Atoms and molecules have three degrees of freedom associated with the translation of their centers of mass about each coordinate axis.
- Translational degrees of freedom in three dimensions yield $$\frac{3}{2}k_BT$$  of energy.

## Rotational degrees of freedom
- Atoms have a negligible amount of rotational energy because their mass is concentrated in the nucleus which is very small (about $$10^{-15}$$ m).
- A linear molecule, has two rotational degrees of freedom. 
- A nonlinear molecule, where the atoms do not lie along a single axis, has three rotational degrees of freedom, because it can rotate around any of three perpendicular axes.
- The rotational degrees of freedom contribute $$k_BT$$  to the energy of linear molecules and $$\frac{3}{2}k_BT$$ to the energy of non-linear molecules.

## Vibrational degrees of freedom
- A diatomic molecule has one molecular vibration mode.
- A linear molecule with *N* atoms has *3N − 5* vibrational modes.
- A non-linear molecule with *N* atoms has *3N − 6* vibrational modes.

| Angle bend                                 | Symmetric stretch                               | Asymmetric stretch                               |
|--------------------------------------------|-------------------------------------------------|--------------------------------------------------|
|![bend]({{ page.root }}/fig/vibr_modes/water_1.gif){: width="100"} |![symmetric]({{ page.root }}/fig/vibr_modes/water_2.gif){: width="100"} |![asymmetric]({{ page.root }}/fig/vibr_modes/water_3.gif){: width="100"}|

- There are two degrees of freedom for each vibrational mode.
- One degree involves the kinetic energy of the moving atoms, and the other involves the potential energy of the springlike chemical bond(s).
- Each vibrational degree of freedom contributes $$k_BT$$ to the energy of a molecule. However, this is valid only when $$k_BT$$ is much bigger than energy spacing between vibrational modes.  
- At low temperature this condition is not satisfied, only a few vibrational states are occupied and the equipartition principle is not typically applicable.

#### Increasing efficiency of thermodynamic sampling.  
By reducing the number of degrees of freedom we can increase thermodynamic sampling efficiency. 

- Force fields remove the electrons’ degrees of freedom by replacing them with atom centered charges. 
- An implicit solvent model eliminates the degrees of freedom associated with the solvent molecules.
- Bond constraints eliminate vibrational degrees of freedom and make possible to use longer time steps.
- Constraints including angles and dihedrals can be also applied. 

#### Reduction of the number of degrees of freedom may lead to artifacts. 
Bond and angle constraints can: 
{: .instructor_notes :}
- slow down dihedral angle transitions [[1]](https://aip.scitation.org/doi/10.1063/1.453488) 
- shift the frequencies of the normal modes in biomolecules [[2]](https://aip.scitation.org/doi/10.1063/1.455654)
- perturb the dynamics of polypeptides [[3]](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.52.6868).

