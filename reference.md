---
layout: reference
root: .
---

## Glossary

FIXME

{:auto_ids}
barostat
:   Pressure control algorithms in molecular dynamics (MD) simulations are commonly referred to
    as barostats and are needed to study isobaric systems.  Barostats work by altering the
    size of the simulation box.  Therefore they can only be used in conjunction with
    [periodic boundary conditions (PBC)](#periodic-boundary-conditions).

ergodicity
:   In statistical mechanics ergodicity describes the principle that studying a single particle 
    averaged over a long time is equivalent to averaging over many particles which are studied
    for a short time.
    See also: ["Ergodicity" on Wikipedia](https://en.wikipedia.org/wiki/Ergodicity)

periodic boundary conditions
:   Periodic boundary conditions (PBC) ...
    See also:
    * ["Periodic boundary conditions" in Wikipedia](https://en.wikipedia.org/wiki/Periodic_boundary_conditions)
    * ["Periodic boundary conditions" on gromacs.org](http://www.gromacs.org/Documentation/Terminology/Periodic_Boundary_Conditions)

thermostat
:   Temperature control algorithms in molecular dynamics (MD) simulations are commonly referred to
    as thermostats and are needed to study isothermal systems.  Thermostats work by altering the
    velocities of particles.


## Bibliography

{% comment %}

Add items to the Bibliography using the following template and add it in the right
location sorted by RefnameYYYY.

The Style is PLOS ONE with separate lines for:
- list of Authors, (line ends with two spaces)
- *Title* (enclosed in single `*` to format as emphasized, line ends with two spaces)
- rest of the citation including the DOI if available.
  You can leave the DOI in doi:xxxxx notation. We can convert them into links using
  this regexp (still untested): `s/ doi:([^ ]+) *$/ [doi:\1](https://doi.org/\1)/`

e.g.:

Refname-YYYY
:   List OF, Author S.  
    *Title of the paper -  PLOS ONE style*  
    PLoS ONE. 2019;1. doi:10.1371/journal.pone.0000000

You can add links to this entry by using on any page:

[RefnameYYYY]({{ page.root }}/reference.html#Refname-YYYY)

{% endcomment %}

{:auto_ids}

Ahmed-2010
:   Ahmed A, Sadus RJ.  
    *Effect of potential truncations and shifts on the solid-liquid phase coexistence of Lennard-Jones fluids.*  
    J Chem Phys. 2010;133: 124515. [doi:10.1063/1.3481102](https://doi.org/doi:10.1063/1.3481102)

Allen-2017
:   Allen MP, Tildesley DJ.  
    *Computer Simulation of Liquids. Second Edition*  
    Oxford University Press; 2017. ISBN: 9780198803195 [doi:10.1093/oso/9780198803195.001.0001](https://doi.org/10.1093/oso/9780198803195.001.0001)

Andersen-1980
:   Andersen HC.  
    *Molecular dynamics simulations at constant pressure and/or temperature.*  
    J Chem Phys. 1980;72: 2384–2393. [doi:10.1063/1.439486](https://doi.org/10.1063/1.439486)

Basconi-2013
:   Basconi JE, Shirts MR.  
    *Effects of Temperature Control Algorithms on Transport Properties and Kinetics in Molecular Dynamics Simulations.*  
    J Chem Theory Comput. 2013;9: 2887–2899.  [doi:10.1021/ct400109a](https://doi.org/10.1021/ct400109a)

Berendsen-1984
:   Berendsen HJC, Postma JPM, van Gunsteren WF, DiNola A, Haak JR.  
    *Molecular dynamics with coupling to an external bath.*  
    J Chem Phys. $abstract.copyright_name.value; 1984;81: 3684–90. [doi:10.1063/1.448118](https://doi.org/10.1063/1.448118)
    
Braun-2019
:   Braun E, Gilmer J, Mayes HB, Mobley DL, Monroe JI, Prasad S, et al.  
    *Best Practices for Foundations in Molecular Simulations [Article v1.0].*  
    Living J Comput Mol Sci. 2019;1: 1–28. [doi:10.33011/livecoms.1.1.5957](https://doi.org/10.33011/livecoms.1.1.5957)

Bussi-2007
:   Bussi G, Donadio D, Parrinello M.  
    *Canonical sampling through velocity rescaling.*  
    J Chem Phys. 2007;126: 014101. [doi:10.1063/1.2408420](https://doi.org/10.1063/1.2408420)

Dauber-Osguthorpe-2018
:   Dauber-Osguthorpe P, Hagler AT.  
    *Biomolecular force fields: where have we been, where are we now, where do we need to go and how do we get there?*  
    Journal of Computer-Aided Molecular Design. 2019;33: 133–203. [doi:10.1007/s10822-018-0111-4](https://doi.org/doi:10.1007/s10822-018-0111-4)

Grosfils-2009
:   Grosfils P, Lutsko JF.  
    *Dependence of the liquid-vapor surface tension on the range of interaction: a test of the law of corresponding states.*  
    J Chem Phys. 2009;130: 054703. [doi:10.1063/1.3072156](https://doi.org/doi:10.1063/1.3072156)

Hagler-2019
:   Hagler AT.  
    *Force field development phase II: Relaxation of physics-based criteria… or inclusion of more rigorous physics into the representation of molecular energetics.*  
    J Comput Aided Mol Des. 2019;33: 205–264. [doi:10.1007/s10822-018-0134-x](https://doi.org/doi:10.1007/s10822-018-0134-x)

Hoover-1985
:   Hoover WG.  
    *Canonical dynamics: Equilibrium phase-space distributions.*  
    Phys Rev A. APS; 1985;31: 1695–1697. [doi:10.1103/PhysRevA.31.1695](https://doi.org/10.1103/PhysRevA.31.1695)

Huang-2014
:   Huang K, García AE.  
    *Effects of truncating van der Waals interactions in lipid bilayer simulations.*  
    Journal of Chemical Physics. 2014;141. [doi:10.1063/1.4893965](https://doi.org/doi:10.1063/1.4893965)

Koopman-2006
:   Koopman EA, Lowe CP.  
*Advantages of a Lowe-Andersen thermostat in molecular dynamics simulations.*  
J Chem Phys. 2006;124: 1–6. [doi:10.1063/1.2198824](https://doi.org/10.1063/1.2198824)

Larsson-2011
:   Larsson P, Hess B, Lindahl E.  
    *Algorithm improvements for molecular dynamics simulations.*  
    Wiley Interdiscip Rev Comput Mol Sci. 2011;1: 93–108. [doi:10.1002/wcms.3](https://doi.org/doi:10.1002/wcms.3)

Lifson-1968
:   Lifson S, Warshel A.  
    *Consistent Force Field for Calculations of Conformations, Vibrational Spectra, and Enthalpies of Cycloalkane and n‐Alkane Molecules.*  
    J Chem Phys. 1968;49: 5116–5129. [doi:10.1063/1.1670007](https://doi.org/doi:10.1063/1.1670007)

Martyna-1992
:   Martyna GJ, Klein ML, Tuckerman M.  
    *Nosé–Hoover chains: The canonical ensemble via continuous dynamics.*  
    J Chem Phys. 1992;97: 2635–2643. [doi:10.1063/1.463940](https://doi.org/10.1063/1.463940)

Nose-1984
:   Nosé S.  
    *A molecular dynamics method for simulations in the canonical ensemble.*  
    Mol Phys. 1984;52: 255–268. [doi:10.1080/00268978400101201](https://doi.org/10.1080/00268978400101201)

Piana-2012
:   Piana S, Lindorff-Larsen K, Dirks RM, Salmon JK, Dror RO, Shaw DE.  
    *Evaluating the Effects of Cutoffs and Treatment of Long-range Electrostatics in Protein Folding Simulations.*  
    PLOS ONE. 2012;7: e39918. [doi:10.1371/journal.pone.0039918](https://doi.org/doi:10.1371/journal.pone.0039918)

Shirts-2013
:   Shirts MR.  
    *Simple Quantitative Tests to Validate Sampling from Thermodynamic Ensembles.*  
    J Chem Theory Comput. 2013;9: 909–926. [doi:10.1021/ct300688p](https://doi.org/10.1021/ct300688p)

Verlet-1967
:   Verlet L. Computer “Experiments” on Classical Fluids. I.  
    *Thermodynamical Properties of Lennard-Jones Molecules.*  
    Phys Rev. 1967;159: 98–103. [doi:10.1103/PhysRev.159.98](https://doi.org/doi:10.1103/PhysRev.159.98)

Winger-2009
:   Winger M, Trzesniak D, Baron R, van Gunsteren WF.  
    *On using a too large integration time step in molecular dynamics simulations of coarse-grained molecular models.*  
    Phys Chem Chem Phys. 2009;11: 1934–1941. [doi:10.1039/b818713d](https://doi.org/doi:10.1039/b818713d)

Wong-ekkabut-2016
:   Wong-ekkabut J, Karttunen M.  
    *The good, the bad and the user in soft matter simulations.*  
    Biochim Biophys Acta - Biomembr. Elsevier B.V.; 2016;1858: 2529–2538. [doi:10.1016/j.bbamem.2016.02.004](https://doi.org/10.1016/j.bbamem.2016.02.004)

Yonetani-2006
:   Yonetani Y.  
    *Liquid water simulation: A critical examination of cutoff length.*  
    Journal of Chemical Physics. 2006;124. [doi:10.1063/1.2198208](https://doi.org/doi:10.1063/1.2198208)
