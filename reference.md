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
  this regexp (still untested): `s/ doi:([^ ]+)/ [\1](http://dx.doi.org/\1)/`

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
    J Chem Phys. 2010;133: 124515. [doi:10.1063/1.3481102](http://dx.doi.org/doi:10.1063/1.3481102)

Allen-2017
:   Allen MP, Tildesley DJ.
    *Computer Simulation of Liquids. Second Edition*
    Oxford University Press; 2017. ISBN: 9780198803195 [doi:10.1093/oso/9780198803195.001.0001](http://dx.doi.org/10.1093/oso/9780198803195.001.0001)

Basconi-2013
:   Basconi JE, Shirts MR.
    *Effects of Temperature Control Algorithms on Transport Properties and Kinetics in Molecular Dynamics Simulations.*
    J Chem Theory Comput. 2013;9: 2887–2899.  [doi:10.1021/ct400109a](http://dx.doi.org/10.1021/ct400109a)

Grosfils-2009
:    Grosfils P, Lutsko JF.
     *Dependence of the liquid-vapor surface tension on the range of interaction: a test of the law of corresponding states.*
     J Chem Phys. 2009;130: 054703. [doi:10.1063/1.3072156](http://dx.doi.org/doi:10.1063/1.3072156)

Huang-2014
:     Huang K, García AE.
      *Effects of truncating van der Waals interactions in lipid bilayer simulations.*
      Journal of Chemical Physics. 2014;141. [doi:10.1063/1.4893965](http://dx.doi.org/doi:10.1063/1.4893965)

<!--
Gowers-2016
:   Gowers RJ, Linke M, Barnoud J, Reddy TJE, Melo MN, Seyler SL, et al.
    *MDAnalysis: A Python Package for the Rapid Analysis of Molecular Dynamics Simulations.*
    Proc 15th Python Sci Conf. 2016; 98–105. Available: http://conference.scipy.org/proceedings/scipy2016/pdfs/oliver_beckstein.pdf -->

Larsson-2011
:   Larsson P, Hess B, Lindahl E.
    *Algorithm improvements for molecular dynamics simulations.*
    Wiley Interdiscip Rev Comput Mol Sci. 2011;1: 93–108. [doi:10.1002/wcms.3](http://dx.doi.org/doi:10.1002/wcms.3)

<!--
Merz-2018
:   Merz PT, Shirts MR.
    *Testing for physical validity in molecular simulations.*
    Huang X, editor. PLoS One. 2018;13: e0202764. [doi:10.1371/journal.pone.0202764](http://dx.doi.org/10.1371/journal.pone.0202764) -->
<!--
Michaud-Agrawal-2011
:    Michaud-Agrawal N, Denning EJ, Woolf TB, Beckstein O.
    *MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.*
    J Comput Chem. 2011;32: 2319–27. [doi:10.1002/jcc.21787](http://dx.doi.org/10.1002/jcc.21787) -->
Piana-2012
:    Piana S, Lindorff-Larsen K, Dirks RM, Salmon JK, Dror RO, Shaw DE.
    *Evaluating the Effects of Cutoffs and Treatment of Long-range Electrostatics in Protein Folding Simulations.*
    PLOS ONE. 2012;7: e39918. [doi:10.1371/journal.pone.0039918](http://dx.doi.org/doi:10.1371/journal.pone.0039918)

Wong-ekkabut-2016
:   Wong-ekkabut J, Karttunen M.
    *The good, the bad and the user in soft matter simulations.*
    Biochim Biophys Acta - Biomembr. Elsevier B.V.; 2016;1858: 2529–2538. [doi:10.1016/j.bbamem.2016.02.004](http://dx.doi.org/10.1016/j.bbamem.2016.02.004)

Yonetani-2006
:    Yonetani Y.
     *Liquid water simulation: A critical examination of cutoff length.* Journal of Chemical Physics. 2006;124. [doi:10.1063/1.2198208](http://dx.doi.org/doi:10.1063/1.2198208)
