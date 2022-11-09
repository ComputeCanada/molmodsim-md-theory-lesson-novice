---
layout: page
title: "Scratch Pad"
permalink: /scratch-pad/
teaching: 0
exercises: 0
objectives:
- The scratch pad contains some rough ideas and unfinished content that is not yet part of the main lesson.
keypoints:
- FIXME
---

## Table of Contents

* Table of Contents
{:toc}


## Flow Charts
### The Global MD Algorithm

<div class="mermaid">
graph TD
    A("Input initial conditions") --> B(fa:fa-spinner)
    B --> C["compute forces"]
    C --> D["update configuration"]
    D --> E[output step]
    E -->|repeat *nsteps* times | B
</div>


### Main phases of a molecular dynamics simulation

<div class="mermaid">
flowchart TB
%%    subgraph work
%%        direction TB
        subgraph systemprep["System Preparation"]
            direction TB
            A1["check input structure"]
            A2["generate topology"]
            A3["define simulation box"]
            A4["add ions and solvent"]
            A1 -->  A2
            A2 -->  A3
            A3 -->  A4
        end

        subgraph minimization["Energy Minimization and Relaxation"]
            direction TB
            B1["Energy Minimization"]
            B1 --> B2["short NVT (possibly w/ position restraints)"]
            B2 --> B3["short NPT"]
        end
     
        subgraph equilibration["Equilibration"]
            direction TB
            C1["Equilibration MD"]
        end

        subgraph production["Production"]
            direction TB
            D1["Production MD and data collection"]
        end

        subgraph analysis["Analysis"]
           direction TB
           E1["Data Analysis"]
        end
    systemprep    --> minimization
    minimization  --> equilibration
    equilibration --> production
    production    --> analysis
%%    end
%%    work --- intent
%%    subgraph intent
%%        direction TB
%%        QA1(["Is the molecule whole? Are there missing atoms or residues?"])
%%        QA2(["Check whether each component as integer charge."])
%%        QA4(["Neutralize system."])
%%    A1 -.-> QA1
%%    A2 -.-> QA2
%%    A4 -.-> QA4
%%        B1 -.-> QB1(["Resolve atomic clashes that lead to high forces."])
%%        B2 -.-> QB2(["Let system relax to target temperature."])
%%        B3 -.-> QB3(["Let system relax to target density."])
%%    end


    
</div>

### The flow of data

<div class="mermaid">
graph LR
    inp1(coordinates)   ---> md{MD ENGINE}
    inp2(topology)      ---> md
    inp3(md parameters) ---> md
    md   --->   out1(log file as plain text)
    md   --->   out2(trajectory of coordinates)
    md   --->   out3(final coordinates)
    md   --->   out4(checkpoint / restart file)
    md   --->   out5(energy components, measurements, ...)
</div>

## MD-parameters

### Time Step

The maximum time step with which MD simulations can be expected to be stable depends on the fastest oscillations of the model.


| conditions                                       | typical maximum time step |
|:-------------------------------------------------|:-------------------------:|
| simulation without bond- or angle-constraints    |     0.5 fs - 1 fs         |
| covalent bonds constrained                       |              2 fs         |
| using virtual interaction sites                  |            ~ 5 fs         |
| using a coarse-grained forcefield (e.g. MARTINI) |        20 - 40 fs         |
{: .container-md  }

Section "Choosing an appropriate timestep" in [Braun-2019]({{ page.root }}/reference.html#braun-2019)


### Cut-Offs

The cut-off scheme should match what was used for parameterization of the forcefield.
It has been shown that using PME for determining long-range electrostatics in many cases improves the results, even if the forcefield has not been parameterized with PME. (FIXME: citations needed)

FIXME: citations needed


### Force Fields

Different force fields have been parameterized with different philosophies and tuned to reproduce different physical properties or with emphasis on different kinds of molecules.

Generally speaking, _different forcefield **cannot** be combined_ unless they are members of the same family and therefore the secondary forcefield has been specifically parameterized to be compatible with the parent-forcefield.

#### A Selection of Common Forcefield Families

* AMBER (Assisted Model Building with Energy Refinement)
    * GAFF (General AMBER force field) for small molecules
    * GLYCAM for carbohydrates
    * LIPID21 for lipids
* CHARMM (Chemistry at Harvard Macromolecular Mechanics)
    * CGenFF (CHARMM General Force Field) for small molecules
* GROMOS (GROningen MOlecular Simulation)
* OPLS (Optimized Potentials for Liquid Simulations)
    * OPLS-UA (OPLS United Atoms forcefield)
    * OPLS-AA (OPLS All Atoms forcefield)
* MARTINI - a general purpose Coarse-Grained Force Field
