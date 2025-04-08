---
title: "Supplemental: Overview of the ML Force Fields"
teaching: 0
exercises: 0
questions:
- "What are the categories of molecular dynamics force fields?"
- "What force fields are available?" 
- "What systems they work best with?"
objectives:
- "Be able to recognize the strengths and weaknesses of different types of force fields."
- "Find out which force fields are available and which systems they are most suitable for."
- "Identify the original papers that introduced force fields."
keypoints:
- "There are different types of force fields designed for different types of simulations."
- "Induction effects are not accounted for by fixed-charge force fields."
- "Using more accurate and diverse target data allows MD force fields to be improved."
---

* Table of Contents
{:toc}


## Machine Learning Force Fields

#### [ANAKIN-ME](https://roitberg.chem.ufl.edu/research/) (Accurate NeurAl networK engINe for Molecular Energies)  
ANI-1 neural networks are trained on molecules up to up to 8 atoms (H, C, N, O)
- [ANI-1x](https://doi.org/10.1063/1.5023802)  Trained on 5 million molecular conformations, ωB97X/6-31G(d). 
- [ANI-1ccx](https://www.nature.com/articles/s41467-019-10827-4) Trained on 500,000 molecular conformations, CCSD(T)/CBS.

ANI-2 neural networks added F, Cl, S.  
- [ANI-2x](https://pubs.acs.org/doi/10.1021/acs.jctc.0c00121)  Trained on 8.9 million molecular conformations, ωB97X/6-31G(d).  

#### [SO3LR](https://chemrxiv.org/engage/chemrxiv/article-details/679bf19781d2151a02991c58)  

Trained on 4 million molecular conformations, PBE0+MBD. SO3krates neural network for semi-local interactions combined with a pairwise force field for short-range repulsion, long-range electrostatics, and dispersion interactions. Scales to 200,000 atoms on a single GPU. Implemented in [JAX-MD](https://jax-md.readthedocs.io/en/main/).

#### [MACE](https://doi.org/10.1063/5.0155322)  
The MACE architecture combines Atomic Cluster Expansion framework with Equivariant Message Passing Neural Networks. Several pre-trained [models](https://mace-docs.readthedocs.io/en/latest/guide/guide.html) are available.  

#### [ViSNet](https://www.nature.com/articles/s41467-023-43720-2)
Equivariant geometry-enhanced graph neural network. Developed by Microsoft Research. A Graph neural Network with Vector-scalar interactive message passing. Avaialble in [Pytorch Geometric](https://pytorch-geometric.readthedocs.io/en/2.5.1/generated/torch_geometric.nn.models.ViSNet.html)  and [GitHub](https://github.com/microsoft/AI2BMD/tree/ViSNet).