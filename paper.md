---
title: 'QXTools: A Julia framework for distributed quantum circuit simulation'
tags:
  - Julia
  - Qunatum computing
  - Simulation
  - Tensor networks
authors:
  - name: Momme Allalen
    affiliation: 1
  - name: David Brayford
    affiliation: 1
  - name: John Brennan
    affiliation: 2
  - name: Myles Doyle
    affiliation: 2
  - name: Kenneth Hanley
    affiliation: 2
  - name: Luigi Iapichino
    affiliation: 2
  - name: Niall Moran^[corresponding author]
    orcid: 0000-0002-2619-5040
    affiliation: 2
#    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Lee O’Riordan
    affiliation: 2
affiliations:
 - name: Leibniz Supercomputing Centre
   index: 1
 - name: Irish Centre for High-End Computing, Ireland
   index: 2
date: 12 August 2021
bibliography: paper.bib

---

# Summary

QXTools is a framework for simulating quantum circuits using tensor network methods.
It is designed to run on large distributed compute clusters and supports GPU
accelerators. The simulation workflow is broken down into a number of stages, each
of which is managed by a special purpose package which can be used independently or
as part of the QXTools framework. A domain specific language (DSL) is used to
represent the simulation as a set of tensor network operations.
This separates the high level index accounting and contraction planning from the
low level implementation of the tensor network operations and makes it easier to
support new hardware and network architectures.

# Statement of need

As quantum processing devices continue to scale and the algorithms and experiments
being run on them grow in complexity, simulations of these systems become much
more computationally demanding. To reduce the turnaround time and
allow larger systems to be simulated it is necessary to move beyond single workstations
and use distributed compute clusters. QXTools provides a flexible, extensible open source
framework for performing these simulations. The use of Julialang
[@Bezanson:2017] makes it easy for the code
to be understood, modified and extended while not sacrificing performance compared to
compiled languages.

# Functionality and design

QXTools consists of a number of Julialang packages available under the
[@JuliaQX] organization and registered in the Julialang package registry.
The QXTools.jl package ties these together to
enable circuit simulation workflows. The individual packages and their
roles are as follows:

- QXTns.jl: Provides data structures representing tensor networks and tensor network
circuits along with functionality for contracting these and keeping track of tensor
indices and hyper-indices.
- QXGraphDecompositions.jl: Provides specialised graph data structures and algorithms
for finding good contraction orderings and edges to slice to decompose computations.
- QXContexts.jl: Provides computation tree data structures to represent computations
and the ability to execute these compute graphs on different hardware platforms.
- QXZoo.jl: Quantum circuit representations and manipulation functionality.
- YaoQX.jl: Enables QXTools to be used as a backend for Yao.jl [@Luo:2020].

# Acknowledgements

This work was financially supported by the PRACE project funded in part by the EU’s
Horizon 2020 Research and Innovation programme (2014-2020) under grant agreement 823767.

# References