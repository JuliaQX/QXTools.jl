---
title: 'QXTools: A Julia framework for distributed quantum circuit simulation'
tags:
  - Julia
  - Qunatum computing
  - Simulation
  - Tensor networks
authors:

  - name: John Brennan
    affiliation: 1
  - name: Lee O’Riordan
    affiliation: 1
  - name: Kenneth Hanley
    affiliation: 1
  - name: Myles Doyle
    affiliation: 1
  - name: Momme Allalen
    affiliation: 2
  - name: David Brayford
    affiliation: 2
  - name: Luigi Iapichino
    affiliation: 2
  - name: Niall Moran^[corresponding author]
    orcid: 0000-0002-2619-5040
    affiliation: 1
affiliations:
 - name: Irish Centre for High-End Computing, Ireland
   index: 1
 - name: Leibniz Supercomputing Centre
   index: 2
date: 12 August 2021
bibliography: paper.bib

---

# Summary

QXTools is a framework for simulating quantum circuits using tensor network methods.
Weak simulation is the primary use case where given a quantum circuit and input state
QXTools will efficiently calculate the probability amplitude of a given output configuration
or set of configurations. Given this ability one can sample from the output distribution using
random sampling approaches. QXTools is intended to be used by researchers interested in
simulating circuits larger than those possible with full wave-function simulators or those
interested in tensor network circuit simulation methods.

QXTools is written in Julia [@Bezanson:2017] and is designed to run on large distributed compute clusters and to support GPU accelerators. The simulation workflow is broken down into a number of stages, each of which is managed by a special purpose package which can be used independently or
as part of the QXTools framework. A domain specific language (DSL) is used to
express the simulation as a set of tensor network operations.
This separates the high level index accounting and contraction planning from the
low level implementation of the tensor network operations and makes it easier to
support new hardware and network architectures.

# Statement of need

As quantum processing devices continue to scale and the algorithms and experiments
being run on them grow in complexity, simulations of these systems become much
more computationally demanding. To reduce the turnaround time and
allow larger systems to be simulated it is necessary to move beyond single workstations
and use distributed compute clusters. QXTools provides a flexible, extensible open source
framework for performing these simulations. The use of Julia
[@Bezanson:2017] makes it easy for the code
to be understood, modified and extended while not sacrificing performance compared to
compiled languages.

# Background

Classical simulation of quantum circuits is essential for debugging and validating the accuracy of quantum computing devices and algorithms. This is a very computationally demanding problem owing to the exponential growth in the state space as the number of qubits is increased. Up until recently it has been possible to simulate the largest prototype universal quantum computers using direct evolution of the quantum state (where the full wave-function is stored in memory (or on disk)) using modest personal computing resources. With the recent emergence of Noisy Intermediate Scale Quantum (NISQ) devices, it has become intractable to use this approach to simulate devices of this size on even the largest supercomputers. These advances have necessitated the development of new circuit simulation methods which can simulate large systems without incurring the memory required to store the full wave-function. Tensor network methods have been demonstrated to achieve state of the art performance in simulating Random Quantum Circuits (RQC) as part of the quantum sumpremacy experiments [@Villalonga:2019]. Despite the impressive results achieved with these methods to date they are not suitable for all types of circuits and in many cases full wave-function methods are preferable. For example for highly entangled circuits, the tensor network representation will consume the same memory as full wave-function methods but will incur additional overhead. However for circuits with moderate entanglement and cases where one is not interested in exact results, but in results up to a particulary fidelity, tensor network approaches can offer significant advantages.

Tensor networks refer to networks of interconnected tensors, the use of which originated in the condensed matter physics and quantum information communities as a means of simulating strongly correlated many body quantum systems. Expressing a quantum circuit as a tensor network is very straightforward and involves replacing each gate with a tensor and registers with sets of interconnected tensors. Single qubit Hadamard gates remain unchanged from their matrix representation, while two qubit gates can be expressed as a single rank 4 tensor or two connected rank 3 tensors. Once the quantum circuit is expressed as a network of interconnected tensors operations can be performed by contracting tensors together. For further details there are many excellent resources on tensor networks and their use in quantum information, see [@Biamonte:2017], [@Bridgeman:2017], [@Wood:2011].

# Functionality and design

QXTools consists of a number of Julia packages available under the
[@JuliaQX] organization and registered in the Julia package registry.
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