```@meta
CurrentModule = QXTools
```

# QXTools

QXTools is a Julia package for simulating quantum circuits using tensor network approaches. It targets large distributed memory clusters with hardware accelerators. It was developed as part of the QuantEx project, one of the individual software projects of WP8 of [PRACE](https://prace-ri.eu/) 6IP.

QXTools ties together a number of other Julia packages which are also part of the QuantEx project. These include QXZoo for generating and manipulating quantum circuits, QXTns for representing and manipulating tensor networks, QXGraphDecompositions which implements a number of graph algorithms for finding good contraction plans and finally QXContexts which is designed to run on large distributed clusters and carry out the computations using input files generated using QXTools.

The design and implementation of QXTools and related packages was inspired by many other frameworks and packages including ITensors, TensorOperations.jl, OMEinsum.jl, Yao.jl, TAL-SH and ExaTN.

## Where to begin

For getting QXTools installed and setup, see the [Getting Started](@ref) section which has instructions on how to install QXTools and some hello world examples.
The "Tutorials" section contains some more in-depth examples and the "Users Guide" has more details of the design and structure of QXTools.

```@contents
Pages = ["users_guide.md"]
Depth = 2
```
