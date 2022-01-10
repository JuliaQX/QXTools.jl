```@meta
CurrentModule = QXTools
```

# QXTools

QXTools is a Julia package for simulating quantum circuits using tensor network approaches. It targets large distributed memory clusters with hardware accelerators. It was developed as part of the QuantEx project, one of the individual software projects of WP8 of [PRACE](https://prace-ri.eu/) 6IP.

QXTools ties together a number of other Julia packages which are also part of the QuantEx project. These include [QXZoo](https://juliaqx.github.io/QXZoo.jl/stable/) for generating and manipulating quantum circuits, [QXTns](https://juliaqx.github.io/QXTns.jl/stable/) for representing and manipulating tensor networks, [QXGraphDecompositions](https://juliaqx.github.io/QXGraphDecompositions.jl/stable/) which implements a number of graph algorithms for finding good contraction plans and finally [QXContexts](https://juliaqx.github.io/QXContexts.jl/stable/) which is designed to run on large distributed clusters and carry out the computations using input files generated using QXTools.

The design and implementation of QXTools and related packages was inspired by many other frameworks and packages including [ITensors](https://itensor.org/), [TensorOperations.jl](https://github.com/Jutho/TensorOperations.jl), [OMEinsum.jl](https://github.com/under-Peter/OMEinsum.jl), [Yao.jl](https://github.com/QuantumBFS/Yao.jl), [TAL-SH](https://github.com/DmitryLyakh/TAL_SH) and [ExaTN](https://github.com/ORNL-QCI/exatn).

## Statement of need

As quantum processing devices continue to scale and the algorithms and experiments being run on them grow in complexity, simulations of these systems become much more computationally demanding. To reduce the turnaround time and allow larger systems to be simulated it is necessary to move beyond single workstations and use distributed compute clusters. QXTools provides a flexible, extensible open source framework for performing these simulations. The use of [Julia](https://epubs.siam.org/doi/10.1137/141000671) makes it easy for the code to be understood, modified and extended while not sacrificing performance compared to compiled languages.

## Where to begin

For getting QXTools installed and setup, see the [Getting Started](@ref) section which has instructions on how to install QXTools and some hello world examples.
The "Tutorials" section contains some more in-depth examples and the "Users Guide" has more details of the design and structure of QXTools.

```@contents
Pages = ["users_guide.md"]
Depth = 2
```
