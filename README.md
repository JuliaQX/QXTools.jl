# QXTools

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaQX.github.io/QXTools.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaQX.github.io/QXTools.jl/dev)
[![Build Status](https://github.com/JuliaQX/QXTools.jl/workflows/CI/badge.svg)](https://github.com/JuliaQX/QXTools.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaQX/QXTools.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaQX/QXTools.jl)

QXTools is a Julia package for simulating quantum circuits using tensor networking approaches. It targets large distributed memory clusters with hardware
accelerators. It was developed as part of the QuantEx project, one of the individual software projects of WP8 of PRACE 6IP.

QXTools depends on a number of other Julia packages developed that were also developed as part of the QuantEx project. These include QXZoo which
is capable of generating and manipulating quantum circuits, QXTns which features data structures and functions for manipulating tensor networks,
QXGraphDecompositions which implements a number of graph algorithms for finding good contraction plans and QXContexts which is designed to run on large distributed
clusters.

The design and implementation of QXTools and related packages was inspired by many other frameworks and packages including ITensors.jl, TensorOperations.jl,
Yao.jl, TAL-SH and ExaTN.

# Installation

QXTools is a Julia package and can be installed using Julia's inbuilt package manager from the Julia REPL using.

```
import Pkg
Pkg.add("QXTools")
```

# Example usage

An example of how QXTools can be used to calculate a set of amplitudes for small GHZ preparation circuit looks like

```
using QXTools
using QXTools.Circuits
using QXTns
using QXGraphDecompositions

# Create ghz circuit
circ = create_ghz_circuit(3)

# Convert the circuit to a tensor network circuit
tnc = convert_to_tnc(circ)

# Find a good contraction plan
plan = quickbb_contraction_plan(tnc)

# Contract the network using this plan to find the given amplitude for different outputs
@show QXTools.single_amplitude(tnc, plan, "000")
@show QXTools.single_amplitude(tnc, plan, "111")
@show QXTools.single_amplitude(tnc, plan, "100")
```

This is only recommended for small test cases. For larger scale runs one can call the `generate_simulation_files`
which will do the conversion to a network, find the contraction plan and create output files describing the required
calculations. For example

```
using QXTools
using QXTools.Circuits

# Create ghz circuit
circ = create_ghz_circuit(3)

generate_simulation_files(circ, 2, "ghz_3", 4)
```

will generate the files:
- `ghz_3.tl`: A DSL file with instructions
- `ghz_3.jld`: A data file with tensors
- `ghz_3.yml`: A parameter file with parameters controlling the simulation

These can be used as input to QXContexts to run the simulation on HPC clusters to calculate the amplitudes for 4 bitstrings sampled uniformly.
For more details and options see the documentation at [docs](https://juliaqx.github.io/QXContexts.jl/dev/).

# Contributing
Contributions from users are welcome and we encourage users to open issues and submit merge/pull requests for any problems or feature requests they have. The
[CONTRIBUTING.md](CONTRIBUTION.md) has further details of the contribution guidelines.


# Building documentation

QXTools.jl using [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/) to generate documentation. To build
the documentation locally run the following from the `docs` folder.

The first time it is will be necessary to instantiate the environment to install dependencies

```
julia --project 'import Pkg; Pkg.instantiate()'
```

and then to build the documentation

```
julia --project make.jl
```

To serve the generated documentation locally use

```
julia --project -e 'using LiveServer; serve(dir="build")'
```

Or with python3 using from the `docs/build` folder using

```
python3 -m http.server
```

The generated documentation should now be viewable locally in a browser at `http://localhost:8000`.