# Getting Started

## Installation

QXSim is a Julia package and can be installed using Julia's inbuilt package manager from the Julia REPL using.

```
import Pkg
Pkg.add("QXSim")
```

## Example usage

An example of how QXSim can be used to calculate a set of amplitudes for small GHZ preparation circuit looks like

```
using QXSim
using QXSim.Circuits
using QXTns
using QXGraph

# Create ghz circuit
circ = create_ghz_circuit(3)

# Convert the circuit to a tensor network circuit
tnc = convert_to_tnc(circ)

# Find a good contraction plan
plan = quickbb_contraction_plan(tnc)

# Contract the network using this plan to find the given amplitude for different outputs
@show QXSim.single_amplitude(tnc, plan, "000")
@show QXSim.single_amplitude(tnc, plan, "111")
@show QXSim.single_amplitude(tnc, plan, "100")
```

This is only recommended for small test cases. For larger scale runs one can call the `generate_simulation_files`
which will do the conversion to a network, find the contraction plan and create output files describing the required
calculations. For example

```
using QXSim
using QXSim.Circuits

# Create ghz circuit
circ = create_ghz_circuit(3)

generate_simulation_files(circ, 2, "ghz_3", 4)
```

will generate the files:
- `ghz_3.qx`: A DSL file with instructions
- `ghz_3.jld2`: A data file with tensors
- `ghz_3.yml`: A parameter file with parameters controlling the simulation

These can be used as input to QXRun to run the simulation on HPC clusters to calculate the amplitudes for 4 bitstrings sampled uniformly.