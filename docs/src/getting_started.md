# Getting Started

## Installation

QXTools is a Julia package and can be installed using Julia's inbuilt package manager from the Julia REPL using.

```
import Pkg
Pkg.add("QXTools")
```

## Example usage

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
plan = flow_cutter_contraction_plan(tnc; time=10)

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

generate_simulation_files(circ; 
                          number_bonds_to_slice=2, 
                          output_prefix="ghz_3",
                          output_method=:uniform
                          num_outputs=4)
```

will generate the files:
- `ghz_3.qx`: A DSL file with instructions
- `ghz_3.jld2`: A data file with tensors
- `ghz_3.yml`: A parameter file with parameters controlling the simulation

These can be used as input to QXContexts to run the simulation on HPC clusters to calculate the amplitudes for 4 bitstrings sampled uniformly.