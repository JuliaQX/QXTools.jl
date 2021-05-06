# User's Guide

QXTools generates output files which provide a description of the computations and data
required to perform a simulation. There are three output files:

- Parameter file: This YAML file provides informations on the sliced edges and dimensions as well
as details of the output amplitudes or sampling method
- DSL file: The DSL file contains instructions describing the tensor operations involved in
performing the contraction
- Data file: The data file contains the numerical values of the initial tensors of the tensor
network

## Parameter File
To distribute the computation over multiple processes and nodes, a number of edges of the tensor
network are sliced. The computation is then divided into multiple smaller computations, each
having a particular assignment of values on the sliced edges. The details of which edges are sliced
and the dimension of these is described in the parameter file. In addition to this, the outputs
to use when contracting the network are also described in the parameter file. An example of a parameter
file is as follows

```
amplitudes:
  - "00"
  - "01"
partitions:
  parameters:
    v2: 2
    v1: 2
```

This tells us there are two edges that will be sliced, each with dimension 2 and that the computation
should be performed by setting the outputs to each of the four bitstrings given in the file. In future
it will be possible to specify parameter of the sampling method to use for selecting which bitstrings
to calculate the amplitude of instead of each being listed.

The above parameter file will result in the contraction described in the DSL file being performed eight times
with different substitutions in each case. The full list of substitutions is

```@eval
using DataFrames
using Latexify
df = let
    df = DataFrame("v1" => Int[], "v2" => Int[], "o1" => String[], "o2" => String[])
    for outs in CartesianIndices((2, 1))
        for vs in CartesianIndices((2, 2))
            push!(df,[ vs[2], vs[1], "o1_$(outs[2])", "o2_$(outs[1])"])
        end
    end
    df
end
mdtable(df,latex=false)
```

The final scalar output of the four different substitutions for each output set are summed to get the
final output amplitude. This means that each of these contractions can be performed independently on
different processes/nodes with a reduction operation performed on the scalars resulting from each contraction.

## DSL Specification

A Domain Specific Language (DSL) has been defined when enables better separation
of concerns between the components. This DSL describes a contraction over a tensor
network to calculate observables of interest in the form of a sequence of instructions
which act on individual tensors. This makes it possible to separate the development
of the high performance distributed tensor network computation code from that of the
higher level contraction planning, circuit and and network manipulation code. For QXTools
DSL files the ".qx" suffix is used for "Tensor Language". We will first show a very
simple example of a QXTools DSL file and then descibe in detail how each of the
instructions work.

### Example DSL file

An example of the DSL generated for the contraction of a two qubit GHZ circuit looks like.

```
# version: 0.3.0
# Determination of contraction plan :
#   Method used : flow cutter
#   Treewidth : 2
#   Time allocated : 30
#   Seed used : -1
#   Returned metadata :
#     1 : c min degree heuristic
#     2 : c status 3 1620299104065
#     3 : c min shortcut heuristic
#     4 : c run with 0.0/0.1/0.2 min balance and node_min_expansion in endless loop with varying seed
#   Hypergraph used : true
#   Hyperedge contraction method : Netcon where possible, min fill heuristic otherwise.
# Slicing :
#   Method used : greedy treewidth deletion
#   Edges sliced : 2
#   Score fucntion used : direct_treewidth
#   Treewidths after slicing consecutive edges : [1, 0]
outputs 2
load t1 data_1
load t2 data_2
load t3 data_3
load t4 data_4
load t5 data_4
view t2_s t2 3 v1
view t3_s t3 1 v1
view t1_s t1 1 v1
view o1_s o1 1 v1
view t5_s t5 1 v2
view t2_s_s t2_s 2 v2
ncon t8 2 t3_s 2 o1_s 2
ncon t9 1,3 t8 1 t5_s 3
ncon I1 2,3 t2_s_s 1,2,3 o2 1
ncon t10 1 t9 1,3 I1 3,1
ncon I2 1 t1_s 1,2 t4 2
ncon t11 0 t10 1 I2 1
save t11 output
```

### DSL Format and Instructions

The DSL file is a regular ASCII text file with one instruction per line. Lines beginning with `#` are comments which are
ignored (except for the first line which contains version information). The first line has a version string
which specifies the format version. Comments following the first line contain metadata about the methods used
to determine the contraction plan used and which edges to slice.

#### Outputs
```
outputs 2
```

This line indicates that the tensor network has two output tensors. These tensors are represented in subsequent DSL
commands by symbols in the format `o{i}`, where `{i}` is the index of the given output.
To contract the tensor network with different values for the output tensors, the data these symbols point to needs to be updated.

#### Load

Load instructions define a new tensor symbol using data from the input file. Here the first argument is the
name of the new symbol and the second argument is the key to find the data in the input [Data File](@ref).
It is normal that multiple tensor symbols have the same data associated with them, e.g. multiple occurrences of the same
gate. In the following example both `t4` and `t5` both used the data labeled `data_4` in the input [Data File][@ref]

```
load t1 data_1
load t2 data_2
load t3 data_3
load t4 data_4
load t5 data_4
```

#### View

A view instruction creates a new tensor symbol by taking a view of an existing tensor.
The first argument is the new symbol name, the second is the symbol of the tensor to take
the view of, the third argument is the dimension and the final value is the particular value
for that dimension. An example view command is as follows

```
view t13_s t13 3 v1
```

The `v1` variable will be replaced by a particular value before execution as described in
the [Parameter File](@ref) section.

#### Contraction

The `ncon` instruction specifies a pairwise contraction of tensors.
An example contraction command is as follows

```
ncon t15 2,3 o1 1 t14_s_s 2,1,3
```

Which indicates that the `o1` tensor should be contracted with the `t14_s_s` tensor to get the `t15` tensor.
The contraction should be performed over the only rank of the `o1` tensor and the second rank of the `t14_s_s` tensor.
Einstein summation notation is used for for specifying which indices to contract over.
This convention uses repeated indices on the right hand side to indicate that those indices should be contracted over.
In the above, the `1` index appears twice on the right hand side which indicates that this index should be contracted over.
For the case where one of the tensors is a scalar, a `0` is used as a placeholder.
For example if `t1` is a scalar tensor and `t2` is a matrix,
the multiplication of the matrix by the scalar can be expressed as

```
ncon t3 1,2 t1 0 t2 1,2
```

Batched contractions (sometimes called hyper-contractions) are also supported. This is where the same index
is repeated on the right hand side and also appears on the left hand side.

#### Save

The save instruction indicates that the given tensor is an output from the contraction and provides a label for this.
The value of this tensor will then be used in a reduction operation and/or written to the output file.


## Data File

The data file is a [JLD](https://github.com/JuliaIO/JLD.jl) file and contains the numerical values of the initial tensors.
Each tensor is stored as a multi-dimensional array with a data label that is referenced in `load` commands of the DSL file.
