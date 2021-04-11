# User's Guide

QXTools generates output files which provide a description of the computations and data
required to perform the simulation. There are three output files:

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
# version: 0.2.0
# Determination of contraction plan :
#   Method used : quickbb
#   Time allocated : 120
#   Ordering used : min_fill
#   Lower bound flag used : false
#   Returned metadata :
#     treewidth : 2
#     is_optimal : true
#     time : 0.000147
#   Hypergraph used : false
#   Hyperedge contraction method : Netcon where possible, min fill heuristic otherwise.
# Slicing :
#   Method used : greedy treewidth deletion
#   Edges sliced : 2
#   Score fucntion used : direct_treewidth
#   Treewidths after slicing consecutive edges : [1, 1]
outputs 2
load t1 data_1
load t2 data_2
load t3 data_3
load t4 data_4
load t5 data_4
view t2_$v1 t2 3 $v1
del t2
view t3_$v1 t3 1 $v1
del t3
view t1_$v1 t1 1 $v1
del t1
view $o1_$v1 $o1 1 $v1
del $o1
ncon I1 2 t3_$v1 2 $o1_$v1 2
del t3_$v1
del $o1_$v1
ncon I2 1 t1_$v1 1,2 t4 2
del t1_$v1
del t4
ncon I3 1,3 t2_$v1 1,2,3 t5 2
del t2_$v1
del t5
ncon I4 2 I3 1,2 $o2 1
del I3
del $o2
ncon t8 2 I1 2 I2 2
del I1
del I2
ncon t9 0 t8 1 I4 1
del t8
del I4
save t9 output
```

### DSL Format and Instructions

The DSL file is a regular ASCII text file with one instruction per line. Lines beginning with `#` are comments which are
ignored (except for the first line which contains version information). The first line should have a version string
which specifies the format version. Comments following the first line contain metadata about the methods used
to determine the contraction plan used and which edges to slice. Symbols prefixed by a dollar sign indicate variables which are take different values
for each iteration of the computation. An explanation of each instruction follows

#### Outputs
```
outputs 2
```

This line indicates that the tensor network has two output tensors. This will create symbols for each of the possible values of
these output tensors. For a quantum circuit, measurement is performed the computational basis which means that the possible values for the output
tensors are the vectors `[1, 0]` and `[0, 1]` corresponding to the $| 0 \rangle$ and $| 1 \rangle$ states respectively. For two outputs the following
symbols and mappings will be defined

```
o1_0 => [1, 0]
o1_1 => [0, 1]
o2_0 => [1, 0]
o2_1 => [0, 1]
```

#### Load

Load instructions define a new tensor symbol using data from the input file. Here the first argument is the
name of the new symbol and the second argument is the key to find the data in the input [Data File](@ref).
Mulitple tensor symbols can have the same data associated with them. In the following example both `t4` and `t5`
will have the same data associated with them.

```
load t1 data_1
load t2 data_2
load t3 data_3
load t4 data_4
load t5 data_4
```

#### View

A view instruction creates a new tensor symbol by taking a view of an existing tensor. The first argument is the new symbol name, the second is the symbol of the tesnor to take the view on, the third argument is the dimension and the final value is the particular value for that dimension. An example view command is as follows

```
view t13_$v1 t13 3 $v1
```

The `$v1` variable will be replaced by a particular value before execution as described in the [Parameter File](@ref) section.

#### Delete

The `del` instructions marks a particular tensor for deletion and indicates that it will not be used for the rest of that contraction.

#### Contraction

The `ncon` instruction specificies a pairwise contraction of tensors. An example contraction command is as follows

```
ncon t15 2,3 $o1 1 t14_$v1_$v2 2,1,3
```

Which indicates that the `$o1` tensor should be contracted with the `t14_$v1_v2` tensor to get the `t15` tensor. The contraction should be performed over the only rank of the `$o1` tensor and the second rank of the `t14_$v1_$v2` tensor. Einstein summation notation is used for for specifying which indices to contract over. This convention uses repeated indices on the right hand side to indicate that those indices should be contracted over. In the above the `1` index appears twice on the right hand side which indicates that this index should be contracted over.
For the case where one of the tensors is a scalar, a `0` is used as a placeholder for the indices.
For example if `t1` is a scalar tensor and `t2` is a matrix,
the multiplication of the matrix by the scalar can be expressed as

```
ncon t3 1,2 t1 0 t2 1,2
```

Batched contractions (sometimes called hyper-contractions) are also supported. This is where the same index
is repeated on the right hand side and also appears on the left hand side.

#### Save

The save instruction indicates that the given tensor is an output from the contraction and provides a label for this. The
value of this tensor will then be used in a reduction operation and/or written to the output file.


## Data File

The data file is a [JLD](https://github.com/JuliaIO/JLD.jl) file and contains the numerical values of the initial tensors. Each tensor is stored as a multi-dimensional array with a data label that is referenced in `load` commands of the DSL file.
