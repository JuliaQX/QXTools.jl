# User's Guide

QXTools generates output files which provide a description of the computations and data
required to perform a simulation. There are three output files:

- Parameter file: This YAML file provides information on the sampling method to use for the simulation
- DSL file: The DSL file contains instructions describing the tensor operations involved in
performing the simulation. Uses `.qx` suffix
- Data file: The data file contains the numerical values of the initial tensors of the tensor network

## Parallel Processing
There are two levels of parallelism that are used to distribute computations over different processing elements.
The highest level of parallism divides the set of output bitstrings among the available processing elements.
The lower level of parallelism breaks each contraction into a sum of smaller contractions by slicing bonds of the tensor network.
Slicing also reduces the total memory requirements and make it possible to perform contractions that would otherwise be intractable on a given system.
However for each additional edge that is sliced, the network must be contracted for all possible values of the edge, leading to exponential growth in the number of contractions with the number of edges that are sliced.
It is thus very important to choose the right edges to slice.
This is a non trivial optimisation problem, but there exist many methods to which can find reasonable sets of edges to slice.
In QXTools, the edges to slice are found when finding the contraciton plan which is handled by the [QXGraphDecompositions](https://github.com/juliaqx/QXGraphDecompositions.jl) package.
Once the edges to select have been identified these appear as `view` commands in the DSL file, explained further below.

## Parameter File
The parameter file parameters for the sampling method used to select output bitstrings.
The following example shows how an explicit list of 3 bitstrings can be provided.

```
output:
  method: List
  params:
    num_samples: 3
    bitstrings:
      - "010"
      - "101"
      - "111"
```

## DSL Specification

A Domain Specific Language (DSL) has been defined when enables better separation of concerns between components.
This DSL describes a contraction over a tensor network to calculate observables of interest in the form of a sequence of instructions which act on individual tensors.
This makes it possible to separate the development of the high performance distributed tensor network computation code from that of the higher level contraction planning, circuit and and network manipulation code.
The ".qx" suffix is used to identify these DSL files.
We will first show a simple example of such a DSL file and then descibe each instruction in detail.

### Example DSL file

An example of the DSL generated for the contraction of a two qubit GHZ circuit looks like.

```
# version: 0.4.0
# Determination of contraction plan:
#   Method used: "flow cutter"
#   Treewidth: 2
#   Time allocated: 2
#   Seed used: -1
#   Returned metadata:
#     1: "c min degree heuristic"
#     2: "c status 3 1627390755700"
#     3: "c min shortcut heuristic"
#     4: "c run with 0.0/0.1/0.2 min balance and node_min_expansion in endless loop with varying seed"
#   Hypergraph used: true
#   Hyperedge contraction method: "Netcon where possible, min fill heuristic otherwise."
# Slicing:
#   Method used: "greedy treewidth deletion"
#   Edges sliced: 2
#   Score fucntion used: direct_treewidth
#   Treewidths after slicing consecutive edges:
#     - 1
#     - 0
#
load t5 data_1 2
view t5_s t5 v2 1 2
load t1 data_2 2,2
view t1_s t1 v1 1 2
load t4 data_1 2
ncon I2 1 t1_s 1,2 t4 2
ncon t8 1,2 t5_s 1 I2 2
load t3 data_4 2
view t3_s t3 v1 1 2
ncon t9 1,2 t8 1,2 t3_s 2
output t6 1 2
view t6_s t6 v1 1 2
ncon t10 1,2 t9 1,2 t6_s 2
load t2 data_3 2,2,2
view t2_s t2 v1 3 2
view t2_s_s t2_s v2 2 2
output t7 2 2
ncon I1 2,3 t2_s_s 1,2,3 t7 1
ncon t11 0 t10 1,2 I1 1,2
save output t11
```

### DSL Format and Instructions

The DSL file is a regular ASCII text file with one instruction per line.
Lines beginning with `#` are comments which are
ignored (except for the first line which contains version information).
The first line has a version string which specifies the format version.
Comments following the first line contain metadata about the methods used to determine the contraction plan used and which edges to slice.
We now go through each of the instructions in order of appearence.

#### Load

Load instructions define a new tensor symbol using data from the input file.
Here the first argument is the name of the new symbol, the second argument is the key to find the data at and the third argument is a comma separated list of dimensions of the tensor.

#### View

The view instruction creates a new tensor symbol by taking a view of an existing tensor.
The first argument is the new symbol name, the second is the symbol of the tensor to take the view of, the third argument is the symbol used to identify the index, the fourth identifies the rank of the index of the tensor to slice and the final argument is the dimension of the bond.
For example the instruction

```
view t5_s t5 v2 1 2
```

will create a new tensor labeled `t5_s` by taking a view on the tensor labeled `t5` by setting the value of the first index to the value labeled by symbol `v2`.
If `t5_s`, `t5` and `v2` were variables this would look something like

```
t5_s = t5[v2]
```

#### Contraction

The `ncon` instruction specifies a pairwise contraction of tensors.
An example contraction command is as follows

```
ncon I2 1 t1_s 1,2 t4 2
```

Which indicates that the `t1_s` tensor should be contracted with the `t4` tensor to get the `I2` tensor.
The contraction should be performed over the second rank of the `t1_s` tensor and first rank of the `t4` tensor.
Einstein summation notation is used for for specifying which indices to contract over.
This convention uses repeated indices on the right hand side to indicate that those indices should be contracted over.
In the above, the `2` index appears twice on the right hand side which indicates that this index should be contracted over.
For the case where one of the tensors is a scalar, a `0` is used as a placeholder.
For example if `t1` is a scalar tensor and `t2` is a matrix,
the multiplication of the matrix by the scalar can be expressed as

```
ncon t3 1,2 t1 0 t2 1,2
```

Batched contractions (sometimes called hyper-contractions) are also supported. This is where the same index is repeated on the right hand side and also appears on the left hand side.


#### Output

The output instruction marks tensors as corresponding to circuit outputs and provides their dimension.
For example, the instruction

```
output t6 1 2
```

specifies that the tensor `t6` will refer to the first output and that it has dimension 2.


#### Save

The save instruction indicates that the given tensor is an output from the contraction and provides a label for this.
The value of this tensor will then be used in a reduction operation and/or written to the output file.


## Data File

The data file is a [JLD2](https://github.com/JuliaIO/JLD2.jl) file and contains the numerical values of the initial tensors.
Each tensor is stored as a multi-dimensional array with a data label that is referenced in `load` commands of the DSL file.
