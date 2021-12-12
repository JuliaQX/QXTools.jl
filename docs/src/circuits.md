# Circuits with QXZoo and YaoQX

In the [introduction](basics.md) tutorial circuits are constructed by explicitly defining each of the gate matrices.
For larger circuits this can get tedious and it is easier to construct circuits using QXZoo or to convert circuits created using [Yao.jl](https://yaoquantum.org/) to TensorNetworkCircuits.

## QXZoo
QXZoo was developed to simplify this process and provides functions for constructing commonly used circuits as well as primatives for building circuits from commonly used gate sets.
QXTools wraps many of the functions from QXZoo to improve integration and provides a simple interface for creating [Greenberger–Horne–Zeilinger](https://arxiv.org/abs/0712.0921) (GHZ), Quantum Fourier Transform (QFT) and Random Quantum Circuits (RQC).

To create a GHZ circuit with 5 qubits one would use

```
using QXTools
circ = create_ghz_circuit(5)
```
This returns a QXZoo Circ struct which contains `num_qubits` and `circ_ops` data members giving the number of qubits and access to the gates respectively.
QXTools provides `gate_qubits` and `gate_matrix` functions to simplify extracting the qubits that the gate acts upon and the matrix representation of the gate.

```
println("Circuit has $(circ.num_qubits) qubits")
for g in circ.circ_ops
    println("Gate acts on qubits $(gate_qubits(g))")
    println("Matrix form $(gate_matrix(g))")
end
```

QXTools provides the `convert_to_tnc` function which will convert a QXZoo circuit to a TensorNetworkCircuit

```
tnc = convert_to_tnc(circ)
```

The `create_qft_circuit` will create a QFT circuit when given the number of qubits.

```
circ = create_qft_circuit(5)
```

Finally the `create_rqc_circuit` creates a RQC circuit and takes arguments specifying the number of rows, columns, depth and the seed to use when sampling unitary gates. One can also specify whether or not to add a final layer of Hadamard which is sometimes ommitted.

```
circ = create_rqc_circuit(4, 4, 16, 42, final_h=true)
```

## YaoQX

The [YaoQX](https://github.com/JuliaQX/YaoQX.jl) package provides a lightweight wrapper which enables circuits constructed using the [Yao](https://yaoquantum.org/) tools to be used with QXTools.
[YaoQX](https://github.com/JuliaQX/YaoQX.jl) can be installed from the Julia package registry with

```
] add YaoQX
```

We also install [YaoBlocks](https://github.com/QuantumBFS/YaoBlocks.jl) to allow us to construct Yao circuits.

```
] add YaoBlocks
```

The following is an example of creating a simple  3 qubit GHZ circuit with YaoBlocks and converting it to a QXTools TensorNetworkCircuit.

```
using YaoQX
using YaoBlocks
n = 3
circ = chain(put(1=>H), chain(map(x -> control(n, x, x+1=>X), 1:n-1)...))
```

YaoQX implements the `convert_to_tnc` function for YaoBlocks circuits so these circuits can be converted to TensorNetworkCircuits in exactly the same way that QXZoo circuits can be with

```
tnc = convert_to_tnc(circ)
```

Note that the same arguments are supported so to convert to a TensorNetworkCircuit with no input or output one can add additional keyword arguments as

```
tnc = convert_to_tnc(circ, no_input=true, no_output=true)
```
