# Introduction

In this tutorial, we will cover the basics of creating a simple quantum circuit, using only primitive operations defined by `QXTools`, and contracting its' associated tensor network to compute a probability amplitude for a bitstrings we might see when the qubits are measured at the end of the circuit. We begin by importing the `QXTools` package:

```
using QXTools
```

## Creating a TensorNetworkCircuit

`QXTools` uses objects called `TensorNetworkCircuit` to represent quantum circuits. A `TensorNetworkCircuit` contains a tensor network, representing the output state of the the circuit's qubits, as one of its field. We can create a `TensorNetworkCircuit` for an empty quantum circuit by calling the `TensorNetworkCircuit` constructor with the number of qubits we want in the circuit. Below we create an empty 3 qubit circuit with no gates.

```
tnc = TensorNetworkCircuit(3)
```

We can add a single qubit gate to our circuit by passing the matrix representation of the gate to the `Base.push!` function. We also need to specify which qubit we want the gate to act on as shown below.  Note, the qubits in the circuit are numbered 1 to N where N is the number of qubits in the circuit.

```
X = [[0., 1.] [1., 0.]]

# Apply an X gate to both the first and third qubits in the circuit.
push!(tnc, [1], X)
push!(tnc, [3], X)
```

Similarly, to add a two qubit gate we can pass the matrix representation of the gate to the `push!` function, along with an array of integers indicating the target qubits we want the gate to act on. Note, for controlled gates, the first qubit in the array indicates the target qubit while the second qubit is taken as the control qubit.

```
# Matrix representation of the controlled not gate.
CX = [[1., 0., 0., 0.] [0., 1., 0., 0.] [0., 0., 0., 1.] [0., 0., 1., 0.]]

# Note, the order of the qubits here is [target, control].
push!(tnc, [1, 2], CX)
```


## Contracting a TensorNetworkCircuit.

After creating a `TensorNetworkCircuit`, we may then want to compute the probability amplitude for a particular bitstring. To compute amplitudes for different measurement outcomes, we first set the initial state of the qubits in the circuit. This can be done using the `add_inputs!` function. For example, calling `add_inputs!(tnc, "111")` will add tensors to the tensor network which set the initial state of the circuit's qubits to the basis state labeled by "111". By default, calling `add_inputs!(tnc)` will set the initial qubit state to the "all zeros" state. 

```
# Set the initial state of the qubits to "000".
add_input!(tnc)
```

Next, we need to calculate a contraction plan for the tensor network. A contraction plan for a `TensorNetworkCircuit` determines the order in which tensors in the network are contracted to compute a probability amplitude. With `QXTools`, we can compute a contraction plan using an algorithm known as FlowCutter for a set amount of time. Below, we run FlowCutter 1 second to find a contraction plan for our quantum circuit.

```
# Find a good contraction plan.
plan = flow_cutter_contraction_plan(tnc; time=1)
```

We can now compute probabiliy amplitudes using the `single_amplitude` function. Calling this function on our `TensorNetworkCircuit` object, its contraction plan and a string, indicating a product state of the qubits, will result in a copy of the circuit's tensor network being contracted according to the given contraction plan and the desired probability amplitude being returned. Note, the characters in the string give the states of the individual qubits in ascending order such that the first character in the string gives the state of the first qubit in the circuit etc. 

```
# Contract the network using this plan to find the given amplitude for different outputs.
@show single_amplitude(tnc, plan, "000")
@show single_amplitude(tnc, plan, "111")

# The `+` and `-` characters are also supported to represent the states ("0" + "1") and ("0" - "1") respectively.
@show single_amplitude(tnc, plan, "1++")
```