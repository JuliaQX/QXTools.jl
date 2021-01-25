import LightGraphs
import ITensors

export convert_to_graph

convert_to_graph(tnc::TensorNetworkCircuit) = convert_to_graph(tnc.tn)

function convert_to_graph(tn::TensorNetwork)
    g = LightGraphs.SimpleGraph(length(tn))
    tensor_map = Dict{Symbol, Int64}([k => i for (i, k) in enumerate(keys(tn))])
    # TODO: add structures to tensor network to avoid N^2 complexity here
    for bond in bonds(tn)
        tensors = tn[bond]
        if length(tensors) == 2
            LightGraphs.add_edge!(g, tensor_map[tensors[1]], tensor_map[tensors[2]])
        end
    end
    g
end