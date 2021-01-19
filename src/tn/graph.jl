using LightGraphs

convert_to_graph(tnc::TensorNetworkCircuit) = convert_to_graph(tnc.tn)

function convert_to_graph(tn::TensorNetwork)
    g = SimpleGraph(length(tn.data))
    # TODO: add structures to tensor network to avoid N^2 complexity here
    for (i, itensor) in enumerate(tn.data[1:end-1])
        for (j,jtensor) in enumerate(tn.data[i+1:end])
            if length(intersect(inds(itensor), inds(jtensor))) > 0
                add_edge!(g, i, i+j)
            end
        end
    end
    g
end