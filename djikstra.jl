# load packages
using Pkg
Pkg.activate("metroflow")
Pkg.add("DataFrames")
Pkg.add("CSV")

using DataFrames
using CSV
using LightGraphs
using DataStructures

# load the data
grapharcs = CSV.read("metroarcs.csv",DataFrame)

# prepare the graph
nodes = unique(vcat(grapharcs.Origin,grapharcs.Destination))
nodeid = Dict([nodes[i] => i for i in eachindex(nodes)])
arcid = Dict((nodeid[grapharcs.Origin[i]],nodeid[grapharcs.Destination[i]])=> i for i in axes(grapharcs,1))

function djisktra(nodes,nodeid,grapharcs)
    println("Computing distances with Djikstras algorithm.")
    nr_nodes = length(nodes)
    distance = fill(Inf,nr_nodes,nr_nodes)
    previous_nodes = zeros(Int64,nr_nodes,nr_nodes)
    for node in nodes
        distance[nodeid[node],nodeid[node]] = 0.0
        queue = PriorityQueue{Int64, Float64}()
        enqueue!(queue, nodeid[node], 0.01)
        while !(isempty(queue))
            origin = dequeue!(queue)
            for arc in axes(grapharcs,1)
                if nodeid[grapharcs.Origin[arc]] == origin
                    destination = nodeid[grapharcs.Destination[arc]]
                    new_distance = distance[nodeid[node],origin] + grapharcs.Traveltime[arc]
                    if new_distance < distance[nodeid[node],destination]
                        distance[nodeid[node],destination] = new_distance
                        previous_nodes[nodeid[node],destination] = origin
                        if !(haskey(queue, destination))
                            enqueue!(queue,destination,new_distance)
                        else
                            delete!(queue, destination)
                            enqueue!(queue,destination,new_distance)
                        end
                    end
                end
            end
        end
    end
    return distance, previous_nodes
end

function create_graph(previous_nodes,grapharcs)
    nr_nodes = size(previous_nodes,1)
    nr_arcs = size(grapharcs,1)
    on_path = zeros(Bool,nr_nodes,nr_nodes,nr_arcs)
    for origin in 1:nr_nodes
        for destination in 1:nr_nodes
            previous = previous_nodes[origin,destination]
            current_destination = destination
            while previous != 0
                on_path[origin,destination,arcid[(previous,current_destination)]] = true
                current_destination = previous
                previous = previous_nodes[origin,current_destination]
            end
        end
    end
    return on_path
end

distance, previous_nodes = djisktra(nodes,nodeid,grapharcs)
on_path = create_graph(previous_nodes,grapharcs)




