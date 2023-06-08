# load packages
using Pkg
Pkg.activate("metroflow")

using DataFrames
using CSV
using DataStructures
using JuMP
using HiGHS
using Plots

# safety factor
safety = 0.9
# minutes in steps
minutes_in_step = 60

# load the data
grapharcs = CSV.read("data_metro/metroarcs.csv", DataFrame)
demand = CSV.read("data_demand/OD_2022-11-30_fullhour.csv", DataFrame)

# prepare the graph
nodes = unique(vcat(grapharcs.origin,grapharcs.destination))
nodeid = Dict([nodes[i] => i for i in eachindex(nodes)])
arcid = Dict((nodeid[grapharcs.origin[i]],nodeid[grapharcs.destination[i]]) => i for i in axes(grapharcs,1))

# prepare the demand
timesteps = sort(unique(demand.dayhour))

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
                if nodeid[grapharcs.origin[arc]] == origin
                    destination = nodeid[grapharcs.destination[arc]]
                    new_distance = distance[nodeid[node],origin] + grapharcs.traveltime[arc]
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

function new_arrivals!(nodeid,demand,step,new_queue)
    for movement in axes(demand,1)
        if demand.dayhour[movement] == step
            new_queue[nodeid[demand.origin[movement]],nodeid[demand.destination[movement]]] += demand.value[movement]
        end
    end
end

function dispatch_queues!(station_moved,queue)
    total_queue = sum(queue,dims=2)
    for o in axes(queue,1)
        for d in axes(queue,2)
            if queue[o,d] > 0
                queue[o,d] = max(0,round(queue[o,d] - station_moved[o] * (queue[o,d]/total_queue[o])))
            end
        end
    end
end

function arc_utilization(grapharcs,station_moved)
    arcs = copy(grapharcs)
    total_queue = sum(queue,dims=2)
    arcs.utilization .= 0.0
    for a in 1:size(grapharcs,1)
        for o in eachindex(station_moved)
            for d in eachindex(station_moved)
                if queue[o,d] > 0
                    if on_path[o,d,a]
                        arcs.utilization[a] += station_moved[o] * (queue[o,d]/total_queue[o])
                    end
                end
            end
        end
    end
    arcs.utilization .= arcs.utilization ./ (arcs.capacity * minutes_in_step)
    return arcs
end


# Optimization

distance, previous_nodes = djisktra(nodes,nodeid,grapharcs)
on_path = create_graph(previous_nodes,grapharcs)
queue = zeros(Float64,length(nodeid),length(nodeid))
allowed_entry = zeros(Int64,length(timesteps),length(nodeid))

for step in timesteps
    println("Starting timestep $step")
    new_arrivals!(nodeid,demand,step,queue)

    if !(3 <= step <= 5)
        inflow_model = Model(HiGHS.Optimizer)
        set_attribute(inflow_model, "presolve", "on")
        set_attribute(inflow_model, "time_limit", 60.0)

        @variable(
            inflow_model, 
            X[eachindex(nodes)].>=0
            )
        @objective(
            inflow_model, Min,
            sum((sum(queue[o,d] for d in eachindex(nodes)) - X[o] for o in eachindex(nodes)).^2)
            )
        @constraint(
            inflow_model,
            [a = 1:size(grapharcs,1)],
            sum(X[o]*(queue[o,d]/(sum(queue[o,e] for e in eachindex(nodes))+1)) for o in eachindex(nodes), d in eachindex(nodes) if on_path[o,d,a] == true) <= grapharcs.capacity[a] * minutes_in_step * safety
            )
        @constraint(
            inflow_model,
            [n = eachindex(nodes)],
            X[n] <= sum(queue[n,d] for d in eachindex(nodes))
            )

        JuMP.optimize!(inflow_model)
        solution_summary(inflow_model)
        station_moved = value.(X)
    else
        station_moved = zeros(Float64,length(nodes))
    end

    station_queue = sum(queue,dims=2)

    dispatch_queues!(station_moved,queue)

    arcs = arc_utilization(grapharcs,station_moved)

    display(bar(station_queue.-station_moved,title="Queue at timestep $step",ylims=(0,25000)))
    #display(bar(station_moved,title="Allowed to enter at timestep $step",ylims=(0,25000)))
    #display(bar(arcs.utilization,title="Arc at timestep $step",ylims=(0,1)))

end

# Simulation

