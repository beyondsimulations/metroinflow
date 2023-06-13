# load packages
using Pkg
Pkg.activate("metroflow")

using DataFrames
using CSV
using DataStructures
using JuMP
using HiGHS
using Plots
using Dates

# safety factor
safety = 0.8

# minutes in steps
minutes_in_step = 15

# prefered dispatch
prefered_dispatch = 100

# max_enter
max_enter = 6000

# max_change in entry (to previous period!)
max_change = 1000

# taxi during night
taxi_night = 750

# daterange
daterange = Date("2022-11-27"):Date("2022-11-30")

# load the data
grapharcs = CSV.read("data_metro/metroarcs.csv", DataFrame)

# prepare the graph
nodes = unique(vcat(grapharcs.origin,grapharcs.destination))
nodeid = Dict([nodes[i] => i for i in eachindex(nodes)])
arcid = Dict((nodeid[grapharcs.origin[i]],nodeid[grapharcs.destination[i]]) => i for i in axes(grapharcs,1))
nodes_clean = replace.(nodes,r"Metro_" => "")
nodes_clean = replace.(nodes_clean,r"SouthBound" => "SB")
nodes_clean = replace.(nodes_clean,r"NorthBound" => "NB")
nodes_clean = replace.(nodes_clean,r"EastBound"  => "EB")
nodes_clean = replace.(nodes_clean,r"WestBound"  => "WB")
nodes_clean = replace.(nodes_clean,r"Platform"  => "Plt")
nodes_clean = replace.(nodes_clean,r"Red"  => "R")
nodes_clean = replace.(nodes_clean,r"Green"  => "G")

function djisktra(
    nodes,
    nodeid,
    grapharcs
    )
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
    return distance, 
    previous_nodes
end

function create_graph(
    previous_nodes,
    grapharcs
    )
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

function new_arrivals!(
    nodeid,
    demand,
    step,
    new_queue
    )
    for movement in axes(demand,1)
        if demand.datetime[movement] == step
            new_queue[nodeid[demand.origin[movement]],nodeid[demand.destination[movement]]] += demand.value[movement]
        end
    end
end

function dispatch_queues!(
    station_moved,
    queue
    )
    total_queue = sum(queue,dims=2)
    for o in axes(queue,1)
        for d in axes(queue,2)
            if queue[o,d] > 0
                queue[o,d] = max(0,round(queue[o,d] - station_moved[o] * (queue[o,d]/total_queue[o])))
            end
        end
    end
end

function arc_utilization(
    grapharcs,
    station_moved,
    queue,
    on_path
    )
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

function metro_model(
    nodes,
    queue,
    grapharcs,
    on_path,
    minutes_in_step,
    safety,
    prefered_dispatch,
    max_enter,
    max_change,
    previous_entry
    )
    inflow_model = Model(HiGHS.Optimizer)
    set_attribute(inflow_model, "presolve", "on")
    set_attribute(inflow_model, "time_limit", 60.0)

    @variable(
        inflow_model, 
        0 .<= X[eachindex(nodes)].<= max_enter
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
        [n = eachindex(nodes); 0 < sum(queue[n,d] for d in eachindex(nodes))],
        X[n] >= prefered_dispatch
        )
    @constraint(
        inflow_model,
        [n = eachindex(nodes)],
        X[n] >= previous_entry[n] - max_change
        )
    @constraint(
        inflow_model,
        [n = eachindex(nodes)],
        X[n] <= previous_entry[n] + max_change
        )
    return inflow_model, X
end


function timestep_optimization(
    date,
    step,
    nodeid,
    nodes,
    demand,
    queue,
    grapharcs,
    on_path,
    minutes_in_step,
    safety,
    prefered_dispatch,
    max_enter,
    max_change,
    previous_entry
    )
    println("Starting timestep $step")
        new_arrivals!(nodeid,demand,step,queue)

        if !(DateTime(date)+Hour(3) <= step <= DateTime(date)+Hour(5))
            inflow_model, X = metro_model(nodes,queue,grapharcs,on_path,minutes_in_step,safety,prefered_dispatch,max_enter,max_change,previous_entry)
            JuMP.optimize!(inflow_model)
            solution_summary(inflow_model)
            station_moved = value.(X)
            arcs = arc_utilization(grapharcs,station_moved,queue,on_path)
            previous_entry .= station_moved
        else
            arcs = arc_utilization(grapharcs,zeros(Int64,length(nodes)),queue,on_path)
            station_moved = fill(taxi_night,length(nodes))
            previous_entry .= max_change
        end

        station_allowed = zeros(Int64,length(station_moved))
        for movement in eachindex(station_moved)
            if 0 < floor(station_moved[movement]) <= prefered_dispatch 
                station_allowed[movement] = prefered_dispatch
            else
                station_allowed[movement] = floor(station_moved[movement])
            end
        end

        station_queue = sum(queue,dims=2)

        dispatch_queues!(station_moved,queue)

        return station_allowed,
        station_moved,
        station_queue,
        arcs,
        previous_entry
    end

function metro_inflow(
    daterange,
    nodeid,
    nodes,
    grapharcs,
    safety,
    minutes_in_step,
    prefered_dispatch,
    max_enter,
    max_change)

    distance, previous_nodes = djisktra(nodes,nodeid,grapharcs)
    on_path = create_graph(previous_nodes,grapharcs)
    queue = zeros(Float64,length(nodeid),length(nodeid))
    previous_entry = zeros(Float64,length(nodeid),length(nodeid))
    results_queues = DataFrame(datetime = DateTime[],station=String[],allowed=Int64[],moved=Int64[],queued=Int64[])
    results_arcs   = DataFrame(datetime = DateTime[],connection=Int64[],line=String[],utilization=Float64[])

    for date in daterange

        demand = CSV.read("data_demand/OD_$(date)v3.csv", DataFrame)

        timesteps = DateTime(date):Minute(minutes_in_step):DateTime(date)+Day(1)-Minute(minutes_in_step)

        for step in timesteps
            
            station_allowed,
            station_moved,
            station_queue,
            arcs = timestep_optimization(
                date,
                step,
                nodeid,
                nodes,
                demand,
                queue,
                grapharcs,
                on_path,
                minutes_in_step,
                safety,
                prefered_dispatch,
                max_enter,
                max_change,
                previous_entry)

                
            for station in eachindex(nodes_clean)
                push!(results_queues,(
                    datetime = date,
                    station=nodes_clean[station],
                    allowed=floor(station_allowed[station]),
                    moved=floor(station_moved[station]),
                    queued=floor(station_queue[station])))
            end

            for arc in 1:size(grapharcs,1)
                push!(results_arcs,(
                    datetime = date,
                    connection=arc,
                    line=grapharcs.category[arc],
                    utilization=arcs.utilization[arc]))
            end

            plot_queue_demand = bar(station_queue,title="Queue at timestep $step",ylims=(0,35000),label="People in Queue",legend=:topright)
                bar!(station_allowed,label="People allowed to Enter",xticks=(1:length(nodes), nodes_clean),xrotation=30,xtickfontsize=4,titlefontsize=12)
            display(plot_queue_demand)

            plot_arc_demand = bar(arcs.utilization,title="Metroarcs at timestep $step",ylims=(0,1.1),label="Utilization",legend=:topright,titlefontsize=12)
                hline!([safety],label="Restriction")
            #display(plot_arc_demand)

        end
    end
    return results_queues, results_arcs
end


# Optimization
results_queues, results_arcs = metro_inflow(
    daterange,
    nodeid,
    nodes,
    grapharcs,
    safety,
    minutes_in_step,
    prefered_dispatch,
    max_enter,
    max_change
    )
# Simulation

