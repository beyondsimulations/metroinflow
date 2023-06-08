# load packages
using Pkg
Pkg.activate("metroflow")

using DataFrames
using CSV
using Dates

for dates in Date("2022-11-27"):Date("2022-11-30")
    demand = CSV.read("data_demand/OD_$(dates)_quarterhour.csv", DataFrame)
    demand.date .= dates
    demand.datetime .= DateTime(dates)
    for add_minutes in axes(demand,1)
        demand.datetime[add_minutes] += Minute(demand.quarterhour[add_minutes] * 15)
    end
    select!(demand,[:date,:datetime,:origin,:destination,:value])
    CSV.write("data_demand/OD_$(dates)", demand)
end

