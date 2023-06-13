# load packages
using Pkg
Pkg.activate("metroflow")

using DataFrames
using CSV
using Dates

for dates in Date("2022-11-27"):Date("2022-11-30")
    demand = CSV.read("data_demand/OD_$(dates)v2.csv", DataFrame)
    for add_minutes in axes(demand,1)
        if demand.origin[add_minutes] in ["Metro_MsheirebGoldPlatformEastBound","Metro_MsheirebRedAndGreenPlatformSouthBound","Metro_MsheirebGoldPlatformWestBound","Metro_MsheirebRedAndGreenPlatformNorthBound"]
            demand.origin[add_minutes] = "Metro_MsheirebConcourse"
        end
        if demand.origin[add_minutes] in ["Metro_AlBiddaGreenPlatformNorthBound","Metro_AlBiddaRedAndGreenPlatformNorthBound","Metro_AlBiddaRedPlatformNorthBound","Metro_AlBiddaRedPlatformSouthBound","Metro_AlBiddaGreenPlatformSouthBound","Metro_AlBiddaRedAndGreenPlatformSouthBound"]
            demand.origin[add_minutes] = "Metro_AlBiddaConcourse"
        end
    end
    select!(demand,[:date,:datetime,:origin,:destination,:value])
    CSV.write("data_demand/OD_$(dates)v3.csv", demand)
end

