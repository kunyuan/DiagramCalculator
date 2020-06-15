
module Observable
include("parameter.jl")
include("grid.jl")
using StaticArrays: MArray
using JLD2, FileIO

const curr = Main.Curr

mutable struct OneBody
    norm::Float # normalization 
    phy::Float # the physcial weight of the normalization diagrams
    # static::MArray{(Order, KGridSize),Float}
    data::Array{Float,3} # order, kgrid, taugrid
    OneBody() = new(1.0e-10, KGridSize, zeros(Float, TauGridSize, KGridSize, Order))
end

function measure(obs::OneBody, weight, factor)
    # @assert isapprox(curr.T[LastTidx], Grid.tau.grid[curr.extTidx]) "Not good"
    if curr.order == 0
        obs.norm += weight * factor
    else
        @assert abs(Grid.K.grid[curr.extKidx] - curr.K[1][1]) < 1.0e-15 "ExtK doesn't match!"
        @assert abs(Grid.tau.grid[curr.extTidx] - curr.T[LastTidx]) < 1.0e-15 "ExtT doesn't match! $(Grid.tau.grid[curr.extTidx]) vs $(curr.T[LastTidx])"
        @assert 0 < curr.extKidx <= KGridSize "K out of range"
        @assert 0 < curr.extTidx <= TauGridSize "Tau out of range"
        obs.data[curr.extTidx, curr.extKidx, curr.order] += weight * factor
    end
end

function save(obs::OneBody)
    filename = "$(name())_pid$(curr.PID).jld2"
    data =
        Dict("PID" => curr.PID, "Norm" => obs.norm, "Data" => obs.data / obs.norm * obs.phy)

    for ki = 1:KGridSize
        println(
            "k=$(Grid.K.grid[ki]): ",
            sum(obs.data[:, ki, 1]) / length(TauGridSize) * Beta / obs.norm * obs.phy,
        )
    end

    # println(obs.data[1,1,1], "norm:", obs.norm)

    # FileIO.save(filename, data, compress = true)
    FileIO.save(filename, data)
end

function name()
    if DiagType == POLAR
        return "polar"
    elseif DiagType == SIGMA
        return "sigma"
    elseif DiagType == DELTA
        return "delta"
    elseif DiagType == GAMMA
        return "gamma"
    else
        throw("Not implemented!")
    end
end

end
