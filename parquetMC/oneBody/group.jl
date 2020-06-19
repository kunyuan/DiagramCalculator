module OneBody

# include("../../parameter.jl")
include("../../utility/constant.jl")
using QuantumStatistics: Green, FastMath
using StaticArrays
const Mom = Main.Mom
include("./update.jl")

struct Group
    diag::Int
    order::Int
    reweight::Float
    innerK::Tuple{Int,Int}
    innerT::Tuple{Int,Int}
    statistics::Array{Float,2} # two dimension array of (k, tau)
    tauGrid::Any
    KGrid::Any
    function Group(diag, order, grids)
        if diag == Diag.∅
            @assert order == 0 "∅ diagram must be order 0!"
            tauGrid = nothing
            KGrid = nothing
        else
            @assert order > 0 "$diag diagram must be order ≥ 0!"
            tauGrid = grids.tau
            if (diag == Diag.Σ) || (diag == Diag.Δ)
                kGrid = grids.fermiK
            elseif diag == Diag.Π
                kGrid = grids.boseK
            else
                error("$diag diagram is not an implemented one-body observable!")
            end
        end

        return new(diag, 0, 1.0, (1, 1), (1, 1), zeros(Float, tauGrid.size, kGrid.size), tauGrid, kGrid)
    end
end

"""
two groups are equal if and only if the diagram and the order are the same
"""
Base.:=(x::Group, y::Group) = (x.diag == y.diag && x.order == y.order)

mutable struct State
    group::Group
    extTidx::Int
    extKidx::Int
    absweight::Float
    T::Array{Float,1}
    K::Array{Mom,1}
    propose::Dict{Symbol,Float} #
    accept::Dict{Symbol,Float}
    visited::Dict{Group,Float}
    # groups::Array{Group, 1}

    """
    MC Variable initialize here
    """
    function State(para, groups, updates)
        LastT = para.order + 2
        LastK = para.order + 1
        varT = rand(LastT) * para.β
        varK = [rand(Mom) * para.Kf for i = 1:LastK]

        curr = groups[1]
        varT[1] = 0.0
        varT[LastT] = curr.tauGrid[curr.extTidx]
        varK[1] = zero(Mom)
        varK[1][1] = curr.KGrid[curr.extKidx]

        state = new(0, curr, 1, 1, 0.0, LastT, varT, varK)
        state.absweight = abs(evaluate(para, state))
        return state
    end
end

"""
    init(parameter, grids)

build Diagram groups and the State from the para.diagrams and other parameters
"""
function init(para, grids)
    # maxOrder=maximum([d[2] for d in diagList])
    lastT = para.order + 1
    groups = [Group(Diag.∅, 0, grids), ]

    for o in 1:para.order
        for diag in para.diagrams
            push!(groups, Group(diag, o, grids))
        end
    end
    state = State(para, groups)
    return state, groups, (propGroup, rejectGroup)
end

function evaluate(para, state)
    group = state.group
    if group.diag == Diag.∅
        return 1.0
    elseif group.diag == Diag.Π
        if group.order == 1
            tau = state.T[state.lastT] - state.T[1]
            q, k = state.K[1], state.K[2]
            gweight = bareFermi(para.β, tau, k, para.μ) * bareFermi(para.β, tau, k + q, para.μ)
            return -para.spin * gweight / (2π)^para.dim
        else
            error("Not implemented")
        end
    else
        error("Not implemented")
    end
end

function measure(para, state)
    group = state.group
    if group.diag == Diag.∅
    else
    end
end

end