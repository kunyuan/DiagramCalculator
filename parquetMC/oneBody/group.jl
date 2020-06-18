module OneBody

# include("../../parameter.jl")
include("../../utility/constant.jl")
using QuantumStatistics: Green, FastMath
using StaticArrays
const Mom=MVector{3, Float}

struct Group
    diag::Symbol
    order::Int
    innerK::Tuple{Int, Int}
    innerT::Tuple{Int, Int}
    statistics::Array{Float, 2} #two dimension array of (k, tau)
    tauGrid::Any
    KGrid::Any
    function Group(diag, order, grids)
        if diag==:∅
            @assert order==0 "∅ diagram must be order 0!"
            return new(diag, 0, (1, 1), (1, 1), zeros(1,1).+1.0e-10, nothing, nothing)
        else
            @assert order>0 "$diag diagram must be order ≥ 0!"
            tauGrid=grids.tau
            if diag==:Σ || diag==:Δ
                kGrid=grids.fermiK
                return new(diag, 0, (1, 1), (1, 1), zeros(Float, tauGrid.size, kGrid.size), tauGrid, kGrid)
            elseif diag==:Π
                kGrid=grids.boseK
                return new(diag, 0, (1, 1), (1, 1), zeros(Float, tauGrid.size, kGrid.size), tauGrid, kGrid)
            else
                error("$diag diagram is not an implemented one-body observable!")
            end
        end
    end
end

mutable struct State
    step::Int
    group::Group
    extTidx::Int
    extKidx::Int
    lastT::Int
    absWeight::Float
    T::Array{Float, 1}
    K::Array{Mom, 1}
    # groups::Array{Group, 1}

    """
    MC Variable initialize here
    """
    function State(para, group)
        LastT=para.order+2
        LastK=para.order+1
        β, Kf=para.β, para.Kf
        varT = [rand() .* β for i = 1:LastT]
        varK = [rand(Mom) .* Kf for i = 1:LastK]

        curr = new(0, group, 1, 1, LastT, 0.0, varT, varK)
        curr.T[1] = 0.0
        curr.K[1] .= zero(Mom)
        if group.diag!=:∅
            curr.T[LastT] = group.tauGrid[curr.extTidx]
            curr.K[1][1] = group.KGrid[curr.extKidx]
        end

        curr.absWeight=abs(evaluate(para, curr))
        return curr
    end
end

"""
    init(parameter, grids)

build Diagram groups and the State from the para.diagrams and other parameters
"""
function init(para, grids)
    # maxOrder=maximum([d[2] for d in diagList])
    lastT=para.order+1
    groups=[Group(:∅, 0, grids), ]

    for o in 1:para.order
        for diag in para.diagrams
            push!(groups, Group(diag, o, grids))
        end
    end
    curr=State(para, groups[1])
    return curr, groups
end

# function measure(group, factor)
    # group.data[0]+=1.0*factor;
# end

# struct Polar<:OneBody
#     function Polar(order, extTidx)
#     end
# end

function evaluate(para, state)
    group=state.group
    if group.diag==:∅
        return 1.0
    elseif group.diag==:Π && group.order==1
        tau=state.T[state.lastT]-state.T[1]
        q, k=state.K[1], state.K[2]
        gweight=bareFermi(para.β, tau, k, para.μ)*bareFermi(para.β, tau, k+q, para.μ)
        return -para.spin*gweight/(2π)^para.dim
    else
        error("Note implemented")
    end
end

include("update.jl")

end