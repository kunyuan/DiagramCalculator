module Markov
include("parameter.jl")
include("grid.jl")
include("observable.jl")
include("utility/utility.jl")
include("diag/propagator.jl")
import .Propagator: green
using Random, StaticArrays, Printf, Dates
import .Vertex4, .Ver4Test, .Polar
const UpdateNum = 6
const INCREASE_ORDER, DECREASE_ORDER, CHANGE_EXTTAU, CHANGE_EXTK, CHANGE_TAU, CHANGE_K =
    1:UpdateNum
const Name = [
    "increase_order",
    "decrease_order",
    "change_ExtTau",
    "change_ExtK",
    "change_Tau",
    "change_K",
]
const Accepted = (@MArray zeros(Float, UpdateNum, Order + 1)) .+ 1.0e-10
const Proposed = (@MArray zeros(Float, UpdateNum, Order + 1)) .+ 1.0e-10

const curr = Main.Curr
const rng = Main.Curr.rng
const varK = Main.Curr.K
const varT = Main.Curr.T

const ver4 = Vector{Vertex4.Ver4}(undef, 0)
const oneBody = Observable.OneBody()
const polar = Vector{Polar.Polarization}(undef, 0)

# function init(_counter, _rng)
function init()
    #######  initialize MC variables  ################################
    # println(oneBody.norm)
    # println(typeof(oneBody))
    ###### initialized diagram trees #######################################

    curr.absWeight = Markov.eval(curr.order)
end

@fastmath function test() end

function eval(order)
    if order == 0
        return 1.0
    else
        Tau = varT[LastTidx] - varT[1]
        @assert 0.0 < Tau < Beta "t out of range"
        @assert abs(Grid.K.grid[curr.extKidx] - curr.K[1][1]) < 1.0e-15 "ExtK doesn't match!"
        @assert abs(Grid.tau.grid[curr.extTidx] - curr.T[LastTidx]) < 1.0e-15 "ExtT doesn't match! $(Grid.tau.grid[curr.extTidx]) vs $(curr.T[LastTidx])"
        @assert 0 < curr.extKidx <= KGridSize "K out of range"
        @assert 0 < curr.extTidx <= TauGridSize "Tau out of range"
        gWeight = green(Tau, varK[2]) * green(-Tau, varK[2] + varK[1])
        return -SPIN * gWeight * PhaseFactor
    end
    # return 
    # order == 0 && return 1.0
    # if DiagType == POLAR
    #     return Polar.eval(polar[order])
    # elseif DiagType == GAMMA
    #     Vertex4.eval(ver4[order], varK[INL], varK[OUTL], varK[INR], varK[OUTR], 5, true)
    #     chanWeight = sum(ver4[order].weight)
    #     return chanWeight[DI] + chanWeight[EX] / SPIN
    # else
    #     throw("Not implemented!")
    # end
end

function measure()

    factor = 1.0 / curr.absWeight / ReWeight[curr.order+1]
    # println("$(curr.order), $(curr.absWeight), $(ReWeight[curr.order + 1])")

    if curr.order == 1
        Tau = curr.T[LastTidx]
        weight =
            -green(Tau, curr.K[2]) * green(-Tau, curr.K[2] + curr.K[1]) * SPIN * PhaseFactor
        @assert abs(weight - eval(curr.order)) < 1.0e-10 "weight wrong"
    end

    @assert isinf(factor) == false "factor is infinite at step $(curr.step)"
    @assert abs(curr.absWeight - abs(eval(curr.order))) < 1.0e-10 "weight is wrong!"

    if DiagType == POLAR || DiagType == SIGMA || DiagType == DELTA
        Observable.measure(oneBody, eval(curr.order), factor)
    end
end

function save()
    println("Saving data to disk ...")
    if DiagType == POLAR || DiagType == SIGMA || DiagType == DELTA
        Observable.save(oneBody)
    end
end

function reweight()
    return
end

@inline propose(name) = (Proposed[name, curr.order+1] += 1)
@inline accept(name) = (Accepted[name, curr.order+1] += 1)

function increaseOrder()
    curr.order == Order && return # already at the highest order
    newOrder = curr.order + 1
    prop = 1.0
    if curr.order == 0
        # create new external Tau
        # curr.extTidx, propT = createExtIdx(TauGridSize)
        # varT[LastTidx] = Grid.tau.grid[curr.extTidx]

        # curr.extKidx, propK = createExtIdx(KGridSize)
        # varK[1][1] = Grid.K.grid[curr.extKidx]
        # prop = propT * propK
    else
        # create new internal Tau
        varT[lastInnerTidx(newOrder)], prop = createTau()
    end
    # newOrder == 1 && println(lastInnerKidx(newOrder))
    # oldK = copy(varK[2])
    # prop *= createK!(varK[lastInnerKidx(newOrder)])
    prop *= createK!(varK[2])
    # @assert norm(oldK) != norm(varK[2]) "K remains the same"

    newAbsWeight = abs(eval(newOrder))
    # println(prop, ", ", newAbsWeight)
    R = prop * newAbsWeight * ReWeight[newOrder+1] / curr.absWeight / ReWeight[curr.order+1]
    propose(INCREASE_ORDER)
    if rand(rng) < R
        accept(INCREASE_ORDER)
        curr.order = newOrder
        curr.absWeight = newAbsWeight
    end
end

function decreaseOrder()
    curr.order == 0 && return
    newOrder = curr.order - 1
    prop = 1.0
    if newOrder == 0
        # remove external Tau 
        # prop *= removeExtIdx(TauGridSize)
        # prop *= removeExtIdx(KGridSize)
    else
        # remove internal Tau
        prop *= removeTau()
    end
    # prop *= removeK(varK[lastInnerKidx(curr.order)])
    prop *= removeK(varK[2])
    newAbsWeight = abs(eval(newOrder))
    R = prop * newAbsWeight * ReWeight[newOrder+1] / curr.absWeight / ReWeight[curr.order+1]
    propose(DECREASE_ORDER)
    if rand(rng) < R
        accept(DECREASE_ORDER)
        curr.order = newOrder
        curr.absWeight = newAbsWeight
    end
end

function changeTau()
    # Proposed[CHANGE_TAU, curr.order + 1] += 1
    # Accepted[CHANGE_TAU, curr.order + 1] += 1
    return
end

function changeK()
    # return
    curr.order == 0 && return
    loopidx = rand(rng, firstInnerKidx():lastInnerKidx(curr.order))
    oldK = copy(varK[loopidx])
    prop = shiftK!(oldK, varK[loopidx])

    propose(CHANGE_K)
    newAbsWeight = abs(eval(curr.order))
    if rand(rng) < prop * newAbsWeight / curr.absWeight
        accept(CHANGE_K)
        curr.absWeight = newAbsWeight
    else
        varK[loopidx] .= oldK
    end
end

function changeExtTau()
    # curr.order == 0 && return
    oldTidx = curr.extTidx
    curr.extTidx, prop = shiftExtIdx(TauGridSize)
    varT[LastTidx] = Grid.tau.grid[curr.extTidx]
    propose(CHANGE_EXTTAU)
    newAbsWeight = abs(eval(curr.order))
    if rand(rng) < prop * newAbsWeight / curr.absWeight
        accept(CHANGE_EXTTAU)
        curr.absWeight = newAbsWeight
    else
        curr.extTidx = oldTidx
        varT[LastTidx] = Grid.tau.grid[curr.extTidx]
    end
end

function changeExtK()
    # curr.order == 0 && return
    oldKidx = curr.extKidx
    prop = 1.0
    if DiagType == POLAR
        curr.extKidx, prop = shiftExtIdx(KGridSize)
        varK[1][1] = Grid.K.grid[curr.extKidx]
    else
        return
    end
    propose(CHANGE_EXTK)
    newAbsWeight = abs(eval(curr.order))
    if rand(rng) < prop * newAbsWeight / curr.absWeight
        accept(CHANGE_EXTK)
        curr.absWeight = newAbsWeight
    else
        curr.extKidx = oldKidx
        varK[1][1] = Grid.K.grid[curr.extKidx]
    end
end


const barbar = "====================================================================================="
const bar = "-------------------------------------------------------------------------------------"

function printStatus()
    # Var.counter += 1
    println(barbar)
    printstyled(Dates.now(), color = :green)
    println("\nStep:", curr.step)
    println(bar)
    for i = 1:UpdateNum
        @printf("%-14s %12s %12s %12s\n", Name[i], "Proposed", "Accepted", "Ratio  ")
        for o = 0:Order
            @printf(
                "  Order%2d:     %12.0f %12.0f %12.6f\n",
                o,
                Proposed[i, o+1],
                Accepted[i, o+1],
                Accepted[i, o+1] / Proposed[i, o+1]
            )
        end
        println(bar)
    end
    println(progressBar(round(curr.step / 1000_000, digits = 2), TotalBlock))
    println()
end

end
