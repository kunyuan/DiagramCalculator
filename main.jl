include("parameter.jl")
include("grid.jl")
# include("utility/utility.jl")

using Random
using StaticArrays

PID=length(ARGS) >0 ? PID = parse(Int, ARGS[1]) : rand(100_000:999_999)
Random.seed!(PID) #initialize rand(), which is already global and by default using MersenneTwister
@assert PID>=0 "RNG seed must be >=0"

macro flush(content)
    return :( $content )
end

@flush printstyled(
    "Rs: $Rs, kF: $Kf, EF: $Ef, β: $β, T/T_F: $(round(1.0 / β / Ef, digits = 4))\nSeed: $PID, TotalBlock: $TotalBlock\n",
    color = :green,
)


mutable struct State
    PID::Int
    step::Int
    order::Int
    scaleidx::Int
    extTidx::Int
    extKidx::Int
    extAngidx::Int
    absWeight::Float
    T::MVector{LastTidx,Float}
    K::SVector{LastKidx,Mom}

    function State(pid)
        varT = @MVector [rand() .* β for i = 1:LastTidx]
        varK = @SVector [rand(Mom) .* Kf for i = 1:LastKidx]
        # varT = rand(rng, LastTidx) .* Beta
        # varK = rand(rng, Mom, LastKidx) .* Kf
        # rand!(rng, varK)
        curr = new(pid, 0, 0, 1, 1, 1, 1, 0.0, varT, varK)
        curr.T[LastTidx] = Grid.tau.grid[curr.extTidx]
        curr.T[1] = 0.0

        if DiagType == GAMMA
            kL, kR = zero(Mom), zero(Mom)
            kL[1] = Kf
            curr.K[OUTL] .= curr.K[INL] .= kL
            θ = acos(Grid.angle.grid[curr.extAngidx])
            kR[1:2] .= [cos(θ), sin(θ)] .* Kf
            curr.K[OUTR] .= curr.K[INR] .= kR
        else
            curr.K[1] .= zero(Mom)
            curr.K[1][1] = Grid.K.grid[curr.extKidx]
        end
        return curr
    end
end

const Curr = State(PID)

# include("markov.jl")
# Markov.init()

################# Test and Benchmark ############################
# Markov.test()

# import BenchmarkTools: @btime
# for _order = 1:Order
#     println("Benchmark Order $_order")
#     @btime Markov.eval(o) samples = 1 evals = 100 setup = (o = $_order)
    # println(sum(Markov.ver4[_order].weight))
# end
# exit()

# println("Start Simulation ...")
# block = 0

# printTimer = StopWatch(PrintTime, Markov.printStatus)
# saveTimer = StopWatch(SaveTime, Markov.save)
# reweightTimer = StopWatch(ReWeightTime, Markov.reweight)
# messageTimer = StopWatch(MessageTime, Markov.save)

# for block = 1:TotalBlock
#     for i = 1:1000_000
#         Curr.step += 1
#         x = rand(Curr.rng)
#         if x < 1.0 / 6.0
#             Markov.increaseOrder()
#         elseif x < 2.0 / 6.0
#             Markov.decreaseOrder()
#         elseif x < 3.0 / 6.0
#             # Markov.changeK()
#         elseif x < 4.0 / 6.0
#             Markov.changeTau()
#         elseif x < 5.0 / 6.0
#             Markov.changeExtK()
#         else
#             Markov.changeExtTau()
#         end
#         # println(Curr.T[1], "-->", Curr.T[LastTidx])

#         i % 8 == 0 && Markov.measure()

#         if i % 1000 == 0
#             check(printTimer)
#             check(saveTimer)
#             check(reweightTimer)
#         end
#     end
# end

# Markov.printStatus()
# Markov.save()
# println("End Simulation. ")
