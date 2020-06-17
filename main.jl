include("./utility/constant.jl")
include("parameter.jl")

using Random
using StaticArrays
using .Parameter

PID = length(ARGS) > 0 ? PID = parse(Int, ARGS[1]) : rand(100_000:999_999)
Random.seed!(PID) #initialize rand(), which is already global and by default using MersenneTwister
@assert PID >= 0 "RNG seed must be >=0"


@flush printstyled(
    "Rs: $(para.Rs), kF: $(para.Kf), EF: $(para.Ef), β: $(para.β), T/T_F: $(round(1.0 / para.β / para.Ef, digits = 4))\nSeed: $PID, TotalBlock: $(mcPara.block)\n",
    color = :green,
)

# println(StaticArrays.SVector{3, Int})
# println(Grid.LogGrid{Float, 64, 2})

# ############  Derived parameters ###################################
# @assert Dim == 2 || Dim == 3 "Dim $Dim is not implemented!"
# @assert length(ReWeight) > Order + 1 "More ReWeight elements are needed!"
# const Kf = Dim == 3 ? (9π / 4.0)^(1.0 / 3) / Rs : √2 / Rs
# const π² = π^2
# const Ef, μ, Nf = Kf^2, Kf^2, Spin * Kf / (4π²)
# const β = _Beta / Ef # rescale the temperature
# const PhaseFactor = 1.0 / (2π)^Dim
# const LastTidx = 2 * (Order + 2)
# const LastKidx = Order + 8
# const TauGridSize, FermiKGridSize, BoseKGridSize, AngGridSize = 128, 64, 64, 32
# const TauScale, FermiKScale, BoseKScale, MaxK = 3Ef, 2 / β^0.5, 2 / β^0.5, 3Kf


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
