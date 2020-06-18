include("parameter.jl")

using Random
using StaticArrays

PID = length(ARGS) > 0 ? PID = parse(Int, ARGS[1]) : rand(100_000:999_999)
Random.seed!(PID) #initialize rand(), which is already global and by default using MersenneTwister
@assert PID >= 0 "RNG seed must be >=0"

@flush println("Seed/PID: $PID")

include("markov.jl")
markov(Para(), MCPara(), grids)


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
