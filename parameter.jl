const GAMMA, SIGMA, POLAR, DELTA = 1, 2, 3, 4
################ Global parameters  ##################################
const DiagType = POLAR
# const DiagType = GAMMA
const Order = 1
const TotalBlock = 101
const beta, Rs, Mass2, Lambda, maxK = 40.0, 1.0, 0.0, 1.0, 3.0
const ReWeight = [1.0, 0.1, 30.0, 1.0, 0.2, 0.1, 0.01, 0.01]
const DIM, SPIN = 3, 2
const BoldG, BoldVer = true, true
const TauGridSize, KGridSize, AngGridSize = 128, 64, 64
const PrintTime, SaveTime, ReWeightTime, MessageTime, CollectTime = 5, 10, 30, 10, 10

############  Derived parameters ###################################
@assert DIM == 2 || DIM == 3 "DIM not implemented!"
@assert length(ReWeight) > Order + 1 "More ReWeight elements are needed!"
const Kf = (DIM == 3) ? (9.0 * pi / 4.0)^(1.0 / 3) / Rs : sqrt(2.0) / Rs
const Ef, Mu, Nf = (Kf^2, Kf^2, Kf / (4.0 * pi^2) * SPIN)
const Beta = beta / Ef # rescale the temperature
const MaxK = maxK * Kf
const PhaseFactor = 1.0 / (2π)^DIM
const LastTidx = 2 * (Order + 2)
const LastKidx = Order + 8

############ Global External Variable Grid #########################
########### Other constants  #####################################
using StaticArrays: MVector
const Float = Float64
const Mom = MVector{3,Float}
const VerWeight = MVector{2,Float}
const IN, OUT = 1, 2
const INL, OUTL, INR, OUTR = 1, 2, 3, 4
const DI, EX = 1, 2
const DOWN, UP = 1, 2
const LEFT, RIGHT = 1, 2
const I, T, U, S, TC, UC = 1, 2, 3, 4, 5, 6

# mutable struct VerWeight <: FieldVector{2,Float}
#     dir::Float
#     ex::Float
# end

########## Global function  #######################################

@inline function lastInnerTidx(order)
    return order
    # if DiagType == SIGMA || DiagType == DELTA
    #     return order - 1
    # elseif DiagType == POLAR
    #     return order
    # else
    #     return order + 1 # GAMMA
    # end
end

@inline function firstInnerKidx()
    if DiagType == GAMMA
        return 5
    elseif DiagType == SIGMA || DiagType == POLAR || DiagType == DELTA
        return 2
    end
end

@inline lastInnerKidx(order) = firstInnerKidx() + order - 1

# import Yeppp
# @inline dot(k, q) = Yeppp.dot(k, q)
# @inline norm(k) = sqrt(dot(k, k))
# @inline squaredNorm(k) = dot(k, k)
@inline squaredNorm(k) = DIM == 3 ? k[1]^2 + k[2]^2 + k[3]^2 : k[1]^2 + k[2]^2
@inline norm(k) = DIM == 3 ? sqrt(k[1]^2 + k[2]^2 + k[3]^2) : sqrt(k[1]^2 + k[2]^2)
@inline dot(k, q) =
    DIM == 3 ? k[1] * q[1] + k[2] * q[2] + k[3] * q[3] : k[1] * q[1] + k[2] * q[2]

# module Weight
# mutable struct VerWeight
#     dir::Float64
#     ex::Float64
#     function VerWeight()
#         new(0.0, 0.0)
#     end
#     function VerWeight(x, y)
#         new(x, y)
#     end
# end

# Base.:+(x::VerWeight, y::VerWeight) = VerWeight(x.dir + y.dir, x.ex + y.ex)
# Base.:-(x::VerWeight, y::VerWeight) = VerWeight(x.dir - y.dir, x.ex - y.ex)
# Base.:*(x::VerWeight, y::Number) = VerWeight(x.dir * y, x.ex * y)
# Base.:/(x::VerWeight, y::Number) = VerWeight(x.dir / y, x.ex / y)
# export VerWeight
# end
