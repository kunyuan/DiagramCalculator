############ Global External Variable Grid #########################
########### Other constants  #####################################
const NORM, SIGMA, POLAR, DELTA, GAMMA_I, GAMMA_T, GAMMA_U, GAMMA_S = 0, 1, 2, 3, 4, 5, 6, 7
using StaticArrays: MVector, SVector
const Float = Float64
const Mom = MVector{3,Float}
const VerWeight = MVector{2,Float}
const IN, OUT = 1, 2
const INL, OUTL, INR, OUTR = 1, 2, 3, 4
const DI, EX = 1, 2
const DOWN, UP = 1, 2
const LEFT, RIGHT = 1, 2
const I, T, U, S, TC, UC = 1, 2, 3, 4, 5, 6

const π = Float(pi)
const π² = π^2

# mutable struct VerWeight <: FieldVector{2,Float}
#     dir::Float
#     ex::Float
# end

@inline squaredNorm(K::MVector{D, Float}) where D = D == 3 ? k[1]^2 + k[2]^2 + k[3]^2 : k[1]^2 + k[2]^2
@inline norm(k::MVector{D, Float}) where D= D == 3 ? sqrt(k[1]^2 + k[2]^2 + k[3]^2) : sqrt(k[1]^2 + k[2]^2)
@inline dot(k::MVector{D, Float}, q::MVector{D, Float}) where D =
    D == 3 ? k[1] * q[1] + k[2] * q[2] + k[3] * q[3] : k[1] * q[1] + k[2] * q[2]

macro flush(content)
    return :($content)
end

