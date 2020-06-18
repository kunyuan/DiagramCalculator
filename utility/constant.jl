############ Global External Variable Grid #########################
########### Other constants  #####################################
const DiagType=[
    :∅  #\emptyset, normalization diagram
    :Σ  #self energy
    :Π  #polarization
    :Δ  #annormal self energy
    :Γ₄ #4-point one-particle-irreducible vertex function
    :Γi #irreduble Γ₄
    :Γt #particle-hole Γ₄
    :Γu #particle-hole-exchange Γ₄
    :Γs #particle-particle Γ₄
    :Γ₆ #6-point one-particle-irreducible vertex function
]

using StaticArrays: MVector
const Float = Float64
# const Mom = MVector{Dim,Float}
const VerWeight = MVector{2,Float}
const IN, OUT = 1, 2
const INL, OUTL, INR, OUTR = 1, 2, 3, 4
const DI, EX = 1, 2
const DOWN, UP = 1, 2
const LEFT, RIGHT = 1, 2
const I, T, U, S, TC, UC = 1, 2, 3, 4, 5, 6

const π² = π^2


# mutable struct VerWeight <: FieldVector{2,Float}
#     dir::Float
#     ex::Float
# end

macro flush(content)
    return :($content)
end

