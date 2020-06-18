############ Global External Variable Grid #########################
########### Other constants  #####################################
const Diag=(
    ∅ =0,  #\emptyset, normalization diagram
    Σ =1,  #self energy
    Π =2,  #polarization
    Δ =3,  #annormal self energy
    Γ₄=4, #4-point one-particle-irreducible vertex function
    Γi=5, #irreduble Γ₄
    Γt=6, #particle-hole Γ₄
    Γu=7, #particle-hole-exchange Γ₄
    Γs=8, #particle-particle Γ₄
    Γ₆=9 #6-point one-particle-irreducible vertex function
)

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

