using Parameters

Dim=3
include("utility/constant.jl")
import QuantumStatistics:Grid
# const QS=QuantumStatistics

_rs=1.0
_kf = (Dim == 3) ? ((9π / 4.0)^(1.0 / 3) / _rs) : (sqrt(2) / _rs)
_beta, _mass2, _lambda2=40.0, 0.0, 1.0


@with_kw struct Para
    dim::Int=Dim
    spin::Int=2
    order::Int=1
    diagrams=[Diag.Π, ]
    β::Float =_beta/_kf^2
    Rs::Float=_rs
    m²::Float=_mass2
    λ²::Float=_lambda2
    Kf::Float = _kf
    Ef::Float = _kf^2
    μ::Float=_kf^2
    boldG::Bool=false
end

@with_kw struct MCPara
    block::Int=101
    reweight::Array{Float, 1}=[1.0, 0.1, 30.0, 1.0, 0.2, 0.1, 0.01, 0.01]
    printTime::Int=5
    saveTime::Int=10
    reweightTime::Int=10
    messageTime::Int=10
end

const para=Para()
const mcPara=MCPara()
const grids=(
    tau=Grid.tau(para.β, 3*para.Kf^2, 128),
    fermiK=Grid.fermiK(para.Kf, 3.0*para.Kf, 0.3*para.Kf, 64),
    boseK=Grid.boseK(para.Kf, 3.0*para.Kf, 0.3*para.Kf, 64)
)

# const groups=[(∅, 0), ]
# for order in 1:para.order
#     for diag in para.diagrams
#         push!(groups, (diag, order))
#     end
# end

@flush printstyled(
    "Rs: $(para.Rs), kF: $(para.Kf), EF: $(para.Ef), β: $(para.β), T/T_F: $(round(1.0 / para.β / para.Ef, digits = 4)), TotalBlock: $(mcPara.block)\n",
    color = :green,
)

