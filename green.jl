module Green

@inline function bareFermi(β::T, τ::T, Ek::T, scale::T = T(0.0)) where {T<:AbstractFloat}
    if τ == T(0.0)
        τ = -prevfloat(T(0.0))
    end
    G = sign(τ)
    if τ < T(0.0)
        τ += β
    end
    @assert T(0.0) < τ <= β "τ must be (0.0, β]"
    x = β * Ek / 2
    y = 2τ / β - 1
    if -T(100.0) < x < T(100.0)
        G *= exp(-x * y) / 2 * cosh(x)
    elseif x >= T(100.0)
        G *= exp(-x * (y + 1))
    else # x<=-100.0
        G *= exp(x * (1 - y))
    end
    return G
end

# function interaction(mom::Mom, verOrder::Int = 0)
#     q2 = squaredNorm(mom)
#     weight = q2 > 1.0e-14 ? -8π / (q2 + Mass2 + Lambda) : 0.0
#     if verOrder > 0
#         weight *= (weight * Lambda / 8π)^verOrder
#     end
#     return weight
# end

# function interaction(
#     kInL::Mom,
#     kOutL::Mom,
#     kInR::Mom,
#     kOutR::Mom,
#     Boxed::Bool,
#     extQ = -1.0,
# )
#     weight = VerWeight(0.0, 0.0)
#     qDi2 = squaredNorm(kInL - kOutL)
#     weight[DI] = -8π / (qDi2 + Mass2 + Lambda)
#     if (DiagType == SIGMA && qDi2 < 1.0e-14) ||
#        (DiagType == POLAR && abs(qDi2 - extQ^2) < 1.0e-14)
#         weight[DI] = 0.0
#     end

#     if Boxed == false
#         qEx2 = squaredNorm(kInL - kOutR)
#         weight[EX] = 8π / (qEx2 + Mass2 + Lambda)
#         if (DiagType == SIGMA && qEx2 < 1.0e-14) ||
#            (DiagType == POLAR && abs(qEx2 - extQ^2) < 1.0e-14)
#             weight.ex = 0.0
#         end
#     else
#         weight[EX] = 0.0
#     end
#     return weight
# end

# counterBubble(K::Mom) =
#     Lambda / (8π * Nf) * green(Beta / 2.0, K) * green(-Beta / 2.0, K)

# angle(K1::Mom, K2::Mom) = dot(K1, K2) / norm(K1) / norm(K2)

end
