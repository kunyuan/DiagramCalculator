module markov
import StaticArrays
using QuantumStatistics
import ..Parameters

mutable struct State{G}
    step::Int
    group::G
    scaleidx::Int
    extTidx::Int
    extKidx::Int
    extAngidx::Int
    absWeight::Float
    T::MVector{LastTidx,Float}
    K::SVector{LastKidx,Mom}
    # groups::Array{Group, 1}

    function State{G}(para, group::G) where {G}
        β, Kf=para.β, para.Kf
        varT = @MVector [rand() .* β for i = 1:LastTidx]
        varK = @SVector [rand(Mom) .* Kf for i = 1:LastKidx]
        # varT = rand(rng, LastTidx) .* Beta
        # varK = rand(rng, Mom, LastKidx) .* Kf
        # rand!(rng, varK)
        curr = new(0, 0, 1, 1, 1, 1, 0.0, varT, varK)
        curr.T[LastTidx] = Grid.tau.grid[curr.extTidx]
        curr.T[1] = 0.0

        if para.DiagType == GAMMA
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

end