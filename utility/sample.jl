"""
All sampling functions have three versions:
1. new* : function to create a new variable, return params: 1) the new variable (through return value or reference) 2) the proposal probability
2. del* : function to remove the old variable, return the proposal probability
3. shift*  : shift the old variable to a new one, return params: 1) the new variable (through return value or reference) 2) the proposal probability
"""

# Before using the sampling functions, call setRand function to reset the random number generator
# function setRand(_rand) 
#     global const rand=_rand
# end

@inline newIdx(size) = rand(1:size), Float(size)
@inline delIdx(oldIdx, size) = Float(1.0) / size
@inline shiftIdx(oldIdx, size) = rand(1:size), Float(1.0)

@inline newTau(β) = rand(Float) * β, β 
@inline delTau(oldTau, β) = Float(1.0) / β
@inline function shiftTau(oldTau, β, TauScale)
    x = rand(Float)
    newTau::Float = 0.0
    if x < 1.0 / 3
        newTau = oldTau + TauScale * (rand(Float) - Float(0.5))
    elseif x < 2.0 / 3
        newTau = -oldTau
    else
        newTau = rand(Float) * β
    end

    if newTau < Float(0.0)
        newTau += β
    elseif newTau > β
        newTau -= β
    end
    return newTau, Float(1.0)
end

@inline function newK!(newK, Kf, Kscale=Kf/2.0)::Float
    ############ Simple Way ########################
    # for i in 1:DIM
    #     newK[i] = Kf * (rand(rng) - 0.5) * 2.0
    # end
    # return (2.0 * Kf)^DIM
    ################################################

    Kamp = Kf + (2.0*rand() - 1.0) * Kscale
    (Kamp <= 0.0) && (return 0.0)
    # Kf-dK<Kamp<Kf+dK 
    ϕ = 2π * rand(Float)
    if length(newK) == 3
        θ = π * rand(Float)
        newK[1] = Kamp * cos(ϕ) * sin(θ)
        newK[2] = Kamp * sin(ϕ) * sin(θ)
        newK[3] = Kamp * cos(θ)
        return 2Kscale * 2π * π * (sin(θ) * Kamp^2)
        # prop density of Kamp in [Kf-dK, Kf+dK), prop density of ϕ
        # prop density of θ, Jacobian
    elseif length(newK)==2
        newK[1] = Kamp * cos(ϕ)
        newK[2] = Kamp * sin(ϕ)
        return 2Kscale * 2π * Kamp
        # prop density of Kamp in [Kf-dK, Kf+dK), prop density of ϕ, Jacobian
    else
        throw(DimensionMismatch("Dimension $(length(newK)) is not supported!"))
    end
end

@inline function delK(oldK, Kf, Kscale=Kf/2.0)::Float
    ############## Simple Way #########################
    # for i in 1:DIM
    #     if abs(oldK[i]) > Kf
    #         return 0.0
    #     end
    # end
    # return 1.0 / (2.0 * Kf)^DIM
    ####################################################

    Kamp = norm(oldK)
    if !(Kf - Kscale < Kamp < Kf + Kscale)
        return 0.0
    end
    # (Kamp < Kf - dK || Kamp > Kf + dK) && return 0.0
    if length(oldK) == 3
        sinθ = sqrt(oldK[1]^2 + oldK[2]^2) / Kamp
        (sinθ < eps(Float)) && (return 0.0)
        return 1.0 / (2Kscale * 2π * π * sinθ * Kamp^2)
    elseif length(oldK)==2 
        return 1.0 / (2Kscale * 2π * Kamp)
    else
        throw(DimensionMismatch("Dimension $(length(newK)) is not supported!"))
    end
end

@inline function shiftK!(oldK, newK, Kf, Kscale=Kf/3.0, λ=1.5)::Float
    x = rand()
    D=length(newK)
    if x < 1.0 / 3
        # dK = β > 1.0 ? Kf / β * 3.0 : Kf
        # newK .= oldK .+ (rand(rng, DIM) .- 0.5) .* dK
        for i in 1:D
            newK[i] = oldK[i] + (rand() - 0.5) * Kscale
        end
        return 1.0
    elseif x < 2.0 / 3
        ratio = 1.0 / λ + rand() * (λ - 1.0 / λ)
        newK .= oldK .* ratio
        if D==3
            return ratio
        elseif D==2
            return 1.0
        else
            throw(DimensionMismatch("Dimension $(length(newK)) is not supported!"))
        end
    else
        newK .= oldK .* (-1.0)
        return 1.0
    end
end