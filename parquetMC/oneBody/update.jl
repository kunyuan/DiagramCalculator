include("../../utility/sample.jl")

function increaseOrder(para, state, groups, rng = Random.GLOBAL_RNG)
    # random select a new group
    new, old = rand(rng, groups), state.group  
    (new.order != old.order + 1) && return

    β, Kf = para.β, para.Kf
    δKnum = new.innerK[2] - old.innerK[2]
    δTnum = new.innerT[2] - old.innerT[2]
    propK, propT = 1.0, 1.0

    for ik in old.innerK[2]:new.innerK[2]
        propK = newK!(state.K[ik], Kf, Kf / 2.0, rng)
    end
    for it in old.innerT[2]:new.innerT[2]
        prop, state.T[new.innerT[2]] = newTau(β, rng)
        propT *= prop
    end

    state.propose[:increaseOrder] += 1
    newAbsWeight = abs(evaluate(para, state))
    R = propT * paraK * newAbsWeight * new.reweight / state.absWeight / old.reweight
    if rand(rng) < R
        state.accept[:increaseOrder] += 1
        state.group = new
        state.absWeight = newAbsWeight
    end
end

function decreaseOrder(para, state, )