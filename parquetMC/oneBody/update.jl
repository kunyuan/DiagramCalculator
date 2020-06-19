include("../../utility/sample.jl")

function changeGroup(para, state, groups, rng = Random.GLOBAL_RNG)
    β, Kf = para.β, para.Kf
    # random select a new group
    new, old = rand(rng, groups), state.group  
    (new == old) && return

    δKnum = new.innerK[2] - old.innerK[2]
    δTnum = new.innerT[2] - old.innerT[2]

    # if there are too many K or T variables to sample,
    # the acceptance ratio will be too low
    (δKnum > 1 || δTnum > 1) && return 

    propK = newK!(state.K[new.innerK[2]], Kf, Kf / 2.0, rng)
    propT, state.T[new.innerT[2]] = newTau(β, rng)

    newAbsWeight = abs(evaluate(para, state))
    R = propT * paraK * newAbsWeight * new.reweight / state.absWeight / old.reweight

    if rand(rng) < R
        state.group = new
        state.absWeight = newAbsWeight
    end
end


function propK(newGroup)
    return
end

function rejectT(newGroup)
    return
end

function propExtK(newGroup)
    return
end

function rejectExtK(newGroup)
    return
end

function propExtT(newGroup)
    return
end

function rejectExtT(newGroup)
    return
end