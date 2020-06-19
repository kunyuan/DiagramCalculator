# include("./oneBody/group.jl")
# using .OneBody

function markov(state, groups, steps, updates, weight = ones(length(updates)), rng = Random.GLOBAL_RNG)
    z = sum(weight)
    N = length(updates)
    prop = [sum(weight[1:iw]) / z for (iw, w) in enumerate(weight)]
    println(prop)
    for s in 1:steps
        x = rand(rng)
        for (i, p) in prop
            if x < p
                # if updates(i)
                # end
            end
        end
        return state, groups
    end