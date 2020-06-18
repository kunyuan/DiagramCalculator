include("./parquetMC/oneBody/group.jl")
using .OneBody

function markov(para, mcpara, grids)
    state, groups=OneBody.init(para, grids)
end