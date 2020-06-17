module OneBodyGroup

include("../../utility/constant.jl")
using QuantumStatistics

struct OneBody
    diag::Int
    order::Int
    extTidx::Int
    innerK::Tuple{Int, Int}
    innerT::Tuple{Int, Int}
    data::Array{Float, 2} #kgrid, taugrid
    # TauGrid::LogGrid{Float, 
    # KGrid::Array{Float, 1}

    # function OneBodyGroup(_diag, _order, _extTidx)
    #     if _diag==NORM
    #         InnerT={1, 1}
    #         InnerK={1, 1}
    #     elseif _diag==POLAR
    #         InnerT={2, }
    #         InnerK={2, 2+_order-1}
    #     end
    #     new(_diag, _order, _extTidx)
    # end
end

struct Norm<:OneBody
    Norm(extTidx)=new(NORM, 0, extTidx, (1, 1), (1, 1), [1.0e-10])
end

function measure(group::Norm, factor)
    group.data[0]+=1.0*factor;
end

# struct Polar<:OneBody
#     function Polar(order, extTidx)
#     end
# end

include("one_body_update.jl")
end