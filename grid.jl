module Grid
using StaticArrays: SVector, MVector

struct Coeff{T<:AbstractFloat}
    bound::SVector{2,T}
    idx::SVector{2,T}
    λ::T
    a::T
    b::T

    function Coeff{T}(bound, idx, λ, dense2sparse::Bool) where {T}
        # bound, idx, λ = Tuple{T,T}(_bound), Tuple{T,T}(_idx), T(_λ)

        if dense2sparse == false
            bound = (bound[2], bound[1])
            idx = (idx[2], idx[1])
            λ = -λ
        end
        # println(bound, ", ", idx)
        _l1, _l2 = T(1.0), exp(λ * (idx[2] - idx[1]))
        b = (bound[2] - bound[1]) / (_l2 - _l1)
        a = (bound[1] * _l2 - bound[2] * _l1) / (_l2 - _l1)
        return new{T}(bound, idx, λ, a, b)
    end
end

function _floor(l::Coeff{T}, x::T) where {T}
    # @tmpassert(bound[0] <= x <= bound[1])
    pos = l.idx[1] + T(1.0) / l.λ * log((x - l.a) / l.b)
    return Base.floor(Int, pos)
end

function _grid(l::Coeff{T}, idx::Int) where {T}
    return l.a + l.b * exp(l.λ * (T(idx) - l.idx[1]))
end

function checkOrder(grid)
    # println($(grid))
    for idx = 2:length(grid)
        @assert grid[idx-1] < grid[idx] "The grid at $idx is not in the increase order: \n$grid"
    end
end

struct LogGrid{T,SIZE,SEG} # create a log grid of the type T with SIZE grids and SEG of segments
    size::Int
    grid::MVector{SIZE,T}
    range::SVector{SEG,UnitRange{Int}} # ends of each segments
    coeff::SVector{SEG,Coeff{T}}

    function LogGrid{T,SIZE,SEG}(coeff, range, isopen) where {T<:AbstractFloat,SIZE,SEG}
        @assert SIZE > 1 "Size must be large than 1"
        # coeff = SVector{SEG,Coeff{T}}(_coeff)
        # range = SVector{SEG,UnitRange{Int}}(_range)
        println(range)
        for ri = 2:SEG
            @assert range[ri-1][end] + 1 == range[ri][1] "ranges must be connected to each other"
        end
        @assert range[1][1] == 1 "ranges should start with the idx 1"
        @assert range[end][end] == SIZE "ranges ends with $(range[end][end]), expected $SIZE"

        grid = []
        for s = 1:SEG
            for idx in range[s]
                push!(grid, _grid(coeff[s], idx))
            end
        end
        if isopen[1] == true
            grid[1] += eps(T)
        end
        if isopen[2] == true
            grid[end] -= eps(T)
        end
        checkOrder(grid)
        return new{T,SIZE,SEG}(SIZE, grid, range, coeff)
    end
end

struct UniformGrid{T,SIZE}
    size::Int
    grid::SVector{SIZE,T}

    function UniformGrid{T,SIZE}(head, tail) where {T<:AbstractFloat,SIZE}
        @assert SIZE > 1 "Size must be large than 1"
        grid = LinRange(T(head), T(tail), SIZE)
        return new{T,SIZE}(SIZE, grid)
    end
end

# cos(theta) grids
@inline function angle(size, type = Float64)
    return UniformGrid{type,size}(-1.0, 1.0)
end

@inline function tau(β, halfLife, size::Integer, type = Float64)
    c1 = Grid.Coeff{type}([0.0, 0.5β], [1.0, 0.5size + 0.5], 1.0 / halfLife, true)
    r1 = 1:Int(0.5size)
    c2 = Grid.Coeff{type}([0.5β, β], [0.5size - 0.5, size], 1.0 / halfLife, false)
    r2 = (Int(0.5size)+1):size
    tau = LogGrid{type,size,2}([c1, c2], [r1, r2], [true, true])
    return tau
end

@inline function fermiK(Kf, maxK, halfLife, size, kFi = Int(floor(0.5size)), type = Float64)
    c1 = Grid.Coeff{type}([0.0, Kf], [1.0, kFi + 1.0], 1.0 / halfLife, false)
    r1 = 1:kFi-1
    c2 = Grid.Coeff{type}([Kf, maxK], [kFi, size], 1.0 / halfLife, true)
    r2 = kFi:size
    K = LogGrid{type,size,2}([c1, c2], [r1, r2], [true, false])
    return K
end

@inline function boseK(
    Kf,
    maxK,
    halfLife,
    size,
    kFi = Int(floor(0.5size)),
    twokFi = Int(floor(2size / 3)),
    type = Float64,
)
    λ = 1.0 / halfLife
    c1 = Grid.Coeff{type}([0.0, Kf], [1.0, kFi + 1.0], λ, true)
    r1 = 1:kFi-1
    c2 = Grid.Coeff{type}([Kf, 2.0 * Kf], [kFi, twokFi + 1.0], λ, false)
    r2 = kFi:twokFi-1
    c3 = Grid.Coeff{type}([2.0 * Kf, maxK], [twokFi, size], λ, true)
    r3 = twokFi:size
    K = Grid.LogGrid{type,size,3}([c1, c2, c3], [r1, r2, r3], [true, false])
    return K
end

end
