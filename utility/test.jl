import Test: @test, @testset
include("constant.jl")

@testset "sample functions" begin
    include("sample.jl")
    Size=64
    idx, prop=newIdx(Size)
    @test prop≈1.0/delIdx(idx, Size)

    β=10.0
    tau, prop=newTau(β)
    @test prop≈1.0/delTau(tau, β)

    Kf=1.0
    K=zero(Mom)
    prop=newK!(K, Kf)
    OldK=copy(K)
    @test prop≈1.0/delK(K, Kf)
    @test K≈OldK #make sure oldK will not be changed

    NewK=zero(Mom)
    #mulitple attempts to test all internal sampling possibilities
    for i in 1:100
        shiftK!(OldK, NewK, Kf)
    end
    @test K≈OldK #make sure shiftK! will not change OldK
end