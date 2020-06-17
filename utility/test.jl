import Test: @test, @testset
include("constant.jl")

@testset "sample functions" begin
    include("sample.jl")
    Size=64
    idx, prop=Sample.newIdx(Size)
    @test prop≈1.0/Sample.delIdx(idx, Size)

    β=10.0
    tau, prop=Sample.newTau(β)
    @test prop≈1.0/Sample.delTau(tau, β)

    Kf=1.0
    K=zero(Mom)
    prop=Sample.newK!(K, Kf)
    OldK=copy(K)
    @test prop≈1.0/Sample.delK(K, Kf)
    @test K≈OldK #make sure oldK will not be changed

    NewK=zero(Mom)
    #mulitple attempts to test all internal sampling possibilities
    for i in 1:100
        Sample.shiftK!(OldK, NewK, Kf)
    end
    @test K≈OldK #make sure shiftK! will not change OldK
end