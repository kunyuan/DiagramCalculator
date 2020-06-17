module Parameter
include("utility/constant.jl")

D=3
rs=1.0
kf = (D == 3) ? ((9π / 4.0)^(1.0 / 3) / rs) : (sqrt(2) / rs)
beta, mass2, lambda2=40.0, 0.0, 1.0

const para=(
    dim=3,
    spin=2,
    order=1,
    groups=[(NORM, 0), (POLAR, 1), ],
    β =beta/kf^2,
    Rs=rs,
    m²=mass2,
    λ²=lambda2,
    Kf = kf,
    Ef = kf^2,
    μ=kf^2,
    boldG=false
)

const mcPara=(
    block=101,
    reweight=[1.0, 0.1, 30.0, 1.0, 0.2, 0.1, 0.01, 0.01],
    printTime=5,
    saveTime=10,
    reweightTime=10,
    messageTime=10
)

const gridPara=(
    tau=(128, para.β, 3*kf^2),
    fermiK=(64, 3.0, 2 / para.β^0.5),
    boseK=(64, 3.0, 2 / para.β^0.5),
    angleSize=32
)

export para, mcPara, gridPara
end