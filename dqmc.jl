include("core.jl")
include("measure.jl")

function main(L::Int,U::Float64,Temp::Float64,Nsweep::Int,Nwarm::Int)
    μ = U / 2
    β = 1/Temp
    Nt = Int(round(β/sqrt(0.06/U)))
    AuxField = rand([1,-1],l.Ns,Nt)
    l = lattice(L,U,μ,Temp,Nt)
    MultBup, MultBdn = initMultBudt(l,AuxField)
    Gup = invoneplus(MultBup[Nt+1])
    Gdn = invoneplus(MultBdn[Nt+1])
    for warm = 1:Nwarm
        sweep!(l,AuxField,Gup,Gdn,MultBup,MultBdn)
    end
    for sweep = 1:Nsweep
        sweep!(l,AuxField,Gup,Gdn,MultBup,MultBdn)
    end
end