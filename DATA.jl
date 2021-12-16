include("core.jl")
include("measure.jl")
using JLD2
using DelimitedFiles
using ProgressBars
using Printf

function main(L::Int,U::Float64,Temp::Float64,Nwarm::Int,Nsweep::Int)
    println("DQMC START!")
    @printf("L = %d U  = %g Temp = %g\n",L,U,Temp)
    resultfile = "DATA/L " * string(L) * " U " * string(U) * " Temp " * string(Temp) * ".dat"
    AuxFieldfile = "AuxField/L " * string(L) * " U " * string(U) * " Temp " * string(Temp) * ".jld2"
    result = zeros(1,5)
    io = open(resultfile, "a")
    μ = U / 2
    β = 1/Temp
    Nt = Int(round(β/sqrt(0.06/U)))
    @printf("Choosing Nt = %d\n",Nt)
    l = lattice(L,U,μ,Temp,Nt)
    println("Initializing ...")
    AuxField = rand([1,-1],l.Ns,Nt)
    MultBup, MultBdn = initMultBudt(l,AuxField)
    Gup = invoneplus(MultBup[Nt+1])
    Gdn = invoneplus(MultBdn[Nt+1])
    println("Warming ...")
    for warm in ProgressBar(1:Nwarm)
        sweep!(l,AuxField,Gup,Gdn,MultBup,MultBdn)
    end
    println("Measure ...")
    for sweep = ProgressBar(1:Nsweep)
        sweep!(l,AuxField,Gup,Gdn,MultBup,MultBdn)
        saveallobser!(Gup,Gdn,l,result)
        writedlm(io, result, '\t')
    end
    close(io)
    jldsave(AuxFieldfile;AuxField)
end

L = 4
Nwarm = 2000
Nsweep = 2000
for U = 6.0:1.0:12.0
    for Temp = 0.25:0.05:0.8
        main(L,U,Temp,Nwarm,Nsweep)
    end
end
end