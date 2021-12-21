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
    result = zeros(Nsweep,5)
    μ = U / 2
    β = 1/Temp
    Nt = Int(round(β/0.05))
    @printf("Choosing Nt = %d\n",Nt)
    l = lattice(L,U,μ,Temp,Nt)
    println("Initializing ...")
    AuxField = rand([1,-1],l.Ns,Nt)
    MultBup, MultBdn = initMultBudt(l,AuxField)
    Gup = greens(MultBup[Nt+1],MultBup[Nt+2])
    Gdn = greens(MultBdn[Nt+1],MultBdn[Nt+2])
    println("Warming ...")
    for warm in ProgressBar(1:Nwarm)
        sweep!(l,AuxField,Gup,Gdn,MultBup,MultBdn)
    end
    println("Measure ...")
    for sweep = ProgressBar(1:Nsweep)
        sweep!(l,AuxField,Gup,Gdn,MultBup,MultBdn)
        saveallobser!(Gup,Gdn,l,result,sweep)
    end
    io = open(resultfile, "a")
    writedlm(io, result, '\t')
    close(io)
    jldsave(AuxFieldfile;AuxField)
end

L = 4
Nwarm = 2000
Utmp = [6.0 7.0 8.0 9.0 10.0 11.0 12.0]
Ttmp = [0.32 0.36 0.40 0.44 0.48 0.54 0.60 0.72 0.84]
for U = 1:7
    for Temp = 1:9
        Nsweep = 24000-2000*(7-U)-1000*(Temp-1)
        main(L,Utmp[U],Ttmp[Temp],Int(Nwarm),Int(Nsweep))
    end
end