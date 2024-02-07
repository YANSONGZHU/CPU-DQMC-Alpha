include("core.jl")
include("measure.jl")
using DelimitedFiles
using ProgressBars
using Printf
using Statistics
using JLD2

function main(L::Int,U::Real,Temp::Real,μ::Real,Ntherm::Int,Nsweep::Int)
    Nt = 81
    Nwrap = 9
    l = lattice(L,U,μ,Temp,Nt,Nwrap,dim=3)
    result = zeros(Nsweep,2)
    MultB = Vector{UDT}(undef, l.NumB)
    AuxField = rand([1,-1],l.Ns,Nt)
    tmp = cache(l,AuxField)
    calcuMultBudt!(MultB,l,tmp)
    G = greens(MultB[l.NumB])
    unequal_green = zeros(Float64,l.Nt,l.Ns,l.Ns)
    for therm = ProgressBar(1:Ntherm)
        sweep!(l,AuxField,G,MultB,tmp)
    end
    for sweep = ProgressBar(1:Nsweep)
        sweep_forward!(l,AuxField,G,MultB,tmp)
        result[sweep,1] = occupy(G,l.Ns)
        result[sweep,2] = doubleoccupy(G,l.Ns)
        accum_unequal_green!(l,unequal_green,G,MultB,tmp)
        calcuMultBudt!(MultB,l,tmp)
    end
    println(sum(result[:,1])/Nsweep)
    println(sum(result[:,2])/Nsweep)
    println(2*sum(result[:,2])/sum(result[:,1]))
    unequal_green = unequal_green ./ Nsweep
    unequal_green_file = "L" * string(L) * " U" * string(U) * " T" * string(Temp) * ".jld2"
    jldsave(unequal_green_file;unequal_green)
end

L = 6
# Ntherm = [6000 8000 10000 12000 14000 16000]
# U = -8.0
# T = [1.0    0.8    1.2    1.6    0.6    2.0]
# μ = [-0.371 -0.331 -0.413 -0.502 -0.292 -0.637]
# U = -4.0
# T = [2.0 1.6 1.2 1.0 0.8 0.6]
# μ = [-0.76 -0.65 -0.55 -0.5 -0.45 -0.41]
# Nsweep = [8000 10000 12000 14000 16000 18000]
U = [-4.0 -6.0 -6.0 -6.0 -6.0 -6.0 -6.0]
T = [0.6 2.0 1.6 1.2 1.0 0.8 0.6]
μ = [-0.41 -0.69 -0.56 -0.47 -0.42 -0.38 -0.34]
Ntherm = [16000 6000 8000 10000 12000 14000 16000]
Nsweep = [18000 8000 10000 12000 14000 16000 18000]
for i = 1:7
    main(L,U[i],T[i],μ[i],Ntherm[i],Nsweep[i])
end