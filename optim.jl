include("core.jl")
include("measure.jl")
using DelimitedFiles
using ProgressBars
using Printf
using Statistics
using JLD2

function main(L::Int,U::Real,Temp::Real,μ::Real,Ntherm::Int,Nsweep::Int)
    Nt = 80
    Nwrap = 80
    l = lattice(L,U,μ,Temp,Nt,Nwrap)
    result = zeros(Nsweep,2)
    AuxField = rand([1,-1],l.Ns,Nt)
    tmp = cache(l,AuxField)
    MultB = initMultBudt(l,AuxField)
    G = greens(MultB[1])
    unequal_green = init_unequal_green(l)
    accum_unequal_green!(l,unequal_green,G,MultB,tmp)
    unequal_green[2]
end

L = 12
Ntherm = 8000
U = -8.0
μ = -0.32
T = 1.0
Nsweep = 16000
main(L,U,T,μ,Ntherm,Nsweep)