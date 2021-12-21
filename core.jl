include("qrudt.jl")
using LinearAlgebra

struct lattice
    L::Int
    Ns::Int
    U::Float64
    μ::Float64
    Temp::Float64
    Nt::Int
    Δτ::Float64
    λ::Float64
    Tmatrix::Matrix{Int}
    expmΔτT::UDT{Float64}
    expmhΔτT::UDT{Float64}
    exphΔτT::UDT{Float64}

    function lattice(L::Int,U::Float64,μ::Float64,Temp::Float64,Nt::Int)
        Ns = L^3
        Δτ = 1/(Temp*Nt)
        λ = acosh(exp(abs(U)*Δτ/2))
        Tmatrix = initT(L)
        expmΔτT = udt(exp(-Δτ * Tmatrix))
        exphmΔτT = udt(exp(-Δτ/2 * Tmatrix))
        exphΔτT = udt(exp(Δτ/2 * Tmatrix))
        new(L,Ns,U,μ,Temp,Nt,Δτ,λ,Tmatrix,expmΔτT,exphmΔτT,exphΔτT)
    end
end

function initT(L::Int)
    Tmatrix = zeros(Int,L^3,L^3)
    index::Int = 1
    for z = 1:L, y = 1:L, x = 1:L
        Tmatrix[index, x == 1 ? index+L-1 : index-1] = -1
        Tmatrix[index, x == L ? index-L+1 : index+1] = -1
        Tmatrix[index, y == 1 ? index+L*(L-1) : index-L] = -1
        Tmatrix[index, y == L ? index-L*(L-1) : index+L] = -1
        Tmatrix[index, z == 1 ? index+L^2*(L-1) : index-L^2] = -1
        Tmatrix[index, z == L ? index-L^2*(L-1) : index+L^2] = -1
        index += 1
    end
    Tmatrix
end

function eTeV(eT::UDT,eV::Vector{Float64})
    inveV = 1 ./ eV
    UDT(eT.U,eT.D .* eV, Diagonal(inveV) * eT.T * Diagonal(eV))
end

function greens(R::UDT,L::UDT)
    DRmaxinv = Diagonal(1 ./ max.(R.D,1))
	DLmaxinv = Diagonal(1 ./ max.(L.D,1))
	DRmin = Diagonal(min.(R.D,1))
	DLmin = Diagonal(min.(L.D,1))
	inv(L.T) * DLmaxinv * inv(DRmaxinv * inv(L.T * R.U) * DLmaxinv +
		DRmin * R.T * L.U * DLmin) * DRmaxinv * inv(R.U)
end

function initMultBudt(l::lattice,AuxField::Matrix{Int})
    MultBup = Vector{UDT}(undef,l.Nt+2)
    MultBdn = Vector{UDT}(undef,l.Nt+2)
    UDTI = udt(Matrix(Diagonal(ones(l.Ns))))
    MultBup[1] = copy(UDTI)
    MultBdn[1] = copy(UDTI)
    MultBup[l.Nt+2] = copy(UDTI)
    MultBdn[l.Nt+2] = copy(UDTI)
    for i = 2:l.Nt+1
        MultBup[i] = udtMult(MultBup[i-1],eTeV(l.expmΔτT,exp.( AuxField[:,i-1]*l.λ .+ (l.μ - l.U/2)*l.Δτ)))
        MultBdn[i] = udtMult(MultBdn[i-1],eTeV(l.expmΔτT,exp.(-AuxField[:,i-1]*l.λ .+ (l.μ - l.U/2)*l.Δτ)))
    end
    MultBup, MultBdn
end

function flip!(slice::Int,l::lattice,AuxField::Matrix{Int},Gup::Matrix{Float64},Gdn::Matrix{Float64})
    γup = exp.(-2*l.λ*AuxField[:,slice]).-1
    γdn = exp.( 2*l.λ*AuxField[:,slice]).-1
    Rup = 0
    Rdn = 0
    gtmp = similar(Gup[1,:])
    @inbounds for site = 1:l.Ns
        Rup = 1+(1-Gup[site,site])*γup[site]
        Rdn = 1+(1-Gdn[site,site])*γdn[site]
        P = Rup * Rdn
        if P > 1 || rand() < P
            AuxField[site,slice] *= -1
            updateg!(site,γup[site]/Rup,Gup,gtmp)
            updateg!(site,γdn[site]/Rdn,Gdn,gtmp)
        end
    end
    nothing
end

function updateg!(site::Int,prop::Float64,g::Matrix{Float64},gtmp::Array{Float64})
    gtmp[:] = -g[site,:]
    gtmp[site] += 1
    @views g[:,:] = g[:,:] - prop * g[:,site] * transpose(gtmp)
    nothing
end

function updateBτ0!(slice::Int,l::lattice,AuxField::Matrix{Int},
    MultBup::Vector{UDT},MultBdn::Vector{UDT},gup::Matrix{Float64},gdn::Matrix{Float64})
    MultBup[l.Nt-slice+2] = udtMult(eTeV(l.expmΔτT,
        exp.(AuxField[:,slice]*l.λ .+ (l.μ - l.U/2)*l.Δτ)), MultBup[l.Nt-slice+3])
    gup[:,:] = greens(MultBup[l.Nt-slice+2],MultBup[l.Nt-slice+1])
    MultBdn[l.Nt-slice+2] = udtMult(eTeV(l.expmΔτT,
        exp.(-AuxField[:,slice]*l.λ .+ (l.μ - l.U/2)*l.Δτ)), MultBdn[l.Nt-slice+3])
    gdn[:,:] = greens(MultBdn[l.Nt-slice+2],MultBdn[l.Nt-slice+1])
    nothing
end

function updateBβτ!(slice::Int,l::lattice,AuxField::Matrix{Int},
    MultBup::Vector{UDT},MultBdn::Vector{UDT},gup::Matrix{Float64},gdn::Matrix{Float64})
    MultBup[l.Nt-slice+2] = udtMult(MultBup[l.Nt-slice+1],
        eTeV(l.expmΔτT, exp.(AuxField[:,slice]*l.λ .+ (l.μ - l.U/2)*l.Δτ)))
    Buptmp = eTeV(l.expmΔτT, exp.(AuxField[:,slice-1]*l.λ .+ (l.μ - l.U/2)*l.Δτ))
    gup[:,:] = greens(MultBup[l.Nt-slice+4],udtMult(MultBup[l.Nt-slice+2],Buptmp))
    MultBdn[l.Nt-slice+2] = udtMult(MultBdn[l.Nt-slice+1],
        eTeV(l.expmΔτT, exp.(-AuxField[:,slice]*l.λ .+ (l.μ - l.U/2)*l.Δτ)))
    Bdntmp = eTeV(l.expmΔτT, exp.(-AuxField[:,slice-1]*l.λ .+ (l.μ - l.U/2)*l.Δτ))
    gdn[:,:] = greens(MultBdn[l.Nt-slice+4],udtMult(MultBdn[l.Nt-slice+2],Bdntmp))
    nothing
end

function sweep!(l::lattice,AuxField::Matrix{Int},Gup::Matrix{Float64},Gdn::Matrix{Float64},
    MultBup::Vector{UDT},MultBdn::Vector{UDT})
    for slice = 1:l.Nt-1
        flip!(slice,l,AuxField,Gup,Gdn)
        updateBτ0!(slice,l,AuxField,MultBup,MultBdn,Gup,Gdn)
    end
    for slice = l.Nt:-1:2
        flip!(slice,l,AuxField,Gup,Gdn)
        updateBβτ!(slice,l,AuxField,MultBup,MultBdn,Gup,Gdn)
    end
end