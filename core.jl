include("qrudt.jl")
using LinearAlgebra
using LoopVectorization

struct lattice
    L::Int
    Ns::Int
    U::Float64
    Temp::Float64
    Nt::Int
    Nwrap::Int
    NumB::Int
    Δτ::Float64
    λ::Float64
    Tmatrix::Matrix{Float64}
    eΔτT::Matrix{Float64}

    function lattice(L::Int,U::Float64,μ::Float64,Temp::Float64,Nt::Int,Nwrap::Int;dim=3)
        Ns = L^dim
        Δτ = 1/(Temp*Nt)
        NumB = Int(Nt/Nwrap)
        λ = acosh(exp(abs(U)*Δτ/2))
        Tmatrix = initT(L,μ,dim)
        eΔτT = exp(-Δτ * Tmatrix)
        new(L,Ns,U,Temp,Nt,Nwrap,NumB,Δτ,λ,Tmatrix,eΔτT)
    end
end

struct cache
    gtmp::Vector{Float64}
    IG::Vector{Float64}
    G::Vector{Float64}
    Btmp::Matrix{Float64}
    B::Vector{Matrix{Float64}}
    Im::Matrix{Float64}
    MultBtmp::Matrix{Float64}
    MultBtmptmp::Matrix{Float64}
    UDTI::UDT
    function cache(l::lattice,AuxField::Matrix{Int})
        gtmp = Vector{Float64}(undef, l.Ns)
        IG = Vector{Float64}(undef, l.Ns)
        G = Vector{Float64}(undef, l.Ns)
        Btmp = Matrix{Float64}(undef, l.Ns, l.Ns)
        B = Vector{Matrix{Float64}}(undef, l.Nt)
        for i = 1:l.Nt
            @views B[i] = l.eΔτT * Diagonal(exp.(AuxField[:,i]*l.λ))
        end
        Im = Matrix(Diagonal(ones(l.Ns)))
        MultBtmp = Matrix{Float64}(undef, l.Ns, l.Ns)
        MultBtmptmp = Matrix{Float64}(undef, l.Ns, l.Ns)
        UDTI = udt!(Matrix(Diagonal(ones(l.Ns))))
        new(gtmp,IG,G, Btmp,B,Im,MultBtmp,MultBtmptmp,UDTI)
    end
end

function initT(L::Int,μ::Real,dim::Int)
    Ns = L^dim
    Tmatrix = zeros(Float64,Ns,Ns)
    for i = 1:Ns
        Tmatrix[i,i] = -μ
    end
    index::Int = 1
    if dim == 3
        for z = 1:L, y = 1:L, x = 1:L
            Tmatrix[index, x == 1 ? index+L-1 : index-1] = -1
            Tmatrix[index, x == L ? index-L+1 : index+1] = -1
            Tmatrix[index, y == 1 ? index+L*(L-1) : index-L] = -1
            Tmatrix[index, y == L ? index-L*(L-1) : index+L] = -1
            Tmatrix[index, z == 1 ? index+L^2*(L-1) : index-L^2] = -1
            Tmatrix[index, z == L ? index-L^2*(L-1) : index+L^2] = -1
            index += 1
        end
    elseif dim == 2
        for y = 1:L, x = 1:L
            Tmatrix[index, x == 1 ? index+L-1 : index-1] = -1
            Tmatrix[index, x == L ? index-L+1 : index+1] = -1
            Tmatrix[index, y == 1 ? index+L*(L-1) : index-L] = -1
            Tmatrix[index, y == L ? index-L*(L-1) : index+L] = -1
            index += 1
        end
    end
    Tmatrix
end

function eTeV(eT::Matrix{Float64},eV::Vector{Float64})
    udt!(eT * Diagonal(eV))
end

function calcuMultBudt!(MultB::Vector{UDT},l::lattice,tmp::cache)
    MultBtmp = I
    slice = l.Nt
    for i = 1:l.Nwrap
        MultBtmp = MultBtmp * tmp.B[slice]
        slice -= 1
    end
    MultB[1] = udt!(MultBtmp)
    for i = 2:l.NumB
        MultBtmp = I
        for j = 1:l.Nwrap
            MultBtmp = MultBtmp * tmp.B[slice]
            slice -= 1
        end
        MultB[i] = udtMult(MultB[i-1],udt!(MultBtmp))
    end
end

function flip!(slice::Int,l::lattice,AuxField::Matrix{Int},G::Matrix{Float64},tmp::cache)
    γ = exp.(-2*l.λ*AuxField[:,slice]).-1
    R = 0
    gtmp = similar(G[1,:])
    @inbounds for site = 1:l.Ns
        R = 1+(1-G[site,site])*γ[site]
        P = R * R * exp.(2*l.λ*AuxField[site,slice])
        if P > 1 || rand() < P
            AuxField[site,slice] *= -1
            update_greens!(tmp, γ[site]/R, G, site)
        end
    end
    nothing
end

function update_greens!(tmp::cache, prop, G, i)
    # calculate (I - G)[:, i:N:end]
    vsub!(tmp.IG, I, G, i)

    # calculate {Δ R⁻¹} * G[i:N:end, :]
    vmul!(tmp.G, prop, G, i)

    # update greens function 
    # G[m, n] -= {(I - G)[m, i:N:end]} {{Δ R⁻¹} * G[i:N:end, n]}
    vsubkron!(G, tmp.IG, tmp.G)

    nothing
end

function vsub!(trg::Vector{Float64}, ::UniformScaling, src::Matrix{Float64}, i::Int)
    @turbo for j in eachindex(trg)
        trg[j] = - src[j, i]
    end
    @inbounds trg[i] += 1.0
    nothing
end

function vmul!(trg::Vector{Float64}, M::Float64, src::Matrix{Float64}, i::Int)
    @turbo for j in eachindex(trg)
        trg[j] = M * src[i, j]
    end
    nothing
end

function vsubkron!(G::Matrix{Float64}, L::Vector{Float64}, R::Vector{Float64})
    @turbo for k in eachindex(L), l in eachindex(R)
        G[k, l] -= L[k] * R[l]
    end  
end

function sweep_forward!(l::lattice,AuxField::Matrix{Int},g::Matrix{Float64},MultB::Vector{UDT},tmp::cache)
    slice = 0
    for i = l.NumB:-1:1
        MultBtmp = I
        for j = 1:l.Nwrap
            slice += 1
            flip!(slice,l,AuxField,g,tmp)
            @views tmp.B[slice][:,:] = l.eΔτT * Diagonal(exp.(AuxField[:,slice]*l.λ))
            MultBtmp = tmp.B[slice] * MultBtmp
            if j == l.Nwrap
                if i == 1
                    MultB[i] = udtMult(udt!(MultBtmp),MultB[i+1])
                    g[:,:] = greens(MultB[i])
                elseif i == l.NumB
                    MultB[i] = udt!(MultBtmp)
                    g[:,:] = greens(MultB[i],MultB[i-1])
                else
                    MultB[i] = udtMult(udt!(MultBtmp),MultB[i+1])
                    g[:,:] = greens(MultB[i],MultB[i-1])
                end
            else
                g[:,:] = tmp.B[slice]*g*inv(tmp.B[slice])
            end
        end
    end
    nothing
end

function accum_unequal_green!(l::lattice,unequal_green::Array{Float64},G::Matrix{Float64},MultB::Vector{UDT},tmp::cache)
    currentg = G*inv(tmp.B[1])
    initgudt = udt(currentg)
    for i = 1:l.NumB
        for j = 1:l.Nwrap
            slice = (i-1)*l.Nwrap + j
            if j == l.Nwrap
                currentg = Matrix(udtMult(MultB[l.NumB-i+1],initgudt))
            else
                currentg = tmp.B[slice]*currentg
            end
            unequal_green[slice,:,:] += currentg
        end
    end
end

function sweep!(l::lattice,AuxField::Matrix{Int},G::Matrix{Float64},MultB::Vector{UDT},tmp::cache)
    sweep_forward!(l,AuxField,G,MultB,tmp)
    calcuMultBudt!(MultB,l,tmp)
end