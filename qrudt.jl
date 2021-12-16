using LinearAlgebra

mutable struct UDT{Type <: Real}
	U::Matrix{Type}
	D::Vector{Type}
	T::Matrix{Type}
end

# iteration for destructuring into components
Base.iterate(S::UDT) = (S.U, Val(:D))
Base.iterate(S::UDT, ::Val{:D}) = (S.D, Val(:T))
Base.iterate(S::UDT, ::Val{:T}) = (S.T, Val(:done))
Base.copy(S::UDT) = UDT(S.U, S.D, S.T)

function udt(A::Matrix{Float64})
    F = qr(A)
    D = diag(F.R)
    UDT(Matrix(F.Q), D, Diagonal(1 ./ D) * F.R)
end

function udtMult(A::UDT{Float64},B::UDT{Float64})
    Mat = A.T * B.U
    lmul!(Diagonal(A.D), Mat)
    rmul!(Mat, Diagonal(B.D))
    F = udt(Mat)
    UDT(A.U * F.U, F.D, F.T * B.T)
end

function invoneplus(F::UDT;u = similar(F.U),t = similar(F.T))
    U, D, T = F
    m = U' / T
    m[diagind(m)] .+= D
    utmp, d, ttmp = udt(m)
    mul!(u, U, utmp)
    mul!(t, ttmp, T)
    return inv(t)*Diagonal(1 ./ d)*u'
end