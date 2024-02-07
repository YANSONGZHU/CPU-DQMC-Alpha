using LoopVectorization

mutable struct UDT{Type <: Real}
	U::Matrix{Type}
	D::Vector{Type}
	T::Matrix{Type}
end

Base.iterate(S::UDT) = (S.U, Val(:D))
Base.iterate(S::UDT, ::Val{:D}) = (S.D, Val(:T))
Base.iterate(S::UDT, ::Val{:T}) = (S.T, Val(:done))
Base.copy(S::UDT) = UDT(S.U, S.D, S.T)
Base.Matrix(F::UDT{Float64}) = (F.U * Diagonal(F.D)) * F.T

function udt!(A::Matrix{Float64})
    U = Matrix{Float64}(undef,size(A))
    D = Vector{Float64}(undef,size(A,1))
    udt_AVX_pivot!(U,D,A)
    UDT(U, D, A)
end

function udt(F::Matrix{Float64})
    A = copy(F)
    udt!(A)
end

function udtMult(A::UDT{Float64},B::UDT{Float64})
    Mat = A.T * B.U
    lmul!(Diagonal(A.D), Mat)
    rmul!(Mat, Diagonal(B.D))
    F = udt(Mat)
    UDT(A.U * F.U, F.D, F.T * B.T)
end

@inline function reflector!(x::Matrix{C}, normu, j=1, n=size(x, 1)) where {C <: Real}
    @inbounds begin
        ξ1 = x[j, j]
        if iszero(normu)
            return zero(ξ1) #zero(ξ1/normu)
        end
        normu = sqrt(normu)
        ν = LinearAlgebra.copysign(normu, real(ξ1))
        ξ1 += ν
        x[j, j] = -ν
        @turbo for i = j+1:n
            x[i, j] /= ξ1
        end
    end
    ξ1/ν
end

function indmaxcolumn(A::Matrix{C}, j=1, n=size(A, 1)) where {C <: Real}
    squared_norm = 0.0
    @turbo for k in j:n
        squared_norm += abs2(A[k, j])
    end
    ii = j
    @inbounds for i in j+1:n
        mi = 0.0
        @turbo for k in j:n
            mi += abs2(A[k, i])
        end
        if abs(mi) > squared_norm
            squared_norm = mi
            ii = i
        end
    end
    return ii, squared_norm
end

# Much higher accuracy, but a bit slower
"""
    udt_AVX_pivot!(
        U::Matrix, D::Vector, T::Matrix[, 
        pivot::Vector, temp::Vector, apply_pivoting = Val(true)
    ])

In-place calculation of a UDT (unitary - diagonal - (upper-)triangular) 
decomposition. The matrix `T` is simultaniously the input matrix that is 
decomposed and the triangular output matrix.

If `apply_pivoting = Val(true)` the `T` matrix will be pivoted such that 
`U * Diagonal(D) * T` matches the input.
If `apply_pivoting = Val(false)` `T` will be a "dirty" upper triangular matrix 
(i.e. with random values elsewhere) which still requires pivoting. 
`rdivp!(A, T, temp, pivot)` is built explicitly for this case - it applies the 
pivoting while calculating `A T^-1`. Warning: BlockDiagonal matrices use 
per-block pivoting.

This assumes correctly sized square matrices as inputs.
"""
function udt_AVX_pivot!(
        U::AbstractArray{C, 2}, 
        D::AbstractArray{C, 1}, 
        input::AbstractArray{C, 2},
        pivot::AbstractArray{Int64, 1} = Vector(UnitRange(1:size(input, 1))),
        temp::AbstractArray{C, 1} = Vector{C}(undef, length(D)),
        apply_pivot::Val = Val(true)
    ) where {C<:Real}
    # Assumptions:
    # - all matrices same size
    # - input can be mutated (input becomes T)

    # @bm "reset pivot" begin
        n = size(input, 1)
        @inbounds for i in 1:n
            pivot[i] = i
        end
    # end

    # @bm "QR decomposition" begin
        @inbounds for j = 1:n
            # Find column with maximum norm in trailing submatrix
            # @bm "get jm" begin
                jm, maxval = indmaxcolumn(input, j, n)
            # end

            # @bm "pivot" begin
                if jm != j
                    # Flip elements in pivoting vector
                    tmpp = pivot[jm]
                    pivot[jm] = pivot[j]
                    pivot[j] = tmpp

                    # Update matrix with
                    @turbo for i = 1:n
                        tmp = input[i,jm]
                        input[i,jm] = input[i,j]
                        input[i,j] = tmp
                    end
                end
            # end

            # Compute reflector of columns j
            # @bm "Reflector" begin
                τj = reflector!(input, maxval, j, n)
                temp[j] = τj
            # end

            # Update trailing submatrix with reflector
            # @bm "apply" begin
                reflectorApply!(input, τj, j, n)
            # end
        end
    # end

    # @bm "Calculate Q" begin
        copyto!(U, I)
        @inbounds begin
            U[n, n] -= temp[n]
            for k = n-1:-1:1
                for j = k:n
                    vBj = U[k,j]
                    @turbo for i = k+1:n
                        vBj += conj(input[i,k]) * U[i,j]
                    end
                    vBj = temp[k]*vBj
                    U[k,j] -= vBj
                    @turbo for i = k+1:n
                        U[i,j] -= input[i,k]*vBj
                    end
                end
            end
        end
        # U done
    # end

    # @bm "Calculate D" begin
        @inbounds for i in 1:n
            # With checkerboard it apparently can happen that the input matrix 
            # takes the form [a[1] * v   a[2] * v   a[3] * v   ...]
            # In this case we should get zeros on the diagonal which would cause 
            # div 0 issues in apply_pivot. To avoid those, we have this ifelse.
            # See #169
            x = abs(input[i, i])
            D[i] = ifelse(x == 0, 1.0, x)
        end
    # end

    # @bm "pivoted zeroed T w/ inv(D)" begin
        _apply_pivot!(input, D, temp, pivot, apply_pivot)
    # end

    nothing
end

function _apply_pivot!(input::Matrix{C}, D, temp, pivot, ::Val{true}) where {C <: Real}
    n = size(input, 1)
    @inbounds for i in 1:n
        d = 1.0 / D[i]
        @inbounds for j in 1:i-1
            temp[pivot[j]] = zero(C)
        end
        @turbo for j in i:n
            temp[pivot[j]] = d * input[i, j]
        end
        @turbo for j in 1:n
            input[i, j] = temp[j]
        end
    end
end

function _apply_pivot!(input::Matrix{C}, D, temp, pivot, ::Val{false}) where {C <: Real}
    n = size(input, 1)
    @inbounds for i in 1:n
        d = 1.0 / D[i]
        @turbo for j in i:n
            input[i, j] = d * input[i, j]
        end
    end
end

@inline function reflectorApply!(x::AbstractVector{<: Real}, τ::Real, A::StridedMatrix{<: Real})
    m, n = size(A)
    @inbounds for j = 1:n
        # dot
        vAj = A[1, j]
        @turbo for i = 2:m
            vAj += conj(x[i]) * A[i, j]
        end

        vAj = conj(τ)*vAj

        # ger
        A[1, j] -= vAj
        @turbo for i = 2:m
            A[i, j] -= x[i]*vAj
        end
    end
    return A
end

@inline function reflectorApply!(M::StridedArray{<: Real}, τ::Real, k::Int, n::Int)
    @inbounds for j = k+1:n
        # dot
        vAj = M[k, j]
        @turbo for i = k+1:n
            vAj += conj(M[i, k]) * M[i, j]
        end

        vAj = conj(τ)*vAj

        # ger
        M[k, j] -= vAj
        @turbo for i = k+1:n
            M[i, j] -= M[i, k]*vAj
        end
    end
    return M
end

function greens(R::UDT,L::UDT)
    DRmaxinv = Diagonal(1 ./ max.(R.D,1))
	DLmaxinv = Diagonal(1 ./ max.(L.D,1))
	DRmin = Diagonal(min.(R.D,1))
	DLmin = Diagonal(min.(L.D,1))
	inv(L.T) * DLmaxinv * inv(DRmaxinv * inv(L.T * R.U) * DLmaxinv +
		DRmin * R.T * L.U * DLmin) * DRmaxinv * inv(R.U)
end

function greens(R::UDT)
    DRmaxinv = Diagonal(1 ./ max.(R.D,1))
	DRmin = Diagonal(min.(R.D,1))
	inv(DRmaxinv * inv(R.U) + DRmin * R.T) * DRmaxinv * inv(R.U)
end