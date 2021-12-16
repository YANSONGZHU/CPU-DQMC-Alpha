function occupy(gup::Matrix{Float64},gdn::Matrix{Float64},Ns::Int)
    occupy = 2 - (sum(diag(gup)) + sum(diag(gdn)))/Ns
end

function doubleoccupy(gup::Matrix{Float64},gdn::Matrix{Float64},Ns::Int)
    doubleoccupy = sum((1 .- diag(gup)).*(1 .- diag(gdn)))/Ns
end

function kinetic(gtildeup::Matrix{Float64},gtildedn::Matrix{Float64},Ns::Int,Tmatrix::Matrix{Int})
    sumk = sum(gtildeup .* Tmatrix) + sum(gtildedn .* Tmatrix)
    sumk / (2 * Ns)
end

function index2xyz(index::Int, L::Int)
	n = index - 1
	xyz = zeros(Int,3)
	for i = 1:3
		xyz[i] = n % L
		n = n ÷ L
	end
	xyz .+ 1
end

function Sπ(gup::Matrix{Float64},gdn::Matrix{Float64},
        gtildeup::Matrix{Float64},gtildedn::Matrix{Float64},l::lattice)
    Sπ = 0
    Q = [pi, pi, pi]
    for i = 1:l.Ns
        for j = 1:l.Ns
            c = real(exp(im*sum(Q.*(index2xyz(i,l.L)-index2xyz(j,l.L)))))
            Sπ += c * (gtildeup[i,i] * gtildeup[j,j] + gtildeup[i,j] * gup[i,j] +
                gtildedn[i,i] * gtildedn[j,j] + gtildedn[i,j] * gdn[i,j] -
                gtildedn[i,i] * gtildeup[j,j] - gtildeup[i,i] * gtildedn[j,j])
        end
    end
    Sπ / l.Ns
end

function nnspincorr(gup::Matrix{Float64},gdn::Matrix{Float64},
        gtildeup::Matrix{Float64},gtildedn::Matrix{Float64},Tmatrix::Matrix{Int})
    Ns = size(gup,1)
    P = 0
    for i = 1:Ns
        for j = 1:Ns
            if Tmatrix[i,j] == -1
                P += (gtildeup[i,i] * gtildedn[j,j]) - 
                    (gtildeup[i,i] * gtildedn[j,j] * gtildedn[i,i]) -
                    (gtildeup[i,i] * gtildedn[j,i] * gdn[j,i]) -
                    (gtildeup[i,i] * gtildedn[j,j] * gtildeup[j,j]) -
                    (gtildeup[i,j] * gup[i,j] * gtildedn[j,j]) +
                    (gtildeup[i,i] * gtildedn[j,j] * gtildedn[i,i] * gtildeup[j,j]) + 
                    (gtildeup[i,i] * gup[i,j] * gtildedn[i,i] * gtildeup[j,j]) +
                    (gtildeup[i,i] * gtildedn[j,i] * gdn[j,i] * gtildeup[j,j]) + 
                    (gtildeup[i,j] * gup[i,j] * gtildedn[j,i] * gdn[j,i])
                P += (gtildedn[i,i] * gtildeup[j,j]) - 
                (gtildedn[i,i] * gtildeup[j,j] * gtildeup[i,i]) -
                (gtildedn[i,i] * gtildeup[j,i] * gup[j,i]) -
                (gtildedn[i,i] * gtildeup[j,j] * gtildedn[j,j]) -
                (gtildedn[i,j] * gdn[i,j] * gtildeup[j,j]) +
                (gtildedn[i,i] * gtildeup[j,j] * gtildeup[i,i] * gtildedn[j,j]) + 
                (gtildedn[i,i] * gdn[i,j] * gtildeup[i,i] * gtildedn[j,j]) +
                (gtildedn[i,i] * gtildeup[j,i] * gup[j,i] * gtildedn[j,j]) + 
                (gtildedn[i,j] * gdn[i,j] * gtildeup[j,i] * gup[j,i])
            end
        end
    end
    avgP = P/6/Ns
end

function saveallobser!(gup::Matrix{Float64},gdn::Matrix{Float64},l::lattice,result::Array{Float64})
    result[1] = occupy(gup,gdn,l.Ns)
    result[2] = doubleoccupy(gup,gdn,l.Ns)
    gtildeup = Diagonal(Vector(ones(l.Ns))) - transpose(gup)
    gtildedn = Diagonal(Vector(ones(l.Ns))) - transpose(gdn)
    result[3] = kinetic(gtildeup,gtildedn,l.Ns,l.Tmatrix)
    result[4] = nnspincorr(gup,gdn,gtildeup,gtildedn,l.Tmatrix)
    result[5] = Sπ(gup,gdn,gtildeup,gtildedn,l)
end