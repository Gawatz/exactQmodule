#
#	Quantum State struct
#
struct QState{T<:Number}
	Coef::AbstractVector{T}
	QBasis::Union{<:AbstractQBasis}
end





#
#	basic function on Quantum States
#

function getCijMatrix(Coef::Vector{<:Number}, N::Int, Np::Int; Basis::Union{AbstractQBasis, Nothing} = nothing)

	CijMatrix = Array{ComplexF64,2}(undef, N, N)

	for i in 1:N
		cr_i = nothing 
		try
			cr_i = loadFermionicCr(i, N, Np)
		catch	
			cr_i = fermionicCrMatrix(i, N, Np; preQBasis=Basis)
			saveFermionicCr(cr_i, i, N, Np)
		end
		for j in i:N
			cr_j = nothing 
			try
				cr_j = loadFermionicCr(j, N, Np)
			catch	
				cr_j = fermionicCrMatrix(j, N, Np; preQBasis=Basis)
				saveFermionicCr(cr_j, j, N, Np)
			end
			
			#cr_j = fermionicCrMatrix(j, N, Np; preQBasis = Basis)
			C = cr_i*transpose(cr_j)
			CijMatrix[i,j] = transpose(conj(Coef))*C*Coef
			if i != j
				CijMatrix[j,i] = conj(CijMatrix[i,j])
			end
		end
	end

	return CijMatrix
end

getCijMatrix(State::QState) = getCijMatrix(State.Coef, State.QBasis.N, State.QBasis.Np-1)


function projectOntoStates(States::Vector{<:Any}, CijMatrix::AbstractArray{<:Number})

	occ = 0
	for x in States
		occ += (transpose(x)*CijMatrix*conj(x))[1,1]
	end

	return occ

end

function entanglementHalfCut(Coef::Vector{<:Number}, N::Int, Np::Int; Basis::Union{AbstractQBasis, Nothing} = nothing)

	if Basis == nothing
	
		Basis = QFermionicBasis(N, Np)
	end

	dim = Int(sum(2.0.^[i for i = 0:Int(N/2-1)]))+1
	CoefMatrix = spzeros(ComplexF64,dim,dim)

	for idx in Basis.OccRep
		
		nleft = idx[1] >> Int(N/2)
		nright = idx[1] - (nleft << Int(N/2))+1
		nleft += 1
		#@show idx[2], idx[1], nleft, nright, Coef[idx[2]]
		#if I want to do symsec I have to count bits here
		

		CoefMatrix[nleft,nright] = Coef[idx[2]]
	
	end

	#@show CoefMatrix

	quality = 0.0 
	k = 10
	svdHandle = []
	svdHandle2 = []
	s = []
	while quality < 0.99
		try 
			svdHandle = svds(CoefMatrix, nsv = k)
			svdHandle2 = svdsolve(CoefMatrix, k)
		catch
			@show "something went wrong"
		end
		quality = sum(svdHandle[1].S.^2.0)
		k += 10
		s = svdHandle[1].S[:]

	end

	@show s
	@show svdHandle2[1]
	@show max(real.((CoefMatrix - sparse(svdHandle[1].U*diagm(0 => s)*transpose(conj(svdHandle[1].V)))).nzval)...)

	@show max(imag.((CoefMatrix - sparse(svdHandle[1].U*diagm(0 => s)*transpose(conj(svdHandle[1].V)))).nzval)...)
	return s, svdHandle[1].U, svdHandle[1].V, quality
end


entanglementHalfCut(State::QState) = entanglementHalfCut(State.Coef, State.QBasis.N, State.QBasis.Np, Basis = State.QBasis)
