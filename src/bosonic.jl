
struct QBosonicBasis <: AbstractQBasis
  	Dim::Int 
	N::Int
	Np::Int
	OccRep::Vector{Vector{Int}}
end

function QBosonicBasis(N::Int, Np::Int)	
	I = 1
	n = zeros(Int,N)
	n[1] = Np
	
	vectors = Vector{Vector{Int}}([])
	while true
		#@show n
		push!(vectors,n[:])
		I += 1
		if n[1] > 0
			n[1] -= 1
			n[2] += 1
		elseif n[1] == 0
			j = findfirst(!iszero, n)
			
			if j == N break end
		
			n[1] = n[j] - 1
			n[j] = 0
			n[j+1] += 1

		end
		#@show n
	end

	return QBosonicBasis(I-1, N, Np, vectors)
end


function bosonicHashfunc(Np::Int, col::Vector{Int})
	hashkey = 0
	particlei = 1
	p_i = findfirst(!iszero,col)
	col = col[p_i:end]	
	while true
		
		#@show col
		
		if col[1] > 0
			col[1] -= 1
		else
 			orbital = findfirst(!iszero,col)
			p_i += orbital-1
			col = col[orbital:end]
			col[1] -= 1
		end

		hashkey += binomial((p_i-1) + particlei -1 , particlei)
	
		
		if particlei == Np break end
 		particlei += 1
	end
	return  hashkey
end


function addParticleBosonic(singleCoef::Vector{Float64}, preState::QState)
 	N = preState.QBasis.N
	Np = preState.QBasis.Np
	
	N >= Np+1 || throw(DomainError("you cannot creat more particles $(Np+1) then sites $N"))

	newQBasis = QBosonicBasis(N,Np+1)
	newCoef = zeros(ComplexF64,newQBasis.Dim)
	preQBasis = preState.QBasis
	
	for i in 1:N
	 	for (colIdx, colOccRep) in enumerate(preQBasis.OccRep)
			#@show colOccRep
	 		colOccRep[i] += 1

			# actualy this bosonicHashfunc gives directly idx in
			# OccRep vector so no dict needed
			rowIdx = bosonicHashfunc(Np+1, colOccRep)+1
			#@show i,rowIdx
			newCoef[rowIdx] +=  singleCoef[i]*preState.Coef[colIdx]*sqrt(colOccRep[i])
		
			#@show newCoef[rowIdx]
			colOccRep[i] -= 1
		end
		@show newCoef
	end
	
	return QState(newCoef, newQBasis)
end


function constructBosonicState(blochWaveFunc::Vector{<:Any})
	@show size(blochWaveFunc[1,:,1])	
	State = QState(blochWaveFunc[1][:,1], QBosonicBasis(size(blochWaveFunc[1])[1],1))
	for i in 2:size(blochWaveFunc)[1]
		State = addParticleBosonic(blochWaveFunc[i][:,1], State)
	end

	return State
end
