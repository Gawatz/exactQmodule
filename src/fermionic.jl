
struct QFermionicBasis <: AbstractQBasis
	Dim::Int
	N::Int
	Np::Int
	OccRep::Vector{Tuple{Int,Int}}
end


function fermionicHashfunc(n::Int)
	i::Int = 1
	I::Int = 0
	while n>0
		idx = findFirstBit(n)
		n -= Int(2^(idx-1))
		
		I += binomial((idx-1),i)
		i += 1
	end
	return I
end

function QFermionicBasis(N::Int, Np::Int)
	OccRep = Vector{Tuple{Int,Int}}()

	n = [0 for i in 1:N]
	n[1:Np] = [1 for i in 1:Np]
	for x in combinations(1:N,Np)
		nidx = 0
		x = N.-x
		nidx = Int(sum(2.0.^(x)))
		#@show nidx	
		
		nReducedBasis = Int(fermionicHashfunc(nidx)+1)
		push!(OccRep,(nidx,nReducedBasis))

	end

	OccRep = OccRep[end:-1:1]
	return QFermionicBasis(size(OccRep)[1], N, Np, OccRep)
end

function fermionicCrMatrix(i::Int, N::Int, Np::Int; preQBasis::Union{AbstractQBasis, Nothing} = nothing)
	N  >= Np+1 || throw(DomainError("you cannot creat more particles $(Np+1) then sites $N"))
	if preQBasis == nothing
		preQBasis = QFermionicBasis(N, Np)
	else	
	 	N == preQBasis.N || throw(DomainError("N $N and number of sites of QBasis mismatch!"))
		Np == preQBasis.Np || throw(DomainError("Np $Np and number of particles of QBasis mismatch!"))
	end

	
	#creationMatrix = spzeros(Float64,newQBasis.Dim, QState.QBasis.Dim)

	rowVec = Vector{Int}([])
	colVec = Vector{Int}([])
	dataVec = Vector{Float64}([])
	for Idx in preQBasis.OccRep
		#@show Idx
		rowIdx = Idx[1] ⊻ (1<<(i-1)) 
		#@show rowIdx
		if rowIdx > Idx[1]
			rowIdx = Int(fermionicHashfunc(rowIdx)+1)
			#@show rowIdx
			push!(rowVec, rowIdx)
			push!(colVec, Idx[2])
			
			
			remainderBits = Idx[1] & Int(sum(2.0.^[j-1 for j = 1:(i-1)]))
			NpreParticles = bitCount(remainderBits)

			c = NpreParticles%2 == 0 ? 1 : -1
	 		push!(dataVec, c)   
	 	end
	end

	creationMatrix = sparse(rowVec, colVec, dataVec, binomial(N,Np+1), preQBasis.Dim)

	return creationMatrix
end

function fermionicCr(State::AbstractVector{<:Number}, preQBasis::QFermionicBasis, i::Int, N::Int, Np::Int)

	#@show State.Coef	
	N == preQBasis.N || throw(DomainError("system size of state does notmatch with N"))
	Np == preQBasis.Np || throw(DomainError("particle numb	er of state does not match with Np"))


	dataVec = spzeros(ComplexF64, binomial(N,Np+1))

	#if sparse I can reduce the for loop to non-zero elements	
	nzind = 0 
	try
		nzind = State.nzind
		#println("sparse")
	catch
		nzind = 1:size(State)[1]
	end

	
	for n in nzind
		Idx = preQBasis.OccRep[n]
		#@show Idx
		rowIdx = Idx[1] ⊻ (1<<(i-1)) 
		#@show rowIdx
		if rowIdx > Idx[1]   # if particle was at i rowIdx will be lower than initial Idx[1]
			rowIdx  = Int(fermionicHashfunc(rowIdx)+1)
			remainderBits = Idx[1] & Int(sum(2.0.^[j-1 for j = 1:(i-1)]))
			NpreParticles = bitCount(remainderBits)
			c = NpreParticles%2 == 0 ? 1 : -1
			dataVec[rowIdx] =  c*State[Idx[2]]
		end
	end
	#
	#State = QState(dataVec,QFermionicBasis(N,Np+1				))
	#
	return dataVec
end

fermionicCr(State::QState, i::Int, N::Int, Np::Int) = fermionicCr(State.Coef, State.QBasis, i, N, Np)


function addParticleFermionic(singleCoef::AbstractVector{<:Number}, preState::QState)
 	N = preState.QBasis.N
	Np = preState.QBasis.Np
	
	N >= Np+1 || throw(DomainError("you cannot creat more particles $(Np+1) then sites $N"))

	newQBasis = QFermionicBasis(N,Np+1)
	newCoef = zeros(ComplexF64,newQBasis.Dim)
	preQBasis = preState.QBasis
	
	for i in 1:N
		
		for Idx in preQBasis.OccRep
			
			rowIdx = Idx[1] ⊻ (1<<(i-1)) 
			
			if rowIdx > Idx[1]
				rowIdx = Int(fermionicHashfunc(rowIdx)+1)
				remainderBits = Idx[1] & Int(sum(2.0.^[j-1 for j = 1:(i-1)]))
				NpreParticles = bitCount(remainderBits)

				c = NpreParticles%2 == 0 ? 1 : -1
				newCoef[rowIdx] += c*singleCoef[i]*preState.Coef[Idx[2]]
			  
			end
		end
	
	
	end
	
	return QState(newCoef, newQBasis)
end

function constructFermionicState(blochWaveFunc::Vector{<:Any})
	@show size(blochWaveFunc[1,:,1])	
	State = QState(blochWaveFunc[1][:,1], QFermionicBasis(size(blochWaveFunc[1])[1],1))
	for i in 2:size(blochWaveFunc)[1]
		State = addParticleFermionic(blochWaveFunc[i][:,1], State)
	end

	return State
end
