
function bitCount(n::Int)
	count = 0 
	while n>0 
		count += n & 1
		n >>= 1
	end
	return count
end

function findFirstBit(n::Int)
	n > 0 || throw(DomainError(" n needs to be larger than 0"))
	
	breakcondition = 0
	searchposition = 1
	idx = 0
	while breakcondition==0
		breakcondition = n & searchposition
		searchposition <<= 1 
		idx += 1
	end
	
	return idx
end


function getBitPos(n::Int)

	Pos = Vector{Int}([])
	while n > 0 
	
		idx = findFirstBit(n)
		push!(Pos,idx)

		n -= 2^(idx-1)

	end

	return Pos
end

"""
    bitSwap(n, i, j)

swaps bit i with j in n, for e.g.

n = 001001 , i = 1, j = 2 results in bitSwap(n,i,j) = 001010

#Arguments
- n: bit-string represented as an integer 
- i: position of the first bit which should be swaped with the second
- j: position of the second bit which should be swaped with the first

"""
function bitSwap(n::Int, i::Int, j::Int)	
	i -= 1
	j -= 1

	bit1 = (n >> i) & 1
	bit2 = (n >> j) & 1

	x = bit1 ⊻ bit2

	x = (x << i) | (x << j)

	return (n ⊻ x)
end



function FermionicBitSwap(n::Int, j::Int, k::Int)
	n_new = bitSwap(n,j,k)
	untouchedBit =  n & n_new	

	#=
	bitParticipating = n ⊻ n_new
	@show bitParticipating
	actionOnN = n & bitParticipating
	@show actionOnN
	actionOnN_new = n_new & bitParticipating
	@show actionOnN_new
	a = getBitPos(actionOnN)[1]
	@show a
	c = Int(sum(2.0.^[x for x in a:L]))
 	@show c
	@show untouchedBit & c
	jumped_C = bitCount(untouchedBit & c)
	@show jumped_C

	a = getBitPos(actionOnN_new)[1]
	@show a
	c = Int(sum(2.0.^[x for x in a:L]))
	@show c
	@show untouchedBit & c
	jumped_Cdagger = bitCount(untouchedBit & c)
	
	@show jumped_Cdagger
	
	return (-1)^(jumped_C+jumped_Cdagger)
	=#
	
	inbetween = j<k ? (j+1:k-1) : (k+1:j-1)
	c = Int(sum(2.0.^[x-1 for x in inbetween]))
	overjumped_fermions = bitCount(untouchedBit & c)
	return n_new, (-1)^overjumped_fermions
end



