module exactQmodule

using Base, SparseArrays, Arpack, LinearAlgebra
using Combinatorics, KrylovKit, JLD

export AbstractQBasis

#from (1)
export bitSwap, bitCount, findFirstBit, getBitPos, FermionicBitSwap 

#from (2)
export QState, shiftStateBasis, getCijMatrix, projectOntoStates, entanglementHalfCut

#from (3)
export QBosonicBasis, addParticleBosonic, constructBosonicState, bosonicHashfunc

#from 4
export fermionicCr, QFermionicBasis, fermionicCrMatrix, addParticleFermionic, constructFermionicState, fermionicHashfunc

#from 5
export loadFermionicCr, saveFermionicCr


abstract type AbstractQBasis end

include("./bitOperations.jl")
include("./Qstate.jl")
include("./bosonic.jl")
include("./fermionic.jl")
include("./preBuild.jl")

end #modlue
