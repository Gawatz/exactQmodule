module exactQmodule

using Base, SparseArrays, Arpack, LinearAlgebra
using Combinatorics, KrylovKit

export AbstractQBasis

#from (1)
export bitSwap, bitCount, findFirstBit, getBitPos, FermionicBitSwap 

#from (2)
export QState, getCijMatrix, projectOntoStates, entanglementHalfCut

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
