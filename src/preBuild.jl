using JLD

currentDir = dirname(@__FILE__)
@show currentDir
@show findlast("/",currentDir)[1]
@show dirname(@__FILE__)
setDir = string(currentDir[1:findlast("/",currentDir)[1]],"preBuild/")

function loadFermionicCr(i::Int, N::Int, Np::Int)
	
	cr = load(string(string(setDir,"Creation/Fermionic/"),"cr_",string(i),"_N",string(N),"_Np",string(Np),".jld"))
	cr = cr["creation"]

	return cr
end

function saveFermionicCr(cr::AbstractArray{<:Number,2},i::Int, N::Int, Np::Int)
	save(string(string(setDir,"Creation/Fermionic/"),"cr_",string(i),"_N",string(N),"_Np",string(Np),".jld"),"creation",cr)
end


#=
function loadFermionicFT()
end

function saveFermionicFT()
end
=#
