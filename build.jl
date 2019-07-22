
function main()
	currentDir = pwd()
	mkpath(string(currentDir,"/preBuild/QBasis/Fermionic"))
	mkpath(string(currentDir,"/preBuild/Creation/Fermionic"))
	mkpath(string(currentDir,"/preBuild/FT/Fermionic"))
	
	return currentDir
end


main()

