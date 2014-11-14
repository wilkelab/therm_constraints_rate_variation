
# function to calculate AICc
AICc = function(m,ndat,kextra=0){
    ktot = m$rank+kextra
    aic =  AIC(m) + 2*kextra
    aicc = aic + 2*ktot*(ktot+1)/(ndat-ktot-1) 
}


read.ddG.foldx = function( pdb.id, chain )
{
    file = paste("../ddG_calculations/foldX/foldx_ddG/", pdb.id, "_", chain, "_foldx_ddG.txt", sep='' )
    if (file.exists(file))
        ddG = read.table( file, header=T )
    else
        ddG = NULL
    ddG
}

read.ddG.rosetta = function( pdb.id, chain )
{
    file = paste("../ddG_calculations/rosetta/rosetta_ddG/", pdb.id, "_", chain, "_rosetta_ddG.txt", sep='' )
    if (file.exists(file))
        ddG = read.table( file, header=T )
    else
        ddG = NULL
    ddG
}
