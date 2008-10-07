.First.lib <- function(lib, pkg){
    library.dynam("bisoreg", pkg, lib)
}

.Last.lib <- function(lib){
    library.dynam.unload("bisoreg", lib)
    ##detach(package:coda)
    ##detach(package:R2WinBUGS)
    ##detach(package:monreg)
    ##detach(package:bootstrap)
}
