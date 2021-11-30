GPF <- list()

.onLoad <- function(libname, pkgname)
  {
    cat(sprintf("Welcome to %s (vers 1.3) !", pkgname))

    GPF$ID <<- initID()
    GPF$UTIL <<- initUTIL()
    GPF$PF <<- initPF()
    GPF$GA <<- initGA()
  }
