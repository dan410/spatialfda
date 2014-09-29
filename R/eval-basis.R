


eval.basis <- function (evalarg, basisobj, Lfdobj = 0, returnMatrix = FALSE) 
{
    if (is.numeric(basisobj) && inherits(evalarg, "basisfd")) {
        temp <- basisobj
        basisobj <- evalarg
        evalarg <- temp
    }
    if (!(is.numeric(evalarg))) 
        stop("Argument EVALARG is not numeric.")
    if (!is.vector(evalarg)) 
        stop("Argument EVALARG is not a vector.")
    if (!(inherits(basisobj, "basisfd"))) 
        stop("Second argument is not a basis object.")
    Lfdobj <- int2Lfd(Lfdobj)
    nderiv <- Lfdobj$nderiv
    bwtlist <- Lfdobj$bwtlist
    basismat <- getbasismatrix(evalarg, basisobj, nderiv, returnMatrix)
    nbasis <- dim(basismat)[2]
    oneb <- matrix(1, 1, nbasis)
    if (nderiv > 0) {
        nonintwrd <- FALSE
        for (j in 1:nderiv) {
            bfd <- bwtlist[[j]]
            bbasis <- bfd$basis
            if (bbasis$type != "constant" || bfd$coefs != 0) 
                nonintwrd <- TRUE
        }
        if (nonintwrd) {
            for (j in 1:nderiv) {
                bfd <- bwtlist[[j]]
                if (!all(c(bfd$coefs) == 0)) {
                  wjarray <- eval.fd(evalarg, bfd, 0, returnMatrix)
                  Dbasismat <- getbasismatrix(evalarg, basisobj, 
                    j - 1, returnMatrix)
                  basismat <- basismat + (wjarray %*% oneb) * 
                    Dbasismat
                }
            }
        }
    }
    if ((!returnMatrix) && (length(dim(basismat)) == 2)) {
        return(as.matrix(basismat))
    }
    return(basismat)
}

