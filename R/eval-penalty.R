

eval.penalty <- function (basisobj, Lfdobj = int2Lfd(0), rng = rangeval) 
{
    if (inherits(basisobj, "fd")) 
        basisobj <- basisobj$basis
    if (inherits(basisobj, "fdPar")) {
        fdobj <- basisobj$fd
        basisobj <- fdobj$basis
    }
    if (!inherits(basisobj, "basisfd")) 
        stop("Argument BASISOBJ is not a functional basis object.")
    rangeval <- basisobj$rangeval
    Lfdobj <- int2Lfd(Lfdobj)
    type <- basisobj$type
    if (type == "bspline") 
        penaltymat <- bsplinepen(basisobj, Lfdobj, rng)
    else if (type == "const") {
        rangeval <- getbasisrange(basisobj)
        penaltymat <- rangeval[2] - rangeval[1]
    }
    else if (type == "expon") 
        penaltymat <- exponpen(basisobj, Lfdobj)
    else if (type == "fourier") 
        penaltymat <- fourierpen(basisobj, Lfdobj)
    else if (type == "monom") 
        penaltymat <- monomialpen(basisobj, Lfdobj)
    else if (type == "polyg") 
        penaltymat <- polygpen(basisobj, Lfdobj)
    else if (type == "power") 
        penaltymat <- powerpen(basisobj, Lfdobj)
    else if (type == "empirical")
        penaltymat <- empiricalpen(basisobj, Lfdobj)
    else stop("Basis type not recognizable, can not find penalty matrix")
    dropind <- basisobj$dropind
    nbasis <- basisobj$nbasis
    if (length(dropind) > 0) {
        index <- 1:nbasis
        for (i in 1:length(dropind)) index <- index(index != 
            dropind[i])
        penaltymat <- penaltymat(index, index)
    }
    penaltymat <- (penaltymat + t(penaltymat))/2
    return(penaltymat)
}

#########################################################
empiricalpen <- function(basisobj, Lfdobj){
	cat("need to create empirical penalty function")
}
