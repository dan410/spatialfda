


create.empirical.basis <- function (basisfuns, rangeval = c(0, 1), params=NULL,
    dropind = NULL, quadvals = NULL, values = NULL, basisvalues = NULL, 
    names = "empirical", axes = NULL) 
{
    type <- "empirical"
    if (!is.list(basisfuns))
    	stop("basisfuns must be a list")
    nbasis <- length(basisfuns)
    if (!is.numeric(rangeval)) 
        stop("rangaval must be numeric;  class(rangeval) = ", 
            class(rangeval))
    if (length(rangeval) < 1) 
        stop("rangeval must be a numeric vector of length 2;  ", 
            "length(rangeval) = 0.")
    if (length(rangeval) == 1) {
        if (rangeval <= 0) 
            stop("rangeval a single value that is not positive:  ", 
                rangeval)
        rangeval <- c(0, rangeval)
    }
    if (length(rangeval) > 2) 
        stop("rangeval must be a vector of length 2;  ", "length(rangeval) = ", 
            length(rangeval))
    if (diff(rangeval) <= 0) 
        stop("rangeval must cover a positive range;  diff(rangeval) = ", 
            diff(rangeval))

    if (length(dropind) == 0) 
        dropind <- NULL
    if (length(dropind) > 0) {
        if (!is.numeric(dropind)) 
            stop("dropind must be numeric;  is ", class(dropind))
        doops <- which((dropind%%1) > 0)
        if (length(doops) > 0) 
            stop("dropind must be integer;  element ", doops[1], 
                " = ", dropind[doops[1]], "; fractional part = ", 
                dropind[doops[1]]%%1)
        doops0 <- which(dropind <= 0)
        if (length(doops0) > 0) 
            stop("dropind must be positive integers;  element ", 
                doops0[1], " = ", dropind[doops0[1]], " is not.")
        doops2 <- which(dropind > nbasis)
        if (length(doops2) > 0) 
            stop("dropind must not exceed nbasis = ", nbasis, 
                ";  dropind[", doops2[1], "] = ", dropind[doops2[1]])
        dropind <- sort(dropind)
        if (length(dropind) > 1) {
            if (min(diff(dropind)) == 0) 
                stop("Multiple index values in DROPIND.")
        }
    }
    params=NULL
    basisobj <- basisfd(type = type, rangeval = rangeval, nbasis = nbasis, 
        params = params, dropind = dropind, quadvals = quadvals, 
        values = values, basisvalues = basisvalues, basisfuns = basisfuns)
    {
        if (length(names) == nbasis) 
            basisobj$names <- names
        else {
            if (length(names) > 1) 
                stop("length(names) = ", length(names), ";  must be either ", 
                  "1 or nbasis = ", nbasis)
            basisobj$names <- paste(names, 0:(nbasis - 1), sep = "")
        }
    }
    if (!is.null(axes)) 
        basisobj$axes <- axes
    basisobj
}

######################################################

empiricalbasis <- function(x, basisfuns, nderiv=0){ 
	x <- as.vector(x)
	n <- length(x)
	nbasis <- length(basisfuns)
	empmat <- matrix(0, nrow = n, ncol = nbasis)
	if (nderiv == 0){
		for (ibasis in seq(nbasis)){
			empmat[,ibasis] <- basisfuns[[ibasis]](x)
		}
	}
	return(empmat)
}
