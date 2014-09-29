eval.fd <- function (evalarg, fdobj, Lfdobj = 0, returnMatrix = FALSE) 
{
    Lfdobj <- int2Lfd(Lfdobj)
    if (inherits(fdobj, "numeric") && inherits(evalarg, "fd")) {
        temp <- fdobj
        fdobj <- evalarg
        evalarg <- temp
    }
    if (!(is.numeric(evalarg))) 
        stop("Argument EVALARG is not numeric.")
    evaldim <- dim(evalarg)
    if (!(length(evaldim) < 3)) 
        stop("Argument EVALARG is not a vector or a matrix.")
    if (!(inherits(fdobj, "fd"))) 
        stop("Argument FD is not a functional data object.")
    basisobj <- fdobj$basis
    nbasis <- basisobj$nbasis
    rangeval <- basisobj$rangeval
    onerow <- rep(1, nbasis)
    temp <- c(evalarg)
    temp <- temp[!(is.na(temp))]
    EPS <- 1e-14
    if (min(temp) < rangeval[1] - EPS || max(temp) > rangeval[2] + 
        EPS) {
        warning(paste("Values in argument EVALARG are outside of permitted range,", 
            "and will be ignored."))
        print(c(rangeval[1] - min(temp), max(temp) - rangeval[2]))
    }
    if (is.vector(evalarg)) {
        n <- length(evalarg)
    }
    else {
        n <- evaldim[1]
    }
    coef <- fdobj$coefs
    coefd <- dim(coef)
    ndim <- length(coefd)
    if (ndim <= 1) 
        nrep <- 1
    else nrep <- coefd[2]
    if (ndim <= 2) 
        nvar <- 1
    else nvar <- coefd[3]
    if (ndim <= 2) 
        evalarray <- matrix(0, n, nrep)
    else evalarray <- array(0, c(n, nrep, nvar))
    if (ndim == 2) 
        dimnames(evalarray) <- list(NULL, dimnames(coef)[[2]])
    if (ndim == 3) 
        dimnames(evalarray) <- list(NULL, dimnames(coef)[[2]], 
            dimnames(coef)[[3]])
    if (is.vector(evalarg)) {
        evalarg[evalarg < rangeval[1] - 1e-10] <- NA
        evalarg[evalarg > rangeval[2] + 1e-10] <- NA
        basismat <- eval.basis(evalarg, basisobj, Lfdobj, returnMatrix)
        if (ndim <= 2) {
            evalarray <- basismat %*% coef
            dimnames(evalarray) <- list(rownames(basismat), colnames(coef))
        }
        else {
            evalarray <- array(0, c(n, nrep, nvar))
            for (ivar in 1:nvar) evalarray[, , ivar] <- basismat %*% 
                coef[, , ivar]
        }
    }
    else {
        for (i in 1:nrep) {
            evalargi <- evalarg[, i]
            if (all(is.na(evalargi))) 
                stop(paste("All values are NA for replication", 
                  i))
            index <- !(is.na(evalargi) | evalargi < rangeval[1] | 
                evalargi > rangeval[2])
            evalargi <- evalargi[index]
            basismat <- eval.basis(evalargi, basisobj, Lfdobj, 
                returnMatrix)
            if (ndim == 2) {
                evalarray[index, i] <- as.vector(basismat %*% 
                  coef[, i])
                evalarray[!(index), i] <- NA
            }
            if (ndim == 3) {
                for (ivar in 1:nvar) {
                  evalarray[index, i, ivar] <- as.vector(basismat %*% 
                    coef[, i, ivar])
                  evalarray[!(index), i, ivar] <- NA
                }
            }
        }
    }
    if ((length(dim(evalarray)) == 2) && !returnMatrix) {
        return(as.matrix(evalarray))
    }
    else return(evalarray)
}
