fdParcheck <- function (fdParobj) 
{
    if (!inherits(fdParobj, "fdPar")) {
        if (inherits(fdParobj, "fd") || inherits(fdParobj, "basisfd")) 
            fdParobj <- fdPar(fdParobj)
        else stop(paste("'fdParobj' is not a functional parameter object,", 
            "not a functional data object, and", "not a basis object."))
    }
    return(fdParobj)
}


##################################################################

fdPar <- function (fdobj = NULL, Lfdobj = NULL, lambda = 0, estimate = TRUE, 
    penmat = NULL) 
{
    if (!inherits(fdobj, "fd")) {
        if (is.null(fdobj)) {
            fdobj = fd()
        }
        else {
            if (inherits(fdobj, "basisfd")) {
                nbasis <- fdobj$nbasis
                dropind <- fdobj$dropind
                coefs <- matrix(0, nbasis - length(dropind), 
                  1)
                fdnames <- list("time", "reps 1", "values")
                if (!is.null(fdobj$names)) {
                  nms <- {
                    if (length(dropind) > 1) 
                      fdobj$names[-dropind]
                    else fdobj$names
                  }
                  dimnames(coefs) <- list(nms, NULL)
                  fdnames[[1]] <- nms
                }
                fdobj <- fd(coefs, fdobj)
            }
            else if (is.numeric(fdobj)) 
                fdobj <- fd(fdobj)
            else stop("First argument is neither a functional data object nor a basis object.")
        }
    }
    {
        if (is.null(Lfdobj)) {
            if (fdobj$basis$type == "fourier") {
                rng <- fdobj$basis$rangeval
                Lfdobj <- vec2Lfd(c(0, (2 * pi/diff(rng))^2, 
                  0), rng)
            }
            else {
                norder <- {
                  if (fdobj$basis$type == "bspline") 
                    norder.bspline(fdobj$basis)
                  else 2
                }
                Lfdobj <- int2Lfd(max(0, norder - 2))
            }
        }
        else Lfdobj <- int2Lfd(Lfdobj)
    }
    if (!inherits(Lfdobj, "Lfd")) 
        stop("'Lfdobj' is not a linear differential operator object.")
    if (!is.numeric(lambda)) 
        stop("Class of LAMBDA is not numeric.")
    if (lambda < 0) 
        stop("LAMBDA is negative.")
    if (!is.logical(estimate)) 
        stop("Class of ESTIMATE is not logical.")
    if (!is.null(penmat)) {
        if (!is.numeric(penmat)) 
            stop("PENMAT is not numeric.")
        penmatsize <- dim(penmat)
        if (any(penmatsize != nbasis)) 
            stop("Dimensions of PENMAT are not correct.")
    }
    fdParobj <- list(fd = fdobj, Lfd = Lfdobj, lambda = lambda, 
        estimate = estimate, penmat = penmat)
    oldClass(fdParobj) <- "fdPar"
    fdParobj
}

#################################################################

fd <- function (coef = NULL, basisobj = NULL, fdnames = NULL) 
{
    if (is.null(coef) && is.null(basisobj)) 
        basisobj <- basisfd()
    if (is.null(coef)) 
        coef <- rep(0, basisobj[["nbasis"]])
    type <- basisobj$type
    {
        if (!is.numeric(coef)) 
            stop("'coef' is not numeric.")
        else if (is.vector(coef)) {
            coef <- as.matrix(coef)
            if (identical(type, "constant")) 
                coef <- t(coef)
            coefd <- dim(coef)
            ndim <- length(coefd)
        }
        else if (is.matrix(coef)) {
            coefd <- dim(coef)
            ndim <- length(coefd)
        }
        else if (is.array(coef)) {
            coefd <- dim(coef)
            ndim <- length(coefd)
        }
        else stop("Type of 'coef' is not correct")
    }
    if (ndim > 3) 
        stop("'coef' not of dimension 1, 2 or 3")
    {
        if (is.null(basisobj)) {
            rc <- range(coef)
            if (diff(rc) == 0) 
                rc <- rc + 0:1
            dimC <- dim(coef)
            nb <- {
                if (is.null(dimC)) 
                  length(coef)
                else dimC[1]
            }
            basisobj <- create.bspline.basis(rc, nbasis = max(4, 
                nb))
            type <- basisobj$type
        }
        else if (!(inherits(basisobj, "basisfd"))) 
            stop("Argument basis must be of basis class")
    }
    nbasis = basisobj$nbasis
    if (coefd[1] != nbasis) 
        stop("First dim. of 'coef' not equal to 'nbasis'.")
    dropind <- basisobj$dropind
    if (length(dropind) > 0) 
        stop("'fd' not yet programmed to handle 'dropind'")
    if (coefd[1] != basisobj$nbasis) 
        stop("Number of coefficients does not match ", "the number of basis functions.")
    if (ndim > 1) 
        nrep <- coefd[2]
    else nrep <- 1
    if (ndim > 2) 
        nvar <- coefd[3]
    else nvar <- 1
    if (is.null(fdnames)) {
        if (ndim == 1) 
            fdnames <- list("time", "reps", "values")
        if (ndim == 2) 
            fdnames <- list("time", paste("reps", as.character(1:nrep)), 
                "values")
        if (ndim == 3) 
            fdnames <- list("time", paste("reps", as.character(1:nrep)), 
                paste("values", as.character(1:nvar)))
        names(fdnames) <- c("args", "reps", "funs")
    }
    if (is.null(dimnames(coef))) {
        dimc <- dim(coef)
        ndim <- length(dimc)
        dnms <- vector("list", ndim)
        if (dimc[1] == length(fdnames[[1]])) 
            dnms[[1]] <- fdnames[[1]]
        if ((ndim > 1) && (dimc[2] == length(fdnames[[2]]))) 
            dnms[[2]] <- fdnames[[2]]
        if ((ndim > 2) && (dimc[3] == length(fdnames[[3]]))) 
            dnms[[3]] <- fdnames[[3]]
        if (!all(sapply(dnms, is.null))) 
            dimnames(coef) <- dnms
    }
    fdobj <- list(coefs = coef, basis = basisobj, fdnames = fdnames)
    oldClass(fdobj) <- "fd"
    fdobj
}

