smooth.basis <- function (argvals = 1:n, y, fdParobj, wtvec = NULL, fdnames = NULL, 
    covariates = NULL, method = "chol", dfscale = 1, returnMatrix = FALSE) 
{
    if (!is.numeric(y)) 
        stop("'y' is not numeric.")
    if (is.vector(y)) 
        y <- as.matrix(y)
    dimy <- dim(y)
    ndy <- length(dimy)
    n <- dimy[1]
    if (!is.numeric(argvals)) 
        stop("'argvals' is not numeric.")
    if (is.vector(argvals)) 
        argvals <- as.matrix(argvals)
    dima <- dim(argvals)
    nda <- length(dima)
    if (ndy < nda) 
        stop("argvals has ", nda, " dimensions  y has only ", 
            ndy)
    if (nda < 3) {
        if (dima[2] == 1) {
            sb2 <- smooth.basis1(argvals, y, fdParobj, wtvec = wtvec, 
                fdnames = fdnames, covariates = covariates, method = method, 
                dfscale = dfscale, returnMatrix = returnMatrix)
        }
        else {
            sb2 <- smooth.basis2(argvals, y = y, fdParobj = fdParobj, 
                wtvec = wtvec, fdnames = fdnames, covariates = covariates, 
                method = method, dfscale = dfscale, returnMatrix = returnMatrix)
        }
        return(sb2)
    }
    if (nda < 4) {
        return(smooth.basis3(argvals, y = y, fdParobj = fdParobj, 
            wtvec = wtvec, fdnames = fdnames, covariates = covariates, 
            method = method, dfscale = dfscale, returnMatrix = returnMatrix))
    }
    else {
        cat("dim(argvals) =", paste(dima, collapse = ", "), "\n")
        cat("dim(y)      = ", paste(dimy, collapse = ", "), "\n")
        stop("Dimensions of argvals do not match those of y")
        return()
    }
}

#'
#'
#'

 smooth.basis1 <- function (argvals = 1:n, y, fdParobj, wtvec = NULL, fdnames = NULL, 
    covariates = NULL, method = "chol", dfscale = 1, returnMatrix = FALSE) 
{
    if (is.vector(y)) {
        y <- matrix(y, length(y), 1)
    }
    dimy <- dim(y)
    n <- dimy[1]
    ycheck <- ycheck(y, n)
    y <- ycheck$y
    y0 <- y
    nrep <- ycheck$ncurve
    nvar <- ycheck$nvar
    ndim <- ycheck$ndim
    ydim <- dim(y)
    if (!is.numeric(argvals)) 
        stop("'argvals' is not numeric.")
    argvals <- as.vector(argvals)
    fdParobj <- fdParcheck(fdParobj)
    fdobj <- fdParobj$fd
    lambda <- fdParobj$lambda
    Lfdobj <- fdParobj$Lfd
    penmat <- fdParobj$penmat
    if (lambda < 0) {
        warning("Value of 'lambda' was negative  0 used instead.")
        lambda <- 0
    }
    wtlist <- wtcheck(n, wtvec)
    wtvec <- wtlist$wtvec
    onewt <- wtlist$onewt
    matwt <- wtlist$matwt
    nderiv <- Lfdobj$nderiv
    basisobj <- fdobj$basis
    nbasis <- basisobj$nbasis
    onebasis <- rep(1, nbasis)
    tnames <- dimnames(y)[[1]]
    if (is.null(tnames)) 
        tnames <- 1:n
    bnames <- fdParobj$fd$basis$names
    if (ndim == 2) {
        coef <- matrix(0, nbasis, nrep)
        ynames <- dimnames(y)[[2]]
        vnames <- "value"
        dimnames(coef) <- list(bnames, ynames)
    }
    if (ndim == 3) {
        coef <- array(0, c(nbasis, nrep, nvar))
        ynames <- dimnames(y)[[2]]
        vnames <- dimnames(y)[[3]]
        dimnames(coef) <- list(bnames, ynames, vnames)
    }
    if (!is.null(covariates)) {
        if (!is.numeric(covariates)) {
            stop(paste("smooth_basis_LS:covariates", "Optional argument COVARIATES is not numeric."))
        }
        if (dim(covariates)[1] != n) {
            stop(paste("smooth_basis_LS:covariates", "Optional argument COVARIATES has incorrect number of rows."))
        }
        q <- dim(covariates)[2]
    }
    else {
        q <- 0
        beta. <- NULL
    }
    basismat <- eval.basis(as.vector(argvals), basisobj, 0, returnMatrix)
    if (method == "chol") {
        if (n > nbasis + q || lambda > 0) {
            if (!is.null(covariates)) {
                ind1 <- 1:n
                ind2 <- (nbasis + 1):(nbasis + q)
                basismat <- as.matrix(basismat)
                basismat <- cbind(basismat, matrix(0, dim(basismat)[1], 
                  q))
                basismat[ind1, ind2] <- covariates
            }
            if (matwt) {
                wtfac <- chol(wtvec)
                basisw <- wtvec %*% basismat
            }
            else {
                rtwtvec <- sqrt(wtvec)
                rtwtmat <- matrix(rtwtvec, n, nrep)
                basisw <- (wtvec %*% matrix(1, 1, nbasis + q)) * 
                  basismat
            }
            Bmat <- t(basisw) %*% basismat
            Bmat0 <- Bmat
            if (ndim < 3) {
                Dmat <- t(basisw) %*% y
            }
            else {
                Dmat <- array(0, c(nbasis + q, nrep, nvar))
                for (ivar in 1:nvar) {
                  Dmat[, , ivar] <- crossprod(basisw, y[, , ivar])
                }
            }
            if (lambda > 0) {
                if (is.null(penmat)) 
                  penmat <- eval.penalty(basisobj, Lfdobj)
                Bnorm <- sqrt(sum(diag(t(Bmat0) %*% Bmat0)))
                pennorm <- sqrt(sum(penmat * penmat))
                condno <- pennorm/Bnorm
                if (lambda * condno > 1e+12) {
                  lambda <- 1e+12/condno
                  warning(paste("lambda reduced to", lambda, 
                    "to prevent overflow"))
                }
                if (!is.null(covariates)) {
                  penmat <- rbind(cbind(penmat, matrix(0, nbasis, 
                    q)), cbind(matrix(0, q, nbasis), matrix(0, 
                    q, q)))
                }
                Bmat <- Bmat0 + lambda * penmat
            }
            else {
                penmat <- NULL
                Bmat <- Bmat0
            }
            Bmat <- (Bmat + t(Bmat))/2
            Lmat <- try(chol(Bmat), silent = TRUE)
            if (class(Lmat) == "try-error") {
                Beig <- eigen(Bmat, symmetric = TRUE)
                BgoodEig <- (Beig$values > 0)
                Brank <- sum(BgoodEig)
                if (Brank < dim(Bmat)[1]) 
                  warning("Matrix of basis function values has rank ", 
                    Brank, " < dim(fdobj$basis)[2] = ", length(BgoodEig), 
                    "  ignoring null space")
                goodVec <- Beig$vectors[, BgoodEig]
                Bmatinv <- (goodVec %*% (Beig$values[BgoodEig] * 
                  t(goodVec)))
            }
            else {
                Lmatinv <- solve(Lmat)
                Bmatinv <- Lmatinv %*% t(Lmatinv)
            }
            if (ndim < 3) {
                coef <- Bmatinv %*% Dmat
                if (!is.null(covariates)) {
                  beta. <- as.matrix(coef[(nbasis + 1):(nbasis + 
                    q), ])
                  coef <- as.matrix(coef[1:nbasis, ])
                }
                else {
                  beta. <- NULL
                }
            }
            else {
                coef <- array(0, c(nbasis, nrep, nvar))
                if (!is.null(covariates)) {
                  beta. <- array(0, c(q, nrep, nvar))
                }
                else {
                  beta. <- NULL
                }
                for (ivar in 1:nvar) {
                  coefi <- Bmatinv %*% Dmat[, , ivar]
                  if (!is.null(covariates)) {
                    beta.[, , ivar] <- coefi[(nbasis + 1):(nbasis + 
                      q), ]
                    coef[, , ivar] <- coefi[1:nbasis, ]
                  }
                  else {
                    coef[, , ivar] <- coefi
                  }
                }
            }
        }
        else {
            if (n == nbasis + q) {
                if (ndim == 2) 
                  coef <- solve(basismat, y)
                else for (ivar in 1:var) coef[1:n, , ivar] <- solve(basismat, 
                  y[, , ivar])
                penmat <- NULL
            }
            else {
                stop("The number of basis functions = ", nbasis + 
                  q, " exceeds ", n, " = the number of points to be smoothed.")
            }
        }
    }
    else {
        if (n > nbasis || lambda > 0) {
            if (!onewt) {
                if (matwt) {
                  wtfac <- chol(wtvec)
                  basismat.aug <- wtfac %*% basismat
                  if (ndim < 3) {
                    y <- wtfac %*% y
                  }
                  else {
                    for (ivar in 1:nvar) {
                      y[, , ivar] <- wtfac %*% y[, , ivar]
                    }
                  }
                }
                else {
                  rtwtvec <- sqrt(wtvec)
                  basismat.aug <- matrix(rtwtvec, n, nbasis) * 
                    basismat
                  if (ndim < 3) {
                    y <- matrix(rtwtvec, n, nrep) * y
                  }
                  else {
                    for (ivar in 1:nvar) {
                      y[, , ivar] <- matrix(rtwtvec, n, nrep) * 
                        y[, , ivar]
                    }
                  }
                }
            }
            else {
                basismat.aug <- basismat
            }
            if (lambda > 0) {
                if (is.null(penmat)) 
                  penmat <- eval.penalty(basisobj, Lfdobj)
                eiglist <- eigen(penmat)
                Dvec <- eiglist$values
                Vmat <- eiglist$vectors
                neiglow <- nbasis - nderiv
                naug <- n + neiglow
                if (Dvec[neiglow] <= 0) {
                  stop(paste("smooth_basis:eig", "Eigenvalue(NBASIS-NDERIV) of penalty matrix ", 
                    "is not positive check penalty matrix."))
                }
                indeig <- 1:neiglow
                penfac <- Vmat[, indeig] %*% diag(sqrt(Dvec[indeig]))
                basismat.aug <- rbind(basismat.aug, sqrt(lambda) * 
                  t(penfac))
                if (ndim < 3) {
                  y <- rbind(y, matrix(0, nbasis - nderiv, nrep))
                }
                else {
                  y <- array(0, c(ydim[1] + nbasis - nderiv, 
                    ydim[2], ydim[3]))
                  y[1:ydim[1], , ] <- y0
                  ind1 <- (1:(nbasis - nderiv)) + ydim[1]
                  for (ivar in 1:nvar) {
                    y[ind1, , ivar] <- matrix(0, nbasis - nderiv, 
                      nrep)
                  }
                }
            }
            else {
                penmat <- NULL
            }
            if (!is.null(covariates)) {
                ind1 <- 1:n
                ind2 <- (nbasis + 1):(nbasis + q)
                basismat.aug <- cbind(basismat.aug, matrix(0, 
                  naug, q))
                if (!onewt) {
                  if (matwt) {
                    basismat.aug[ind1, ind2] <- wtfac %*% covariates
                  }
                  else {
                    wtfac <- matrix(rtwtvec, n, q)
                    basismat.aug[ind1, ind2] <- wtfac * covariates
                  }
                }
                else {
                  basismat.aug[ind1, ind2] <- covariates
                }
                penmat <- rbind(cbind(penmat, matrix(0, nbasis, 
                  q)), cbind(matrix(0, q, nbasis), matrix(0, 
                  q, q)))
            }
            qr <- qr(basismat.aug)
            if (ndim < 3) {
                coef <- qr.coef(qr, y)
                if (!is.null(covariates)) {
                  beta. <- coef[ind2, ]
                  coef <- coef[1:nbasis, ]
                }
                else {
                  beta. <- NULL
                }
            }
            else {
                coef <- array(0, c(nbasis, nrep, nvar))
                if (!is.null(covariates)) {
                  beta. <- array(0, c(q, nrep, nvar))
                }
                else {
                  beta. <- NULL
                }
                for (ivar in 1:nvar) {
                  coefi <- qr.coef(qr, y[, , ivar])
                  if (!is.null(covariates)) {
                    beta.[, , ivar] <- coefi[ind2, ]
                    coef[, , ivar] <- coefi[1:nbasis, ]
                  }
                  else {
                    coef[, , ivar] <- coefi
                  }
                }
            }
        }
        else {
            if (n == nbasis + q) {
                {
                  if (ndim == 2) 
                    coef <- solve(basismat, y)
                  else for (ivar in 1:var) coef[, , ivar] <- solve(basismat, 
                    y[, , ivar])
                }
                penmat <- NULL
            }
            else {
                stop(paste("The number of basis functions = ", 
                  nbasis, " exceeds ", n, " = the number of points to be smoothed.  "))
            }
        }
    }
    if (onewt) {
        temp <- t(basismat) %*% basismat
        if (lambda > 0) {
            temp <- temp + lambda * penmat
        }
        L <- chol(temp)
        MapFac <- solve(t(L), t(basismat))
        y2cMap <- solve(L, MapFac)
    }
    else {
        if (matwt) {
            temp <- t(basismat) %*% wtvec %*% basismat
        }
        else {
            temp <- t(basismat) %*% (as.vector(wtvec) * basismat)
        }
        if (lambda > 0) {
            temp <- temp + lambda * penmat
        }
        L <- chol((temp + t(temp))/2)
        MapFac <- solve(t(L), t(basismat))
        if (matwt) {
            y2cMap <- solve(L, MapFac %*% wtvec)
        }
        else {
            y2cMap <- solve(L, MapFac * rep(as.vector(wtvec), 
                e = nrow(MapFac)))
        }
    }
    df. <- sum(diag(y2cMap %*% basismat))
    if (ndim < 3) {
        yhat <- basismat[, 1:nbasis] %*% coef
        SSE <- sum((y[1:n, ] - yhat)^2)
        if (is.null(ynames)) 
            ynames <- dimnames(yhat)[[2]]
    }
    else {
        SSE <- 0
        yhat <- array(0, c(n, nrep, nvar))
        dimnames(yhat) <- list(dimnames(basismat)[[1]], dimnames(coef)[[2]], 
            dimnames(coef)[[3]])
        for (ivar in 1:nvar) {
            yhat[, , ivar] <- basismat[, 1:nbasis] %*% coef[, 
                , ivar]
            SSE <- SSE + sum((y[1:n, , ivar] - yhat[, , ivar])^2)
        }
        if (is.null(ynames)) 
            ynames <- dimnames(yhat)[[2]]
        if (is.null(vnames)) 
            vnames <- dimnames(yhat)[[2]]
    }
    if (is.null(ynames)) 
        ynames <- paste("rep", 1:nrep, sep = "")
    if (is.null(vnames)) 
        vnames <- paste("value", 1:nvar, sep = "")
    if (!is.null(df.) && df. < n) {
        if (ndim < 3) {
            gcv <- rep(0, nrep)
            for (i in 1:nrep) {
                SSEi <- sum((y[1:n, i] - yhat[, i])^2)
                gcv[i] <- (SSEi/n)/((n - df.)/n)^2
            }
            if (ndim > 1) 
                names(gcv) <- ynames
        }
        else {
            gcv <- matrix(0, nrep, nvar)
            for (ivar in 1:nvar) {
                for (i in 1:nrep) {
                  SSEi <- sum((y[1:n, i, ivar] - yhat[, i, ivar])^2)
                  gcv[i, ivar] <- (SSEi/n)/((n - df.)/n)^2
                }
            }
            dimnames(gcv) <- list(ynames, vnames)
        }
    }
    else {
        gcv <- NULL
    }
    if (is.null(fdnames)) {
        fdnames <- list(time = tnames, reps = ynames, values = vnames)
    }
    if (ndim < 3) {
        coef <- as.matrix(coef)
        fdobj <- fd(coef[1:nbasis, ], basisobj, fdnames)
    }
    else {
        fdobj <- fd(coef[1:nbasis, , ], basisobj, fdnames)
    }
    if (!is.null(penmat) && !is.null(covariates)) 
        penmat <- penmat[1:nbasis, 1:nbasis]
    smoothlist <- list(fd = fdobj, df = df., gcv = gcv, beta = beta., 
        SSE = SSE, penmat = penmat, y2cMap = y2cMap, argvals = argvals, 
        y = y0)
    class(smoothlist) <- "fdSmooth"
    return(smoothlist)
}

#'
#'
#'
#'

smooth.basis2 <- function (argvals = matrix(1:n, n, N), y, fdParobj, wtvec = NULL, 
    fdnames = NULL, covariates = NULL, method = "chol", dfscale = 1, 
    returnMatrix = FALSE) 
{
    dimy <- dim(y)
    ndy <- length(dimy)
    n <- dimy[1]
    N <- dimy[2]
    ynames <- dimnames(y)
    argNames <- dimnames(argvals)
    if (ndy < 3) {
        sb1 <- smooth.basis1(argvals[, 1], y = y[, 1], fdParobj = fdParobj, 
            wtvec = wtvec, fdnames = fdnames, covariates = covariates, 
            method = method, dfscale = dfscale, returnMatrix = returnMatrix)
        dimc1 <- dim(sb1$fd$coefs)
        dimc <- c(dimc1[1], dimy[-1])
        coefs <- array(NA, dim = dimc)
        c1names <- dimnames(sb1$fd$coefs)
        cNames <- vector("list", 2)
        if (!is.null(c1names[[1]])) 
            cNames[[1]] <- c1names[[1]]
        if (!is.null(ynames[[2]])) 
            cNames[[2]] <- ynames[[2]]
        dimnames(coefs) <- cNames
        coefs[, 1] <- sb1$fd$coefs
        if (!is.null(covariates)) {
            q <- dim(covariates)[2]
            beta. <- matrix(0, q, dimy[2])
            beta.[, 1] <- sb1$beta
        }
        else {
            beta. <- NULL
        }
        for (i in seq(2, length = dimy[2] - 1)) {
            sbi <- smooth.basis1(argvals[, i], y = y[, i], fdParobj = fdParobj, 
                wtvec = wtvec, fdnames = fdnames, covariates = covariates, 
                method = method, dfscale = dfscale, returnMatrix = returnMatrix)
            coefs[, i] <- sbi$fd$coefs
            if (!is.null(covariates)) {
                beta.[, i] <- sbi$beta
            }
        }
        if (is.null(fdnames)) {
            fdnames <- sb1$fdnames
            if (is.null(fdnames)) 
                fdnames <- list(time = NULL, reps = NULL, values = "value")
            valueChk <- ((length(fdnames$values) == 1) && (fdnames$values == 
                "value") && (length(fdnames$reps) == 1) && (!is.null(ynames[[2]])))
            if (valueChk) 
                fdnames$values <- fdnames$reps
            if (!is.null(ynames[[2]])) 
                fdnames[[2]] <- ynames[[2]]
        }
    }
    else {
        sb1 <- smooth.basis1(argvals[, 1], y = y[, 1, ], fdParobj = fdParobj, 
            wtvec = wtvec, fdnames = fdnames, covariates = covariates, 
            method = method, dfscale = dfscale, returnMatrix = returnMatrix)
        coef1 <- sb1$fd$coefs
        dimc1 <- dim(coef1)
        dimc <- c(dimc1[1], dimy[-1])
        coefs <- array(NA, dim = dimc)
        yNames <- dimnames(y)
        c1Names <- dimnames(coef1)
        cNames <- vector("list", 3)
        if (!is.null(c1Names[[1]])) 
            cNames[[1]] <- c1Names[[1]]
        if (!is.null(yNames[[2]])) 
            cNames[[2]] <- yNames[[2]]
        if (is.null(c1Names[[2]])) {
            if (!is.null(yNames[[3]])) 
                cNames[[3]] <- yNames[[3]]
        }
        else {
            cNames[[3]] <- c1Names[[2]]
        }
        dimnames(coefs) <- cNames
        coefs[, 1, ] <- coef1
        if (!is.null(covariates)) {
            q <- dim(covariates)[2]
            beta. <- array(0, c(q, dimy[2], dimy[3]))
            beta.[, , 1] <- sb1$beta
        }
        else {
            beta. <- NULL
        }
        for (i in seq(2, length = dimy[2] - 1)) {
            sbi <- smooth.basis1(argvals[, i], y = y[, i, ], 
                fdParobj = fdParobj, wtvec = wtvec, fdnames = fdnames, 
                covariates = covariates, method = method, dfscale = dfscale)
            coefs[, i, ] <- sbi$fd$coefs
            if (!is.null(covariates)) {
                beta.[, , i] <- sbi$beta
            }
            else {
                beta. <- NULL
            }
        }
        if (is.null(fdnames)) {
            fdnames <- sb1$fdnames
            if (is.null(fdnames)) {
                fdnames <- list(time = NULL, reps = NULL, values = NULL)
                if (!is.null(argNames[[1]])) {
                  fdnames[[1]] <- argNames[[1]]
                }
                else {
                  fdnames[[1]] <- ynames[[1]]
                }
                if (!is.null(ynames[[2]])) 
                  fdnames[[2]] <- ynames[[2]]
                if (!is.null(ynames[[3]])) 
                  fdnames[[3]] <- ynames[[3]]
            }
        }
    }
    sb <- sb1
    sb$beta <- beta.
    sb$fd$coefs <- coefs
    sb$fd$fdnames <- fdnames
    sb
}

#'
#'
#'
#'
smooth.basis3 <- function (argvals = array(1:n, c(n, N, M)), y, fdParobj, wtvec = NULL, 
    fdnames = NULL, covariates = NULL, method = "chol", dfscale = 1, 
    returnMatrix = FALSE) 
{
    dimy <- dim(y)
    ndy <- length(dimy)
    n <- dimy[1]
    N <- dimy[2]
    M <- dimy[3]
    if (ndy < 3) 
        stop("length(dim(y)) must be 3  is ", ndy)
    if (any(dima != dimy)) {
        stop("dim(argvals) = ", paste(dima, collapse = ", "), 
            " != dim(y) = ", paste(dimy, collapse = ", "))
    }
    dima <- dim(argvals)
    nda <- length(dima)
    if (nda < 3) 
        stop("length(dim(argvals)) must be 3  is ", nda)
    sb1 <- smooth.basis2(argvals[, , 1], y = y[, , 1], fdParobj = fdParobj, 
        wtvec = wtvec, fdnames = fdnames, covariates = covariates, 
        method = method, dfscale = dfscale, returnMatrix = returnMatrix)
    coef1 <- sb1$fd$coefs
    dimc1 <- dim(coef1)
    dimc <- c(dimc1[1], dimy[-1])
    coefs <- array(NA, dim = dimc)
    argNames <- dimnames(argvals)
    yNames <- dimnames(y)
    c1Names <- dimnames(coef1)
    cNames <- vector("list", 3)
    if (!is.null(c1Names[[1]])) 
        cNames[[1]] <- c1Names[[1]]
    if (!is.null(yNames[[2]])) 
        cNames[[2]] <- yNames[[2]]
    if (!is.null(yNames[[3]])) 
        cNames[[3]] <- yNames[[3]]
    dimnames(coefs) <- cNames
    if (!is.null(covariates)) {
        q <- dim(covariates)[2]
        beta. <- array(0, c(q, dimy[2], dimy[3]))
        beta.[, , 1] <- sb1$beta
    }
    else {
        beta. <- NULL
    }
    for (i in seq(2, length = dimy[3] - 1)) {
        sbi <- smooth.basis2(argvals[, , i], y = y[, , i], fdParobj = fdParobj, 
            wtvec = wtvec, fdnames = fdnames, covariates = covariates, 
            method = method, dfscale = dfscale, returnMatrix = returnMatrix)
        coefs[, , i] <- sbi$fd$coefs
        if (!is.null(covariates)) {
            beta.[, , i] <- sbi$beta
        }
    }
    if (is.null(fdnames)) {
        fdnames <- list(time = NULL, reps = NULL, values = NULL)
        if (!is.null(yNames[[1]])) {
            fdnames[[1]] <- yNames[[1]]
        }
        else {
            if (!is.null(argNames[[1]])) 
                fdnames[[1]] <- argNames[[1]]
        }
        if (!is.null(yNames[[2]])) {
            fdnames[[2]] <- yNames[[2]]
        }
        else {
            if (!is.null(argNames[[2]])) 
                fdnames[[2]] <- argNames[[2]]
        }
        if (!is.null(yNames[[3]])) {
            fdnames[[3]] <- yNames[[3]]
        }
        else {
            if (!is.null(argNames[[3]])) 
                fdnames[[3]] <- argNames[[3]]
        }
    }
    sb <- sb1
    sb$fd$coefs <- coefs
    sb$fd$fdnames <- fdnames
    sb$beta <- beta.
    sb
}
############################################################

ycheck <- function (y, n) 
{
    if (is.vector(y)) 
        y <- as.matrix(y)
    if (!inherits(y, "matrix") && !inherits(y, "array")) 
        stop("Y is not of class matrix or class array.")
    ydim = dim(y)
    if (ydim[1] != n) 
        stop("Y is not the same length as ARGVALS.")
    ndim = length(ydim)
    if (ndim == 2) {
        ncurve = ydim[2]
        nvar = 1
    }
    if (ndim == 3) {
        ncurve = ydim[2]
        nvar = ydim[3]
    }
    if (ndim > 3) 
        stop("Second argument must not have more than 3 dimensions")
    return(list(y = y, ncurve = ncurve, nvar = nvar, ndim = ndim))
}

###############################################################

wtcheck <- function (n, wtvec = NULL) 
{
    if (n != round(n)) 
        stop("n is not an integer.")
    if (n < 1) 
        stop("n is less than 1.")
    if (!is.null(wtvec)) {
        dimw = dim(as.matrix(wtvec))
        if (any(is.na(as.vector(wtvec)))) 
            stop("WTVEC has NA values.")
        if (all(dimw == n)) {
            onewt = FALSE
            matwt = TRUE
            wteig = eigen(wtvec)$values
            if (any(is.complex(wteig))) 
                stop("Weight matrix has complex eigenvalues.")
            if (min(wteig) <= 0) 
                stop("Weight matrix is not positive definite.")
        }
        else {
            if ((length(dimw) > 1 && dimw[1] > 1 && dimw[2] > 
                1) || length(dimw) > 2) {
                stop("WTVEC is neither a vector nor a matrix of order n.")
            }
            wtvec = as.matrix(wtvec)
            if (length(wtvec) == 1) {
                wtvec = wtvec * matrix(1, n, 1)
            }
            if (length(wtvec) != n) {
                stop("WTVEC of wrong length")
            }
            if (min(wtvec) <= 0) 
                stop("Values in WTVEC are not positive.")
            onewt = FALSE
            matwt = FALSE
        }
    }
    else {
        wtvec = matrix(1, n, 1)
        onewt = TRUE
        matwt = FALSE
    }
    return(list(wtvec = wtvec, onewt = onewt, matwt = matwt))
}




