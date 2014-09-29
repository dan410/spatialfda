

pca.fd <- function (fdobj, nharm = 2, harmfdPar = fdPar(fdobj), centerfns = TRUE) 
{
    if (!(inherits(fdobj, "fd"))) 
        stop("Argument FD  not a functional data object.")
    meanfd <- mean.fd(fdobj)
    if (centerfns) 
        fdobj <- center.fd(fdobj)
    coef <- fdobj$coefs
    coefd <- dim(coef)
    ndim <- length(coefd)
    nrep <- coefd[2]
    coefnames <- dimnames(coef)
    if (nrep < 2) 
        stop("PCA not possible without replications.")
    basisobj <- fdobj$basis
    nbasis <- basisobj$nbasis
    type <- basisobj$type
    harmbasis <- harmfdPar$fd$basis
    nhbasis <- harmbasis$nbasis
    Lfdobj <- harmfdPar$Lfd
    lambda <- harmfdPar$lambda
    if (ndim == 3) {
        nvar <- coefd[3]
        ctemp <- matrix(0, nvar * nbasis, nrep)
        for (j in 1:nvar) {
            index <- 1:nbasis + (j - 1) * nbasis
            ctemp[index, ] <- coef[, , j]
        }
    }
    else {
        nvar <- 1
        ctemp <- coef
    }
    Lmat <- eval.penalty(harmbasis, 0)
    if (lambda > 0) {
        Rmat <- eval.penalty(harmbasis, Lfdobj)
        Lmat <- Lmat + lambda * Rmat
    }
    Lmat <- (Lmat + t(Lmat))/2
    Mmat <- chol(Lmat)
    Mmatinv <- solve(Mmat)
    Wmat <- crossprod(t(ctemp))/nrep
    Jmat = inprod(harmbasis, basisobj)
    MIJW = crossprod(Mmatinv, Jmat)
    if (nvar == 1) {
        Cmat = MIJW %*% Wmat %*% t(MIJW)
    }
    else {
        Cmat = matrix(0, nvar * nhbasis, nvar * nhbasis)
        for (i in 1:nvar) {
            indexi <- 1:nbasis + (i - 1) * nbasis
            for (j in 1:nvar) {
                indexj <- 1:nbasis + (j - 1) * nbasis
                Cmat[indexi, indexj] <- MIJW %*% Wmat[indexi, 
                  indexj] %*% t(MIJW)
            }
        }
    }
    Cmat <- (Cmat + t(Cmat))/2
    result <- eigen(Cmat)
    eigvalc <- result$values
    eigvecc <- as.matrix(result$vectors[, 1:nharm])
    sumvecc <- apply(eigvecc, 2, sum)
    eigvecc[, sumvecc < 0] <- -eigvecc[, sumvecc < 0]
    varprop <- eigvalc[1:nharm]/sum(eigvalc)
    if (nvar == 1) {
        harmcoef <- Mmatinv %*% eigvecc
    }
    else {
        harmcoef <- array(0, c(nbasis, nharm, nvar))
        for (j in 1:nvar) {
            index <- 1:nbasis + (j - 1) * nbasis
            temp <- eigvecc[index, ]
            harmcoef[, , j] <- Mmatinv %*% temp
        }
    }
    harmnames <- rep("", nharm)
    for (i in 1:nharm) harmnames[i] <- paste("PC", i, sep = "")
    if (length(coefd) == 2) 
        harmnames <- list(coefnames[[1]], harmnames, "values")
    if (length(coefd) == 3) 
        harmnames <- list(coefnames[[1]], harmnames, coefnames[[3]])
    harmfd <- fd(harmcoef, harmbasis, harmnames)
    if (nvar == 1) {
        harmscr <- inprod(fdobj, harmfd)
    }
    else {
        harmscr <- array(0, c(nrep, nharm, nvar))
        coefarray <- fdobj$coefs
        harmcoefarray <- harmfd$coefs
        for (j in 1:nvar) {
            fdobjj <- fd(as.matrix(coefarray[, , j]), basisobj)
            harmfdj <- fd(as.matrix(harmcoefarray[, , j]), basisobj)
            harmscr[, , j] <- inprod(fdobjj, harmfdj)
        }
    }
    pcafd <- list(harmfd, eigvalc, harmscr, varprop, meanfd)
    class(pcafd) <- "pca.fd"
    names(pcafd) <- c("harmonics", "values", "scores", "varprop", 
        "meanfd")
    return(pcafd)
}