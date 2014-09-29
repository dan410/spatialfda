

Lfd <- function (nderiv = 0, bwtlist = vector("list", 0)) 
{
    if (!is.numeric(nderiv)) 
        stop("Order of operator is not numeric.")
    if (nderiv != round(nderiv)) 
        stop("Order of operator is not an integer.")
    if (nderiv < 0) 
        stop("Order of operator is negative.")
    if (!inherits(bwtlist, "list") && !inherits(bwtlist, "fd") && 
        !is.null(bwtlist) && !missing(bwtlist)) 
        stop("BWTLIST is neither a LIST or a FD object")
    if (is.null(bwtlist)) {
        bwtlist <- vector("list", nderiv)
        if (nderiv > 0) {
            conbasis <- create.constant.basis()
            for (j in 1:nderiv) bwtlist[[j]] <- fd(0, conbasis)
        }
    }
    if (inherits(bwtlist, "fd")) 
        bwtlist <- fd2list(bwtlist)
    nbwt <- length(bwtlist)
    if (nbwt != nderiv & nbwt != nderiv + 1) 
        stop("The size of bwtlist inconsistent with NDERIV.")
    if (nderiv > 0) {
        rangevec <- c(0, 1)
        for (j in 1:nbwt) {
            bfdj <- bwtlist[[j]]
            if (inherits(bfdj, "fdPar")) {
                bfdj <- bfdj$fd
                bwtlist[[j]] <- bfdj
            }
            if (!inherits(bfdj, "fd") && !inherits(bfdj, "integer")) 
                stop(paste("An element of BWTLIST contains something other ", 
                  " than an fd object or an integer"))
            if (inherits(bfdj, "fd")) {
                bbasis <- bfdj$basis
                rangevec <- bbasis$rangeval
            }
            else {
                if (length(bfdj) == 1) {
                  bwtfd <- fd(bfdj, conbasis)
                  bwtlist[[j]] <- bwtfd
                }
                else stop("An element of BWTLIST contains a more than one integer.")
            }
        }
        for (j in 1:nbwt) {
            bfdj <- bwtlist[[j]]
            if (inherits(bfdj, "fdPar")) 
                bfdj <- bfdj$fd
            bbasis <- bfdj$basis
            btype <- bbasis$type
            if (!btype == "const") {
                brange = bbasis$rangeval
                if (any(rangevec != brange)) 
                  stop("Ranges are not compatible.")
            }
        }
    }
    Lfd.call <- match.call()
    Lfdobj <- list(call = Lfd.call, nderiv = nderiv, bwtlist = bwtlist)
    oldClass(Lfdobj) <- "Lfd"
    Lfdobj
}