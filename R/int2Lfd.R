

int2Lfd <- function (m = 0) 
{
    if (inherits(m, "Lfd")) {
        Lfdobj <- m
        return(Lfdobj)
    }
    if (!is.numeric(m)) 
        stop("Argument not numeric and not a linear differential operator.")
    if (length(m) != 1) 
        stop("Argument is not a scalar.")
    if (round(m) != m) 
        stop("Argument is not an integer.")
    if (m < 0) 
        stop("Argument is negative.")
    if (m == 0) {
        bwtlist <- NULL
    }
    else {
        basisobj <- create.constant.basis(c(0, 1))
        bwtlist <- vector("list", m)
        for (j in 1:m) bwtlist[[j]] <- fd(0, basisobj)
    }
    Lfdobj <- Lfd(m, bwtlist)
    return(Lfdobj)
}
