

fdlabels <- function (fdnames, nrep, nvar) 
{
    if (!inherits(fdnames, "list")) 
        stop("Argument fdnames is not a list object.")
    if (length(fdnames) != 3) 
        stop("Argument fdnames is not of length 3.")
    xlabel = fdnames[[1]]
    if (length(xlabel) > 1 || is.null(xlabel)) 
        xlabel = names(fdnames)[1]
    if (!is.character(xlabel)) 
        xlabel = ""
    ylabel = fdnames[[3]]
    if ((nvar > 1 && length(ylabel) == nvar) || is.null(ylabel)) 
        ylabel = names(fdnames)[3]
    if (length(ylabel) > 1) {
        if (inherits(ylabel, "character")) 
            ylabel = ylabel[1]
        else {
            if (inherits(ylabel, "list")) 
                ylabel = ylabel[[1]]
            else ylabel = ""
        }
    }
    if (!is.character(ylabel)) 
        ylabel = ""
    if (length(fdnames[[2]]) == nrep) {
        casenames = as.character(fdnames[[2]])
    }
    else {
        casenames = NULL
    }
    if (length(fdnames[[3]]) == nvar) {
        varnames = as.character(fdnames[[3]])
    }
    else {
        varnames = NULL
    }
    return(list(xlabel = xlabel, ylabel = ylabel, casenames = casenames, 
        varnames = varnames))
}
