


basisfd <- function (type, rangeval, nbasis, params, dropind = vector("list", 
    0), quadvals = vector("list", 0), values = vector("list", 
    0), basisvalues = vector("list", 0), basisfuns=NULL) 
{
    if (nargs() == 0) {
        type <- "bspline"
        rangeval <- c(0, 1)
        nbasis <- 2
        params <- vector("list", 0)
        dropind <- vector("list", 0)
        quadvals <- vector("list", 0)
        values <- vector("list", 0)
        basisvalues <- vector("list", 0)
        basisobj <- list(type = type, rangeval = rangeval, nbasis = nbasis, 
            params = params, dropind = dropind, quadvals = quadvals, 
            values = values, basisvalues = basisvalues)
        oldClass(basisobj) <- "basisfd"
        return(basisobj)
    }
    if (class(type) == "basisfd") {
        basisobj <- type
        return(basisobj)
    }
    if (type == "bspline" || type == "Bspline" || type == "spline" || 
        type == "Bsp" || type == "bsp") {
        type = "bspline"
    }
    else if (type == "con" || type == "const" || type == "constant") {
        type = "const"
    }
    else if (type == "exp" || type == "expon" || type == "exponential") {
        type = "expon"
    }
    else if (type == "Fourier" || type == "fourier" || type == 
        "Fou" || type == "fou") {
        type = "fourier"
    }
    else if (type == "mon" || type == "monom" || type == "monomial") {
        type = "monom"
    }
    else if (type == "polyg" || type == "polygon" || type == 
        "polygonal") {
        type = "polyg"
    }
    else if (type == "poly" || type == "pol" || type == "polynom" || 
        type == "polynomial") {
        type = "polynom"
    }
    else if (type == "pow" || type == "power") {
        type = "power"
    }
    else if( type == "empirical" || type == "empir"){
    	type = "empirical"
    }
    else {
        type = "unknown"
    }
    if (type == "unknown") {
        stop("'type' unrecognizable.")
    }
    if (missing(quadvals)) 
        quadvals <- vector("list", 0)
    else if (!(length(quadvals) == 0 || is.null(quadvals))) {
        nquad <- dim(quadvals)[1]
        ncol <- dim(quadvals)[2]
        if ((nquad == 2) && (ncol > 2)) {
            quadvals <- t(quadvals)
            nquad <- dim(quadvals)[1]
            ncol <- dim(quadvals)[2]
        }
        if (nquad < 2) 
            stop("Less than two quadrature points are supplied.")
        if (ncol != 2) 
            stop("'quadvals' does not have two columns.")
    }
    if (!(length(values) == 0 || missing(values) || is.null(values))) {
        n <- dim(values)[1]
        k <- dim(values)[2]
        if (n != nquad) 
            stop(paste("Number of rows in 'values' not equal to number of", 
                "quadrature points."))
        if (k != nbasis) 
            stop(paste("Number of columns in 'values' not equal to number of", 
                "basis functions."))
    }
    else values <- vector("list", 0)
    if (!(length(basisvalues) == 0 || missing(basisvalues) || 
        !is.null(basisvalues))) {
        if (!is.list(basisvalues)) 
            stop("BASISVALUES is not a list object.")
        sizevec <- dim(basisvalues)
        if (length(sizevec) != 2) 
            stop("BASISVALUES is not 2-dimensional.")
        for (i in 1:sizevec[1]) {
            if (length(basisvalues[[i, 1]]) != dim(basisvalues[[i, 
                2]])[1]) 
                stop(paste("Number of argument values not equal number", 
                  "of values."))
        }
    }
    else basisvalues <- vector("list", 0)
    if (missing(dropind)) 
        dropind <- vector("list", 0)
    if (type == "fourier") {
        paramvec <- rangeval[2] - rangeval[1]
        period <- params[1]
        if (period <= 0) 
            stop("Period must be positive for (a Fourier basis")
        params <- period
        if ((2 * floor(nbasis/2)) == nbasis) 
            nbasis <- nbasis + 1
    }
    else if (type == "bspline") {
        if (!missing(params)) {
            nparams <- length(params)
            if (nparams > 0) {
                if (params[1] <= rangeval[1]) 
                  stop("Smallest value in BREAKS not within RANGEVAL")
                if (params[nparams] >= rangeval[2]) 
                  stop("Largest value in BREAKS not within RANGEVAL")
            }
        }
    }
    else if (type == "expon") {
        if (length(params) != nbasis) 
            stop("No. of parameters not equal to no. of basis fns for (exponential basisobj$")
    }
    else if (type == "polyg") {
        if (length(params) != nbasis) 
            stop("No. of parameters not equal to no. of basis fns for (polygonal basisobj$")
    }
    else if (type == "power") {
        if (length(params) != nbasis) 
            stop("No. of parameters not equal to no. of basis fns for (power basisobj$")
    }
    else if (type == "const") {
        params <- 0
    }
    else if (type == "monom") {
        if (length(params) != nbasis) 
            stop("No. of parameters not equal to no. of basis fns for (monomial basisobj$")
    }
    else if (type == "polynom") {
        if (length(params) > 1) 
            stop("More than one parameter for (a polynomial basisobj$")
    }
    else if (type == "empirical"){
    	
    }
    else stop("Unrecognizable basis")
    obj.call <- match.call()
    if (type == "empirical"){
    	basisobj <- list(call = obj.call, type = type, rangeval = rangeval, 
        	nbasis = nbasis, params = params, dropind = dropind, 
        	quadvals = quadvals, values = values, basisvalues = basisvalues, basisfuns = basisfuns)
    }
    else{
    	basisobj <- list(call = obj.call, type = type, rangeval = rangeval, 
        	nbasis = nbasis, params = params, dropind = dropind, 
        	quadvals = quadvals, values = values, basisvalues = basisvalues)
    }
    oldClass(basisobj) <- "basisfd"
    basisobj
}
##################################################