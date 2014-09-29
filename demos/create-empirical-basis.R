

# create a list of functions to use as a basis
funlist <- list(f1 = function(x) sin(2*pi*x),
				f2 = function(x) x,
				f3 = function(x) x^2)

# create a basis object 
mybasis <- create.empirical.basis(basisfuns = funlist, rangeval=c(0,1))
plot(mybasis)

# Currently only regression smoothing is implemented, which does not use a smoothing penalty. 
fit.emp <- smooth.basis(argvals = argvals, y = y, fdParobj = emp.basis)$fd

# your basis

chem1 <- function(x){
	eval.fd(x, fit) # fit is your fitted objected 
}











































