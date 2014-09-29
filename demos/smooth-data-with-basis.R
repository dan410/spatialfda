###################################################################
# project a spatial functional data set onto an empirical basis
###################################################################

sourceDir("R") # spatialfda
library(plyr)
library(sfdasim)
library(sseigfun)
data("sfdat") # example data set in the sseigfun package
simdat <- sfdat

nfuns <- length(unique(simdat$ID))
y <- matrix(simdat$X, ncol = nfuns, byrow = FALSE)
argvals <- simdat[simdat$ID == 1, ]$Time

# use the plyr package to extract the locations
locs <- ddply(simdat, .(ID), function(x){x[1, c("locs.Var1", "locs.Var2")]})
locs <- locs[,-1]

# create a list of functions to use as a basis
funlist <- list(f1 = function(x) cos(pi*x),
				f2 = function(x) (-1/4)*cos(2*pi*x),
				f3 = function(x) (1/9)*cos(3*pi*x))
				
funlist <- list(f1 = function(x) cos(pi*x),
				f2 = function(x) cos(2*pi*x),
				f3 = function(x) cos(3*pi*x))
				
# create a list from estimated eigenfunctions
cfit <- estimate_cov_function(simdat, nbasis=3)
eigen.fit <- estimate_eigenfunctions(cfit)
funlist <- list(f1 = function(x) eval_ef(x, eigen.fit,1),
				f2 = function(x) eval_ef(x, eigen.fit,2),
				f3 = function(x) eval_ef(x, eigen.fit,3))

# create a basis object 
emp.basis <- create.empirical.basis(basisfuns = funlist, rangeval=c(0,1))
plot(emp.basis)

# or use fixed basis
b.basis <- create.bspline.basis(rangeval=c(0,1), nbasis=6)
plot(b.basis)

fit.emp <- smooth.basis(argvals = argvals, y = y, fdParobj = emp.basis)$fd
fit.b <- smooth.basis(argvals = argvals, y = y, fdParobj = b.basis)$fd
op <- par(mfrow=c(1,2))
plot(fit.emp, ylim=c(-1,1))
plot(fit.b, ylim = c(-1,1))
par(op)

coef.fields.emp <- t(fit.emp$coefs)
coef.fields.b <- t(fit.b$coefs)

# empirical basis
dat.geo.emp <- list()
for(i in 1:emp.basis$nbasis){
	dat.geo.emp[[i]] <- as.geodata(cbind(locs, coef.fields.emp), coords.col=1:2, data.col=(2 + i))
}
variograms.emp <- llply(dat.geo.emp, variog, estimator.type="modulus", max.dist = 0.6, messages=FALSE)
variograms.fit.emp <- llply(variograms.emp, variofit, ini.cov.pars = c(.04,0.2), cov.model="exponential", fix.nugget=FALSE, nugget = 0.0)

# bspline basis
dat.geo.b <- list()
for(i in 1:b.basis$nbasis){
	dat.geo.b[[i]] <- as.geodata(cbind(locs, coef.fields.b), coords.col=1:2, data.col=(2 + i))
}
variograms.b <- llply(dat.geo.b, variog, estimator.type="modulus", max.dist = 0.6, messages=FALSE)
variograms.fit.b <- llply(variograms.b, variofit, ini.cov.pars = c(.04,0.2), cov.model="exponential", fix.nugget=FALSE, nugget = 0.0)



ids <- sample(1:nlocs, 2)
ids %in% samp.id
pred.locs <- mycurves$LOCS[ids ,]
preds.emp <- list()
for( i in 1:emp.basis$nbasis){
	preds.emp[[i]] <- krige.conv(dat.geo.emp[[i]], locations = pred.locs, krige = krige.control(obj.model = variograms.fit.emp[[i]]))
}

preds.b <- list()
for( i in 1:b.basis$nbasis){
	preds.b[[i]] <- krige.conv(dat.geo.b[[i]], locations = pred.locs, krige = krige.control(obj.model = variograms.fit.b[[i]]))
}

coef.emp <- ldply(preds.emp, function(x){x$predict})
coef.b <- ldply(preds.b, function(x){x$predict})


op <- par(mfrow=c(1,2))
emp.obj <- fd(coef=coef.emp, emp.basis)
b.obj <- fd(coef=coef.b, b.basis)
plot_curves(coef=mycurves$COEF[ids,], basis.fns=mycurves$basis.fns, ylim=c(-1,1)); abline(h=0, lty=2)
plot(emp.obj, ylim=c(-1,1), add=TRUE, col="blue")
plot_curves(coef=mycurves$COEF[ids,], basis.fns=mycurves$basis.fns, ylim=c(-1,1)); abline(h=0, lty=2)
plot(b.obj, ylim=c(-1,1), add=TRUE, col="blue")
par(op)

plot(locs.obs)
points(mycurves$LOCS[ids,], col="red", pch=16)


junk <- fd(coef=t(j), basisobj=emp.basis)
plot(junk, ylim=c(-1,1))















