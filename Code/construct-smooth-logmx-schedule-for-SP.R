#.......................................................
# modified version of 
#   construct mortality_info.R 
# from
#   https://www.demographic-research.org/volumes/vol42/17/files/demographic-research.42-17.zip
# or  
#   https://www.overleaf.com/project/5c8fa0afe7632a543424e2fa
#   
# creates a smooth version of a discrete mortality schedule over a 
# fine age grid, and extrapolates to very high ages (300!)
#.......................................................

library(splines)
library(Matrix)

## construct a smoothed version of log mortality rates ----

deltax     = .10
max_age    = 300

logmx_raw      = log(Dx_SP/Nx_SP)  # treat as avgs over [0,1)...[99,100)

xgrid = seq(from=deltax/2, to=max_age-deltax/2, by=deltax)

n = length(xgrid) # (large) number of params

W = outer( 0:99, xgrid, function(a,b) {deltax*(floor(b)==a)})

D1 = diag(0, n-1, n)
diag(D1) = -1
diag(D1[,-1]) = +1

K = as.matrix( crossprod( Matrix(D1)))

##quad programming problem ----
## min a'Ka subj to W a=logmx_raw 

Amat = rbind( cbind( 2*K,  t(W)),
              cbind( W  , diag(0, length(logmx_raw))))
bvec = c( rep(0,n), logmx_raw)

a = solve( Amat, bvec)[1:n]

logmx_nonparam = a

# smooth the nonparametric rates with bsplines ----

B = bs(xgrid, knots=c(0,1,10,20,seq(40,140,20)))

logmx_spline = B %*% solve(crossprod(B)) %*% t(B) %*% logmx_nonparam

# add Kannisto extrapolations ages 100+ ----
# using parameters fit from the spline rates at ages (90,100)

ikeep   = which(xgrid <= 100)
iregr   = which(xgrid > 90 & xgrid < 100)   # which xgrid values for regression
ipred   = which(xgrid > 100)                # which xgrid values get extrapolated rates
ianchor = max(iregr)                        # which (x,y) value is the point from which to extrap 

# Kannisto fit and prediction
logit = function(p) {log(p/(1-p))}

ystar = logit(exp(logmx_spline[ianchor]))
xstar = xgrid[ianchor]

yy = logit( exp(logmx_spline[iregr]) ) - ystar
xx = xgrid[iregr] - xstar

b = sum(xx*yy)/sum(xx*xx)
a = exp( ystar  -b*xstar ) 

numer       = a * exp(b * xgrid[ipred])

logmx_fine =  c( logmx_spline[ikeep],
                          log( numer/(1+numer) ))  

#... extra for this paper: accumulate back to single-year
# grid for ages 0,1,...119

W2 = outer( age, xgrid, function(a,b) {deltax*(floor(b)==a)})

logmx_smooth        = W2 %*% logmx_fine
names(logmx_smooth) = 0:119
