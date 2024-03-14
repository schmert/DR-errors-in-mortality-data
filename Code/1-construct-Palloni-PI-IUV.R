# construct the Palloni matrix from scratch
# with a new linear spline on the PROB scale
# that has a zero value at age 40
# and matches the value alogit(a + b*60) at age 60

library(tidyverse)
library(here)

rm(list=ls()) # clear memory

# calculate thetas for each age ----
age       = 0:99
true_age  = 0:110
rept_age  = 0:120
old_age   = 55:110

# transition ages from L to H: cubic spline f(y) such that
# f(L)  = f'(L) = 0
# f(H)  = prob. value from Palloni et al. equation
# f'(H) = prob. slope from Palloni et al. equation

L = 40
H = min(old_age)

A = rbind( c(  1,  L,  L^2,   L^3),
           c(  0,  1,  2*L, 3*L^2),
           c(  1,  H,  H^2,   H^3),
           c(  0,  1,  2*H, 3*H^2)
          )

logit  = function(p) log(p/(1-p))
alogit = function(z) 1/(1+exp(-z))

# understatement probs by age
theta_u        = 0 * true_age
names(theta_u) = true_age

theta_u[paste(old_age)] = alogit(-1.846 + 0.002 * old_age)

PH   = theta_u[paste(H)]
coef = solve(A, c(0, 0, PH, .002*PH*(1-PH) ) )

theta_u[ paste(L:H) ] = coef[1] + coef[2]*(L:H) + 
                        coef[3]*(L:H)^2 + coef[4]*(L:H)^3

# overstatement
theta_o        = 0 * true_age
names(theta_o) = true_age

theta_o[paste(old_age)] = alogit(-2.127 + 0.014 * old_age)

PH   = theta_o[paste(H)]
coef = solve(A, c(0, 0, PH, .014*PH*(1-PH)) )

theta_o[ paste(L:H)] = coef[1] + coef[2]*(L:H) + 
                       coef[3]*(L:H)^2 + coef[4]*(L:H)^3


# at ages 40+, calc off-diag elements of Pi ----
# prob of incorrect age reporting, from 
# Table 2, p. 408 of Palloni et al.
# 
# I've corrected so that they sum to one, but
# kept proportions

rho_u = c(.510, .128, .091, .052, .041, .035, .028, .026, .013, .060) %>% 
  prop.table()

rho_o = c(.621, .191, .079, .040, .023, .015, .009, .007, .005, .009) %>% 
  prop.table()

PI = matrix(0, length(rept_age), length(true_age),
            dimnames=list(rept_age,true_age))

# diagonal
for (y in true_age) {
  i = paste(y)
  PI[i,i] = 1 - theta_u[i] - theta_o[i]
}

# understatement probs
for (y in 40:max(true_age)) {
  j = paste(y)
  i = paste(y - (1:10))
  PI[i,j] = theta_u[j] * rho_u
}

# overstatement probs
for (y in 40:max(true_age)) {
  j = paste(y)
  i = paste(y + (1:10))
  PI[i,j] = theta_o[j] * rho_o
}

# collapse on 110+ for rows (won't have any
# important empirical consequences)

PI = rbind(PI[paste(0:109),],
           colSums(PI[paste(110:120),]))

rownames(PI) = 0:110

# separate matrix PI into I + U + V components

I = diag(nrow(PI))

U               = 0*PI
U[upper.tri(U)] = PI[upper.tri(PI)]
diag(U)         = -colSums(U) 

V               = 0*PI
V[lower.tri(V)] = PI[lower.tri(PI)]
diag(V)         = -colSums(V) 

# write the PI, I, U, and V matrices to a file ----

save( PI, I, U, V, file=here('Data','Palloni-PI-IUV.Rdata'))








