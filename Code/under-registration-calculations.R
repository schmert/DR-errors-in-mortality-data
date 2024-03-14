# calculate the under-registration derivatives for
# plot

# Theoretical calculation with derivs, assuming
# P=C=Q=V=I, N=n, D=d at starting point

e0_effect = function(y,kc=1,kv=1) {
  Ty = mean( M$Tx[paste(y+0:1)] ) 
  eff = -Ty * (kc-kv) * M$mx[paste(y)]
  return(eff)
}

mx_effect = function(y,kc=1,kv=1) {
  eff = (kc-kv) * M$mx[paste(y)]
  return(eff)
}

e0_empirical_effect = function(y,kc=1, kv=1) {
  eps = .01
  P = C = Q = V = diag(length(age))
  dimnames(P) = list(age,age)
  dimnames(C) = list(age,age)
  dimnames(Q) = list(age,age)
  dimnames(V) = list(age,age)
  
  iy = paste(y)
  C[iy,iy] = C[iy,iy] - kc*eps
  V[iy,iy] = V[iy,iy] - kv*eps
  
  nn = P %*% C %*% Nx
  dd = Q %*% V %*% Dx
  
  M2 = mortality_info(dd,nn)
  
  tmp = (M2$e0 - M$e0)/eps
  
  return(tmp)
} 

sel_ages = 0:98

df = tibble( y=sel_ages) %>% 
       mutate(Ceff = map_dbl(y, e0_effect, kc=.01, kv=  0),
              Deff = map_dbl(y, e0_effect, kc=  0, kv=.01),
              Ceff_emp = map_dbl(y, e0_empirical_effect, kc=.01, kv=  0),
              Deff_emp = map_dbl(y, e0_empirical_effect, kc=  0, kv=.01),
              Ceff_mx  = map_dbl(y,mx_effect, kc=.01, kv=0),
              Deff_mx  = map_dbl(y,mx_effect, kc= 0,  kv=0.01))
