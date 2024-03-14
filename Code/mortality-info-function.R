mortality_info = function(Dx,Nx) {
  
  age       = seq(Dx)-1  # 0,1,2,... 
  mx        = Dx/Nx
  names(mx) = age
  logmx     = log(mx)
  Hx        = c(0, cumsum(mx))
  lx        = exp(-Hx)
  names(lx) = age
  
  closed_Lx = (head(lx,-1) + tail(lx,-1))/2
  open_Lx   = last(lx)/last(mx)
  
  Tx = c(closed_Lx,open_Lx) %>% rev() %>% cumsum() %>% rev()
  ex = Tx/lx
  
  delta   = 1/2 * (head(Tx,-1)+tail(Tx,-1)) /Nx 
  epsilon = -mx * delta
  
  return(list(
    age=age,
    Dx=Dx, Nx=Nx, mx=mx, logmx=logmx, lx=lx,
    Tx=Tx, ex=ex, e0=ex['0'], 
    delta=delta, epsilon=epsilon
  ))
}