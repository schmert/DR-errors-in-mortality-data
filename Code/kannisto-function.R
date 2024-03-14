# the Kannisto() function extends a logmx schedule to
# very old ages by fitting a regression to the logits
# of high ages

Kannisto = function(logmx, 
                    ages        = 0:89,        
                    fit_ages    = 70:89,
                    extrap_ages = 90:119) {
  
  fit_mx     = exp(logmx[ages %in% fit_ages])
  fit_logits = log(fit_mx / (1-fit_mx))
  logit_coef = coef( lm( fit_logits ~ fit_ages) )
  loga       = logit_coef[1]
  b          = logit_coef[2]
  
  extrap_logits = loga + b*extrap_ages
  extrap_logmx  = log( 1 / (1+exp(-extrap_logits)) )
  
  tmp        = c(logmx, extrap_logmx)
  names(tmp) = c(ages,extrap_ages)
  
  return(tmp)
  
} # Kannisto
