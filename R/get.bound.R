get.bound <- function(t,c,kappa.max){
  
  ### Checks
  if(!( is.numeric(t) & is.numeric(c) & is.numeric(kappa.max) )) stop('c and kappa.max should be numeric')

  if(t<0 | t>1) stop('The threshold t should be between s1 and s2, hence between 0 and 1')
  
  if(c<0) stop('c should be non-negative')
  if(c>1) stop('c should be small')
  
  if(kappa.max<=0) stop('kappa.max should be positive')
  
  ### Compute bound for the nr of false positives, evaluated at t
  B = floor((t+c)/kappa.max)
  
  B
}
 