mFDP.equiv <- function(Ts, delta, gamma=0.05){
  
  #checks
  if(gamma<0 | gamma>=1) stop('gamma should be at least 0 and less than 1')
  if(missing(delta)) stop('specify delta')
  if(delta<=0) stop('delta should be positive')
  
  Ts = abs(Ts) #only the absolute values are required
  
  Ts_small = Ts[which(Ts<delta)] 
  Ts_large = Ts[ which(Ts >= delta) ]
  JumpsR = delta - Ts_small 
  JumpsRmin = Ts_large - delta
  Jumps = c(JumpsR,JumpsRmin) #all t>=0 where FDP_tilde has a jump
  Jumps_so = sort(Jumps)     #sorted jumping points
  nJumps = length(Jumps_so)
  
  FDPtilde_Jumps = numeric(nJumps)
  
  for(i in 1:nJumps){
    t= Jumps_so[i]
    R_t = sum(Ts_small < delta-t)
    bound_t = sum(Ts_large > delta+t)
    if(R_t==0){  #Can only happen for equivalence testing
      FDPtilde_Jumps[i] = 0  
    } else{
      FDPtilde_Jumps[i] = bound_t / max(R_t,1)  #compute FDPtilde at its jumping points
    }
  }
  
  if(min(FDPtilde_Jumps)>gamma) rejectsome = FALSE   else   rejectsome = TRUE
  
  if(rejectsome){
    
    if(max(FDPtilde_Jumps) > gamma ){
      index_s = max(which(FDPtilde_Jumps>gamma)) 
      s_plus = Jumps_so[ min(index_s, nJumps-1) +1 ]
    }
    
    if(max(FDPtilde_Jumps) <= gamma ){
      s_plus=0
    }
    
    rejections = which(Ts < delta - s_plus)
  } else{
    rejections = numeric(0)
  }
  
  rejections #indices of the rejected hypotheses
}



