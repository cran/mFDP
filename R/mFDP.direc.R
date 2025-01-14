mFDP.direc <- function(Ts, delta=0, gamma=0.05){
  
  if(gamma<0 | gamma>=1) stop('gamma should be at least 0 and less than 1')
  
  
  Ts_above_delta = Ts[which(Ts>delta)]
  Ts_below_minusdelta = Ts[ which(Ts < -delta) ]
  JumpsR = Ts_above_delta - delta
  JumpsRmin = -delta - Ts_below_minusdelta 
  Jumps = c(JumpsR,JumpsRmin) #all t>=0 where FDP_tilde has a jump
  Jumps_so = sort(Jumps)     #sorted jumping points
  nJumps = length(Jumps_so)
  
  FDPtilde_Jumps = numeric(nJumps)
  
  for(i in 1:nJumps){
    t= Jumps_so[i]
    R_t = sum(Ts_above_delta > delta+t)
    bound_t = sum(Ts_below_minusdelta < -delta-t)

    FDPtilde_Jumps[i] = bound_t / max(R_t,1)  #compute FDPtilde at its jumping points
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
    
    rejections = which( Ts > delta + s_plus)
  } else{
    rejections = numeric(0)
  }
  
  rejections
}
   