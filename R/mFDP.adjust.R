mFDP.adjust <- function(P, c="1/(2m)", s1=0, s2=0.1){
  
  ### Checks
  if(s1<0) stop('s1 should be non-negative')
  if(s2>1) stop('s2 should be at most 1')
  if(s2<=s1) stop('s2 should be larger than s1')
  
  if( !((is.numeric(c))|(c=="1/(2m)")|(c=="1/m"))) stop('c should be numeric or "1/(2m)" or "1/m"')
  
  if(!(length(P)>1))  stop('P should be a vector')
  if(max(P)>1 | min(P)<0) stop('The p-values in P should be between 0 and 1')
  
  m=length(P)
  if(c=="1/(2m)") c=1/(2*m)
  if(c=="1/m") c=1/m
  
  if(c<0) stop('c should be non-negative')
  if(c>1) stop('c should be small')
 
  
  ###Compute confidence envelope (i.e., kappa.max)
  P.unsorted=P
  P=sort(P)
  
  mincutoff=s1; maxcutoff=s2
  nr.large = sum(P>=1-maxcutoff)  #nr of p-values larger than 1-maxcutoff
  nr.toolarge = sum(P>1-mincutoff)

  if(sum(P>=1-mincutoff)==0){
    kappa0=Inf
  }else{
    kappa0=(mincutoff+c)/( sum(P>=1-mincutoff))
  }
  
  #for every relevant i, compute kappa_i:
  if(nr.large>nr.toolarge){ # i.e. if there is at least 1 p-value for which 1-p is in [mincutoff, maxcutoff]
    kappas.i = numeric(m)+Inf
    for(i in (m-nr.large+1):(m-nr.toolarge)) {   #only use the p-values p for which 1-p is in [mincutoff, maxcutoff]
      kappas.i[i] = (1-P[i]+c)/( sum(P>=P[i]))   
    }
    kappa.max = min(kappa0,min(kappas.i))   #(paradoxically, kappa.max equals this minimum)
  }
  if(nr.large==nr.toolarge) kappa.max=kappa0
  if(nr.large==0) kappa.max = Inf


  
  ### Compute adjusted p-values
  Padj <- numeric(m)
  r <- sum(P<=s2)
  if(r<m) Padj[(r+1):m] <- Inf 
  if(r>0){
    if(s1<=P[r]){
      Padj[r] <- floor((P[r]+c)/(kappa.max)) / r
    } else{
      Padj[r] <- floor((s1+c)/(kappa.max)) / sum(P<=s1)
    }
    
    l <- r-1
    continue <- FALSE
    if(l>0) if(s1<=P[l]) continue <- TRUE
    while(continue==TRUE){
      Padj[l] <- min( Padj[l+1], floor((P[l]+c)/(kappa.max))/l )
      l <- l-1
      if(l==0){
        continue <- FALSE
      } else{
        if(P[l]<s1)  continue <- FALSE
      }
    }
    if(0<l) Padj[1:l] <- min( Padj[l+1], floor((s1+c)/(kappa.max))/(sum(P<=s1)))
  }
  
  
  Padj.unsorted <- numeric(m)
  Padj.unsorted[order(P.unsorted)] <- Padj
  
  Padj.unsorted
}
 
   