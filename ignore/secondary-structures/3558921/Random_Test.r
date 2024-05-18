

###################################################################################################
############## The function for Monte Carlo randomization test of modified MTE model###############
###################################################################################################
#cod: The index of K-means groups (e.g. cod= c(1,1,2,1,2,2......))
#res: The dataset(mxn)consists of the residuals of size spectra from linear fit.
#     m is the number of observation & n is the number of size class
#x: The size class (log body volume or weight) of size spectra corresponding to each column of res
#n: The times of Monte Carlo simulation
Random_Test=function(cod,res,x,n=1000){
  library(MASS)
  lcr=list(qua=NULL,cub=NULL);bpr=list(qua=NULL,cub=NULL)
  for(k in 1:n){
    lcrr=bprr=NULL
    for(i in 1:ng){
      # resampled size spectral dataset
      res.r=res[sample(which(cod!=i),tb[i]),]
      x.t=rep(x,each=nrow(res.r))

      # Generate null predictions (SP-> STL; Null predicted size-TL relationship)
      res.rn=matrix(c(res.r),ncol=1)
      X=cbind(rep(1),x.t)
      P=cbind(x.t^2,x.t^3)
      Q=-(diag(1,nrow(X))- X%*%solve(t(X)%*%X)%*%t(X))%*%P
      bprr=rbind(bprr,t(ginv(t(Q)%*%Q)%*%t(Q)%*%res.rn))

      # Generate null observations (STL-> SP; Null size-deviation function)
      los.t=lm( c(res.r)~I(x.t)+I(x.t^2)+I(x.t^3) )
      lcrr=rbind(lcrr,los.t$coefficients[3:4])
    }
    # The collection of null predictions (qua: the 2nd order, cub: the 3rd order)
    bpr$qua=cbind(bpr$qua,bprr[,1]);  bpr$cub=cbind(bpr$cub,bprr[,2])
    # The collection of null observation (qua: the 2nd order, cub: the 3rd order)
    lcr$qua=cbind(lcr$qua,lcrr[,1]);  lcr$cub=cbind(lcr$cub,lcrr[,2])

  }
  return(list(lcr,bpr))
}