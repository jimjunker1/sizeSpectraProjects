########################################################################## 
##  Prediction 1: Predict size-TL relationship by ISD (ISD->STL)#
########################################################################## 
prediction1=function(x,res,m=3){                      # Prediction function
  # x: size categories in ISD (log);
  # res: The residuals of ISD from power-law fitting
  # m: the order of polynomial fitting
  library(MASS)
  X=cbind(rep(1),x)                                   # Power-law terms
  n= length(x)                                        # Sample size
  P=NULL;for(i in 2:m){P=cbind(P,x^i)}                # Higher order terms
  nQ=-(diag(1,nrow(X))-X%*%solve(t(X)%*%X)%*%t(X))%*%P# The covariates for the multiple regression
  bp=c(ginv(t(nQ)%*%nQ)%*%t(nQ)%*%res)                # Regression coefficients for the multiple regression
  mse.t=sum((res-nQ%*%bp)^2)/(n-(m+1))                # Mean square error (MSE) of the predictive fitting
  up11=bp[1]+qt(0.975,n-(m+1))*(mse.t*(t(c(1,0))%*%ginv(t(nQ)%*%nQ)%*%c(1,0)))^0.5# Compute 95%CI for the estimated coefficients
  up12=bp[2]+qt(0.975,n-(m+1))*(mse.t*(t(c(0,1))%*%ginv(t(nQ)%*%nQ)%*%c(0,1)))^0.5#
  lp11=bp[1]-qt(0.975,n-(m+1))*(mse.t*(t(c(1,0))%*%ginv(t(nQ)%*%nQ)%*%c(1,0)))^0.5#
  lp12=bp[2]-qt(0.975,n-(m+1))*(mse.t*(t(c(0,1))%*%ginv(t(nQ)%*%nQ)%*%c(0,1)))^0.5#
  cfd=rbind(sort(c(up11,lp11)),sort(c(up12,lp12)))
  colnames(cfd)=c("2.5%","97.5%");names(bp)=rownames(cfd)=c("2nd order","3rd order")
  return(list(Predicted_coefficients_P1=bp,Confidence_intervals_P1=cfd))# Return the predicted higher order coefficients of size-TL relationship with 95%CIs
}



##################################################################################################################### 
## Prediction2: Predict the deviation of ISD from power-law distribution by size-TL relationship (STL->ISD)#  
##################################################################################################################### 
prediction2=function(rang=c(4,29),sizeTL,nb=1000){                    # Prediction function
# rang: the size range of ISD in log space 
# sizeTL: The data set consist of size categories and corresponding trophic level
# nb: number of resampling in wild bootstrap
  lmtl=lm(sizeTL[,2]~I(sizeTL[,1])+I(sizeTL[,1]^2)+I(sizeTL[,1]^3))   # Empirical size-TL relationship
  cost=lmtl$coefficients[-1]                                          # The polynomial coefficients of the empirical size-TL relationship
  cofb=NULL
  for(j in 1:nb){                                                     # Wild bootstrap with nb times resampling
    str.t= lmtl$fitted+rnorm(nrow(sizeTL))*residuals(lmtl)
    wco=lm(str.t~sizeTL[,1]+I(sizeTL[,1]^2)+I(sizeTL[,1]^3))$coefficients
    cofb=rbind(cofb,simNBSS(cost=wco[-1]))
  }
  wci=t(apply(cofb,2,sort)[c(25,975),])                               # The results of wild bootstrap
  colnames(wci)=c("2.5%","97.5%");rownames(wci)=c("2nd order","3rd order")
  return(list(Predicted_coefficients=simNBSS(cost=cost),Confidence_intervals=wci))# Return the predicted coefficients of size-deviation relationship with 95%CI
}

simNBSS=function(rang=c(4,29),cost){
# rang: the size range of ISD in log space 
# cost: the polynomial coefficents of size-TL relationship (exclude intercept)
  x.c=seq(rang[1],rang[2],0.01)                       # Size categories
  m=length(cost)                                      # Polynomial order
  stlc=0;for(i in 1:m){stlc=stlc+cost[i]*x.c^(i-1)}   # Predictor from size-TL relationship
  y=30+(-(stlc)-0.75)*x.c                             # ISD simulation
  rd=residuals(lm(y~x.c))                             # Simulated deviation of ISD from linear fitting
  X=NULL;for(i in 1:m){X=cbind(X,x.c^i)}              
  sdiv=lm(rd~.,data=data.frame(X))                    # The Simulated relationship between body size and the deviation of ISDs 
  cofs=sdiv$coefficients[3:4];names(cofs)=c("2nd order","3rd order")# The polynomial coefficients of the simulated size-deviation relationship
  return(cofs)
}



