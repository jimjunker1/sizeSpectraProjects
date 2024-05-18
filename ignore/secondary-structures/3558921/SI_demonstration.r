
#######Loading data #################################################################################################################################
source(file.choose())                                       # Loading prediction function: SIp1.r 
source(file.choose())                                       # Loading prediction function: SIp2.r 
sizeTL=read.table(file.choose(),sep=",",header=T)           # Loading size-TL data: sizeTL.csv
ssize=read.table(file.choose(),sep=",",header=T)            # Loading log individual size distribution(ISD) data: ssize.csv


#######Empirical size-TL relationship################################################################################################################
plot(sizeTL[,2]~sizeTL[,1],xlab="log Body size",
ylab="Trophic level",main="Size-TL relationship")           # Plot empirical size_TL relationship
lmtl=lm(sizeTL[,2]~I(sizeTL[,1])+I(sizeTL[,1]^2)+I(sizeTL[,1]^3)) # Fit size-TL relationship by polynomial function (cubbic)
tlc=lmtl$coefficients                                       # The estimates of the polynomial coefficients in empirical sizeTL relationship
tlci=confint(lmtl)                                          # The confidence interval of the polynomial coefficients 
x.c=seq(0,30,0.1)
lines(c(tlc[1]+tlc[2]*x.c+tlc[3]*x.c^2+tlc[4]*x.c^3)~x.c)   # The fitted size-TL relationship


#######Empirical deviation of ISD from power-law fitting #####################################################################################
windows()
plot(ssize[,2]~ssize[,1],xlab="log Body size",
ylab="Abundance",main="ISD")                      # Plot empirical log ISD
ress=residuals(lm(ssize[,2]~ssize[,1]))                     # Extract the deviation of log ISD from the power-law fitting 
windows()
plot(ress~ssize[,1],xlab="log Body size",,ylab="Residual",
main="Secondary structure of ISD" )       # Plot relationship between size & deviation 
silm=lm(ress~I(ssize[,1])+I(ssize[,1]^2)+I(ssize[,1]^3))    # Empirical secondary strucutre: Fitting the residuals of log ISD from power-law fit by a polynomial function of body size 
silc=silm$coefficients                                      # The estimates of the polynomial coefficients in the secondary structure of ISD
silci=confint(silm)                                         # The confidence interval of the polynomial coefficients
lines(c(silc[1]+silc[2]*ssize[,1]+silc[3]*ssize[,1]^2+silc[4]*ssize[,1]^3)~ssize[,1])# the fitted polynomial secondary structure of log ISD
names(silc)[3:4]=names(tlc)[3:4]=rownames(silci)[3:4]=rownames(tlci)[3:4]=c("2nd order","3rd order");

#######The comparison between empirical observation and model prediction#############################################################################
  ####### Prediction 1: ISD--> size-TL relationship (ISD->STL)########################################################################################
(p1=prediction1(x=ssize[,1],res=ress))                      # Execute the prediction1 function (output the predictions with 95%C.I.)
(empirical1=list(Empirical_coefficients_P1=tlc[3:4],Confidence_interval_P1=tlci[3:4,]))  # Empirical coefficients in size-TL relationship 
  ########Plot the comparison between observation (solid curcle)and prediction(solid triangle)#######################################################
windows()
par(mfcol=c(1,2))
plot(c(tlc[3])~c(1),ylim=c(0.07,0.17),xlim=c(-1,5),pch=16,ylab="",
xlab="",main="2nd order" ,axes=F,cex=2)
points(c(p1$Predicted_coefficients_P1[1])~c(3),pch=17,cex=2)
box();axis(2)
arrows(x0=1,x1=1,y0=tlci[3,1],y1=tlci[3,2],length=0.05,angle=90)
arrows(x0=1,x1=1,y0=tlci[3,2],y1=tlci[3,1],length=0.05,angle=90)
arrows(x0=3,x1=3,y0=p1$Confidence_intervals_P1[1,1],y1=p1$Confidence_intervals_P1[1,2],length=0.05,angle=90)
arrows(x0=3,x1=3,y0=p1$Confidence_intervals_P1[1,2],y1=p1$Confidence_intervals_P1[1,1],length=0.05,angle=90)
legend(0.25,0.09,c("Observation","Prediction"),pch=c(16:17))
plot(c(tlc[4])~c(1),ylim=c(-0.0028,-0.0015),xlim=c(-1,5),pch=16,ylab="",
xlab="",main="3rd order" ,axes=F,cex=2)
points(c(p1$Predicted_coefficients_P1[2])~c(3),pch=17,cex=2)
box();axis(2)
arrows(x0=1,x1=1,y0=tlci[4,1],y1=tlci[4,2],length=0.05,angle=90)
arrows(x0=1,x1=1,y0=tlci[4,2],y1=tlci[4,1],length=0.05,angle=90)
arrows(x0=3,x1=3,y0=p1$Confidence_intervals_P1[2,1],y1=p1$Confidence_intervals_P1[2,2],length=0.05,angle=90)
arrows(x0=3,x1=3,y0=p1$Confidence_intervals_P1[2,2],y1=p1$Confidence_intervals_P1[2,1],length=0.05,angle=90)

  ####### Prediction 2: size-TL relationship--> ISD (STL->ISD)########################################################################################
(p2= prediction2(sizeTL=sizeTL))                             # Execute the prediction2 function (output the predictions with 95%C.I.)
(empirical2=list(Empirical_coefficients_P2=silc[3:4],Confidence_interval_P2=silci[3:4,]))# Empirical coefficients in polynomial secondary structure of ISD (from ISD) 
  ########Plot the comparison between observation (solid curcle)and prediction(solid triangle)#######################################################
windows()
par(mfcol=c(1,2))
plot(c(silc[3])~c(1),ylim=c(-0.20,-0.08),xlim=c(-1,5),pch=16,ylab="",
xlab="",main="2nd order" ,axes=F,cex=2)
points(c(p2$Predicted_coefficients[1])~c(3),pch=17,cex=2)
box();axis(2)
arrows(x0=1,x1=1,y0=silci[3,1],y1=silci[3,2],length=0.05,angle=90)
arrows(x0=1,x1=1,y0=silci[3,2],y1=silci[3,1],length=0.05,angle=90)
arrows(x0=3,x1=3,y0=p2$Confidence_intervals[1,1],y1=p2$Confidence_intervals[1,2],length=0.05,angle=90)
arrows(x0=3,x1=3,y0=p2$Confidence_intervals[1,2],y1=p2$Confidence_intervals[1,1],length=0.05,angle=90)
legend(0.25,-0.18,c("Observation","Prediction"),pch=c(16:17))
plot(c(silc[4])~c(1),ylim=c(0.0014,0.0029),xlim=c(-1,5),pch=16,ylab="",
xlab="",main="3rd order" ,axes=F,cex=2)
points(c(p2$Predicted_coefficients[2])~c(3),pch=17,cex=2)
box();axis(2)
arrows(x0=1,x1=1,y0=silci[4,1],y1=silci[4,2],length=0.05,angle=90)
arrows(x0=1,x1=1,y0=silci[4,2],y1=silci[4,1],length=0.05,angle=90)
arrows(x0=3,x1=3,y0=p2$Confidence_intervals[2,1],y1=p2$Confidence_intervals[2,2],length=0.05,angle=90)
arrows(x0=3,x1=3,y0=p2$Confidence_intervals[2,2],y1=p2$Confidence_intervals[2,1],length=0.05,angle=90)




