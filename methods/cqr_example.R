###Conformalized quantile regression


###Function to split the data into training and calibration sets

library(ggplot2)
library(ggpubr)
library(cowplot)
library(quantreg)
library(matrixcalc)
library(hdm)

########################################
###Set up random data
########################################
set.seed(1)
n<-1000
X<-rnorm(n,1,3)
U<-rnorm(n,0,1)
Y<-2+3*X+U
datXY<-cbind.data.frame(X,Y)
###data point to predict

library(matrixcalc)
sing<-TRUE

###check for singularity of training data 
##########################################
while (sing==TRUE){
  ind <- sample(n,n,replace=FALSE)
  y   <- Y[ind]
  x   <- X[ind]
  x   <- as.matrix(x)
  
  obj.uneven <- aux.uneven(y,x)
  
  Xa <- cbind(obj.uneven$X0)
  Ya<- obj.uneven$Y0
  X1a <- cbind(obj.uneven$X1)
  Y1a <- obj.uneven$Y1
  
  sing <- is.singular.matrix(t(Xa)%*%Xa)
  
}





#######################################################################################
##########################################################################################
############################################################################################
####################CONFORMALIZED QUANTILE REGRESSION#######################################
############################################################################################
#############################################################################################
alpha.sig<-0.1

Xnew<-0.67
Ynew<-2.7

point<-c(Xnew,Ynew)
newDat<-rbind.data.frame(datXY,point)
names(newDat)<-c("X","Y")

train<-cbind.data.frame(Xa,Ya)
calib<-cbind.data.frame(X1a,Y1a)


a<-ggplot(data=newDat,aes(x=X,y=Y))+geom_point()+
  geom_vline(xintercept=0.67,col="red",linetype="dashed",linewidth=1)+
  geom_hline(yintercept=2.7,col="blue",linetype="dashed",linewidth=1)

a1<-a+annotate(x=0.67,y=+Inf,label="Xtest:0.67",vjust=2,geom="label")
a2<-a1+annotate(x=+Inf,y=2.7,label="Ytest:2.7",hjust=1.25,geom="label")
a3<-a2+annotate(x=+Inf,y=+Inf,label="Entire Data",vjust=5,hjust=1,
                geom="label")
  

beta.lo <- rq.fit.br(cbind(1,Xa),Ya,tau=(alpha.sig/2))$coefficients
beta.hi <- rq.fit.br(cbind(1,Xa),Ya,tau=(1-alpha.sig/2))$coefficients
beta.50 <- rq.fit.br(cbind(1,Xa),Ya,tau=0.5)$coefficients

#############################################################################
###############################################################################

####Plot of quantiles on training data
cq<-ggplot2::ggplot(train,aes(x=Xa,y=Ya))+geom_point(col="lightgreen")+
  xlab("X(training data)")+ylab("Y(training data)")+
  geom_abline(intercept=beta.lo[1],slope=beta.lo[2],colour="orange",lwd=1.5,linetype="dashed")+
  geom_abline(intercept=beta.hi[1],slope=beta.hi[2],colour="darkred",lwd=1.5,linetype="dashed")+
  geom_abline(intercept=beta.50[1],slope=beta.50[2],colour="purple",lwd=1.5,linetype="dashed")+
  annotate(x=0,y=5,label="beta.hi",geom="label",hjust=1,fill="pink")+
  annotate(x=0,y=0,label="beta.lo",geom="label",hjust=-0.15,fill="orange")+
  annotate(x=5,y=17,label="beta.50",geom="label",fill="purple",hjust=-0.15)
####################################################################################  
  
tq.lo  <- cbind(1,X1a)%*%beta.lo
tq.hi  <- cbind(1,X1a)%*%beta.hi
tq.50  <- cbind(1,X1a)%*%beta.50
  
qsr <- t(apply(cbind(tq.lo,tq.50,tq.hi),1,FUN=sort))
  
q.lo <- qsr[,1]
q.50 <- qsr[,2]
q.hi <- qsr[,3]  

#######################################################################################

a0<- 0.9
a1<-1.1
highlight<- data.frame(start=a0, end=a1, group=seq_along(a0))

b0<-1.9
b1<-2.1
highlight2<- data.frame(start=b0, end=b1, group=seq_along(b0))

c0<-2.9
c1<-3.1
highlight3<- data.frame(start=c0, end=c1, group=seq_along(c0))

###plot of different observations for 3 different quantiles
qsrt<-data.frame(t(qsr))
qsrt$no<-1:3
cq2<-ggplot2::ggplot(qsrt,aes(x=no,y=X1))+geom_line(col="#CCFF00",size=2)+
  geom_point(col="orange",size=3)+
  geom_line(aes(x=no,y=X50),col="#FF6699",size=2)+geom_point(aes(x=no,y=X50),col="#993300",size=3)+
  geom_line(aes(x=no,y=X100),col="orange",size=2)+geom_point(aes(x=no,y=X100),col="#FF0033",size=3)+
     geom_point(aes(x=no,y=X150),col="#3300FF",size=3)+
  geom_point(aes(x=no,y=X200),col="#FF6633",size=3)+
geom_line(aes(x=no,y=X250),col="yellow",size=2)+geom_point(aes(x=no,y=X250),col="#339900",size=3)+
 geom_line(aes(x=no,y=X300),col="lightgreen",size=2)+geom_point(aes(x=no,y=X300),col="#336600",size=3)+
 geom_line(aes(x=no,y=X350),col="lightblue",size=2)+geom_point(aes(x=no,y=X350),col="#3399FF",size=3)+
  geom_point(aes(x=no,y=X400),col="#CC0066",size=3)+
  geom_line(aes(x=no,y=X450),col="darkgrey",size=2)+geom_point(aes(x=no,y=X450),col="#003300",size=3)+
  geom_line(aes(x=no,y=X500),col="#CC00CC",size=2)+geom_point(aes(x=no,y=X500),col="#FF3333",size=3)+
  xlab("Number of different quantiles")+
  ylab("Sorted coefficients")+
  geom_rect(data=highlight,inherit.aes=FALSE, 
            aes(xmin=a0, xmax=a1, ymin=-20,
                ymax=30, group=group), color="transparent", fill="magenta", alpha=0.1)+
  geom_rect(data=highlight2,inherit.aes=FALSE, 
            aes(xmin=b0, xmax=b1, ymin=-20,
                ymax=30, group=group), color="transparent", fill="purple", alpha=0.1)+
  geom_rect(data=highlight3,inherit.aes=FALSE, 
            aes(xmin=c0, xmax=c1, ymin=-20,
                ymax=30, group=group), color="transparent", fill="aquamarine", alpha=0.35)+
  annotate(x=1,y=25,label=expression(Q[low]),geom="label")+
  annotate(x=2,y=25,label=expression(Q[50]),geom="label")+
  annotate(x=3,y=25,label=expression(Q[high]),geom="label")


#######################################################################################
######################################################################################
############################################################################################

###Conformity scores computation
Eo.vec <- Em.vec <- Er.vec <- rep(NA,length(Y1a))
for (t in 1:(length(Y1a))){
    Eo.vec[t]   <- max(q.lo[t]-Y1a[t],Y1a[t]-q.hi[t])
    Em.vec[t]   <- max((q.lo[t]-Y1a[t])/(q.50[t]-q.lo[t]),(Y1a[t]-q.hi[t])/(q.hi[t]-q.50[t]))
    Er.vec[t]   <- max((q.lo[t]-Y1a[t])/(q.hi[t]-q.lo[t]),(Y1a[t]-q.hi[t])/(q.hi[t]-q.lo[t]))
  }
  
  k     <- ceiling((1-alpha.sig)*(1+length(Y1a)))
  Q.Eo  <- sort(Eo.vec)[k]
  Q.Em  <- sort(Em.vec)[k]
  Q.Er  <- sort(Er.vec)[k]
  
  tq.test.lo <- cbind(1,Xnew)%*%beta.lo
  tq.test.50 <- cbind(1,Xnew)%*%beta.50
  tq.test.hi <- cbind(1,Xnew)%*%beta.hi
  
  qs.test <- t(apply(cbind(tq.test.lo,tq.test.50,tq.test.hi),1,sort))
  
  q.test.lo <- qs.test[,1]
  q.test.50 <- qs.test[,2]
  q.test.hi <- qs.test[,3] 
  
  lb.o  <- q.test.lo - Q.Eo
  ub.o  <- q.test.hi + Q.Eo
  lb.m  <- q.test.lo - Q.Em * (q.test.50-q.test.lo) 
  ub.m  <- q.test.hi + Q.Em * (q.test.hi-q.test.50)  
  lb.r  <- q.test.lo - Q.Er * (q.test.hi-q.test.lo)  
  ub.r  <- q.test.hi + Q.Er * (q.test.hi-q.test.lo) 
  
  cov.o <- (Ynew<=ub.o & Ynew>=lb.o)
  cov.m <- (Ynew<=ub.m & Ynew>=lb.m)
  cov.r <- (Ynew<=ub.r & Ynew>=lb.r)
  
  leng.o <- ub.o-lb.o
  leng.m <- ub.m-lb.m
  leng.r <- ub.r-lb.r
  
##################################################################################
###################################################################################
qg<-data.frame(qs.test)
names(qg)<-c("q.lo","q.50","q.hi")
qgt<-data.frame(t(qg))
qgt$no<-1:3
names(qgt)<-c("qg","no")
qgt$label<-c("qlo","q50","qhi")
#######################################################################################
######################################################################################
#######################################################################################

###Conformity scores for CQR
ho<-ggplot2::ggplot(qgt,aes(x=no,y=qg,colour=label))+geom_point(size=3)+
   xlab("Number of values of quantile")+ylab("Q.Eo & coefficients for test set")+
   geom_hline(yintercept =lb.o,col="darkred",lwd=1)+
   geom_hline(yintercept=ub.o,col="darkred",lwd=1)+
  geom_segment(aes(x=2.5,xend=2.5,y=lb.o,yend=ub.o),col="purple",lwd=1,linetype="dashed")+
   geom_segment(aes(x=1,xend=1,y=lb.o,yend=qgt[1,1]),col="hotpink",lwd=1,linetype="dashed")+
   geom_segment(aes(x=3,xend=3,y=qgt[3,1],yend=ub.o),col="hotpink",lwd=1,linetype="dashed")+
   annotate(x=+Inf,y=lb.o,label="Lower bound",geom="label",hjust=1.3,vjust=1.05)+
   annotate(x=+Inf,y=ub.o,label="Upper bound",geom="label",hjust=1.3,vjust=-0.1)+
  geom_curve(
    aes(x = 2, y = (lb.o+qgt[1,1])/2, xend = 1, yend = (lb.o+qgt[1,1])/2   ),angle=180,curvature=0.5,
    arrow = arrow(length = unit(0.03, "npc"))
  )+
  geom_curve(
    aes(x = 2, y = (ub.o+qgt[3,1])/2, xend = 3, yend = (ub.o+qgt[3,1])/2    ),angle=180,curvature=0.5,
    arrow = arrow(length = unit(0.03, "npc"))
  )+
  annotate(x=2,y=(lb.o+qgt[1,1])/2,label="Q.Eo",geom="label")+
  annotate(x=2,y=(ub.o+qgt[3,1])/2,label="Q.Eo",geom="label")+
  annotate(x=2.55,y=4.5,label=expression(Length[CQR]:3.715477),geom="label",size=8/.pt)+
  scale_x_continuous(breaks= c(0,1,2,3,4))
  

ho  
###################################################################################################
###################################################################################################
###Conformity scores for CQR-m

hm<-ggplot2::ggplot(qgt,aes(x=no,y=qg,colour=label))+geom_point(size=3)+
  xlab("Number of values of quantile")+ylab("Q.Eo & coefficients for test set")+
  geom_hline(yintercept =lb.m,col="darkred",lwd=1)+
  geom_hline(yintercept=ub.m,col="darkred",lwd=1)+
  geom_segment(aes(x=2.5,xend=2.5,y=lb.o,yend=ub.o),col="purple",lwd=1,linetype="dashed")+
  geom_segment(aes(x=1,xend=1,y=lb.m,yend=qgt[1,1]),col="hotpink",lwd=1,linetype="dashed")+
  geom_segment(aes(x=3,xend=3,y=qgt[3,1],yend=ub.m),col="hotpink",lwd=1,linetype="dashed")+
  annotate(x=+Inf,y=lb.m,label="Lower bound",geom="label",hjust=1.3,vjust=1.05)+
  annotate(x=+Inf,y=ub.m,label="Upper bound",geom="label",hjust=1.3,vjust=-0.1)+
  geom_curve(
    aes(x = 2, y = (lb.m+qgt[1,1])/2, xend = 1, yend = (lb.m+qgt[1,1])/2   ),angle=180,curvature=0.5,
    arrow = arrow(length = unit(0.03, "npc"))
  )+
  geom_curve(
    aes(x = 2, y = (ub.m+qgt[3,1])/2, xend = 3, yend = (ub.m+qgt[3,1])/2    ),angle=180,curvature=0.5,
    arrow = arrow(length = unit(0.03, "npc"))
  )+
  annotate(x=2,y=(lb.m+qgt[1,1])/2,label="Q.Eo",geom="label")+
  annotate(x=2,y=(ub.m+qgt[3,1])/2,label="Q.Eo",geom="label")+
  annotate(x=2.55,y=4.5,label=expression(Length[CQR-m]:3.708842),geom="label",size=8/.pt)+
  scale_x_continuous(breaks= c(0,1,2,3,4))
  
hm 

##########################################################################################################
##########################################################################################################

###Confority scores for CQR-r
hr<-ggplot2::ggplot(qgt,aes(x=no,y=qg,colour=label))+geom_point(size=3)+
  xlab("Number of values of quantile")+ylab("Q.Eo & coefficients for test set")+
  geom_hline(yintercept =lb.r,col="darkred",lwd=1)+
  geom_hline(yintercept=ub.r,col="darkred",lwd=1)+
  geom_segment(aes(x=2.5,xend=2.5,y=lb.o,yend=ub.o),col="purple",lwd=1,linetype="dashed")+
  geom_segment(aes(x=1,xend=1,y=lb.r,yend=qgt[1,1]),col="hotpink",lwd=1,linetype="dashed")+
  geom_segment(aes(x=3,xend=3,y=qgt[3,1],yend=ub.r),col="hotpink",lwd=1,linetype="dashed")+
  annotate(x=+Inf,y=lb.r,label="Lower bound",geom="label",hjust=1.3,vjust=1.05)+
  annotate(x=+Inf,y=ub.r,label="Upper bound",geom="label",hjust=1.3,vjust=-0.1)+
  geom_curve(
    aes(x = 2, y = (lb.r+qgt[1,1])/2, xend = 1, yend = (lb.r+qgt[1,1])/2   ),angle=180,curvature=0.5,
    arrow = arrow(length = unit(0.03, "npc"))
  )+
  geom_curve(
    aes(x = 2, y = (ub.r+qgt[3,1])/2, xend = 3, yend = (ub.r+qgt[3,1])/2    ),angle=180,curvature=0.5,
    arrow = arrow(length = unit(0.03, "npc"))
  )+
  annotate(x=2,y=(lb.r+qgt[1,1])/2,label="Q.Eo",geom="label")+
  annotate(x=2,y=(ub.r+qgt[3,1])/2,label="Q.Eo",geom="label")+
  annotate(x=2.55,y=4.5,label=expression(Length[CQR-r]:3.713439),geom="label",size=8/.pt)+
  scale_x_continuous(breaks= c(0,1,2,3,4))

hr

############################################################################################################
############################################################################################################

######CQR#################################################################

###Prediction  intervals:CQR
hfo<-a+geom_hline(yintercept = lb.o,col="darkgreen",size=1.5)+
  geom_hline(yintercept=ub.o,col="darkgreen",size=1.5)+
  annotate(x=0.67,y=40,label="Xtest:0.67",geom="label",size=9/.pt)+
  annotate(x=-5,y=2.7,label="Ytest:2.7",vjust=-0.1,geom="label",size=9/.pt)+
  annotate(x=10,y=ub.o,label="Upper bound:5.7874",geom="label",vjust=0.1,size=9/.pt)+
  annotate(x=10,y=lb.o,label="Lower bound:2.0719",geom="label",size=9/.pt,
           vjust=1)
hfo

#################################################################################

#####CQR-m-prediction intervals

hfm<-a+geom_hline(yintercept = lb.m,col="darkgreen",size=1.5)+
  geom_hline(yintercept=ub.m,col="darkgreen",size=1.5)+
  annotate(x=0.67,y=40,label="Xtest:0.67",geom="label",size=9/.pt)+
  annotate(x=-5,y=2.7,label="Ytest:2.7",vjust=-0.1,geom="label",size=9/.pt)+
  annotate(x=10,y=ub.m,label="Upper bound:5.7836",geom="label",vjust=0.1,size=9/.pt)+
  annotate(x=10,y=lb.m,label="Lower bound:2.0747",geom="label",size=9/.pt,
           vjust=1)
hfm

#################################################################################

#####CQR-r-prediction intervals

hfr<-a+geom_hline(yintercept = lb.r,col="darkgreen",size=1.5)+
  geom_hline(yintercept=ub.r,col="darkgreen",size=1.5)+
  annotate(x=0.67,y=40,label="Xtest:0.67",geom="label",size=9/.pt)+
  annotate(x=-5,y=2.7,label="Ytest:2.7",vjust=-0.1,geom="label",size=9/.pt)+
  annotate(x=10,y=ub.r,label="Upper bound:5.7863",geom="label",vjust=0.1,size=9/.pt)+
  annotate(x=10,y=lb.r,label="Lower bound:2.0729",geom="label",size=9/.pt,
           vjust=1)
hfr

##################################################################################
###################################################################################





