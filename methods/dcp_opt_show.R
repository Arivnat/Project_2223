#############################################################
####Distributional conformal prediction:opt#################
##############################################################

####################################################################################
# Paper: Distributional conformal prediction
# Authors: V. Chernozhukov, K. Wuthrich, Y. Zhu,PNAS,118(48)
##Relevant github:
#https://github.com/kwuthrich/Replication_DCP
###################################################################################
####The above github is source material for the functions for the different 
###methods-DCP/DCP-QR,DCP-opt,CQR,CQR-m,CQR-r and CP-reg
#The functions for the methods, the function to create train,calibration
##and test data, the loop to check for singularity of the matrix as well 
##as well as the function for quantile bins have been taken from the github
##and used for this thesis.
#The rest of the code(data generation,setting up replications, additional plots,descriptive 
#analysis) has been written by me
########################################################################################
#######################################################################################
#########################################################
###loading required packages
library(quantreg)
library(matrixcalc)
library(tidyverse)
library(ggplot2)

#############################################################
set.seed(1)

aux.uneven <- function(Y,X){
  T01     <- length(Y)
  T01even <- floor(T01/2)*2
  
  Yeven <- Y[(T01-T01even+1):T01]
  Xeven <- cbind(X[(T01-T01even+1):T01,])
  
  ind0 <- 1:(T01even/2)
  ind1 <- (T01even/2+1):T01even
  
  Y0 <- Yeven[ind0]
  X0 <- cbind(Xeven[ind0,])
  Y1 <- Yeven[ind1]
  X1 <- cbind(Xeven[ind1,])
  
  return(list(Y0=Y0,X0=X0,Y1=Y1,X1=X1))
  
}

########################################################################
##########################################################################

#############################################################################
#############################################################################
###Set up random data
#############################################################################
set.seed(1)
n<-1000
X<-rnorm(n,1,3)
U<-rnorm(n,0,1)
Y<-2+3*X+U
datXY<-cbind.data.frame(X,Y)
###data point to predict
Xnew<-0.67
Ynew<-2.7

library(matrixcalc)
sing<-TRUE

#######checking if matrix is singular
####################################################################

while (sing==TRUE){
  set.seed<-1
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
######################################################################
XXX <- rbind(X1a,Xnew)
YYY <- c(Y1a,Ynew)

train<-cbind.data.frame(Xa,Ya)
calib<-cbind.data.frame(X1a,Y1a)

#######################################################################
point<-c(Xnew,Ynew)
newDat<-rbind.data.frame(datXY,point)
names(newDat)<-c("X","Y")


a<-ggplot(data=newDat,aes(x=X,y=Y))+geom_point()+
  geom_vline(xintercept=0.67,col="red",linetype="dashed",linewidth=1)+
  geom_hline(yintercept=2.7,col="blue",linetype="dashed",linewidth=1)

a1<-a+annotate(x=0.67,y=+Inf,label="Xtest:0.67",vjust=2,geom="label")
a2<-a1+annotate(x=+Inf,y=2.7,label="Ytest:2.7",hjust=1.25,geom="label")
a3<-a2+annotate(x=+Inf,y=+Inf,label="Entire Data",vjust=5,hjust=1,
                geom="label")

######################################################################
####Visualizing the training phase-quantile regression on training data
######################################################################

taus    <- seq(0.001,0.999,length=200)  
beta.qr <- matrix(NA,dim(Xa)[2]+1,length(taus))
  for (t in 1:length(taus)){
    beta.qr[,t] <- rq.fit.br(cbind(1,Xa),Ya,tau=taus[t])$coefficients
  }

#######################################################################
#######################################################################

go<-data.frame(beta.qr[,c(20,40,60,80,100,120,140,160,180,200)])
go<-data.frame(t(go))

b3<-ggplot2::ggplot(train,aes(x=Xa,y=Ya))+geom_point(col="darkgreen")+
  xlab("X(training data)")+ylab("Y(training data)")+
  geom_abline(intercept=go[1,1],slope=go[1,2])+
  geom_abline(intercept=go[2,1],slope=go[2,2])+
  geom_abline(intercept=go[3,1],slope=go[3,2])+
  geom_abline(intercept=go[4,1],slope=go[4,2])+
  geom_abline(intercept=go[5,1],slope=go[5,2])+
  geom_abline(intercept=go[6,1],slope=go[6,2])+
  geom_abline(intercept=go[7,1],slope=go[7,2])+
  geom_abline(intercept=go[8,1],slope=go[8,2])+
  geom_abline(intercept=go[9,1],slope=go[9,2])+
  geom_abline(intercept=go[10,1],slope=go[10,2])

#######################################################################

###Visualizing estimated cdf for calibration and test data

########################################################################

###covariates(Xcalib,Xnew)'betaqr
tQ.yx   <- cbind(1,XXX)%*%beta.qr

##Sorting estimate row-wise(observation-wise) to ensure estimated cdf is monotonic
################################################################################

Q.yx    <- t(apply(tQ.yx,1,FUN=sort))  ##501 200

###Computed ranks
u.hat <- rowMeans((Q.yx <= matrix(YYY,length(YYY),dim(beta.qr)[2])))
  
########################################################################

q<-data.frame(Q.yx) 
qt<-data.frame(t(q))  ##200 501
qt<-data.frame(qt)
qt$no<-1:200

library(latex2exp)  
f<-ggplot2::ggplot(qt,mapping=aes(no,qt[,20]))+geom_point()+
    geom_point(aes(x=no,y=qt[,40]),col="#CC0000")+
    geom_point(aes(x=no,y=qt[,60]),col="#CC9900")+
    geom_point(aes(x=no,y=qt[,80]),col="#0066CC")+
    geom_point(aes(x=no,y=qt[,100]),col="#336600")+
    geom_point(aes(x=no,y=qt[,120]),col="#00FF00")+
    geom_point(aes(x=no,y=qt[,140]),col="#330099")+
    geom_point(aes(x=no,y=qt[,160]),col="#0099CC")+
    geom_point(aes(x=no,y=qt[,180]),col="#FF3300")+
    geom_point(aes(x=no,y=qt[,200]),col="#FF33CC")+
    xlab("Total number of unique values of quantile")+
    ylab(TeX(r'($(X_{calib,test}$'$\hat{\beta}$) for each  value of quantile-sorted by observation)'))
             
f
##################################################################################
##################################################################################

qt2<-qt[,-502]
yg<-matrix(t(matrix(YYY,length(YYY),length(taus))),nrow=200,ncol=501)
yg<-data.frame(yg)

g<-ggplot2::ggplot(qt,mapping=aes(no,qt[,20]))+geom_point(col="#333333")+
  geom_hline(yintercept=yg[1,20], lwd=1,linetype="dashed")+
  geom_point(aes(x=no,y=qt[,40]),col="#CC0000")+
  geom_hline(yintercept=yg[1,40],lwd=1,col="#990000",linetype="dashed")+
 # geom_point(aes(x=no,y=qt[,60]),col="#CC9900")+
  #geom_hline(yintercept=yg[1,60],lwd=1,col="#993300",linetype="dashed")+
  # geom_point(aes(x=no,y=qt[,80]),col="#0066CC")+
  #  geom_hline(yintercept=yg[1,80],lwd=1,col="#003366")+
  #geom_point(aes(x=no,y=qt[,100]),col="#336600")+
  #geom_hline(yintercept=yg[1,100],lwd=1,col="#003300")+
  geom_point(aes(x=no,y=qt[,120]),col="#00FF00")+
  geom_hline(yintercept=yg[1,120],lwd=1,col="#006600",linetype="dashed")+
  #geom_point(aes(x=no,y=qt[,140]),col="#330099")+
  #geom_hline(yintercept=yg[1,140],lwd=1,col="#003366",linetype="dashed")+
  geom_point(aes(x=no,y=qt[,160]),col="#0099CC")+
  geom_hline(yintercept=yg[1,160],lwd=1,col="#3399FF",linetype="dashed")+
  geom_point(aes(x=no,y=qt[,180]),col="#FF3300")+
  geom_hline(yintercept=yg[1,180],lwd=1,col="#660000",linetype="dashed")+
  geom_point(aes(x=no,y=qt[,501]),col="#FF33CC")+
  geom_hline(yintercept=yg[1,501],lwd=1,col="#9900FF",linetype="dashed")+
  xlab("Number of values of quantile")+
  ylab(TeX(r'(Sorted $X_{calib,test}$'$\hat{\beta}$ for each  value of quantile & corresponding response values)'))


#########################################################################################
g1<-g+annotate(x=75,y=yg[1,20],label="Y20",geom="label",
               size=8/.pt)+
  annotate(x=100,y=yg[1,40],label="Y40",geom="label",
           size=8/.pt)+
  #annotate(x=+Inf,y=yg[1,60],label=expression(Y[60]),hjust=30,vjust=0.25,geom="label",
   #        size=12/.pt)+
  annotate(x=150,y=yg[1,120],label="Y120",geom="label",
           size=8/.pt)+
  #annotate(x=+Inf,y=yg[1,140],label=expression(Y[140]),hjust=24,geom="label",
   #        size=12/.pt)+
  annotate(x=150,y=yg[1,160],label="Y160",geom="label",
           size=8/.pt)+
  annotate(x=150,y=yg[1,180],label="Y180",geom="label",
           size=8/.pt)+
  annotate(x=150,y=yg[1,501],label="Ynew",geom="label",
           size=8/.pt)+
  annotate(x=+Inf,y=qt[198,20],label=TeX("${X}^{T}$$\\hat{\\beta}_{train}$:20th obs"),hjust=1.5,geom="label",
           size=8/.pt)+
  annotate(x=+Inf,y=qt[198,40],label=TeX("${X}^{T}$$\\hat{\\beta}_{train}$:40th obs"),hjust=2,
           vjust=1.5,geom="label",
           size=8/.pt)+
  #annotate(x=+Inf,y=qt[198,60],label=TeX("${X}^{T}$$\\hat{\\beta}_{train}$:60th obs"),hjust=1.5,geom="label",
   #        size=8/.pt)+
  annotate(x=+Inf,y=qt[198,120],label=TeX("${X}^{T}$$\\hat{\\beta}_{train}$:120th obs"),hjust=1.5,geom="label",
           size=8/.pt)+
 # annotate(x=+Inf,y=qt[198,140],label=TeX("${X}^{T}$$\\hat{\\beta}_{train}$:140th obs"),hjust=1.5,geom="label",
  #         size=8/.pt)+
  annotate(x=+Inf,y=qt[198,160],label=TeX("${X}^{T}$$\\hat{\\beta}_{train}$:160th obs"),hjust=2,
           vjust=2,geom="label",
           size=8/.pt)+
  annotate(x=+Inf,y=qt[198,180],label=TeX("${X}^{T}$$\\hat{\\beta}_{train}$:180th obs"),hjust=1,
           vjust=0.5,geom="label",
           size=8/.pt)+
  annotate(x=+Inf,y=qt[198,501],label=TeX("${X}^{T}$$\\hat{\\beta}_{train}$:$X_{new}$"),hjust=1.5,geom="label",
           size=8/.pt)

g1
#########################################################################################
#########################################################################################

##qt 200 502 
qt2<-qt[,-502]  ##501
qt2<=yg
qt3<-qt2[,c(20,40,120,140,160,501)]
qt3$no<-1:200
yg2<-yg[,c(20,40,120,140,160,501)]
names(yg2)<-c("Y20","Y40","Y120","Y140","Y160","Y501")

qt3$Y20<-yg2[,1]
qt3$Y40<-yg2[,2]
qt3$Y120<-yg2[,3]
qt3$Y140<-yg2[,4]
qt3$Y160<-yg2[,5]
qt3$Y501<-yg2[,6]

qt4<-qt2[,c(20,40,120,140,160,501)]
S<-colSums(qt4<=yg2)

##############################################
### X.19  X.39 X.119 X.139 X.159  Xnew 
#####141   176     2   100    82    21

##############################################################################
################################################################################

g3<-ggplot2::ggplot(qt3,mapping=aes(x=no,y=X.19))+geom_line(col="#333333")+
  geom_hline(yintercept=qt3$Y20,lwd=1,linetype="dashed")+
  geom_ribbon(data=subset(qt3,X.19<=qt3$X.19[141]),
              aes(ymin=X.19,ymax=Y20),fill="grey")+
  
  geom_line(aes(x=no,y=X.39),col="#CC0000")+
  geom_hline(yintercept=qt3$Y40,lwd=1,col="#990000",linetype="dashed")+
  geom_ribbon(data=subset(qt3,X.39<=qt3$X.39[176]),
              aes(ymin=X.39,ymax=Y40),fill="pink")+
  
  geom_line(aes(x=no,y=X.119),col="#00FF00")+
  geom_hline(yintercept=qt3$Y120,lwd=1,col="#006600",linetype="dashed")+
  geom_ribbon(data=subset(qt3,X.119<=qt3$X.119[2]),
              aes(ymin=X.119,ymax=Y120),fill="yellow")+
  
  
  geom_line(aes(x=no,y=X.159),col="#0099CC")+
  geom_hline(yintercept=qt3$Y160,lwd=1,col="#3399FF",linetype="dashed")+
  geom_ribbon(data=subset(qt3,X.159<=qt3$X.159[82]),
              aes(ymin=X.159,ymax=Y160),fill="lightblue")+
  
  geom_line(aes(x=no,y=Xnew),col="#CC66FF")+
  geom_hline(yintercept=qt3$Y501,lwd=1,col="#3300FF",linetype="dashed")+
  geom_ribbon(data=subset(qt3,0<=Xnew & Xnew<=qt3$Xnew[21]),
              aes(ymin=Xnew,ymax=Ynew),fill="#FF66CC")+
  
  xlab("Number of values of quantile")+
  ylab(TeX(r'(Sorted $X_{calib,test}$'$\hat{\beta}$ for each  value of quantile & corresponding response values)'))

g4<-g3+annotate(x=200,y=yg[1,20],label="Y20",geom="label",
                size=8/.pt)+
  annotate(x=150,y=yg[1,40]+1,label="Y40",geom="label",
           size=8/.pt)+
  annotate(x=150,y=yg[1,120],label="Y120",geom="label",
           size=8/.pt)+
  annotate(x=150,y=yg[1,160],label="Y160",geom="label",
           size=8/.pt)+
  annotate(x=150,y=yg[1,501],label="Ynew",geom="label",
           size=8/.pt)+
  geom_curve(
    aes(x = 50, y = (yg[1,40]-0.25), xend = 75, yend = -5   ),angle=90,curvature=0.5,
    arrow = arrow(length = unit(0.03, "npc")))+
  annotate(x=75,y=-3.5,label=TeX("\\hat{U}:\\frac{#(X<=Y40)}{200}"),geom="label",
           size=8/.pt)+
  geom_curve(
    aes(x = 20, y = (yg[1,160])-0.05, xend = 50, yend = -12   ),angle=90,curvature=0.5,
    arrow = arrow(length = unit(0.03, "npc")))+
  annotate(x=70,y=-13 ,label=TeX("\\hat{U}:\\frac{#(X<=Y160)}{200}"),geom="label",
           size=8/.pt)+
  geom_curve(
    aes(x = 50, y =yg[1,20]-0.05 , xend = 80, yend = -20   ),angle=90,curvature=0.5,
    arrow = arrow(length = unit(0.03, "npc")))+
  annotate(x=100,y=-20 ,label=TeX("\\hat{U}:\\frac{#(X<=Y20)}{200}"),geom="label",
           size=8/.pt)+
  geom_curve(
    aes(x = 15, y =yg[1,501]-0.1 , xend = 35, yend = 0   ),angle=90,curvature=0.5,
    arrow = arrow(length = unit(0.03, "npc")))+
  annotate(x=52,y=0 ,label=TeX("\\hat{U}:\\frac{#(X<=Ynew)}{200}"),geom="label",
           size=8/.pt)

g4
##################################################################
##################################################################

####################################################################
###Computing upper and lower bound of prediction interval 
#####################################################################
alpha.sig<-0.1
bhat <- rep(NA,length(YYY))
b.grid <- taus[taus<=alpha.sig]  ##20 values 

###
sum(taus<=alpha.sig) ##20

l<-list()
  for (t in 1:length(YYY)){
    leng <- rep(NA,length(b.grid))
    leng.test <- rep(NA,length(b.grid))
    for (b in 1:length(b.grid)){
      Q.yx.u <- approx(x=taus,y=Q.yx[t,],xout=(b.grid[b]+1-alpha.sig),rule=2)$y
      leng[b] <- Q.yx.u -Q.yx[t,b]
    }
    bhat[t] <- b.grid[which.min(leng)]
  }
  
####################################################################################
#####################################################################################
r1<-Q.yx[1,1:20]

leng <- rep(NA,length(b.grid))
leng.test <- rep(NA,length(b.grid))
for (b in 1:length(b.grid)){
  Q.yx.u <- approx(x=taus,y=Q.yx[1,],xout=(b.grid[b]+1-alpha.sig),rule=2)$y
  leng[b] <- Q.yx.u -Q.yx[1,b]
}
Q.yx.u
m <- b.grid[which.min(leng)]  

leng <- rep(NA,length(b.grid))
leng.test <- rep(NA,length(b.grid))
for (b in 1:length(b.grid)){
  Q.yx.u20 <- approx(x=taus,y=Q.yx[20,],xout=(b.grid[b]+1-alpha.sig),rule=2)$y
  leng[b] <- Q.yx.u20 -Q.yx[1,b]
}
Q.yx.u20



B<-cbind.data.frame(taus,Q.yx[1,],1:200,Q.yx[20,])
names(B)<-c("taus","Qyx1","no","Qyx20")
highlight<- data.frame(start=1, end=20, group=seq_along(20))

n<-ggplot2::ggplot(B,mapping=aes(x=no,y=Qyx1))+geom_point()+
   geom_line(B,mapping=aes(x=no,y=taus),col="red",lwd=1)+
  annotate(x=150,y=taus[150],label="Quantiles",geom="label",
           size=8/.pt,hjust=-0.05)+
   geom_hline(yintercept=Q.yx.u,col="blue",lwd=1)+
   annotate(x=20.5,y=(Q.yx.u+min(B$Qyx1))/1.75,label=TeX("$\\hat{b}1$;min.length:1st obs"),geom="label",
            size=8/.pt,hjust=-0.05)+
  annotate(x=150,y=Q.yx.u,label="Upper bound: 1st obs",geom="label",
           size=8/.pt)+
  annotate(x=125,y=B$Qyx1[125],label=TeX("Sorted ${X}^{T}$$\\beta$:1st obs"),geom="label",
           size=8/.pt)+
  geom_rect(data=highlight,inherit.aes=FALSE, 
  aes(xmin=1, xmax=20, ymin=-20,
                ymax=10, group=group), color="transparent", fill="magenta", alpha=0.25)+
  geom_segment(aes(x=20,xend=20,y=Q.yx[1,20],yend=Q.yx.u),col="darkred",lwd=1,linetype="dashed")+
  geom_point(B,mapping=aes(x=no,y=Qyx20),col="#3399FF")+
  annotate(x=125,y=B$Qyx20[125],label=TeX("Sorted ${X}^{T}$$\\beta$:20th obs"),geom="label",
           size=8/.pt)+
  geom_hline(yintercept=Q.yx.u20,col="#339900",lwd=1)+
  annotate(x=150,y=Q.yx.u20,label="Upper bound: 20th obs",geom="label",
           size=8/.pt)+
  geom_segment(aes(x=20,xend=20,y=Q.yx[20,20],yend=Q.yx.u20),col="darkred",lwd=1,linetype="dashed")+
  annotate(x=20.5,y=(Q.yx.u20+min(B$Qyx20))/2.25,label=TeX("$\\hat{b}20$;min.length:20th obs"),geom="label",
           size=8/.pt,hjust=-0.05)+
  xlab("Number of values of quantile")+
  ylab(TeX(r'(Sorted $X_{calib,test}$'$\hat{\beta}$ and quantiles)'))
           

n
#################################################################################
#################################################################################
D<-cbind.data.frame(taus,bhat[501],1:200,0.9/2)
names(D)<-c("taus","bhat","no","sig")
library(latex2exp)

h<-ggplot2::ggplot(D,mapping=aes(x=no,y=taus))+geom_line(col="purple",lwd=1)+
   geom_hline(yintercept = D$bhat,col="darkred",lwd=1)+
   geom_hline(yintercept = D$sig,col="#FF9933",lwd=1)+
   geom_hline(yintercept = D$sig+D$bhat,col="#FF0066",lwd=1)+
  geom_segment(aes(x=200,xend=200,y=sig+bhat,yend=taus[200]),col="#660000",lwd=1,linetype="dashed")+
  annotate(x=184,y=(D$sig+D$bhat+taus[200])/2-0.02,
  label=TeX("$\\tau-\\hat{b}_{Xtest}-\\frac{(1-\\alpha)}{2}$"),geom="label",size=8/.pt)+
  annotate(x=132,y=D$sig+D$bhat,label=TeX("$\\hat{b}_{Xtest}+\\frac{(1-\\alpha)}{2}$"),
           geom="label",size=8/.pt)+
  annotate(x=184,y=D$sig,label=TeX("$\\frac{(1-\\alpha)}{2}$"),geom="label",size=8/.pt)+
  annotate(x=125,y=D$bhat,label=TeX("$\\hat{b}_{Xtest}$"),geom="label",size=8/.pt)+
  xlab("Number of values of quantile")+
  ylab(TeX(r'($\tau,\hat{b}_{Xtest},\frac{1-\alpha}{2}$)'))


  
################################################################################
################################################################################
   
ind.test <- (length(Y1a)+1):length(YYY)
  
cs.opt <- abs(u.hat-bhat-(1-alpha.sig)/2)

#################################################################################
#################################################################################
k           <- ceiling((1-alpha.sig)*(1+length(Y1a)))
threshold   <- sort(cs.opt[-ind.test])[k]

cov.opt   <- (cs.opt[ind.test] <= threshold)

################################################################################
################################################################################
ind.test<-501
leng.opt <- NULL
for (t in ind.test){
   ci.grid <- abs(taus - bhat[t]-(1-alpha.sig)/2)
   ci <- Q.yx[t,(ci.grid <= threshold)]
   ub <- max(ci)
   lb <- min(ci)
   leng.opt <- c(leng.opt,ub-lb)
  }
leng.opt
ind.test
  
###################################################################################
##############Visualizing upper and lower bound of prediction intervals as 
###############well as its construction###########################################
###################################################################################

D2<-cbind.data.frame(ci.grid,threshold,Q.yx[501,],1:200)
names(D2)<-c("grid","threshold","Qyx50","no")
D2$grid<=D2$threshold  ##9---195

h2<-ggplot2::ggplot(D2,aes(x=no,y=grid))+geom_point()+
    geom_line(D2,mapping=aes(x=no,y=Qyx50),col="orange",size=1)+
  geom_line(D2,mapping=aes(x=no,y=threshold),col="darkred",size=1)+
  #geom_vline(xintercept=9,col="magenta")+
  #geom_vline(xintercept=195,col="magenta" )+
  geom_segment(aes(x=9,xend=9,y=threshold,yend=Qyx50[9],col="magenta"))+
  geom_segment(aes(x=195,xend=195,y=threshold,yend=Qyx50[195],col="magenta"))+
  annotate(x=28,y=(D2$Qyx50[9]+D2$threshold[9])/2,label="Lower Bound",geom="label",size=8/.pt)+
  annotate(x=176,y=(D2$Qyx50[195]+D2$threshold[195])/2,label="Lower Bound",geom="label",size=8/.pt)+
  geom_segment(aes(x=9,xend=195,y=2,yend=2,col="red"))+
  annotate(x=100,y=2,label="Length of prediction interval",geom="label",size=8/.pt)+
  xlab("Number of values of quantile")+
  ylab(TeX("Sorted $\\X_{test}$'$\\beta$, threshold, grid"))+
  annotate(x=100,y=threshold,label="threshold",geom="label",size=8/.pt)+
  annotate(x=150,y=D2$grid[150],label="Grid",geom="label",size=8/.pt)+
  annotate(x=150,y=D2$Qyx50[150],label=TeX("Sorted $\\{X}^{T}_{test}$$\\beta$"),geom="label",size=8/.pt)
######################################################################################  
######################################################################################

###########Constructed prediction interval
#####################################################################################
#####################################################################################

a4<-a+geom_hline(yintercept = lb,col="darkgreen",size=1.5)+
  geom_hline(yintercept=ub,col="darkgreen",size=1.5)+
  annotate(x=0.67,y=40,label="Xtest:0.67",geom="label",size=9/.pt)+
  annotate(x=-5,y=2.7,label="Ytest:2.7",vjust=-0.1,geom="label",size=9/.pt)+
  annotate(x=9,y=ub,label=TeX("$\\Upper bound_{DCP-opt}$:5.712972"),geom="label",vjust=0.1,size=9/.pt)+
  annotate(x=9,y=lb,label=TeX("$\\Lower bound_{DCP-opt}$:2.230669"),geom="label",size=9/.pt,
           vjust=1)
a4
  





    