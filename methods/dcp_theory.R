##Distributional Conformal prediction

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

###################################################################
###Function to split the data into training and calibration sets###
###################################################################


library(ggplot2)
library(ggpubr)
library(cowplot)

library(matrixcalc)
library(quantreg)

devtools::install_github('stefano-meschiari/latex2exp')
install.packages("latex2exp")
library(latex2exp)


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

while (sing==TRUE){
  set.seed(1)
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




#############################################################################
############################################################################

#dcp.qr <- function(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig){



###################################################################################
###Step-1: Splitting data into training and calibration sets 
###################################################################################
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

train<-cbind.data.frame(Xa,Ya)
calib<-cbind.data.frame(X1a,Y1a)

b<-ggplot(train,aes(x=Xa,y=Ya))+geom_point(col="darkgreen")+
  xlab("X(training data)")+ylab("Y(training data)")
c<-ggplot(calib,aes(x=X1a,y=Y1a))+geom_point(col="orange")+
  xlab("X(calibration data)")+ylab("Y(calibration data)")



b1<-b+annotate(x=+Inf,y=+Inf,label="Training Data",vjust=5,hjust=1,geom="label")
c1<-c+annotate(x=+Inf,y=+Inf,label="Calibration Data",vjust=5,hjust=1,geom="label")

plot_grid(
  a3, b1,c1,
  align = "h", axis = "tb",
  nrow = 1, rel_widths = c(1, 1)
)


###Step-2: Performing quantile regression on training data 
######################################################################
library(quantreg)
taus    <- seq(0.001,0.999,length=200)  

beta.qr <- matrix(NA,dim(Xa)[2]+1,length(taus))
for (t in 1:length(taus)){
  beta.qr[,t] <- rq.fit.br(cbind(1,Xa),Ya,tau=taus[t])$coefficients
}

g<-data.frame(beta.qr[,c(20,40,60,80,100,120,140,160,180,200)])
g<-data.frame(t(g))

b3<-ggplot2::ggplot(train,aes(x=Xa,y=Ya))+geom_point(col="darkgreen")+
  xlab("X(training data)")+ylab("Y(training data)")+
  geom_abline(intercept=g[1,1],slope=g[1,2])+
  geom_abline(intercept=g[2,1],slope=g[2,2])+
  geom_abline(intercept=g[3,1],slope=g[3,2])+
  geom_abline(intercept=g[4,1],slope=g[4,2])+
  geom_abline(intercept=g[5,1],slope=g[5,2])+
  geom_abline(intercept=g[6,1],slope=g[6,2])+
  geom_abline(intercept=g[7,1],slope=g[7,2])+
  geom_abline(intercept=g[8,1],slope=g[8,2])+
  geom_abline(intercept=g[9,1],slope=g[9,2])+
  geom_abline(intercept=g[10,1],slope=g[10,2])
b3


###Step 3: Using estimates from quantile regression on training data 
##to multiply with the data on the calibration set
###############################################################################

tQ.yx   <- cbind(1,X1a)%*%beta.qr

Q.yx    <- t(apply(tQ.yx,1,FUN=sort))

q<-data.frame(Q.yx) 
qt<-data.frame(t(q))
qt<-data.frame(qt)
qt$no<-1:200

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
  ylab(TeX(r'(Sorted  ${\X}^{T}$${\beta}$ per observation for each unique value of quantile)'))

f

################################################################################
qt2<-qt[,-501]
yg<-matrix(t(matrix(Y1a,length(Y1a),length(taus))),nrow=200,ncol=500)
yg<-data.frame(yg)

g<-ggplot2::ggplot(qt,mapping=aes(no,qt[,20]))+geom_point(col="#333333")+
  geom_hline(yintercept=yg[1,20], lwd=1,linetype="dashed")+
  geom_point(aes(x=no,y=qt[,40]),col="#CC0000")+
  geom_hline(yintercept=yg[1,40],lwd=1,col="#990000",linetype="dashed")+
  geom_point(aes(x=no,y=qt[,60]),col="#CC9900")+
  geom_hline(yintercept=yg[1,60],lwd=1,col="#993300",linetype="dashed")+
  # geom_point(aes(x=no,y=qt[,80]),col="#0066CC")+
  #  geom_hline(yintercept=yg[1,80],lwd=1,col="#003366")+
  #geom_point(aes(x=no,y=qt[,100]),col="#336600")+
  #geom_hline(yintercept=yg[1,100],lwd=1,col="#003300")+
  geom_point(aes(x=no,y=qt[,120]),col="#00FF00")+
  geom_hline(yintercept=yg[1,120],lwd=1,col="#006600",linetype="dashed")+
  geom_point(aes(x=no,y=qt[,140]),col="#330099")+
  geom_hline(yintercept=yg[1,140],lwd=1,col="#003366",linetype="dashed")+
  geom_point(aes(x=no,y=qt[,160]),col="#0099CC")+
  geom_hline(yintercept=yg[1,160],lwd=1,col="#3399FF",linetype="dashed")+
  geom_point(aes(x=no,y=qt[,180]),col="#FF3300")+
  geom_hline(yintercept=yg[1,180],lwd=1,col="#660000",linetype="dashed")+
  #geom_point(aes(x=no,y=qt[,200]),col="#FF33CC")+
  #geom_hline(yintercept=yg[1,200],lwd=1,col="#9900FF")+
  xlab("Number of values of quantile")+
  ylab(TeX(r'(Sorted ${\X}^{T}$$\beta$ for each observation,corresponding response values)'))

g1<-g+annotate(x=+Inf,y=yg[1,20],label=expression(Y[20]),hjust=30,vjust=-0.05,geom="label",
               size=12/.pt)+
  annotate(x=175,y=yg[1,40],label=expression(Y[40]),hjust=31,vjust=-0.05,geom="label",
           size=12/.pt)+
  annotate(x=+Inf,y=yg[1,60],label=expression(Y[60]),hjust=30,vjust=0.25,geom="label",
           size=12/.pt)+
  annotate(x=+Inf,y=yg[1,120],label=expression(Y[120]),hjust=24,geom="label",
           size=12/.pt)+
  annotate(x=+Inf,y=yg[1,140],label=expression(Y[140]),hjust=24,geom="label",
           size=12/.pt)+
  annotate(x=+Inf,y=yg[1,160],label=expression(Y[160]),hjust=24,geom="label",
           size=12/.pt)+
  annotate(x=+Inf,y=yg[1,180],label=expression(Y[180]),hjust=28,vjust=-0.25,geom="label",
           size=12/.pt)+
  annotate(x=175,y=qt[175,20],label=TeX("${\\X}_{20}^{T}$$\\beta$"),vjust=-0.1,geom="label",
           size=12/.pt)+
  annotate(x=175,y=qt[175,40],label=TeX("${\\X}_{40}^{T}$$\\beta$"),geom="label",
           size=12/.pt,vjust=1.2)+
  annotate(x=75,y=qt[75,60],label=TeX("${\\X}_{60}^{T}$$\\beta$"),geom="label",
           size=12/.pt)+
  annotate(x=75,y=qt[75,120],label=TeX("${\\X}_{120}^{T}$$\\beta$"),geom="label",
           size=12/.pt,vjust=1.2)+
  annotate(x=175,y=qt[175,140],label=TeX("${\\X}_{140}^{T}$$\\beta$"),hjust=2,geom="label",
           size=12/.pt)+
  annotate(x=175,y=qt[175,160],label=TeX("${\\X}_{160}^{T}$$\\beta$"),hjust=2,geom="label",
           size=12/.pt)+
  annotate(x=160,y=qt[160,180],label=TeX("${\\X}_{180}^{T}$$\\beta$"),geom="label",
           size=12/.pt)

g1

plot_grid(b3,f)


##
##taus[20];taus[40];taus[60];taus[120];taus[140];taus[160];taus[180]
################################################################################
################################################################################
#################################################################################


u.hat   <- rowMeans((Q.yx <= matrix(Y1a,length(Y1a),length(taus))))
cs      <- abs(u.hat-0.5)
#################################################################################
#################################################################################

qt2<-qt[,-501]
qt2<=yg
qt3<-qt2[,c(20,40,120,140,160)]
qt3$no<-1:200
yg2<-yg[,c(20,40,120,140,160)]
names(yg2)<-c("Y20","Y40","Y120","Y140","Y160")
yg2$no<-1:200
library(dplyr)
yg2<-yg2%>%select(no,Y20,Y40,Y120,Y140,Y160)

dfjoin <- left_join(qt3,yg2, by = "no")

###left join 
qt4<-qt3[,-6]
yg3<-yg2[,-1]
C<-colSums(qt4<=yg3)  ##gives maximum index upto which 
###sorted coefficients are less than corresponding y value 
C

#X20  X40 X120 X140 X160 
#141  176    2  100   82 

qt3$Y20<-yg2[,2]
qt3$Y40<-yg2[,3]
qt3$Y120<-yg2[,4]
qt3$Y140<-yg2[,5]
qt3$Y160<-yg2[,6]

#X20  X40 X120 X140 X160 
#10    67  72   87   96

g3<-ggplot2::ggplot(qt3,mapping=aes(x=no,y=X20))+geom_line(col="#333333")+
  geom_hline(yintercept=qt3$Y20,lwd=1,linetype="dashed")+
  geom_ribbon(data=subset(qt3,X20<=qt3$X20[10]),
              aes(ymin=X20,ymax=Y20),fill="grey")+
  geom_line(aes(x=no,y=X40),col="#CC0000")+
  geom_hline(yintercept=qt3$Y40,lwd=1,col="#990000",linetype="dashed")+
  geom_ribbon(data=subset(qt3,X40<=qt3$X40[67]),
              aes(ymin=X40,ymax=Y40),fill="pink")+
  geom_line(aes(x=no,y=X120),col="#00FF00")+
  geom_hline(yintercept=qt3$Y120,lwd=1,col="#006600",linetype="dashed")+
  geom_ribbon(data=subset(qt3,X120<=qt3$X120[72]),
              aes(ymin=X120,ymax=Y120),fill="yellow")+
  geom_line(aes(x=no,y=X140),col="#330099")+
  geom_hline(yintercept=qt3$Y140,lwd=1,col="#003366",linetype="dashed")+
  geom_ribbon(data=subset(qt3,X140>=0 & X140<=qt3$X140[87]),
              aes(ymin=X140,ymax=Y140),fill="purple")+
  geom_line(aes(x=no,y=X160),col="#0099CC")+
  geom_hline(yintercept=qt3$Y160,lwd=1,col="#3399FF",linetype="dashed")+
  geom_ribbon(data=subset(qt3,0<=X160 & X160<=qt3$X160[96]),
              aes(ymin=X160,ymax=Y160),fill="lightblue")+
  xlab("Number of values of quantile")+
  ylab(TeX(r'(Sorted ${\X}^{T}$$\beta$ by observation,response values)'))

g4<-g3+annotate(x=+Inf,y=yg[1,20],label="Y20",hjust=2,vjust=-0.05,geom="label",
                size=12/.pt)+
  annotate(x=+Inf,y=yg[1,40],label="Y40",hjust=2,vjust=-0.05,geom="label",
           size=12/.pt)+
  annotate(x=+Inf,y=yg[1,120],label="Y120",hjust=4,geom="label",
           size=12/.pt)+
  annotate(x=+Inf,y=yg[1,140],label="Y140",hjust=1,geom="label",
           size=12/.pt)+
  annotate(x=125,y=yg[1,160],label="Y160",geom="label",
           size=12/.pt)
#####################################################################################
#####################################################################################
#####################################################################################
plot_grid(g1,g4)


tQ.test   <- cbind(1,Xnew)%*%beta.qr
Q.test    <- t(apply(tQ.test,1,FUN=sort))
####################################################################################
####################################################################################
qtest<-data.frame(Q.test)
qtestt<-data.frame(t(qtest))
qtestt$no<-1:200

ft<-ggplot2::ggplot(qtestt,mapping=aes(no,t.qtest.))+geom_point()+
  xlab("Total number of unique values of quantile")+
  ylab(TeX(r'(Sorted ${\X}^{T}$$\beta$ for test data)'))

ft

####################################################################################
####################################################################################


qtestt2<-data.frame(qtestt[,-2])
names(qtestt2)<-c("Qtest")
yt<-matrix(t(matrix(Ynew,length(Ynew),length(taus))),nrow=200,ncol=1)
yt<-data.frame(yt) 


ft2<-ggplot2::ggplot(qtestt,mapping=aes(no,t.qtest.))+geom_point(col="#333333")+
  geom_hline(yintercept=yt[1,1], lwd=1,linetype="dashed")+
  xlab("Total number of unique values of quantile")+
  ylab(TeX(r'(Sorted ${\X}^{T}$$\beta$ for test data)'))

ft2


####################################################################################
####################################################################################
####################################################################################
u.test    <- rowMeans((Q.test  <= matrix(Ynew,length(Ynew),length(taus))))
cs.test   <- abs(u.test-0.5)


alpha.sig<-0.1

k<- ceiling((1-alpha.sig)*(1+length(Y1a)))
threshold <- sort(cs)[k]

####################################################################################
####################################################################################
qtestt$t.qtest.t<=yt
yt$no<-1:200
library(dplyr)
yt<-yt%>%select(no,yt)

dfjoin2 <- left_join(qtestt,yt, by = "no")

C2<-colSums(qtestt2<=yt[,2])  ##gives maximum index upto which 
###sorted coefficients are less than corresponding y value 
C2
#15

qtestt$y<-yt[,2]

gt<-ggplot2::ggplot(qtestt,mapping=aes(x=no,y=t.qtest.))+geom_line(col="#333333")+
  geom_hline(yintercept=yt[1,2],lwd=1,linetype="dashed")+
  geom_ribbon(data=subset(qtestt,t.qtest.<=qtestt$t.qtest.[15]),
              aes(ymin=t.qtest.,ymax=y),fill="grey")+
  xlab("Total number of unique values of quantile")+
  ylab(TeX(r'(Sorted ${\X}_{test}^{T}$$\beta$)'))


gt2<-gt+annotate(x=+Inf,y=yt[1,2],label="Ynew",hjust=2,vjust=-0.05,geom="label",
                 size=12/.pt)
gt2

#################################################################################
#################################################################################

ci.grid   <- abs(taus - 0.5)
cov.qr    <- (cs.test <= threshold)
D<-data.frame(ci.grid)
D$t<-threshold
D$no<-1:200
D$c<-cs.test
names(D)<-c("grid","threshold","no","test")
##################################################################################
##################################################################################
ct<-ggplot2::ggplot(D,aes(x=no,y=ci.grid))+geom_point()+
  geom_hline(yintercept = D$threshold[1],col="darkred",lwd=1)+
  geom_vline(xintercept = 6,col="magenta",lwd=1)+
  geom_vline(xintercept = 195,col="magenta",lwd=1)+
  
  geom_hline(yintercept = D$test[1],col="green",lwd=1)+
  xlab("Total number of values of quantiles")+
  ylab("Set of grid values and threshold")
ct2<-ct+annotate(x=+Inf,y=threshold,label="Threshold",hjust=2,vjust=-0.05,geom="label",
                 size=12/.pt)+
  annotate(x=135,y=D$grid[120],label="Grid values",geom="label",size=12/.pt,
           hjust=1.5)+
  annotate(x=125,y=D$test[1],label="Test Point",hjust=2,vjust=-0.05,geom="label",
           size=12/.pt)+
  annotate(x=175,y=0.1,label="Upper Bound:Grid<=Threshold",geom="label",
           size=9/.pt)+
  annotate(x=25,y=0.1,label="Lower Bound:Grid<=Threshold",geom="label",
           size=9/.pt)
  
plot_grid(ft2,gt2)
#####################################################################################
####################################################################################

lb <- ub <- rep(NA,length(Ynew))
for (i in 1:length(Ynew)){
  ci     <- Q.test[i,(ci.grid <= threshold)]
  ub[i]  <- max(ci)
  lb[i]  <- min(ci)
}

leng.qr <- ub-lb
leng.qr[which(leng.qr==-Inf)] <- NA

###################################################################################
####################################################################################


lb <- ub <- rep(NA,length(Ynew))
for (i in 1:length(Ynew)){
  ci2     <- qtest[i,(ci.grid <= threshold)]
  ub[i]  <- max(ci2)
  lb[i]  <- min(ci2)
}
###X6-X195
##########################################################################
######################################################################

a0<- 6
a1<-195
highlight<- data.frame(start=a0, end=a1, group=seq_along(a0))

ft3<-ggplot2::ggplot(qtestt,mapping=aes(no,t.qtest.))+geom_point(col="#333333")+
  geom_rect(data=highlight,inherit.aes=FALSE, 
            aes(xmin=a0, xmax=a1, ymin=min(qtestt$t.qtest.),
                ymax=max(qtestt$t.qtest.), group=group), color="transparent", fill="magenta", alpha=0.25)+
  xlab("Number of values of quantiles")+
  ylab(TeX(r'(Sorted ${\X}_{test}^{T}$$\beta$)'))+
  geom_point(aes(x=6,y=lb),col="maroon",size=5)+
  geom_point(aes(x=195,y=ub),col="maroon",size=5)+
  annotate(x=10,y=lb-0.25,label="Lower bound",geom="label",size=12/.pt)+
  annotate(x=193,y=ub+0.25,label="Upper bound",geom="label",size=12/.pt,
           hjust=0.65)

ft3

plot_grid(ct2,ft3,nrow=2)


###########################################################################
###########################################################################
a4<-a+geom_hline(yintercept = lb,col="darkgreen",size=1.5)+
  geom_hline(yintercept=ub,col="darkgreen",size=1.5)+
  annotate(x=0.67,y=40,label="Xtest:0.67",geom="label",size=11/.pt)+
  annotate(x=-5,y=2.7,label="Ytest:2.7",vjust=-0.1,geom="label",size=11/.pt)+
  annotate(x=10,y=ub,label="Upper bound:6.1836",geom="label",vjust=0.1,size=11/.pt)+
  annotate(x=10,y=lb,label="Lower bound:2.2348",geom="label",size=11/.pt,
           vjust=1)
a4

#############################################################################