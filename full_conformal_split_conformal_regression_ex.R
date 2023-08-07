#################################################################################################
###SIMPLE EXAMPLE: FULL CONFORMAL PREDICTION METHOD##############################################
##################################################################################################
install.packages(ggplot2)
library(ggplot2)
install.packages("ggpubr")
library(ggpubr)
install.packages("cowplot")
library(cowplot)
##Setting up random data 
set.seed(1)
n<-1000
X<-rnorm(n,1,3)
U<-rnorm(n,0,1)
Y<-2+3*X+U
datXY<-cbind.data.frame(X,Y)
###data point to predict
Xnew<-0.67

###set of candidate y-values for test value,X
candidate_y<-seq(from=min(Y),to=max(Y),length=200)
alpha<-0.1

####Full conformal prediction#####################################
###################################################################
listc<-NULL
for (i in 1:length(candidate_y)){
  point<-c(Xnew,candidate_y[i])
  newDat<-rbind.data.frame(datXY,point)
  names(newDat)<-c("X","Y")
  mod<-lm(Y~X,data=newDat)   ###Step 1:linear regression
  resid<-abs(mod$residuals)    ###Step 2: Absolute value of residuals
  threshold<-resid[length(resid)]  
  pi0<-ifelse(resid<=threshold,1,0) ##Step 3: proportion of residuals<=residuals for Xnew
  piy<-sum(pi0)/(nrow(newDat))
  threshold2<-ceiling((1-alpha)*(nrow(newDat))) ##Step 4: threshold of (n+1)(1-alpha)
  piyp<-piy*(nrow(newDat))
  if (piyp<=threshold2)
  {listc<-rbind(listc,candidate_y[i])}        ###Step 5: choosing candidate y values<= threshold
  
  else {next}
  
}

v<-range(listc)

############################################################################
############################################################################
##plotting results##########################################################

###Demonstrating-Step 1: linear regression-using candidate value y1
point<-c(Xnew,candidate_y[1])
newDat<-rbind.data.frame(datXY,point)
names(newDat)<-c("X","Y")
mod<-lm(Y~X,data=newDat)
a<-ggplot(data=newDat,aes(x=X,y=Y))+geom_point()+geom_smooth(method="lm",se=FALSE)+
     geom_vline(xintercept=0.67,col="red",linetype="dashed",size=1)

a1<-a+annotate(x=0.67,y=+Inf,label="Xtest:0.67",vjust=2,geom="label")

##step2: Absolute value of residuals
z<-1:1001
resid<-abs(mod$residuals)
f<-cbind.data.frame(z,resid)
names(f)<-c("z","resid")
lastval<-f[1001,2]

b<-ggplot(data=f,aes(x=z,y=resid))+geom_point()+
  geom_vline(xintercept=1001,col="red",linetype="dashed",size=1)+
  xlab("Index")+ylab("Absolute Value of Residuals")


b1<-b+annotate(x=101,y=+Inf,label="R(Xtest,y)",vjust=2,hjust=-6.5,geom="label")

###Step-3
listp<-NULL
for (i in 1:length(candidate_y)){
  point<-c(Xnew,candidate_y[i])
  newDat<-rbind.data.frame(datXY,point)
  names(newDat)<-c("X","Y")
  mod<-lm(Y~X,data=newDat)   ###Step 1:linear regression
  resid<-abs(mod$residuals)    ###Step 2: Absolute value of residuals
  threshold<-resid[length(resid)]  
  pi0<-ifelse(resid<=threshold,1,0) ##Step 3: proportion of residuals<=residuals for Xnew
  piy<-sum(pi0)/(nrow(newDat))
  threshold2<-ceiling((1-alpha)*(nrow(newDat))) ##Step 4: threshold of (n+1)(1-alpha)
  piyp<-piy*(nrow(newDat))
  listp<-rbind(listp,piyp)
  
}

th<-ceiling((1-alpha)*(nrow(newDat)))
t<-rep(th,200)

g<-cbind.data.frame(t,listp,1:200)
names(g)<-c("Threshold","pi","Index")

thc<-g$Threshold
c<-ggplot(data=g,aes(x=Index,y=pi))+geom_point()+
  geom_hline(yintercept=thc,col="red",linetype="dashed",size=1)+
  ylab("pi*size of Augmented Data")

c1<-c+annotate(x=+Inf,y=thc,label="(1-alpha)(1+Data Size)=901",hjust=1,geom="label")

###############################################################################
####Prediction intervals#######################################################
j<-range(listc)

d<-ggplot(data=newDat,aes(x=X,y=Y))+geom_point()+geom_smooth(method="lm",se=FALSE)+
  geom_vline(xintercept=0.67,col="red",linetype="dashed",size=1)+
  geom_hline(yintercept=j[1],col="green",linetype="dashed",size=1)+
  geom_hline(yintercept=j[2],col="green",linetype="dashed",size=1)
  

d1<-d+annotate(x=0.67,y=+Inf,label="Xtest:0.67",vjust=2,geom="label")+
  annotate(x=+Inf,y=j[1],label="Lower tail:2.36",hjust=1.8,geom="label")+
  annotate(x=+Inf,y=j[2],label="Upper tail:5.58",hjust=1.5,geom="label")
  

plot_grid(
  a1,b1,c1,d1,
  align = "h", axis = "tb",
  nrow = 2, rel_widths = c(2, 2)
)


############################################################################################
###########SPLIT CONFORMAL PREDICTION-SIMPLE EXAMPLE########################################
############################################################################################
set.seed(2)
alpha<-0.1
ind<-sample(nrow(datXY),floor(nrow(datXY)/2),replace=FALSE)##split into training and calibration sets
Ytrain<-data.frame(datXY[ind,2])
Xtrain<-data.frame(datXY[ind,1])
train<-cbind.data.frame(Xtrain,Ytrain)
names(train)<-c("X","Y")
mod<-lm(Y~X,data=train)
Ycalib<-datXY[-ind,2]
Xcalib<-datXY[-ind,1]
calib<-cbind.data.frame(Xcalib,Ycalib)
names(calib)<-c("X","Y")
resid_calib<-abs(calib[,2]-predict(mod,calib)) ##Conformity scores on calibration set
yhat<-predict(mod,newdata=data.frame(X=Xnew))
emp<-ceiling((1+nrow(calib))*(1-alpha)) ##threshold or k
##find empth empirical quantile
resid_calib_sorted<-sort(resid_calib) ##sorted conformity scores on calibration set
emp_qtile<-resid_calib[emp] ##finding the kth empirical quantile 
pred_interval<-c(yhat-emp_qtile,yhat+emp_qtile)  ##prediction interval

############################################################################################
############################################################################################
##Step 1: split into training and calibration data 
##Step 2: regression on training data 
a<-ggplot(datXY,aes(x=X,y=Y))+geom_point()+
  geom_vline(xintercept=0.67,col="red",linetype="dashed",size=1)
b<-ggplot(train,aes(x=X,y=Y))+geom_point(col="blue")+geom_smooth(method="lm",se=FALSE)
c<-ggplot(calib,aes(x=X,y=Y))+geom_point(col="orange")

a1<-a+annotate(x=+Inf,y=+Inf,label="Entire Data",vjust=3,hjust=-7,geom="label")+
  annotate(x=0.67,y=+Inf,label="Xtest:0.67",vjust=2,geom="label")

b1<-b+annotate(x=+Inf,y=+Inf,label="Training Data",vjust=5,hjust=1,geom="label")
c1<-c+annotate(x=+Inf,y=+Inf,label="Calibration Data",vjust=5,hjust=1,geom="label")

plot_grid(
  a1, b1,c1,
  align = "h", axis = "tb",
  nrow = 1, rel_widths = c(1, 1)
)

###step 3: conformity scores
###ranking of conformity scores to 
###obtain kth smallest empirical quantile 
###where k=(1-alpha)*(1+size(calibration set))
f<-cbind.data.frame(resid_calib,1:500,yhat)
g<-cbind.data.frame(resid_calib_sorted,1:500,yhat)
names(f)<-c("resid","conseq","yhat")
names(g)<-names(f)
d<-ggplot(f,aes(x=conseq,y=resid))+geom_point(col="orange")+
  xlab("Index from 1 to 500")+ylab("Conformity Scores(unsorted)")
e<-ggplot(g,aes(x=conseq,y=resid))+geom_point(col="purple")+
  geom_vline(xintercept=emp,col="red",linetype="dashed",size=1)+
  geom_hline(yintercept=emp_qtile,col="red",linetype="dashed",size=1)+
  xlab("Index from 1 to 500")+ylab("Conformity scores")

e1<-e+annotate(x=emp,y=+Inf,label="kth position",vjust=2,geom="label")+
  annotate(x=+Inf,y=emp_qtile,label="kth empirical quantile",hjust=1.8,geom="label")+
  annotate(x=200,y=0.7,label="Sorted conformity scores(calibration data)",geom="label")

###Step 4: obtaining (1-alpha) prediction intervals using
##predicted y and kth empirical quantilea 
k<-ggplot(datXY,aes(x=X,y=Y))+geom_point()+geom_smooth(method="lm",se=FALSE)+
  geom_vline(xintercept=0.67,col="red",linetype="dashed",size=1)+
  geom_hline(yintercept=yhat,col="red",linetype="dashed",size=1)+
  geom_hline(yintercept=3.885431,col="green",linetype="dashed",size=1)+
  geom_hline(yintercept=4.148118 ,col="yellow",linetype="dashed",size=1)

k

############make dataset small in R to 
#######################zoom in on upper and lower tails of prediction interval
library(dplyr)
library(tidyverse)

dat2<-datXY%>%filter(X>=-2.5 & X<=2.5)%>%filter(Y>0& Y<6)
k<-ggplot(dat2,aes(x=X,y=Y))+geom_point()+
  geom_vline(xintercept=0.67,col="red",linetype="dashed",size=1)+
  geom_hline(yintercept=yhat,col="red",linetype="dashed",size=1)+
  geom_hline(yintercept=pred_interval[1],col="green",linetype="dashed",size=1)+
  geom_hline(yintercept=pred_interval[2],col="green",linetype="dashed",size=1)

#a+annotate(x=0.67,y=+Inf,label="Xtest:0.67",vjust=2,geom="label")
k1<-k+annotate(x=0.67,y=+Inf,label="Xtest:0.67",vjust=5,hjust=1,geom="label")+
  annotate(x=+Inf,y=yhat,label="Predicted y:4.017",hjust=2,geom="label")+
  annotate(x=+Inf,y=pred_interval[1],label="Lower tail:3.89",hjust=5,vjust=1,geom="label")+
  annotate(x=+Inf,y=pred_interval[2],label="Upper tail:4.15",hjust=5,vjust=0.05,geom="label")+
  annotate(x=0,y=1,label="Upper tail:Predicted y+kth empirical quantile",geom="label")+
  annotate(x=0,y=0.75,label="Lower tail:Predicted y-kth empirical quantile",geom="label")



