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
########################################################
###loading required packages




install.packages("VGAM")
library(VGAM)
library(openxlsx)
library(lmtest)
library(quantreg)
library(hdm)
library(matrixcalc)
library(cowplot)
devtools::install_github("thomasp85/patchwork")
library(patchwork) 
library(handyplots)
library(tidyverse)
library(ggplot2)

boston<-read.xlsx("boston.xlsx")
head(boston)
names(boston)<-c("crim","zn","indus","chas","nox",
                 "rm","age","dis","rad","tax",
                 "ptratio","b","lstat","medv")


#CRIM     per capita crime rate by town
#ZN       proportion of residential land zoned for lots over 25,000 sq.ft.
#INDUS    proportion of non-retail business acres per town
#CHAS     Charles River dummy variable (= 1 if tract bounds river; 0 otherwise)
#NOX      nitric oxides concentration (parts per 10 million)
#RM       average number of rooms per dwelling
#AGE      proportion of owner-occupied units built prior to 1940
#DIS      weighted distances to five Boston employment centres
#RAD      index of accessibility to radial highways
#TAX      full-value property-tax rate per $10,000
#PTRATIO  pupil-teacher ratio by town
#B        1000(Bk - 0.63)^2 where Bk is the proportion of blacks by town
#LSTAT    % lower status of the population
#MEDV     Median value of owner-occupied homes in $1000's

x<- model.matrix(~-1 + (crim+zn+indus+chas+nox+
                          rm+age+dis+rad+tax+
                          ptratio+b+lstat),data=boston)
y <-boston$medv






###eda-linear regression and tobit
###Uncorrected versions:
lfit<-lm(medv~crim+zn+indus+chas +nox+rm+age+dis+rad+tax+ptratio+b+lstat,data=boston)
l<-summary(lfit)
lmtest::bptest(lfit)

lfit2<-lm(log(medv)~crim+zn+indus+chas +nox+rm+age+dis+rad+tax+ptratio+b+lstat,data=boston)
l2<-summary(lfit2)
lmtest::bptest(lfit2)



resplot(lfit,highlight.outliers=TRUE)

########################################################################################
###Tobit regression####################################################################
########################################################################################
bt <- vglm(medv~crim+zn+indus+chas +nox+rm+age+dis+rad+tax+ptratio+b+lstat , tobit(Upper = 50), data = boston)

fit <- fitted(bt)
residuals <- resid(bt, type = "response")
qqnorm(residuals)
qqline(residuals,col="red")




library(quantreg)


###function for prediction interval
yogen<-function(t,x,y){
  T.ho <- floor(length(y)*0.2)
  cov.mat <- leng.mat <- X.test.mat <- NULL
  Y.test.mat<-Y0.mat<-Y1.mat<-NULL
  X0.mat<-X1.mat<-NULL
  start<-Sys.time()
  for (i in 1:t){
    start<-Sys.time()
    # Define training and holdout samples (while loop is to avoid singular design)
    sing <- TRUE
    while (sing==TRUE){
      set.seed(i)
      
      ind <- sample(length(y),length(y),replace=FALSE)
      Y   <- y[ind]
      X   <- x[ind,]
      
      ind.test  <- (length(Y)-T.ho+1):length(Y)
      Y.test    <- Y[ind.test]
      X.test    <- cbind(X[ind.test,])
      
      obj.uneven <- aux.uneven(Y[-ind.test],cbind(X[-ind.test,]))
      
      X0 <- cbind(obj.uneven$X0)
      Y0 <- obj.uneven$Y0
      X1 <- cbind(obj.uneven$X1)
      Y1 <- obj.uneven$Y1
      
      sing  <- is.singular.matrix(t(X0)%*%X0)
    }
    
    # Choosing grids for QR and DR (we found that DR often works better with a coarser grid)
    taus    <-  seq(0.001,0.999,length=200)  ##change taus to avoid infinite values 
    #ys      <- quantile(unique(c(Y0,Y1)),seq(0.001,0.999,length=200))
    
    # Applying the difference conformal prediction methods
    res.qr    <- dcp.qr(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)
    res.opt   <- dcp.opt(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)
    res.cqr   <- cqr(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
    res.reg   <- cp.reg(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
    # res.loc   <- cp.loc(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
    
    # Results
    X.test.mat <- rbind(X.test.mat,X.test)
    Y.test.mat <- rbind(Y.test.mat,Y.test)
    X0.mat <- rbind(X0.mat,X0)
    Y0.mat <- rbind(Y0.mat,Y0)
    X1.mat <- rbind(X1.mat,X1)
    Y1.mat <- rbind(Y1.mat,Y1)
    
    
    cov.mat.temp  <- cbind(res.qr$cov.qr,res.opt$cov.opt,res.cqr$cov.o,res.cqr$cov.m,res.cqr$cov.r,res.reg$cov.reg)
    leng.mat.temp <- cbind(res.qr$leng.qr,res.opt$leng.opt,res.cqr$leng.o,res.cqr$leng.m,res.cqr$leng.r,res.reg$leng.reg)
    
    cov.mat   <- rbind(cov.mat,cov.mat.temp)
    leng.mat  <- rbind(leng.mat,leng.mat.temp)
    
  }
  j<-Sys.time()-start
  
  l<-list(cov.mat,leng.mat,
          X.test.mat,Y.test.mat,X0.mat,Y0.mat,
          X1.mat,Y1.mat,j)
  return(l)
}


b<-yogen(100,x,y)
#####################################################################################
###Saving results
#####################################################################################
library(openxlsx)
write.xlsx(b,"bost100.xlsx")

COV<-data.frame(b[1])
LENG<-data.frame(b[2])

ac<-colMeans(COV,na.rm=TRUE)
al<-colMeans(LENG,na.rm=TRUE)


##############################################sd

ps<-function(t){
  
  COV<-data.frame(t[1])
  LENG<-data.frame(t[2])
  
  ac<-colMeans(COV,na.rm=TRUE)
  al<-colMeans(LENG,na.rm=TRUE)
  pred.cov.qr       <- predict(glm(COV[,1]~as.matrix(data.frame(t[3])),family=binomial(link="logit")),type="response")
  pred.cov.opt      <- predict(glm(COV[,2]~as.matrix(data.frame(t[3])),family=binomial(link="logit")),type="response")
  pred.cov.cqr.o    <- predict(glm(COV[,3]~as.matrix(data.frame(t[3])),family=binomial(link="logit")),type="response")
  pred.cov.cqr.m    <- predict(glm(COV[,4]~as.matrix(data.frame(t[3])),family=binomial(link="logit")),type="response")
  pred.cov.cqr.r    <- predict(glm(COV[,5]~as.matrix(data.frame(t[3])),family=binomial(link="logit")),type="response")
  pred.cov.reg      <- predict(glm(COV[,6]~as.matrix(data.frame(t[3])),family=binomial(link="logit")),type="response")
  
  pred.cov <- cbind(pred.cov.qr,pred.cov.opt,pred.cov.cqr.o,pred.cov.cqr.m,pred.cov.cqr.r,pred.cov.reg)
  
  res.cov.cond <- apply(pred.cov,2,sd)
  
  k<-list(ac,al,res.cov.cond)
  
  return(k)
}

ps(b)

#################################################################################################
#################################################################################################
dat<-cbind.data.frame(ps(b)[1],ps(b)[2],ps(b)[3])
names(dat)<-c("Avg.Coverage","Avg.Length","SD")
xtable(dat,digits=4)


################################################################################################################
####################################################################################################################
########################################################################################################################

###Histogram

COV<-data.frame(b[1])

pred.cov.qr       <- predict(glm(COV[,1]~as.matrix(data.frame(b[3])),family=binomial(link="logit")),type="response")
pred.cov.opt      <- predict(glm(COV[,2]~as.matrix(data.frame(b[3])),family=binomial(link="logit")),type="response")
pred.cov.cqr.o    <- predict(glm(COV[,3]~as.matrix(data.frame(b[3])),family=binomial(link="logit")),type="response")
pred.cov.cqr.m    <- predict(glm(COV[,4]~as.matrix(data.frame(b[3])),family=binomial(link="logit")),type="response")
pred.cov.cqr.r    <- predict(glm(COV[,5]~as.matrix(data.frame(b[3])),family=binomial(link="logit")),type="response")
pred.cov.reg      <- predict(glm(COV[,6]~as.matrix(data.frame(b[3])),family=binomial(link="logit")),type="response")


H<-cbind.data.frame(pred.cov.qr,pred.cov.opt,pred.cov.cqr.o,pred.cov.cqr.m,pred.cov.cqr.r,pred.cov.reg)
H$breaks<-seq(0,1,length.out=dim(H)[1])
names(H)<-c("DCP","DCPopt","CQR","CQRm","CQRr","CPreg","breaks")


f1<-ggplot2::ggplot(H,aes(x=DCP))+geom_histogram(color="darkblue", fill="lightblue")+
  geom_vline(xintercept = 0.90,col="red",lwd=1)+
  annotate(x=0.90,y=2500,label="0.90",geom="label")+
  ylab("Count")

f2<-ggplot2::ggplot(H,aes(x=DCPopt))+geom_histogram(color="darkblue", fill="lightblue")+
  geom_vline(xintercept = 0.90,col="red",lwd=1)+
  annotate(x=0.90,y=2500,label="0.90",geom="label")+
  ylab("Count")+xlab("DCP-opt")

f3<-ggplot2::ggplot(H,aes(x=CQR))+geom_histogram(color="darkblue", fill="lightblue")+
  geom_vline(xintercept = 0.90,col="red",lwd=1)+
  annotate(x=0.90,y=2500,label="0.90",geom="label")+
  ylab("Count")+xlab("CQR")

f4<-ggplot2::ggplot(H,aes(x=CQRm))+geom_histogram(color="darkblue", fill="lightblue")+
  geom_vline(xintercept = 0.90,col="red",lwd=1)+
  annotate(x=0.90,y=2500,label="0.90",geom="label")+
  ylab("Count")+xlab("CQR-m")

f5<-ggplot2::ggplot(H,aes(x=CQRr))+geom_histogram(color="darkblue", fill="lightblue")+
  geom_vline(xintercept = 0.90,col="red",lwd=1)+
  annotate(x=0.90,y=2500,label="0.90",geom="label")+
  ylab("Count")+xlab("CQR-r")

f6<-ggplot2::ggplot(H,aes(x=CPreg))+geom_histogram(color="darkblue", fill="lightblue")+
  geom_vline(xintercept = 0.90,col="red",lwd=1)+
  annotate(x=0.90,y=2500,label="0.90",geom="label")+
  ylab("Count")+xlab("CP-reg")


plot_grid(f1,f2,f3,f4,f5,f6) 

###################################################################################
###################################################################################
###Alternative specification#######################################################

boston$mv<-log(boston$medv)
boston$n2<-(boston$nox)^2
boston$r2<-(boston$rm)^2
boston$ld<-log(boston$dis)
boston$lr<-log(boston$rad)
boston$ls<-log(boston$lstat)

x2<-model.matrix(~-1 + (crim+chas+nox+n2+
                          r2+age+ld+lr+tax+
                          ptratio+b+ls+rm+indus+zn),data=boston)
y2<-boston$mv


b2<-yogen(100,x2,y2)

#####################################################################################
#####################################################################################


write.xlsx(b2,"bost_2_100.xlsx")

COV2<-data.frame(b2[1])
LENG2<-data.frame(b2[2])

ac2<-colMeans(COV2,na.rm=TRUE)
al2<-colMeans(LENG2,na.rm=TRUE)


##############################################sd

ps2<-function(t){
  
  COV2<-data.frame(t[1])
  LENG2<-data.frame(t[2])
  
  ac<-colMeans(COV2,na.rm=TRUE)
  al<-colMeans(LENG2,na.rm=TRUE)
  pred.cov.qr       <- predict(glm(COV[,1]~as.matrix(data.frame(t[3])),family=binomial(link="logit")),type="response")
  pred.cov.opt      <- predict(glm(COV[,2]~as.matrix(data.frame(t[3])),family=binomial(link="logit")),type="response")
  pred.cov.cqr.o    <- predict(glm(COV[,3]~as.matrix(data.frame(t[3])),family=binomial(link="logit")),type="response")
  pred.cov.cqr.m    <- predict(glm(COV[,4]~as.matrix(data.frame(t[3])),family=binomial(link="logit")),type="response")
  pred.cov.cqr.r    <- predict(glm(COV[,5]~as.matrix(data.frame(t[3])),family=binomial(link="logit")),type="response")
  pred.cov.reg      <- predict(glm(COV[,6]~as.matrix(data.frame(t[3])),family=binomial(link="logit")),type="response")
  
  pred.cov <- cbind(pred.cov.qr,pred.cov.opt,pred.cov.cqr.o,pred.cov.cqr.m,pred.cov.cqr.r,pred.cov.reg)
  
  res.cov.cond <- apply(pred.cov,2,sd)
  
  k<-list(ac,al,res.cov.cond)
  
  return(k)
}

ps2(b2)

#################################################################################################
#################################################################################################
dat<-cbind.data.frame(ps2(b2)[1],ps2(b2)[2],ps2(b2)[3])
names(dat)<-c("Avg.Coverage","Avg.Length","SD")
xtable(dat,digits=4)


################################################################################################################
####################################################################################################################
########################################################################################################################

###Histogram

COV<-data.frame(b2[1])

pred.cov.qr       <- predict(glm(COV[,1]~as.matrix(data.frame(b2[3])),family=binomial(link="logit")),type="response")
pred.cov.opt      <- predict(glm(COV[,2]~as.matrix(data.frame(b2[3])),family=binomial(link="logit")),type="response")
pred.cov.cqr.o    <- predict(glm(COV[,3]~as.matrix(data.frame(b2[3])),family=binomial(link="logit")),type="response")
pred.cov.cqr.m    <- predict(glm(COV[,4]~as.matrix(data.frame(b2[3])),family=binomial(link="logit")),type="response")
pred.cov.cqr.r    <- predict(glm(COV[,5]~as.matrix(data.frame(b2[3])),family=binomial(link="logit")),type="response")
pred.cov.reg      <- predict(glm(COV[,6]~as.matrix(data.frame(b2[3])),family=binomial(link="logit")),type="response")


H<-cbind.data.frame(pred.cov.qr,pred.cov.opt,pred.cov.cqr.o,pred.cov.cqr.m,pred.cov.cqr.r,pred.cov.reg)
H$breaks<-seq(0,1,length.out=dim(H)[1])
names(H)<-c("DCP","DCPopt","CQR","CQRm","CQRr","CPreg","breaks")


f1<-ggplot2::ggplot(H,aes(x=DCP))+geom_histogram(color="darkblue", fill="lightblue")+
  geom_vline(xintercept = 0.90,col="red",lwd=1)+
  annotate(x=0.90,y=4000,label="0.90",geom="label")+
  ylab("Count")

f2<-ggplot2::ggplot(H,aes(x=DCPopt))+geom_histogram(color="darkblue", fill="lightblue")+
  geom_vline(xintercept = 0.90,col="red",lwd=1)+
  annotate(x=0.90,y=4000,label="0.90",geom="label")+
  ylab("Count")+xlab("DCP-opt")

f3<-ggplot2::ggplot(H,aes(x=CQR))+geom_histogram(color="darkblue", fill="lightblue")+
  geom_vline(xintercept = 0.90,col="red",lwd=1)+
  annotate(x=0.90,y=4000,label="0.90",geom="label")+
  ylab("Count")+xlab("CQR")

f4<-ggplot2::ggplot(H,aes(x=CQRm))+geom_histogram(color="darkblue", fill="lightblue")+
  geom_vline(xintercept = 0.90,col="red",lwd=1)+
  annotate(x=0.90,y=4000,label="0.90",geom="label")+
  ylab("Count")+xlab("CQR-m")

f5<-ggplot2::ggplot(H,aes(x=CQRr))+geom_histogram(color="darkblue", fill="lightblue")+
  geom_vline(xintercept = 0.90,col="red",lwd=1)+
  annotate(x=0.90,y=4000,label="0.90",geom="label")+
  ylab("Count")+xlab("CQR-r")

f6<-ggplot2::ggplot(H,aes(x=CPreg))+geom_histogram(color="darkblue", fill="lightblue")+
  geom_vline(xintercept = 0.90,col="red",lwd=1)+
  annotate(x=0.90,y=4000,label="0.90",geom="label")+
  ylab("Count")+xlab("CP-reg")

 
plot_grid(f1,f2,f3,f4,f5,f6) 































