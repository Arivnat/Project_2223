###SIMULATION-CASE 1

library(xtable)
library(hdm)
library(matrixcalc)
library(hdm)
library(quantreg)
library(openxlsx)
library(readxl)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyverse)

source("functions_final.R") ###loads all relevant algorithms and functions

###All functions under functions_final.R as well as 
##functions used after data generation under algsine function are 
##from:
#Chernozhukov, V.,Wuthrich,K. &  Zhu,Y.(2021).Distributional Conformal Prediction. PNAS,118 (48).
##https://doi.org/10.1073/pnas.2107794118

##The relevant code taken from the authors is in their github account:
#https://github.com/kwuthrich/Replication_DCP 

##########################################################################
###########################################################################

###General function for data generation and then obtaining length,coverage and 
##standard deviation of prediction intervals

algsine<-function(n){
  alpha.sig<-0.1
  T.ho <- floor(n*0.2) ##20% holdout sample for test set 
  

  cov.mat <- leng.mat <- X.test.mat <- NULL
  cov.mat2 <- leng.mat2 <- X.test.mat2 <- NULL
  cov.mat3 <- leng.mat3 <- X.test.mat3 <- NULL
  cov.mat4 <- leng.mat4 <- X.test.mat4 <- NULL
  X0.mat<-Y0.mat<-X1.mat<-Y1.mat<-NULL
  X02.mat<-Y02.mat<-X12.mat<-Y12.mat<-NULL
  X03.mat<-Y03.mat<-X13.mat<-Y13.mat<-NULL
  X04.mat<-Y04.mat<-X14.mat<-Y14.mat<-NULL
  Y.test.mat<-Y.test.mat2<-Y.test.mat3<-Y.test.mat4<-NULL
  
 ####data generation for 4 specifications 
  for (i in 1:1000){
    set.seed(i)
    x1<-runif(n,-2,2)
    x2<-rnorm(n,0,1)
    eps<-rnorm(n,0,1)
    y1<-x1+2*exp(-16*(x1^2))+0.4*eps
    y2<-sin(2*x1)+2*(exp(-16*(x1^2)))+0.3*eps
    y3<-0.3*exp(-4*(x1+1)^2)+0.7*exp(-16*(x1-1)^2)+0.1*eps
    y4<-0.4*x2+0.15*eps
   
###checking of matrix is singular or not-necessary for 
###quantile regression by quantreg package in R
    sing <- TRUE
    sing2<-TRUE
    sing3<-TRUE
    sing4<-TRUE
    while (sing==TRUE){
      while (sing2==TRUE){
        while (sing3==TRUE){
          while (sing4==TRUE){
        
        ind <- sample(n,n,replace=FALSE)
        Y   <- y1[ind]
        Y2   <- y2[ind]
        Y3   <- y3[ind]
        Y4   <- y4[ind]
        Xa   <- x1[ind]
        Xb   <- x2[ind]
        Xa<-as.matrix(Xa)
        Xb<-as.matrix(Xb)
        Y<-as.matrix(Y)
        Y2<-as.matrix(Y2)
        Y3<-as.matrix(Y3)
        Y4<-as.matrix(Y4)

        ####creation of training,calibration and test sets         
        ind.test  <- (length(Y)-T.ho+1):length(Y)
        Y.test    <- Y[ind.test]
        X.test    <- cbind(Xa[ind.test,])
        
        ind.test2  <- (length(Y2)-T.ho+1):length(Y2)
        Y.test2    <- Y2[ind.test2]
        X.test2    <- cbind(Xa[ind.test2,])
        
        ind.test3  <- (length(Y3)-T.ho+1):length(Y3)
        Y.test3    <- Y3[ind.test3]
        X.test3    <- cbind(Xa[ind.test3,])
        
        ind.test4  <- (length(Y4)-T.ho+1):length(Y4)
        Y.test4    <- Y4[ind.test4]
        X.test4    <- cbind(Xb[ind.test4,])
        
        obj.uneven <- aux.uneven(Y[-ind.test],cbind(Xa[-ind.test,]))
        obj.uneven2 <- aux.uneven(Y2[-ind.test2],cbind(Xa[-ind.test2,]))
        obj.uneven3 <- aux.uneven(Y3[-ind.test3],cbind(Xa[-ind.test3,]))
        obj.uneven4 <- aux.uneven(Y4[-ind.test4],cbind(Xb[-ind.test4,]))
        
        X0 <- cbind(obj.uneven$X0)
        Y0 <- obj.uneven$Y0
        X1 <- cbind(obj.uneven$X1)
        Y1 <- obj.uneven$Y1
        
        X02 <- cbind(obj.uneven2$X0)
        Y02 <- obj.uneven2$Y0
        X12 <- cbind(obj.uneven2$X1)
        Y12 <- obj.uneven2$Y1
        
        X03 <- cbind(obj.uneven3$X0)
        Y03 <- obj.uneven3$Y0
        X13 <- cbind(obj.uneven3$X1)
        Y13 <- obj.uneven3$Y1
        
        X04 <- cbind(obj.uneven4$X0)
        Y04 <- obj.uneven4$Y0
        X14 <- cbind(obj.uneven4$X1)
        Y14 <- obj.uneven4$Y1
        
        sing  <- is.singular.matrix(t(X0)%*%X0)
        sing2  <- is.singular.matrix(t(X02)%*%X02)
        sing3  <- is.singular.matrix(t(X03)%*%X03)
        sing4  <- is.singular.matrix(t(X04)%*%X04)
        
      }
        }
      }
    }
    
    ###Implementation of methods
    taus    <- seq(0.001,0.999,length=200)
  
    
    # Applying the difference conformal prediction methods
    res.qr    <- dcp.qr(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)
    res.qro<-dcp.opt(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)
    res.cqr   <- cqr(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
    res.reg   <- cp.reg(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
    
    res.qr2    <- dcp.qr(Y02,X02,Y12,X12,Y.test2,X.test2,taus,alpha.sig)
    res.qro2<-dcp.opt(Y02,X02,Y12,X12,Y.test2,X.test2,taus,alpha.sig)
    res.cqr2   <- cqr(Y02,X02,Y12,X12,Y.test2,X.test2,alpha.sig)
    res.reg2   <- cp.reg(Y02,X02,Y12,X12,Y.test2,X.test2,alpha.sig)
    
    res.qr3    <- dcp.qr(Y03,X03,Y13,X13,Y.test3,X.test3,taus,alpha.sig)
    res.qro3<-dcp.opt(Y03,X03,Y13,X13,Y.test3,X.test3,taus,alpha.sig)
    res.cqr3   <- cqr(Y03,X03,Y13,X13,Y.test3,X.test3,alpha.sig)
    res.reg3   <- cp.reg(Y03,X03,Y13,X13,Y.test3,X.test3,alpha.sig)
    
    res.qr4    <- dcp.qr(Y04,X04,Y14,X14,Y.test4,X.test4,taus,alpha.sig)
    res.qro4<-dcp.opt(Y04,X04,Y14,X14,Y.test4,X.test4,taus,alpha.sig)
    res.cqr4   <- cqr(Y04,X04,Y14,X14,Y.test4,X.test4,alpha.sig)
    res.reg4   <- cp.reg(Y04,X04,Y14,X14,Y.test4,X.test4,alpha.sig)
    
    X.test.mat <- rbind(X.test.mat,X.test)
    X.test.mat2 <- rbind(X.test.mat2,X.test2)
    X.test.mat3 <- rbind(X.test.mat3,X.test3)
    X.test.mat4 <- rbind(X.test.mat4,X.test4)
    
    Y.test.mat <- rbind(Y.test.mat,Y.test)
    Y.test.mat2 <- rbind(Y.test.mat2,Y.test2)
    Y.test.mat3 <- rbind(Y.test.mat3,Y.test3)
    Y.test.mat4 <- rbind(Y.test.mat4,Y.test4)
    
    X0.mat <- rbind(X0.mat,X0)
    X1.mat <- rbind(X1.mat,X1)
    Y0.mat <- rbind(Y0.mat,Y0)
    Y1.mat <- rbind(Y1.mat,Y1)
    
    X02.mat <- rbind(X02.mat,X02)
    X12.mat <- rbind(X12.mat,X12)
    Y02.mat <- rbind(Y02.mat,Y02)
    Y12.mat <- rbind(Y12.mat,Y12)
    
    X03.mat <- rbind(X03.mat,X03)
    X13.mat <- rbind(X13.mat,X13)
    Y03.mat <- rbind(Y03.mat,Y03)
    Y13.mat <- rbind(Y13.mat,Y13)
    
    X04.mat <- rbind(X04.mat,X04)
    X14.mat <- rbind(X14.mat,X14)
    Y04.mat <- rbind(Y04.mat,Y04)
    Y14.mat <- rbind(Y14.mat,Y14)
    
    cov.mat.temp  <- cbind(res.qr$cov.qr,res.qro$cov.opt,res.cqr$cov.o,res.cqr$cov.m,res.cqr$cov.r,res.reg$cov.reg)
    leng.mat.temp <- cbind(res.qr$leng.qr,res.qro$leng.opt,res.cqr$leng.o,res.cqr$leng.m,res.cqr$leng.r,res.reg$leng.reg)
    
    cov.mat.temp2  <- cbind(res.qr2$cov.qr,res.qro2$cov.opt,res.cqr2$cov.o,res.cqr2$cov.m,res.cqr2$cov.r,res.reg2$cov.reg)
    leng.mat.temp2 <- cbind(res.qr2$leng.qr,res.qro2$leng.opt,res.cqr2$leng.o,res.cqr2$leng.m,res.cqr2$leng.r,res.reg2$leng.reg)
    
    cov.mat.temp3  <- cbind(res.qr3$cov.qr,res.qro3$cov.opt,res.cqr3$cov.o,res.cqr3$cov.m,res.cqr3$cov.r,res.reg3$cov.reg)
    leng.mat.temp3 <- cbind(res.qr3$leng.qr,res.qro3$leng.opt,res.cqr3$leng.o,res.cqr3$leng.m,res.cqr3$leng.r,res.reg3$leng.reg)
    
    cov.mat.temp4  <- cbind(res.qr4$cov.qr,res.qro4$cov.opt,res.cqr4$cov.o,res.cqr4$cov.m,res.cqr4$cov.r,res.reg4$cov.reg)
    leng.mat.temp4 <- cbind(res.qr4$leng.qr,res.qro4$leng.opt,res.cqr4$leng.o,res.cqr4$leng.m,res.cqr4$leng.r,res.reg4$leng.reg)
    
    cov.mat   <- rbind(cov.mat,cov.mat.temp)
    leng.mat  <- rbind(leng.mat,leng.mat.temp)
    
    cov.mat2   <- rbind(cov.mat2,cov.mat.temp2)
    leng.mat2  <- rbind(leng.mat2,leng.mat.temp2)
    
    cov.mat3   <- rbind(cov.mat3,cov.mat.temp3)
    leng.mat3  <- rbind(leng.mat3,leng.mat.temp3)
    
    cov.mat4   <- rbind(cov.mat4,cov.mat.temp4)
    leng.mat4  <- rbind(leng.mat4,leng.mat.temp4)
    
  }
  ###Computing average coverage and average length
  cov <- colMeans(cov.mat,na.rm=TRUE)
  cov2 <- colMeans(cov.mat2,na.rm=TRUE)
  cov3 <- colMeans(cov.mat3,na.rm=TRUE)
  cov4 <- colMeans(cov.mat4,na.rm=TRUE)
  leng<-colMeans(leng.mat,na.rm=TRUE)
  leng2<-colMeans(leng.mat2,na.rm=TRUE)
  leng3<-colMeans(leng.mat3,na.rm=TRUE)
  leng4<-colMeans(leng.mat4,na.rm=TRUE)
  
  ###standard deviation of predicted conditional coverage probability
  
  pred.cov.qr       <- predict(glm(cov.mat[,1]~X.test.mat,family=binomial(link="logit")),type="response")
  pred.cov.qr.opt   <- predict(glm(cov.mat[,2]~X.test.mat,family=binomial(link="logit")),type="response")
  pred.cov.cqr.o    <- predict(glm(cov.mat[,3]~X.test.mat,family=binomial(link="logit")),type="response")
  pred.cov.cqr.m    <- predict(glm(cov.mat[,4]~X.test.mat,family=binomial(link="logit")),type="response")
  pred.cov.cqr.r    <- predict(glm(cov.mat[,5]~X.test.mat,family=binomial(link="logit")),type="response")
  pred.cov.reg      <- predict(glm(cov.mat[,6]~X.test.mat,family=binomial(link="logit")),type="response")
  
  pred.cov.qr2       <- predict(glm(cov.mat2[,1]~X.test.mat2,family=binomial(link="logit")),type="response")
  pred.cov.qr2.opt  <- predict(glm(cov.mat2[,2]~X.test.mat2,family=binomial(link="logit")),type="response")
  pred.cov.cqr2.o    <- predict(glm(cov.mat2[,3]~X.test.mat2,family=binomial(link="logit")),type="response")
  pred.cov.cqr2.m    <- predict(glm(cov.mat2[,4]~X.test.mat2,family=binomial(link="logit")),type="response")
  pred.cov.cqr2.r    <- predict(glm(cov.mat2[,5]~X.test.mat2,family=binomial(link="logit")),type="response")
  pred.cov.reg2      <- predict(glm(cov.mat2[,6]~X.test.mat2,family=binomial(link="logit")),type="response")
  
  pred.cov.qr3       <- predict(glm(cov.mat3[,1]~X.test.mat3,family=binomial(link="logit")),type="response")
  pred.cov.qr3.opt  <- predict(glm(cov.mat3[,2]~X.test.mat3,family=binomial(link="logit")),type="response")
  pred.cov.cqr3.o    <- predict(glm(cov.mat3[,3]~X.test.mat3,family=binomial(link="logit")),type="response")
  pred.cov.cqr3.m    <- predict(glm(cov.mat3[,4]~X.test.mat3,family=binomial(link="logit")),type="response")
  pred.cov.cqr3.r    <- predict(glm(cov.mat3[,5]~X.test.mat3,family=binomial(link="logit")),type="response")
  pred.cov.reg3      <- predict(glm(cov.mat3[,6]~X.test.mat3,family=binomial(link="logit")),type="response")
  
  pred.cov.qr4       <- predict(glm(cov.mat4[,1]~X.test.mat4,family=binomial(link="logit")),type="response")
  pred.cov.qr4.opt  <- predict(glm(cov.mat4[,2]~X.test.mat4,family=binomial(link="logit")),type="response")
  pred.cov.cqr4.o    <- predict(glm(cov.mat4[,3]~X.test.mat4,family=binomial(link="logit")),type="response")
  pred.cov.cqr4.m    <- predict(glm(cov.mat4[,4]~X.test.mat4,family=binomial(link="logit")),type="response")
  pred.cov.cqr4.r    <- predict(glm(cov.mat4[,5]~X.test.mat4,family=binomial(link="logit")),type="response")
  pred.cov.reg4      <- predict(glm(cov.mat4[,6]~X.test.mat4,family=binomial(link="logit")),type="response")
  
  
  
  pred.cov <- cbind(pred.cov.qr,pred.cov.qr.opt,pred.cov.cqr.o,pred.cov.cqr.m,pred.cov.cqr.r,pred.cov.reg)
  
  pred.cov2 <- cbind(pred.cov.qr2,pred.cov.qr2.opt,pred.cov.cqr2.o,pred.cov.cqr2.m,pred.cov.cqr2.r,pred.cov.reg2)
  
  pred.cov3 <- cbind(pred.cov.qr3,pred.cov.qr3.opt,pred.cov.cqr3.o,pred.cov.cqr3.m,pred.cov.cqr3.r,pred.cov.reg3)
  
  pred.cov4 <- cbind(pred.cov.qr4,pred.cov.qr4.opt,pred.cov.cqr4.o,pred.cov.cqr4.m,pred.cov.cqr4.r,pred.cov.reg4)
  
  res.cov.cond <- apply(pred.cov,2,sd)
  
  res.cov.cond2 <- apply(pred.cov2,2,sd)
  
  res.cov.cond3 <- apply(pred.cov3,2,sd)
  
  res.cov.cond4 <- apply(pred.cov4,2,sd)
  
  l<-list(cov,cov2,cov3,cov4,res.cov.cond,res.cov.cond2,res.cov.cond3,res.cov.cond4,leng.mat,leng.mat2,leng.mat3,leng.mat4,
          cov.mat,cov.mat2,cov.mat3,cov.mat4,pred.cov.qr,pred.cov.qr.opt,pred.cov.cqr.o,pred.cov.cqr.m,pred.cov.cqr.r,pred.cov.reg,
          pred.cov.qr2,pred.cov.qr2.opt,pred.cov.cqr2.o,pred.cov.cqr2.m,pred.cov.cqr2.r,pred.cov.reg2,
          pred.cov.qr3,pred.cov.qr3.opt,pred.cov.cqr3.o,pred.cov.cqr3.m,pred.cov.cqr3.r,pred.cov.reg3,
          pred.cov.qr4,pred.cov.qr4.opt,pred.cov.cqr4.o,pred.cov.cqr4.m,pred.cov.cqr4.r,pred.cov.reg4,
          X0.mat,Y0.mat,X1.mat,Y1.mat,
          X02.mat,Y02.mat,X12.mat,Y12.mat,
          X03.mat,Y03.mat,X13.mat,Y13.mat,
          X04.mat,Y04.mat,X14.mat,Y14.mat,
          X.test.mat,X.test.mat2,X.test.mat3,X.test.mat4,
          Y.test.mat,Y.test.mat2,Y.test.mat3,Y.test.mat4)
  return(l)
}


###Running algsine function on different sample sizes
###Samples of sizes 50,70,100,200,300,500 and 1000
########################################################

start<-Sys.time()
s<-algsine(50)  #4.17575 mins
Sys.time()-start

library(readxl)
k00<- lapply(5:8, function(i) read_excel("n0.xlsx", sheet = i))
d<-data.frame(rbind(t(data.frame(k00[1])),t(data.frame(k00[2])),t(data.frame(k00[3])),t(data.frame(k00[4]))))
names(d)<-c("pred.cov.qr","pred.cov.qr.opt","pred.cov.cqr.o","pred.cov.cqr.m",  
            "pred.cov.cqr.r","pred.cov.reg")
dt<-t(d)

start<-Sys.time()
k0<-algsine(70)  ###3.502714 mins
Sys.time()-start

start<-Sys.time()
k<-algsine(100)  ###4.551998 mins
Sys.time()-start

start<-Sys.time()
k2<-algsine(200) ###7.73741 mins 
Sys.time()-start

start<-Sys.time()
k3<-algsine(300) ###12.25624 mins 
Sys.time()-start

start<-Sys.time()
k4<-algsine(500)  ###26.75069 mins 
Sys.time()-start

start<-Sys.time()
k5<-algsine(1000)  ##38.72706 mins  
Sys.time()-start
#################################################
##################################################
##Avg.Length for all 4 specifications
avgl<-function(d){
  l1<-colMeans(data.frame(d[9]),na.rm=TRUE)
  l2<-colMeans(data.frame(d[10]),na.rm=TRUE)
  l3<-colMeans(data.frame(d[11]),na.rm=TRUE)
  l4<-colMeans(data.frame(d[12]),na.rm=TRUE)
  dat<-data.frame(rbind(l1,l2,l3,l4))
  dat$spec<-c("1","2","3","4")
  names(dat)<-c("DCP-QR","DCP-opt","CQR","CQR-m","CQR-r","CP-reg","specification")
  return(dat)
}
################################################
##Avg.coverage for all 4 specifications
avgc<-function(d){
  l1<-colMeans(data.frame(d[13]),na.rm=TRUE)
  l2<-colMeans(data.frame(d[14]),na.rm=TRUE)
  l3<-colMeans(data.frame(d[15]),na.rm=TRUE)
  l4<-colMeans(data.frame(d[16]),na.rm=TRUE)
  dat<-data.frame(rbind(l1,l2,l3,l4))
  dat$spec<-c("1","2","3","4")
  names(dat)<-c("DCP-QR","DCP-opt","CQR","CQR-m","CQR-r","CP-reg","specification")
  return(dat)
}
###################################################
###################################################
#####################################################
a<-avgl(s)
b<-avgl(k0)
c<-avgl(k)
d<-avgl(k2)
e<-avgl(k3)
f<-avgl(k4)
g<-avgl(k5)


a1<-avgc(s)
b1<-avgc(k0)
c1<-avgc(k)
d1<-avgc(k2)
e1<-avgc(k3)
f1<-avgc(k4)
g1<-avgc(k5)


###Table for average length
al<-rbind.data.frame(a,b,c,d,e,f,g)
xtable(al,digits=3)

##Table for average coverage
ac<-rbind.data.frame(a1,b1,c1,d1,e1,f1,g1)

###Standard deviation of predicted coverage probability
sd<-function(d){
  l1<-data.frame(d[5])
  l2<-data.frame(d[6])
  l3<-data.frame(d[7])
  l4<-data.frame(d[8])
  dat<-cbind.data.frame(l1,l2,l3,l4)
  names(dat)<-c("1","2","3","4")
  return(dat)
}


###Table for standard deviation
a2<-sd(s)
dt<-data.frame(dt)
names(dt)<-c("1","2","3","4")
b2<-sd(k0)
c2<-sd(k)
d2<-sd(k2)
e2<-sd(k3)
f2<-sd(k4)
g2<-sd(k5)
dsd<-rbind.data.frame(t(a2),t(b2),t(c2),t(d2),t(e2),t(f2),t(g2))


###scientific notation format function for very small standard deviation values 
formatn<-function(x){
  #if (x>0.09){y=x}
  #else {y=formatC(x,format="e",digits=2)}
  y<-ifelse(abs(x)>=0.01,round(x,3),formatC(x,format="e",digits=3))
  return(y)
}
dsd[]<-lapply(dsd,formatn)
dsd$type<-rep(c(1,2,3,4),7)
xtable(dsd)

###############QUANTILE BINS#########################################

###coverage is in positions 13,14,15,16 for list obtained from algsine()
###length is in positions 9,10,11,12
##X.test.mat,X.test.mat2,X.test.mat3,X.test.mat4 positions: 57,58,59,60

####Drawing quantile bins for conditional coverage and length 

#https://statisticsglobe.com/add-common-legend-to-combined-ggplot2-plots-in-r/#example-2-add-shared-legend-to-ggplot2-plots-using-gridextra-package
###user-defined function for created shared legend for grid.arrange() function under 
###gridExtra
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

qbc<-function(d,i,j,num.seg){
  X1<-as.vector(unlist(d[i]))
  cov1<-data.frame(d[j])
  cov.cond1<-binning(X1,cov1,num.seg)$cond
  cov.cond1<-data.frame(cov.cond1)
  c1<-cbind.data.frame(1:num.seg,cov.cond1)
  n1<-length(X1)/(1000*0.2)
  spec1<-1
  
  X2<-as.vector(unlist(d[i+1]))
  cov2<-data.frame(d[j+1])
  cov.cond2<-binning(X2,cov2,num.seg)$cond
  cov.cond2<-data.frame(cov.cond2)
  c2<-cbind.data.frame(1:num.seg,cov.cond2)
  n2<-length(X2)/(1000*0.2)
  spec2<-2
  
  X3<-as.vector(unlist(d[i+2]))
  cov3<-data.frame(d[j+2])
  cov.cond3<-binning(X2,cov3,num.seg)$cond
  cov.cond3<-data.frame(cov.cond3)
  c3<-cbind.data.frame(1:num.seg,cov.cond3)
  n3<-length(X3)/(1000*0.2)
  spec3<-3
  
  X4<-as.vector(unlist(d[i+3]))
  cov4<-data.frame(d[j+3])
  cov.cond4<-binning(X4,cov4,num.seg)$cond
  cov.cond4<-data.frame(cov.cond4)
  c4<-cbind.data.frame(1:num.seg,cov.cond4)
  n4<-length(X4)/(1000*0.2)
  spec4<-4
  
  names(c1)<-names(c2)<-names(c3)<-names(c4)<-c("Bin","DCP-QR","DCP-opt","CQR","CQR-m","CQR-r","CP-reg")
  
  c11<- melt(c1 ,  id.vars = 'Bin', variable.name = 'Method')
  c22<- melt(c2 ,  id.vars = 'Bin', variable.name = 'Method')
  c33<- melt(c3 ,  id.vars = 'Bin', variable.name = 'Method')
  c44<- melt(c4 ,  id.vars = 'Bin', variable.name = 'Method')
  
  p1<-ggplot(c11, aes(Bin, value)) +
    geom_line(aes(colour = Method),linetype="dashed",linewidth=1.25)+
    geom_hline(yintercept = 0.9,linewidth=1)+scale_colour_brewer(palette="Dark2") +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    ggtitle(paste("Sample Size:",n1,",Specification:",spec1))+xlab("Quantile Bin")+
    ylab("Conditional Coverage")+theme(legend.position = "none")
  
  p2<-ggplot(c22, aes(Bin, value)) +
    geom_line(aes(colour = Method),linetype="dashed",linewidth=1.25)+
    geom_hline(yintercept = 0.9,linewidth=1)+scale_colour_brewer(palette="Dark2") +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    ggtitle(paste("Sample Size:",n2,",Specification:",spec2))+xlab("Quantile Bin")+
    ylab("Conditional Coverage")+theme(legend.position = "none")
  
  p3<-ggplot(c33, aes(Bin, value)) +
    geom_line(aes(colour = Method),linetype="dashed",linewidth=1.25)+
    geom_hline(yintercept = 0.9,linewidth=1)+scale_colour_brewer(palette="Dark2") +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    ggtitle(paste("Sample Size:",n3,",Specification:",spec3))+xlab("Quantile Bin")+
    ylab("Conditional Coverage")+theme(legend.position = "none")
  
  p4<-ggplot(c44, aes(Bin, value)) +
    geom_line(aes(colour = Method),linetype="dashed",linewidth=1.25)+
    geom_hline(yintercept = 0.9,linewidth=1)+scale_colour_brewer(palette="Dark2") +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    ggtitle(paste("Sample Size:",n4,",Specification:",spec4))+xlab("Quantile Bin")+
    ylab("Conditional Coverage")+theme(legend.position = "none")
  
  p4_legend <-ggplot(c44, aes(Bin, value)) +
    geom_line(aes(colour = Method),linetype="dashed",linewidth=1.25)+
    geom_hline(yintercept = 0.9,linewidth=1)+scale_colour_brewer(palette="Dark2") +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    ggtitle(paste("Sample Size:",n4,",Specification:",spec4))+xlab("Quantile Bin")+
    ylab("Conditional Coverage")+theme(legend.position = "right")
  
  shared_legend <- extract_legend(p4_legend)
  
  grid.arrange(arrangeGrob(p1,p2,p3,p4,ncol=4),
               shared_legend)
  
  k<-list(p1,p2,p3,p4,c1,c2,c3,c4)
  return(k)
}

####Getting the results
a0<-qbc(s,57,13,10)
a<-qbc(k0,57,13,10)  
b<-qbc(k,57,13,10)
c<-qbc(k2,57,13,10)
d<-qbc(k3,57,13,10)
e<-qbc(k4,57,13,10)
f<-qbc(k5,57,13,10)


###conditional length
qbl<-function(d,i,j,num.seg){  ##i=1
  X1<-as.vector(unlist(d[i]))
  leng1<-data.frame(d[j])
  leng.cond1<-binning(X1,leng1,num.seg)$cond
  leng.cond1<-data.frame(leng.cond1)
  c1<-cbind.data.frame(1:num.seg,leng.cond1)
  n1<-length(X1)/(1000*0.2)
  spec1<-1
  
  X2<-as.vector(unlist(d[i+1]))
  leng2<-data.frame(d[j+1])
  leng.cond2<-binning(X2,leng2,num.seg)$cond
  leng.cond2<-data.frame(leng.cond2)
  c2<-cbind.data.frame(1:num.seg,leng.cond2)
  n2<-length(X2)/(1000*0.2)
  spec2<-2
  
  X3<-as.vector(unlist(d[i+2]))
  leng3<-data.frame(d[j+2])
  leng.cond3<-binning(X2,leng3,num.seg)$cond
  leng.cond3<-data.frame(leng.cond3)
  c3<-cbind.data.frame(1:num.seg,leng.cond3)
  n3<-length(X3)/(1000*0.2)
  spec3<-3
  
  X4<-as.vector(unlist(d[i+3]))
  leng4<-data.frame(d[j+3])
  leng.cond4<-binning(X4,leng4,num.seg)$cond
  leng.cond4<-data.frame(leng.cond4)
  c4<-cbind.data.frame(1:num.seg,leng.cond4)
  n4<-length(X4)/(1000*0.2)
  spec4<-4
  
  names(c1)<-names(c2)<-names(c3)<-names(c4)<-c("Bin","DCP-QR","DCP-opt","CQR","CQR-m","CQR-r","CP-reg")
  
  c11<- melt(c1 ,  id.vars = 'Bin', variable.name = 'Method')
  c22<- melt(c2 ,  id.vars = 'Bin', variable.name = 'Method')
  c33<- melt(c3 ,  id.vars = 'Bin', variable.name = 'Method')
  c44<- melt(c4 ,  id.vars = 'Bin', variable.name = 'Method')
  
  p1<-ggplot(c11, aes(Bin, value)) +
    geom_line(aes(colour = Method),linewidth=1.25)+
    scale_colour_brewer(palette="Dark2") +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    ggtitle(paste("Sample Size:",n1,",Specification:",spec1))+xlab("Quantile Bin")+
    ylab("Conditional Length")+theme(legend.position = "none")
  
  p2<-ggplot(c22, aes(Bin, value)) +
    geom_line(aes(colour = Method),linewidth=1.25)+
    scale_colour_brewer(palette="Dark2") +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    ggtitle(paste("Sample Size:",n2,",Specification:",spec2))+xlab("Quantile Bin")+
    ylab("Conditional Length")+theme(legend.position = "none")
  
  p3<-ggplot(c33, aes(Bin, value)) +
    geom_line(aes(colour = Method),linewidth=1.25)+
    scale_colour_brewer(palette="Dark2") +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    ggtitle(paste("Sample Size:",n3,",Specification:",spec3))+xlab("Quantile Bin")+
    ylab("Conditional Length")+theme(legend.position = "none")
  
  p4<-ggplot(c44, aes(Bin, value)) +
    geom_line(aes(colour = Method),linewidth=1.25)+
    scale_colour_brewer(palette="Dark2") +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    ggtitle(paste("Sample Size:",n4,",Specification:",spec4))+xlab("Quantile Bin")+
    ylab("Conditional Length")+theme(legend.position = "none")
  
  p4_legend <-ggplot(c44, aes(Bin, value)) +
    geom_line(aes(colour = Method),linewidth=1.25)+
    scale_colour_brewer(palette="Dark2") +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    ggtitle(paste("Sample Size:",n4,",Specification:",spec4))+xlab("Quantile Bin")+
    ylab("Conditional Length")+theme(legend.position = "right")
  
  shared_legend <- extract_legend(p4_legend)
  
  grid.arrange(arrangeGrob(p1,p2,p3,p4,ncol=4),
               shared_legend)
  
  k<-list(p1,p2,p3,p4,c1,c2,c3,c4)
  return(k)
}


p0<-qbl(s,57,9,10)
p<-qbl(k0,57,9,10)  
q<-qbl(k,57,9,10)
r<-qbl(k2,57,9,10)
s<-qbl(k3,57,9,10)
t<-qbl(k4,57,9,10)
u<-qbl(k5,57,9,10)

p0
p
q
r
s
t
u

write.xlsx(s,"s.xlsx")
write.xlsx(k0,"k0.xlsx")
write.xlsx(k,"k.xlsx")
write.xlsx(k2,"k2.xlsx")
write.xlsx(k3,"k3.xlsx")
write.xlsx(k4,"k4.xlsx")
write.xlsx(k5,"k5.xlsx")


##Avg Coverage
#c<-cbind.data.frame(k[1],k[2],k[3],k[4])
#names(c)<-c("c1","c2","c3","c4")

##sd 
#r<-cbind.data.frame(k[5],k[6],k[7],k[8])
#names(r)<-c("r1","r2","r3","r4")

##Merging all 
#all<-cbind.data.frame(l1,l2,l3,l4,c,r)

##rounding digits and using scientific notation for extremely small values
##just to help compare methods for small differences
#formatn<-function(x){
  #if (x>0.09){y=x}
  #else {y=formatC(x,format="e",digits=2)}
 # y<-ifelse(abs(x)>=0.01,round(x,3),formatC(x,format="e",digits=2))
  #return(y)
#}
#all[]<-lapply(all,formatn)


##taking to latex
#print(xtable(all,include.rownames=F))


###Avg.coverage for each iteration
#cov1<-data.frame(k[13])
#cov2<-data.frame(k[14])
#cov3<-data.frame(k[15])
#cov4<-data.frame(k[16])

##Empirical coverage
#empcov<-function(n,m){
#sv<-NULL
#for (i in 1:1000){
 #    subval<-m[n*(i-1)+1:n*i,]
  #   subval_avg<-colMeans(subval,na.rm=TRUE)
   #  sv<-rbind(sv,subval_avg)
#}
 #    return(sv)
#}
##for correct coverage, avg.coverage over each iteration
##should be around (1-alpha) and average of these averages 
##should be around (1-alpha).
#a1<-empcov(20,cov1)  ##for test sample size of 20
#a2<-empcov(20,cov2)  ##for test sample size of 20
#a3<-empcov(20,cov3)  ##for test sample size of 20
#a4<-empcov(20,cov4)  ##for test sample size of 20

#aa<-data.frame(colMeans(a1,na.rm=TRUE))
#ab<-data.frame(colMeans(a2,na.rm=TRUE))
#ac<-data.frame(colMeans(a3,na.rm=TRUE))
#ad<-data.frame(colMeans(a4,na.rm=TRUE))
#A<-cbind.data.frame(aa,ab,ac,ad)
#print(xtable(A,digits=4))


#p = seq(0,1, length=100)
#l<-function(n,alpha){
 # u<-floor((n+1)*alpha)
  #return(u)
#}
##parameters of beta: Beta(n+1-l,l)


#hist()
#plot(p, dbeta(p, 111-10, 10), type='l')
#c_bar
#colMeans(a1,na.rm=TRUE)

##saving results for methods in Python 
#X0<-data.frame(k[41])
#Y0a<-data.frame(k[42])
#Y0<-data.frame(as.vector(t(Y0a)))
#X1<-data.frame(k[43])
#Y1a<-data.frame(k[44])
#Y1<-data.frame(as.vector(t(Y1a)))
#X02<-data.frame(k[45])
#Y02a<-data.frame(k[46])
#Y02<-data.frame(as.vector(t(Y02a)))
#X12<-data.frame(k[47])
#Y12a<-data.frame(k[48])
#Y12<-data.frame(as.vector(t(Y12a)))
#X03<-data.frame(k[49])
#Y03a<-data.frame(k[50])
#Y03<-data.frame(as.vector(t(Y03a)))
#X13<-data.frame(k[51])
#Y13a<-data.frame(k[52])
#Y13<-data.frame(as.vector(t(Y13a)))
#X04<-data.frame(k[53])
#Y04a<-data.frame(k[54])
#Y04<-data.frame(as.vector(t(Y04a)))
#X14<-data.frame(k[55])
#Y14a<-data.frame(k[56])
#Y14<-data.frame(as.vector(t(Y14a)))
#Xtest<-data.frame(k[57])
#Xtest2<-data.frame(k[58])
#Xtest3<-data.frame(k[59])
#Xtest4<-data.frame(k[60])
#Ytesta<-data.frame(k[61])
#Ytest<-data.frame(as.vector(t(Ytesta)))
#Ytestb<-data.frame(k[62])
#Ytest2<-data.frame(as.vector(t(Ytestb)))
#Ytestc<-data.frame(k[63])
#Ytest3<-data.frame(as.vector(t(Ytestc)))
#Ytestd<-data.frame(k[64])
#Ytest4<-data.frame(as.vector(t(Ytestd)))

#library(writexl)
##Example-1
#names(X0)<-c("X")
#names(Y0)<-c("Y")
#t100<-write_xlsx(cbind.data.frame(X0,Y0),"t100.xlsx")
#names(X1)<-c("X1")
#names(Y1)<-c("Y1")
#v100<-write_xlsx(cbind.data.frame(X1,Y1),"v100.xlsx")
#names(Xtest)<-c("Xtest")
#names(Ytest)<-c("Ytest")
#ts100<-write_xlsx(cbind.data.frame(Xtest,Ytest),"ts100.xlsx")

##Example-2
#names(X02)<-c("X02")
#names(Y02)<-c("Y02")
#t2100<-write_xlsx(cbind.data.frame(X02,Y02),"t2100.xlsx")
#names(X12)<-c("X12")
#names(Y12)<-c("Y12")
#v2100<-write_xlsx(cbind.data.frame(X12,Y12),"v2100.xlsx")
#names(Xtest2)<-c("Xtest2")
#names(Ytest2)<-c("Ytest2")
#ts2100<-write_xlsx(cbind.data.frame(Xtest2,Ytest2),"ts2100.xlsx")

##Example-3
#names(X03)<-c("X03")
#names(Y03)<-c("Y03")
#t3100<-write_xlsx(cbind.data.frame(X03,Y03),"t3100.xlsx")
#names(X13)<-c("X13")
#names(Y13)<-c("Y13")
#v3100<-write_xlsx(cbind.data.frame(X13,Y13),"v3100.xlsx")
#names(Xtest3)<-c("Xtest3")
#names(Ytest3)<-c("Ytest3")
#ts3100<-write_xlsx(cbind.data.frame(Xtest3,Ytest3),"ts3100.xlsx")

##Example-4
#names(X04)<-c("X04")
#names(Y04)<-c("Y04")
#t4100<-write_xlsx(cbind.data.frame(X04,Y04),"t4100.xlsx")
#names(X14)<-c("X14")
#names(Y14)<-c("Y14")
#v4100<-write_xlsx(cbind.data.frame(X14,Y14),"v4100.xlsx")
#names(Xtest4)<-c("Xtest4")
#names(Ytest4)<-c("Ytest4")
#ts4100<-write_xlsx(cbind.data.frame(Xtest4,Ytest4),"ts4100.xlsx")


#k1<-algsine(200) ##Test size:0.2*200=40 per iteration
##get empirical coverage for k1,k2,k3,k4 and k5 as well
#covb1<-data.frame(k1[13])
#covb2<-data.frame(k1[14])
#covb3<-data.frame(k1[15])
#covb4<-data.frame(k1[16])



#b1<-empcov(40,covb1)
#b2<-empcov(40,covb2)
#b3<-empcov(40,covb3)
#b4<-empcov(40,covb4)
##cbar
#ba<-data.frame(colMeans(b1,na.rm=TRUE))
#bb<-data.frame(colMeans(b2,na.rm=TRUE))
#bc<-data.frame(colMeans(b3,na.rm=TRUE))
#bd<-data.frame(colMeans(b4,na.rm=TRUE))
#B<-cbind.data.frame(ba,bb,bc,bd)
#print(xtable(B,row.names=FALSE,digits=4))
##histogram with c_bar and corresponding beta distribution
#k2<-algsine(300)  ##test size: 300*0.2=60
#covc1<-data.frame(k2[13])
#covc2<-data.frame(k2[14])
#covc3<-data.frame(k2[15])
#covc4<-data.frame(k2[16])

#c1<-empcov(60,covc1)
#c2<-empcov(60,covc2)
#c3<-empcov(60,covc3)
#c4<-empcov(60,covc4)
#colMeans(c1,na.rm=TRUE)
#colMeans(c2,na.rm=TRUE)
#colMeans(c3,na.rm=TRUE)
#colMeans(c4,na.rm=TRUE)
#ca<-data.frame(colMeans(c1,na.rm=TRUE))
#cb<-data.frame(colMeans(c2,na.rm=TRUE))
#cc<-data.frame(colMeans(c3,na.rm=TRUE))
#cd<-data.frame(colMeans(c4,na.rm=TRUE))
#C<-cbind.data.frame(ca,cb,cc,cd)
#print(xtable(B,row.names=FALSE,digits=4))
##histogram with c_bar and corresponding beta distribution
k3<-algsine(500)  ##500*0.2
covd1<-data.frame(k3[13])
covd2<-data.frame(k3[14])
covd3<-data.frame(k3[15])
covd4<-data.frame(k3[16])

d1<-empcov(100,covd1)
d2<-empcov(100,covd2)
d3<-empcov(100,covd3)
d4<-empcov(100,covd4)
colMeans(d1,na.rm=TRUE)
colMeans(d2,na.rm=TRUE)
colMeans(d3,na.rm=TRUE)
colMeans(d4,na.rm=TRUE)
da<-data.frame(colMeans(d1,na.rm=TRUE))
db<-data.frame(colMeans(d2,na.rm=TRUE))
dc<-data.frame(colMeans(d3,na.rm=TRUE))
dd<-data.frame(colMeans(d4,na.rm=TRUE))
D<-cbind.data.frame(da,db,dc,dd)
print(xtable(D,row.names=FALSE,digits=4))


k4<-algsine(700)  ##700*0.2
cove1<-data.frame(k4[13])
cove2<-data.frame(k4[14])
cove3<-data.frame(k4[15])
cove4<-data.frame(k4[16])
e1<-empcov(140,cove1)
e2<-empcov(140,cove2)
e3<-empcov(140,cove3)
e4<-empcov(140,cove4)
colMeans(e1,na.rm=TRUE)
colMeans(e2,na.rm=TRUE)
colMeans(e3,na.rm=TRUE)
colMeans(e4,na.rm=TRUE)
ea<-data.frame(colMeans(e1,na.rm=TRUE))
eb<-data.frame(colMeans(e2,na.rm=TRUE))
ec<-data.frame(colMeans(e3,na.rm=TRUE))
ed<-data.frame(colMeans(e4,na.rm=TRUE))
E<-cbind.data.frame(ea,eb,ec,ed)
print(xtable(E,row.names=FALSE,digits=4))

##histogram with c_bar and corresponding beta distribution
k5<-algsine(1000)  ##100*0.2
covf1<-data.frame(k5[13])
covf2<-data.frame(k5[14])
covf3<-data.frame(k5[15])
covf4<-data.frame(k5[16])
f1<-empcov(200,covf1)
f2<-empcov(200,covf2)
f3<-empcov(200,covf3)
f4<-empcov(200,covf4)
colMeans(f1,na.rm=TRUE)
colMeans(f2,na.rm=TRUE)
colMeans(f3,na.rm=TRUE)
colMeans(f4,na.rm=TRUE)
fa<-data.frame(colMeans(f1,na.rm=TRUE))
fb<-data.frame(colMeans(f2,na.rm=TRUE))
fc<-data.frame(colMeans(f3,na.rm=TRUE))
fd<-data.frame(colMeans(f4,na.rm=TRUE))
F<-cbind.data.frame(fa,fb,fc,fd)
print(xtable(F,row.names=FALSE,digits=4))
##since functions are simple, we can also do sample size of 10,000
k6<-algsine(10000)  ##100*0.2
covg1<-data.frame(k6[13])
covg2<-data.frame(k6[14])
covg3<-data.frame(k6[15])
covg4<-data.frame(k6[16])
g1<-empcov(2000,covg1)
g2<-empcov(2000,covg2)
g3<-empcov(2000,covg3)
g4<-empcov(2000,covg4)
colMeans(g1,na.rm=TRUE)
colMeans(g2,na.rm=TRUE)
colMeans(g3,na.rm=TRUE)
colMeans(g4,na.rm=TRUE)
ga<-data.frame(colMeans(g1,na.rm=TRUE))
gb<-data.frame(colMeans(g2,na.rm=TRUE))
gc<-data.frame(colMeans(g3,na.rm=TRUE))
gd<-data.frame(colMeans(g4,na.rm=TRUE))
G<-cbind.data.frame(ga,gb,gc,gd)
print(xtable(G,row.names=FALSE,digits=4))

k7<-algsine(5000)  ##100*0.2
covk1<-data.frame(k7[13])
covk2<-data.frame(k7[14])
covk3<-data.frame(k7[15])
covk4<-data.frame(k7[16])
k1<-empcov(1000,covk1)
k2<-empcov(1000,covk2)
k3<-empcov(1000,covk3)
k4<-empcov(1000,covk4)
colMeans(k1,na.rm=TRUE)
colMeans(k2,na.rm=TRUE)
colMeans(k3,na.rm=TRUE)
colMeans(k4,na.rm=TRUE)
ka<-data.frame(colMeans(k1,na.rm=TRUE))
kb<-data.frame(colMeans(k2,na.rm=TRUE))
kc<-data.frame(colMeans(k3,na.rm=TRUE))
kd<-data.frame(colMeans(k4,na.rm=TRUE))
K<-cbind.data.frame(ka,kb,kc,kd)
print(xtable(K,row.names=FALSE,digits=4))

sink("k.txt")
k
sink()

#r<- read.delim("k.txt")
# Separate elements by one or more whitepace
#w <- strsplit(r, "[[0:99]]+")