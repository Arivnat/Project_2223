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
#######################################################
###loading required packages













##############################################################################
#########SIMULATION EXAMPLE FROM KIRANOV ET.AL(2019)##########################
################################################################################


library(ggplot2)
library(quantreg)
library(matrixcalc)
library(hdm)
library(writexl)
library(readxl)
library(stringr)
library(readr)
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
############################################################################################
###########################################################################################
############################################################################################

###Function for data generation
###of homoskedastic and heteroskedastic errors
data<-function(n,s){
  Xmat<-r<-r2<-yh<-yho<-res<-NULL
  for (i in 1:s){
    set.seed(i)
    x<-matrix(runif(n=n*100, min=0, max=1), nrow=n)
    xd<-data.frame(x)
    beta<-c(rep(1,5),rep(0,95))
    xb<-as.vector(sqrt(1+(beta%*%t(x))^2))
    for (k in 1:n){
    res[k]<-rnorm(1,0,xb[k])
    }
    res2<-rnorm(n,0,1)
    t2<-as.matrix(res)
    t3<-as.matrix(res2)
    y<-2*sin(pi*beta%*%t(x))+(pi*beta%*%t(x))+res  ##heteroskedastic 
    y2<-2*sin(pi*beta%*%t(x))+(pi*beta%*%t(x))+res2  ###homoskedastic
    y<-as.matrix(as.vector(y))
    y2<-as.matrix(as.vector(y2))
    Xmat<-rbind(Xmat,xd)
    r<-rbind(r,t2)
    r2<-rbind(r2,t3)
    yh<-rbind(yh,y)
    yho<-rbind(yho,y2)
  }
  l<-list(Xmat,yh,yho,t2)
  return(l)
}
##700 data points-200 iterations
a<-data(700,200)
x<-data.frame(a[1])   
y1<-data.frame(a[2])  ##140000 1
y2<-data.frame(a[3])  ##140000 1
names(y1)<-names(y2)<-names(x)<-NULL   ##y1[1,]


##index generating function for train and calibration data modify:
###modification of aux.uneven function in Chernozhukov.et.al(2021)
ax <- function(Y,X){
  T01     <- dim(Y)[1]
  T01even <- floor(T01/2)*2
  
  Yeven <- Y[(T01-T01even+1):T01,]
  Xeven <- cbind(X[(T01-T01even+1):T01,])
  
  ind0 <- 1:(T01even/2)
  ind1 <- (T01even/2+1):T01even
  
  Y0 <- Yeven[ind0,]
  X0 <- cbind(Xeven[ind0,])
  Y1 <- Yeven[ind1,]
  X1 <- cbind(Xeven[ind1,])
  
  return(list(Y0=Y0,X0=X0,Y1=Y1,X1=X1))
  
}
alpha.sig<-0.1
T.ho <- floor(700*0.2) ##20% holdout sample  

##Use only 70,000 (100 iterations)
#reps<-function(n,nsims,alpha.sig){

###Function for implementation of methods 
cpo<-function(t,n,x,y,z){
cov.mat <- leng.mat <- X.test.mat <- NULL
cov.mat2 <- leng.mat2 <- X.test.mat2 <- NULL

sx0<-sy0<-sx1<-sy1<-stx<-sty<-NULL
sx02<-sy02<-sx12<-sy12<-stx2<-sty2<-NULL

start<-Sys.time()
for (i in 1:t){
  ya<-y1[n*(i-1)+1:n*i,]
  yb<-y2[n*(i-1)+1:n*i,]
  xo<-x[n*(i-1)+1:n*i,]
  sing <- TRUE
  sing2<-TRUE
  while (sing==TRUE){
    
  while (sing2==TRUE){
      set.seed(i)
      ind <- sample(n,n,replace=FALSE)
      Y   <- ya[ind,]
      Y2   <- y2[ind,]
      X   <- xo[ind,]
      
      ind.test  <- (dim(Y)[1]-T.ho+1):dim(Y)[1]
      Y.test    <- Y[ind.test,]
      X.test    <- cbind(X[ind.test,])
      
      ind.test2  <- (dim(Y2)[1]-T.ho+1):dim(Y2)[1]
      Y.test2    <- Y2[ind.test2,]
      X.test2    <- cbind(X[ind.test2,])
      
      obj.uneven <- ax(Y[-ind.test,],cbind(X[-ind.test,]))
      obj.uneven2 <- ax(Y2[-ind.test2,],cbind(X[-ind.test2,]))
      
      X0 <- cbind(obj.uneven$X0)
      Y0 <- obj.uneven$Y0
      X1 <- cbind(obj.uneven$X1)
      Y1 <- obj.uneven$Y1
      
      X02 <- cbind(obj.uneven2$X0)
      Y02 <- obj.uneven2$Y0
      X12 <- cbind(obj.uneven2$X1)
      Y12 <- obj.uneven2$Y1
      X0<-as.matrix(X0)
      X02<-as.matrix(X02)
      X1<-as.matrix(X1)
      X12<-as.matrix(X12)
      Y0<-as.matrix(Y0)
      Y02<-as.matrix(Y02)
      Y1<-as.matrix(Y1)
      Y12<-as.matrix(Y12)
      X.test<-as.matrix(X.test)
      X.test2<-as.matrix(X.test2)
      Y.test<-as.matrix(Y.test)
      Y.test2<-as.matrix(Y.test2)
      
      sing  <- is.singular.matrix(t(X0)%*%X0)
      sing2  <- is.singular.matrix(t(X02)%*%X02)
      
    }
  }
  
  taus    <- seq(0.001,0.999,length=200)
  #ys      <- quantile(unique(c(Y0,Y1)),seq(0.001,0.999,length=200))
  #ys2      <- quantile(unique(c(Y02,Y12)),seq(0.001,0.999,length=200))
  #Y0<-as.matrix(Y0)
  #Y1<-as.matrix(Y1)
  #Y02<-as.matrix(Y02)
  #Y12<-as.matrix(Y12)
  #Y.test<-as.matrix(Y.test)
  #Y.test2<-as.matrix(Y.test2)
  
  
  # Applying the difference conformal prediction methods
  res.qr    <- dcp.qr(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)
  res.qro<-dcp.opt(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)
  res.cqr   <- cqr(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
  res.reg   <- cp.reg(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
  
  res.qr2    <- dcp.qr(Y02,X02,Y12,X12,Y.test2,X.test2,taus,alpha.sig)
  res.qro2<-dcp.opt(Y02,X02,Y12,X12,Y.test2,X.test2,taus,alpha.sig)
  res.cqr2   <- cqr(Y02,X02,Y12,X12,Y.test2,X.test2,alpha.sig)
  res.reg2   <- cp.reg(Y02,X02,Y12,X12,Y.test2,X.test2,alpha.sig)
  
  X.test.mat <- rbind(X.test.mat,X.test)
  X.test.mat2 <- rbind(X.test.mat2,X.test2)
  
  cov.mat.temp  <- cbind(res.qr$cov.qr,res.qro$cov.opt,res.cqr$cov.o,res.cqr$cov.m,res.cqr$cov.r,res.reg$cov.reg)
  leng.mat.temp <- cbind(res.qr$leng.qr,res.qro$leng.opt,res.cqr$leng.o,res.cqr$leng.m,res.cqr$leng.r,res.reg$leng.reg)
  
  cov.mat.temp2  <- cbind(res.qr2$cov.qr,res.qro2$cov.opt,res.cqr2$cov.o,res.cqr2$cov.m,res.cqr2$cov.r,res.reg2$cov.reg)
  leng.mat.temp2 <- cbind(res.qr2$leng.qr,res.qro2$leng.opt,res.cqr2$leng.o,res.cqr2$leng.m,res.cqr2$leng.r,res.reg2$leng.reg)
  
  cov.mat   <- rbind(cov.mat,cov.mat.temp)
  leng.mat  <- rbind(leng.mat,leng.mat.temp)
  
  cov.mat2   <- rbind(cov.mat2,cov.mat.temp2)
  leng.mat2  <- rbind(leng.mat2,leng.mat.temp2)
  
  
  sx0<-rbind(sx0,X0)
  sy0<-rbind(sy0,Y0)
  sx1<-rbind(sx1,X1)
  sy1<-rbind(sy1,Y1)
  stx<-rbind(stx,X.test)
  sty<-rbind(sty,Y.test)
  sx02<-rbind(sx02,X02)
  sy02<-rbind(sy02,Y02)
  sx12<-rbind(sx12,X12)
  sy12<-rbind(sy12,Y12)
  stx2<-rbind(stx2,X.test2)
  sty2<-rbind(sty2,Y.test2)
}
j=Sys.time()-start  

###heteroskdastic case
cov1<-colMeans(cov.mat)  ##unconditional coverage 
leng1<-colMeans(leng.mat,na.rm=TRUE) ##Average length

##homo
cov2<-colMeans(cov.mat2)  ##unconditional coverage 
leng2<-colMeans(leng.mat2,na.rm=TRUE) ##Average length

pred.cov.qr       <- predict(glm(cov.mat[,1]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.qr.opt   <- predict(glm(cov.mat[,2]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.cqr.o    <- predict(glm(cov.mat[,3]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.cqr.m    <- predict(glm(cov.mat[,4]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.cqr.r    <- predict(glm(cov.mat[,5]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.reg      <- predict(glm(cov.mat[,6]~X.test.mat,family=binomial(link="logit")),type="response")

pred.cov.qr2       <- predict(glm(cov.mat2[,1]~X.test.mat2,family=binomial(link="logit")),type="response")
pred.cov.qr2.opt  <- predict(glm(cov.mat2[,2]~X.test.mat2,family=binomial(link="logit")),type="response")
pred.cov.cqr2.o    <- predict(glm(cov.mat2[,3]~X.test.mat2,family=binomial(link="logit")),type="response")
pred.cov.cqr2.m    <- predict(glm(cov.mat[,4]~X.test.mat2,family=binomial(link="logit")),type="response")
pred.cov.cqr2.r    <- predict(glm(cov.mat[,5]~X.test.mat2,family=binomial(link="logit")),type="response")
pred.cov.reg2      <- predict(glm(cov.mat2[,6]~X.test.mat2,family=binomial(link="logit")),type="response")

pred.cov <- cbind(pred.cov.qr,pred.cov.qr.opt,pred.cov.cqr.o,pred.cov.cqr.m,pred.cov.cqr.r,pred.cov.reg)

pred.cov2 <- cbind(pred.cov.qr2,pred.cov.qr2.opt,pred.cov.cqr2.o,pred.cov.cqr2.m,pred.cov.cqr2.r,pred.cov.reg2)

# for comparison sqrt(colMeans((pred.cov-0.9)^2))*100
res.cov.cond <- apply(pred.cov,2,sd)


res.cov.cond2 <- apply(pred.cov2,2,sd)

h<-list(cov1,cov2,leng1,leng2,res.cov.cond,res.cov.cond2,
        cov.mat,cov.mat2,leng.mat,leng.mat2,
        sx0,sy0,sx1,sy1,stx,sty,
        sx02,sy02,sx12,sy12,stx2,sty2,j)
return(h)
}

##sample size of 700 with 50 iterations
start<-Sys.time()
sin700<-cpo(50,700,x,y1,y2)  ###14.05808 mins
Sys.time()-start
write.xlsx(sin700,"sin700a.xlsx")

##############1000-100
a<-sin700[3]  ###avg.length-heteroskedastic case
b<-sin700[1]  ##avg.coverage-heteroskedastic case
c<-sin700[5]  ##standard deviation of conditional coverage probability-heteroskedastic case

a1<-sin700[4] ##avg length-homoskedastic case
b1<-sin700[2] ##avg.coverage-homoskedastic case
c1<-sin700[6]##standard deviation of conditional coverage probability-homoskedastic case

j<-cbind.data.frame(a,a1,b,b1,c,c1)
names(j)<-c("a","a1","b","b1","c","c1")
j<-round(j,4)


###1000-100
a<-data(1000,200)
x<-data.frame(a[1])
y1<-data.frame(a[2])
y2<-data.frame(a[3]) 
names(y1)<-names(y2)<-names(x)<-NULL 

alpha.sig<-0.1
T.ho <- floor(1000*0.2) ##20% holdout sample  ##try 3/4-1/4 split as well 

start<-Sys.time()
sin1000<-cpo(100,1000,x,y1,y2) ##19.42501 mins
Sys.time()-start 


##Set up a function that returns all results
##g is lsit returned from alg fitting function,n is sample size in it, it is # of iteration rounds 
res<-function(g,n,it){
  a<-g[3]
  b<-g[1]
  c<-g[5]
  
  a1<-g[4]
  b1<-g[2]
  c1<-g[6]
  
  j<-cbind.data.frame(a,a1,b,b1,c,c1)
  names(j)<-c("a","a1","b","b1","c","c1")
  j<-round(j,4)
  
  ###2000-200
  a<-data(2000,200)
  x<-data.frame(a[1])
  y1<-data.frame(a[2])
  y2<-data.frame(a[3]) 
  names(y1)<-names(y2)<-names(x)<-NULL 
  
  alpha.sig<-0.1
  T.ho <- floor(2000*0.2) ##20% holdout sample  ##try 3/4-1/4 split as well 
  
  start<-Sys.time()
  sin2000<-cpo(100,2000,x,y1,y2)
  Sys.time()-start  ###50.17938 mins
  write.xlsx(sin2000,"sin2000b.xlsx")
  

###############################################################################
  
  sin700_ds<-lapply(1:6, function(i) read_excel("sin700a.xlsx", sheet = i,col_names=FALSE))
  sin1000_ds<-lapply(1:6, function(i) read_excel("sin1000.xlsx", sheet = i,col_names=FALSE))
  sin2000_ds<-lapply(1:6, function(i) read_excel("sin2000b.xlsx", sheet = i,col_names=FALSE))
  
  ###Create table for avg coverage,avg length and sd of conditional coverage probability 
  a<-cbind.data.frame(sin700_ds[1],sin700_ds[2],sin700_ds[3],sin700_ds[4],sin700_ds[5],sin700_ds[6])
  b<-cbind.data.frame(sin1000_ds[1],sin1000_ds[2],sin1000_ds[3],sin1000_ds[4],sin1000_ds[5],sin1000_ds[6])
  c<-cbind.data.frame(sin2000_ds[1],sin2000_ds[2],sin2000_ds[3],sin2000_ds[4],sin2000_ds[5],sin2000_ds[6])
  
  d<-rbind.data.frame(a,b,c)
  names(d)<-c("Cov-1","Cov-2","Leng-1","Leng-2","SD-1","SD-2")
  xtable(d,digits=3)

  ####load data directly from saved xlsx worksheets
  setwd("C:/Users/RimJhim/Documents/doc")
  
  ####PER ITERATION AVERAGE COVERAGE PLOT
  s700<-lapply(7:8, function(i) read_excel("sin700b.xlsx", sheet = i))
  s1000<-lapply(7:8, function(i) read_excel("sin1000.xlsx", sheet = i))
  s2000<-lapply(7:8, function(i) read_excel("sin2000b.xlsx", sheet = i))
  
  ###Heteroskedastic case
  g700<-data.frame(s700[1])
  g1000<-data.frame(s1000[1])
  g2000<-data.frame(s2000[1])
  
  n0<-nrow(g700)/50
  n01<-nrow(g1000)/100
  n02<-nrow(g2000)/100
 
  g1<-data.frame(aggregate(g700, list(rep(1:(nrow(g700) %/% n0 + 1), each = n0, len = nrow(g700))), mean)[-1])
  g2<-data.frame(aggregate(g1000, list(rep(1:(nrow(g1000) %/% n01 + 1), each = n01, len = nrow(g1000))), mean)[-1])
  g3<-data.frame(aggregate(g2000, list(rep(1:(nrow(g2000) %/% n02 + 1), each = n02, len = nrow(g2000))), mean)[-1])
  
  ###Homoskedastic case
  g700b<-data.frame(s700[2])
  g1000b<-data.frame(s1000[2])
  g2000b<-data.frame(s2000[2])
 
 
  g4<-data.frame(aggregate(g700b, list(rep(1:(nrow(g700b) %/% n0 + 1), each = n0, len = nrow(g700b))), mean)[-1])
  g5<-data.frame(aggregate(g1000b, list(rep(1:(nrow(g1000b) %/% n01 + 1), each = n01, len = nrow(g1000b))), mean)[-1])
  g6<-data.frame(aggregate(g2000b, list(rep(1:(nrow(g2000b) %/% n02 + 1), each = n02, len = nrow(g2000b))), mean)[-1])
 
  g1$sample<-as.character(nrow(g700)/(50*0.2))
  g2$sample<-as.character(nrow(g1000)/(100*0.2))
  g3$sample<-as.character(nrow(g2000)/(100*0.2))
  g4$sample<-as.character(nrow(g700b)/(50*0.2))
  g5$sample<-as.character(nrow(g1000b)/(100*0.2))
  g6$sample<-as.character(nrow(g2000b)/(100*0.2))
  names(g1)<-names(g2)<-names(g3)<-names(g4)<-names(g5)<-names(g6)<-c("DCP-QR","DCP-opt",
                                                                      "CQR","CQR-m","CQR-r","CP-reg","Sample")
 
  
  ######################################################################################################################
  ########################################BOXPLOT FOR COVERAGE#########################################################
  
  plotsinc<-function(d1,d2,d3,d4,d5,d6){
    dat<-rbind.data.frame(d1,d2,d3)
    dat2<-rbind.data.frame(d4,d5,d6)
    dat$Sample<-factor(dat$Sample,levels=c("700","1000","2000"))
    dat2$Sample<-factor(dat2$Sample,levels=c("700","1000","2000"))
    datm<-melt(dat,id.vars="Sample")
    dat2m<-melt(dat2,id.vars="Sample")
    
    p1<-ggplot(datm,aes(x=Sample,y=value,col=variable))+facet_wrap(~variable)+
        geom_boxplot()+
        geom_hline(yintercept=0.9,col="darkgrey",linewidth=1)+
        stat_summary(fun = mean, geom = "point", col = "darkred") +  
        stat_summary(fun = mean, geom = "text", col = "red",     
                   vjust = 1.5, aes(label = paste(round(..y.., digits = 3))))+
      #theme_classic()+
        xlab("Sample Size")+ylab("Avg.Coverage")
    
    p2<-ggplot(dat2m,aes(x=Sample,y=value,col=variable))+facet_wrap(~variable)+
      geom_boxplot()+
      geom_hline(yintercept=0.9,col="darkgrey",linewidth=1)+
      stat_summary(fun = mean, geom = "point", col = "darkred") +  
      stat_summary(fun = mean, geom = "text", col = "red",     
                   vjust = 1.5, aes(label = paste(round(..y.., digits = 3))))+
      xlab("Sample Size")+ylab("Avg.Coverage")
    
    k<-list(p1,p2)
    return(k)
    
  }
  plotsinc(g1,g2,g3,g4,g5,g6)[1]
  plotsinc(g1,g2,g3,g4,g5,g6)[2]
#############################################################################
####AVERAGE LENGTH
  
  s700l<-lapply(9:10, function(i) read_excel("sin700a.xlsx", sheet = i))
  s1000l<-lapply(9:10, function(i) read_excel("sin1000.xlsx", sheet = i))
  s2000l<-lapply(9:10, function(i) read_excel("sin2000b.xlsx", sheet = i))
  
  ###Heteroskedastic case
  g700<-data.frame(s700l[1])
  g1000<-data.frame(s1000l[1])
  g2000<-data.frame(s2000l[1])
  
  n0<-nrow(g700)/50
  n01<-nrow(g1000)/100
  n02<-nrow(g2000)/100
  
  g1<-data.frame(aggregate(g700, list(rep(1:(nrow(g700) %/% n0 + 1), each = n0, len = nrow(g700))), mean)[-1])
  g2<-data.frame(aggregate(g1000, list(rep(1:(nrow(g1000) %/% n01 + 1), each = n01, len = nrow(g1000))), mean)[-1])
  g3<-data.frame(aggregate(g2000, list(rep(1:(nrow(g2000) %/% n02 + 1), each = n02, len = nrow(g2000))), mean)[-1])
  
  ###Homoskedastic case
  g700b<-data.frame(s700l[2])
  g1000b<-data.frame(s1000l[2])
  g2000b<-data.frame(s2000l[2])
  
  
  g4<-data.frame(aggregate(g700b, list(rep(1:(nrow(g700b) %/% n0 + 1), each = n0, len = nrow(g700b))), mean)[-1])
  g5<-data.frame(aggregate(g1000b, list(rep(1:(nrow(g1000b) %/% n01 + 1), each = n01, len = nrow(g1000b))), mean)[-1])
  g6<-data.frame(aggregate(g2000b, list(rep(1:(nrow(g2000b) %/% n02 + 1), each = n02, len = nrow(g2000b))), mean)[-1])
  
  g1$sample<-as.character(nrow(g700)/(50*0.2))
  g2$sample<-as.character(nrow(g1000)/(100*0.2))
  g3$sample<-as.character(nrow(g2000)/(100*0.2))
  g4$sample<-as.character(nrow(g700b)/(50*0.2))
  g5$sample<-as.character(nrow(g1000b)/(100*0.2))
  g6$sample<-as.character(nrow(g2000b)/(100*0.2))
  names(g1)<-names(g2)<-names(g3)<-names(g4)<-names(g5)<-names(g6)<-c("DCP-QR","DCP-opt",
                                                                      "CQR","CQR-m","CQR-r","CP-reg","Sample")
  
  
  ######################################################################################################################
  ########################################BOXPLOT FOR COVERAGE#########################################################
  
  plotsinl<-function(d1,d2,d3,d4,d5,d6){
    dat<-rbind.data.frame(d1,d2,d3)
    dat2<-rbind.data.frame(d4,d5,d6)
    dat$Sample<-factor(dat$Sample,levels=c("700","1000","2000"))
    dat2$Sample<-factor(dat2$Sample,levels=c("700","1000","2000"))
    datm<-melt(dat,id.vars="Sample")
    dat2m<-melt(dat2,id.vars="Sample")
    
    p1<-ggplot(datm,aes(x=Sample,y=value,col=variable))+facet_wrap(~variable)+
      geom_boxplot()+
      stat_summary(fun = mean, geom = "point", col = "darkred") +  
      stat_summary(fun = mean, geom = "text", col = "red",     
                   vjust = 1.5, aes(label = paste(round(..y.., digits = 3))))+
      xlab("Sample Size")+ylab("Avg.Length")
    
    p2<-ggplot(dat2m,aes(x=Sample,y=value,col=variable))+facet_wrap(~variable)+
      geom_boxplot()+
      stat_summary(fun = mean, geom = "point", col = "darkred") +  
      stat_summary(fun = mean, geom = "text", col = "red",     
                   vjust = 1.5, aes(label = paste(round(..y.., digits = 3))))+
      xlab("Sample Size")+ylab("Avg.Length")
    
    k<-list(p1,p2)
    return(k)
    
  }
  plotsinl(g1,g2,g3,g4,g5,g6)[1]
  plotsinl(g1,g2,g3,g4,g5,g6)[2]
######################################################################################################
######STANDARD DEVIATION OF PREDICTED CONDITIONAL COVERAGE PROBABILITY#############################

  
  
  
  
  
  
  
  
  
#############################SINE WITH OUTLIERS###############################################################
  ###############################################################################################################
  data_out<-function(n,s){
    Xmat<-r<-r2<-yh<-yho<-res<-kk<-NULL
    for (i in 1:s){
      set.seed(i)
      x<-matrix(runif(n=n*100, min=0, max=1), nrow=n)
      xd<-data.frame(x)
      beta<-t(as.matrix(c(rep(1,5),rep(0,95))))
      xb<-as.vector(sqrt(1+(beta%*%t(x))^2))
      for (k in 1:n){
        res[k]<-rnorm(1,0,xb[k])
      }
      res2<-rnorm(n,0,1)
      t2<-as.matrix(res)
      t3<-as.matrix(res2)
      v<-as.matrix(rnorm(n,1,2))
      for (j in 1:n){
        if (j%%2==0){kk[j]=0.05}
        else {kk[j]=0.9}
      }
      kk<-as.matrix(kk)
      c<-floor((i+1)/2)
      cator<-ifelse(x[,c]>kk,1,0)
      cator<-as.matrix(cator)
      y<-t(2*sin(pi*beta%*%t(x)))+t((pi*beta%*%t(x)))+t2+100*(cator)  ##heteroskedastic 
      y2<-t(2*sin(pi*beta%*%t(x)))+t((pi*beta%*%t(x)))+t3+100*(cator) ##homoskedastic
      y<-as.matrix(as.vector(y))
      y2<-as.matrix(as.vector(y2))
      Xmat<-rbind(Xmat,xd)
      r<-rbind(r,t2)
      r2<-rbind(r2,t3)
      yh<-rbind(yh,y)
      yho<-rbind(yho,y2)
    }
    l<-list(Xmat,yh,yho,t2,t3,kk)
    return(l)
  }
  
  a<-data_out(700,200) ##get 50 iterations only
  b<-data_out(1000,200)  ##100
  c<-data_out(2000,200) ##100
  d<-data_out(5000,200)  ##50
  
  
  ##a
  xa<-data.frame(a[1])   ###140000 100
  y1a<-data.frame(a[2])  ##140000 1
  y2a<-data.frame(a[3])  ##140000 1
  names(y1a)<-names(y2a)<-names(xa)<-NULL   ##y1[1,]
  
  T.ho <- floor(700*0.2)
  sin700_o<-cpo(50,700,xa,y1a,y2a)
  
  
  ##b
  xb<-data.frame(b[1])   ###200000 100
  y1b<-data.frame(b[2])  ##200000 1
  y2b<-data.frame(b[3])  ##200000 1
  names(y1b)<-names(y2b)<-names(xb)<-NULL   ##y1[1,]
  
  T.ho <- floor(1000*0.2)
  sin1000_o<-cpo(100,1000,xb,y1b,y2b)
  
  ##c
  x<-data.frame(c[1])   ###200000 100
  y1<-data.frame(c[2])  ##200000 1
  y2<-data.frame(c[3])  ##200000 1
  names(y1)<-names(y2)<-names(x)<-NULL   ##y1[1,]
  
  T.ho <- floor(2000*0.2)
  sin2000_o<-cpo(100,2000,x,y1,y2)    ###Time difference of 1.115348 hours
  
  ##Stop at 2000 and now get results 
  
  write.xlsx(sin700_o,"sin700_o.xlsx")
  write.xlsx(sin1000_o,"sin1000_o.xlsx")
  write.xlsx(sin2000_o,"sin2000_o.xlsx")
  
  ###import data and do computations
  sin700_ds<-lapply(1:6, function(i) read_excel("sin700b.xlsx", sheet = i,col_names=FALSE))
  sin1000_ds<-lapply(1:6, function(i) read_excel("sin1000.xlsx", sheet = i,col_names=FALSE))
  sin2000_ds<-lapply(1:6, function(i) read_excel("sin2000b.xlsx", sheet = i,col_names=FALSE))
  
  ###Create table for avg coverage,avg length and sd of conditional coverage probability 
  a<-cbind.data.frame(sin700_ds[1],sin700_ds[2],sin700_ds[3],sin700_ds[4],sin700_ds[5],sin700_ds[6])
  b<-cbind.data.frame(sin1000_ds[1],sin1000_ds[2],sin1000_ds[3],sin1000_ds[4],sin1000_ds[5],sin1000_ds[6])
  c<-cbind.data.frame(sin2000_ds[1],sin2000_ds[2],sin2000_ds[3],sin2000_ds[4],sin2000_ds[5],sin2000_ds[6])
  
  d<-rbind.data.frame(a,b,c)
  names(d)<-c("Cov-1","Cov-2","Leng-1","Leng-2","SD-1","SD-2")
  xtable(d,digits=3)
  
################################################################################################
##############################################################################################
  
###import data and do computations
sin700_ods<-lapply(1:6, function(i) read_excel("sin700_o.xlsx", sheet = i,col_names=FALSE))
sin1000_ods<-lapply(1:6, function(i) read_excel("sin1000_o.xlsx", sheet = i,col_names=FALSE))
sin2000_ods<-lapply(1:6, function(i) read_excel("sin2000_o.xlsx", sheet = i,col_names=FALSE))

###Create table for avg coverage,avg length and sd of conditional coverage probability 
a<-cbind.data.frame(sin700_ods[1],sin700_ods[2],sin700_ods[3],sin700_ods[4],sin700_ods[5],sin700_ods[6])
b<-cbind.data.frame(sin1000_ods[1],sin1000_ods[2],sin1000_ods[3],sin1000_ods[4],sin1000_ods[5],sin1000_ods[6])
c<-cbind.data.frame(sin2000_ods[1],sin2000_ods[2],sin2000_ods[3],sin2000_ods[4],sin2000_ods[5],sin2000_ods[6])
  
d<-rbind.data.frame(a,b,c)
names(d)<-c("Cov-1","Cov-2","Leng-1","Leng-2","SD-1","SD-2")
xtable(d,digits=3)

#####################################################################################################
#########BOXPLOTS FOR AVERAGE COVERAGE AND LENGTH#############################################################
s700<-lapply(7:8, function(i) read_excel("sin700_o.xlsx", sheet = i))
s1000<-lapply(7:8, function(i) read_excel("sin1000_o.xlsx", sheet = i))
s2000<-lapply(7:8, function(i) read_excel("sin2000_o.xlsx", sheet = i))

###Heteroskedastic case
g700<-data.frame(s700[1])
g1000<-data.frame(s1000[1])
g2000<-data.frame(s2000[1])

n0<-nrow(g700)/50
n01<-nrow(g1000)/100
n02<-nrow(g2000)/100

g1<-data.frame(aggregate(g700, list(rep(1:(nrow(g700) %/% n0 + 1), each = n0, len = nrow(g700))), mean)[-1])
g2<-data.frame(aggregate(g1000, list(rep(1:(nrow(g1000) %/% n01 + 1), each = n01, len = nrow(g1000))), mean)[-1])
g3<-data.frame(aggregate(g2000, list(rep(1:(nrow(g2000) %/% n02 + 1), each = n02, len = nrow(g2000))), mean)[-1])

###Homoskedastic case
g700b<-data.frame(s700[2])
g1000b<-data.frame(s1000[2])
g2000b<-data.frame(s2000[2])


g4<-data.frame(aggregate(g700b, list(rep(1:(nrow(g700b) %/% n0 + 1), each = n0, len = nrow(g700b))), mean)[-1])
g5<-data.frame(aggregate(g1000b, list(rep(1:(nrow(g1000b) %/% n01 + 1), each = n01, len = nrow(g1000b))), mean)[-1])
g6<-data.frame(aggregate(g2000b, list(rep(1:(nrow(g2000b) %/% n02 + 1), each = n02, len = nrow(g2000b))), mean)[-1])

g1$sample<-as.character(nrow(g700)/(50*0.2))
g2$sample<-as.character(nrow(g1000)/(100*0.2))
g3$sample<-as.character(nrow(g2000)/(100*0.2))
g4$sample<-as.character(nrow(g700b)/(50*0.2))
g5$sample<-as.character(nrow(g1000b)/(100*0.2))
g6$sample<-as.character(nrow(g2000b)/(100*0.2))
names(g1)<-names(g2)<-names(g3)<-names(g4)<-names(g5)<-names(g6)<-c("DCP-QR","DCP-opt",
                                                                    "CQR","CQR-m","CQR-r","CP-reg","Sample")


######################################################################################################################
########################################BOXPLOT FOR COVERAGE#########################################################

plotsinc<-function(d1,d2,d3,d4,d5,d6){
  dat<-rbind.data.frame(d1,d2,d3)
  dat2<-rbind.data.frame(d4,d5,d6)
  dat$Sample<-factor(dat$Sample,levels=c("700","1000","2000"))
  dat2$Sample<-factor(dat2$Sample,levels=c("700","1000","2000"))
  datm<-melt(dat,id.vars="Sample")
  dat2m<-melt(dat2,id.vars="Sample")
  
  p1<-ggplot(datm,aes(x=Sample,y=value,col=variable))+facet_wrap(~variable)+
    geom_boxplot()+
    geom_hline(yintercept=0.9,col="darkgrey",linewidth=1)+
    stat_summary(fun = mean, geom = "point", col = "darkred") +  
    stat_summary(fun = mean, geom = "text", col = "red",     
                 vjust = 1.5, aes(label = paste(round(..y.., digits = 3))))+
    xlab("Sample Size")+ylab("Avg.Coverage")
  
  p2<-ggplot(dat2m,aes(x=Sample,y=value,col=variable))+facet_wrap(~variable)+
    geom_boxplot()+
    geom_hline(yintercept=0.9,col="darkgrey",linewidth=1)+
    stat_summary(fun = mean, geom = "point", col = "darkred") +  
    stat_summary(fun = mean, geom = "text", col = "red",     
                 vjust = 1.5, aes(label = paste(round(..y.., digits = 3))))+
    xlab("Sample Size")+ylab("Avg.Coverage")
  
  k<-list(p1,p2)
  return(k)
  
}
plotsinc(g1,g2,g3,g4,g5,g6)[1]  
plotsinc(g1,g2,g3,g4,g5,g6)[2]  
  
####################################################################################
##########################AVERAGE LENGTH #############################################
s700l<-lapply(9:10, function(i) read_excel("sin700_o.xlsx", sheet = i))
s1000l<-lapply(9:10, function(i) read_excel("sin1000_o.xlsx", sheet = i))
s2000l<-lapply(9:10, function(i) read_excel("sin2000_o.xlsx", sheet = i))

###Heteroskedastic case
g700<-data.frame(s700l[1])
g1000<-data.frame(s1000l[1])
g2000<-data.frame(s2000l[1])

n0<-nrow(g700)/50
n01<-nrow(g1000)/100
n02<-nrow(g2000)/100

g1<-data.frame(aggregate(g700, list(rep(1:(nrow(g700) %/% n0 + 1), each = n0, len = nrow(g700))), mean)[-1])
g2<-data.frame(aggregate(g1000, list(rep(1:(nrow(g1000) %/% n01 + 1), each = n01, len = nrow(g1000))), mean)[-1])
g3<-data.frame(aggregate(g2000, list(rep(1:(nrow(g2000) %/% n02 + 1), each = n02, len = nrow(g2000))), mean)[-1])

###Homoskedastic case
g700b<-data.frame(s700l[2])
g1000b<-data.frame(s1000l[2])
g2000b<-data.frame(s2000l[2])


g4<-data.frame(aggregate(g700b, list(rep(1:(nrow(g700b) %/% n0 + 1), each = n0, len = nrow(g700b))), mean)[-1])
g5<-data.frame(aggregate(g1000b, list(rep(1:(nrow(g1000b) %/% n01 + 1), each = n01, len = nrow(g1000b))), mean)[-1])
g6<-data.frame(aggregate(g2000b, list(rep(1:(nrow(g2000b) %/% n02 + 1), each = n02, len = nrow(g2000b))), mean)[-1])

g1$sample<-as.character(nrow(g700)/(50*0.2))
g2$sample<-as.character(nrow(g1000)/(100*0.2))
g3$sample<-as.character(nrow(g2000)/(100*0.2))
g4$sample<-as.character(nrow(g700b)/(50*0.2))
g5$sample<-as.character(nrow(g1000b)/(100*0.2))
g6$sample<-as.character(nrow(g2000b)/(100*0.2))
names(g1)<-names(g2)<-names(g3)<-names(g4)<-names(g5)<-names(g6)<-c("DCP-QR","DCP-opt",
                                                                    "CQR","CQR-m","CQR-r","CP-reg","Sample")


######################################################################################################################
########################################BOXPLOT FOR COVERAGE#########################################################

plotsinl<-function(d1,d2,d3,d4,d5,d6){
  dat<-rbind.data.frame(d1,d2,d3)
  dat2<-rbind.data.frame(d4,d5,d6)
  dat$Sample<-factor(dat$Sample,levels=c("700","1000","2000"))
  dat2$Sample<-factor(dat2$Sample,levels=c("700","1000","2000"))
  datm<-melt(dat,id.vars="Sample")
  dat2m<-melt(dat2,id.vars="Sample")
  
  p1<-ggplot(datm,aes(x=Sample,y=value,col=variable))+facet_wrap(~variable)+
    geom_boxplot()+
    stat_summary(fun = mean, geom = "point", col = "red") +  
    stat_summary(fun = mean, geom = "text", col = "darkred",     
                 vjust = 1.65, aes(label = paste(round(..y.., digits = 3))))+
    xlab("Sample Size")+ylab("Avg.Length")
  
  p2<-ggplot(dat2m,aes(x=Sample,y=value,col=variable))+facet_wrap(~variable)+
    geom_boxplot()+
    stat_summary(fun = mean, geom = "point", col = "red") +  
    stat_summary(fun = mean, geom = "text", col = "darkred",     
                 vjust = -1.5, aes(label = paste(round(..y.., digits = 3))))+
    xlab("Sample Size")+ylab("Avg.Length")
  
  k<-list(p1,p2)
  return(k)
  
}
plotsinl(g1,g2,g3,g4,g5,g6)[1]  
plotsinl(g1,g2,g3,g4,g5,g6)[2]  













































##Export data to excel and then data to python 
write.xlsx(data.frame(sin700[11]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X0.xlsx")
write.xlsx(data.frame(sin700[12]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y0.xlsx")
write.xlsx(data.frame(sin700[13]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X1.xlsx")
write.xlsx(data.frame(sin700[14]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y1.xlsx")
write.xlsx(data.frame(sin700[15]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Xt.xlsx")
write.xlsx(data.frame(sin700[16]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Yt.xlsx")
write.xlsx(data.frame(sin700[17]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X02.xlsx")
write.xlsx(data.frame(sin700[18]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y02.xlsx")
write.xlsx(data.frame(sin700[19]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X12.xlsx")
write.xlsx(data.frame(sin700[20]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y12.xlsx")
write.xlsx(data.frame(sin700[21]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Xt2.xlsx")
write.xlsx(data.frame(sin700[22]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Yt2.xlsx")


###Now I need to get empirical coverage and grand avg from cov.mat,cov.mat2, go to python as well 


###########EMPIRICAL COVERAGE FOR CONDITIONAL VALIDITY#########################
    covm1<-data.frame(g[7])
    covm2<-data.frame(g[8])
    case1<-empcov(n,covm1,it)  ##matrix array 
    case2<-empcov(n,covm2,it) 
    gavg<-data.frame(colMeans(case1,na.rm=TRUE))
    names(gavg)<-c("gavg")
    gavg2<-data.frame(colMeans(case2,na.rm=TRUE))
    names(gavg2)<-c("gavg2")
    both<-cbind.data.frame(gavg,gavg2)
    both<-round(both,4)
    
    w<-list(j,case1,case2,both)
    return(w)
}
sine1000_res<-res(sin1000,200,50)
##take results to Latex
data.frame(sine1000_res[1])
data.frame(sine1000_res[4])

###Exporting to excel
write.xlsx(data.frame(sin1000[11]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X0_000.xlsx")
write.xlsx(data.frame(sin1000[12]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y0_000.xlsx")
write.xlsx(data.frame(sin1000[13]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X1_000.xlsx")
write.xlsx(data.frame(sin1000[14]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y1_000.xlsx")
write.xlsx(data.frame(sin1000[15]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Xt_000.xlsx")
write.xlsx(data.frame(sin1000[16]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Yt_000.xlsx")
write.xlsx(data.frame(sin1000[17]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X02_000.xlsx")
write.xlsx(data.frame(sin1000[18]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y02_000.xlsx")
write.xlsx(data.frame(sin1000[19]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X12_000.xlsx")
write.xlsx(data.frame(sin1000[20]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y12_000.xlsx")
write.xlsx(data.frame(sin1000[21]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Xt2_000.xlsx")
write.xlsx(data.frame(sin1000[22]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Yt2_000.xlsx")



##Case 1:Heteroskedasticity, Case 2: Homoskedasticity
case10001<-empcov(200,covm1,50)  ##matrix array 
case2<-empcov(700,covm2,50) 

##Visualization-need more algorithms from python-rf,nn,cqr-rf,cqr-nn,...
boxplot(case1[,2],case1[,3],case1[,4],case1[,5],case2[,6],names=c("DCP-opt","CQR","CQR-m","CQR-r","CP-reg"))
abline(h=0.9,col="red",lwd=2)  ##Except for CP-reg, median number of computed averages show undercoverage-most of CQR-m
###show undercoverage. Variability is lower for CP-reg and it shows least undercoverage. Most obs lie in a small interval around 
##0.9
##However sample size is small. Angeloupoulos suggest at least 1000 calibration obs. 

boxplot(case2[,2],case2[,3],case1[,4],case1[,5],case2[,6],names=c("DCP-opt","CQR","CQR-m","CQR-r","CP-reg"))
abline(h=0.9,col="red",lwd=2)

###read in data from excel:
###Sample size:1000,100 iterations 
s1000<-lapply(7:8, function(i) read_excel("sin1000.xlsx", sheet = i))
het1000<-data.frame(s1000[1])  ###20000 6
hom1000<-data.frame(s1000[2])  ###20000 6



####PER ITERATION AVERAGE COVERAGE PLOT
s700<-lapply(7:8, function(i) read_excel("sin700b.xlsx", sheet = i))
s1000<-lapply(7:8, function(i) read_excel("sin1000.xlsx", sheet = i))
s2000<-lapply(7:8, function(i) read_excel("sin2000b.xlsx", sheet = i))


  g700<-data.frame(s700[1])
  g1000<-data.frame(s1000[1])
  g2000<-data.frame(s2000[1])
  
  n0<-nrow(g700)/50
  n01<-nrow(g1000)/100
  n02<-nrow(g2000)/100
  n03<-nrow(gd)/1000
  n04<-nrow(ge)/1000
  n05<-nrow(gf)/1000
  g1<-data.frame(aggregate(g, list(rep(1:(nrow(g) %/% n0 + 1), each = n0, len = nrow(g))), mean)[-1])
  g2<-data.frame(aggregate(gb, list(rep(1:(nrow(gb) %/% n01 + 1), each = n01, len = nrow(gb))), mean)[-1])
  g3<-data.frame(aggregate(gc, list(rep(1:(nrow(gc) %/% n02 + 1), each = n02, len = nrow(gc))), mean)[-1])
  g4<-data.frame(aggregate(gd, list(rep(1:(nrow(gd) %/% n03 + 1), each = n03, len = nrow(gd))), mean)[-1])
  g5<-data.frame(aggregate(ge, list(rep(1:(nrow(ge) %/% n04 + 1), each = n04, len = nrow(ge))), mean)[-1])
  g6<-data.frame(aggregate(gf, list(rep(1:(nrow(gf) %/% n05 + 1), each = n05, len = nrow(gf))), mean)[-1])
  g1$type<-g2$type<-g3$type<-g4$type<-g5$type<-g6$type<-i
  g1$sample<-as.character(nrow(g)/(1000*0.2))
  g2$sample<-as.character(nrow(gb)/(1000*0.2))
  g3$sample<-as.character(nrow(gc)/(1000*0.2))
  g4$sample<-as.character(nrow(gd)/(1000*0.2))
  g5$sample<-as.character(nrow(ge)/(1000*0.2))
  g6$sample<-as.character(nrow(gf)/(1000*0.2))
  names(g1)<-names(g2)<-names(g3)<-names(g4)<-names(g5)<-names(g6)<-c("DCP-QR","DCP-opt",
                                                                      "CQR","CQR-m","CQR-r","CP-reg","i","Sample")
  assign(paste("el50",i,sep=""),g1)
  assign(paste("el70",i,sep=""),g2)
  assign(paste("el100",i,sep=""),g3)
  assign(paste("el200",i,sep=""),g4)
  assign(paste("el300",i,sep=""),g5)
  assign(paste("elfive",i,sep=""),g6)
}


##########2000-200
##Data generation
##############


sin2000_res<-res(sin2000,400,100)
sin2000_res[1]

write.xlsx(data.frame(sin2000[11]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X0_2000.xlsx")
write.xlsx(data.frame(sin2000[12]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y0_2000.xlsx")
write.xlsx(data.frame(sin2000[13]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X1_2000.xlsx")
write.xlsx(data.frame(sin2000[14]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y1_2000.xlsx")
write.xlsx(data.frame(sin2000[15]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Xt_2000.xlsx")
write.xlsx(data.frame(sin2000[16]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Yt_2000.xlsx")
write.xlsx(data.frame(sin2000[17]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X02_2000.xlsx")
write.xlsx(data.frame(sin2000[18]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y02_2000.xlsx")
write.xlsx(data.frame(sin2000[19]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X12_2000.xlsx")
write.xlsx(data.frame(sin2000[20]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y12_2000.xlsx")
write.xlsx(data.frame(sin2000[21]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Xt2_2000.xlsx")
write.xlsx(data.frame(sin2000[22]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Yt2_2000.xlsx")

dat<-cbind.data.frame(as.matrix(cov1),as.matrix(cov2),as.matrix(leng1),as.matrix(leng2),as.matrix(res.cov.cond),as.matrix(res.cov.cond2))
names(dat)<-c("cov_het","cov_hom","leng_het","leng_hom","Sd_het","Sd_hom")
formatn<-function(x){
  #if (x>0.09){y=x}
  #else {y=formatC(x,format="e",digits=2)}
  y<-ifelse(abs(x)>=0.01,round(x,4),formatC(x,format="e",digits=2))
  return(y)
}
dat[]<-lapply(dat,formatn)
xtable(dat)

##Use separate script for generating data for exporting
####5000-50
##############
a<-data(5000,200) 
x<-data.frame(a[1])
y1<-data.frame(a[2])
y2<-data.frame(a[3]) 
names(y1)<-names(y2)<-names(x)<-NULL 

alpha.sig<-0.1
T.ho <- floor(5000*0.2) ##20% holdout sample  ##try 3/4-1/4 split as well 

sin5000<-cpo(50,5000,x,y1,y2)

sin5000_res<-res(sin5000,1000,50)
sin5000_res[1]

write.xlsx(data.frame(sin5000[11]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X0_5000.xlsx")
write.xlsx(data.frame(sin5000[12]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y0_5000.xlsx")
write.xlsx(data.frame(sin5000[13]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X1_5000.xlsx")
write.xlsx(data.frame(sin5000[14]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y1_5000.xlsx")
write.xlsx(data.frame(sin5000[15]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Xt_5000.xlsx")
write.xlsx(data.frame(sin5000[16]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Yt_5000.xlsx")
write.xlsx(data.frame(sin5000[17]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X02_5000.xlsx")
write.xlsx(data.frame(sin5000[18]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y02_5000.xlsx")
write.xlsx(data.frame(sin5000[19]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X12_5000.xlsx")
write.xlsx(data.frame(sin5000[20]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y12_5000.xlsx")
write.xlsx(data.frame(sin5000[21]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Xt2_5000.xlsx")
write.xlsx(data.frame(sin5000[22]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Yt2_5000.xlsx")


############smaller samples-50,70




dat2<-cbind.data.frame(as.matrix(cov),as.matrix(covh),as.matrix(leng),as.matrix(lengh),as.matrix(res.cov.cond.h),as.matrix(res.cov.cond2.het))
names(dat2)<-c("cov","covh","leng","lengh","Sd_cond_het","Sd_cov_cond2_het")
formatn<-function(x){
  #if (x>0.09){y=x}
  #else {y=formatC(x,format="e",digits=2)}
  y<-ifelse(abs(x)>=0.01,round(x,4),formatC(x,format="e",digits=2))
  return(y)
}
dat2[]<-lapply(dat2,formatn)
xtable(dat)

#library(writexl)
#write_xlsx(dat,"res_2000_100_sin.xlsx")
#write_xlsx(dat2,"res_5000_50_sin.xlsx")

#write_xlsx(data.frame(sx0),"sx0.xlsx")
#write_xlsx(data.frame(sy0),"sy0.xlsx")
#write_xlsx(data.frame(stx),"stx.xlsx")
#write_xlsx(data.frame(sty),"sty.xlsx")
#write_xlsx(data.frame(sx1),"sx1.xlsx")
#write_xlsx(data.frame(sy1),"sy1.xlsx")

#write_xlsx(data.frame(sx02),"sx02.xlsx")
#write_xlsx(data.frame(sy02),"sy02.xlsx")
#write_xlsx(data.frame(sx12),"sx12.xlsx")
#write_xlsx(data.frame(sy12),"sy12.xlsx")
#write_xlsx(data.frame(stx2),"stx2.xlsx")
#write_xlsx(data.frame(sty2),"sty2.xlsx")

##empirical coverage
empcov<-function(n,m,s){
  sv<-NULL
  for (i in 1:s){
    subval<-m[n*(i-1)+1:n*i,]
    subval_avg<-colMeans(subval,na.rm=TRUE)
    sv<-rbind(sv,subval_avg)
  }
  return(sv)
}
h<-data.frame(cov.mat)
w<-empcov(1000,h,50) ##2000 is training-validation size ,1000 for testing
write_xlsx(data.frame(w),"w_5000_50.xlsx")

gavg<-data.frame(colMeans(w,na.rm=TRUE))
names(gavg)<-c("gavg")
write_xlsx(data.frame(gavg),"gavg_5000_50.xlsx")

######################################################################
##create outliers
data_out<-function(n,s){
  Xmat<-r<-r2<-yh<-yho<-res<-kk<-NULL
  for (i in 1:s){
    set.seed(i)
    x<-matrix(runif(n=n*100, min=0, max=1), nrow=n)
    xd<-data.frame(x)
    beta<-t(as.matrix(c(rep(1,5),rep(0,95))))
    xb<-as.vector(sqrt(1+(beta%*%t(x))^2))
    for (k in 1:n){
      res[k]<-rnorm(1,0,xb[k])
    }
    res2<-rnorm(n,0,1)
    t2<-as.matrix(res)
    t3<-as.matrix(res2)
    v<-as.matrix(rnorm(n,1,2))
    for (j in 1:n){
      if (j%%2==0){kk[j]=0.05}
      else {kk[j]=0.9}
    }
    kk<-as.matrix(kk)
    c<-floor((i+1)/2)
    cator<-ifelse(x[,c]>kk,1,0)
    cator<-as.matrix(cator)
    y<-t(2*sin(pi*beta%*%t(x)))+t((pi*beta%*%t(x)))+t2+100*(cator)  ##heteroskedastic 
    y2<-t(2*sin(pi*beta%*%t(x)))+t((pi*beta%*%t(x)))+t3+100*(cator) ##homoskedastic
    y<-as.matrix(as.vector(y))
    y2<-as.matrix(as.vector(y2))
    Xmat<-rbind(Xmat,xd)
    r<-rbind(r,t2)
    r2<-rbind(r2,t3)
    yh<-rbind(yh,y)
    yho<-rbind(yho,y2)
  }
  l<-list(Xmat,yh,yho,t2,t3,kk)
  return(l)
}

a<-data_out(700,200) ##get 50 iterations only
b<-data_out(1000,200)  ##100
c<-data_out(2000,200) ##100
d<-data_out(5000,200)  ##50


##a
xa<-data.frame(a[1])   ###140000 100
y1a<-data.frame(a[2])  ##140000 1
y2a<-data.frame(a[3])  ##140000 1
names(y1a)<-names(y2a)<-names(xa)<-NULL   ##y1[1,]

T.ho <- floor(700*0.2)
sin700_o<-cpo(50,700,xa,y1a,y2a)


##b
xb<-data.frame(b[1])   ###200000 100
y1b<-data.frame(b[2])  ##200000 1
y2b<-data.frame(b[3])  ##200000 1
names(y1b)<-names(y2b)<-names(xb)<-NULL   ##y1[1,]

T.ho <- floor(1000*0.2)
sin1000_o<-cpo(100,1000,xb,y1b,y2b)

##c
x<-data.frame(c[1])   ###200000 100
y1<-data.frame(c[2])  ##200000 1
y2<-data.frame(c[3])  ##200000 1
names(y1)<-names(y2)<-names(x)<-NULL   ##y1[1,]

T.ho <- floor(2000*0.2)
sin2000_o<-cpo(100,2000,x,y1,y2)    ###Time difference of 1.115348 hours

##Stop at 2000 and now get results 

write.xlsx(sin700_o,"sin700_o.xlsx")
write.xlsx(sin1000_o,"sin1000_o.xlsx")
write.xlsx(sin2000_o,"sin2000_o.xlsx")

h<-list(cov1,cov2,leng1,leng2,res.cov.cond,res.cov.cond2,
        cov.mat,cov.mat2,leng.mat,leng.mat2,
        sx0,sy0,sx1,sy1,stx,sty,
        sx02,sy02,sx12,sy12,stx2,sty2,j)

func<-leng(x,y,z){  ##not working
  a<-cbind.data.frame(x,y,z)
  names(a)<-c("700","1000","2000")
  a<-as.matrix(a)
  a<-data.frame(t(a))
  a1<-data.frame(x=seq_along(a[,1]),a)
  a1<-melt(a1,id.vars="x")
  a1$sample<-rep(c(700,1000,2000),3)
  p<-ggplot(a1, aes(x = sample, y = value, color = variable)) +
    geom_line(size=1.25)+geom_point()+
    xlab("Sample Size") +
    ylab("Avg.Length")+
    labs(xlab="Sample Size",ylab="Avg.Length")+
    theme(legend.title = element_blank(),
          legend.position = "bottom") +
    scale_colour_discrete(labels = c('DCP-QR','DCP-opt','CQR','CQR-m','CQR-r','CP-reg'))
  return(p)
}
coverhet<-cbind.data.frame(sin700_o[1],sin1000_o[1],sin2000_o[1])
names(coverhet)<-c("Cov_700","Cov_1000","Cov_2000")
coverhet<-as.matrix(coverhet)
coverhet<-data.frame(t(coverhet))
#xax<-c(700,1000,2000)
#coverhom<-cbind.data.frame(coverhom,xax)
names(coverhet)<-c("X1","X2","X3","X4","X5","X6")


coverhom<-cbind.data.frame(sin700_o[2],sin1000_o[2],sin2000_o[2])
names(coverhom)<-c("Cov_700h","Cov_1000h","Cov_2000h")
coverhom<-as.matrix(coverhom)
coverhom<-data.frame(t(coverhom))
#xax<-c(700,1000,2000)
#coverhom<-cbind.data.frame(coverhom,xax)
names(coverhom)<-c("X1","X2","X3","X4","X5","X6")

library(reshape)
c1<- data.frame(x = seq_along(coverhom[, 1]),coverhom)
c1 <- melt(c1, id.vars = "x")
c1$sample<-rep(c(700,1000,2000),3)

p1<-ggplot(c1, aes(x = sample, y = value, color = variable)) +
  geom_line(size=1.25)+geom_point()+
  xlab("Sample Size") +
  ylab("Avg.Coverage")+
 labs(xlab="Sample Size",ylab="Avg.Coverage")+
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  scale_colour_discrete(labels = c('DCP-QR','DCP-opt','CQR','CQR-m','CQR-r','CP-reg'))

c2<- data.frame(x = seq_along(coverhet[, 1]),coverhet)
c2 <- melt(c2, id.vars = "x")
c2$sample<-rep(c(700,1000,2000),3)

########################################################################
p2<-ggplot(c2, aes(x = sample, y = value, color = variable)) +
  geom_line(size=1.25)+geom_point()+
  xlab("Sample Size") +
  ylab("Avg.Coverage")+
  labs(xlab="Sample Size",ylab="Avg.Coverage")+
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  scale_colour_discrete(labels = c('DCP-QR','DCP-opt','CQR','CQR-m','CQR-r','CP-reg'))
  =
https://stackoverflow.com/questions/1249548/side-by-side-plots-with-ggplot2
##################################################################################


install.packages("cowplot")
library(cowplot)
plot_grid(p1,p2, labels = "AUTO")

install.packages("gridExtra")

##homoskedasticity
lenghom<-cbind.data.frame(sin700_o[4],sin1000_o[4],sin2000_o[4])
names(lenghom)<-c("700h","1000h","2000h")
lenghom<-as.matrix(lenghom)
lenghom<-data.frame(t(lenghom))
#xax<-c(700,1000,2000)
#coverhom<-cbind.data.frame(coverhom,xax)
names(lenghom)<-c("X1","X2","X3","X4","X5","X6")

library(reshape)
c1<- data.frame(x = seq_along(lenghom[, 1]),lenghom)
c1 <- melt(c1, id.vars = "x")
c1$sample<-rep(c(700,1000,2000),3)

p1<-ggplot(c1, aes(x = sample, y = value, color = variable)) +
  geom_line(size=1.25)+geom_point()+
  xlab("Sample Size") +
  ylab("Avg.Length")+
  labs(xlab="Sample Size",ylab="Avg.Length")+
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  scale_colour_discrete(labels = c('DCP-QR','DCP-opt','CQR','CQR-m','CQR-r','CP-reg'))

##heteroskedasticity
lenghet<-cbind.data.frame(sin700_o[3],sin1000_o[3],sin2000_o[3])
names(lenghet)<-c("700h","1000h","2000h")
lenghet<-as.matrix(lenghet)
lenghet<-data.frame(t(lenghet))
#xax<-c(700,1000,2000)
#coverhom<-cbind.data.frame(coverhom,xax)
names(lenghet)<-c("X1","X2","X3","X4","X5","X6")

library(reshape)
c1<- data.frame(x = seq_along(lenghet[, 1]),lenghet)
c1 <- melt(c1, id.vars = "x")
c1$sample<-rep(c(700,1000,2000),3)

p2<-ggplot(c1, aes(x = sample, y = value, color = variable)) +
  geom_line(size=1.25)+geom_point()+
  xlab("Sample Size") +
  ylab("Avg.Length")+
  labs(xlab="Sample Size",ylab="Avg.Length")+
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  scale_colour_discrete(labels = c('DCP-QR','DCP-opt','CQR','CQR-m','CQR-r','CP-reg'))

#################################################################################################
##############################################################################################
library(writexl)
###exporting for python
write_xlsx(data.frame(sin700_o[11]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X0_700_o.xlsx")
write_xlsx(data.frame(sin700_o[12]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y0_700_o.xlsx")
write_xlsx(data.frame(sin700_o[13]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X1_700_o.xlsx")
write_xlsx(data.frame(sin700_o[14]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y1_700_o.xlsx")
write_xlsx(data.frame(sin700_o[15]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Xt_700_o.xlsx")
write_xlsx(data.frame(sin700_o[16]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Yt_700_o.xlsx")
write_xlsx(data.frame(sin700_o[17]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X02_700_o.xlsx")
write_xlsx(data.frame(sin700_o[18]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y02_700_o.xlsx")
write_xlsx(data.frame(sin700_o[19]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X12_700_o.xlsx")
write_xlsx(data.frame(sin700_o[20]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y12_700_o.xlsx")
write_xlsx(data.frame(sin700_o[21]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Xt2_700_o.xlsx")
write_xlsx(data.frame(sin700_o[22]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Yt2_700_o.xlsx")

write_xlsx(data.frame(sin1000_o[11]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X0_1000_o.xlsx")
write_xlsx(data.frame(sin1000_o[12]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y0_1000_o.xlsx")
write_xlsx(data.frame(sin1000_o[13]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X1_1000_o.xlsx")
write_xlsx(data.frame(sin1000_o[14]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y1_1000_o.xlsx")
write_xlsx(data.frame(sin1000_o[15]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Xt_1000_o.xlsx")
write_xlsx(data.frame(sin1000_o[16]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Yt_1000_o.xlsx")
write_xlsx(data.frame(sin1000_o[17]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X02_1000_o.xlsx")
write_xlsx(data.frame(sin1000_o[18]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y02_1000_o.xlsx")
write_xlsx(data.frame(sin1000_o[19]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X12_1000_o.xlsx")
write_xlsx(data.frame(sin1000_o[20]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y12_1000_o.xlsx")
write_xlsx(data.frame(sin1000_o[21]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Xt2_1000_o.xlsx")
write_xlsx(data.frame(sin1000_o[22]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Yt2_1000_o.xlsx")

write_xlsx(data.frame(sin2000_o[11]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X0_2000_o.xlsx")
write_xlsx(data.frame(sin2000_o[12]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y0_2000_o.xlsx")
write_xlsx(data.frame(sin2000_o[13]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X1_2000_o.xlsx")
write_xlsx(data.frame(sin2000_o[14]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y1_2000_o.xlsx")
write_xlsx(data.frame(sin2000_o[15]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Xt_2000_o.xlsx")
write_xlsx(data.frame(sin2000_o[16]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Yt_2000_o.xlsx")
write_xlsx(data.frame(sin2000_o[17]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X02_2000_o.xlsx")
write_xlsx(data.frame(sin2000_o[18]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y02_2000_o.xlsx")
write_xlsx(data.frame(sin2000_o[19]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/X12_2000_o.xlsx")
write_xlsx(data.frame(sin2000_o[20]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Y12_2000_o.xlsx")
write_xlsx(data.frame(sin2000_o[21]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Xt2_2000_o.xlsx")
write_xlsx(data.frame(sin2000_o[22]),"C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/Yt2_2000_o.xlsx")







###Avg.coverage and length by sample size and 
##boxplot of avg coverage and length for outliers

###700
##heteroskedastic case
cov700het<-data.frame(sin700_o[8])
##get per iteration average coverage-s is no of iterations,m is the matrix used,n is test set size
empcov<-function(n,m,s){
  sv<-NULL
  for (i in 1:s){
    subval<-m[n*(i-1)+1:n*i,]
    subval_avg<-colMeans(subval,na.rm=TRUE)
    sv<-rbind(sv,subval_avg)
  }
  return(sv)
}

cov_it_700_het<-empcov(140,cov700het,50)
##loading cp-rf and cp-lw from Python by Excel 
library(openxlsx)
cov_het_700_cps<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s700_cov_cps.xlsx")
cov_het_700_cps<-cov_het_700_cps[,-1]
cov_het_700_lw<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s700_cov_lwo.xlsx")
cov_het_700_lw<-cov_het_700_lw[,-1]

cps700het<-data.frame(apply(cov_het_700_cps,1,function(x){sum(x)/140}))
lw700het<-data.frame(apply(cov_het_700_lw,1,function(x){sum(x)/140}))

cov700het<-cbind.data.frame(cov_it_700_het_truncate,cps700het,lw700het)
names(cov700het)<-c("DCP-opt","CQR","CQR-m","CQR-r","CP-reg","CP-RF","CP-LW")

library(ggplot2)
cov700het<- melt(cov700het)

theme_set(theme_gray(base_size = 20))
p<- ggplot(data = cov700het, aes(x = variable, y = value))
p + geom_boxplot() + labs(x = "Method", y = "Average Coverage/iteration")+geom_hline(yintercept=0.9,
                                                                             col="red",size=0.75)

###1000 iterations
cov1000het<-data.frame(sin1000_o[8])
cov_it_1000_het<-empcov(200,cov1000het,100) ##100 6
#cov_it_1000_het<-na.omit(cov_it_1000_het) ##100 201
cov_het_1000_cps<-read.xlsx("C:/Users/RimJhim/Documents/s1000_cov_cps.xlsx")
cov_het_1000_cps<-cov_het_1000_cps[,-1]  ##100 200
cov_het_1000_lw<-read.xlsx("C:/Users/RimJhim/Documents/s1000_cov_lwo.xlsx")
cov_het_1000_lw<-cov_het_1000_lw[,-1]  ##100 140

##find column means of imported files--they are being made and so wait
cps1000het<-data.frame(apply(cov_het_1000_cps,1,function(x){sum(x)/200}))
lw1000het<-data.frame(apply(cov_het_1000_lw,1,function(x){sum(x)/200}))

cov1000het_av<-cbind.data.frame(cov_it_1000_het,cps1000het,lw1000het)
names(cov1000het_av)<-c("DCP","DCP-opt","CQR","CQR-m","CQR-r","CP-reg","CP-RF","CP-LW")
cov1000het_av<-cov1000het_av[,-1]

library(ggplot2)
cov1000het_av<- melt(cov1000het_av)

theme_set(theme_gray(base_size = 20))
p1<- ggplot(data = cov1000het_av, aes(x = variable, y = value))
p1 + geom_boxplot() + labs(x = "Method", y = "Average Coverage/iteration")+geom_hline(yintercept=0.9,
                                                                                     col="red",size=0.75)

############################################################################################################
############################################################################################################
plotfunc<-function(k,dat,dat2,dat3,n,it){
  x<-data.frame(dat[k])
  x<-empcov(n,x,it)
  #dat2<-read.xlsx("")
  x1<-data.frame(apply(dat2,1,function(x){sum(x)/n}))
  x2<-data.frame(apply(dat3,1,function(x){sum(x)/n}))
  dat4<-cbind.data.frame(x,x1,x2)
  names(dat4)<-c("DCP","DCP-opt","CQR","CQR-m","CQR-r","CP-reg","CP-RF","CP-LW")
  dat4<-dat4[,which(apply(dat4,2,var)!=0)] ##removing columns with constant values
  dat5<-melt(dat4)
  theme_set(theme_gray(base_size = 20))
  plotg<-ggplot(data = dat5, aes(x = variable, y = value))
  plotg2<-plotg+geom_boxplot()+labs(x = "Method", y = "Average Coverage/iteration")+geom_hline(yintercept=0.9,col="red",size=0.75)
  return(plotg2)
}
###implementation
plotfunc(8,sin1000_o,cov_het_1000_cps,cov_het_1000_lw,200,100)

##Loading the xlsx files and apply function


#########################################################################################
####Heteroskedastic case-
##Sammple size 700

cov_het_700_cps<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s700_cov_cps.xlsx")
cov_het_700_cps<-cov_het_700_cps[,-1]
cov_het_700_lw<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s700_cov_lwo.xlsx")
cov_het_700_lw<-cov_het_700_lw[,-1]

het700c<-plotfunc(7,sin700_o,cov_het_700_cps,cov_het_700_lw,140,50)

##Sample size of 1000
cov_het_1000_cps<-read.xlsx("C:/Users/RimJhim/Documents/s1000_cov_cps.xlsx")
cov_het_1000_cps<-cov_het_700_cps[,-1]
cov_het_1000_lw<-read.xlsx("C:/Users/RimJhim/Documents/s1000_cov_lwo.xlsx")
cov_het_1000_lw<-cov_het_700_lw[,-1]

het1000c<-plotfunc(7,sin1000_o,cov_het_1000_cps,cov_het_1000_lw,200,100)

######Sample size of 2000
cov_het_2000_cps<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s2000_cov_cps.xlsx")
cov_het_2000_cps<-cov_het_2000_cps[,-1]
cov_het_2000_lw<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s2000_cov_lwo.xlsx")
cov_het_2000_lw<-cov_het_2000_lw[,-1]

het2000c<-plotfunc(7,sin2000_o,cov_het_2000_cps,cov_het_2000_lw,400,100)

#########################################################################################
###Homoskedasticity

##Sammple size 700

cov_hom_700_cps<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s700_cov_cpsh.xlsx")
cov_hom_700_cps<-cov_hom_700_cps[,-1]
cov_hom_700_lw<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s700_cov_lwoh.xlsx")
cov_hom_700_lw<-cov_hom_700_lw[,-1]

hom700c<-plotfunc(8,sin1000_o,cov_hom_700_cps,cov_hom_700_lw,140,50)

###Sample size 1000

cov_hom_1000_cps<-read.xlsx("C:/Users/RimJhim/Documents/s1000_cov_cpsh.xlsx")
cov_hom_1000_cps<-cov_hom_1000_cps[,-1]
cov_hom_1000_lw<-read.xlsx("C:/Users/RimJhim/Documents/s1000_cov_lwoh.xlsx")
cov_hom_1000_lw<-cov_hom_1000_lw[,-1]

hom1000c<-plotfunc(8,sin1000_o,cov_hom_1000_cps,cov_hom_1000_lw,200,100)

###Sample size of 2000
cov_hom_2000_cps<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s2000_cov_cpsh.xlsx")
cov_hom_2000_cps<-cov_hom_2000_cps[,-1]
cov_hom_2000_lw<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s2000_cov_lwoh.xlsx")
cov_hom_2000_lw<-cov_hom_2000_lw[,-1]

hom2000c<-plotfunc(8,sin2000_o,cov_hom_2000_cps,cov_hom_2000_lw,400,100)

#########################################################################################
###Length-700
leng_het_700_cps<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s700_leng_cps.xlsx")
leng_het_700_cps<-leng_het_700_cps[,-1]  #
leng_het_700_lw<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s700_leng_lwo.xlsx")
leng_het_700_lw<-leng_het_700_lw[,-1]  ##

het700l<-plotfunc2(9,sin700_o,leng_het_700_cps,leng_het_700_lw,140,50)

###1000-het
leng_het_1000_cps<-read.xlsx("C:/Users/RimJhim/Documents/s1000_leng_cps.xlsx")
leng_het_1000_cps<-leng_het_1000_cps[,-1]  #
leng_het_1000_lw<-read.xlsx("C:/Users/RimJhim/Documents/s1000_leng_lwo.xlsx")
leng_het_1000_lw<-leng_het_1000_lw[,-1]  ##

het1000l<-plotfunc2(9,sin1000_o,leng_het_1000_cps,leng_het_1000_lw,200,100)

###2000-het
leng_het_2000_cps<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s2000_leng_cps.xlsx")
leng_het_2000_cps<-leng_het_2000_cps[,-1]  #
leng_het_2000_lw<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s2000_leng_lwo.xlsx")
leng_het_2000_lw<-leng_het_2000_lw[,-1]  ##

het2000l<-plotfunc2(9,sin2000_o,leng_het_2000_cps,leng_het_2000_lw,400,100)

#########Length-homoskedastic#######################################
leng_hom_700_cps<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s700_leng_cpsh.xlsx")
leng_hom_700_cps<-leng_hom_700_cps[,-1]  #
leng_hom_700_lw<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s700_leng_lwoh.xlsx")
leng_hom_700_lw<-leng_hom_700_lw[,-1]  ##

hom700l<-plotfunc2(10,sin700_o,leng_hom_700_cps,leng_hom_700_lw,140,50)

###1000-het
leng_hom_1000_cps<-read.xlsx("C:/Users/RimJhim/Documents/s1000_leng_cpsh.xlsx")
leng_hom_1000_cps<-leng_hom_1000_cps[,-1]  #
leng_hom_1000_lw<-read.xlsx("C:/Users/RimJhim/Documents/s1000_leng_lwoh.xlsx")
leng_hom_1000_lw<-leng_hom_1000_lw[,-1]  ##

hom1000l<-plotfunc2(10,sin1000_o,leng_hom_1000_cps,leng_hom_1000_lw,200,100)

###2000-het
leng_hom_2000_cps<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s2000_leng_cpsh.xlsx")
leng_hom_2000_cps<-leng_hom_2000_cps[,-1]  #
leng_hom_2000_lw<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s2000_leng_lwoh.xlsx")
leng_hom_2000_lw<-leng_hom_2000_lw[,-1]  ##

hom2000l<-plotfunc2(10,sin2000_o,leng_hom_2000_cps,leng_hom_2000_lw,400,100)

figure <- ggarrange(het700c,het1000c,het2000c,
                    hom700c,hom1000c,hom2000c,
                    het700l,het1000l,het2000l,
                    hom700l,hom1000l,hom2000l,
                    labels = c("700","1000","2000","700","1000","2000",
                               "700","1000","2000","700","1000","2000"),
                    ncol = 3, nrow = 4)
figure




####2000 iteration

cov2000het<-data.frame(sin2000_o[8])
cov_it_2000_het<-empcov(400,cov2000het,100)  ##100 6 

cov_het_2000_cps<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s2000_cov_cps.xlsx")
cov_het_2000_cps<-cov_het_2000_cps[,-1]  ##100 400
cov_het_2000_lw<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s2000_cov_lwo.xlsx")
cov_het_2000_lw<-cov_het_2000_lw[,-1]  ##100 400

cps2000het<-data.frame(apply(cov_het_2000_cps,1,function(x){sum(x)/400}))
lw2000het<-data.frame(apply(cov_het_2000_lw,1,function(x){sum(x)/400}))

cov2000het<-cbind.data.frame(cov_it_2000_het,cps2000het,lw2000het)
names(cov2000het)<-c("DCP","DCP-opt","CQR","CQR-m","CQR-r","CP-reg","CP-RF","CP-LW")

cov2000het<- melt(cov2000het)
theme_set(theme_gray(base_size = 20))
p2<- ggplot(data = cov2000het, aes(x = variable, y = value))
p2 + geom_boxplot() + labs(x = "Method", y = "Average Coverage/iteration")+geom_hline(yintercept=0.9,
                                                                                     col="red",size=0.75)

###################################################################
##Do it for the length for het cases
leng_het_1000_cps<-read.xlsx("C:/Users/RimJhim/Documents/s1000_leng_cps.xlsx")
leng_het_1000_cps<-leng_het_1000_cps[,-1]  ##100 200
leng_het_1000_lw<-read.xlsx("C:/Users/RimJhim/Documents/s1000_leng_lwo.xlsx")
leng_het_1000_lw<-leng_het_1000_lw[,-1]  ##100 140


###plot function for length
plotfunc2<-function(k,dat,dat2,dat3,n,it){
  x<-data.frame(dat[k])
  x<-empcov(n,x,it)
  #dat2<-read.xlsx("")
  x1<-data.frame(apply(dat2,1,function(x){sum(x)/n}))
  x2<-data.frame(apply(dat3,1,function(x){sum(x)/n}))
  dat4<-cbind.data.frame(x,x1,x2)
  names(dat4)<-c("DCP","DCP-opt","CQR","CQR-m","CQR-r","CP-reg","CP-RF","CP-LW")
  dat4<-dat4[,which(apply(dat4,2,var)!=0)] ##removing columns with constant values
  dat5<-melt(dat4)
  theme_set(theme_gray(base_size = 20))
  plotg<-ggplot(data = dat5, aes(x = variable, y = value))
  plotg2<-plotg+geom_boxplot()+labs(x = "Method", y = "Average Length/iteration")
  return(plotg2)
}

###sample size of 1000
plotfunc2(10,sin1000_o,leng_het_1000_cps,leng_het_1000_lw,200,100)

###for sample size of 700 and 2000

leng_het_700_cps<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s700_leng_cps.xlsx")
leng_het_700_cps<-leng_het_700_cps[,-1]  ##100 400
leng_het_700_lw<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s700_leng_lwo.xlsx")
leng_het_700_lw<-leng_het_700_lw[,-1]  ##100 400


###sample size of 700
plotfunc2(10,sin700_o,leng_het_700_cps,leng_het_700_lw,140,50)

###sample size of 2000
leng_het_2000_cps<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s2000_leng_cps.xlsx")
leng_het_2000_cps<-leng_het_2000_cps[,-1]  ##
leng_het_2000_lw<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s2000_leng_lwo.xlsx")
leng_het_2000_lw<-leng_het_2000_lw[,-1]  ##

###sample size of 2000
plotfunc2(10,sin2000_o,leng_het_2000_cps,leng_het_2000_lw,400,100)



###Repeat for homoskedastic case

##700 
cov_hom_700_cps<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s700_cov_cpsh.xlsx")
cov_hom_700_cps<-cov_hom_700_cps[,-1]  ##
cov_hom_700_lw<-read.xlsx("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/sinecurve/results/s700_cov_lwoh.xlsx")
cov_hom_700_lw<-cov_hom_700_lw[,-1]

plotfunc(8,sin700_o,cov_hom_700_cps,cov_hom_700_lw,200,100)


###############################################
library(openxlsx)
write.xlsx(sin700,"sin700.xlsx")
write.xlsx(sin1000,""
cov1,cov2,leng1,leng2,res.cov.cond,res.cov.cond2,
cov.mat,cov.mat2,leng.mat,leng.mat2,
sx0,sy0,sx1,sy1,stx,sty,
sx02,sy02,sx12,sy12,stx2,sty2,j



cov1,cov2,leng1,leng2,res.cov.cond,res.cov.cond2,
cov.mat,cov.mat2,leng.mat,leng.mat2,
sx0,sy0,sx1,sy1,stx,sty,
sx02,sy02,sx12,sy12,stx2,sty2,j