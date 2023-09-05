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

#######################################################################
##(1-alpha) PREDICTION INTERVAL FOR CRITICAL TEMPERATURE: alpha=0.05
######################################################################
library(openxlsx)
library(ggplot2)
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
library(cowplot)

###All functions under functions_final.R as well as 
##functions used after data generation are 
##from:
#Chernozhukov, V.,Wuthrich,K. &  Zhu,Y.(2021).Distributional Conformal Prediction. PNAS,118 (48).
##https://doi.org/10.1073/pnas.2107794118

##The relevant code taken from the authors is in their github account:
#https://github.com/kwuthrich/Replication_DCP 


##################################################################
#####LOADING DATA#################################################
####################################################################
train<-read.csv("train.csv",header=TRUE,sep=",")  
###Data source: https://archive.ics.uci.edu/dataset/464/superconductivty+data

dim(train)   ####21263    82

y<-train$critical_temp  ###response variable 

#####model matrix of covariates
formd <- as.formula(paste("~-1+",paste(names(train)[1:81],collapse="+")))
mm<-model.matrix(formd,data=train)

###Changing alpha.sig to 0.05
alpha.sig<-0.05
T.ho <- floor(length(y)*0.2)  ###80% training and calibration-20% testing 

cov.mat <- leng.mat <- X.test.mat <- NULL
X0.mat<-Y0.mat<-Y.test.mat<-X1.mat<-Y1.mat<-NULL

source("functions_final.R")

start<-Sys.time()
for (i in 1:50){
  
 
  set.seed(i)
  ####checking if matrix is singular or not-required for quantile regression
  sing <- TRUE
  while (sing==TRUE){
    
    ind <- sample(length(y),length(y),replace=FALSE)
    Y   <- y[ind]
    X   <- mm[ind,]
    
    ind.test  <- (length(Y)-T.ho+1):length(Y)
    Y.test    <- Y[ind.test]
    X.test    <- cbind(X[ind.test,])
    
    obj.uneven <- aux.uneven(Y[-ind.test],cbind(X[-ind.test,]))
    
    X0 <- cbind(obj.uneven$X0)
    Y0 <- obj.uneven$Y0
    X1 <- cbind(obj.uneven$X1)
    Y1 <- obj.uneven$Y1
    Y0<-as.matrix(Y0)
    Y1<-as.matrix(Y1)
    
    sing  <- is.singular.matrix(t(X0)%*%X0)
  }
  
  ###grid of probability values for DCP-QR and DCP-opt
  taus    <- seq(0.05,0.95,length=200)  
  
  # Conformal prediction methods 
  res.qr    <- dcp.qr(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)
  res.opt   <- dcp.opt(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)
  res.cqr   <- cqr(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
  res.reg   <- cp.reg(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
  
  
  # Results
  X.test.mat <- rbind(X.test.mat,X.test) ##test set-covariates
  X0.mat <- rbind(X0.mat,X0)  ##train set-covariates
  Y0.mat <- rbind(Y0.mat,Y0)  ##train set-response
  X1.mat <- rbind(X1.mat,X1)  ##test set-covariates
  Y1.mat <- rbind(Y1.mat,Y1)  ##test set-response 
  Y.test.mat <- rbind(Y.test.mat,Y.test) ##test set-response
  
  ###storing coverage and length values from each iteration and method
  cov.mat.temp  <- cbind(res.qr$cov.qr,res.opt$cov.opt,res.cqr$cov.o,res.cqr$cov.m,res.cqr$cov.r,res.reg$cov.reg)
  leng.mat.temp <- cbind(res.qr$leng.qr,res.opt$leng.opt,res.cqr$leng.o,res.cqr$leng.m,res.cqr$leng.r,res.reg$leng.reg)
  
  ####binding the coverage and length values together 
  cov.mat   <- rbind(cov.mat,cov.mat.temp)
  leng.mat  <- rbind(leng.mat,leng.mat.temp)
  
}
Sys.time()-start  ###5.501 hours


###Computing column means to find average coverage and average length over data generated over all 50 iterations
colMeans(leng.mat,na.rm=TRUE)  ##212600 6
colMeans(cov.mat,na.rm=TRUE)  ##212600 6 


# Estimate predicted coverage fro logistic regression
pred.cov.qr       <- predict(glm(cov.mat[,1]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.opt      <- predict(glm(cov.mat[,2]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.cqr.o    <- predict(glm(cov.mat[,3]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.cqr.m    <- predict(glm(cov.mat[,4]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.cqr.r    <- predict(glm(cov.mat[,5]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.reg      <- predict(glm(cov.mat[,6]~X.test.mat,family=binomial(link="logit")),type="response")


# Dispersion of conditional coverage

pred.cov <- cbind(pred.cov.qr,pred.cov.opt,pred.cov.cqr.o,pred.cov.cqr.m,pred.cov.cqr.r,pred.cov.reg)
res.cov.cond <- apply(pred.cov,2,sd)

####Preparing coverage, length and standard deviation results from latex
##########################################################################################
c<-colMeans(cov.mat,na.rm=TRUE)
l<-colMeans(leng.mat,na.rm=TRUE)

s<-data.frame(c,l,res.cov.cond)
names(s)<-c("Coverage","Length","SD")
rownames(s)<-c("DCP-QR","DCP-opt","CQR","CQR-m","CQR-r","CP-reg")
xtable(s,digits=4)

################################################################################
################################################################################
#################################################################################
#################################################################################
###Visualizing the spread of predicted coverage 

################################################################################
alg<-c("DCP-QR","DCP-opt","CQR","CQR-m","CQR-r","CP-reg")
pred.cov<-data.frame(pred.cov)


###Source: https://www.statology.org/mode-in-r/
find_mode<-function(x){
    u<-unique(x)
    tab<-tabulate(match(x,u))
    u[tab==max(tab)]
}

####function to render histograms by ggplot
prep<-function(z,alg,i){
  a<-ggplot(data=pred.cov,aes(z))+geom_histogram()+
  geom_vline(xintercept=0.95,col="red")+xlab("Predicted coverage probability")+
  ylab("Frequency")
  a1<-a+annotate(x=find_mode(z),y=150000,label=alg[i],geom="label")
  return(a1)
}

prep(pred.cov.qr,alg,1)
prep(pred.cov.opt,alg,2)
prep(pred.cov.cqr.o,alg,3)
prep(pred.cov.cqr.m,alg,4)
prep(pred.cov.cqr.r,alg,5)
prep(pred.cov.reg,alg,6)


#################################################################################


write.csv(data.frame(X0.mat),"X0.csv")

write.csv(data.frame(X1.mat),"X1.csv")

write.csv(data.frame(X.test.mat),"X.test.mat.csv")

write.csv(data.frame(Y0.mat),"Y0.csv")

write.csv(data.frame(Y1.mat),"Y1c.csv")

write.csv(data.frame(Y.test.mat),"Ytest.csv")

cov<-data.frame(cov.mat)


write.csv(cov,"COV_CT.csv")
write.csv(data.frame(leng.mat),"LENG_CT.csv")

##################################################################################
###################################################################################
###################################################################################
####Boxplots of per iteration average coverage 
#####################################################################################
#####################################################################################

####Compute avg.coverage and avg.length 
cov<-data.frame(cov.mat)
leng<-data.frame(leng.mat)

n0<-nrow(cov)/50

c1<-data.frame(aggregate(cov, list(rep(1:(nrow(cov) %/% n0 + 1), each = n0, len = nrow(cov))), mean)[-1])
l1<-data.frame(aggregate(leng, list(rep(1:(nrow(leng) %/% n0 + 1), each = n0, len = nrow(leng))), mean)[-1])

names(c1)<-names(l1)<-c("DCPQR","DCPopt","CQR","CQRm","CQRr","CPreg")

c1$x<-l1$x<-"Method"

   
c2<-melt(c1,id.vars = "x")


b<-ggplot(c2,aes(variable,value,col=variable))+geom_boxplot()+stat_summary(fun = mean)+
  stat_summary(fun = mean, geom = "text", col = "black",  
               vjust = 1.5, aes(label = paste(round(..y.., digits = 3))))+
  geom_hline(yintercept=0.95,col="darkred",linewidth=1)+
  xlab("Method")+ylab("Average Coverage")

b

l2<-melt(l1,id.vars = "x")


b1<-ggplot(l2,aes(variable,value,col=variable))+geom_boxplot()+stat_summary(fun = mean)+
  stat_summary(fun = mean, geom = "text", col = "black",  
               vjust = 1.5, aes(label = paste(round(..y.., digits = 3))))+
  xlab("Method")+ylab("Average Length")

b1


grid.arrange(b,b1)  


