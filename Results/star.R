##Star
library(tidyverse)
library(dplyr)


install.packages("mlmRev")
library(mlmRev)
library(quantreg)
library(hdm)
library(matrixcalc)
library(ggplot2)
library(readxl)

library(cowplot)
library(patchwork)



data(star)
star<-na.omit(star)
dim(star)  ##22815 18
names(star)
unique(star$gr)  ##K 1 2 3
unique(star$cltype)  ##small reg   reg+A
unique(star$hdeg)  ##BS/BA     MS/MA/MEd Ed.S      MA+       Ed.D/Ph.D
unique(star$schtype)  ##inner  suburb rural  urban 
unique(star$sx)  ##M F
unique(star$eth) ##B W O A I H
##yrs is numeric version of grade-so drop grade
unique(star$ses) ##F N
unique(star$trace) ##B W A
unique(star$clad)  ##1    APPR PROB NOT  3    2    PEND
unique(star$ses)


##data processing

##cltype
star$type1<-ifelse(star$schtype=="small",1,0)
star$type2<-ifelse(star$schtype=="reg",1,0)

##hdeg
star$ted1<-ifelse(star$hdeg=="BS/BA",1,0)
star$ted2<-ifelse(star$hdeg=="MS/MA/MEd",1,0)
star$ted3<-ifelse(star$hdeg=="Ed.S",1,0)
star$ted4<-ifelse(star$hdeg=="MA+",1,0)

##schtype
star$school1<-ifelse(star$schtype=="inner",1,0)
star$school2<-ifelse(star$schtype=="suburb",1,0)
star$school3<-ifelse(star$schtype=="rural",1,0)

##sx
star$gender<-ifelse(star$sx=="F",1,0)

##eth
star%>%count(eth)
star$srace<-ifelse(star$eth=="W",1,0)

##trace
star%>%count(trace)
star$tr<-ifelse(star$trace=="W",1,0)

##clad
star$c1<-ifelse(star$clad=="1",1,0)
star$c2<-ifelse(star$clad=="APPR",1,0)
star$c3<-ifelse(star$clad=="PROB",1,0)
star$c4<-ifelse(star$clad=="NOT",1,0)
star$c5<-ifelse(star$clad=="3",1,0)
star$c6<-ifelse(star$clad=="2",1,0)

##ses
star$ss<-ifelse(star$ses=="F",1,0)

s2<-cbind.data.frame(star$type1,star$type2,star$ted1,star$ted2,star$ted3,star$ted4,
                     star$school1,star$school2,star$school3,star$gender,star$srace,star$tr,
                     star$c1,star$c2,star$c3,star$c4,star$c5,star$c6,star$ss,star$yrs,star$exp,star$read,star$math)
names(s2)<-c("type1","type2","ted1","ted2","ted3","ted4","school1","school2","school3",
             "gender","srace","tr","c1","c2","c3","c4","c5","c6","ss","yrs","exp",
             "read","math")
s2$grade<-s2$read+s2$math
drop <- c("math","read")
s2<-s2[ , !(names(s2) %in% drop)]
s2

##PREDICTIONS
x<- model.matrix(~1 +type1+type2+ted1+ted2+ted3+ted4+school1+school2+school3+gender+srace+tr+c1+     
                   c2+c3+c4+c5+c6+ss+yrs+exp,data=s2)
x<- x[, which(apply(x, 2, var) != 0)]

y<-s2$grade

alpha.sig<-0.1
T.ho <- floor(length(y)*0.1)

cov.mat <- leng.mat <- X.test.mat <- NULL
X0.mat<-Y0.mat<-X1.mat<-Y1.mat<-Y.test.mat<-NULL

start<- Sys.time()


for (i in 1:50){
  
  #print(r)
  
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
  taus    <- seq(0.001,0.999,length=200)  ##change taus to avoid infinite values 
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
  Y0.mat <- rbind(Y0.mat,Y0)
  Y1.mat <- rbind(Y1.mat,Y1)
  X0.mat <- rbind(X0.mat,X0)
  X1.mat <- rbind(X1.mat,X1)
  
  cov.mat.temp  <- cbind(res.qr$cov.qr,res.opt$cov.opt,res.cqr$cov.o,res.cqr$cov.m,res.cqr$cov.r,res.reg$cov.reg)
  leng.mat.temp <- cbind(res.qr$leng.qr,res.opt$leng.opt,res.cqr$leng.o,res.cqr$leng.m,res.cqr$leng.r,res.reg$leng.reg)
  
  cov.mat   <- rbind(cov.mat,cov.mat.temp)
  leng.mat  <- rbind(leng.mat,leng.mat.temp)
  
}
Sys.time()-start
##Time difference of 1.793022 hours
k<-colMeans(leng.mat,na.rm=TRUE)
h<-colMeans(cov.mat,na.rm=TRUE)


# Estimate predicted coverage
pred.cov.qr       <- predict(glm(cov.mat[,1]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.opt      <- predict(glm(cov.mat[,2]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.cqr.o    <- predict(glm(cov.mat[,3]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.cqr.m    <- predict(glm(cov.mat[,4]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.cqr.r    <- predict(glm(cov.mat[,5]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.reg      <- predict(glm(cov.mat[,6]~X.test.mat,family=binomial(link="logit")),type="response")
#pred.cov.loc      <- predict(glm(cov.mat[,8]~X.test.mat,family=binomial(link="logit")),type="response")

# Dispersion of conditional coverage

pred.cov <- cbind(pred.cov.qr,pred.cov.opt,pred.cov.cqr.o,pred.cov.cqr.m,pred.cov.cqr.r,pred.cov.reg)
# for comparison sqrt(colMeans((pred.cov-0.9)^2))*100
res.cov.cond <- apply(pred.cov,2,sd)

d<-cbind.data.frame(k,h,res.cov.cond)
write_xlsx(d,"star_5iter_res.xlsx")
write_xlsx(data.frame(cov.mat),"star_5iter_cov.xlsx")
write_xlsx(data.frame(leng.mat),"star_5iter_leng.xlsx")

##Also get the excel files of the train,test dfs

write_xlsx(data.frame(X0.mat),"X0_star.xlsx")
write_xlsx(data.frame(Y0.mat),"Y0_star.xlsx")
write_xlsx(data.frame(X1.mat),"X1_star.xlsx")
write_xlsx(data.frame(Y1.mat),"Y1_star.xlsx")
write_xlsx(data.frame(X.test.mat),"Xtest.xlsx")
write_xlsx(data.frame(Y.test.mat),"Ytest.xlsx")

#####



##50 iterations

C<-read.xlsx("star_5iter_cov.xlsx")
L<-read.xlsx("star_5iter_leng.xlsx")

k<-colMeans(C,na.rm=TRUE)
h<-colMeans(L,na.rm=TRUE)

Xt<-read.xlsx("Xtest.xlsx")

pred.cov.qr       <- predict(glm(C[,1]~as.matrix(Xt),family=binomial(link="logit")),type="response")
pred.cov.opt      <- predict(glm(C[,2]~as.matrix(Xt),family=binomial(link="logit")),type="response")
pred.cov.cqr.o    <- predict(glm(C[,3]~as.matrix(Xt),family=binomial(link="logit")),type="response")
pred.cov.cqr.m    <- predict(glm(C[,4]~as.matrix(Xt),family=binomial(link="logit")),type="response")
pred.cov.cqr.r    <- predict(glm(C[,5]~as.matrix(Xt),family=binomial(link="logit")),type="response")
pred.cov.reg      <- predict(glm(C[,6]~as.matrix(Xt),family=binomial(link="logit")),type="response")

pred.cov <- cbind(pred.cov.qr,pred.cov.opt,pred.cov.cqr.o,pred.cov.cqr.m,pred.cov.cqr.r,pred.cov.reg)

res.cov.cond <- apply(pred.cov,2,sd)

d<-cbind.data.frame(k,h,res.cov.cond)
rownames(d)<-c("DCP","DCPopt","CQR","CQRm","CQRr","CPreg")
names(d)<-c("Avg Length","Avg.Coverage","SD")
formatn<-function(x){
  #if (x>0.09){y=x}
  #else {y=formatC(x,format="e",digits=2)}
  y<-ifelse(abs(x)>=0.01,round(x,4),formatC(x,format="e",digits=2))
  return(y)
}
d[]<-lapply(d,formatn)
xtable(d,digits=4)

###################################################################
####################################################################



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

########################################################################################



###EDA
par(mfrow=c(4,2))
boxplot(grade~gender,data=s2)
boxplot(grade~srace,data=s2,xlab="student's race")
boxplot(grade~tr,data=s2,xlab="teacher's race")
star$grade<-star$math+star$read
boxplot(grade~schtype,data=star,xlab="school type")
boxplot(grade~cltype,data=star,xlab="class type")
boxplot(grade~hdeg,data=star,xlab="education level of teacher")
boxplot(grade~ses,data=star,xlab="socioeconomic status of student")
boxplot(grade~gr,data=star,xlab="year of schooling")

mod<-lm(grade~type1+type2+ted1+ted2+ted3+ted4+school1+school2+school3+gender+srace+tr+c1+     
          c2+c3+c4+c5+c6+ss+yrs+exp,data=s2)
hist(mod$residuals,xlab="Residuals")
library(lmtest)
lmtest::bptest(mod)

#library(car)
#ncvTest(mod)


####################################################################################################
#######################################################################################################
###model with interactions
#####################################################################################################
x<- model.matrix(~-1 +(type1+type2+ted1+ted2+ted3+ted4+school1+school2+school3+gender+srace+tr+c1+     
                         c2+c3+c4+c5+c6+ss+yrs+exp)^2,data=s2)

x<- x[, which(apply(x, 2, var) != 0)]

y<-s2$grade


alpha.sig<-0.1
T.ho <- floor(length(y)*0.1)

cov.mat <- leng.mat <- X.test.mat <- NULL
X0.mat<-Y0.mat<-X1.mat<-Y1.mat<-Y.test.mat<-NULL

start<- Sys.time()


for (i in 1:50){
  
  #print(r)
  
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
  taus    <- seq(0.001,0.999,length=200)  ##change taus to avoid infinite values 
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
  Y0.mat <- rbind(Y0.mat,Y0)
  Y1.mat <- rbind(Y1.mat,Y1)
  X0.mat <- rbind(X0.mat,X0)
  X1.mat <- rbind(X1.mat,X1)
  
  cov.mat.temp  <- cbind(res.qr$cov.qr,res.opt$cov.opt,res.cqr$cov.o,res.cqr$cov.m,res.cqr$cov.r,res.reg$cov.reg)
  leng.mat.temp <- cbind(res.qr$leng.qr,res.opt$leng.opt,res.cqr$leng.o,res.cqr$leng.m,res.cqr$leng.r,res.reg$leng.reg)
  
  cov.mat   <- rbind(cov.mat,cov.mat.temp)
  leng.mat  <- rbind(leng.mat,leng.mat.temp)
  
}
Sys.time()-start
k<-colMeans(leng.mat,na.rm=TRUE)
h<-colMeans(cov.mat,na.rm=TRUE)


# Estimate predicted coverage
pred.cov.qr       <- predict(glm(cov.mat[,1]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.opt      <- predict(glm(cov.mat[,2]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.cqr.o    <- predict(glm(cov.mat[,3]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.cqr.m    <- predict(glm(cov.mat[,4]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.cqr.r    <- predict(glm(cov.mat[,5]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.reg      <- predict(glm(cov.mat[,6]~X.test.mat,family=binomial(link="logit")),type="response")
#pred.cov.loc      <- predict(glm(cov.mat[,8]~X.test.mat,family=binomial(link="logit")),type="response")

# Dispersion of conditional coverage

pred.cov <- cbind(pred.cov.qr,pred.cov.opt,pred.cov.cqr.o,pred.cov.cqr.m,pred.cov.cqr.r,pred.cov.reg)
# for comparison sqrt(colMeans((pred.cov-0.9)^2))*100
res.cov.cond <- apply(pred.cov,2,sd)
