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


##conduct 
library(openxlsx)
train<-read.csv("train.csv",header=TRUE,sep=",")
m<-read.csv("unique_m.csv",header=TRUE,sep=",")
dim(train)
dim(m)
y<-train$critical_temp
##

######Visualization
#########################################################################
hist(train$critical_temp,xlab="Critical Temperature")
boxplot(train$critical_temp,ylab="Critical Temperature")

#########################################################################
l<-lm(critical_temp~.,data=train)

install.packages("handyplots")
library(handyplots)

resplot(l,highlight.outliers=TRUE)

formd <- as.formula(paste("~-1+",paste(names(train)[1:81],collapse="+")))
mm<-model.matrix(formd,data=train)

###
alpha.sig<-0.1
T.ho <- floor(length(y)*0.2)

cov.mat <- leng.mat <- X.test.mat <- NULL
X0.mat<-Y0.mat<-Y.test.mat<-X1.mat<-Y1.mat<-NULL

for (i in 1:30){
  
  #print(r)
  
  # Define training and holdout samples (while loop is to avoid singular design)
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
  X0.mat <- rbind(X0.mat,X0)
  Y0.mat <- rbind(Y0.mat,Y0)
  X1.mat <- rbind(X1.mat,X1)
  Y1.mat <- rbind(Y1.mat,Y1)
  Y.test.mat <- rbind(Y.test.mat,Y.test)
  
  cov.mat.temp  <- cbind(res.qr$cov.qr,res.opt$cov.opt,res.cqr$cov.o,res.cqr$cov.m,res.cqr$cov.r,res.reg$cov.reg)
  leng.mat.temp <- cbind(res.qr$leng.qr,res.opt$leng.opt,res.cqr$leng.o,res.cqr$leng.m,res.cqr$leng.r,res.reg$leng.reg)
  
  cov.mat   <- rbind(cov.mat,cov.mat.temp)
  leng.mat  <- rbind(leng.mat,leng.mat.temp)
  
}
colMeans(leng.mat,na.rm=TRUE)
colMeans(cov.mat,na.rm=TRUE)


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

#pred.cov.qr   pred.cov.opt pred.cov.cqr.o pred.cov.cqr.m pred.cov.cqr.r   pred.cov.reg 
#0.01155703     0.04235030     0.01175060     0.01257260     0.01277881     0.11665529

###
setwd("D:/")
#setwd("C:/Users/RimJhim/Downloads/distribution_conformal_prediction/cqr_try/crittemp")
library(readxl)
library(writexl)

write_xlsx(data.frame(cov.mat),"COV_CT90.xlsx")
write_xlsx(data.frame(leng.mat),"LENG_CT90.xlsx")

####Saving results to Excel 
##################################################################
cov<-data.frame(read_excel("COV_CT.xlsx"))
leng<-data.frame(read_excel("LENG_CT.xlsx"))


####################################################################
###Preparing results for Latex
####################################################################
C<-colMeans(cov.mat)
L<-colMeans(leng.mat)
S<-res.cov.cond

#################################################################################
#################################################################################

lt<-cbind.data.frame(C,L,S)
rownames(lt)<-c("DCP","DCP-opt","CQR","CQR-m","CQR-r","CP-reg")
names(lt)<-c("Avg.Coverage","Avg.Length","Standard Deviation")

xtable(lt,digits=4)

##50 repeats 
n0<-nrow(cov)/30
g1<-data.frame(aggregate(cov, list(rep(1:(nrow(cov) %/% n0 + 1), each = n0, len = nrow(cov))), mean)[-1])


###Do same for length
l1<-data.frame(aggregate(leng, list(rep(1:(nrow(leng) %/% n0 + 1), each = n0, len = nrow(leng))), mean)[-1])

###boxplot
names(g1)<-names(l1)<-c("DCP-QR","DCP-opt","CQR","CQR-m","CQR-r","CP-reg")

c<-colMeans(cov,na.rm=TRUE)
l<-colMeans(leng,na.rm=TRUE)

Xt<-read_excel("Xtc.xlsx")
pred.cov.qr       <- predict(glm(cov[,1]~Xt,family=binomial(link="logit")),type="response")
pred.cov.opt      <- predict(glm(cov[,2]~Xt,family=binomial(link="logit")),type="response")
pred.cov.cqr.o    <- predict(glm(cov[,3]~Xt,family=binomial(link="logit")),type="response")
pred.cov.cqr.m    <- predict(glm(cov[,4]~Xt,family=binomial(link="logit")),type="response")
pred.cov.cqr.r    <- predict(glm(cov[,5~Xt,family=binomial(link="logit")),type="response")
pred.cov.reg      <- predict(glm(cov[,6]~Xt,family=binomial(link="logit")),type="response")


########################################################################################################
########################################################################################################
####95%####################################################################################################

alpha.sig<-0.05
T.ho <- floor(length(y)*0.2)

cov.mat <- leng.mat <- X.test.mat <- NULL
X0.mat<-Y0.mat<-Y.test.mat<-X1.mat<-Y1.mat<-NULL

for (i in 1:50){
  
  #print(r)
  
  # Define training and holdout samples (while loop is to avoid singular design)
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
  X0.mat <- rbind(X0.mat,X0)
  Y0.mat <- rbind(Y0.mat,Y0)
  X1.mat <- rbind(X1.mat,X1)
  Y1.mat <- rbind(Y1.mat,Y1)
  Y.test.mat <- rbind(Y.test.mat,Y.test)
  
  cov.mat.temp  <- cbind(res.qr$cov.qr,res.opt$cov.opt,res.cqr$cov.o,res.cqr$cov.m,res.cqr$cov.r,res.reg$cov.reg)
  leng.mat.temp <- cbind(res.qr$leng.qr,res.opt$leng.opt,res.cqr$leng.o,res.cqr$leng.m,res.cqr$leng.r,res.reg$leng.reg)
  
  cov.mat   <- rbind(cov.mat,cov.mat.temp)
  leng.mat  <- rbind(leng.mat,leng.mat.temp)
  
}

colMeans(leng.mat,na.rm=TRUE)
colMeans(cov.mat,na.rm=TRUE)

##colMeans(leng.mat,na.rm=TRUE)
#[1] 55.80296 54.02805 55.80966 56.12485 56.03888 70.78978
# colMeans(cov.mat,na.rm=TRUE)
# 0.9560066 0.9518156 0.9504421 0.9503951 0.9504233 0.9502117

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

#res.cov.cond
#pred.cov.qr   pred.cov.opt pred.cov.cqr.o pred.cov.cqr.m pred.cov.cqr.r   pred.cov.reg 
#0.006602592    0.019411681    0.006989461    0.009251765    0.009318539    0.069205963 

write_xlsx(data.frame(cov.mat),"COV_CT95.xlsx")
write_xlsx(data.frame(leng.mat),"LENG_CT95.xlsx")

C<-colMeans(cov.mat)
L<-colMeans(leng.mat)
S<-res.cov.cond

#################################################################################
#################################################################################

lt<-cbind.data.frame(C,L,S)
rownames(lt)<-c("DCP","DCP-opt","CQR","CQR-m","CQR-r","CP-reg")
names(lt)<-c("Avg.Coverage","Avg.Length","Standard Deviation")

xtable(lt,digits=4)

\begin{table}[ht]
\centering
\begin{tabular}{rrrr}
\hline
& Avg.Coverage & Avg.Length & Standard Deviation \\ 
\hline
DCP & 0.9560 & 55.8030 & 0.0066 \\ 
DCP-opt & 0.9518 & 54.0280 & 0.0194 \\ 
CQR & 0.9504 & 55.8097 & 0.0070 \\ 
CQR-m & 0.9504 & 56.1249 & 0.0093 \\ 
CQR-r & 0.9504 & 56.0389 & 0.0093 \\ 
CP-reg & 0.9502 & 70.7898 & 0.0692 \\ 
\hline
\end{tabular}
\end{table}


#########Histograms of standard deviation of conditional coverage probability
#####################################################################################
H<-cbind.data.frame(pred.cov.qr,pred.cov.opt,pred.cov.cqr.o,pred.cov.cqr.m,pred.cov.cqr.r,pred.cov.reg)
H$breaks<-seq(0,1,length.out=dim(H)[1])
names(H)<-c("DCP","DCPopt","CQR","CQRm","CQRr","CPreg","breaks")


f1<-ggplot2::ggplot(H,aes(x=DCP))+geom_histogram(color="darkblue", fill="lightblue")+
      geom_vline(xintercept = 0.95,col="red",lwd=1)+
      annotate(x=0.95,y=110000,label="0.95",geom="label")+
      ylab("Count")

f2<-ggplot2::ggplot(H,aes(x=DCPopt))+geom_histogram(color="darkblue", fill="lightblue")+
  geom_vline(xintercept = 0.95,col="red",lwd=1)+
  annotate(x=0.95,y=110000,label="0.95",geom="label")+
  ylab("Count")+xlab("DCP-opt")

f3<-ggplot2::ggplot(H,aes(x=CQR))+geom_histogram(color="darkblue", fill="lightblue")+
  geom_vline(xintercept = 0.95,col="red",lwd=1)+
  annotate(x=0.95,y=110000,label="0.95",geom="label")+
  ylab("Count")+xlab("CQR")

f4<-ggplot2::ggplot(H,aes(x=CQRm))+geom_histogram(color="darkblue", fill="lightblue")+
  geom_vline(xintercept = 0.95,col="red",lwd=1)+
  annotate(x=0.95,y=110000,label="0.95",geom="label")+
  ylab("Count")+xlab("CQR-m")

f5<-ggplot2::ggplot(H,aes(x=CQRr))+geom_histogram(color="darkblue", fill="lightblue")+
  geom_vline(xintercept = 0.95,col="red",lwd=1)+
  annotate(x=0.95,y=110000,label="0.95",geom="label")+
  ylab("Count")+xlab("CQR-r")

f6<-ggplot2::ggplot(H,aes(x=CPreg))+geom_histogram(color="darkblue", fill="lightblue")+
  geom_vline(xintercept = 0.95,col="red",lwd=1)+
  annotate(x=0.95,y=110000,label="0.95",geom="label")+
  ylab("Count")+xlab("CP-reg")
  
library(cowplot)
devtools::install_github("thomasp85/patchwork")
library(patchwork) 
plot_grid(f1,f2,f3,f4,f5,f6) 
