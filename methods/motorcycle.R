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
#####################MOTORCYCLE DATA###################################
###loading required packages




install.packages("MASS")
library(MASS)
data(mcycle)  ##times in -axis, accel in y-axis

y<-mcycle$accel
x<-as.matrix(mcycle$times)  ##25% testing, 75% training
##small dataset 

mcycle<-cbind.data.frame(x,y)
names(mcycle)<-c("time","accelertion")

ggplot(mcycle,aes(x=time,y=accelertion))+geom_point()+
  xlab("Time")+ylab("Acceleration")

yogen<-function(t,x,y,p){
  T.ho <- floor(length(y)*p)
  cov.mat <- leng.mat <- X.test.mat <- NULL
  Y.test.mat<-Y0.mat<-Y1.mat<-NULL
  X0.mat<-X1.mat<-NULL
  start<-Sys.time()
  for (i in 1:t){
    
    # Define training and holdout samples (while loop is to avoid singular design)
    sing <- TRUE
    while (sing==TRUE){
      set.seed(i)
      
      ind <- sample(length(y),length(y),replace=FALSE)
      Y   <- y[ind]
      X   <- as.matrix(x[ind,])
      
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
    X0.mat <- rbind(X0.mat,X0)
    Y0.mat <- rbind(Y0.mat,Y0)
    X1.mat <- rbind(X1.mat,X1)
    Y1.mat <- rbind(Y1.mat,Y1)
    
    
    cov.mat.temp  <- cbind(res.qr$cov.qr,res.opt$cov.opt,res.cqr$cov.o,res.cqr$cov.m,res.cqr$cov.r,res.reg$cov.reg)
    leng.mat.temp <- cbind(res.qr$leng.qr,res.opt$leng.opt,res.cqr$leng.o,res.cqr$leng.m,res.cqr$leng.r,res.reg$leng.reg)
    
    cov.mat   <- rbind(cov.mat,cov.mat.temp)
    leng.mat  <- rbind(leng.mat,leng.mat.temp)
    
  }
  l<-colMeans(leng.mat,na.rm=TRUE)
  c<-colMeans(cov.mat,na.rm=TRUE)
  
  
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
  j=Sys.time()-start
  
  dat<-cbind.data.frame(as.matrix(l),as.matrix(c),as.matrix(res.cov.cond))
  
  l<-list(dat,cov.mat,leng.mat,
          X.test.mat,Y.test.mat,X0.mat,Y0.mat,
          X1.mat,Y1.mat,j)
  return(l)
}

alpha.sig<-0.1

library(writexl)
h<-yogen(1000,x,y,0.2)
#h1<-yogen(2000,x,y,0.2)

X<-cbind.data.frame(data.frame(h[6]),data.frame(h[8]))
write_xlsx(X,"xtxv_mcycle2.xlsx")
write_xlsx(data.frame(h[4]),"xtest_mcycle.xlsx")
write_xlsx(data.frame(h[7]),"ytrain_mcycle2.xlsx")
write_xlsx(data.frame(h[9]),"yvalid_mcycle2.xlsx")
write_xlsx(data.frame(h[5]),"ytest_mcycle2.xlsx")

###coverage
write_xlsx(data.frame(h[2]),"covmat_2.xlsx")
##Length
write_xlsx(data.frame(h[3]),"lengmat_2.xlsx")

##results
write_xlsx(data.frame(h[1]),"res_2.xlsx")

install.packages("ggplot2")
library(ggplot2)





k<-yogen(1000,x,y,0.25)
#1<-yogen(2000,x,y,0.25)

X2<-cbind.data.frame(data.frame(k[6]),data.frame(k[8]))
write_xlsx(X2,"xtxv_mcycle22.xlsx")
write_xlsx(data.frame(k[4]),"xtest_mcycle22.xlsx")
write_xlsx(data.frame(k[7]),"ytrain_mcycle22.xlsx")
write_xlsx(data.frame(k[9]),"yvalid_mcycle22.xlsx")
write_xlsx(data.frame(k[5]),"ytest_mcycle22.xlsx")

###coverage
write_xlsx(data.frame(k[2]),"covmat_25.xlsx")
##Length
write_xlsx(data.frame(k[3]),"lengmat_25.xlsx")

##results
write_xlsx(data.frame(k[1]),"res_25.xlsx")






v<-yogen(2000,x,y,0.3)

X2<-cbind.data.frame(data.frame(v[6]),data.frame(v[8]))
write_xlsx(X2,"xtxv_mcycle23.xlsx")
write_xlsx(data.frame(v[4]),"xtest_mcycle23.xlsx")
write_xlsx(data.frame(v[7]),"ytrain_mcycle23.xlsx")
write_xlsx(data.frame(v[9]),"yvalid_mcycle23.xlsx")
write_xlsx(data.frame(v[5]),"ytest_mcycle23.xlsx")

###coverage
write_xlsx(data.frame(v[2]),"covmat_30.xlsx")
##Length
write_xlsx(data.frame(v[3]),"lengmat_30.xlsx")

##results
write_xlsx(data.frame(v[1]),"res_30.xlsx")

w<-yogen(1000,x,y,0.35)
#w1<-yogen(2000,x,y,0.35)
X2<-cbind.data.frame(data.frame(w[6]),data.frame(w[8]))
write_xlsx(X2,"xtxv_mcycle24.xlsx")
write_xlsx(data.frame(w[4]),"xtest_mcycle24.xlsx")
write_xlsx(data.frame(w[7]),"ytrain_mcycle24.xlsx")
write_xlsx(data.frame(w[9]),"yvalid_mcycle24.xlsx")
write_xlsx(data.frame(w[5]),"ytest_mcycle24.xlsx")

###coverage
write_xlsx(data.frame(w[2]),"covmat_35.xlsx")
##Length
write_xlsx(data.frame(w[3]),"lengmat_35.xlsx")

##results
write_xlsx(data.frame(w[1]),"res_35.xlsx")
##change split to have more test data 25%-30%

##now also get computation time
h[10] ##Time difference of 1.652863 mins
k[10] ##Time difference of 1.696985 mins
v[10] ##Time difference of 1.75321 mins
w[10] ##Time difference of 1.768663 mins


######################################################################
###Preparing dataframe for latex
D<-rbind.data.frame(data.frame(h[1]),data.frame(k[1]),data.frame(v[1]),data.frame(w[1]))
names(D)<-c("Length","Coverage","SD")
D

######################################################################
###Obtaining quantile bins############################################
######################################################################

num.seg <- 20
Xth<-as.vector(unlist(h[4]))
Xtk<-as.vector(unlist(k[4]))
Xtv<-as.vector(unlist(v[4]))
Xtw<-as.vector(unlist(w[4]))

ch<-as.matrix(data.frame(h[2]))
ck<-as.matrix(data.frame(k[2]))
cv<-as.matrix(data.frame(v[2]))
cw<-as.matrix(data.frame(w[2]))

lh<-as.matrix(data.frame(h[3]))
lk<-as.matrix(data.frame(k[3]))
lv<-as.matrix(data.frame(v[3]))
lw<-as.matrix(data.frame(w[3]))

cov.cond.h  <- binning(Xth,ch,num.seg)$cond
leng.cond.h <- binning(Xth,lh,num.seg)$cond

cov.cond.k  <- binning(Xtk,ck,num.seg)$cond
leng.cond.k <- binning(Xtk,lk,num.seg)$cond

cov.cond.v  <- binning(Xtv,cv,num.seg)$cond
leng.cond.v <- binning(Xtv,lv,num.seg)$cond

cov.cond.w  <- binning(Xtw,cw,num.seg)$cond
leng.cond.w <- binning(Xtw,lw,num.seg)$cond

#################################################################

QB<-function(dat){
d<-data.frame(dat)
d$seg<-1:20
names(d)<-c("DCP","DCPopt","CQR","CQRm","CQRr","CPreg","segment")
c<-ggplot2::ggplot(d,aes(x=segment,y=DCP))+geom_line(aes(x=segment,y=DCP,
    colour="DCP"),lwd=1,linetype="dashed")+
    geom_line(d,mapping=aes(x=segment,y=DCPopt,colour="DCPopt"),lwd=1,linetype="dashed")+
    geom_line(d,mapping=aes(x=segment,y=CQR,colour="CQR"),lwd=1,linetype="dashed")+
    geom_line(d,mapping=aes(x=segment,y=CQRm,colour="CQRm"),lwd=1,linetype="dashed")+
    geom_line(d,mapping=aes(x=segment,y=CQRr,colour="CQRr"),lwd=1,linetype="dashed")+
  geom_line(d,mapping=aes(x=segment,y=CPreg,colour="CPreg"),lwd=1,linetype="dashed")+
  xlab("Bin")+ylab("Conditional Coverage")+
  scale_colour_manual(name="Legend",
                      values=c(DCP="red", DCPopt="blue", CQR="darkgreen",
                               CQRm="orange",CQRr="magenta",CPreg="black"))+
  geom_hline(yintercept = 0.9,col="darkgrey",lwd=1)+
  annotate(x=5,y=0.9,label="0.9",geom="label")
 return(c)
}

QB2<-function(dat){
  d<-data.frame(dat)
  d$seg<-1:20
  names(d)<-c("DCP","DCPopt","CQR","CQRm","CQRr","CPreg","segment")
  c<-ggplot2::ggplot(d,aes(x=segment,y=DCP))+geom_line(aes(x=segment,y=DCP,
                                                           colour="DCP"),lwd=1,linetype="dashed")+
    geom_line(d,mapping=aes(x=segment,y=DCPopt,colour="DCPopt"),lwd=1,linetype="dashed")+
    geom_line(d,mapping=aes(x=segment,y=CQR,colour="CQR"),lwd=1,linetype="dashed")+
    geom_line(d,mapping=aes(x=segment,y=CQRm,colour="CQRm"),lwd=1,linetype="dashed")+
    geom_line(d,mapping=aes(x=segment,y=CQRr,colour="CQRr"),lwd=1,linetype="dashed")+
    geom_line(d,mapping=aes(x=segment,y=CPreg,colour="CPreg"),lwd=1,linetype="dashed")+
    xlab("Bin")+ylab("Conditional Length")+
    scale_colour_manual(name="Legend",
                        values=c(DCP="red", DCPopt="blue", CQR="darkgreen",
                                 CQRm="orange",CQRr="magenta",CPreg="black"))
  return(c)
}
  
QB(cov.cond.h)
QB(cov.cond.k)
QB(cov.cond.v)
QB(cov.cond.w)

QB2(leng.cond.h)
QB2(leng.cond.k)
QB2(leng.cond.v)
QB2(leng.cond.w)

