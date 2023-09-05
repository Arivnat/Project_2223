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
#####################EQUATIONS-43,45 AND 46####################################
###loading required packages








###Necessary libraries
###################################################

install.packages("readr", repos=c("http://rstudio.org/_packages",   "http://cran.rstudio.com"))
library(readr)
library(cowplot)
library(ggplot2)
library(quantreg)
library(hdm)
library(matrixcalc)





###########################################################
###############Function#######################################

sinerun<-function(n){
  alpha.sig<-0.1
  T.ho <- floor(n*0.2) ##20% holdout sample
  
  cov.mat <- leng.mat <- X.test.mat <- NULL
  cov.mat2 <- leng.mat2 <- X.test.mat2 <- NULL
  cov.mat3 <- leng.mat3 <- X.test.mat3 <- NULL
  
  for (i in 1:1000){
    set.seed(i)
    x<-runif(n,-1.5,1.5)
    f<-((x-1)^2)*(x+1)
    sig<-sqrt(0.25+abs(x))
    g<-2*sqrt(abs(x)+0.5)*as.numeric(x>=-0.05)
    #g2<-2*sqrt(abs(x+0.5))*as.numeric(x>=-0.05)
    y<-0.5*rnorm(n,f-g,sig)+0.5*rnorm(n,f+g,sig)
    e<-rnorm(n,1,2)
    y2<-0.5*rnorm(n,f-g,sig)+0.5*rnorm(n,f+g,sig)+e
    y3<-0.5*rnorm(n,f-g,sig)+0.5*rnorm(n,f+g,sig)+e*abs(x)
    x<-as.matrix(x)
    y<-as.matrix(y)
    y2<-as.matrix(y2)
    y3<-as.matrix(y3)
    
    
    sing <- TRUE
    sing2<-TRUE
    sing3<-TRUE
    while (sing==TRUE){
      while (sing2==TRUE){
        while (sing3==TRUE){
          
          ind <- sample(n,n,replace=FALSE)
          Y   <- y[ind]
          Y2   <- y2[ind]
          Y3<-y3[ind]
          X   <- x[ind,]
          X<-as.matrix(X)
          
          ind.test  <- (length(Y)-T.ho+1):length(Y)
          Y.test    <- Y[ind.test]
          X.test    <- cbind(X[ind.test,])
          
          ind.test2  <- (length(Y2)-T.ho+1):length(Y2)
          Y.test2    <- Y2[ind.test2]
          X.test2    <- cbind(X[ind.test2,])
          
          ind.test3  <- (length(Y3)-T.ho+1):length(Y3)
          Y.test3    <- Y3[ind.test3]
          X.test3    <- cbind(X[ind.test3,])
          
          obj.uneven <- aux.uneven(Y[-ind.test],cbind(X[-ind.test,]))
          obj.uneven2 <- aux.uneven(Y2[-ind.test2],cbind(X[-ind.test2,]))
          obj.uneven3 <- aux.uneven(Y3[-ind.test3],cbind(X[-ind.test3,]))
          
          Xa <- cbind(obj.uneven$X0)
          Ya <- obj.uneven$Y0
          X1a <- cbind(obj.uneven$X1)
          Y1a <- obj.uneven$Y1
          
          X02 <- cbind(obj.uneven2$X0)
          Y02 <- obj.uneven2$Y0
          X12 <- cbind(obj.uneven2$X1)
          Y12 <- obj.uneven2$Y1
          
          X03 <- cbind(obj.uneven3$X0)
          Y03 <- obj.uneven3$Y0
          X13 <- cbind(obj.uneven3$X1)
          Y13 <- obj.uneven3$Y1
          sing  <- is.singular.matrix(t(Xa)%*%Xa)
          sing2  <- is.singular.matrix(t(X02)%*%X02)
          sing3 <- is.singular.matrix(t(X03)%*%X03)
          
        }
      }
    }
    
    taus    <- seq(0.001,0.999,length=200)
    
    # Applying the difference conformal prediction methods
    res.qr    <- dcp.qr(Ya,Xa,Y1a,X1a,Y.test,X.test,taus,alpha.sig)
    res.qro<-dcp.opt(Ya,Xa,Y1a,X1a,Y.test,X.test,taus,alpha.sig)
    res.cqr   <- cqr(Ya,Xa,Y1a,X1a,Y.test,X.test,alpha.sig)
    res.reg   <- cp.reg(Ya,Xa,Y1a,X1a,Y.test,X.test,alpha.sig)
    
    res.qr2    <- dcp.qr(Y02,X02,Y12,X12,Y.test2,X.test2,taus,alpha.sig)
    res.qro2<-dcp.opt(Y02,X02,Y12,X12,Y.test2,X.test2,taus,alpha.sig)
    res.cqr2   <- cqr(Y02,X02,Y12,X12,Y.test2,X.test2,alpha.sig)
    res.reg2   <- cp.reg(Y02,X02,Y12,X12,Y.test2,X.test2,alpha.sig)
    
    res.qr3    <- dcp.qr(Y03,X03,Y13,X13,Y.test3,X.test3,taus,alpha.sig)
    res.qro3<-dcp.opt(Y03,X03,Y13,X13,Y.test3,X.test3,taus,alpha.sig)
    res.cqr3   <- cqr(Y03,X03,Y13,X13,Y.test3,X.test3,alpha.sig)
    res.reg3   <- cp.reg(Y03,X03,Y13,X13,Y.test3,X.test3,alpha.sig)
    
    X.test.mat <- rbind(X.test.mat,X.test)
    X.test.mat2 <- rbind(X.test.mat2,X.test2)
    X.test.mat3 <- rbind(X.test.mat3,X.test3)
    
    cov.mat.temp  <- cbind(res.qr$cov.qr,res.qro$cov.opt,res.cqr$cov.o,res.cqr$cov.m,res.cqr$cov.r,res.reg$cov.reg)
    leng.mat.temp <- cbind(res.qr$leng.qr,res.qro$leng.opt,res.cqr$leng.o,res.cqr$leng.m,res.cqr$leng.r,res.reg$leng.reg)
    
    cov.mat.temp2  <- cbind(res.qr2$cov.qr,res.qro2$cov.opt,res.cqr2$cov.o,res.cqr2$cov.m,res.cqr2$cov.r,res.reg2$cov.reg)
    leng.mat.temp2 <- cbind(res.qr2$leng.qr,res.qro2$leng.opt,res.cqr2$leng.o,res.cqr2$leng.m,res.cqr2$leng.r,res.reg2$leng.reg)
    
    cov.mat.temp3  <- cbind(res.qr3$cov.qr,res.qro3$cov.opt,res.cqr3$cov.o,res.cqr3$cov.m,res.cqr3$cov.r,res.reg3$cov.reg)
    leng.mat.temp3 <- cbind(res.qr3$leng.qr,res.qro3$leng.opt,res.cqr3$leng.o,res.cqr3$leng.m,res.cqr3$leng.r,res.reg3$leng.reg)
    
    cov.mat   <- rbind(cov.mat,cov.mat.temp)
    leng.mat   <- rbind(leng.mat,leng.mat.temp)
    
    
    cov.mat2   <- rbind(cov.mat2,cov.mat.temp2)
    leng.mat2  <- rbind(leng.mat2,leng.mat.temp2)
    
    cov.mat3   <- rbind(cov.mat3,cov.mat.temp3)
    leng.mat3  <- rbind(leng.mat3,leng.mat.temp3)
    
  }
  cov1 <- colMeans(cov.mat,na.rm=TRUE)
  cov2 <- colMeans(cov.mat2,na.rm=TRUE)
  cov3 <- colMeans(cov.mat3,na.rm=TRUE)
  leng1<-colMeans(leng.mat,na.rm=TRUE)
  leng2<-colMeans(leng.mat2,na.rm=TRUE)
  leng3<-colMeans(leng.mat3,na.rm=TRUE)
  
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
  
  pred.cov.qr3       <- predict(glm(cov.mat3[,1]~X.test.mat3,family=binomial(link="logit")),type="response")
  pred.cov.qr3.opt  <- predict(glm(cov.mat3[,2]~X.test.mat3,family=binomial(link="logit")),type="response")
  pred.cov.cqr3.o    <- predict(glm(cov.mat3[,3]~X.test.mat3,family=binomial(link="logit")),type="response")
  pred.cov.cqr3.m    <- predict(glm(cov.mat3[,4]~X.test.mat3,family=binomial(link="logit")),type="response")
  pred.cov.cqr3.r    <- predict(glm(cov.mat3[,5]~X.test.mat3,family=binomial(link="logit")),type="response")
  pred.cov.reg3      <- predict(glm(cov.mat3[,6]~X.test.mat3,family=binomial(link="logit")),type="response")
  
  
  pred.cov <- cbind(pred.cov.qr,pred.cov.qr.opt,pred.cov.cqr.o,pred.cov.cqr.m,pred.cov.cqr.r,pred.cov.reg)
  
  pred.cov2 <- cbind(pred.cov.qr2,pred.cov.qr2.opt,pred.cov.cqr2.o,pred.cov.cqr2.m,pred.cov.cqr2.r,pred.cov.reg2)
  
  pred.cov3 <- cbind(pred.cov.qr3,pred.cov.qr3.opt,pred.cov.cqr3.o,pred.cov.cqr3.m,pred.cov.cqr3.r,pred.cov.reg3)
  
  res.cov.cond <- apply(pred.cov,2,sd)
  
  res.cov.cond2 <- apply(pred.cov2,2,sd)
  
  res.cov.cond3 <- apply(pred.cov3,2,sd)
  
  l<-list(cov1,cov2,cov3,leng1,leng2,leng3,res.cov.cond,res.cov.cond2,res.cov.cond3,leng.mat,leng.mat2,leng.mat3,cov.mat,
          cov.mat2,cov.mat3,pred.cov,pred.cov2,pred.cov3,X.test.mat,X.test.mat2,X.test.mat3)
  return(l)
}

start<-Sys.time()
m<-sinerun(50)  ###2.368462 mins
Sys.time()-start

start<-Sys.time()
n<-sinerun(70)  ##2.964876 mins
Sys.time()-start

start<-Sys.time()
o<-sinerun(100)  ###4.17443 mins
Sys.time()-start

#start<-Sys.time()
#b<-sinerun(200)
#Sys.time()-start


start<-Sys.time()
p<-sinerun(300)  ###9.028689 mins
Sys.time()-start

start<-Sys.time()
q<-sinerun(500) ###24.13591 mins
Sys.time()-start


bind<-function(a,b,c){
  k<-cbind.data.frame(a,b,c)
  names(k)<-c("1","2","3")
  return(k)
}

formatn<-function(x){
  #if (x>0.09){y=x}
  #else {y=formatC(x,format="e",digits=2)}
  y<-ifelse(abs(x)>=0.01,round(x,4),formatC(x,format="e",digits=2))
  return(y)
}


comb<-function(df1,df2,df3,df4,df5){
  cov<-rbind.data.frame(bind(df1[1],df1[2],df1[3]),
                        bind(df2[1],df2[2],df2[3]),
                        bind(df3[1],df3[2],df3[3]),
                        bind(df4[1],df4[2],df4[3]),
                        bind(df5[1],df5[2],df5[3]))
  names(cov)<-c("c1","c2","c3")
  cov<-round(cov,3)
  
  leng<-rbind.data.frame(bind(df1[4],df1[5],df1[6]),
                         bind(df2[4],df2[5],df2[6]),
                         bind(df3[4],df3[5],df3[6]),
                         bind(df4[4],df4[5],df4[6]),
                         bind(df5[4],df5[5],df5[6]))
  names(leng)<-c("l1","l2","l3")
  leng<-round(leng,3)
  
  sd<-rbind.data.frame(bind(df1[7],df1[8],df1[9]),
                       bind(df2[7],df2[8],df2[9]),
                       bind(df3[7],df3[8],df3[9]),
                       bind(df4[7],df4[8],df4[9]),
                       bind(df5[7],df5[8],df5[9]))
  names(sd)<-c("(sd1)","(sd2)","(sd3)")
  sd[]<-lapply(sd,formatn)
  
  a<-cbind.data.frame(cov,leng,sd)
  a<-as.matrix(a)
  rownames(a)<-rep(c("DCP-QR","DCP-opt","CQR","CQR-m","CQR-r","CP-reg"),5)
  a2<-as.data.frame(a)
  l<-list(cov,leng,sd,a2)
  return(l)
  
}
xtable(data.frame(comb(m,n,o,p,q)[4]))

##############do plots###########################
n<-1000
set.seed(21)
x<-runif(n,-1.5,1.5)
f<-((x-1)^2)*(x+1)
sig<-sqrt(0.25+abs(x))
g<-2*sqrt(abs(x)+0.5)*as.numeric(x>=-0.05)
#g2<-2*sqrt(abs(x+0.5))*as.numeric(x>=-0.05)
y<-0.5*rnorm(n,f-g,sig)+0.5*rnorm(n,f+g,sig)
e<-rnorm(n,1,2)
y2<-0.5*rnorm(n,f-g,sig)+0.5*rnorm(n,f+g,sig)+e
y3<-0.5*rnorm(n,f-g,sig)+0.5*rnorm(n,f+g,sig)+e*abs(x)

g1<-ggplot(data=NULL,aes(x,y))+geom_point()+xlab("x")+
  ylab("y")
g2<-ggplot(data=NULL,aes(x,y2))+geom_point()+xlab("x")+
  ylab("y")
g3<-ggplot(data=NULL,aes(x,y3))+geom_point()+xlab("x")+
  ylab("y")

library(gridExtra)
grid.arrange(g1,g2,g3,ncol=3)

###################################################################
##################### QUANTILE BINS########################################

###FUNCTION FOR AVG.COVERAGE##################################################

###Obtain Xtest matrix from sincostandatagen folder
library(readxl)

s50xt<-data.frame(read_excel("s50.xlsx",sheet=5))
s50xt2<-data.frame(read_excel("s50.xlsx",sheet=11))
s50xt3<-data.frame(read_excel("s50.xlsx",sheet=17))

s70xt<-data.frame(read_excel("s70.xlsx",sheet=5))
s70xt2<-data.frame(read_excel("s70.xlsx",sheet=11))
s70xt3<-data.frame(read_excel("s70.xlsx",sheet=17))

s100xt<-data.frame(read_excel("s100.xlsx",sheet=5))
s100xt2<-data.frame(read_excel("s100.xlsx",sheet=11))
s100xt3<-data.frame(read_excel("s100.xlsx",sheet=17))

s300xt<-data.frame(read_excel("s300.xlsx",sheet=5))
s300xt2<-data.frame(read_excel("s300.xlsx",sheet=11))
s300xt3<-data.frame(read_excel("s300.xlsx",sheet=17))

s500xt<-data.frame(read_excel("s500.xlsx",sheet=5))
s500xt2<-data.frame(read_excel("s500.xlsx",sheet=11))
s500xt3<-data.frame(read_excel("s500.xlsx",sheet=17))


#################################################################################
################################################################################
####################FUNCTION FOR CONDITIONAL COVERAGE###############################
######################################################################################
qbc<-function(d,c,i,num.seg){  
  X<-as.numeric(unlist(d, use.names=FALSE))
  cov<-data.frame(c[i])
  cov.cond<-binning(X,cov,num.seg)$cond
  ###plotting cov.co
  
  cov.cond<-data.frame(cov.cond)
  c0<-cbind.data.frame(1:num.seg,cov.cond)
  
  names(c0)<-c("Bin","DCP-QR","DCP-opt","CQR","CQR-m","CQR-r","CP-reg")
  c2 <- melt(c0 ,  id.vars = 'Bin', variable.name = 'Method')
  
  n<-deparse(substitute(d))
  n2<-parse_number(n)
  #n3<-ifelse(n2=="2","20",n2)
  
  p<-ggplot(c2, aes(Bin, value)) +
    geom_line(aes(colour = Method),linetype="dashed",linewidth=1.25)+
    geom_hline(yintercept = 0.9,linewidth=1)+scale_colour_brewer(palette="Dark2") +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    #ggtitle(paste("Sample Size:",n2))+
    xlab("Quantile Bin")+
    ylab("Conditional Coverage")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
  k<-list(p,cov.cond)
  return(p)
}

###Sample-50
q0<-qbc(m[19],m,13,10) ###specification 1
r0<-qbc(m[20],m,14,10) ###specification 2
s0<-qbc(m[21],m,15,10)###specification 3


###Sample-70
q1<-qbc(n[19],n,13,10) ###specification 1
r1<-qbc(n[20],n,14,10) ###specification 2
s1<-qbc(n[21],n,15,10)###specification 3


###Sample-100
q2<-qbc(o[19],o,13,10) ###specification 1
r2<-qbc(o[20],o,14,10) ###specification 2
s2<-qbc(o[21],o,15,10)###specification 3

###Sample-300
q3<-qbc(p[19],p,13,10) ###specification 1
r3<-qbc(p[20],p,14,10) ###specification 2
s3<-qbc(p[21],p,15,10)###specification 3

####Sample-500
q4<-qbc(q[19],q,13,10) ###specification 1
r4<-qbc(q[20],q,14,10) ###specification 2
s4<-qbc(q[21],q,15,10)###specification 3

plot_grid(q0,r0,s0,nrow=1)
plot_grid(q1,r1,s1,nrow=1)
plot_grid(q2,r2,s2,nrow=1)
plot_grid(q3,r3,s3,nrow=1)
plot_grid(q4,r4,s4,nrow=1)

#################################################################################
##FUNCTION FOR CONDITIONAL LENGTH#########################################################
qbl<-function(d,l,i,num.seg){  
  X<-as.numeric(unlist(d, use.names=FALSE))
  leng<-data.frame(l[i])
  leng.cond <- binning(X,leng,num.seg)$cond
  leng.cond<-data.frame(leng.cond)
  c<-cbind.data.frame(1:num.seg,leng.cond)
  names(c)<-c("Bin","DCP-QR","DCP-opt","CQR","CQR-m","CQR-r","CP-reg")
  c2 <- melt(c ,  id.vars = 'Bin', variable.name = 'Method')
  n<-deparse(substitute(d))
  n2<-parse_number(n)
  #n3<-ifelse(n2=="2","20",n2)
  p<-ggplot(c2, aes(Bin, value)) +
    geom_line(aes(colour = Method),linewidth=1.25)+
    scale_colour_brewer(palette="Dark2") +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    # ggtitle(paste("Sample Size:",n2))+
    xlab("Quantile Bin")+
    ylab("Conditional Length")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
  k<-list(p,leng.cond)
  return(p)
}

##################################################################
################################################################################
q0l<-qbl(m[19],m,10,10)
r0l<-qbl(m[20],m,11,10)
s0l<-qbl(m[21],m,12,10)


q1l<-qbl(n[19],n,10,10)
r1l<-qbl(n[20],n,11,10)
s1l<-qbl(n[21],n,12,10)

q2l<-qbl(o[19],o,10,10)
r2l<-qbl(o[20],o,11,10)
s2l<-qbl(o[21],o,12,10)

q3l<-qbl(p[19],p,10,10)
r3l<-qbl(p[20],p,11,10)
s3l<-qbl(p[21],p,12,10)

q4l<-qbl(q[19],q,10,10)
r4l<-qbl(q[20],q,11,10)
s4l<-qbl(q[21],q,12,10)


plot_grid(q0l,r0l,s0l,nrow=1)
plot_grid(q1l,r1l,s1l,nrow=1)
plot_grid(q2l,r2l,s2l,nrow=1)
plot_grid(q3l,r3l,s3l,nrow=1)
plot_grid(q4l,r4l,s4l,nrow=1)


#q0l<-qbl(s50xt,m[10],10)[1]
#r0l<-qbl(s50xt2,m[11],10)[1]
#s0l<-qbl(s50xt3,m[12],10)[1]


#q1l<-qbl(s70xt,n[10],10)[1]
#r1l<-qbl(s70xt2,n[11],10)[1]
#s1l<-qbl(s70xt3,n[12],10)[1]

#q2l<-qbl(s100xt,o[10],10)[1]
#2l<-qbl(s100xt2,o[11],10)[1]
#s2l<-qbl(s100xt3,o[12],10)[1]

#q3l<-qbl(s300xt,p[10],10)[1]
#r3l<-qbl(s300xt2,p[11],10)[1]
#s3l<-qbl(s300xt3,p[12],10)[1]

#q4l<-qbl(s500xt,q[10],10)[1]
#r4l<-qbl(s500xt2,q[11],10)[1]
#s4l<-qbl(s500xt3,q[12],10)[1]
################################################################################
################################################################################

