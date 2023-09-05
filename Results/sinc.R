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



library(rdetools)
library(xtable)
sine<-function(n){
  alpha.sig<-0.1
  T.ho <- floor(n*0.2)
  cov.mat <- leng.mat <- X.test.mat <- NULL
  cov.mat2 <- leng.mat2 <- X.test.mat2 <- NULL
  cov.mat0 <- leng.mat0 <- X.test.mat0 <- NULL
  cov.mat002 <- leng.mat002 <- X.test.mat002 <- NULL
  
  for (i in 1:1000){
    set.seed(i)
    x1<-runif(n,-3,3)
    x2<-rnorm(n,1,2)
    y0<-sinc(x1)+rnorm(n,0,1)
    y02<-sinc(x1)+rcauchy(n,1,3)
    y<-sinc(x1)+2*x2+rnorm(n,0,1)
    y2<-sinc(x1)+2*x2+rcauchy(n,1,3)
    x<-cbind.data.frame(x1,x2)
    x<-as.matrix(x)
    
    sing<-TRUE
    sing2<-TRUE
    sing3<-TRUE
    sing4<-TRUE
    while (sing==TRUE){
      while (sing2==TRUE){
        while (sing3==TRUE){
          while (sing4==TRUE){
      ind <- sample(n,n,replace=FALSE)
      Y0<-y0[ind]
      Y002<-y02[ind]
      Y   <- y[ind]
      Y2   <- y2[ind]
      X   <- x[ind,]
      X<-as.matrix(X)
      X10<-as.matrix(x1)
      X10t<-as.matrix(X10[ind,])
      
      ind.test0  <- (length(Y0)-T.ho+1):length(Y0)
      Y.test0    <- Y0[ind.test0]
      X.test0    <- cbind(X10t[ind.test0,])
      
      ind.test002  <- (length(Y002)-T.ho+1):length(Y002)
      Y.test002    <- Y002[ind.test002]
      X.test002    <- cbind(X10t[ind.test002,])
      
      ind.test  <- (length(Y)-T.ho+1):length(Y)
      Y.test    <- Y[ind.test]
      X.test    <- cbind(X[ind.test,])
      
      ind.test2  <- (length(Y2)-T.ho+1):length(Y2)
      Y.test2    <- Y2[ind.test2]
      X.test2    <- cbind(X[ind.test2,])
      
      obj.uneven0 <- aux.uneven(Y0[-ind.test0],cbind(X10t[-ind.test0,]))
      obj.uneven002 <- aux.uneven(Y002[-ind.test002],cbind(X10t[-ind.test002,]))
      obj.uneven <- aux.uneven(Y[-ind.test],cbind(X[-ind.test,]))
      obj.uneven2 <- aux.uneven(Y2[-ind.test2],cbind(X[-ind.test2,]))
      
      X0a <- cbind(obj.uneven0$X0)
      Y0a <- obj.uneven0$Y0
      X01a <- cbind(obj.uneven0$X1)
      Y01a <- obj.uneven0$Y1
      
      X02a <- cbind(obj.uneven002$X0)
      Y02a <- obj.uneven002$Y0
      X012a <- cbind(obj.uneven002$X1)
      Y012a <- obj.uneven002$Y1
      
      Xa <- cbind(obj.uneven$X0)
      Ya <- obj.uneven$Y0
      X1a <- cbind(obj.uneven$X1)
      Y1a <- obj.uneven$Y1
      
      X02 <- cbind(obj.uneven2$X0)
      Y02 <- obj.uneven2$Y0
      X12 <- cbind(obj.uneven2$X1)
      Y12 <- obj.uneven2$Y1
      
      sing  <- is.singular.matrix(t(X0a)%*%X0a)
      sing2  <- is.singular.matrix(t(X02a)%*%X02a)
      sing3  <- is.singular.matrix(t(Xa)%*%Xa)
      sing4  <- is.singular.matrix(t(X02)%*%X02)
      }      
        }
      }
    }

  taus    <- seq(0.001,0.999,length=200)
  
  # Applying the difference conformal prediction methods
  res.qr0    <- dcp.qr(Y0a,X0a,Y01a,X01a,Y.test0,X.test0,taus,alpha.sig)
  res.qro0<-dcp.opt(Y0a,X0a,Y01a,X01a,Y.test0,X.test0,taus,alpha.sig)
  res.cqr0   <- cqr(Y0a,X0a,Y01a,X01a,Y.test0,X.test0,alpha.sig)
  res.reg0   <- cp.reg(Y0a,X0a,Y01a,X01a,Y.test0,X.test0,alpha.sig)
  
  res.qr002    <- dcp.qr(Y02a,X02a,Y012a,X012a,Y.test002,X.test002,taus,alpha.sig)
  res.qro002<-dcp.opt(Y02a,X02a,Y012a,X012a,Y.test002,X.test002,taus,alpha.sig)
  res.cqr002   <- cqr(Y02a,X02a,Y012a,X012a,Y.test002,X.test002,alpha.sig)
  res.reg002   <- cp.reg(Y02a,X02a,Y012a,X012a,Y.test002,X.test002,alpha.sig)
  
  res.qr    <- dcp.qr(Ya,Xa,Y1a,X1a,Y.test,X.test,taus,alpha.sig)
  res.qro<-dcp.opt(Ya,Xa,Y1a,X1a,Y.test,X.test,taus,alpha.sig)
  res.cqr   <- cqr(Ya,Xa,Y1a,X1a,Y.test,X.test,alpha.sig)
  res.reg   <- cp.reg(Ya,Xa,Y1a,X1a,Y.test,X.test,alpha.sig)
  
  res.qr2    <- dcp.qr(Y02,X02,Y12,X12,Y.test2,X.test2,taus,alpha.sig)
  res.qro2<-dcp.opt(Y02,X02,Y12,X12,Y.test2,X.test2,taus,alpha.sig)
  res.cqr2   <- cqr(Y02,X02,Y12,X12,Y.test2,X.test2,alpha.sig)
  res.reg2   <- cp.reg(Y02,X02,Y12,X12,Y.test2,X.test2,alpha.sig)
  
  X.test.mat0 <- rbind(X.test.mat0,X.test0)
  X.test.mat002 <- rbind(X.test.mat002,X.test002)
  
  X.test.mat <- rbind(X.test.mat,X.test)
  X.test.mat2 <- rbind(X.test.mat2,X.test2)
  
  cov.mat.temp0  <- cbind(res.qr0$cov.qr,res.qro0$cov.opt,res.cqr0$cov.o,res.cqr0$cov.m,res.cqr0$cov.r,res.reg0$cov.reg)
  leng.mat.temp0 <- cbind(res.qr0$leng.qr,res.qro0$leng.opt,res.cqr0$leng.o,res.cqr0$leng.m,res.cqr0$leng.r,res.reg0$leng.reg)
  
  cov.mat.temp002  <- cbind(res.qr002$cov.qr,res.qro002$cov.opt,res.cqr002$cov.o,res.cqr002$cov.m,res.cqr002$cov.r,res.reg002$cov.reg)
  leng.mat.temp002 <- cbind(res.qr002$leng.qr,res.qro002$leng.opt,res.cqr002$leng.o,res.cqr002$leng.m,res.cqr002$leng.r,res.reg002$leng.reg)
 
  cov.mat.temp  <- cbind(res.qr$cov.qr,res.qro$cov.opt,res.cqr$cov.o,res.cqr$cov.m,res.cqr$cov.r,res.reg$cov.reg)
  leng.mat.temp <- cbind(res.qr$leng.qr,res.qro$leng.opt,res.cqr$leng.o,res.cqr$leng.m,res.cqr$leng.r,res.reg$leng.reg)
  
  cov.mat.temp2  <- cbind(res.qr2$cov.qr,res.qro2$cov.opt,res.cqr2$cov.o,res.cqr2$cov.m,res.cqr2$cov.r,res.reg2$cov.reg)
  leng.mat.temp2 <- cbind(res.qr2$leng.qr,res.qro2$leng.opt,res.cqr2$leng.o,res.cqr2$leng.m,res.cqr2$leng.r,res.reg2$leng.reg)
  
  cov.mat0   <- rbind(cov.mat0,cov.mat.temp0)
  leng.mat0  <- rbind(leng.mat0,leng.mat.temp0)
  
  cov.mat002   <- rbind(cov.mat002,cov.mat.temp002)
  leng.mat002  <- rbind(leng.mat002,leng.mat.temp002)
  
  cov.mat   <- rbind(cov.mat,cov.mat.temp)
  leng.mat  <- rbind(leng.mat,leng.mat.temp)
  
  cov.mat2   <- rbind(cov.mat2,cov.mat.temp2)
  leng.mat2  <- rbind(leng.mat2,leng.mat.temp2)
  
  
  }
  
  l<-list(cov.mat0,cov.mat002,cov.mat,cov.mat2,
          leng.mat0,leng.mat002,leng.mat,leng.mat2,
          X.test.mat0,X.test.mat002,X.test.mat,X.test.mat2)
  return(l)
}

start<-Sys.time()
s00<-sine(50)  ##3.664296 mins
Sys.time()-start

##Time difference of 

start<-Sys.time()
s0<-sine(70)
Sys.time()-start
#Time difference of 3.547034 mins
#######################################
start<-Sys.time()
s<-sine(100)  ###4.542645 mins
Sys.time()-start 

start<-Sys.time()
s2<-sine(200)   ###9.829671 mins
Sys.time()-start  

start<-Sys.time()
s3<-sine(300)  ###11.22557 mins
Sys.time()-start  

start<-Sys.time()
s4<-sine(500)  ##21.71108 mins
Sys.time()-start

start<-Sys.time()
s5<-sine(700)  ###Time difference of 14.71807 mins
Sys.time()-start

start<-Sys.time()
s6<-sine(1000)
Sys.time()-start   ###Time difference of 17.74256 mins

##########################################################################
###getting avg coverage and avg length by method
#######################################################################
##s00,s0,s,s2,s3,s4
coverage<-function(d){
  cov1<-colMeans(data.frame(d[1]),na.rm=TRUE)
  cov2<-colMeans(data.frame(d[2]),na.rm=TRUE)
  cov3<-colMeans(data.frame(d[3]),na.rm=TRUE)
  cov4<-colMeans(data.frame(d[4]),na.rm=TRUE)
  C<-cbind.data.frame(cov1,cov2,cov3,cov4)
  return(C)
}

cov<-rbind.data.frame(coverage(s00),coverage(s0),coverage(s),coverage(s2),coverage(s3),coverage(s4))

len<-function(d){
  len1<-colMeans(data.frame(d[5]),na.rm=TRUE)
  len2<-colMeans(data.frame(d[6]),na.rm=TRUE)
  len3<-colMeans(data.frame(d[7]),na.rm=TRUE)
  len4<-colMeans(data.frame(d[8]),na.rm=TRUE)
  L<-cbind.data.frame(len1,len2,len3,len4)
  return(L)
}

leng<-rbind.data.frame(len(s00),len(s0),len(s),len(s2),len(s3),len(s4))


sd<-function(d){
  SD<-NULL
  for (i in 1:4){
  d1<-as.matrix(data.frame(d[i]))
  e<-as.matrix(data.frame(d[i+8]))
  pred.cov.qr       <- predict(glm(d1[,1]~e,family=binomial(link="logit")),type="response")
  pred.cov.qr.opt   <- predict(glm(d1[,2]~e,family=binomial(link="logit")),type="response")
  pred.cov.cqr.o    <- predict(glm(d1[,3]~e,family=binomial(link="logit")),type="response")
  pred.cov.cqr.m    <- predict(glm(d1[,4]~e,family=binomial(link="logit")),type="response")
  pred.cov.cqr.r    <- predict(glm(d1[,5]~e,family=binomial(link="logit")),type="response")
  pred.cov.reg      <- predict(glm(d1[,6]~e,family=binomial(link="logit")),type="response")
  p<- cbind(pred.cov.qr,pred.cov.qr.opt,pred.cov.cqr.o,pred.cov.cqr.m,pred.cov.cqr.r,pred.cov.reg)
  res.cov.cond <- apply(p,2,function(x){sqrt(var(x))})
  SDtemp<-res.cov.cond
  SD<-rbind(SD,SDtemp)
  }
 return(SD)
  
}

ss<-rbind.data.frame(t(sd(s00)),t(sd(s0)),t(sd(s)),t(sd(s2)),t(sd(s3)),t(sd(s4)))

CLS<-cbind.data.frame(cov,leng,ss)

###Apply format for latex
formatn<-function(x){
  #if (x>0.09){y=x}
  #else {y=formatC(x,format="e",digits=2)}
  y<-ifelse(abs(x)>=0.01,round(x,4),formatC(x,format="e",digits=2))
  return(y)
}
CLS[]<-lapply(CLS,formatn)
xtable(CLS,digits=4)

formatn<-function(x){
  #if (x>0.09){y=x}
  #else {y=formatC(x,format="e",digits=2)}
  y<-ifelse(abs(x)>=0.01,round(x,4),formatC(x,format="e",digits=2))
  return(y)
}

#######################
library(openxlsx)
write.xlsx(s00,"s00.xlsx")
write.xlsx(s0,"s0.xlsx")
write.xlsx(s,"s.xlsx")
write.xlsx(s2,"s2.xlsx")
write.xlsx(s3,"s3.xlsx")
write.xlsx(s4,"s4.xlsx")

write.xlsx(CLS,"CLS.xlsx")
#########################################

############QUANTILE BINS#######################

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
  
  
  names(c1)<-names(c2)<-c("Bin","DCP-QR","DCP-opt","CQR","CQR-m","CQR-r","CP-reg")
  
  c11<- melt(c1 ,  id.vars = 'Bin', variable.name = 'Method')
  c22<- melt(c2 ,  id.vars = 'Bin', variable.name = 'Method')
 
  
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

  
  p2_legend <-ggplot(c22, aes(Bin, value)) +
    geom_line(aes(colour = Method),linetype="dashed",linewidth=1.25)+
    geom_hline(yintercept = 0.9,linewidth=1)+scale_colour_brewer(palette="Dark2") +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    ggtitle(paste("Sample Size:",n2,",Specification:",spec2))+xlab("Quantile Bin")+
    ylab("Conditional Coverage")+theme(legend.position = "right")
  
  shared_legend <- extract_legend(p2_legend)
  
  grid.arrange(arrangeGrob(p1,p2,ncol=2),
               shared_legend)
  
  k<-list(p1,p2,c1,c2)
  return(k)
}


qbc(s00,9,1,10)
qbc(s0,9,1,10)
qbc(s,9,1,10)
qbc(s2,9,1,10)
qbc(s3,9,1,10)
qbc(s4,9,1,10)

###################################################################################
#######################AVERAGE LENGTH###############################################
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
  
  
  names(c1)<-names(c2)<-c("Bin","DCP-QR","DCP-opt","CQR","CQR-m","CQR-r","CP-reg")
  
  c11<- melt(c1 ,  id.vars = 'Bin', variable.name = 'Method')
  c22<- melt(c2 ,  id.vars = 'Bin', variable.name = 'Method')
 
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
  
  p2_legend <-ggplot(c22, aes(Bin, value)) +
    geom_line(aes(colour = Method),linewidth=1.25)+
    scale_colour_brewer(palette="Dark2") +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    ggtitle(paste("Sample Size:",n2,",Specification:",spec2))+xlab("Quantile Bin")+
    ylab("Conditional Length")+theme(legend.position = "right")
  
  shared_legend <- extract_legend(p2_legend)
  
  grid.arrange(arrangeGrob(p1,p2,ncol=2),
               shared_legend)
  
  k<-list(p1,p2,c1,c2)
  return(k)
}

qbl(s00,9,5,10)
qbl(s0,9,5,10)
qbl(s,9,5,10)
qbl(s2,9,5,10)
qbl(s3,9,5,10)
qbl(s4,9,5,10)


###############################################################
################VISUALIZING THE DATA ##############################

set.seed(500)
n<-300
t<-runif(n,-3,3)
v<-rnorm(n,1,2)
w<-sinc(x1)+rnorm(n,0,1)
w2<-sinc(x1)+rcauchy(n,1,3)
dz<-sinc(x1)+2*x2+rnorm(n,0,1)
dz2<-sinc(x1)+2*x2+rcauchy(n,1,3)

ggplot(data=NULL,aes(t,w))+geom_point()+xlab("x1")+ylab("y")
ggplot(data=NULL,aes(t,w2))+geom_point()+xlab("x1")+ylab("y")
plot_ly(x=t, y=v, z=dz, type="scatter3d", mode="markers", color=dz)
plot_ly(x=t, y=v, z=dz2, type="scatter3d", mode="markers", color=dz2)



  avgcov<-colMeans(data.frame(dat[k]),na.rm=TRUE)
  avgleng<-colMeans(data.frame(dat[l]),na.rm=TRUE)
  sdcov<-data.frame(dat[s])
  res<-cbind.data.frame(avgcov,avgleng,sdcov)
  names(res)<-c("Avg.Coverage","Avg. Length","SD")
  res[]<-lapply(res,formatn)
  return(res)
}

samp50a<-results(s00,9,7,5)
samp50b<-results(s00,10,8,6)

xtable(samp50a)
xtable(samp50b)

xtable(t(samp50a))
xtable(t(samp50b))

samp70a<-results(s0,9,7,5)
samp70b<-results(s0,10,8,6)

xtable(samp70a)
xtable(samp70b)

xtable(t(samp70a))
xtable(t(samp70b))

samp100a<-results(s,9,7,5)
samp100b<-results(s,10,8,6)

xtable(samp100a)
xtable(samp100b)

xtable(t(samp100a))
xtable(t(samp100b))

samp200a<-results(s2,9,7,5)
samp200b<-results(s2,10,8,6)

xtable(samp200a)
xtable(samp200b)

xtable(t(samp200a))
xtable(t(samp200b))

samp300a<-results(s3,9,7,5)
samp300b<-results(s3,10,8,6)

xtable(samp300a)
xtable(samp300b)

xtable(t(samp300a))
xtable(t(samp300b))

samp500a<-results(s4,9,7,5)
samp500b<-results(s4,10,8,6)

xtable(samp500a)
xtable(samp500b)

xtable(t(samp500a))
xtable(t(samp500b))

samp700a<-results(s5,9,7,5)
samp700b<-results(s5,10,8,6)

xtable(t(samp700a))
xtable(t(samp700b))

samp1000a<-results(s6,9,7,5)
samp1000b<-results(s6,10,8,6)

xtable(t(samp1000a))
xtable(t(samp1000b))

##########which list elements to use?
###for first case use: 9,7,5
##for second case use: 10,8,6
###get lengmat,covmat to get per iteration results
###################################################################
#cov,cov2,leng,leng2,res.cov.cond,res.cov.cond2,leng.mat,leng.mat2,cov.mat,
#cov.mat2,pred.cov.qr,pred.cov.qr.opt,pred.cov.cqr.o,pred.cov.cqr.m,pred.cov.cqr.r,pred.cov.reg,
#pred.cov.qr2,pred.cov.qr2.opt,pred.cov.cqr2.o,pred.cov.cqr2.m,pred.cov.cqr2.r,pred.cov.reg2

##s00-leng 7,8, cov 9,10 
empcov<-function(n,m,s){
  sv<-NULL
  for (i in 1:s){
    subval<-m[n*(i-1)+1:n*i,]
    subval_avg<-colMeans(subval,na.rm=TRUE)
    sv<-rbind(sv,subval_avg)
  }
  return(sv)
}

###find empirical avg coverage per iteration

###get files from Python 
plotfunc<-function(k,dat,dat2,dat3,n,it){
  x<-data.frame(dat[k])
  x<-empcov(n,x,it)
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

##Use it to create boxplots for avg. coverage and avg.length per iteration


sink("s.txt")
s
sink()

sink("s2.txt")
s2
sink()

sink("s3.txt")
s3
sink()

sink("s4.txt")
s4
sink()

sink("s5.txt")
s5
sink()

sink("s6.txt")
s6
sink()


##Where are these sink files?
##Also need to find location for sinc in python
##Need to export files to Excel and from Excel to Python


##You also need the X0,X1 and Xtest files for Python

sinep<-function(n){
  alpha.sig<-0.1
  T.ho <- floor(n*0.2)
  cov.mat <- leng.mat <- X.test.mat <- X0.mat<-Y0.mat<-X1.mat<-Y1.mat<-Y.test.mat<-NULL
  cov.mat2 <- leng.mat2 <- X.test.mat2 <- X02.mat<-Y02.mat<-X12.mat<-Y12.mat<-Y.test.mat2<-NULL
  
  for (i in 1:1000){
    set.seed(i)
    x1<-runif(n,-3,3)
    x2<-rnorm(n,1,2)
    y<-sinc(x1)+2*x2+rnorm(n,0,1)
    y2<-sinc(x1)+2*x2+rcauchy(n,1,3)
    x<-cbind.data.frame(x1,x2)
    x<-as.matrix(x)
    
    sing<-TRUE
    sing2<-TRUE
    while (sing==TRUE){
      while (sing2==TRUE){
        ind <- sample(n,n,replace=FALSE)
        Y   <- y[ind]
        Y2   <- y2[ind]
        X   <- x[ind,]
        X<-as.matrix(X)
        
        ind.test  <- (length(Y)-T.ho+1):length(Y)
        Y.test    <- Y[ind.test]
        X.test    <- cbind(X[ind.test,])
        
        ind.test2  <- (length(Y2)-T.ho+1):length(Y2)
        Y.test2    <- Y2[ind.test2]
        X.test2    <- cbind(X[ind.test2,])
        
        
        obj.uneven <- aux.uneven(Y[-ind.test],cbind(X[-ind.test,]))
        obj.uneven2 <- aux.uneven(Y2[-ind.test2],cbind(X[-ind.test2,]))
        
        X0 <- cbind(obj.uneven$X0)
        Y0 <- obj.uneven$Y0
        X1 <- cbind(obj.uneven$X1)
        Y1 <- obj.uneven$Y1
        
        X02 <- cbind(obj.uneven2$X0)
        Y02 <- obj.uneven2$Y0
        X12 <- cbind(obj.uneven2$X1)
        Y12 <- obj.uneven2$Y1
        
        sing  <- is.singular.matrix(t(X0)%*%X0)
        sing2  <- is.singular.matrix(t(X02)%*%X02)
      }      
    }
    
    
    X.test.mat <- rbind(X.test.mat,X.test)
    X.test.mat2 <- rbind(X.test.mat2,X.test2)
    
    X0.mat <- rbind(X0.mat,X0)
    X1.mat <- rbind(X1.mat,X1)
    
    Y0.mat <- rbind(Y0.mat,Y0)
    Y1.mat<- rbind(Y1.mat,Y1)
    
    X02.mat <- rbind(X02.mat,X02)
    X12.mat <- rbind(X12.mat,X12)
    
    Y02.mat <- rbind(Y02.mat,Y02)
    Y12.mat<- rbind(Y12.mat,Y12)
    
    Y.test.mat <- rbind(Y.test.mat,Y.test)
    Y.test.mat2 <- rbind(Y.test.mat2,Y.test2)
    
  }
  
  
  l<-list(X0.mat,Y0.mat,X1.mat,Y1.mat,X02.mat,Y02.mat,X12.mat,Y12.mat,
      X.test.mat,X.test.mat2,Y.test.mat,Y.test.mat2)
  return(l)
}

a<-sinep(50)
b<-sinep(70)
c<-sinep(100)
d<-sinep(200)
e<-sinep(300)
f<-sinep(500)
g<-sinep(700)
h<-sinep(1000)


###Exporting to Excel for Python 
eaX0<-data.frame(a[1]) ###20000 2
eaY0<-data.frame(a[2])
eaX1<-data.frame(a[3])
eaY1<-data.frame(a[4])
eaX02<-data.frame(a[5])
eaY02<-data.frame(a[6])
eaX12<-data.frame(a[7])
eaY12<-data.frame(a[8])
eaXt<-data.frame(a[9])
eaXt2<-data.frame(a[10])
eaYt<-data.frame(a[11])
eaYt2<-data.frame(a[12])
###########################
list_exa <- list("X0" = eaX0, "Y0" = eaY0,"X1"=eaX1,"Y1"=eaY1, "X02"=eaX02,
                 "Y02"=eaY02,"X12"=eaX12,"Y12"=eaY12,
                 "Xt"=eaXt,"Xt2"=eaXt2,"Yt"=eaYt,"Yt2"=eaYt2)
write.xlsx(list_exa, file = "Ex1.xlsx")

exportfunc<-function(item){
  list_ex<-list()
  for (i in 1:12){
        k=paste0(i)
        h=list(k=data.frame(item[i]))
  list_ex=append(list_ex,h)
  #names(list_ex)<-names(list_exa)
  #j<-write.xlsx(list_ex, file ="Exa.xlsx" )
  }
 
  return(list_ex)
}
a1<-exportfunc(a)
names(a1)<-names(list_exa)
write.xlsx(a1, file = "Ex50.xlsx")

b1<-exportfunc(b)
names(b1)<-names(list_exa)
write.xlsx(b1, file = "Ex70.xlsx")



#exportxl<-function(item){
 # listn=c(50,70,100,200,300,500,700,1000)
  #for (i in listn){
   # item1<-exportfunc(item)
    #names(item1)<-names(list_exa)
    #j<-write.xlsx(item1,file=paste0("ex",i,"t.xlsx"))
    
#}
#return(j)
#}

c1<-exportfunc(c)
names(c1)<-names(list_exa)
write.xlsx(b1, file = "Ex100.xlsx")

d1<-exportfunc(d)
names(d1)<-names(list_exa)
write.xlsx(d1, file = "Ex200.xlsx")

e1<-exportfunc(e)
names(e1)<-names(list_exa)
write.xlsx(e1, file = "Ex300.xlsx")

f1<-exportfunc(f)
names(f1)<-names(list_exa)
write.xlsx(f1, file = "Ex500.xlsx")

g1<-exportfunc(g)
names(g1)<-names(list_exa)
write.xlsx(g1, file = "Ex700.xlsx")

h1<-exportfunc(h)
names(h1)<-names(list_exa)
write.xlsx(h1, file = "Ex1000.xlsx")


####################################################
##getting results from R output



