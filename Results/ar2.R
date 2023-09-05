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

#############################################################
################AR(p=2)##########################################

install.packages("ggthemes")
library(ggthemes)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(quantreg)
library(matrixcalc)
library(hdm)
au2<-function(n,j,k,t,a,b){
  alpha.sig<-0.1
  cov.mat<-leng.mat<-X.test.mat<-Y.test.mat<-NULL
  X0.mat<-Y0.mat<-X1.mat<-Y1.mat<-NULL
  
  for (i in 1:t){
    g<-seq(1,n,1)
    p<-0.8*n
    indt<-g[p+1:n]
    indt<-indt[!is.na(indt)]
    indp<-g[-indt]
    
    set.seed(i)
    x1<-rnorm(n,2,6)
    x2<-rnorm(n,1.5,3.25)
    err<-15*arima.sim(n,model=list(ar=c(j,k)),innov=runif(n,a,b))
    y<-4.5+1.5*x1+2.5*x2+err
    y<-data.frame(y)
    x<-cbind.data.frame(x1,x2)
    
    Y.test<-y[indt,]
    X.test<-cbind(x[indt,])
    obj.uneven<-aux.uneven(y[indp,],cbind(x[indp,]))
    Xa<-cbind(obj.uneven$X0)
    Ya<-obj.uneven$Y0
    X1a<-cbind(obj.uneven$X1)
    Y1a<-obj.uneven$Y1
    Xa<-as.matrix(Xa)
    Ya<-as.matrix(Ya)
    X1a<-as.matrix(X1a)
    Y1a<-as.matrix(Y1a)
    X.test<-as.matrix(X.test)
    Y.test<-as.matrix(Y.test)
    
    taus    <- seq(0.001,0.999,length=200)  ##change taus to avoid infinite values 
    #ys      <- quantile(unique(c(Y0,Y1)),seq(0.001,0.999,length=200))
    
    # Applying the difference conformal prediction methods
    res.qr    <- dcp.qr(Ya,Xa,Y1a,X1a,Y.test,X.test,taus,alpha.sig)
    res.opt   <- dcp.opt(Ya,Xa,Y1a,X1a,Y.test,X.test,taus,alpha.sig)
    res.cqr   <- cqr(Ya,Xa,Y1a,X1a,Y.test,X.test,alpha.sig)
    res.reg   <- cp.reg(Ya,Xa,Y1a,X1a,Y.test,X.test,alpha.sig)
    
    # Results
    X.test.mat <- rbind(X.test.mat,X.test)
    Y.test.mat <- rbind(Y.test.mat,Y.test)
    X0.mat <- rbind(X0.mat,Xa)
    Y0.mat <- rbind(Y0.mat,Ya)
    X1.mat <- rbind(X1.mat,X1a)
    Y1.mat <- rbind(Y1.mat,Y1a)
    
    
    cov.mat.temp  <- cbind(res.qr$cov.qr,res.opt$cov.opt,res.cqr$cov.o,res.cqr$cov.m,res.cqr$cov.r,res.reg$cov.reg)
    leng.mat.temp <- cbind(res.qr$leng.qr,res.opt$leng.opt,res.cqr$leng.o,res.cqr$leng.m,res.cqr$leng.r,res.reg$leng.reg)
    
    cov.mat   <- rbind(cov.mat,cov.mat.temp)
    leng.mat  <- rbind(leng.mat,leng.mat.temp)
    
  }
  
  l<-list(cov.mat,leng.mat,X0.mat,Y0.mat,X1.mat,Y1.mat,X.test.mat,Y.test.mat,j,k)
  return(l)
}

autor<-list()
u<-list(0.005,0.1,0.3,0.5,0.7,0.9)
u2<--0.1  ###-0.1 0.02 0.09
  
 # s<-Sys.time()-start
  
  alu2<-function(u,v,n,t,a,b){
    start<-Sys.time()
    for (i in 1:length(u)){
      autor[[i]]<-au2(n,u[[i]],u2,t,a,b)
    }
    j<-Sys.time()-start
    c1<-colMeans(data.frame(autor[[1]][1]),na.rm=TRUE)
    c2<-colMeans(data.frame(autor[[2]][1]),na.rm=TRUE)
    c3<-colMeans(data.frame(autor[[3]][1]),na.rm=TRUE)
    c4<-colMeans(data.frame(autor[[4]][1]),na.rm=TRUE)
    c5<-colMeans(data.frame(autor[[5]][1]),na.rm=TRUE)
    c6<-colMeans(data.frame(autor[[6]][1]),na.rm=TRUE)
    
    C<-data.frame(rbind(c1,c2,c3,c4,c5,c6))
    C$AR<-c(0.005,0.1,0.3,0.5,0.7,0.9)
    names(C)<-c("DCP-QR","DCP-opt","CQR","CQR-m","CQR-r","CP-reg","AR(p)")
    
    l1<-colMeans(data.frame(autor[[1]][2]),na.rm=TRUE)
    l2<-colMeans(data.frame(autor[[2]][2]),na.rm=TRUE)
    l3<-colMeans(data.frame(autor[[3]][2]),na.rm=TRUE)
    l4<-colMeans(data.frame(autor[[4]][2]),na.rm=TRUE)
    l5<-colMeans(data.frame(autor[[5]][2]),na.rm=TRUE)
    l6<-colMeans(data.frame(autor[[6]][2]),na.rm=TRUE)
    
    L<-data.frame(rbind(l1,l2,l3,l4,l5,l6))
    L$AR<-c(0.005,0.1,0.3,0.5,0.7,0.9)
    names(L)<-c("DCP-QR","DCP-opt","CQR","CQR-m","CQR-r","CP-reg","AR(p)")
    
    #q<-dataframe(cbind(autor[[1]][9],autor[[2]][9],autor[[3]][9],
     #        autor[[4]][9],autor[[5]][9],autor[[6]][9]))
    s<-data.frame(rep(v,6))
    k<-list(C,L,s,j)
    return(k)
  }

autor<-list()  
ar1<-alu2(u,u2,50,1000,0,1)  

autor<-list()
ar2<-alu2(u,u2,70,1000,0,1)  

autor<-list()
ar3<-alu2(u,u2,100,500,0,1)  

autor<-list()
ar4<-alu2(u,u2,300,300,0,1)  

autor<-list()
ar5<-alu2(u,u2,500,300,0,1)  

write.xlsx(ar1,"ar1.xlsx")
write.xlsx(ar2,"ar2.xlsx")
write.xlsx(ar3,"ar3.xlsx")
write.xlsx(ar4,"ar4.xlsx")
write.xlsx(ar5,"ar5.xlsx")
#####################################################
u3<-0.02
autor<-list()  
ar12<-alu2(u,u3,50,1000,0,1)  

autor<-list()
ar22<-alu2(u,u3,70,1000,0,1)  

autor<-list()
ar32<-alu2(u,u3,100,500,0,1)  

autor<-list()
ar42<-alu2(u,u3,300,300,0,1)  

autor<-list()
ar52<-alu2(u,u3,500,300,0,1) 

#############################################################
u4<-0.09
autor<-list()  
ar13<-alu2(u,u4,50,1000,0,1)  

autor<-list()
ar23<-alu2(u,u4,70,1000,0,1)  

autor<-list()
ar33<-alu2(u,u4,100,500,0,1)  

autor<-list()
ar43<-alu2(u,u4,300,300,0,1)  

autor<-list()
ar53<-alu2(u,u4,500,300,0,1) 

#################################################################
#######################################################################
u2<--0.1
autor<-list()  
ar1b<-alu2(u,u2,50,1000,-1.5,1.5)  

autor<-list()
ar2b<-alu2(u,u2,70,1000,-1.5,1.5)  

autor<-list()
ar3b<-alu2(u,u2,100,500,-1.5,1.5)  

autor<-list()
ar4b<-alu2(u,u2,300,300,-1.5,1.5)  

autor<-list()
ar5b<-alu2(u,u2,500,300,-1.5,1.5)

#############################################################
#########################################################################
u3<-0.02

autor<-list()  
ar1c<-alu2(u,u3,50,1000,-1.5,1.5)  

autor<-list()
ar2c<-alu2(u,u3,70,1000,-1.5,1.5)  

autor<-list()
ar3c<-alu2(u,u3,100,500,-1.5,1.5)  

autor<-list()
ar4c<-alu2(u,u3,300,300,-1.5,1.5)  

autor<-list()
ar5c<-alu2(u,u3,500,300,-1.5,1.5)

#############################################################
#######################################################################
u4<-0.09
autor<-list()  
ar1d<-alu2(u,u4,50,1000,-1.5,1.5)  

autor<-list()
ar2d<-alu2(u,u4,70,1000,-1.5,1.5)  

autor<-list()
ar3d<-alu2(u,u4,100,500,-1.5,1.5)  

autor<-list()
ar4d<-alu2(u,u4,300,300,-1.5,1.5)  

autor<-list()
ar5d<-alu2(u,u4,500,300,-1.5,1.5)

##########################################################################
#################################################################################


####Preparig for Latex ########################################
#########################################################################

plotf<-function(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o){
  t1<-round(cbind.data.frame(a[1],a[2],rep(50,6),rep(as.vector(unlist(a[3]))[1],6)),3)
  t2<-round(cbind.data.frame(b[1],b[2],rep(70,6),rep(as.vector(unlist(b[3]))[1],6)),3)
  t3<-round(cbind.data.frame(c[1],c[2],rep(100,6),rep(as.vector(unlist(c[3]))[1],6)),3)
  t4<-round(cbind.data.frame(d[1],d[2],rep(300,6),rep(as.vector(unlist(d[3]))[1],6)),3)
  t5<-round(cbind.data.frame(e[1],e[2],rep(500,6),rep(as.vector(unlist(e[3]))[1],6)),3)
  
  t6<-round(cbind.data.frame(f[1],f[2],rep(50,6),rep(as.vector(unlist(f[3]))[1],6)),3)
  t7<-round(cbind.data.frame(g[1],g[2],rep(70,6),rep(as.vector(unlist(g[3]))[1],6)),3)
  t8<-round(cbind.data.frame(h[1],h[2],rep(100,6),rep(as.vector(unlist(h[3]))[1],6)),3)
  t9<-round(cbind.data.frame(i[1],i[2],rep(300,6),rep(as.vector(unlist(i[3]))[1],6)),3)
  t10<-round(cbind.data.frame(j[1],j[2],rep(500,6),rep(as.vector(unlist(j[3]))[1],6)),3)
  
  t11<-round(cbind.data.frame(k[1],k[2],rep(50,6),rep(as.vector(unlist(k[3]))[1],6)),3)
  t12<-round(cbind.data.frame(l[1],l[2],rep(70,6),rep(as.vector(unlist(l[3]))[1],6)),3)
  t13<-round(cbind.data.frame(m[1],m[2],rep(100,6),rep(as.vector(unlist(m[3]))[1],6)),3)
  t14<-round(cbind.data.frame(n[1],n[2],rep(300,6),rep(as.vector(unlist(n[3]))[1],6)),3)
  t15<-round(cbind.data.frame(o[1],o[2],rep(500,6),rep(as.vector(unlist(o[3]))[1],6)),3)
  
  names(t1)<-names(t2)<-names(t3)<-names(t4)<-names(t5)<-c("Cov1","Cov2","Cov3","Cov4","Cov5","Cov6",
                                                           "ar","L1","L2","L3","L4","L5","L6",
                                                           "s1","Sample","s3")
  
  names(t6)<-names(t7)<-names(t8)<-names(t9)<-names(t10)<-c("Cov1","Cov2","Cov3","Cov4","Cov5","Cov6",
                                                            "ar","L1","L2","L3","L4","L5","L6",
                                                            "s1","Sample","s3")
  
  names(t11)<-names(t12)<-names(t13)<-names(t14)<-names(t15)<-c("Cov1","Cov2","Cov3","Cov4","Cov5","Cov6",
                                                                "ar","L1","L2","L3","L4","L5","L6",
                                                                "s1","Sample","s3")
  
  T<-rbind.data.frame(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,
                      t11,t12,t13,t14,t15)
  T2<-T%>%mutate(group=paste(s1,s3,sep=","))
  T2$Sample<-as.character(T2$Sample)
  
  T2$Sample<-as.character(T2$Sample)
  T2<-T2[!(T2$Sample=="10"),]
  
  T2a<-T2[,-c(7:14,16)]
  names(T2a)<-c("DCP-QR","DCP-opt","CQR","CQR-m",
                "CQR-r","CP-reg","Sample","g")
  
  T3<-melt(T2a,id.vars=c("Sample","g"))
  T3$Sample<-as.numeric(T3$Sample)
  
  g<-ggplot(T3,aes(Sample,value,col=variable))+geom_line(linewidth=1.25)+
    facet_wrap(~g)+scale_y_continuous(breaks = seq(0,1, by = 0.05))+ 
    geom_hline(yintercept=0.9,col="black",linewidth=1)+
    xlab("Sample Size")+
    ylab("Average Coverage")+geom_point()+
    scale_color_manual(values=c("#cf597e","#eeb479","#9ccb86",
                        "#1b98e0","#ca7dcc","lightblue"))
  
  g1<-g+theme_base()
  #hc()+scale_colour_hc()
  
  T2b<-T2[,-c(1:7,14,16)]
  names(T2b)<-c("DCP-QR","DCP-opt","CQR","CQR-m",
                "CQR-r","CP-reg","Sample","g")
  T3<-melt(T2b,id.vars=c("Sample","g"))
  T3$Sample<-as.numeric(T3$Sample)
  
  g2<-ggplot(T3,aes(Sample,value,col=variable))+geom_line(linewidth=1)+
    facet_wrap(~g)+
    #scale_y_continuous(breaks=seq(0,75,by=5))+
    xlab("Sample Size")+
    ylab("Average Length")+geom_point()+
    scale_color_manual(values=c("#cf597e","#eeb479","#9ccb86",
                                         "#1b98e0","#ca7dcc","lightblue"))
                                         
  
 g12<-g2+theme_base()
  #hc()+scale_colour_hc()
  
  k<-list(g1,g12)
  return(k)
}

#####Uniform(0,1)
plotf(ar1,ar2,ar3,ar4,ar5,ar12,ar22,ar32,ar42,ar52,ar13,ar23,ar33,ar43,ar53)
write.xlsx(ar12,"ar12.xlsx")
write.xlsx(ar22,"ar22.xlsx")
write.xlsx(ar32,"ar32.xlsx")
write.xlsx(ar42,"ar42.xlsx")
write.xlsx(ar52,"ar52.xlsx")
write.xlsx(ar13,"ar13.xlsx")
write.xlsx(ar23,"ar23.xlsx")
write.xlsx(ar33,"ar33.xlsx")
write.xlsx(ar43,"ar43.xlsx")
write.xlsx(ar53,"ar53.xlsx")

#####Uniform(-1.5,1.5)  
plotf(ar1b,ar2b,ar3b,ar4b,ar5b,ar1c,ar2c,ar3c,ar4c,ar5c,ar1d,ar2d,ar3d,ar4d,ar5d)
write.xlsx(ar12,"ar1b.xlsx")
write.xlsx(ar22,"ar2b.xlsx")
write.xlsx(ar32,"ar3b.xlsx")
write.xlsx(ar42,"ar4b.xlsx")
write.xlsx(ar52,"ar5b.xlsx")
write.xlsx(ar13,"ar1c.xlsx")
write.xlsx(ar23,"ar2c.xlsx")
write.xlsx(ar33,"ar3c.xlsx")
write.xlsx(ar43,"ar4c.xlsx")
write.xlsx(ar53,"ar5c.xlsx")
write.xlsx(ar13,"ar1d.xlsx")
write.xlsx(ar23,"ar2d.xlsx")
write.xlsx(ar33,"ar3d.xlsx")
write.xlsx(ar43,"ar4d.xlsx")
write.xlsx(ar53,"ar5d.xlsx")

################################################################################

###Avg-Coverage-Unif(0,1)
a<-colMeans(data.frame(ar1[1]),na.rm=TRUE)
b<-colMeans(data.frame(ar2[1]),na.rm=TRUE)
c<-colMeans(data.frame(ar3[1]),na.rm=TRUE)
d<-colMeans(data.frame(ar4[1]),na.rm=TRUE)
e<-colMeans(data.frame(ar5[1]),na.rm=TRUE)
f<-colMeans(data.frame(ar12[1]),na.rm=TRUE)
g<-colMeans(data.frame(ar22[1]),na.rm=TRUE)
h<-colMeans(data.frame(ar32[1]),na.rm=TRUE)
i<-colMeans(data.frame(ar42[1]),na.rm=TRUE)
j<-colMeans(data.frame(ar52[1]),na.rm=TRUE)
k<-colMeans(data.frame(ar13[1]),na.rm=TRUE)
l<-colMeans(data.frame(ar23[1]),na.rm=TRUE)
m<-colMeans(data.frame(ar33[1]),na.rm=TRUE)
n<-colMeans(data.frame(ar43[1]),na.rm=TRUE)
o<-colMeans(data.frame(ar53[1]),na.rm=TRUE)

A<-rbind.data.frame(a,b,c,d,e)
A$Sample<-c("50","70","100","300","500")
names(A)<-c("DCP","DCPopt","CQR","CQRm","CQRr","CPreg","ARp","Sample")

plot(A$Sample,A$DCP,type="l",xlab="n",ylab="Avg.Coverage",col="red",lwd=2,ylim=c(0.85,1)) 
lines(A$Sample,A$DCPopt,col="blue",lwd=2)
lines(A$Sample,A$CQR,col="darkgreen",lwd=2)
lines(A$Sample,A$CQRm,col="orange",lwd=2)
lines(A$Sample,A$CQRr,col="purple",lwd=2)
lines(A$Sample,A$CPreg,lwd=2,col="brown")
abline(h=0.9,lty="dashed",lwd=2)
legend(390,1,legend=c("DCP","DCP-opt","CQR","CQR-m","CQR-r","CP-reg"),
       fill=c("red","blue","darkgreen","orange","purple","brown"))

####################################################################################
#######################Avg.Coverage:-1.5,1.5#########################################
#####################################################################################






  
  




