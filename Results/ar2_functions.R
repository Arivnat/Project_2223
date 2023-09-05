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
#####################CAUCHY DISTRIBUTION FOR ERROR####################################
###loading required packages


###ar2 functions ##############################\
##################################################

aec2<-function(n,j,k,t,c,d){
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
    err<-15*arima.sim(n,model=list(ar=c(j,k)),innov=rcauchy(n,c,d))
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
  
  l<-list(cov.mat,leng.mat,X0.mat,Y0.mat,X1.mat,Y1.mat,X.test.mat,Y.test.mat)
  return(l)
}

###########################################
##########################################################################################

autor<-list()
u<-list(0.005,0.1,0.3,0.5,0.7,0.9)
u2<--0.1  ###-0.1 0.02 0.09


ca2<-function(u,v,n,t,c,d){
  start<-Sys.time()
  for (i in 1:length(u)){
    autor[[i]]<-aec2(n,u[[i]],u2,t,c,d)
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

#######################################################################
################################################################################
autor<-list()
j1<-ca2(u,u2,50,1000,-1,3)

autor<-list()
j2<-ca2(u,u2,70,1000,-1,3)

autor<-list()
j3<-ca2(u,u2,100,500,-1,3)

autor<-list()
j4<-ca2(u,u2,300,300,-1,3)

autor<-list()
j5<-ca2(u,u2,500,300,-1,3)


#########################################################
###################################################################################
u3<-0.02
autor<-list()
j12<-ca2(u,u3,50,1000,-1,3)

autor<-list()
j22<-ca2(u,u3,70,1000,-1,3)

autor<-list()
j32<-ca2(u,u3,100,500,-1,3)

autor<-list()
j42<-ca2(u,u3,300,300,-1,3)

autor<-list()
j52<-ca2(u,u3,500,300,-1,3)
###########################################################################
#################################################################################

u4<-0.09
autor<-list()
j13<-ca2(u,u4,50,1000,-1,3)

autor<-list()
j23<-ca2(u,u4,70,1000,-1,3)

autor<-list()
j33<-ca2(u,u4,100,500,-1,3)

autor<-list()
j43<-ca2(u,u4,300,300,-1,3)

autor<-list()
j53<-ca2(u,u4,500,300,-1,3)

##################################################################################
####################################################################################
#####################################################################################
library(writexl)
setwd("C:/Users/Alhmdulilla/Documents/ts_res")

###Coverage
write_xlsx(data.frame(j1[1]))
write_xlsx(data.frame(j2[1]))
write_xlsx(data.frame(j3[1]))
write_xlsx(data.frame(j4[1]))
write_xlsx(data.frame(j5[1]))
write_xlsx(data.frame(j12[1]))
write_xlsx(data.frame(j22[1]))
write_xlsx(data.frame(j32[1]))
write_xlsx(data.frame(j42[1]))
write_xlsx(data.frame(j52[1]))
write_xlsx(data.frame(j13[1]))
write_xlsx(data.frame(j23[1]))
write_xlsx(data.frame(j33[1]))
write_xlsx(data.frame(j43[1]))
write_xlsx(data.frame(j53[1]))

###Length
write_xlsx(data.frame(j1[2]))
write_xlsx(data.frame(j2[2]))
write_xlsx(data.frame(j3[2]))
write_xlsx(data.frame(j4[2]))
write_xlsx(data.frame(j5[2]))
write_xlsx(data.frame(j12[2]))
write_xlsx(data.frame(j22[2]))
write_xlsx(data.frame(j32[2]))
write_xlsx(data.frame(j42[2]))
write_xlsx(data.frame(j52[2]))
write_xlsx(data.frame(j13[2]))
write_xlsx(data.frame(j23[2]))
write_xlsx(data.frame(j33[2]))
write_xlsx(data.frame(j43[2]))
write_xlsx(data.frame(j53[2]))

###Computation time
write_xlsx(data.frame(j1[4]))
write_xlsx(data.frame(j2[4]))
write_xlsx(data.frame(j3[4]))
write_xlsx(data.frame(j4[4]))
write_xlsx(data.frame(j5[4]))
write_xlsx(data.frame(j12[4]))
write_xlsx(data.frame(j22[4]))
write_xlsx(data.frame(j32[4]))
write_xlsx(data.frame(j42[4]))
write_xlsx(data.frame(j52[4]))
write_xlsx(data.frame(j13[4]))
write_xlsx(data.frame(j23[4]))
write_xlsx(data.frame(j33[4]))
write_xlsx(data.frame(j43[4]))
write_xlsx(data.frame(j53[4]))


######################################################
######################################################
######################################################

####cauchy (0,1)
u2<--0.1

autor<-list()
k1<-ca2(u,u2,50,1000,0,1)

autor<-list()
k2<-ca2(u,u2,70,1000,0,1)

autor<-list()
k3<-ca2(u,u2,100,500,0,1)

autor<-list()
k4<-ca2(u,u2,300,300,0,1)

autor<-list()
k5<-ca2(u,u2,500,300,0,1)

################################################
################################################

u3<-0.02

autor<-list()
k12<-ca2(u,u3,50,1000,0,1)

autor<-list()
k22<-ca2(u,u3,70,1000,0,1)

autor<-list()
k32<-ca2(u,u3,100,500,0,1)

autor<-list()
k42<-ca2(u,u3,300,300,0,1)

autor<-list()
k52<-ca2(u,u3,500,300,0,1)

###########################################################
##########################################################

u4<-0.09
autor<-list()
k13<-ca2(u,u4,50,1000,0,1)

autor<-list()
k23<-ca2(u,u4,70,1000,0,1)

autor<-list()
k33<-ca2(u,u4,100,500,0,1)

autor<-list()
k43<-ca2(u,u4,300,300,0,1)

autor<-list()
k53<-ca2(u,u4,500,300,0,1)

##################################################################
########################################################################

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
                "CQR-r","CP-reg","Sample","group")
  
  T3<-melt(T2a,id.vars=c("Sample","group"))
  T3$Sample<-as.numeric(T3$Sample)
  
  g<-ggplot(T3,aes(Sample,value,col=variable))+geom_line(linewidth=1.25)+
    facet_wrap(~group)+scale_y_continuous(breaks = seq(0,1, by = 0.05))+ 
    geom_hline(yintercept=0.9,col="black",linewidth=1)+
    xlab("Sample Size")+
    ylab("Average Coverage")
  
  g1<-g+theme_base()
  
  T2b<-T2[,-c(1:7,14,16)]
  names(T2b)<-c("DCP-QR","DCP-opt","CQR","CQR-m",
                "CQR-r","CP-reg","Sample","group")
  T3<-melt(T2b,id.vars=c("Sample","group"))
  T3$Sample<-as.numeric(T3$Sample)
  
  g2<-ggplot(T3,aes(Sample,value,col=variable))+geom_line(linewidth=1)+
    facet_wrap(~group)+
    #scale_y_continuous(breaks=seq(0,75,by=5))+
    xlab("Sample Size")+
    ylab("Average Length")
  
  g12<-g2+theme_base()
  
  k<-list(g1,g12)
  return(k)
}

plotf(j1,j2,j3,j4,j5,j12,j22,j32,j42,j52,j13,j23,j33,j43,j53)


plotf(k1,k2,k3,k4,k5,k12,k22,k32,k42,k52,k13,k23,k33,k43,k53)


###############################################################
###################################################################
########################CAUCHY(2,5)#################################

####cauchy 
u2<--0.1

autor<-list()
s1<-ca2(u,u2,50,1000,2,5)

autor<-list()
s2<-ca2(u,u2,70,1000,2,5)

autor<-list()
s3<-ca2(u,u2,100,500,2,5)

autor<-list()
s4<-ca2(u,u2,300,300,2,5)

autor<-list()
s5<-ca2(u,u2,500,300,2,5)

################################################
################################################

u3<-0.02

autor<-list()
s12<-ca2(u,u3,50,1000,2,5)

autor<-list()
s22<-ca2(u,u3,70,1000,2,5)

autor<-list()
s32<-ca2(u,u3,100,500,2,5)

autor<-list()
s42<-ca2(u,u3,300,300,2,5)

autor<-list()
s52<-ca2(u,u3,500,300,2,5)

###########################################################
##########################################################

u4<-0.09

autor<-list()
s13<-ca2(u,u4,50,1000,2,5)

autor<-list()
s23<-ca2(u,u4,70,1000,2,5)

autor<-list()
s33<-ca2(u,u4,100,500,2,5)

autor<-list()
s43<-ca2(u,u4,300,300,2,5)

autor<-list()
s53<-ca2(u,u4,500,300,2,5)

plotf(s1,s2,s3,s4,s5,s12,s22,s32,s42,s52,s13,s23,s33,s43,s53)

####################################################################
#####################################################################
#################NORMAL#############################################
nor<-function(n,j,k,sd,t){
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
    err<-15*arima.sim(n,model=list(ar=c(j,k)),sd=sqrt(sd))
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
    
    taus    <- seq(0.05,0.95,length=200)  ##change taus to avoid infinite values 
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
  
  l<-list(cov.mat,leng.mat,X0.mat,Y0.mat,X1.mat,Y1.mat,X.test.mat,Y.test.mat)
  return(l)
}

autor<-list()
nor2<-function(u,v,n,sd,t){
  start<-Sys.time()
  for (i in 1:length(u)){
    autor[[i]]<-nor(n,u[[i]],u2,sd,t)
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

###############################################################################
##################################################################################
u2<--0.1

autor<-list()
t1<-nor2(u,u2,50,0.695,1000)

autor<-list()
t2<-nor2(u,u2,70,0.695,1000)

autor<-list()
t3<-nor2(u,u2,100,0.695,500)

autor<-list()
t4<-nor2(u,u2,300,0.695,300)

autor<-list()
t5<-nor2(u,u2,500,0.695,300)

#######################################################################################
##########################################################################################
u3<-0.02

autor<-list()
t12<-nor2(u,u3,50,0.695,1000)

autor<-list()
t22<-nor2(u,u3,70,0.695,1000)

autor<-list()
t32<-nor2(u,u3,100,0.695,500)

autor<-list()
t42<-nor2(u,u3,300,0.695,300)

autor<-list()
t52<-nor2(u,u3,500,0.695,300)

####################################################################################
#########################################################################################
u4<-0.09

autor<-list()
t13<-nor2(u,u4,50,0.695,1000)

autor<-list()
t23<-nor2(u,u4,70,0.695,1000)

autor<-list()
t33<-nor2(u,u4,100,0.695,500)

autor<-list()
t43<-nor2(u,u4,300,0.695,300)

autor<-list()
t53<-nor2(u,u4,500,0.695,300)

plotf(t1,t2,t3,t4,t5,t12,t22,t32,t42,t52,t13,t23,t33,t43,t53)

###################################################################################
#####################################################################################

##########################################################################
##################################################################################
u2<--0.1

autor<-list()
z1<-nor2(u,u2,50,1,1000)

autor<-list()
z2<-nor2(u,u2,70,1,1000)

autor<-list()
z3<-nor2(u,u2,100,1,500)

autor<-list()
z4<-nor2(u,u2,300,1,300)

autor<-list()
z5<-nor2(u,u2,500,1,300)

#######################################################################################
##########################################################################################
u3<-0.02

autor<-list()
z12<-nor2(u,u3,50,1,1000)

autor<-list()
z22<-nor2(u,u3,70,1,1000)

autor<-list()
z32<-nor2(u,u3,100,1,500)

autor<-list()
z42<-nor2(u,u3,300,1,300)

autor<-list()
z52<-nor2(u,u3,500,1,300)

####################################################################################
#########################################################################################
u4<-0.09

autor<-list()
z13<-nor2(u,u4,50,1,1000)

autor<-list()
z23<-nor2(u,u4,70,1,1000)

autor<-list()
z33<-nor2(u,u4,100,1,500)

autor<-list()
z43<-nor2(u,u4,300,1,300)

autor<-list()
z53<-nor2(u,u4,500,1,300)

plotf(z1,z2,z3,z4,z5,z12,z22,z32,z42,z52,z13,z23,z33,z43,z53)

