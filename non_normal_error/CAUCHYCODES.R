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

library(xtable)
library(hdm)
library(matrixcalc)
library(quantreg)
library(writexl)
library(oclust)


################simulation using cauchy distribution
other<-function(n){
  alpha.sig<-0.1
  T.ho <- floor(n*0.2) ##20% holdout sample
  
  #reps<-function(n,nsims,lpha.sig){
  cov.mat <- leng.mat <- X.test.mat <-NULL
  cov.mat2 <- leng.mat2<- X.test.mat2 <-NULL
  cov.mat3 <- leng.mat3 <- X.test.mat3 <- NULL
  cov.mat4 <- leng.mat4 <- X.test.mat4 <- NULL
  cov.mat5 <- leng.mat5 <- X.test.mat5 <- NULL
  cov.mat6 <- leng.mat6 <- X.test.mat6 <- NULL
  cov.mat7 <- leng.mat7 <- X.test.mat7 <- NULL
  cov.mat8 <- leng.mat8 <- X.test.mat8 <- NULL
  cov.mat9 <- leng.mat9 <- X.test.mat9 <- NULL
  cov.mat10 <- leng.mat10<- X.test.mat10 <- NULL
  cov.mat11 <- leng.mat11<- X.test.mat11 <- NULL
  cov.mat12<- leng.mat12<- X.test.mat12 <- NULL
  cov.mat13<- leng.mat13<- X.test.mat13<- NULL
  cov.mat14<- leng.mat14 <- X.test.mat14 <- NULL
  cov.mat15<- leng.mat15<- X.test.mat15 <- NULL
  
  for (i in 1:1000){
    set.seed(i)
    x1<-rnorm(n,2,3)
    x2<-rnorm(n,4,2)
    res<-rcauchy(n,1,3)
    eps<-rcauchy(n,0,1)
    eps2<-rcauchy(n,-1,3)
    x<-cbind.data.frame(x1,x2)
    x<-as.matrix(x)
    y<-2+3.5*x1+2.5*x2+res
    y1<-2+3.5*x1+2.5*x2+res*sqrt(abs(x1-x2))
    y2<-2+3.5*x1+2.5*x2+res*(exp(x1))  ##outliers 
    y3<-2+3.5*x1+2.5*x2+eps 
    y4<-2+3.5*x1+2.5*x2+eps*sqrt(abs(x1-x2))
    y5<-2+3.5*x1+2.5*x2+eps*(exp(0.5*(x1+x2))) 
    y6<-2+3.5*x1+2.5*x2+eps*(exp(0.05*(x1+x2))) 
    y7<-2+3.5*x1+2.5*x2+eps2 
    y8<-2+3.5*x1+2.5*x2+eps2*sqrt(abs(x1-x2))
    y9<-2+3.5*x1+2.5*x2+eps2*(exp(0.5*(x1+x2))) 
    y10<-2+3.5*x1+2.5*x2+eps2*(exp(0.05*(x1+x2)))
    eps3<-runif(n,0,1)
    y11<-2+3.5*x1+2.5*x2+eps3
    y12<-2+3.5*x1+2.5*x2+eps3*sqrt(abs(x1-x2))
    y13<-2+3.5*x1+2.5*x2+eps3*(exp(0.5*(x1+x2)))
    y14<-2+3.5*x1+2.5*x2+eps3*(exp(0.05*(x1+x2)))
    eps3<-runif(n,1,10)
    y14<-2+3.5*x1+2.5*x2+eps3
    y15<-2+3.5*x1+2.5*x2+eps3*sqrt(abs(x1-x2))
    #y16<-2+3.5*x1+2.5*x2+eps3*(exp(0.5*(x1+x2)))#
    #y17<-2+3.5*x1+2.5*x2+eps3*(exp(0.05*(x1+x2)))#
    
    
    sing <- sing2<-sing3<-sing4<-sing5<-sing6<-sing7<-sing8<-sing9<-sing10<-TRUE
    sing11 <- sing12<-sing13<-sing14<-sing15<-TRUE
    while (sing==TRUE){
      while (sing2==TRUE){
        while (sing3==TRUE){
          while (sing4==TRUE){
            while (sing5==TRUE){
              while (sing6==TRUE){
                while (sing7==TRUE){
                  while (sing8==TRUE){
                    while (sing9==TRUE){
                      while (sing10==TRUE){
                        while (sing11==TRUE){
                          while (sing12==TRUE){
                            while (sing13==TRUE){
                              while (sing14==TRUE){
                                while (sing15==TRUE){
                                  
                                  ind <- sample(length(y),length(y),replace=FALSE)
                                  Y   <- y[ind]
                                  Y2   <- y2[ind]
                                  Y3   <- y3[ind]
                                  Y4   <- y4[ind]
                                  Y5   <- y5[ind]
                                  Y6   <- y6[ind]
                                  Y7   <- y7[ind]
                                  Y8   <- y8[ind]
                                  Y9   <- y9[ind]
                                  Y10   <- y10[ind]
                                  Y11   <- y11[ind]
                                  Y12  <- y12[ind]
                                  Y13  <- y13[ind]
                                  Y14   <- y14[ind]
                                  Y15  <- y15[ind]
                                  X   <- x[ind,]
                                  X   <- as.matrix(X)
                                  
                                  ind.test  <- (length(Y)-T.ho+1):length(Y)
                                  Y.test    <- Y[ind.test]
                                  X.test    <- cbind(X[ind.test,])
                                  
                                  ind.test2  <- (length(Y2)-T.ho+1):length(Y2)
                                  Y.test2    <- Y2[ind.test2]
                                  X.test2    <- cbind(X[ind.test2,])
                                  
                                  ind.test3  <- (length(Y3)-T.ho+1):length(Y3)
                                  Y.test3    <- Y3[ind.test3]
                                  X.test3    <- cbind(X[ind.test3,])
                                  
                                  
                                  ind.test4  <- (length(Y4)-T.ho+1):length(Y4)
                                  Y.test4  <- Y4[ind.test2]
                                  X.test4    <- cbind(X[ind.test4,])
                                  
                                  ind.test5  <- (length(Y5)-T.ho+1):length(Y5)
                                  Y.test5    <- Y5[ind.test2]
                                  X.test5    <- cbind(X[ind.test5,])
                                  
                                  ind.test6  <- (length(Y6)-T.ho+1):length(Y6)
                                  Y.test6    <- Y6[ind.test6]
                                  X.test6    <- cbind(X[ind.test6,])
                                  
                                  
                                  ind.test7  <- (length(Y7)-T.ho+1):length(Y7)
                                  Y.test7    <- Y7[ind.test7]
                                  X.test7    <- cbind(X[ind.test7,])
                                  
                                  ind.test8  <- (length(Y8)-T.ho+1):length(Y8)
                                  Y.test8    <- Y8[ind.test8]
                                  X.test8    <- cbind(X[ind.test8,])
                                  
                                  ind.test9  <- (length(Y9)-T.ho+1):length(Y9)
                                  Y.test9   <- Y9[ind.test9]
                                  X.test9   <- cbind(X[ind.test9,])
                                  
                                  ind.test10  <- (length(Y10)-T.ho+1):length(Y10)
                                  Y.test10   <- Y10[ind.test10]
                                  X.test10    <- cbind(X[ind.test10,])
                                  
                                  ind.test11  <- (length(Y11)-T.ho+1):length(Y11)
                                  Y.test11    <- Y11[ind.test11]
                                  X.test11    <- cbind(X[ind.test11,])
                                  
                                  ind.test12  <- (length(Y12)-T.ho+1):length(Y12)
                                  Y.test12   <- Y12[ind.test12]
                                  X.test12    <- cbind(X[ind.test12,])
                                  
                                  ind.test13  <- (length(Y13)-T.ho+1):length(Y13)
                                  Y.test13    <- Y13[ind.test13]
                                  X.test13    <- cbind(X[ind.test13,])
                                  
                                  ind.test14  <- (length(Y14)-T.ho+1):length(Y14)
                                  Y.test14   <- Y14[ind.test14]
                                  X.test14   <- cbind(X[ind.test14,])
                                  
                                  ind.test15  <- (length(Y15)-T.ho+1):length(Y15)
                                  Y.test15    <- Y15[ind.test15]
                                  X.test15    <- cbind(X[ind.test15,])
                                  
                                  obj.uneven <- aux.uneven(Y[-ind.test],cbind(X[-ind.test,]))
                                  obj.uneven2 <- aux.uneven(Y2[-ind.test2],cbind(X[-ind.test2,]))
                                  obj.uneven3 <- aux.uneven(Y3[-ind.test3],cbind(X[-ind.test3,]))
                                  obj.uneven4 <- aux.uneven(Y4[-ind.test4],cbind(X[-ind.test4,]))
                                  obj.uneven5 <- aux.uneven(Y5[-ind.test5],cbind(X[-ind.test5,]))
                                  obj.uneven6 <- aux.uneven(Y6[-ind.test6],cbind(X[-ind.test6,]))
                                  obj.uneven7 <- aux.uneven(Y7[-ind.test7],cbind(X[-ind.test7,]))
                                  obj.uneven8 <- aux.uneven(Y8[-ind.test8],cbind(X[-ind.test8,]))
                                  obj.uneven9 <- aux.uneven(Y9[-ind.test9],cbind(X[-ind.test9,]))
                                  obj.uneven10 <- aux.uneven(Y10[-ind.test10],cbind(X[-ind.test10,]))
                                  obj.uneven11 <- aux.uneven(Y11[-ind.test11],cbind(X[-ind.test11,]))
                                  obj.uneven12 <- aux.uneven(Y12[-ind.test12],cbind(X[-ind.test12,]))
                                  obj.uneven13 <- aux.uneven(Y13[-ind.test13],cbind(X[-ind.test13,]))
                                  obj.uneven14 <- aux.uneven(Y14[-ind.test14],cbind(X[-ind.test14,]))
                                  obj.uneven15 <- aux.uneven(Y15[-ind.test15],cbind(X[-ind.test15,]))
                                  
                                  X0 <- cbind(obj.uneven$X0)
                                  Ya<- obj.uneven$Y0
                                  X1 <- cbind(obj.uneven$X1)
                                  Y1 <- obj.uneven$Y1
                                  
                                  X02 <- cbind(obj.uneven2$X0)
                                  Y02 <- obj.uneven2$Y0
                                  X12 <- cbind(obj.uneven2$X1)
                                  Y12 <- obj.uneven2$Y1
                                  
                                  X03 <- cbind(obj.uneven3$X0)
                                  Y03 <- obj.uneven3$Y0
                                  X13<- cbind(obj.uneven3$X1)
                                  Y13 <- obj.uneven3$Y1
                                  
                                  X04 <- cbind(obj.uneven4$X0)
                                  Y04 <- obj.uneven4$Y0
                                  X14<- cbind(obj.uneven4$X1)
                                  Y14 <- obj.uneven4$Y1
                                  
                                  X05 <- cbind(obj.uneven5$X0)
                                  Y05 <- obj.uneven5$Y0
                                  X15<- cbind(obj.uneven5$X1)
                                  Y15 <- obj.uneven5$Y1
                                  
                                  X06 <- cbind(obj.uneven6$X0)
                                  Y06 <- obj.uneven6$Y0
                                  X16<- cbind(obj.uneven6$X1)
                                  Y16 <- obj.uneven6$Y1
                                  
                                  X07 <- cbind(obj.uneven7$X0)
                                  Y07 <- obj.uneven7$Y0
                                  X17<- cbind(obj.uneven7$X1)
                                  Y17 <- obj.uneven7$Y1
                                  
                                  X08 <- cbind(obj.uneven8$X0)
                                  Y08 <- obj.uneven8$Y0
                                  X18<- cbind(obj.uneven8$X1)
                                  Y18 <- obj.uneven8$Y1
                                  
                                  X09 <- cbind(obj.uneven9$X0)
                                  Y09 <- obj.uneven9$Y0
                                  X19<- cbind(obj.uneven9$X1)
                                  Y19 <- obj.uneven9$Y1
                                  
                                  X10 <- cbind(obj.uneven10$X0)
                                  Y10 <- obj.uneven10$Y0
                                  X110<- cbind(obj.uneven10$X1)
                                  Y110 <- obj.uneven10$Y1
                                  
                                  X11 <- cbind(obj.uneven11$X0)
                                  Y11 <- obj.uneven11$Y0
                                  X111<- cbind(obj.uneven11$X1)
                                  Y111 <- obj.uneven11$Y1
                                  
                                  X12 <- cbind(obj.uneven12$X0)
                                  Y12 <- obj.uneven12$Y0
                                  X112<- cbind(obj.uneven12$X1)
                                  Y112 <- obj.uneven12$Y1
                                  
                                  X13 <- cbind(obj.uneven13$X0)
                                  Y13 <- obj.uneven13$Y0
                                  X113<- cbind(obj.uneven13$X1)
                                  Y113 <- obj.uneven13$Y1
                                  
                                  X14 <- cbind(obj.uneven14$X0)
                                  Y14 <- obj.uneven14$Y0
                                  X114<- cbind(obj.uneven14$X1)
                                  Y114 <- obj.uneven14$Y1
                                  
                                  X15 <- cbind(obj.uneven15$X0)
                                  Y15 <- obj.uneven15$Y0
                                  X115<- cbind(obj.uneven15$X1)
                                  Y115 <- obj.uneven15$Y1
                                  
                                  sing  <- is.singular.matrix(t(X0)%*%X0)
                                  sing2  <- is.singular.matrix(t(X02)%*%X02)
                                  sing3  <- is.singular.matrix(t(X03)%*%X03)
                                  sing4  <- is.singular.matrix(t(X04)%*%X04)
                                  sing5  <- is.singular.matrix(t(X05)%*%X05)
                                  sing6  <- is.singular.matrix(t(X06)%*%X06)
                                  sing7  <- is.singular.matrix(t(X07)%*%X07)
                                  sing8  <- is.singular.matrix(t(X08)%*%X08)
                                  sing9  <- is.singular.matrix(t(X09)%*%X09)
                                  sing10  <- is.singular.matrix(t(X10)%*%X10)
                                  sing11 <- is.singular.matrix(t(X11)%*%X11)
                                  sing12 <- is.singular.matrix(t(X12)%*%X12)
                                  sing13 <- is.singular.matrix(t(X13)%*%X13)
                                  sing14 <- is.singular.matrix(t(X14)%*%X14)
                                  sing15  <- is.singular.matrix(t(X15)%*%X15)
                                  
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    
    
    taus    <- seq(0.001,0.999,length=200)
    # ys      <- quantile(unique(c(Y0,Y1)),seq(0.001,0.999,length=200))
    #  ys2      <- quantile(unique(c(Y02,Y12)),seq(0.001,0.999,length=200))
    
    # Applying the difference conformal prediction methods
    res.qr    <- dcp.qr(Ya,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)
    res.qro<-dcp.opt(Ya,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)
    res.cqr   <- cqr(Ya,X0,Y1,X1,Y.test,X.test,alpha.sig)
    res.reg   <- cp.reg(Ya,X0,Y1,X1,Y.test,X.test,alpha.sig)
    
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
    res.reg4   <- cp.reg(Y04,X04,Y14,X14,Y.test2,X.test4,alpha.sig)
    
    res.qr5    <- dcp.qr(Y05,X05,Y15,X15,Y.test5,X.test5,taus,alpha.sig)
    res.qro5<-dcp.opt(Y05,X05,Y15,X15,Y.test5,X.test5,taus,alpha.sig)
    res.cqr5   <- cqr(Y05,X05,Y15,X15,Y.test5,X.test5,alpha.sig)
    res.reg5   <- cp.reg(Y05,X05,Y15,X15,Y.test5,X.test5,alpha.sig)
    
    res.qr6    <- dcp.qr(Y06,X06,Y16,X16,Y.test6,X.test6,taus,alpha.sig)
    res.qro6<-dcp.opt(Y06,X06,Y16,X16,Y.test6,X.test6,taus,alpha.sig)
    res.cqr6   <- cqr(Y06,X06,Y16,X16,Y.test6,X.test6,alpha.sig)
    res.reg6   <- cp.reg(Y06,X06,Y16,X16,Y.test6,X.test6,alpha.sig)
    
    res.qr7    <- dcp.qr(Y07,X07,Y17,X17,Y.test7,X.test7,taus,alpha.sig)
    res.qro7<-dcp.opt(Y07,X07,Y17,X17,Y.test7,X.test7,taus,alpha.sig)
    res.cqr7   <- cqr(Y07,X07,Y17,X17,Y.test7,X.test7,alpha.sig)
    res.reg7   <- cp.reg(Y07,X07,Y17,X17,Y.test7,X.test7,alpha.sig)
    
    res.qr8    <- dcp.qr(Y08,X08,Y18,X18,Y.test8,X.test8,taus,alpha.sig)
    res.qro8<-dcp.opt(Y08,X08,Y18,X18,Y.test8,X.test8,taus,alpha.sig)
    res.cqr8   <- cqr(Y08,X08,Y18,X18,Y.test8,X.test8,alpha.sig)
    res.reg8   <- cp.reg(Y08,X08,Y18,X18,Y.test8,X.test8,alpha.sig)
    
    res.qr9    <- dcp.qr(Y09,X09,Y19,X19,Y.test9,X.test9,taus,alpha.sig)
    res.qro9<-dcp.opt(Y09,X09,Y19,X19,Y.test9,X.test9,taus,alpha.sig)
    res.cqr9   <- cqr(Y09,X09,Y19,X19,Y.test9,X.test9,alpha.sig)
    res.reg9   <- cp.reg(Y09,X09,Y19,X19,Y.test9,X.test9,alpha.sig)
    
    res.qr10    <- dcp.qr(Y10,X10,Y110,X110,Y.test10,X.test10,taus,alpha.sig)
    res.qro10<-dcp.opt(Y10,X10,Y110,X110,Y.test10,X.test10,taus,alpha.sig)
    res.cqr10  <- cqr(Y10,X10,Y110,X110,Y.test10,X.test10,alpha.sig)
    res.reg10  <- cp.reg(Y10,X10,Y110,X110,Y.test10,X.test10,alpha.sig)
    
    res.qr11    <- dcp.qr(Y11,X11,Y111,X111,Y.test11,X.test11,taus,alpha.sig)
    res.qro11<-dcp.opt(Y11,X11,Y111,X111,Y.test11,X.test11,taus,alpha.sig)
    res.cqr11  <- cqr(Y11,X11,Y111,X111,Y.test11,X.test11,alpha.sig)
    res.reg11  <- cp.reg(Y11,X11,Y111,X111,Y.test11,X.test11,alpha.sig)
    
    res.qr12    <- dcp.qr(Y12,X12,Y112,X112,Y.test12,X.test12,taus,alpha.sig)
    res.qro12<-dcp.opt(Y12,X12,Y112,X112,Y.test12,X.test12,taus,alpha.sig)
    res.cqr12  <- cqr(Y12,X12,Y112,X112,Y.test12,X.test12,alpha.sig)
    res.reg12  <- cp.reg(Y12,X12,Y112,X112,Y.test12,X.test12,alpha.sig)
    
    res.qr13    <- dcp.qr(Y13,X13,Y113,X113,Y.test13,X.test13,taus,alpha.sig)
    res.qro13<-dcp.opt(Y13,X13,Y113,X113,Y.test13,X.test13,taus,alpha.sig)
    res.cqr13  <- cqr(Y13,X13,Y113,X113,Y.test13,X.test13,alpha.sig)
    res.reg13  <- cp.reg(Y13,X13,Y113,X113,Y.test13,X.test13,alpha.sig)
    
    res.qr14    <- dcp.qr(Y14,X14,Y114,X114,Y.test14,X.test14,taus,alpha.sig)
    res.qro14<-dcp.opt(Y14,X14,Y114,X114,Y.test14,X.test14,taus,alpha.sig)
    res.cqr14  <- cqr(Y14,X14,Y114,X114,Y.test14,X.test14,alpha.sig)
    res.reg14  <- cp.reg(Y14,X14,Y114,X114,Y.test14,X.test14,alpha.sig)
    
    res.qr15    <- dcp.qr(Y15,X15,Y115,X115,Y.test15,X.test15,taus,alpha.sig)
    res.qro15<-dcp.opt(Y15,X15,Y115,X115,Y.test15,X.test15,taus,alpha.sig)
    res.cqr15  <- cqr(Y15,X15,Y115,X115,Y.test15,X.test15,alpha.sig)
    res.reg15  <- cp.reg(Y15,X15,Y115,X115,Y.test15,X.test15,alpha.sig)
    
    X.test.mat <- rbind(X.test.mat,X.test)
    X.test.mat2 <- rbind(X.test.mat2,X.test2)
    X.test.mat3 <- rbind(X.test.mat3,X.test3)
    X.test.mat4 <- rbind(X.test.mat4,X.test4)
    X.test.mat5 <- rbind(X.test.mat5,X.test5)
    X.test.mat6 <- rbind(X.test.mat6,X.test6)
    X.test.mat7 <- rbind(X.test.mat7,X.test7)
    X.test.mat8 <- rbind(X.test.mat8,X.test8)
    X.test.mat9<- rbind(X.test.mat9,X.test9)
    X.test.mat10 <- rbind(X.test.mat10,X.test10)
    X.test.mat11<- rbind(X.test.mat11,X.test11)
    X.test.mat12<- rbind(X.test.mat12,X.test12)
    X.test.mat13<- rbind(X.test.mat13,X.test13)
    X.test.mat14<- rbind(X.test.mat14,X.test14)
    X.test.mat15<- rbind(X.test.mat15,X.test15)
    
    cov.mat.temp  <- cbind(res.qr$cov.qr,res.qro$cov.opt,res.cqr$cov.o,res.cqr$cov.m,res.cqr$cov.r,res.reg$cov.reg)
    leng.mat.temp <- cbind(res.qr$leng.qr,res.qro$leng.opt,res.cqr$leng.o,res.cqr$leng.m,res.cqr$leng.r,res.reg$leng.reg)
    
    cov.mat.temp2  <- cbind(res.qr2$cov.qr,res.qro2$cov.opt,res.cqr2$cov.o,res.cqr2$cov.m,res.cqr2$cov.r,res.reg2$cov.reg)
    leng.mat.temp2 <- cbind(res.qr2$leng.qr,res.qro2$leng.opt,res.cqr2$leng.o,res.cqr2$leng.m,res.cqr2$leng.r,res.reg2$leng.reg)
    
    cov.mat.temp3  <- cbind(res.qr3$cov.qr,res.qro3$cov.opt,res.cqr3$cov.o,res.cqr3$cov.m,res.cqr3$cov.r,res.reg3$cov.reg)
    leng.mat.temp3 <- cbind(res.qr3$leng.qr,res.qro3$leng.opt,res.cqr3$leng.o,res.cqr3$leng.m,res.cqr3$leng.r,res.reg3$leng.reg)
    
    cov.mat.temp4 <- cbind(res.qr4$cov.qr,res.qro4$cov.opt,res.cqr4$cov.o,res.cqr4$cov.m,res.cqr4$cov.r,res.reg4$cov.reg)
    leng.mat.temp4 <- cbind(res.qr4$leng.qr,res.qro4$leng.opt,res.cqr4$leng.o,res.cqr4$leng.m,res.cqr4$leng.r,res.reg4$leng.reg)
    
    cov.mat.temp5 <- cbind(res.qr5$cov.qr,res.qro5$cov.opt,res.cqr5$cov.o,res.cqr5$cov.m,res.cqr5$cov.r,res.reg5$cov.reg)
    leng.mat.temp5 <- cbind(res.qr5$leng.qr,res.qro5$leng.opt,res.cqr5$leng.o,res.cqr5$leng.m,res.cqr5$leng.r,res.reg5$leng.reg)
    
    cov.mat.temp6 <- cbind(res.qr6$cov.qr,res.qro6$cov.opt,res.cqr6$cov.o,res.cqr6$cov.m,res.cqr6$cov.r,res.reg6$cov.reg)
    leng.mat.temp6 <- cbind(res.qr6$leng.qr,res.qro6$leng.opt,res.cqr6$leng.o,res.cqr6$leng.m,res.cqr6$leng.r,res.reg6$leng.reg)
    
    cov.mat.temp7 <- cbind(res.qr7$cov.qr,res.qro7$cov.opt,res.cqr7$cov.o,res.cqr7$cov.m,res.cqr7$cov.r,res.reg7$cov.reg)
    leng.mat.temp7 <- cbind(res.qr7$leng.qr,res.qro7$leng.opt,res.cqr7$leng.o,res.cqr7$leng.m,res.cqr7$leng.r,res.reg7$leng.reg)
    
    cov.mat.temp8 <- cbind(res.qr8$cov.qr,res.qro8$cov.opt,res.cqr8$cov.o,res.cqr8$cov.m,res.cqr8$cov.r,res.reg8$cov.reg)
    leng.mat.temp8 <- cbind(res.qr8$leng.qr,res.qro8$leng.opt,res.cqr8$leng.o,res.cqr8$leng.m,res.cqr8$leng.r,res.reg8$leng.reg)
    
    
    cov.mat.temp9 <- cbind(res.qr9$cov.qr,res.qro9$cov.opt,res.cqr9$cov.o,res.cqr9$cov.m,res.cqr9$cov.r,res.reg9$cov.reg)
    leng.mat.temp9 <- cbind(res.qr9$leng.qr,res.qro9$leng.opt,res.cqr9$leng.o,res.cqr9$leng.m,res.cqr9$leng.r,res.reg9$leng.reg)
    
    cov.mat.temp10 <- cbind(res.qr10$cov.qr,res.qro10$cov.opt,res.cqr10$cov.o,res.cqr10$cov.m,res.cqr10$cov.r,res.reg10$cov.reg)
    leng.mat.temp10 <- cbind(res.qr10$leng.qr,res.qro10$leng.opt,res.cqr10$leng.o,res.cqr10$leng.m,res.cqr10$leng.r,res.reg10$leng.reg)
    
    cov.mat.temp11<- cbind(res.qr11$cov.qr,res.qro11$cov.opt,res.cqr11$cov.o,res.cqr11$cov.m,res.cqr11$cov.r,res.reg11$cov.reg)
    leng.mat.temp11 <- cbind(res.qr11$leng.qr,res.qro11$leng.opt,res.cqr11$leng.o,res.cqr11$leng.m,res.cqr11$leng.r,res.reg11$leng.reg)
    
    cov.mat.temp12 <- cbind(res.qr12$cov.qr,res.qro12$cov.opt,res.cqr12$cov.o,res.cqr12$cov.m,res.cqr12$cov.r,res.reg12$cov.reg)
    leng.mat.temp12 <- cbind(res.qr12$leng.qr,res.qro12$leng.opt,res.cqr12$leng.o,res.cqr12$leng.m,res.cqr12$leng.r,res.reg12$leng.reg)
    
    cov.mat.temp13 <- cbind(res.qr13$cov.qr,res.qro13$cov.opt,res.cqr13$cov.o,res.cqr13$cov.m,res.cqr13$cov.r,res.reg13$cov.reg)
    leng.mat.temp13 <- cbind(res.qr13$leng.qr,res.qro13$leng.opt,res.cqr13$leng.o,res.cqr13$leng.m,res.cqr13$leng.r,res.reg13$leng.reg)
    
    cov.mat.temp14 <- cbind(res.qr14$cov.qr,res.qro14$cov.opt,res.cqr14$cov.o,res.cqr14$cov.m,res.cqr14$cov.r,res.reg14$cov.reg)
    leng.mat.temp14 <- cbind(res.qr14$leng.qr,res.qro14$leng.opt,res.cqr14$leng.o,res.cqr14$leng.m,res.cqr14$leng.r,res.reg14$leng.reg)
    
    cov.mat.temp15 <- cbind(res.qr15$cov.qr,res.qro15$cov.opt,res.cqr15$cov.o,res.cqr15$cov.m,res.cqr15$cov.r,res.reg15$cov.reg)
    leng.mat.temp15 <- cbind(res.qr15$leng.qr,res.qro15$leng.opt,res.cqr15$leng.o,res.cqr15$leng.m,res.cqr15$leng.r,res.reg15$leng.reg)
    
    cov.mat   <- rbind(cov.mat,cov.mat.temp)
    leng.mat  <- rbind(leng.mat,leng.mat.temp)
    
    cov.mat2   <- rbind(cov.mat2,cov.mat.temp2)
    leng.mat2  <- rbind(leng.mat2,leng.mat.temp2)
    
    cov.mat3   <- rbind(cov.mat3,cov.mat.temp3)
    leng.mat3 <- rbind(leng.mat3,leng.mat.temp3)
    
    cov.mat4   <- rbind(cov.mat4,cov.mat.temp4)
    leng.mat4  <- rbind(leng.mat4,leng.mat.temp4)
    
    cov.mat5   <- rbind(cov.mat5,cov.mat.temp5)
    leng.mat5  <- rbind(leng.mat5,leng.mat.temp5)
    
    cov.mat6   <- rbind(cov.mat6,cov.mat.temp6)
    leng.mat6  <- rbind(leng.mat6,leng.mat.temp6)
    
    cov.mat7   <- rbind(cov.mat7,cov.mat.temp7)
    leng.mat7  <- rbind(leng.mat7,leng.mat.temp7)
    
    cov.mat8   <- rbind(cov.mat8,cov.mat.temp8)
    leng.mat8  <- rbind(leng.mat8,leng.mat.temp8)
    
    cov.mat9   <- rbind(cov.mat9,cov.mat.temp9)
    leng.mat9  <- rbind(leng.mat9,leng.mat.temp9)
    
    cov.mat10   <- rbind(cov.mat10,cov.mat.temp10)
    leng.mat10  <- rbind(leng.mat10,leng.mat.temp10)
    
    cov.mat11   <- rbind(cov.mat11,cov.mat.temp11)
    leng.mat11  <- rbind(leng.mat11,leng.mat.temp11)
    
    cov.mat12 <- rbind(cov.mat12,cov.mat.temp12)
    leng.mat12  <- rbind(leng.mat12,leng.mat.temp12)
    
    cov.mat13 <- rbind(cov.mat13,cov.mat.temp13)
    leng.mat13  <- rbind(leng.mat13,leng.mat.temp13)
    
    cov.mat14 <- rbind(cov.mat14,cov.mat.temp14)
    leng.mat14  <- rbind(leng.mat14,leng.mat.temp14)
    
    cov.mat15 <- rbind(cov.mat15,cov.mat.temp15)
    leng.mat15  <- rbind(leng.mat15,leng.mat.temp15)
  }
  
  
  
  
  
  
  l<-list(leng.mat,leng.mat2,leng.mat3,leng.mat4,leng.mat5,leng.mat6,leng.mat7,leng.mat8,leng.mat9,leng.mat10,leng.mat11,
          leng.mat12,leng.mat13,leng.mat14,leng.mat15,
          cov.mat,cov.mat2,cov.mat3,cov.mat4,cov.mat5,cov.mat6,cov.mat7,cov.mat8,cov.mat9,cov.mat10,cov.mat11,cov.mat12,
          cov.mat13,cov.mat14,cov.mat15,X.test.mat,X.test.mat2,X.test.mat3,X.test.mat4,X.test.mat5,X.test.mat6,
          X.test.mat7,X.test.mat8, X.test.mat9,
          X.test.mat10,X.test.mat11,X.test.mat12,X.test.mat13,X.test.mat14,X.test.mat15)
  return(l)
}


##sample size-50
start<-Sys.time()
o00<-other(50)
Sys.time()-start

##sample size-70
start<-Sys.time()
o0<-other(70)
Sys.time()-start

##sample size-100
start<-Sys.time()
o1<-other(100)
Sys.time()-start

##sample size-200
start<-Sys.time()
o2<-other(200)
Sys.time()-start

#start<-Sys.time()
#o3<-other(300)
#Sys.time()-start

###sample size-500
start<-Sys.time()
o4<-other(500)
Sys.time()-start

#######################################################

####Function to compute average coverage
avgc<-function(md){   
  avgcov<-NULL
  for (i in 16:30){
    md2<-data.frame(md[i])
    ac<-colMeans(md2,na.rm=TRUE)
    avgcov<-rbind(avgcov,ac)
  }
  row1<-c("(1)","(2)","(3)","(4)","(5)","(6)","(7)","(8)",
          "(9)","(10)","(11)","(12)","(13)","(14)","(15)")
  avgcov<-cbind.data.frame(row1,avgcov)
  names(avgcov)<-c("Type","DCP-QR","DCP-opt","CQR","CQR-m","CQR-r","CP-reg")
  #a<-xtable(avgcov)
  return(avgcov)
}

datc<-rbind.data.frame(q,r,s,t,v)
rownames(datc)<-NULL
##latex
xtable(datc,digits=3)

q<-avgc(o00)
r<-avgc(o0)
s<-avgc(o1)
t<-avgc(o2)
v<-avgc(o4)

datn<-rbind.data.frame(q,r,s,t,v)

##latex
xtable(datn,digits=3)


###function for computing average length
avgl<-function(md){   
  avglen<-NULL
  for (i in 1:15){
    md2<-data.frame(md[i])
    ac<-colMeans(md2,na.rm=TRUE)
    avglen<-rbind(avglen,ac)
  }
  row1<-c("(1)","(2)","(3)","(4)","(5)","(6)","(7)","(8)",
          "(9)","(10)","(11)","(12)","(13)","(14)","(15)")
  avglen<-cbind.data.frame(row1,avglen)
  names(avglen)<-c("Type","DCP-QR","DCP-opt","CQR","CQR-m","CQR-r","CP-reg")
  #a<-xtable(avglen)
  return(avglen)
}
q0<-avgl(o00)
r0<-avgl(o0)
s0<-avgl(o1)
t0<-avgl(o2)
v0<-avgl(o4)


################################################################
#########################################################################

datl<-rbind.data.frame(q0,r0,s0,t0,v0)
rownames(datl)<-NULL
xtable(datl,digits=3)

###function for computing standard deviation

sd1<-function(md){
  sdm<-NULL
  for (i in 16:30){
    cov<-data.frame(md[i])
    xmat<-as.matrix(data.frame(md[i+15]))
    pred.cov.qr       <- predict(glm(cov[,1]~xmat,family=binomial(link="logit")),type="response")
    pred.cov.opt      <- predict(glm(cov[,2]~xmat,family=binomial(link="logit")),type="response")
    pred.cov.cqr.o    <- predict(glm(cov[,3]~xmat,family=binomial(link="logit")),type="response")
    pred.cov.cqr.m    <- predict(glm(cov[,4]~xmat,family=binomial(link="logit")),type="response")
    pred.cov.cqr.r    <- predict(glm(cov[,5]~xmat,family=binomial(link="logit")),type="response")
    pred.cov.reg      <- predict(glm(cov[,6]~xmat,family=binomial(link="logit")),type="response")
    #pred.cov.loc      <- predict(glm(cov.mat[,8]~X.test.mat,family=binomial(link="logit")),type="response")
    pred.cov<-cbind.data.frame(pred.cov.qr,pred.cov.opt,pred.cov.cqr.o,pred.cov.cqr.m,pred.cov.cqr.r,pred.cov.reg)
    res.cov.cond <- apply(pred.cov,2,sd)
    sdm<-rbind(sdm,res.cov.cond)
  }
  sdm<-data.frame(sdm)
  names(sdm)<-c("DCP-QR","DCP-opt","CQR","CQR-m","CQR-r","CP-reg")
  row1<-c("(1)","(2)","(3)","(4)","(5)","(6)","(7)","(8)",
          "(9)","(10)","(11)","(12)","(13)","(14)","(15)")
  sdm1<-cbind.data.frame(row1,sdm)
  names(sdm1)<-c("Type","DCP-QR","DCP-opt","CQR","CQR-m","CQR-r","CP-reg")
  return(sdm1)
}

a<-sd1(o00)
b<-sd1(o1)
c<-sd1(o2)
d<-sd1(o3)
f<-sd1(o4)

sdl<-rbind.data.frame(a,b,c,d,f)
rownames(sdl)<-NULL
formatn<-function(x){
  y<-ifelse(abs(x)>=0.01,round(x,4),formatC(x,format="e",digits=2))
  return(y)
}

sdl[]<-lapply(sdl,formatn)
##latex
xtable(sdl,digits=5)

