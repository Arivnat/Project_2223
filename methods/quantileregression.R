library(quantreg)
library(ggplot2)
library(latex2exp)

set.seed(12345)
x<-runif(n,1,10)
eps<-rnorm(n,0,4*(abs(x))^2)
y<--10+6.25*x+eps

data<-cbind.data.frame(x,y)

low<-rq(y~x,data=data,tau=0.05)
med<-rq(y~x,data=data,tau=0.5)
hi<-rq(y~x,data=data,tau=0.95)

ggplot(data, aes(x,y)) +
  geom_point(col="orange") + 
  geom_abline(intercept=low$coefficients[1], slope=low$coefficients[2],col="red",lwd=1)+
  geom_abline(intercept=med$coefficients[1], slope=med$coefficients[2])+
  geom_abline(intercept=hi$coefficients[1], slope=hi$coefficients[2],col="green",lwd=1)+
  geom_smooth(method="lm",se=FALSE)+
  xlab("X")+ylab("Y")+
  annotate(x=7.5,y=-400,label=TeX("$\\hat{Q}_{Y|X}(0.05)$"),geom="label")+
  annotate(x=7.5,y=400,label=TeX("$\\hat{Q}_{Y|X}(0.95)$"),geom="label")+
  annotate(x=9,y=200,label=TeX("$\\hat{Q}_{Y|X}(0.5)$"),geom="label")+
  annotate(x=7.5,y=-15,label="Line of best fit(given in blue)",geom="label")+
  geom_curve(
    aes(x = 9, y = 200, xend = 10.5, yend = 50 ),angle=90,curvature=0.5,
    arrow = arrow(length = unit(0.03, "npc"))
  )

######Pinball loss
x<--5:5
h<-function(x,alpha,alpha2,alpha3){
  f<-ifelse(x>=0,x*alpha,abs(x)*(1-alpha))
  g<-ifelse(x>=0,x*alpha2,abs(x)*(1-alpha2))
  h<-ifelse(x>=0,x*alpha3,abs(x)*(1-alpha3))
  j<-cbind.data.frame(x,f,g,h)
  z<-ggplot(j,mapping=aes(x,f))+geom_line(lwd=1)+
     geom_line(j,mapping=aes(x,g),col="red",lwd=1)+
     geom_line(j,mapping=aes(x,h),col="green",lwd=1)+
    xlab("X")+ylab(TeX(r'($\Q_{\tau}$(x)))'))+
    annotate(x=2.5,y=0.5*2.5,label=TeX("$\\tau$:0.5"),geom="label")+
    annotate(x=3,y=0.75*2.5,label=TeX("$\\tau$:0.75"),geom="label")+
    annotate(x=2,y=0.9*2.5,label=TeX("$\\tau$:0.9"),geom="label")
    
  return(z)
}
  
  
  
  
