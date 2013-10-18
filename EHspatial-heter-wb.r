###############################
library("mgcv")
library("spatstat")
#include required library
########
setwd('F:/DataW/lg-data/composition')
dir()
xm0a <- read.csv('lgdat2013.csv')
xm0 <- subset(xm0a,sp=='蚬木' & is.na(bra)==T )
if(any(xm0$dbh < 1)){xm0$dbh[xm0$dbh<1]=1}
w00　<- duplicated(paste(xm0$x,xm0$y))
xm0$x[w00] <- xm0$x[w00] + runif(sum(w00),0,0.1)
head(xm0)
###
tiff('DBH rank--abundance.tiff',
     width = 3000, height = 2600,res=600,compression = "lzw")
par( mar=c(5.5,6,0.5,0.5 ),pty="m",mex=0.6) 
 barplot(tab <- table(cut(xm0$dbh,c(1,3,9,16,28, Inf))), ylim=c(0,900),
         ylab="个体数 abundance", xlab="胸径等级 DBH rank (cm)",cex.lab=1.2)
 text(c(1-0.3, 2, 3+0.1, 4+0.3, 5+0.5,6+0.6),tab+40, labels=tab,cex=1.2)
dev.off()
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####
alpha=0.05
envelope=FALSE
#set the windows
win=owin(c(0,plotdim[1]),c(0,plotdim[2]))
#species abundance
n=table(xm0a$sp)
#richness
S=length(n)
# species name
sp=unique(xm0a$sp)
#total area
A0=plotdim[1]*plotdim[2]
#sample area vector
A=exp(seq(1,log(A0),length.out=nx))
#define quadrat shape
qsp=plotdim[2]/plotdim[1]
y=sqrt(qsp*A)
x=y/qsp
#
covr0=Topogra.tran()
names(covr0)
covr <- covr0[c(1,2,3,5,8)] 

#########################################################################
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
xm1 <- subset(xm0, dbh>=28)
#
xm2 <- ppp(x=xm1$x,y=xm1$y, window=win)
xm3.0 <-envelope(xm2,pcf,nsim=11,nrank=5)
xm3.1 <- as.data.frame(xm3.0)$obs
#
data.sp2<- data.frame()
xm.obs2 <- data.frame()
for(i in 1:49){
  sp1 = Inhom.Po(xm2,covr=covr,win=win,sig.t=sig.t)
data.sp1 <- as.data.frame(sp1); data.sp1$marks <- i
data.sp2 <-  rbind(data.sp2, data.sp1)
xm2.1 <- ppp(x=sp1$x,y=sp1$y, window=win)
xm3 <-envelope(xm2.1,pcf,nsim=11,nrank=5)
xm.obs1 <- as.data.frame(xm3); xm.obs1$theo <- i
xm.obs2 <- rbind(xm.obs2, xm.obs1)
}
 xm.wide <- xm.wide01 <- reshape(xm.obs2, v.names = "obs", idvar = "r",
                timevar = "theo", direction = "wide")
######### 
tiff('Inhomogeneous Poisson EH 28+cm.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 

plot(xm3.0,legend=FALSE,main="",#xlab="",ylab="", 
     xlim=c(0,80), ylim=c(0,25), main="",bty="l")
matpoints(x=xm.wide$r,y=xm.wide[,-(1:3)], type = "l",lty=1, col=3)
points(x=xm.wide$r,y=xm3.1,col=2,cex=0.5)

legend("top", "Inhomogeneous Poisson Process",lwd=3,
       cex=1.3,text.col=3,col=3,
       inset =0.01,bty="n")
text(50,15,"蚬木 > 28cm")
dev.off()
###############·······················
#
data.sp2<- data.frame()
xm.obs2 <- data.frame()
for(i in 1:49){
  sp1=Inhom.Th(xm2,covr=covr,win=win,sig.t=sig.t,P=alpha)
  data.sp1 <- as.data.frame(sp1); data.sp1$marks <- i
  data.sp2 <-  rbind(data.sp2, data.sp1)
  xm2.1 <- ppp(x=sp1$x,y=sp1$y, window=win)
  xm3 <-envelope(xm2.1,pcf,nsim=11,nrank=5)
  xm.obs1 <- as.data.frame(xm3); xm.obs1$theo <- i
  xm.obs2 <- rbind(xm.obs2, xm.obs1)
}
xm.wide <- xm.wide02 <- reshape(xm.obs2, v.names = "obs", idvar = "r",
                                timevar = "theo", direction = "wide")
 
##########
tiff('Inhomogeneous Thomas EH 28+ cm.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 

plot(xm3.0,legend=FALSE,main="",#xlab="",ylab="", 
     xlim=c(0,80), ylim=c(0,25), main="",bty="l")
matpoints(x=xm.wide$r,y=xm.wide[,-(1:3)], type = "l",lty=1, col=4)
points(x=xm.wide$r,y=xm3.1,col=2,cex=0.5)

legend("top", "Inhomogeneous Thomas Process",lwd=3,
       cex=1.3,text.col=4,col=4,
       inset =0.01,bty="n")
text(50,15,"蚬木 > 28cm")
dev.off()
################################################ 
##########################################################################



#########################################################################
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
xm1 <- subset(xm0,  dbh>=16 & dbh<28 )
#
xm2 <- ppp(x=xm1$x,y=xm1$y, window=win)
xm3.0 <-envelope(xm2,pcf,nsim=11,nrank=5)
xm3.1 <- as.data.frame(xm3.0)$obs
data.sp2<- data.frame()
xm.obs2 <- data.frame()
for(i in 1:49){
  sp1 = Inhom.Po(xm2,covr=covr,win=win,sig.t=sig.t)
  data.sp1 <- as.data.frame(sp1); data.sp1$marks <- i
  data.sp2 <-  rbind(data.sp2, data.sp1)
  xm2.1 <- ppp(x=sp1$x,y=sp1$y, window=win)
  xm3 <-envelope(xm2.1,pcf,nsim=11,nrank=5)
  xm.obs1 <- as.data.frame(xm3); xm.obs1$theo <- i
  xm.obs2 <- rbind(xm.obs2, xm.obs1)
}
xm.wide <- xm.wide11 <- reshape(xm.obs2, v.names = "obs", idvar = "r",
                                timevar = "theo", direction = "wide")
###############

tiff('Inhomogeneous Poisson EH  dbh 16 -- 28 cm.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 

plot(xm3.0,legend=FALSE,main="", 
     xlim=c(0,80), ylim=c(0,25), main="",bty="l")
matpoints(x=xm.wide$r,y=xm.wide[,-(1:3)], type = "l",lty=1, col=3)
points(x=xm.wide$r,y=xm3.1,col=2,cex=0.5)

legend("top", "Inhomogeneous Poisson Process",lwd=3,
       cex=1.3,text.col=3,col=3,
       inset =0.01,bty="n")
text(50,15," dbh>=16 & dbh<28 ")
dev.off()
##########
#
data.sp2<- data.frame()
xm.obs2 <- data.frame()
for(i in 1:49){
  sp1=Inhom.Th(xm2,covr=covr,win=win,sig.t=sig.t,P=alpha)
  data.sp1 <- as.data.frame(sp1); data.sp1$marks <- i
  data.sp2 <-  rbind(data.sp2, data.sp1)
  xm2.1 <- ppp(x=sp1$x,y=sp1$y, window=win)
  xm3 <-envelope(xm2.1,pcf,nsim=11,nrank=5)
  xm.obs1 <- as.data.frame(xm3); xm.obs1$theo <- i
  xm.obs2 <- rbind(xm.obs2, xm.obs1)
}
xm.wide <- xm.wide12 <- reshape(xm.obs2, v.names = "obs", idvar = "r",
                                timevar = "theo", direction = "wide")
##############################################3
tiff('Inhomogeneous Thomas EH  dbh 16 -- 28  cm.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 

plot(xm3.0,legend=FALSE,main="", 
     xlim=c(0,80), ylim=c(0,25), main="",bty="l")
matpoints(x=xm.wide$r,y=xm.wide[,-(1:3)], type = "l",lty=1, col=4)
points(x=xm.wide$r,y=xm3.1,col=2,cex=0.5)

legend("top", "Inhomogeneous Thomas Process",lwd=3,
       cex=1.3,text.col=4,col=4,
       inset =0.01,bty="n")
text(50,15," dbh>=16 & dbh<28 ")
dev.off()
################################################



###################################################################################################################################################
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
xm1 <- subset(xm0,  dbh>=9 & dbh<16 )
#
xm2 <- ppp(x=xm1$x,y=xm1$y, window=win)
xm3.0 <-envelope(xm2,pcf,nsim=11,nrank=5)
xm3.1 <- as.data.frame(xm3.0)$obs
data.sp2<- data.frame()
xm.obs2 <- data.frame()
for(i in 1:49){
  sp1 = Inhom.Po(xm2,covr=covr,win=win,sig.t=sig.t)
  data.sp1 <- as.data.frame(sp1); data.sp1$marks <- i
  data.sp2 <-  rbind(data.sp2, data.sp1)
  xm2.1 <- ppp(x=sp1$x,y=sp1$y, window=win)
  xm3 <-envelope(xm2.1,pcf,nsim=11,nrank=5)
  xm.obs1 <- as.data.frame(xm3); xm.obs1$theo <- i
  xm.obs2 <- rbind(xm.obs2, xm.obs1)
}
xm.wide <- xm.wide011 <- reshape(xm.obs2, v.names = "obs", idvar = "r",
                                timevar = "theo", direction = "wide")
###############

tiff('Inhomogeneous Poisson EH dbh 9 -- 16cm.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 
plot(xm3.0,legend=FALSE,main="",
     xlim=c(0,80), ylim=c(0,25), main="",bty="l")
matpoints(x=xm.wide$r,y=xm.wide[,-(1:3)], type = "l",lty=1, col=3)
points(x=xm.wide$r,y=xm3.1,col=2,cex=0.5)

legend("top", "Inhomogeneous Poisson Process",lwd=3,
       cex=1.3,text.col=3,col=3,
       inset =0.01,bty="n")
text(50,15,"dbh>=9 & dbh<16 cm")
dev.off()
##########
#
data.sp2<- data.frame()
xm.obs2 <- data.frame()
for(i in 1:49){
  sp1=Inhom.Th(xm2,covr=covr,win=win,sig.t=sig.t,P=alpha)
  data.sp1 <- as.data.frame(sp1); data.sp1$marks <- i
  data.sp2 <-  rbind(data.sp2, data.sp1)
  xm2.1 <- ppp(x=sp1$x,y=sp1$y, window=win)
  xm3 <-envelope(xm2.1,pcf,nsim=11,nrank=5)
  xm.obs1 <- as.data.frame(xm3); xm.obs1$theo <- i
  xm.obs2 <- rbind(xm.obs2, xm.obs1)
}
xm.wide <- xm.wide012 <- reshape(xm.obs2, v.names = "obs", idvar = "r",
                                 timevar = "theo", direction = "wide")
#############################################
tiff('Inhomogeneous Thomas EH dbh 9 -- 16 cm.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 
plot(xm3.0,legend=FALSE,main="",
     xlim=c(0,80), ylim=c(0,25), main="",bty="l")
matpoints(x=xm.wide$r,y=xm.wide[,-(1:3)], type = "l",lty=1, col=4)
points(x=xm.wide$r,y=xm3.1,col=2,cex=0.5)
legend("top", "Inhomogeneous Thomas Process",lwd=3,
       cex=1.3,text.col=4,col=4,
       inset =0.01,bty="n")
text(50,15,"dbh>=9 & dbh<16 cm")
dev.off()
################################################ 


###################################################################################################################################################
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
xm1 <- subset(xm0, dbh>=3 & dbh<9)
#
xm2 <- ppp(x=xm1$x,y=xm1$y, window=win)
xm3.0 <-envelope(xm2,pcf,nsim=11,nrank=5)
xm3.1 <- as.data.frame(xm3.0)$obs
data.sp2<- data.frame()
xm.obs2 <- data.frame()
for(i in 1:49){
  sp1 = Inhom.Po(xm2,covr=covr,win=win,sig.t=sig.t)
  data.sp1 <- as.data.frame(sp1); data.sp1$marks <- i
  data.sp2 <-  rbind(data.sp2, data.sp1)
  xm2.1 <- ppp(x=sp1$x,y=sp1$y, window=win)
  xm3 <-envelope(xm2.1,pcf,nsim=11,nrank=5)
  xm.obs1 <- as.data.frame(xm3); xm.obs1$theo <- i
  xm.obs2 <- rbind(xm.obs2, xm.obs1)
}
xm.wide <- xm.wide21 <- reshape(xm.obs2, v.names = "obs", idvar = "r",
                                timevar = "theo", direction = "wide")
###############

tiff('Inhomogeneous Poisson EH db 3 -- 9 cm.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 
plot(xm3.0,legend=FALSE,main="",
     xlim=c(0,80), ylim=c(0,25), main="",bty="l")
matpoints(x=xm.wide$r,y=xm.wide[,-(1:3)], type = "l",lty=1, col=3)
points(x=xm.wide$r,y=xm3.1,col=2,cex=0.5)
legend("top", "Inhomogeneous Poisson Process",lwd=3,
       cex=1.3,text.col=3,col=3,
       inset =0.01,bty="n")
text(50,15,"dbh>=3 & dbh<9")
dev.off()
##########
#
data.sp2<- data.frame()
xm.obs2 <- data.frame()
for(i in 1:49){
  sp1=Inhom.Th(xm2,covr=covr,win=win,sig.t=sig.t,P=alpha)
  data.sp1 <- as.data.frame(sp1); data.sp1$marks <- i
  data.sp2 <-  rbind(data.sp2, data.sp1)
  xm2.1 <- ppp(x=sp1$x,y=sp1$y, window=win)
  xm3 <-envelope(xm2.1,pcf,nsim=11,nrank=5)
  xm.obs1 <- as.data.frame(xm3); xm.obs1$theo <- i
  xm.obs2 <- rbind(xm.obs2, xm.obs1)
}
xm.wide <- xm.wide22 <- reshape(xm.obs2, v.names = "obs", idvar = "r",
                                timevar = "theo", direction = "wide")
##############################################3
tiff('Inhomogeneous Thomas EH dbh 3 -- 9 cm.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 

plot(xm3.0,legend=FALSE,main="",#xlab="",ylab="", 
     xlim=c(0,80), ylim=c(0,25), main="",bty="l")
matpoints(x=xm.wide$r,y=xm.wide[,-(1:3)], type = "l",lty=1, col=4)
points(x=xm.wide$r,y=xm3.1,col=2,cex=0.5)

legend("top", "Inhomogeneous Thomas Process",lwd=3,
       cex=1.3,text.col=4,col=4,
       inset =0.01,bty="n")
text(50,15,"dbh>=3 & dbh<9 cm")
dev.off()
################################################ 



###################################################################################################################################################
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
xm1 <- subset(xm0, dbh<3)
#
xm2 <- ppp(x=xm1$x,y=xm1$y, window=win)
xm3.0 <-envelope(xm2,pcf,nsim=11,nrank=5)
xm3.1 <- as.data.frame(xm3.0)$obs
data.sp2<- data.frame()
xm.obs2 <- data.frame()
for(i in 1:49){
  sp1 = Inhom.Po(xm2,covr=covr,win=win,sig.t=sig.t)
  data.sp1 <- as.data.frame(sp1); data.sp1$marks <- i
  data.sp2 <-  rbind(data.sp2, data.sp1)
  xm2.1 <- ppp(x=sp1$x,y=sp1$y, window=win)
  xm3 <-envelope(xm2.1,pcf,nsim=11,nrank=5)
  xm.obs1 <- as.data.frame(xm3); xm.obs1$theo <- i
  xm.obs2 <- rbind(xm.obs2, xm.obs1)
}
xm.wide <- xm.wide31 <- reshape(xm.obs2, v.names = "obs", idvar = "r",
                                timevar = "theo", direction = "wide")
###############
tiff('Inhomogeneous Poisson EH dbh 0 -- 3 cm.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 

plot(xm3.0,legend=FALSE,main="",#xlab="",ylab="", 
     xlim=c(0,80), ylim=c(0,25), main="",bty="l")
matpoints(x=xm.wide$r,y=xm.wide[,-(1:3)], type = "l",lty=1, col=3)
points(x=xm.wide$r,y=xm3.1,col=2,cex=0.5)

legend("top", "Inhomogeneous Poisson Process",lwd=3,
       cex=1.3,text.col=3,col=3,
       inset =0.01,bty="n")
text(50,15,"dbh<3 cm")
dev.off()
##########
#
data.sp2<- data.frame()
xm.obs2 <- data.frame()
for(i in 1:49){
  sp1=Inhom.Th(xm2,covr=covr,win=win,sig.t=sig.t,P=alpha)
  data.sp1 <- as.data.frame(sp1); data.sp1$marks <- i
  data.sp2 <-  rbind(data.sp2, data.sp1)
  xm2.1 <- ppp(x=sp1$x,y=sp1$y, window=win)
  xm3 <-envelope(xm2.1,pcf,nsim=11,nrank=5)
  xm.obs1 <- as.data.frame(xm3); xm.obs1$theo <- i
  xm.obs2 <- rbind(xm.obs2, xm.obs1)
}
xm.wide <- xm.wide32 <- reshape(xm.obs2, v.names = "obs", idvar = "r",
                                timevar = "theo", direction = "wide")
##############################################3

tiff('Inhomogeneous Thomas EH dbh 0 -- 3 cm.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 

plot(xm3.0,legend=FALSE,main="",#xlab="",ylab="", 
     xlim=c(0,80), ylim=c(0,25), main="",bty="l")
matpoints(x=xm.wide$r,y=xm.wide[,-(1:3)], type = "l",lty=1, col=4)
points(x=xm.wide$r,y=xm3.1,col=2,cex=0.5)

legend("top", "Inhomogeneous Thomas Process",lwd=3,
       cex=1.3,text.col=4,col=4,
       inset =0.01,bty="n")
text(50,15," dbh<3 cm")
dev.off()
################################################ 



###################################################################################################################################################
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
xm1 <- subset(xm0, dbh>0 )
#
xm2 <- ppp(x=xm1$x,y=xm1$y, window=win)
xm3.0 <-envelope(xm2,pcf,nsim=11,nrank=5)
xm3.1 <- as.data.frame(xm3.0)$obs
data.sp2<- data.frame()
xm.obs2 <- data.frame()
for(i in 1:49){
  sp1 = Inhom.Po(xm2,covr=covr,win=win,sig.t=sig.t)
  data.sp1 <- as.data.frame(sp1); data.sp1$marks <- i
  data.sp2 <-  rbind(data.sp2, data.sp1)
  xm2.1 <- ppp(x=sp1$x,y=sp1$y, window=win)
  xm3 <-envelope(xm2.1,pcf,nsim=11,nrank=5)
  xm.obs1 <- as.data.frame(xm3); xm.obs1$theo <- i
  xm.obs2 <- rbind(xm.obs2, xm.obs1)
}
dat.in.po <- data.sp2
write.csv(dat.in.po, 'sample.Inhom.Po.csv')

xm.wide <- xm.wide41 <- reshape(xm.obs2, v.names = "obs", idvar = "r",
                                timevar = "theo", direction = "wide")
###############

tiff('Inhomogeneous Poisson EH all 0+ cm.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 
plot(xm3.0,legend=FALSE,main="",
     xlim=c(0,80), ylim=c(0,25), main="",bty="l")
matpoints(x=xm.wide$r,y=xm.wide[,-(1:3)], type = "l",lty=1, col=3)
points(x=xm.wide$r,y=xm3.1,col=2,cex=0.5)

legend("top", "Inhomogeneous Poisson Process",lwd=3,
       cex=1.3,text.col=3,col=3,
       inset =0.01,bty="n")
text(50,15,"dbh>0 cm")
dev.off()
##########
#
data.sp2<- data.frame()
xm.obs2 <- data.frame()
for(i in 1:49){
  sp1=Inhom.Th(xm2,covr=covr,win=win,sig.t=sig.t,P=alpha)
  data.sp1 <- as.data.frame(sp1); data.sp1$marks <- i
  data.sp2 <-  rbind(data.sp2, data.sp1)
  xm2.1 <- ppp(x=sp1$x,y=sp1$y, window=win)
  xm3 <-envelope(xm2.1,pcf,nsim=11,nrank=5)
  xm.obs1 <- as.data.frame(xm3); xm.obs1$theo <- i
  xm.obs2 <- rbind(xm.obs2, xm.obs1)
}
dat.in.th <- data.sp2
write.csv(dat.in.th, 'sample.Inhom.Po.csv')
xm.wide <- xm.wide42 <- reshape(xm.obs2, v.names = "obs", idvar = "r",
                                timevar = "theo", direction = "wide")
##############################################3
tiff('Inhomogeneous Thomas EH all 0+ cm.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 

plot(xm3.0,legend=FALSE,main="", 
     xlim=c(0,80), ylim=c(0,25), main="",bty="l")
matpoints(x=xm.wide$r,y=xm.wide[,-(1:3)], type = "l",lty=1, col=4)
points(x=xm.wide$r,y=xm3.1,col=2,cex=0.5)

legend("top", "Inhomogeneous Thomas Process",lwd=3,
       cex=1.3,text.col=4,col=4,
       inset =0.01,bty="n")
text(50,15,"dbh>0 cm")
dev.off()
################################################ 
###################################################################################################################################################
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 