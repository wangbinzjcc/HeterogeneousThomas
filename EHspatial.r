###############################
library("mgcv")
library("spatstat")

########
setwd('F:/DataW/lg-data/composition')
dir()
xm0 <- read.csv('lgdat2013.csv')
xm0 <- subset(xm0,sp=='蚬木' & is.na(bra)==T )
if(any(xm0$dbh < 1)){xm0$dbh[xm0$dbh<1]=1}
w00　<- duplicated(paste(xm0$x,xm0$y))
xm0$x[w00] <- xm0$x[w00] + runif(sum(w00),0,0.1)
head(xm0)
###
tiff('DBH rank--abundance.tiff',
     width = 3000, height = 2600,res=600,compression = "lzw")
par( mar=c(5.5,6,0.5,0.5 ),pty="m",mex=0.6) 
 barplot(tab <- table(cut(xm0$dbh,c(1,3,10,20, Inf))), ylim=c(0,900),
         ylab="个体数 abundance", xlab="胸径等级 DBH rank (cm)",cex.lab=1.2)
 text(c(1-0.3, 2, 3+0.1, 4+0.3, 5+0.5),tab+40, labels=tab,cex=1.2)
dev.off()
 

tiff('DBH rank-- bare area.tiff',
     width = 3000, height = 2600,res=600,compression = "lzw")
par( mar=c(5.5,6,0.5,0.5 ),pty="m",mex=0.6) 
barplot(tab <- tapply(pi*xm0$dbh^2/4,cut(xm0$dbh, c(1,3,10,20, Inf)), sum), 
        ylim=c(0,98000),
        ylab=expression('断面积之和  basal area ' * ( cm^2 ) ), 
        xlab="胸径等级 DBH rank (cm)",cex.lab=1.2)
text(c(1-0.3, 2, 3+0.1, 4+0.3, 5+0.5),tab+2500, labels=round(tab), cex=1.2)
dev.off()
#
#####
#
xm1 <- subset(xm0, dbh>0)
xm2<-ppp(xm1$x,xm1$y,marks=as.factor(xm1$sp),c(0,500),c(0,300))
xm3a<-envelope(xm2,pcf,nsim=199,nrank=5)

#
xm1 <- subset(xm0, dbh<3)
xm2<-ppp(xm1$x,xm1$y,marks=as.factor(xm1$sp),c(0,500),c(0,300))
xm3<-envelope(xm2,pcf, nsim=199, nrank=5)

#
xm1=subset(xm0, dbh>=3 & dbh <9)
xm2<-ppp(xm1$x,xm1$y,marks=as.factor(xm1$sp),c(0,500),c(0,300))
xm3.0<-envelope(xm2,pcf,nsim=199,nrank=5)

#
xm1=subset(xm0, dbh>=9 & dbh<16)
xm2<-ppp(xm1$x,xm1$y,marks=as.factor(xm1$sp),c(0,500),c(0,300))
xm3.1<-envelope(xm2,pcf,nsim=199,nrank=5)

#
xm1=subset(xm0, dbh>=16 & dbh<28 )
xm2<-ppp(xm1$x,xm1$y,marks=as.factor(xm1$sp),c(0,500),c(0,300))
xm3.2<-envelope(xm2,pcf,nsim=599,nrank=5)

#
xm1=subset(xm0, dbh>=28)
xm2<-ppp(xm1$x,xm1$y,marks=as.factor(xm1$sp),c(0,500),c(0,300))
xm3.3<-envelope(xm2,pcf,nsim=199,nrank=5)
##
#
##
tiff('xm.0.inf.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 
   plot(xm3a,legend=FALSE,xlab="",ylab="", 
        xlim=c(0,80), ylim=c(0,25), main="",bty="l")
legend("top", "蚬木（全部）", inset =0.01,bty="n")
dev.off()
##
tiff('xm.0.3.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 
   plot(xm3,legend=FALSE,xlab="",ylab="", 
     xlim=c(0,80),ylim=c(0,25), main="",bty="l")
   legend("top", "蚬木 1-3 cm", inset =0.01,bty="n")
dev.off()
#########
tiff('xm.3.9.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 
   plot(xm3.0,legend=FALSE,xlab="",ylab="", 
     xlim=c(0,80), ylim=c(0,25), main="",bty="l")
legend("top", "蚬木 3-9 cm", inset =0.01,bty="n")
dev.off()
#########
tiff('xm.9.16.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 
   plot(xm3.1,legend=FALSE,xlab="",ylab="",
        ylim=c(0,25),xlim=c(0,80), main="",bty="l")
    legend("top", "蚬木 9-16 cm", inset =0.01,bty="n")
dev.off()
##########

tiff('xm.16.28.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 
  plot(xm3.2,legend=FALSE,xlab="",ylab="",
       ylim=c(0,25), xlim=c(0,80),main="",bty="l")
    legend("top", "蚬木 16-28 cm", inset =0.01,bty="n")
dev.off()
#########

tiff('xm.28.inf.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 
plot(xm3.3,legend=FALSE,xlab="",ylab="", 
     ylim=c(0,25), xlim=c(0,80), main="",bty="l")
    legend("top", "蚬木 > 28cm", inset =0.01,bty="n")
dev.off()
#########################
#
#########################
#
data<-as.data.frame(xm3a)
data <-na.omit(data)
sum(data$obs-data$hi > 0)
write.csv(data,file="xm.0.inf.csv")
#
data<-as.data.frame(xm3)
data <-na.omit(data)
sum(data$obs-data$hi > 0)
write.csv(data,file="xm.0.3.csv")
#

data<-as.data.frame(xm3.0)
data <-na.omit(data)
sum(data$obs-data$hi > 0)
write.csv(data,file="xm.3.9.csv")

#
data<-as.data.frame(xm3.1)
data <-na.omit(data)
sum(data$obs-data$hi > 0)
write.csv(data,file="xm.9.16.csv")

#
data<-as.data.frame(xm3.2)
data <-na.omit(data)
sum(data$obs-data$hi > 0)
write.csv(data,file="xm.16.28.csv")
#
data<-as.data.frame(xm3.3)
data <-na.omit(data)
sum(data$obs-data$hi > 0)
write.csv(data,file="xm.28.csv")

##############################


#####################################################
c00 <- cut(xm0$dbh,c(0,3,9,16,28,Inf))
m00 <- match(c00,sort(unique(c00)))
sp0 <- c('A','B',"C","D",'E')[m00]
#############################

xm3<-ppp(xm0$x,xm0$y,marks=as.factor(sp0),c(0,500),c(0,300))

#############################################
xianbi<-envelope(xm3,fun=pcfcross,nsim=199,i="A",j="B",nrank=5)

plot(xianbi,legend=FALSE,xlab="",ylab="",main="",bty="l")
legend("top", "蚬木A-B", inset =0.01,bty="n")
#
#
xianbi.ac <- envelope(xm3,fun=pcfcross,nsim=199,i="A",j="C",nrank=5)

plot(xianbi.ac,legend=FALSE,xlab="",ylab="",main="",bty="l")
legend("top", "蚬木A-C", inset =0.01,bty="n")
#
#
xianbi.ad<-envelope(xm3,fun=pcfcross,nsim=199,i="A",j="D",nrank=5)

plot(xianbi.ad,legend=FALSE,xlab="",ylab="",main="",bty="l")
 legend("top", "蚬木A-D", inset =0.01,bty="n")
#
#
xianbi.ae<-envelope(xm3,fun=pcfcross,nsim=199,i="A",j="E",nrank=5)

plot(xianbi.ae,legend=FALSE,xlab="",ylab="",main="",bty="l")
legend("top", "蚬木A-E", inset =0.01,bty="n")
#
#
xianbi<-envelope(xm3,fun=pcfcross,nsim=199,i="A",j="B",nrank=5)

plot(xianbi,legend=FALSE,xlab="",ylab="",main="",bty="l")
legend("top", "蚬木A-B", inset =0.01,bty="n")

##########################################

#
xianbi.bc <- envelope(xm3,fun=pcfcross,nsim=199,i="B",j="C",nrank=5)

plot(xianbi.bc,legend=FALSE,xlab="",ylab="",main="",bty="l")
legend("top", "蚬木B-C", inset =0.01,bty="n")
#
#
xianbi.bd<-envelope(xm3,fun=pcfcross,nsim=199,i="B",j="D",nrank=5)

plot(xianbi.bd,legend=FALSE,xlab="",ylab="",main="",bty="l")
legend("top", "蚬木B-D", inset =0.01,bty="n")
#
#
xianbi.be<-envelope(xm3,fun=pcfcross,nsim=199,i="B",j="E",nrank=5)

plot(xianbi.be,legend=FALSE,xlab="",ylab="",main="",bty="l")
legend("top", "蚬木B-E", inset =0.01,bty="n")
#
#
##########################################
 
#
xianbi.cd<-envelope(xm3,fun=pcfcross,nsim=199,i="C",j="D",nrank=5)

plot(xianbi.cd,legend=FALSE,xlab="",ylab="",main="",bty="l")
legend("top", "蚬木C-D", inset =0.01,bty="n")
#
#
xianbi.ce<-envelope(xm3,fun=pcfcross,nsim=199,i="C",j="E",nrank=5)

plot(xianbi.ce,legend=FALSE,xlab="",ylab="",main="",bty="l")
legend("top", "蚬木C-E", inset =0.01,bty="n")
##########


xianbi.de<-envelope(xm3,fun=pcfcross,nsim=199,i="D",j="E",nrank=5)

plot(xianbi.de,legend=FALSE,xlab="",ylab="",main="",bty="l")
legend("top", "蚬木D-E", inset =0.01,bty="n")
#


