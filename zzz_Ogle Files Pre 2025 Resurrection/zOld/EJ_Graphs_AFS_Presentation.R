# Data entry and cleaning
library(NCStats)
mysis <- read.csv("Mysis.csv",head=T)
mysis$c.wt <- mysis$c.wt*1000                                # convert mg to ug
# mysis$tx <- relevel(mysis$tx,"Fresh")                        # make fresh the reference group
mysis1 <- mysis[mysis$include,]                              # excludes values not to be included
mysis2 <- mysis1[mysis1$use.lw==T,]                          # data frame of just individuals for LW analysis
mysis2$logwt <- log10(mysis2$c.wt)                           # transform weight to log10
mysis2$loglen <- log10(mysis2$len)                           # transform length to log10
mysis2$treat <- mysis2$tx:mysis2$stage                       # create an overall treatment variale (makes life easier at times)

mysis2.juv <- mysis2[mysis2$stage=="Juv",]
mysis2.grav <- mysis2[mysis2$stage=="Grav",]
mysis2.grav$tx <- factor(mysis2.grav$tx,levels=c("Fresh","8SBF"))
mysis2.ad <- mysis2[mysis2$stage=="Male" | mysis2$stage=="Fem",]
mysis2.ad$stage <- factor(mysis2.ad$stage)

# declare some labels
xlbl1 <- "Standard Length (mm)"
xlbl2 <- "log(Standard Length (mm))"
ylbl1 <- expression(paste("Weight (",mu,"g)"))
ylbl2 <- expression(paste("log(Weight (",mu,"g))"))

# Male and Female
windows(8,8); par(mar=c(4,4,1,1),mgp=c(2.5,0.75,0),cex=1.5,cex.lab=1.25)
attach(mysis2.ad)
tx <- factor(tx,levels=c("8SBF","8BF","Fresh"))
lm5ad1 <- lm(logwt~loglen*tx*stage)
plot.fit(lm5ad1,xlab=xlbl2,ylab=ylbl2,main="",pts=c(1,19),clrs=c("white","white","black"),xlim=c(0.9,1.25),ylim=c(-0.03,1.28))
legend(x="topleft",legend=c("","","","","Fresh:Fem","Fresh:Male",),pch=c(1,19,1,19,1,19),col=c("white","white","white","white","black","black"),lty=c(1,2,1,2,1,2))
tx <- factor(tx,levels=c("8SBF","Fresh","8BF"))
lm5ad2 <- lm(logwt~loglen*tx*stage)
plot.fit(lm5ad2,xlab=xlbl2,ylab=ylbl2,main="",pts=c(1,19),clrs=c("white","gray","red"),xlim=c(0.9,1.25),ylim=c(-0.03,1.28))
legend(x="topleft",legend=c("","","Fresh:Fem","Fresh:Male","8BF:Fem","8BF:Male"),pch=c(1,19,1,19,1,19),col=c("white","white","gray","gray","red","red"),lty=c(1,2,1,2,1,2))
tx <- factor(tx,levels=c("Fresh","8BF","8SBF"))
lm5ad3 <- lm(logwt~loglen*tx*stage)
plot.fit(lm5ad3,xlab=xlbl2,ylab=ylbl2,main="",pts=c(1,19),clrs=c("lightgray","gray","blue"),xlim=c(0.9,1.25),ylim=c(-0.03,1.28))
legend(x="topleft",legend=c("Fresh:Fem","Fresh:Male","8BF:Fem","8BF:Male","8SBF:Fem","8SBF:Male"),pch=c(1,19,1,19,1,19),col=c("lightgray","lightgray","gray","gray","blue","blue"),lty=c(1,2,1,2,1,2))
plot.fit(lm5ad3,xlab=xlbl2,ylab=ylbl2,main="",pts=c(1,19),clrs=c("black","red","blue"),legend="topleft")

lm5ad4 <- lm(logwt~loglen+tx)
plot.fit(lm5ad4,xlab=xlbl2,ylab=ylbl2,main="",pts=c(1,19),clrs=c("black","red","blue"),legend="topleft",xlim=c(0.9,1.25),ylim=c(-0.03,1.28),plot.pts=F)
detach(mysis2.ad)
rm(tx)

# GRAVIDS
attach(mysis2.grav)
tx <- factor(tx,levels=c("8SBF","Fresh"))
lm5grav1 <- lm(logwt~loglen*tx)
plot.fit(lm5grav1,legend="topleft",xlab=xlbl2,ylab=ylbl2,main="",pts=19,clrs=c("white","black"),xlim=c(0.9,1.25),ylim=c(-0.03,1.28))
plot.fit(lm5grav1,legend="topleft",xlab=xlbl2,ylab=ylbl2,main="",pts=19,clrs=c("blue","black"),xlim=c(0.9,1.25),ylim=c(-0.03,1.28))

lm5grav2 <- lm(logwt~loglen+tx)
plot.fit(lm5grav2,legend="topleft",xlab=xlbl2,ylab=ylbl2,main="",pts=19,clrs=c("blue","black"),xlim=c(0.9,1.25),ylim=c(-0.03,1.28),plot.pts=F)
detach(mysis2.grav)
rm(tx)


# JUVENILES
attach(mysis2.juv)
tx <- factor(tx,levels=c("8SBF","8BF","Fresh"))
lm5juv1 <- lm(logwt~loglen*tx)
plot.fit(lm5juv1,xlab=xlbl2,ylab=ylbl2,main="",pts=19,clrs=c("white","white","black"))
legend(x="topleft",legend=c("","","Fresh"),pch=19,col=c("white","white","black"),lty=1)
tx <- factor(tx,levels=c("8SBF","Fresh","8BF"))
lm5juv2 <- lm(logwt~loglen*tx)
plot.fit(lm5juv2,xlab=xlbl2,ylab=ylbl2,main="",pts=19,clrs=c("white","gray","red"))
legend(x="topleft",legend=c("","Fresh","8BF"),pch=19,col=c("white","gray","red"),lty=1)
tx <- factor(tx,levels=c("Fresh","8BF","8SBF"))
lm5juv3 <- lm(logwt~loglen*tx)
plot.fit(lm5juv3,xlab=xlbl2,ylab=ylbl2,main="",pts=19,clrs=c("lightgray","gray","blue"))
legend(x="topleft",legend=c("Fresh","8BF","8SBF"),pch=19,col=c("lightgray","gray","blue"),lty=1)
plot.fit(lm5juv3,xlab=xlbl2,ylab=ylbl2,main="",pts=19,clrs=c("black","red","blue"),legend="topleft")

lm5juv4 <- lm(logwt~loglen+tx)
plot.fit(lm5juv4,xlab=xlbl2,ylab=ylbl2,main="",pts=19,clrs=c("black","red","blue"),legend="topleft",plot.pts=F)
detach(mysis2.juv)
rm(tx)
