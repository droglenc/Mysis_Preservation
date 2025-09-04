

> str(mysis2)
'data.frame':   413 obs. of  28 variables:
 $ include.L : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
 $ include.LW: logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
 $ purpose   : Factor w/ 3 levels "L","L,LW","LW": 3 3 3 3 3 3 3 3 3 3 ...
 $ use.lw    : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
 $ use.l     : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
 $ tech      : Factor w/ 2 levels "EJ","Kate": 2 2 2 2 2 2 2 2 2 2 ...
 $ batch     : Factor w/ 3 levels "1","2","NA": 3 3 3 3 3 3 3 3 3 3 ...
 $ id        : num  1 2 3 4 5 6 8 9 10 11 ...
 $ tin.id    : Factor w/ 191 levels "501","502","503",..: 191 191 191 191 191 191 191 191 191 191 ...
 $ photo.id  : Factor w/ 190 levels "001-3","002-3",..: 189 189 189 189 189 189 189 189 189 189 ...
 $ date      : Factor w/ 6 levels "38995","39000",..: 6 6 6 6 6 6 6 6 6 6 ...
 $ time      : Factor w/ 5 levels "0.125","0.135416666666667",..: 5 5 5 5 5 5 5 5 5 5 ...
 $ station   : Factor w/ 6 levels "171","Grid 1409",..: 4 4 4 4 4 4 4 4 4 4 ...
 $ tx        : Factor w/ 3 levels "8BF","8SBF","Fresh": 3 3 3 3 3 3 2 3 3 3 ...
 $ stage     : Factor w/ 4 levels "Fem","Grav","Juv",..: 1 1 1 4 4 1 1 1 2 2 ...
 $ len.fr.bin: num  13 15 13 13 13 17 NaN 17 15 13 ...
 $ len.fr    : num  14.2 16.1 13.6 13.1 13.2 17 NaN 17.1 15.3 13.7 ...
 $ len.p2w   : Factor w/ 128 levels "10.1","10.2",..: 128 128 128 128 128 128 33 128 128 128 ...
 $ len.p2m   : Factor w/ 44 levels "10.2","10.4",..: 44 44 44 44 44 44 44 44 44 44 ...
 $ raw.wt    : num  0.0035 0.0083 0.0037 0.0036 0.0025 ...
 $ c.wt      : num  3.96 8.76 4.17 4.06 2.97 ...
 $ len       : num  14.2 16.1 13.6 13.1 13.2 17 13.3 17.1 15.3 13.7 ...
 $ len.bin   : num  13 15 13 13 13 17 13 17 15 13 ...
 $ notes     : Factor w/ 29 levels "changed purpose to LW (from L,LW) b/c no fresh length was recorded",..: 14 14 14 14 14 14 14 14 14 14 ...
 $ V25       : logi  NA NA NA NA NA NA ...
 $ logwt     : num  0.598 0.943 0.620 0.609 0.472 ...
 $ loglen    : num  1.15 1.21 1.13 1.12 1.12 ...
 $ treat     : Factor w/ 12 levels "8BF:Fem","8BF:Grav",..: 9 9 9 12 12 9 5 9 10 10 ...

windows(5,5); par(mar=c(3.5,3.5,1,1),mgp=c(2,0.75,0))

# Data entry and cleaning
library(NCStats)
library(xlsReadWrite)                                        # load package to make reading and writing from Excel easier
mysis <- read.xls("Mysis_Data_Master.xls",sheet="Mysis")     # must change directory to where data file is
mysis$c.wt <- mysis$c.wt*1000                                # convert g to mg
mysis1 <- mysis[mysis$include.LW,]                           # excludes values not to be included
mysis2 <- mysis1[mysis1$use.lw==T,]                          # data frame of just individuals for LW analysis
mysis2$logwt <- log10(mysis2$c.wt)                           # transform weight to log10
mysis2$loglen <- log10(mysis2$len)                           # transform length to log10
mysis2$treat <- mysis2$tx:mysis2$stage                       # create an overall treatment variale (makes life easier at times)
attach(mysis2)

# declare some labels
xlbl1 <- "Standard Length (mm)"
xlbl2 <- "log(Standard Length (mm))"
ylbl1 <- "Dry Weight (mg)"
ylbl2 <- "log(Dry Weight (mg))"

# EDA type stuff -- none of this is in the manuscript
ftable(xtabs(~tx+stage+len.bin))                                                   # get sample sizes -- note that gravid-8BF are missing
ftable(xtabs(~tx+stage))
tapply(len,stage:tx,Summary,na.rm=T)                                               # sample summaries
plot(c.wt~len,pch=as.numeric(tx),col=as.numeric(stage),xlab=xlbl1,ylab=ylbl1)      # plot relationship on original scale
legend(x="topleft",legend=levels(tx:stage),col=rep(seq(1,4),3),pch=rep(seq(1,4),each=4))
plot(logwt~loglen,pch=as.numeric(tx),col=as.numeric(stage),xlab=xlbl2,ylab=ylbl2)  # plot relationship on log-log scale -- note increased variance at low values
legend(x="topleft",legend=levels(tx:stage),col=rep(seq(1,4),3),pch=rep(seq(1,4),each=4))

# Permanently excludes juveniles
detach(mysis2)
mysis3 <- mysis2[mysis2$stage!="Juv",] 
mysis3$stage <- factor(mysis3$stage)                         # re-factors so juvenile is not in list of factors
mysis3$treat <- factor(mysis3$treat)                         # re-factors to remove unused treatments
attach(mysis3)
ftable(xtabs(~tx+stage+len.bin))                             # get sample sizes -- double-check that juveniles have been removed -- not in ms
ftable(xtabs(~tx+stage))                                     # sample size -- in ms
tapply(len,stage:tx,Summary,na.rm=T)                         # summary statistics -- in ms
tapply(c.wt,stage:tx,Summary,na.rm=T)                        # summary statistics -- in ms
plot(logwt~loglen,pch=as.numeric(tx),col=as.numeric(stage),xlab=xlbl2,ylab=ylbl2)          # plot relationship on log-log scale -- not in ms
legend(x="topleft",legend=levels(tx:stage),col=rep(seq(1,3),3),pch=rep(seq(1,3),each=3))

# Temporarily excludes gravids
detach(mysis3)
mysis3.nograv <- mysis3[mysis3$stage!="Grav",]               # excludes gravids
mysis3.nograv$stage <- factor(mysis3.nograv$stage)           # re-factors so gravid is not in list of factors
mysis3.nograv$treat <- factor(mysis3.nograv$treat)           # re-factors to remove unused treatments
attach(mysis3.nograv)
ftable(xtabs(~tx+stage+len.bin))                             # get sample sizes -- double-check that gravids have been removed -- not in ms

# FIRST LM -- Analysis Male and Females among 8BF, 8SBF, and Fresh
lm1 <- lm(logwt~loglen*tx*stage)
Anova(lm1)                                                   # slope differences due to treatments -- in ms
comp.slopes(lm1)                                             # in ms
plot.fit(lm1,leg="topleft",main="")

ad.test(lm1$residuals)
ncv.test(lm1)
plot.resid(lm1,main="")
hist(lm1$residuals,main="")

detach(mysis3.nograv)
mysis3.nograv.nofresh <- mysis3.nograv[mysis3.nograv$tx!="Fresh",]
mysis3.nograv.nofresh$tx <- factor(mysis3.nograv.nofresh$tx)
attach(mysis3.nograv.nofresh)
lm1a <- lm(logwt~loglen+tx)                                 # compare intercepts for 8BF and 8SBF combined for M & F -- in ms
comp.intercepts(lm1a)                                       # in ms

# SECOND LM -- Analysis Male, Female, Gravids among 8SBF and Fresh
detach(mysis3.nograv.nofresh)
mysis3.no8BF <- mysis3[mysis3$tx!="8BF",]
mysis3.no8BF$stage <- factor(mysis3.no8BF$stage)            # re-factors so gravid is not in list of factors
mysis3.no8BF$tx <- factor(mysis3.no8BF$tx)                  # re-factors to remove unused treatments
mysis3.no8BF$treat <- factor(mysis3.no8BF$treat)            # re-factors to remove unused treatments
attach(mysis3.no8BF)
ftable(xtabs(~tx+stage+len.bin))

lm2 <- lm(logwt~loglen*tx*stage)
Anova(lm2)                                                 # slope differences due to treatments, intercept due to stage -- in ms
comp.slopes(lm2)                                           # in ms
plot.fit(lm2,leg="topleft",main="")

ad.test(lm2$residuals)
ncv.test(lm2)
plot.resid(lm2,main="")
hist(lm2$residuals,main="")

lm2a <- lm(logwt~loglen+stage)
Anova(lm2a)                                                # in ms
comp.intercepts(lm2a)                                      # in ms

detach(mysis3.no8BF)
# Finding regression relationships for different groups
mysis3.fresh <- mysis3[mysis3$tx=="Fresh",]                                    # just the fresh tx
mysis3.fresh$stage1 <- factor(ifelse(mysis3.fresh$stage=="Grav","Grav","MF"))  # combine males and females
mysis3.fresh$tx <- factor(mysis3.fresh$tx)
attach(mysis3.fresh)
lm3.fresh <- lm(logwt~loglen+stage1)
summary(lm3.fresh)                                         # in ms
plot.fit(lm3.fresh,main="",leg="topleft")                  # not in ms
detach(mysis3.fresh)

mysis3.nfresh <- mysis3[mysis3$tx!="Fresh",]                                    # not fresh
mysis3.nfresh$stage1 <- factor(ifelse(mysis3.nfresh$stage=="Grav","Grav","MF")) # combine males and females
mysis3.fresh$tx <- factor(mysis3.fresh$tx)
mysis3.nfresh$treat1 <- factor(mysis3.nfresh$stage1:mysis3.nfresh$tx)
attach(mysis3.nfresh)
lm3.nfresh <- lm(logwt~loglen+treat1)
summary(lm3.nfresh)                                        # in ms
comp.intercepts(lm3.nfresh,0)                              # not in ms
plot.fit(lm3.nfresh,main="",leg="topleft")                 # not in ms
detach(mysis3.nfresh)


# THE PLOT OF INDIVIDUAL POINTS and REGRESSION LINES -- FIGURE 2
mysis4 <- mysis3
mysis4$stage1 <- factor(ifelse(mysis4$stage=="Grav","Grav","MF"))
attach(mysis4)
windows(5,5); par(mar=c(3,3,1,1),mgp=c(1.75,0.75,0))
plot(c.wt~len,log="xy",col="white",axes=FALSE,xlab=xlbl1,ylab=ylbl1,xlim=c(8,20),ylim=c(0.8,20))                      # plot base
axis(2,at=c(0.8,0.9,seq(1,10,by=1),20),labels=F)                                                                      # y-axis ticks
axis(2,at=c(1,5,10,20),labels=T)                                                                                      # y-axis labels
axis(1,at=c(8:20),labels=F)                                                                                           # x-axis ticks
axis(1,at=c(8,10,15,20),labels=T)                                                                                     # x-axis labels
data <- subset(mysis4,stage1=="MF" & tx=="Fresh"); points(data$c.wt~data$len,col="black",pch=19)
data <- subset(mysis4,stage1=="Grav" & tx=="Fresh"); points(data$c.wt~data$len,col="black",pch=17,cex=1.25)
data <- subset(mysis4,stage1=="MF" & tx=="8BF"); points(data$c.wt~data$len,col="black",pch=0)
data <- subset(mysis4,stage1=="MF" & tx=="8SBF"); points(data$c.wt~data$len,col="black",pch=1)
data <- subset(mysis4,stage1=="Grav" & tx=="8SBF"); points(data$c.wt~data$len,col="black",pch=2)

lines(c(8.1,19.5),10^c(2.7659*log10(8.1)-2.4304,2.7659*log10(19.5)-2.4304),lwd=3,col="black",lty=1)
lines(c(12.9,19.5),10^c(2.7659*log10(12.9)-2.3719,2.7659*log10(19.5)-2.3719),lwd=3,col="black",lty=3)
lines(c(8.9,19.5),10^c(3.0658*log10(8.9)-2.7178,3.0658*log10(19.5)-2.7178),lwd=3,col="gray",lty=2)
lines(c(8.2,19.5),10^c(3.0658*log10(8.2)-2.5819,3.0658*log10(19.5)-2.5819),lwd=3,col="gray",lty=1)
lines(c(13.7,19.5),10^c(3.0658*log10(13.7)-2.5579,3.0658*log10(19.5)-2.5579),lwd=3,col="gray",lty=3)
legend(x=14,y=1.5,legend=c("Fresh,MF","Fresh,Grav"),lty=c(1,3),lwd=3,inset=0.02)
legend(x=14,y=1.5,legend=c("",""),pch=c(19,17),pt.cex=1.25,inset=0.02,bty="n")
legend("topleft",legend=c("8BF,MF","8SBF,MF","8SBF,Grav"),lty=c(2,1,3),lwd=3,col="gray",inset=0.02)
legend("topleft",legend=c("","",""),pch=c(0,1,2),pt.cex=1.25,inset=0.02,bty="n")



# THE PLOT OF INDIVIDUAL POINTS and REGRESSION LINES -- THIS HAS TWO PLOTS & SCALED SO AXES SHOW SAME ORDER OF MAGNITUDE CHANGE.
mysis4 <- mysis3
mysis4$stage1 <- factor(ifelse(mysis4$stage=="Grav","Grav","MF"))
attach(mysis4)
windows(8,15); par(mar=c(3,3,1,1),mgp=c(1.75,0.75,0),mfrow=c(1,2))
plot(c.wt~len,log="xy",col="white",axes=FALSE,xlab=xlbl1,ylab=ylbl1,xlim=c(8,20))
axis(2,at=c(seq(0.1,0.9,by=0.1),seq(1,10,by=1),20),labels=F)
axis(2,at=c(0.1,1,10,20),labels=T)
axis(1,at=c(8:20),labels=F)                                                                                           # x-axis ticks
axis(1,at=c(8,10,15,20),labels=T)
data <- subset(mysis4,stage1=="MF" & tx=="Fresh"); points(data$c.wt~data$len,col="black",pch=19)
data <- subset(mysis4,stage1=="Grav" & tx=="Fresh"); points(data$c.wt~data$len,col="black",pch=17,cex=1.25)
data <- subset(mysis4,stage1=="MF" & tx=="8BF"); points(data$c.wt~data$len,col="black",pch=0)
data <- subset(mysis4,stage1=="MF" & tx=="8SBF"); points(data$c.wt~data$len,col="black",pch=1)
data <- subset(mysis4,stage1=="Grav" & tx=="8SBF"); points(data$c.wt~data$len,col="black",pch=2)
legend("bottomright",legend=c("Fresh,MF","Fresh,Grav","8BF,MF","8SBF,MF","8SBF,Grav"),
       pch=c(19,17,0,1,2),pt.cex=1.25,inset=0.02,cex=0.8)

plot(c.wt~len,log="xy",col="white",axes=FALSE,xlab=xlbl1,ylab=ylbl1,xlim=c(8,20))
axis(2,at=c(seq(0.1,0.9,by=0.1),seq(1,10,by=1),20),labels=F)
axis(2,at=c(0.1,1,10,20),labels=T)
axis(1,at=c(8:20),labels=F)                                                                                           # x-axis ticks
axis(1,at=c(8,10,15,20),labels=T)
lines(c(8.1,19.5),10^c(2.7659*log10(8.1)-2.4304,2.7659*log10(19.5)-2.4304),lwd=3,col="black",lty=1)
lines(c(12.9,19.5),10^c(2.7659*log10(12.9)-2.3719,2.7659*log10(19.5)-2.3719),lwd=3,col="black",lty=3)
lines(c(8.9,19.5),10^c(3.0658*log10(8.9)-2.7178,3.0658*log10(19.5)-2.7178),lwd=3,col="gray",lty=2)
lines(c(8.2,19.5),10^c(3.0658*log10(8.2)-2.5819,3.0658*log10(19.5)-2.5819),lwd=3,col="gray",lty=1)
lines(c(13.7,19.5),10^c(3.0658*log10(13.7)-2.5579,3.0658*log10(19.5)-2.5579),lwd=3,col="gray",lty=3)
legend("bottomright",legend=c("Fresh,MF","Fresh,Grav","8BF,MF","8SBF,MF","8SBF,Grav"),
       col=c("black","black","gray","gray","gray"),lty=c(1,3,2,1,3),lwd=3,inset=0.02,cex=0.8)
