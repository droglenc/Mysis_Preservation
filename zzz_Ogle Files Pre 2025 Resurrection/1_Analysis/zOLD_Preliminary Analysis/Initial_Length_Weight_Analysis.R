> str(mysis2)
`data.frame':   414 obs. of  26 variables:
 $ include   : logi  TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE ...
 $ purpose   : Factor w/ 3 levels "L","L,LW","LW": 2 3 3 3 3 3 3 2 3 3 ...
 $ use.lw    : logi  TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE ...
 $ use.l     : logi   TRUE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE ...
 $ tech      : Factor w/ 2 levels "EJ","Kate": 1 2 2 1 2 2 2 1 2 2 ...
 $ tx        : Factor w/ 3 levels "8BF","8SBF","Fresh": 2 2 3 3 3 3 2 2 2 3 ...
 $ stage     : Factor w/ 4 levels "Fem","Grav","Juv",..: 3 3 3 3 3 3 3 3 3 3 ...
 $ c.wt      : num  0.78 0.63 0.53 0.55 0.27 0.27 0.67 1.01 1.06 0.47 ...
 $ len       : num  4 4 4.8 5 5 5.1 5.3 5.4 5.6 5.6 ...      <-- this is the length to use
 $ len.bin   : int  3 3 3 5 5 5 5 5 5 5 ...
 $ logwt     : num  -0.108 -0.201 -0.276 -0.260 -0.569 ...
 $ loglen    : num  0.602 0.602 0.681 0.699 0.699 ...

# Data entry and cleaning
library(NCStats)
mysis <- read.csv("Mysis.csv",head=T)
mysis$c.wt <- mysis$c.wt*1000                                # convert mg to ug
mysis1 <- mysis[mysis$include,]                              # excludes values not to be included
mysis2 <- mysis1[mysis1$use.lw==T,]                          # data frame of just individuals for LW analysis
mysis2$logwt <- log10(mysis2$c.wt)                           # transform weight to log10
mysis2$loglen <- log10(mysis2$len)                           # transform length to log10
mysis2$treat <- mysis2$tx:mysis2$stage                       # create an overall treatment variale (makes life easier at times)
attach(mysis2)

# declare some labels
xlbl1 <- "Standard Length (mm)"
xlbl2 <- "log(Standard Length (mm))"
ylbl1 <- expression(paste("Weight (",mu,"g)"))
ylbl2 <- expression(paste("log(Weight (",mu,"g))"))

# EDA type stuff
ftable(xtabs(~tx+stage+len.bin))                             # get sample sizes -- note that gravid-8BF are missing
tapply(len,stage:tx,Summary,na.rm=T)                         # sample summaries
plot(c.wt~len,pch=as.numeric(tx),col=as.numeric(stage),xlab=xlbl1,ylab=ylbl1)      # plot relationship on original scale
legend(x="topleft",legend=levels(tx:stage),col=rep(seq(1,4),3),pch=rep(seq(1,4),each=4))
plot(logwt~loglen,pch=as.numeric(tx),col=as.numeric(stage),xlab=xlbl2,ylab=ylbl2)  # plot relationship on log-log scale -- not increased variance at low values
legend(x="topleft",legend=levels(tx:stage),col=rep(seq(1,4),3),pch=rep(seq(1,4),each=4))

#=============================================================================================================================================
# Examined relationships for non-gravids -- all three treatments
detach(mysis2)
mysis2.nograv <- mysis2[mysis2$stage!="Grav",]               # excludes gravids
mysis2.nograv$stage <- factor(mysis2.nograv$stage)           # re-factors so gravid is not in list of factors
mysis2.nograv$treat <- factor(mysis2.nograv$treat)           # re-factors so gravid is not in list of factors
attach(mysis2.nograv)
ftable(xtabs(~tx+stage+len.bin))                             # get sample sizes -- double-check that gravids have been removed

lm1 <- lm(logwt~loglen*tx*stage)
Anova(lm1)                                                   # slope differences observed in three-way interaction

lm2 <- lm(logwt~loglen*treat)                                # treat as a one-way IVR to determine where slope differences are
Anova(lm2) 
comp.slopes(lm2)                                             # post-hoc analysis on slopes
plot.fit(lm2,legend="topleft")


# all of the slope differences were related to 8SBF:Juv (NOT different than Fresh:Juv, 8BF:Male, 8BF:Juv)
# concerned that the two smallest individuals (with seemingly large weights) may have depressed the slope (8SBF:Juv hadd shallowest slope)
# tried removing these individuals to determine their effect
detach(mysis2.nograv)
mysis2.nograv2 <- mysis2.nograv[-c(81,348),]
attach(mysis2.nograv2)
lm3 <- lm(logwt~loglen*tx*stage)
Anova(lm3) 

# the three-way interaction disappeared without those two individuals
# slopes still differed by stage, but not by treatment
lm4 <- lm(logwt~loglen+tx+stage+loglen:stage)                # reduced to significant factors
Anova(lm4)                                                   # slopes still differ by stage
summary(lm4)                                                 # get J-F and M-F slope p-values
stage1 <- relevel(stage,"Juv")
lm4a <- lm(logwt~loglen+tx+stage1+loglen:stage1)
summary(lm4a)                                                # get J-M slope p-value
fdr.control(c(0.00154,0.81717,0.00492))                      # FDR control for p-values for J-F,M-F,M-J
plot.fit(lm4a,legend="topleft")

# slope for juveniles are significantly shallower than both males and females -- slope for males and females is not different
# reduced to analysis of males and females, analysis of juveniles, and analysis of gravids
# put 2 small juveniles removed back in for now
detach(mysis2.nograv2)
mysis2.juv <- mysis2[mysis2$stage=="Juv",]
mysis2.grav <- mysis2[mysis2$stage=="Grav",]
mysis2.grav$tx <- factor(mysis2.grav$tx)
mysis2.ad <- mysis2[mysis2$stage=="Male" | mysis2$stage=="Fem",]
mysis2.ad$stage <- factor(mysis2.ad$stage)

#=============================================================================================================================================
# male and female analysis -- all three treatments
attach(mysis2.ad)
lm5ad <- lm(logwt~loglen*tx*stage)
Anova(lm5ad)
# slopes are different among treatments -- probably because variability was reduced when juvs were removed
lm5ada <- lm(logwt~loglen+tx+stage+loglen:tx)
Anova(lm5ada)
summary(lm5ada)
tx1 <- relevel(tx,"8SBF")
lm5adb <- lm(logwt~loglen+tx1+stage+loglen:tx1)
summary(lm5adb)
fdr.control(c(0.7529,0.0450,0.0207))                         # FDR control for p-values for 8SBF-8BF,FRESH-8BF,FRESH-8SBF
plot.fit(lm5adb,legend="topleft")
#q-values are 0.06 indicating only a weak significance ....
lm5adc <- lm(logwt~loglen+tx+stage)
Anova(lm5adc)
plot.fit(lm5adc,legend="topleft")
# no stage effect
lm5add <- lm(logwt~loglen+tx)
Anova(lm5add)
plot.fit(lm5add,legend="topleft")

#=============================================================================================================================================
# juvenile analysis -- all three treatments
detach(mysis2.ad)
attach(mysis2.juv)
lm5juv <- lm(logwt~loglen*tx)
Anova(lm5juv)
comp.intercepts(lm5juv)
plot.fit(lm5juv,legend="topleft")

#=============================================================================================================================================
# gravid analysis -- all three treatments
detach(mysis2.juv)
attach(mysis2.grav)
lm5grav <- lm(logwt~loglen*tx)
Anova(lm5grv)
comp.intercepts(lm5grav)
plot.fit(lm5grav,legend="topleft")
detach(mysis2.grav)



#=============================================================================================================================================
# all but juveniles -- 8SBF and Fresh treatments
mysis3 <- mysis2[mysis2$stage!="Juv",]
mysis4 <- mysis3[mysis3$tx!="8BF",]
mysis4$tx <- factor(mysis4$tx)
mysis4$stage <- factor(mysis4$stage)
attach(mysis4)
ftable(xtabs(~tx+stage+len.bin))                             # get sample sizes -- double-check
lm6 <- lm(logwt~loglen*tx*stage)
Anova(lm6)                                                   # only suggests a slope difference among treatments (not stages)
lm6a <- lm(logwt~loglen+tx+stage+loglen:tx)
plot.fit(lm6,legend="topleft")


#q-values are 0.06 indicating only a weak significance ....
lm5adc <- lm(logwt~loglen+tx+stage)
Anova(lm5adc)
plot.fit(lm5adc,legend="topleft")
# no stage effect
lm5add <- lm(logwt~loglen+tx)
Anova(lm5add)
plot.fit(lm5add,legend="topleft")








# This is old stuff --------------------------------------------------------------------------------
library(NCStats)
mysis.orig <- read.table("LW_Mysis4.txt",head=T)              # reads raw data file
mysis.orig$wt <- mysis.orig$wt*1000                           # convert mg to ug
mysis.1 <- mysis.orig[-which(is.na(mysis.orig$wt)),]          # removes indivs with missing wts
mysis.2 <- mysis.1[-which(mysis.1$len>5 & mysis.1$len<7),]    # removes indivs between 5 & 7 mm
mysis.3 <- mysis.2[-c(17,23,41,171),]                         # removed outliers from dataframe
mysis.3$logwt <- log10(mysis.3$wt)
mysis.3$loglen <- log10(mysis.3$len)

# THE PLOT OF INDIVIDUAL POINTS -- FIGURE 1.
attach(mysis.3)
windows(8,8); par(mar=c(3,3,1,1),mgp=c(1.75,0.75,0),mfrow=c(1,2))
plot(wt~len,log="xy",col="white",axes=FALSE,xlab="Standard Length (mm)",ylab="Dry Weight (mg)",xlim=c(2,20))
axis(2,at=c(seq(0.1,0.9,by=0.1),seq(1,10,by=1),20),labels=F)
axis(2,at=c(0.1,1,10,20),labels=T)
axis(1,at=c(2:10,20),labels=F)
axis(1,at=c(5,10,15,20),labels=T)
data <- mysis.3[mysis.3$stage=="Juv" & mysis.3$tx=="Fr",]; points(data$wt~data$len,col="black",pch=5)
data <- mysis.3[mysis.3$stage=="Juv" & mysis.3$tx=="Pr",]; points(data$wt~data$len,col="black",pch=18,cex=1.25)
data <- mysis.3[mysis.3$stage=="Male" & mysis.3$tx=="Fr",]; points(data$wt~data$len,col="black",pch=1,cex=1.25)
data <- mysis.3[mysis.3$stage=="Male" & mysis.3$tx=="Pr",]; points(data$wt~data$len,col="black",pch=16,cex=1.25)
data <- mysis.3[mysis.3$stage=="Fem" & mysis.3$tx=="Fr",]; points(data$wt~data$len,col="black",pch=2)
data <- mysis.3[mysis.3$stage=="Fem" & mysis.3$tx=="Pr",]; points(data$wt~data$len,col="black",pch=17,cex=1.25)
data <- mysis.3[mysis.3$stage=="Grav" & mysis.3$tx=="Fr",]; points(data$wt~data$len,col="black",pch=0)
data <- mysis.3[mysis.3$stage=="Grav" & mysis.3$tx=="Pr",]; points(data$wt~data$len,col="black",pch=15,cex=1.25)
legend("topleft",legend=c("Fr,Juvenile","Pr,Juvenile","Fr,Male","Pr,Male","Fr,Female","Pr,Female","Fr,Gravid","Pr,Gravid"),
       pch=c(5,18,1,16,2,17,0,15),pt.cex=c(1,1.25,1.25,1.25,1,1.25,1,1.25),inset=0.02)
lines(c(10,20,20,10,10),c(2,2,20,20,2),lty=3,lwd=2)

plot(wt~len,log="xy",col="white",axes=FALSE,xlab="Standard Length (mm)",ylab="Dry Weight (mg)",xlim=c(10,20),ylim=c(2,20))
axis(2,at=seq(2,20,by=1),labels=F)
axis(2,at=c(2,10,20),labels=T)
axis(1,at=10:20,labels=F)
axis(1,at=c(10,15,20),labels=T)
data <- mysis.3[mysis.3$stage=="Male" & mysis.3$tx=="Fr",]; points(data$wt~data$len,col="black",pch=1,cex=1.25)
data <- mysis.3[mysis.3$stage=="Male" & mysis.3$tx=="Pr",]; points(data$wt~data$len,col="black",pch=16,cex=1.25)
data <- mysis.3[mysis.3$stage=="Fem" & mysis.3$tx=="Fr",]; points(data$wt~data$len,col="black",pch=2)
data <- mysis.3[mysis.3$stage=="Fem" & mysis.3$tx=="Pr",]; points(data$wt~data$len,col="black",pch=17,cex=1.25)
data <- mysis.3[mysis.3$stage=="Grav" & mysis.3$tx=="Fr",]; points(data$wt~data$len,col="black",pch=0)
data <- mysis.3[mysis.3$stage=="Grav" & mysis.3$tx=="Pr",]; points(data$wt~data$len,col="black",pch=15,cex=1.25)

# THE SAMPLE CHARACTERISTICS -- TABLE 1.
tapply(len,stage:tx,length)
tapply(len,stage:tx,summary)
tapply(len,stage:tx,sd)

# THE OVERALL MODEL FIT
lm1 <- lm(logwt~loglen*tx*stage)
Anova(lm1)

lm1a <- lm(logwt~loglen+tx+stage+loglen:stage)                 # took out insignificant terms
Anova(lm1a)
summary(lm1a)                                                  # find differences in slope from female line
stage1 <- relevel(stage,"Grav")
lm1a1 <- lm(logwt~loglen+tx+stage1+loglen:stage1)
summary(lm1a1)                                                 # find differences in slope from gravid line
stage1 <- relevel(stage,"Juv")
lm1a2 <- lm(logwt~loglen+tx+stage1+loglen:stage1)
summary(lm1a2)                                                 # find differences in slope from Juvenile line
fdr.control(c(0.2210,0.0162,0.2307,0.8282,0.0717,0.0033))      # FDR control for p-values for G-F,J-F,M-F,G-J,G-M,M-J

# JUST JUVENILES
  # Separate into juveniles and non-juveniles
  detach(mysis.3)
  mysis3.juv <- mysis.3[mysis.3$stage=="Juv",]
  mysis3.njuv <- mysis.3[mysis.3$stage!="Juv",]
  mysis3.njuv$stage <- factor(mysis3.njuv$stage)               # removes "Juv" stage which is not needed

 # analysis of juveniles
 attach(mysis3.juv)
 lm1.juv <- lm(logwt~loglen*tx)
 Anova(lm1.juv)
 lm1a.juv <- lm(logwt~loglen+tx)                               # took out insignificant terms
 summary(lm1a.juv)
 confint(lm1a.juv)

# NON-JUVENILES
 detach(mysis3.juv)
 attach(mysis3.njuv)
 lm1.njuv <- lm(logwt~loglen*tx*stage)
 Anova(lm1.njuv)
 lm1a.njuv <- lm(logwt~loglen+tx+stage+tx:stage)               # took out insignificant terms
 Anova(lm1a.njuv)
 
 detach(mysis3.njuv)
 mysis3.njuv$group <- factor(mysis3.njuv$tx:mysis3.njuv$stage) # convert to a one-way ANOVA (combine two factors into one)
 attach(mysis3.njuv)
 lm2.njuv <- lm(logwt~loglen*group)
 Anova(lm2.njuv)
 
 lm2a.njuv <- lm(logwt~loglen+group)                           # took out insignificant terms
 Anova(lm2a.njuv)
 comp.intercepts(lm2a.njuv)
 
# THE PLOT OF INDIVIDUAL REGRESSIONS -- Figure 2
 detach(mysis3.njuv)
 mysis3.frjuv <- mysis.3[mysis.3$stage=="Juv" & mysis.3$tx=="Fr",]
 mysis3.frmale <- mysis.3[mysis.3$stage=="Male" & mysis.3$tx=="Fr",]
 mysis3.frfem <- mysis.3[mysis.3$stage=="Fem" & mysis.3$tx=="Fr",]
 mysis3.frgrav <- mysis.3[mysis.3$stage=="Grav" & mysis.3$tx=="Fr",]
 mysis3.prjuv <- mysis.3[mysis.3$stage=="Juv" & mysis.3$tx=="Pr",]
 mysis3.prmale <- mysis.3[mysis.3$stage=="Male" & mysis.3$tx=="Pr",]
 mysis3.prfem <- mysis.3[mysis.3$stage=="Fem" & mysis.3$tx=="Pr",]
 mysis3.prgrav <- mysis.3[mysis.3$stage=="Grav" & mysis.3$tx=="Pr",]

 lm.frjuv <- lm(logwt~loglen,data=mysis3.frjuv)
 lm.frmale <- lm(logwt~loglen,data=mysis3.frmale)
 lm.frfem <- lm(logwt~loglen,data=mysis3.frfem)
 lm.frgrav <- lm(logwt~loglen,data=mysis3.frgrav)
 lm.prjuv <- lm(logwt~loglen,data=mysis3.prjuv)
 lm.prmale <- lm(logwt~loglen,data=mysis3.prmale)
 lm.prfem <- lm(logwt~loglen,data=mysis3.prfem)
 lm.prgrav <- lm(logwt~loglen,data=mysis3.prgrav)

 windows(8,8); par(mar=c(3,3,1,1),mgp=c(1.75,0.75,0),mfrow=c(1,2))
 plot(wt~len,data=mysis.3,log="xy",col="white",axes=FALSE,xlab="Standard Length (mm)",ylab="Dry Weight (mg)",xlim=c(2,20))
 axis(2,at=c(seq(0.1,0.9,by=0.1),seq(1,10,by=1),20),labels=F)
 axis(2,at=c(0.1,1,10,20),labels=T)
 axis(1,at=c(2:10,20),labels=F)
 axis(1,at=c(5,10,15,20),labels=T)

 lines(c(3.7,8.6),10^predict(lm.frjuv,data.frame(loglen=log10(c(3.7,8.6)))),lwd=3,col="gray",lty=4)
 lines(c(3.3,8.6),10^predict(lm.prjuv,data.frame(loglen=log10(c(3.3,8.6)))),lwd=3,col="black",lty=4)
 lines(c(8.1,15.8),10^predict(lm.frmale,data.frame(loglen=log10(c(8.1,15.8)))),lwd=3,col="gray",lty=1)
 lines(c(12.4,15.1),10^predict(lm.prmale,data.frame(loglen=log10(c(12.4,15.1)))),lwd=3,col="black",lty=1)
 lines(c(8.6,18.3),10^predict(lm.frfem,data.frame(loglen=log10(c(8.6,18.3)))),lwd=3,col="gray",lty=2)
 lines(c(8.2,18.0),10^predict(lm.prfem,data.frame(loglen=log10(c(8.2,18.0)))),lwd=3,col="black",lty=2)
 lines(c(12.9,17.9),10^predict(lm.frgrav,data.frame(loglen=log10(c(12.9,17.9)))),lwd=3,col="gray",lty=3)
 lines(c(13.7,17.9),10^predict(lm.prgrav,data.frame(loglen=log10(c(13.7,17.9)))),lwd=3,col="black",lty=3)

 legend("topleft",legend=c("Fr,Juvenile","Pr,Juvenile","Fr,Male","Pr,Male","Fr,Female","Pr,Female","Fr,Gravid","Pr,Gravid"),
       lty=c(4,4,1,1,2,2,3,3),col=c("gray","black","gray","black","gray","black","gray","black"),lwd=2,inset=0.02)
 lines(c(10,20,20,10,10),c(2,2,20,20,2),lty=3,lwd=2)

 plot(wt~len,data=mysis.3,log="xy",col="white",axes=FALSE,xlab="Standard Length (mm)",ylab="Dry Weight (mg)",xlim=c(10,20),ylim=c(2,20))
 axis(2,at=seq(2,20,by=1),labels=F)
 axis(2,at=c(2,10,20),labels=T)
 axis(1,at=10:20,labels=F)
 axis(1,at=c(10,15,20),labels=T)
 lines(c(8.1,15.8),10^predict(lm.frmale,data.frame(loglen=log10(c(8.1,15.8)))),lwd=3,col="gray",lty=1)
 lines(c(12.4,15.1),10^predict(lm.prmale,data.frame(loglen=log10(c(12.4,15.1)))),lwd=3,col="black",lty=1)
 lines(c(8.6,18.3),10^predict(lm.frfem,data.frame(loglen=log10(c(8.6,18.3)))),lwd=3,col="gray",lty=2)
 lines(c(8.2,18.0),10^predict(lm.prfem,data.frame(loglen=log10(c(8.2,18.0)))),lwd=3,col="black",lty=2)
 lines(c(12.9,17.9),10^predict(lm.frgrav,data.frame(loglen=log10(c(12.9,17.9)))),lwd=3,col="gray",lty=3)
 lines(c(13.7,17.9),10^predict(lm.prgrav,data.frame(loglen=log10(c(13.7,17.9)))),lwd=3,col="black",lty=3) 

 
# THE TABLE OF MODEL VALUES
 summary(lm.frjuv); confint(lm.frjuv)
 summary(lm.frmale); confint(lm.frmale)
 summary(lm.frfem); confint(lm.frfem)
 summary(lm.frgrav); confint(lm.frgrav)
 summary(lm.prjuv); confint(lm.prjuv)
 summary(lm.prmale); confint(lm.prmale)
 summary(lm.prfem); confint(lm.prfem)
 summary(lm.prgrav); confint(lm.prgrav) 
 
# OVERALL REGRESSION -- NO STAGE EFFECT, JUST PRESERVATION
lm.oall <- lm(logwt~loglen*tx)
Anova(lm.oall)
lm.oall2 <- lm(logwt~loglen+tx)
summary(lm.oall2)
confint(lm.oall2)

 detach(mysis.3)
 mysis3.fr <- mysis.3[mysis.3$tx=="Fr",]
 mysis3.pr <- mysis.3[mysis.3$tx=="Pr",]
 
 lm.fr <- lm(logwt~loglen,data=mysis3.fr)
 summary(lm.fr); confint(lm.fr)
 lm.pr <- lm(logwt~loglen,data=mysis3.pr)
 summary(lm.pr); confint(lm.r) 
 

 
 
