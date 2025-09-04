# The structure for reference purposes
> str(mysis.len)
'data.frame':   188 obs. of  28 variables:
 $ include   : logi  TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE ...
 $ purpose   : Factor w/ 3 levels "L","L,LW","LW": 1 1 1 1 1 1 1 1 1 1 ...
 $ use.lw    : logi  FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE ...
 $ use.l     : logi  TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE ...
 $ tech      : Factor w/ 2 levels "EJ","Kate": 2 2 2 2 2 2 2 2 2 2 ...
 $ batch     : int  NA NA NA NA NA NA NA NA NA NA ...
 $ id        : int  303 306 309 311 321 327 331 333 334 337 ...
 $ tin.id    : int  NA NA NA NA NA NA NA NA NA NA ...
 $ photo.id  : Factor w/ 137 levels "523fresh_male",..: NA NA NA NA NA NA NA NA NA NA ...
 $ date      : Factor w/ 5 levels "10/10/2006","10/11/2006",..: NA NA NA NA NA NA NA NA NA NA ...
 $ time      : Factor w/ 4 levels "19:26","21:17",..: NA NA NA NA NA NA NA NA NA NA ...
 $ station   : Factor w/ 5 levels "171","Grid 1409",..: NA NA NA NA NA NA NA NA NA NA ...
 $ tx        : Factor w/ 2 levels "8BF","8SBF": 2 2 2 2 2 2 2 2 2 2 ...
 $ stage     : Factor w/ 4 levels "Fem","Grav","Juv",..: 1 1 1 1 4 1 4 1 1 4 ...
 $ len.fr.bin: int  11 11 9 15 13 11 13 11 9 11 ...
 $ len.fr    : num  11 12.4 10.9 16.5 13.9 11.7 14.2 12.2 10 11.4 ...
 $ len.p2w   : num  10.7 11.8 9.9 15 13 10.6 12.1 11.6 9.7 11.3 ...
 $ len.p2m   : num  10.6 11.8 9.9 15 13 10.6 12.1 11.6 9.7 11.3 ...
 $ loss.p2w  : num  -0.3 -0.6 -1.0 -1.5 -0.9 ...
 $ loss.p2m  : num  -0.4 -0.6 -1.0 -1.5 -0.9 ...
 $ ploss.p2w : num  -0.0273 -0.0484 -0.0917 -0.0909 -0.0647 ...
 $ ploss.p2m : num  -0.0364 -0.0484 -0.0917 -0.0909 -0.0647 ...

# load needed libraries
library(NCStats)

# Data entry, cleaning, and variable creation
mysis <- read.csv("Mysis.csv",head=T)
mysis1 <- mysis[mysis$include,]                            # excludes values not to be included
mysis.len <- mysis1[mysis1$use.l,]                         # just a data frame for length analysis
mysis.len$tx <- factor(mysis.len$tx)                       # this is needed to remove the "Fresh" level which is not included in this analysis.
mysis.len$stage <- factor(mysis.len$stage)                 # this is needed to remove the "Gravid" level which is not included in this analysis.
mysis.len$group <- mysis.len$stage:mysis.len$tx            # creates a group variable that is a combination of the stage and treatment.
mysis.len$loss.p2w <- mysis.len$len.p2w-mysis.len$len.fr   # computes loss in length for 2 weeks preservation (negative means loss)
mysis.len$loss.p2m <- mysis.len$len.p2m-mysis.len$len.fr   # same but for 2 months preservation
mysis.len$ploss.p2w <- mysis.len$loss.p2w/mysis.len$len.fr # computes proportional loss in length for 2 weeks perservation
mysis.len$ploss.p2m <- mysis.len$loss.p2m/mysis.len$len.fr # same but for 2 months preservation
attach(mysis.len)

# Sample Summaries
ftable(xtabs(~tx+stage))                                   # get sample sizes
ftable(xtabs(~tx+stage+len.fr.bin))
tapply(len.fr,stage:tx,Summary,na.rm=T)                    # sample summaries

# initial model fits on absolute loss
lm1 <- lm(loss.p2w~len.fr*stage*tx)
Anova(lm1)                                                 # slopes differ by treatment and stage -- combined to a one-way ANOVA below
plot.fit(lm1,legend="topleft")

lm1a <- lm(loss.p2w~len.fr*group)
Anova(lm1a)
plot.fit(lm1a,legend="topleft")
comp.slopes(lm1a)

# initial model fits on proportional loss
lm2 <- lm(ploss.p2w~len.fr*stage*tx)
Anova(lm2)
plot.fit(lm2,legend="topleft")

lm2a <- lm(ploss.p2w~len.fr+stage+tx+len.fr:stage)
Anova(lm2a)
plot.fit(lm2a,legend="topleft")

   # removed two outliers -- #18 & 22
detach(mysis.len)
mysis.len2 <- mysis.len[-c(18,22),]
attach(mysis.len2)
lm2b <- lm(ploss.p2w~len.fr*stage*tx)
Anova(lm2b)
plot.fit(lm2b,legend="topleft")

  # tried without juveniles
detach(mysis.len2)
mysis.len3 <- mysis.len[mysis.len$stage!="Juv",]
mysis.len3$stage <- factor(mysis.len3$stage)
mysis.len3$group <- factor(mysis.len3$group)
attach(mysis.len3)

lm2c <- lm(loss.p2w~len.fr*stage*tx)
Anova(lm2c)
plot.fit(lm2c,legend="topleft")

lm2d <- lm(ploss.p2w~len.fr*stage*tx)
Anova(lm2d)
plot.fit(lm2d,legend="topleft")


lm3 <- lm(loss.p2w~len.fr*group)
Anova(lm3)
plot.fit(lm3)

   # tried without juveniles and all individuals less than 9.5
detach(mysis.len3)
mysis.len4 <- mysis.len3[mysis.len3$len.fr>=10,]
attach(mysis.len4)
lm4 <- lm(loss.p2w~len.fr*stage*tx)
Anova(lm4)
plot.fit(lm4,legend="topleft")



lm1a <- lm(loss.p2w~len.fr*group)
Anova(lm1a)
plot.fit(lm1a,legend="topleft")
comp.slopes(lm1a)
