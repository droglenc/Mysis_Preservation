library(NCStats)
mysis.orig <- read.table("Mysis3.txt",head=T)           # reads data file
mysis.all <- mysis.orig[-which(is.na(mysis.orig$wt)),]  # removes indivs with missing wts

attach(mysis.all)
lm1 <- lm(wt~len*tx*stage)                              # full model
plot.fit(lm1,legend=T)                                  # observed data with model fit -- click on plot to place legend
detach(mysis.all)

mysis1 <- mysis.all[-201,]                              # remove an obvious outliers
attach(mysis1)
lm2 <- lm(wt~len*tx*stage)
plot.fit(lm2,legend=T) 
trans.chooser(lm2,starty=0.001)                         # explore transformations
detach(mysis1)

mysis.juv <- mysis1[mysis1$stage=="Juv",]               # get just juveniles
mysis.juv$stage <- factor(mysis.juv$stage)              # resets levels for stage
mysis.njuv <- mysis1[mysis1$stage!="Juv",]              # get just non-juveniles
mysis.njuv$stage <- factor(mysis.njuv$stage)            # resets levels for stage

attach(mysis.njuv)                                      # explore non-juveniles
lm3 <- lm(wt~len*tx*stage)
plot.fit(lm3,legend=T) 
trans.chooser(lm3,starty=0.001)
Anova(lm3)
detach(mysis.njuv)

mysis.njuv1 <- mysis.njuv[mysis.njuv$len>11,]           # explore non-juveniles that are longer than 11 mm (compare to above analysis)
attach(mysis.njuv1)
lm4 <- lm(wt~len*tx*stage)
plot.fit(lm4,legend=T) 
trans.chooser(lm4,starty=0.001)
Anova(lm4)
detach(mysis.njuv1)

attach(mysis.juv)                                       # explore juveniles
lm5 <- lm(wt~len*tx)
plot.fit(lm5,legend=T)
trans.chooser(lm5,starty=0.0005)
detach(mysis.juv)
