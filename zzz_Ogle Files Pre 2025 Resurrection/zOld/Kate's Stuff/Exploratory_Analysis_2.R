library(NCStats)
mysis.orig <- read.table("Mysis3.txt",head=T)                 # reads raw data file
mysis.orig$wt <- mysis.orig$wt*1000                           # convert mg to ug
mysis.1 <- mysis.orig[-which(is.na(mysis.orig$wt)),]          # removes indivs with missing wts
mysis.2 <- mysis.1[-which(mysis.1$len>5 & mysis.1$len<7),]    # removes indivs between 5 & 7 mm
plot(mysis.2$wt~mysis.2$len)
identify(mysis.2$wt~mysis.2$len)                              # identified outliers on raw plot
mysis.2[c(17,23,41,171),]                                     # info on outliers -- rows in command, IDs shown below
#     ID stage tx  len      wt
#19   19   Fem Fr  9.6 -0.0008
#40   40   Juv Fr  5.0 -0.0002
#62   62   Juv Pr  7.9  0.0000
#208 201   Fem Pr 11.6  0.0138

mysis.3 <- mysis.2[-c(17,23,41,171),]                         # removed outliers from dataframe

# OVERALL MODEL -- NOT BY STAGES ##############
# Fit LW model with just preservation effect, no stage effect
mysis.3$lnwt <- log(mysis.3$wt)
mysis.3$lnlen <- log(mysis.3$len)
attach(mysis.3)
lm1 <- lm(wt~len*tx)
plot.fit(lm1,legend=T) 
trans.chooser(lm1)                   # log-log transformation is not perfect -- some heteroscedasticity & some non-normality

lm2 <- lm(lnwt~lnlen*tx, data=mysis.3)
plot.fit(lm2,legend=T) 
plot.resid(lm2)
Anova(lm2)                           # slopes are the same, intercepts are different
summary(lm2)                         # use this to find separate equations
detach(mysis.3)

# Separate into fresh and preserved samples to get individual fits
mysis3.fr <- mysis.3[mysis.3$tx=="Fr",]
mysis3.pr <- mysis.3[mysis.3$tx=="Pr",]
lm2.fr <- lm(lnwt~lnlen,data=mysis3.fr)
summary(lm2.fr)
confint(lm2.fr)
lm2.pr <- lm(lnwt~lnlen,data=mysis3.pr)
summary(lm2.pr)
confint(lm2.pr)

# OVERALL MODEL -- INCLUDING STAGES ################
attach(mysis.3)
lm3 <- lm(wt~len*tx*stage)
plot.resid(lm3)
ad.test(lm3$residuals)
hist(lm3$residuals)                   # log-log transformation is not perfect -- some heteroscedasticity & some non-normality

lm4 <- lm(lnwt~lnlen*tx*stage)
plot.resid(lm4)
ad.test(lm4$residuals)
hist(lm4$residuals)
Anova(lm4)

lm5 <- lm(lnwt~lnlen+tx+stage+lnlen:stage)
Anova(lm5)
summary(lm5)                                   # looks like juveniles have different slope -- need to test this with mult comp

  # Separate into juveniles and non-juveniles
  detach(mysis.3)
  mysis3.juv <- mysis.3[mysis.3$stage=="Juv",]
  mysis3.njuv <- mysis.3[mysis.3$stage!="Juv",]
  mysis3.njuv$stage <- factor(mysis3.njuv$stage)  # removes "Juv" stage which is not needed

  # Compare juveniles by treatment -- same slope, different intercepts
  attach(mysis3.juv)
  lm4.juv <- lm(lnwt~lnlen*tx)
  plot.resid(lm4.juv)
  ad.test(lm4.juv$residuals)
  hist(lm4.juv$residuals)
  Anova(lm4.juv)

  lm4a.juv <- lm(lnwt~lnlen+tx)
  Anova(lm4a.juv)
  summary(lm4a.juv)
  
    # Separate JUVENILES into fresh and preserved
    detach(mysis3.juv)
    mysis3juv.fr <- mysis3.juv[mysis3.juv$tx=="Fr",]
    mysis3juv.pr <- mysis3.juv[mysis3.juv$tx=="Pr",]
  
    lm4juv.fr <- lm(lnwt~lnlen,data=mysis3juv.fr)
    summary(lm4juv.fr)
    confint(lm4juv.fr)
    lm4juv.pr <- lm(lnwt~lnlen,data=mysis3juv.pr)
    summary(lm4juv.pr)
    confint(lm4juv.pr)

  # Compare non-juveniles by stage and treatment
  attach(mysis3.njuv)
  lm4.njuv <- lm(lnwt~lnlen*tx*stage)
  plot.resid(lm4.njuv)
  ad.test(lm4.njuv$residuals)
  hist(lm4.njuv$residuals)
  Anova(lm4.njuv)
