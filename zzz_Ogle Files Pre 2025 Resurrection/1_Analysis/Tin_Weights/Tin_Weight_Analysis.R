`data.frame':   30 obs. of  6 variables:
 $ batch   : int  1 1 1 1 1 1 1 1 1 1 ...
 $ tin     : int  1 2 3 4 5 6 7 8 9 10 ...
 $ tx      : Factor w/ 2 levels "1h","24h": 2 2 2 2 2 2 2 2 2 2 ...
 $ fresh.wt: num  1.13 1.13 1.13 1.13 1.14 ...
 $ driet.wt: num  1.13 1.13 1.13 1.13 1.14 ...
 $ wt.loss : num  4e-04 7e-04 7e-04 8e-04 8e-04 9e-04 9e-04 8e-04 8e-04 9e-04 ...

library(NCStats)
tins <- read.table("Tins.txt",head=T)
tins$fbatch <- factor(tins$batch)
attach(tins)
boxplot(wt.loss~batch)
boxplot(wt.loss~tx)

# ANALYSIS BY BATCHES
# Comparing tin weight losses across batches -- all significantly different
tinbatch.lm <- lm(wt.loss~fbatch)
Levenes(tinbatch.lm)
ad.test(tinbatch.lm$residuals)
plot.resid(tinbatch.lm,cex.main=0.8,loess=F)
outlier.test(tinbatch.lm)
anova(tinbatch.lm)
plot.means(tinbatch.lm,n.label=F,xlab="Batch",ylab="Weight Loss (mg)")
TukeyHSD(tinbatch.lm)
add.sig.letters(tinbatch.lm,c("a","b","c"),pos=c(4,4,4),col=c("red","red","red"))
tapply(wt.loss,fbatch,mean)

# Comparing initial tin weights across batches -- No significant difference
tininit.lm <- lm(fresh.wt~fbatch)
Levenes(tininit.lm)
ad.test(tininit.lm$residuals)
plot.resid(tininit.lm,cex.main=0.8,loess=F)
outlier.test(tininit.lm)
anova(tininit.lm)
plot.means(tininit.lm,n.label=F,xlab="Batch",ylab="Initial Weight (mg)")
TukeyHSD(tininit.lm)

# Comparing dried tin weights across batches -- Slight difference between 3 and 1
tindried.lm <- lm(driet.wt~fbatch)
Levenes(tindried.lm)
ad.test(tindried.lm$residuals)
plot.resid(tindried.lm,cex.main=0.8,loess=F)
outlier.test(tindried.lm)
anova(tindried.lm)
plot.means(tindried.lm,n.label=F,xlab="Batch",ylab="Dried Weight (mg)")
TukeyHSD(tindried.lm)
add.sig.letters(tindried.lm,c("a","ab","b"),pos=c(4,4,4),col=c("red","red","red"))


# ANALYSIS BY TREATMENTS
# Comparing tin weight losses across batches -- Strongly significantly different
tinloss1.lm <- lm(wt.loss~tx)
Levenes(tinloss1.lm)
ad.test(tinloss1.lm$residuals)
plot.resid(tinloss1.lm,cex.main=0.8,loess=F)
outlier.test(tinloss1.lm)
anova(tinloss1.lm)
plot.means(tinloss1.lm,n.label=F,xlab="Treatment",ylab="Weight Loss (mg)")
TukeyHSD(tinloss1.lm)
add.sig.letters(tinloss1.lm,c("a","b"),pos=c(4,4),col=c("red","red"))
tapply(wt.loss,tx,mean)

# Comparing initial tin weights across treatments --Slight difference
tininit1.lm <- lm(fresh.wt~tx)
Levenes(tininit1.lm)
ad.test(tininit1.lm$residuals)
plot.resid(tininit1.lm,cex.main=0.8,loess=F)
outlier.test(tininit1.lm)
anova(tininit1.lm)
plot.means(tininit1.lm,n.label=F,xlab="Treatment",ylab="Initial Weight (mg)")
TukeyHSD(tininit1.lm)

# Comparing dried tin weights across treatments -- Moderate difference
tindried1.lm <- lm(driet.wt~tx)
Levenes(tindried1.lm)
ad.test(tindried1.lm$residuals)
plot.resid(tindried1.lm,cex.main=0.8,loess=F)
outlier.test(tindried1.lm)
anova(tindried1.lm)
plot.means(tindried1.lm,n.label=F,xlab="Treatment",ylab="Dried Weight (mg)")
TukeyHSD(tindried1.lm)
add.sig.letters(tindried1.lm,c("a","b"),pos=c(4,4),col=c("red","red"))
tapply(driet.wt,tx,mean)


######
# Tin weight loss correction method 1 -- [(DM+DT)-FT]+(FT-DT)
# add 0.00077 to batch 1 (mean weight loss of batch 1)
# add 0.00048 to batch 2 (mean weight loss of batch 2)
# add 0.000465 to batch 3 (difference in 1 h loss (mean of batch 3) and 24 h loss (i.e., mean of batch 1 and 2))


#####
# Tin weight loss correction method 2 -- (DM+DT) - DTbar
# assume constant 24-h dry tin weight of 1.131565 (there was no statistical difference between batch 1 and 2 final dry weights)
# subtract this from each dry Mysis weight + tin value to get a dry Mysis weight
