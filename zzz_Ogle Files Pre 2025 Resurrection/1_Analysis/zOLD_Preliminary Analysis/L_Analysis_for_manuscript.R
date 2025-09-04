library(NCStats)
library(FSA)

# Reads first raw data file -- each row is an individual with 3 measurements
Lmysis1 <- read.table("L_Mysis1.txt",head=T)
Lmysis1$diff.p2w <- Lmysis1$len.p2w - Lmysis1$len.fr
Lmysis1$diff.p2m <- Lmysis1$len.p2m - Lmysis1$len.fr
Lmysis1$diff.pr  <- Lmysis1$len.p2m - Lmysis1$len.p2w
Lmysis1$chng.p2w <- ifelse(Lmysis1$diff.p2w<0,"Decrease",ifelse(Lmysis1$diff.p2w==0,"Same","Increase"))
Lmysis1$chng.p2m <- ifelse(Lmysis1$diff.p2m<0,"Decrease",ifelse(Lmysis1$diff.p2m==0,"Same","Increase"))
Lmysis1$chng.pr <- ifelse(Lmysis1$diff.pr<0,"Decrease",ifelse(Lmysis1$diff.pr==0,"Same","Increase"))
Lmysis1.pdiff.p2w <- Lmysis1$diff.p2w/Lmysis1$len.fr
Lmysis1.pdiff.p2m <- Lmysis1$diff.p2m/Lmysis1$len.fr
Lmysis1.pdiff.pr <- Lmysis1$diff.pr/Lmysis1$len.fr
Lmysis1 <- lencat(Lmysis1,"len.fr",2,3,vname="LCat.fr")  # creates length category variable -- 2-mm bins starting at 3 mm

# Reads second raw data file -- each row is an individual with 2 measurements -- 2 week and 2 month separated
Lmysis2 <- read.table("L_Mysis2.txt",head=T)
Lmysis2$diff <- Lmysis2$len.pr - Lmysis2$len.fr
Lmysis2$pdiff <- Lmysis2$diff/Lmysis2$len.fr

# EDA, not used in manuscript
  # length frequency of sample -- EDA, not used in manuscript
  table(Lmysis1$stage,Lmysis1$LCat.fr)

  # plot of all three measurements by ID
  attach(Lmysis1)
  plot(len.fr~id,pch=19)
  points(len.p2w~id,pch=19,col="red")
  points(len.p2m~id,pch=19,col="blue")
  
  # table of how fish changed due to preservation
  prop.table(table(chng.p2w))
  prop.table(table(chng.p2m))
  prop.table(table(chng.pr))

  # plots of change versus length
  plot(diff.p2w~len.fr,pch=19,cex=1.5)
  points(diff.p2m~len.fr,pch=19,col="gray")
  
  plot(diff.pr~len.frm)
  abline(h=0)
  detach(Lmysis1)

  # Summary statistics
  summary(diff.p2w)
  summary(diff.p2m)
  

# Model exploreation, not used in manuscript
  attach(Lmysis2)
  lm1 <- lm(diff~len.fr*stage*tx)
  plot.fit(lm1)
  plot.resid(lm1)
  
  lm2 <- lm(pdiff~len.fr*stage*tx)
  plot.fit(lm2)
  plot.resid(lm2)
