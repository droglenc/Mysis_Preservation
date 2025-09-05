# ===== Setup =====
# ----- Load packages
library(tidyverse)
library(gt)
library(emmeans)

# ----- Set some options
options(show.signif.stars=FALSE)
theme_set(theme_bw())

# ----- Set some constants
stage.lvls1 <- c("Male","Fem")          # to match order in ms
stage.lvls2 <- c("Male","Fem","Grav")   # to match order in ms
tx.lvls2 <- c("Fresh","8BF","8SBF")     # to match order in ms

# ===== Prepare Main data.frame =====
#  o convert tx to factor
#  o convert stage to factor
#  o remove variables not used in this analysis (for conciseness)
mysis <- readxl::read_xls("Mysis_Data_Master.xls",sheet="Mysis",
                          na=c("","NA","na")) |>
  mutate(tx=factor(tx),
         stage=factor(stage)) |>
  select(-purpose,-(tech:station),-len.p2m,-raw.wt,-notes)
#  o checks
str(mysis)
dim(mysis)
levels(mysis$tx)
levels(mysis$stage)




# ===== Length Analysis =====
# ----- Data.frame for length analysis
#  o include only individuals to be used for length analysis
#  o exclude values that were problematic for some reason (see notes)
#  o remove variables not used further (for conciseness)
#  o refactor tx, which will remove "Fresh" as an option
#  o refactor stage, which will remove "Grav" as an option
#  o create group that is combo of the stage and treatment
#  o compute loss in len for 2 weeks preservation (negative means loss)
#  o compute proportional loss in len for 2 weeks preservation
#  o categorize "loss" type
mysisL <- mysis |>
  filter(include.L) |>
  filter(use.l) |>
  select(-(include.L:use.l),-c.wt) |>
  mutate(tx=factor(tx),
         stage=factor(stage),
         group=tx:stage,
         loss.p2w=len.p2w-len.fr,
         ploss.p2w=loss.p2w/len.fr,
         sloss.p2w=case_when(
           loss.p2w<0 ~ "lost",
           loss.p2w>0 ~ "gained",
           TRUE ~ "same"
         ),
         sloss.p2w=factor(sloss.p2w,levels=c("lost","same","gained")))
#  CHECK
str(mysisL)
dim(mysisL)
levels(mysisL$tx)
levels(mysisL$stage)
levels(mysisL$group)

#  o Remove juveniles - decided this from len-wt work (i.e., precision
#    of wt measurements not adequate)
#  o need to re-factor to remove juvenile levels
mysisL <- mysisL |>
  filter(stage!="Juv") |>
  mutate(stage=factor(stage,levels=stage.lvls1),
         group=tx:stage)
#  CHECK ... the rows should match 65 in ms results
dim(mysisL)
levels(mysisL$tx)
levels(mysisL$stage)
levels(mysisL$group)


# ----- Results for Table 1
mysisL |>
  group_by(group) |>
  summarize(n=n(),
            minSL=min(len.fr),
            maxSL=max(len.fr),
            meanSL=mean(len.fr),
            sdSL=sd(len.fr),
            minChng=min(loss.p2w),
            maxChng=max(loss.p2w),
            meanChng=mean(loss.p2w),
            sdChng=sd(loss.p2w)) |>
  as.data.frame() |>
  gt() |>
  fmt_number(columns=c(3:5,7,8),decimals=1) |>
  fmt_number(columns=c(6,9,10),decimals=2)


# ----- Loss category summary ... First sentence in results
( tmp <- xtabs(~sloss.p2w,data=mysisL) )
round(tmp/sum(tmp)*100,1)


# ----- Loss summary and relationships ... Second sentence in results
t.test(loss.p2w~1,data=mysisL)
lenloss.lm <- lm(loss.p2w~len.fr*stage*tx,data=mysisL)
car::Anova(lenloss.lm)





# ===== Length-Weight Analysis I =====
# ----- Data.frame for length-weight analysis
#  o include only individuals to be used for length-weight analysis
#  o exclude values that were problematic for some reason (see notes)
#  o remove variables not used further (for conciseness)
#  o refactor tx, just to be sure
#  o refactor stage, just to be sure
#  o create group that is combo of the stage and treatment
#  o convert g to mg
#  o transform wt and len to log10
mysisLW <- mysis |>
  filter(include.LW) |>
  filter(use.lw) |>
  select(-(include.L:use.l)) |>
  mutate(tx=factor(tx),
         stage=factor(stage),
         group=tx:stage,
         c.wt=c.wt*1000,
         logwt=log10(c.wt),
         loglen=log10(len))
#  o CHECK
dim(mysisLW)
levels(mysisLW$tx)
levels(mysisLW$stage)
levels(mysisLW$group)

#  o remove juveniles as weights seem inaccurate
#  o need to re-factor to remove juvenile levels
mysisLW2 <- mysisLW |>
  filter(stage!="Juv") |>
  mutate(tx=factor(tx,levels=tx.lvls2),
         stage=factor(stage,levels=stage.lvls2),
         group=tx:stage)
#  o CHECK ... the rows should match 316 in ms results
dim(mysisLW2)
levels(mysisLW2$tx)
levels(mysisLW2$stage)
levels(mysisLW2$group)


# ----- Results for Table 2
mysisLW2 |>
  group_by(group) |>
  summarize(n=n(),
            minSL=min(len),
            maxSL=max(len),
            meanSL=mean(len),
            sdSL=sd(len),
            minDW=min(c.wt),
            maxDW=max(c.wt),
            meanDW=mean(c.wt),
            sdDW=sd(c.wt)) |>
  as.data.frame() |>
  gt() |>
  fmt_number(columns=c(3:5,7:9),decimals=1) |>
  fmt_number(columns=c(6,10),decimals=2)


# ----- First LW model analysis ... excludes gravid females b/c none in 8BF
#  o Remove gravid individuals
#  o Refactor stage (to remove grav)
mysisLW2_nograv <- mysisLW2 |>
  filter(stage!="Grav") |>
  mutate(stage=factor(stage))

#  o Plot to get a general visual ... not in manuscript
#    x Looks like no stage effect, slopes are likely parallel, at least 8SBF
#      intercept is greater than the other two tx
ggplot(data=mysisLW2_nograv,
       mapping=aes(y=logwt,x=loglen,color=tx,shape=stage)) +
  geom_smooth(method="lm",alpha=0.2) +
  geom_point()

#  o Fit the model & extract ANOVA results
#    x Three-way interaction not significant
#    x Two-way interactions with stage not significant
#    x stage main effect not significant
LW.lm1 <- lm(logwt~loglen*tx*stage,data=mysisLW2_nograv)
car::Anova(LW.lm1)

#  o Remove non-significant interaction terms, fit model again, check ANOVA,
#    compare slopes
#    ! emtrends() is a more robust way to compare slopes then compSlopes(), see
#      here ... https://fishr-core-team.github.io/fishR/blog/posts/2021-5-11_compSlopes-replacement/
#    ! compSlopes() has since been removed from the FSA package (I used it from
#      FSAmisc package below for comparison ... will have to install from GitHub
#      to run on your machine ... https://github.com/droglenc/FSAmisc).
#    ! emtrends() and compSlopes() provide different specific results but
#      qualitatively the same conclusion
#    ! unfortunately this conclusion seems different then what is noted in the
#      manuscript, as indicated by Bianca's question in the results
#    x This is one of those unfortunate statistical results ... the ANOVA
#      indicates a difference in slopes, but the post hoc test does not find a
#      difference due to the multiple comparison correction
LW.lm1a <- lm(logwt~loglen*tx,data=mysisLW2_nograv)
car::Anova(LW.lm1a)
LW.et1a <- emtrends(LW.lm1a,specs=pairwise~tx,var="loglen")
LW.ets1a <- summary(LW.et1a,infer=TRUE)
LW.ets1a$contrasts
LW.ets1a$emtrends
# FSAmisc::compSlopes(LW.lm1a)  # matches original analysis results

#  o Remove interaction (assume parallel lines), fit model, compare intercepts
#    ! emmeans() is more robust then compIntercepts() ... see here ...
#      https://fishr-core-team.github.io/fishR/blog/posts/2021-5-12_compIntercepts-replacement/
#    ! This will be different than what is in the manuscript, because the
#      conclusion about slopes was different
#    x Intercepts (so mean logDW at all logSL b/c parallel) differ
LW.lm1b <- lm(logwt~loglen+tx,data=mysisLW2_nograv)
LW.em1b <- emmeans(LW.lm1b,specs=pairwise~tx)
LW.ems1b <- summary(LW.em1b,infer=TRUE)
LW.ems1b$contrasts


# ----- Second LW model analysis ... excludes 8BF b/c no gravid females
#  o Remove 8BF tx individuals
#  o Refactor tx (to remove 8SBF)
mysisLW2_no8BF <- mysisLW2 |>
  filter(tx!="8BF") |>
  mutate(tx=factor(tx))

#  o Plot to get a general visual ... not in manuscript
#    x Looks like 8SBF may have steeper slope, no slope difference among stages
#      within tx, 8SBF generally heavier in the range of data, Gravids may be
#      heavier (esp for Fresh)
ggplot(data=mysisLW2_no8BF,
       mapping=aes(y=logwt,x=loglen,color=group,shape=group)) +
  geom_smooth(method="lm",alpha=0.2,se=FALSE) +
  geom_point()

#  o Fit the model & extract ANOVA results
#    x Three-way interaction not significant
#    x Two-way interactions with stage not significant
#    x Two-way interaction with tx significant ... diff slope by tx
LW.lm2 <- lm(logwt~loglen*tx*stage,data=mysisLW2_no8BF)
car::Anova(LW.lm2)

#  o Remove non-significant interaction terms, fit model again, check ANOVA,
#    compare slopes
#    x Sig difference in slope between Fresh and 8SBF tx, regardless of stage
#    # 8SBF tx is steeper
LW.lm2a <- lm(logwt~loglen*tx,data=mysisLW2_no8BF)
car::Anova(LW.lm2a)
LW.et2a <- emtrends(LW.lm2a,specs=pairwise~tx,var="loglen")
LW.ets2a <- summary(LW.et2a,infer=TRUE)
LW.ets2a$contrasts
LW.ets2a$emtrends
# FSAmisc::compSlopes(LW.lm2a)  # matches original analysis results

#  o Compare intercepts within 8SBF tx (assume same slope from above)
#    x No diff in intercepts ... same mean log DW for all log SL
mysisLW2_just8SBF <- mysisLW2 |>
  filter(tx=="8SBF") |>
  mutate(tx=factor(tx))
LW.lm2a1 <- lm(logwt~loglen+stage,data=mysisLW2_just8SBF)
car::Anova(LW.lm2a1)

#  o Compare intercepts within Fresh tx (assume same slope from above)
#    x Diff in intercepts ... Gravids greater that both Male and Female
mysisLW2_justFresh <- mysisLW2 |>
  filter(tx=="Fresh") |>
  mutate(tx=factor(tx))
LW.lm2a2 <- lm(logwt~loglen+stage,data=mysisLW2_justFresh)
car::Anova(LW.lm2a2)
LW.em2a2 <- emmeans(LW.lm2a2,specs=pairwise~stage)
LW.ems2a2 <- summary(LW.em2a2,infer=TRUE)
LW.ems2a2$contrasts




# ===== Length-Weight Analysis I =====
#  ! use mysisLW2 from above ... do not separate data.frames to deal with
#    no gravids in 8SBF ... instead see if lm() etc. can handle that

#  o fit model, assess ANOVA
#  x no three-way interaction or two-way with stage are significant; thus, no
#    slope difference due to stage (among other things)
#  x sig two-way of loglen and tx; thus, slope diff due to treatment
LW.lmX <- lm(logwt~loglen*tx*stage,data=mysisLW2)
car::Anova(LW.lmX)

#  o re-fit model with nonsignificant terms removed, recheck ANOVA
#  x still sig loglen and tx interaction, but barely
LW.lmXa <- lm(logwt~loglen+tx+stage+loglen:tx,data=mysisLW2)
car::Anova(LW.lmXa)

#  o compare slopes among treatments (averaged across stages)
#  x after correcting for multiple comparisons, it appears that no slopes differ
LW.etXa <- emtrends(LW.lmXa,specs=pairwise~tx,var="loglen")
LW.etsXa <- summary(LW.etXa,infer=TRUE)
LW.etsXa$contrasts

#  o fit model assuming no slopes differ, see if intercepts differ by tx or stage
#  x sig diff intercepts by tx and stage
LW.lmXb <- lm(logwt~loglen+tx+stage,data=mysisLW2)
car::Anova(LW.lmXb)

#  o see how tx intercepts differ
#  x intercepts differ among all three pairs of stages (all have unique ints)
LW.emXb <- emmeans(LW.lmXb,specs=pairwise~tx)
LW.emsXb <- summary(LW.emXb,infer=TRUE)
LW.emsXb$contrasts

#  o see how stage intercepts differ
#  x intercepts for Grav are greater than for male and female, which are same
LW.emXc <- emmeans(LW.lmXb,specs=pairwise~stage)
LW.emsXc <- summary(LW.emXc,infer=TRUE)
LW.emsXc$contrasts

#  X So, there would be 6 models males & females combined for all three tx, and
#    gravids for all three tx. However, we did not sample any gravids for 8BF,
#    so we would only have three models. All five models have the same slope
#    (are parallel) but different intercepts (once M&F are combined)

# ----- Final regression models
#   o combine males and females into one category in a new stage2 variable
#   o combine tx and stage2 into a new group2 variable
mysisLW2 <- mysisLW2 |>
  mutate(stage2=ifelse(stage=="Grav","Grav","MF"),
         stage2=factor(stage2,levels=c("MF","Grav")),
         group2=tx:stage2)

LW.final <- lm(logwt~loglen+group2,data=mysisLW2)
coef(LW.final)

ggplot(data=mysisLW2,
       mapping=aes(y=logwt,x=loglen,color=group2,shape=group2)) +
  geom_smooth(method="lm",alpha=0.2,se=FALSE) +
  geom_point()
