library(car)
library(MASS)
library(ggplot2)
#
# Statistical analysis of the influence of ER, ST, TEA, AI and NRC [explanation below] on the number of
# campaigns
#
# List of Acronyms
#
# Country:____________________________________________________ Country
# Year ORV began:_____________________________________________ OIFS
# Cases at start of ORV:______________________________________ NRC
# Cases in 2010:______________________________________________ C_2010
# Last rabies case (except bats, humans and imported cases)___ LAST_CASE_Y
# Reduction of rabies:________________________________________ R_R
# ORV campaigns to eliminate rabies:__________________________ ORV_2_C100
# ORV campaigns to control* rabies (90%):_____________________ ORV_2_C90
# Size of territory (km2):____________________________________ ST
# Border length with endemic areas (km):______________________ BL
# Vaccinated area (km2):______________________________________ SVA
# Proportion of territory vaccinated:_________________________ TEA
# Mean Area Index:____________________________________________ AI
# campaigns until 2010:_______________________________________ Y <== dependent variable
# ER (ORV finished / not finished):___________________________ ER
#
data <- read.table(file="/Users/tselhorst/Desktop/OIF/data.txt", header=TRUE, sep="\t")
attach(data)
data
#
# Summary description of data
#
tapply(Y, ER, summary)
tapply(ST, ER, summary)
tapply(TEA, ER, summary)
tapply(AI, ER, summary)
tapply(NRC, ER, summary)
tapply(BL, ER, summary)
tapply(SVA, ER, summary)
#
# Correlatons and correlation matrix plot
#
library(Hmisc)
d <- data.frame(AI, TEA, BL, log(SVA), log(ST), log(NRC))
df <- data.matrix(d)
rcorr(df, type="spearman")
#
# corr matrix plot
#
library(lattice)
# log(ST) is correlated with TEA
scatterplotMatrix(~Y+AI+TEA+BL+log(SVA)+log(ST)+log(NRC) | ER, reg.line=FALSE, smoother=FALSE)
#
# Perform stepwise main effects nb regression model
#
# Parameters of the model are:
#   Y   : number of OIF campaigns
#   ER  : eradication status (finished | not finished)
#   TEA : territory ever affected (proportion)
#   AI  : area index
#   SVA : size of vaccination areq <== element within the calc of AI
#
#
model.nb <- glm.nb(Y ~ (log(SVA) + ER + TEA + AI), data = data)
model.nb.step <- stepAIC(model.nb, direction="backward")
model.nb.step$anova
summary(model.nb.step)
with(model.nb.step, cbind(res.deviance = deviance, df=df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))
#
# include second oder interactions to see whether the better fit better
model.nb.step.2 <- stepAIC(model.nb.step,~.^2, direction="both")
model.nb.step.2$anova
summary(model.nb.step.2)
with(model.nb.step.2, cbind(res.deviance = deviance, df=df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))
#
# why not use the poisson model?
#
# setup a poission regression model with is nested in the negative binomial model and assuming that the disperion parameter is
# constant.
#
model.poi <- glm(Y ~ (log(SVA) + ER + TEA + AI), family = poisson, data = data)
summary(model.poi)
model.poi.step = stepAIC(model.poi, direction="backward")
model.poi.step$anova
summary(model.poi.step)
#
# include second order interactions to see whether the model performs better
model.poi.step.2 <- stepAIC(model.poi.step,~.^2, direction="both")
model.poi.step.2$anova
summary(model.poi.step.2)
with(model.poi.step.2, cbind(res.deviance = deviance, df=df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))
#
# compare negbin with poi
#
chi <- 2 * (logLik(model.nb.step.2) - logLik(model.poi.step.2))
chi
pchisq(chi, df = 1, lower.tail = FALSE)
#
##############################
# Compare all AIC copied from the results screen 
# negBin      AIC: 132.78  
# negBin.log  AIC: 127.16 <===
# poi         AIC: 146.26
# poi.log     AIC: 130.71
##############################
#
# Compare the negative binomial to the poission model
#  This is to answer question 8 of reviewer 2 
#  The negbin is superior to the poi model.
#  The difference is significant!
#
chi <- 2 * (logLik(nbFit.log) - logLik(poiFit.log))
chi
pchisq(chi, df = 1, lower.tail = FALSE)
#
#
# Goodness of fit of negbin
# Conclusion: The model fits the data
#
with(nbFit.log, cbind(res.deviance = deviance, df=df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))
#
# include and test 2nd order interactions for nbFit model 
#
stepNb <- stepAIC(nbFit.log,~.^2, direction="both")
stepNb$anova
summary(stepNb)
#
# Now we have a model with lost of interactions and lots of dicussions :-)
# AIC: 115.14
# 
#
# goodness of fit
# the model fits the data
#
with(stepNb, cbind(res.deviance = deviance, df=df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))
#
# Confidence intervals for the coefficients of the model
#
est <- cbind(Estimate=coef(stepNb), confint(stepNb))
est
#
# Predicted counts for ER
#
predictedER <- data.frame(AI=mean(AI), TEA = mean(TEA), ER=factor(1:2, level=1:2, labels=levels(ER)))
predictedER <- cbind(predictedER, predict(stepNb, predictedER, type="response", se.fit=TRUE))
predictedERconfint <- within(predictedER, {
  LL <- fit-1.96*se.fit
  UL <- fit+1.96*se.fit
})
predictedERconfint
#
#######################
# Grafics have to be reshaped because ER is now out!
#  BUT I WAIT FOR THE NEW DATA !
# code below will not run
#######################
#
# Grafiken
# Einfluss des Area Index
#
#require(ggplot2)
#newdataAI <- data.frame(AI = rep(seq(from=min(data$AI), to = max(data$AI),length.out=100),2), ER=factor(rep(1:2, each=100), levels=1:2, labels = levels(data$ER)), TEA = mean(data$TEA))
#newdataAI <- cbind(newdataAI, predict(stepNb, newdataAI, type="response", se.fit=TRUE))
#newdataAI <- within(newdataAI, {
#	LL <- fit-1.96*se.fit
#	UL <- fit+1.96*se.fit
#})
#ggplot(newdataAI, aes(AI, fit))+geom_ribbon(aes(ymin=LL, ymax=UL, fill=ER),alpha=0.2)+geom_line(aes(color=ER),size=1) + labs(x="Area Index", y= "Predicted Number of OIF campaigns")
#
# Grafiken
# Einfluss des Territory ever affected
#
#newdataTEA <- data.frame(TEA = rep(seq(from=min(data$TEA), to = max(data$TEA),length.out=100),2), ER=factor(rep(1:2, each=100), levels=1:2, labels = levels(data$ER)), AI = mean(data$AI))
#newdataTEA <- cbind(newdataTEA, predict(stepNb, newdataTEA, type="response", se.fit=TRUE))
#newdataTEA <- within(newdataTEA, {
#	LL <- fit-1.96*se.fit
#	UL <- fit+1.96*se.fit
#})
#ggplot(newdataTEA, aes(TEA, fit))+geom_ribbon(aes(ymin=LL, ymax=UL, fill=ER),alpha=0.2)+geom_line(aes(color=ER),size=1) + labs(x="Territory ever affected (proportion)", y= "Predichted Number of OIF campaigns")


