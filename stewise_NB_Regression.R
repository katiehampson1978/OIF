library(car)
library(MASS)
library(ggplot2)
#
# negative binomial regression model for the analysis of the influence of ER, ST, TEA, AI and NRC [explanation below]
# on the number campaigns conducted 
#####################################
# data will change !
#####################################
#
# get the data
# Turkey already excluded
#
data <- read.table(file="/Users/tselhorst/Desktop/Table4paper_cc_woTurkey.txt", header=TRUE, sep="\t")
attach(data)
data
#
# description of the data
#
tapply(Y, ER, summary)
tapply(ST, ER, summary)
tapply(TEA, ER, summary)
tapply(AI, ER, summary)
tapply(NRC, ER, summary)
#
# Correlaton matrix plot
#
library(lattice)
scatterplotMatrix(~Y+log(ST)+TEA+AI+log(NRC) | ER, reg.line=FALSE, smoother=FALSE)
#
# This is the answer to question (8) of reviewer 2
#  He stated that TEA is calculated from SVA and ST
#  SVA is now excluded from the model.
#
# Variables ST and NRC are log transformed!!
#
# Correlations between influential variables
# AI is NOT significantly correlated to log.ST, TEA. 
# There is a tendency that AI is negatively correlated to log.NRC
# 
#
library(Hmisc)
d <- data.frame(log(ST),TEA,AI,log(NRC))
df <- data.matrix(d)
rcorr(df, type="spearman")
#
# Perform main effects nb regression model
#
# Parameters of the model are:
#   Y   : number of OIF campaigns
#   ER  : eradication status (finished | not finished)
#   ST  : size of territory (km^2)
#   TEA : territory ever affected (km^2)
#   AI  : area index
#   NRC : number of rabies cases at the beginning of ORV
#
# negative binomial model with and without log transformed data
#  log transformed data produce a better fit
#
summary(nbFit <- glm.nb(Y ~ (ER + ST + TEA + AI + NRC), data = data))
summary(nbFit.log <- glm.nb(Y ~ (ER + log(ST) + TEA + AI + log(NRC)), data = data))
#
# analyse the effect of ER
#
nbFit_wo_ER.log <- update(nbFit.log, . ~ . -ER)
anova(nbFit.log, nbFit_wo_ER.log)
#
# why not use the poisson model?
#
# setup a poission regression model with is nested in the negative binomial model and assuming that the disperion parameter is
# constant.
#
summary(poiFit <- glm(Y ~ (ER + ST + TEA + AI + NRC), family = poisson, data = data))
summary(poiFit.log <- glm(Y ~ (ER + log(ST) + TEA + AI + log(NRC)), family = poisson, data = data))
##############################
# Compare all AIC copied from the results screen
# negBin      AIC: 200.55  
# negBin.log  AIC: 197.44 <===
# poi         AIC: 246.15
# poi.log     AIC: 230.04
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
# Conclusion: ER (eradication status [finsihed / not finished]) is out!
# AIC: 196.28
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
require(ggplot2)
newdataAI <- data.frame(AI = rep(seq(from=min(data$AI), to = max(data$AI),length.out=100),2), ER=factor(rep(1:2, each=100), levels=1:2, labels = levels(data$ER)), TEA = mean(data$TEA))
newdataAI <- cbind(newdataAI, predict(stepNb, newdataAI, type="response", se.fit=TRUE))
newdataAI <- within(newdataAI, {
	LL <- fit-1.96*se.fit
	UL <- fit+1.96*se.fit
})
ggplot(newdataAI, aes(AI, fit))+geom_ribbon(aes(ymin=LL, ymax=UL, fill=ER),alpha=0.2)+geom_line(aes(color=ER),size=1) + labs(x="Area Index", y= "Predicted Number of OIF campaigns")
#
# Grafiken
# Einfluss des Territory ever affected
#
newdataTEA <- data.frame(TEA = rep(seq(from=min(data$TEA), to = max(data$TEA),length.out=100),2), ER=factor(rep(1:2, each=100), levels=1:2, labels = levels(data$ER)), AI = mean(data$AI))
newdataTEA <- cbind(newdataTEA, predict(stepNb, newdataTEA, type="response", se.fit=TRUE))
newdataTEA <- within(newdataTEA, {
	LL <- fit-1.96*se.fit
	UL <- fit+1.96*se.fit
})
ggplot(newdataTEA, aes(TEA, fit))+geom_ribbon(aes(ymin=LL, ymax=UL, fill=ER),alpha=0.2)+geom_line(aes(color=ER),size=1) + labs(x="Territory ever affected (proportion)", y= "Predichted Number of OIF campaigns")


