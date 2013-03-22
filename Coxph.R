#
# Cox proportional hazard analysis of the OIF data
#
library(survival)
library(KMsurv)
library(MASS)
library(car) # for logit transformation
#
# read the data
data <- read.table(file="/Users/tselhorst/Desktop/OIF/data.txt", header=TRUE, sep="\t")
#Include 'years of ORV experience' and 'rabies incidence - cases per km2'
data$YE <- data$OIFS-1978 
data$RI <- data$NRC/data$SVA 
data
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
# Years of ORV experience:____________________________________ YE
# Rabies incidence (cases per km2):___________________________ RI

#
#Examine correlations between variables 
#Correlation matrix plot
library(lattice)
setwd("/Users/KHampson/collaborators/Eradication endgames/Phil_Trans_Drafts/Mueller")
pdf("CorrelationMatrix.pdf", 10, 10)
scatterplotMatrix(~Y+log(ST)+TEA+AI+log(NRC)+experience | ER, reg.line=FALSE, smoother=FALSE)
dev.off()

d <- data.frame(log(ST),TEA,AI,log(NRC),experience,BL,log(SVA), RI)
df <- data.matrix(d)
rcorr(df, type="spearman")
#Cannot include: log(ST), log(NRC), experience, BL 
#Can include: TEA, AI, log(SVA), RI
#

#
# reshape the data so that we have a column Censor conditioned on ER which 
# tells us wheter the OIF is finishes or not. 
# When OIF is finished Censor is TRUE which tells me that the event was 
# really observed
data4surv <- transform(data, Censor=ifelse(ER=="finished", TRUE, FALSE))
data4surv$experience <- data$OIFS-1978
data4surv
attach(data4surv)
#
# create surv object with censored data
#
so <- Surv(Y, Censor, type='right')
so
#Create new dependent variable (YC) to look at time to control
#censor for control
Censor2=ifelse(!is.na(ORV_2_C90), TRUE, FALSE)
YC=replace(ORV_2_C90, which(is.na(ORV_2_C90)), Y[which(is.na(ORV_2_C90))])
s2 <- Surv(YC, Censor2, type='right')
s2

#
# null model
null.model <- coxph(so~1)
summary(null.model)
# we alread know form negative binomial regression that onyl AI, TEA and log(SVA)
# are not correlated and can be used in the model. So let's start with 
# this one.
##Can include: TEA, AI, log(SVA), RI
s.model <- coxph(so~AI+TEA+log(SVA)+RI)
summary(s.model)
survfit(s.model)
#
# fit Coxph model in both directions. This time I want to restrict model
# parameters within interval [0,1]. This is achived with the logit transformation.
# interactions of order 2 are allowed
#
Scope = list(upper=~logit(AI)+logit(TEA)+log(SVA)+RI+.^2,lower=~1)
phm_e=coxph(so~1, method="efron")
phm_f=stepAIC(phm_e,Scope,direction="both")
summary(phm_f)
plot(survfit(phm_f))
survfit(phm_f)
## log odds increase of 1 in TEA results in >50% (57%) decrease in hazard (to eliminate infection) - increasing time to elimination 
## log odds increase of 1 in AI results in almost 3X (290%) increase in hazard -reducing time to elimination 
# model explains 57% variance (R2)
#
# plots for different figures of AI (min and max)
OIF.Ai <- data.frame(AI=c(min(AI),max(AI)), TEA=rep(median(TEA), 2), SVA=rep(median(SVA), 2), RI=rep(median(RI), 2))
plot(survfit(phm_f, newdata=OIF.Ai), conf.int=TRUE, lty=c(1,2), col=c("red", "blue"),
     xlab="number of Campaigns", ylab="p(eradication not finished)")
legend(locator(1), legend=c('AI=0.016', 'AI=1.00'), lty=c(1,2), col=c("red", "blue"))
summary(survfit(phm_f, newdata=OIF.Ai))
#
# plots for diferent figures of TEA (min and max)
OIF.Tea <- data.frame(TEA=c(min(TEA),max(TEA)),AI=rep(median(AI),2), SVA=rep(median(SVA), 2), RI=rep(median(RI), 2))
plot(survfit(phm_f, newdata=OIF.Tea), conf.int=TRUE, lty=c(1,2), col=c("red", "blue"),
     xlab="number of Campaigns", ylab="p(eradication not finished)")
legend(locator(1), legend=c('TEA=0.016', 'TEA=1.00'), lty=c(1,2), col=c("red", "blue"))
summary(survfit(phm_f, newdata=OIF.Tea))



###########################################################
#
# Control - null model 
null.model2 <- coxph(s2~1)
summary(null.model2)
s.model2 <- coxph(s2~AI+TEA+log(SVA)+RI)
summary(s.model2)
survfit(s.model2)
#        
# fit Coxph model in both directions. This time I want to restrict model
# parameters within interval [0,1]. This is achived with the logit transformation.
# interactions of order 2 are allowed
#
phm_e2=coxph(s2~1, method="efron")
phm_f2=stepAIC(phm_e2,Scope,direction="both")
summary(phm_f2)
plot(survfit(phm_f2))
survfit(phm_f2, show.rmean=TRUE)
## increase of 1 in RI results in (1-2.07E-18)% decrease in hazard (to eliminate infection) - increasing time to elimination 
#

# plots for different ranges of RI (min and max)
OIF.RI <- data.frame(TEA=rep(median(TEA),2),AI=rep(median(AI),2), SVA=rep(median(SVA), 2), RI=c(min(RI), max(RI)))
plot(survfit(phm_f2, newdata=OIF.RI), conf.int=TRUE, lty=c(1,2),
     xlab="number of Campaigns", ylab="p(not yet controlled)")
legend(locator(1), legend=c(paste('RI=', round(min(RI),3)), paste('RI=', round(max(RI),3))), lty=c(1,2))
#

#Compare elimination versus control
plot(survfit(phm_f), xlim=c(0,70), conf.int=TRUE)
lines(survfit(phm_f2), col="red", conf.int=TRUE)
