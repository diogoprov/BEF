#Analysis of Amy Downing's experimental data to incorporate PD as a predictor of 
#ecosystem functioning. Analysis based on linear mixed effect models.
#PhD chapter of Diogo Provete. Last modified: JAN 2014.
#(c) Diogo B. Provete | dbprovete@gmail.com

library(lme4)
library(bbmle)
library(MuMIn)
## Data input
###########

dados <- read.csv("data_diogo.csv", header = T, sep = ";")
dados=dados[,-15]
head(dados)
attach(dados)

dados$PD2<-dados$PD/1000
source("p.values.lmer.R")

#---Exploratory Plots
plot(prod~as.factor(PD2), xlab="Faith's PD", ylab="Productivity")
plot(prod~as.factor(div), xlab="Species Richness (Treatments)", ylab="Productivity")
boxplot(prod~as.factor(div)+time)
plot(prod~PD2, xlab="Faith's PD", ylab="Productivity")
boxplot(prod~as.factor(time), xlab="Time", ylab="Productivity")

## models
###########
#--Productivity
##Establishing the best random effects
model1 <- lmer(prod~PD2 + time + factor(div) + (1+time+PD2|tank), verbose=T, data=dados)
model2 <- lmer(prod~PD2 + time + factor(div) + (1+time|tank),data=dados)
model3 <- lmer(prod~PD2 + time + factor(div) + (1|tank),data=dados)#best model
model31 <- lmer(prod~PD2 + time + (1|tank)+(1|time),data=dados)#best model
model32 <- lmer(prod~PD2 + time + factor(div) + (1|tank/time),data=dados)#best model
model4 <- lmer(prod~PD2 + time + factor(div) + (1+PD2|tank), data=dados)
anova(model1, model2, model3, model4)
(ICtab(model1, model2, model3, model4,model31,
									 type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE, nobs=504))

##Establishing the best fixed effects model
ic<-function(dados){
mymodel2 <- lmer(prod~PD* time  + (1|tank)+(1|time),data=dados)
mymodel3 <- lmer(prod ~ PD + time + (1|tank)+(1|time),data=dados)
#mymodel04<-lmer(prod ~ PD2*div + (1|tank),data=dados)#best model
#mymodel24<-lmer(prod ~ PD2*factor(div) + time+ (1|tank),data=dados)
#mymodel4<-lmer(prod ~ PD2 + factor(div) + (1|tank),data=dados)
#mymodel5 <- lmer(prod ~ factor(div) + time + (1|tank),data=dados)
mymodel6<-lmer(prod ~ comp + (1|tank)+(1|time),data=dados)
#mymodel7 <- lmer(prod ~ factor(div) + (1|tank),data=dados)
mymodel8<-lmer(prod ~ time + (1|tank)+(1|time),data=dados)
mymodel9<-lmer(prod ~ 1 + (1|tank)+(1|time),data=dados)#"null" model only with intercept

modelos<-list(mymodel2=mymodel2, mymodel6=mymodel6, mymodel8=mymodel8, mymodel3=mymodel3)
tab1<-do.call(rbind, lapply(modelos, r.squaredGLMM))
#Model selection table for the models
tab2<-ICtab(mymodel2, mymodel3, mymodel6, mymodel8,mymodel9, 
									 type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE,nobs=84)
return(list(R2=tab1, AICTab=tab2))
}

lapply(dados1,ic)

#Calculating conditional pseudo R2 for each model


#r.squaredGLMM(mymodel04)
#r.squaredGLMM(mymodel4)
r.squaredGLMM(mymodel2)
#r.squaredGLMM(mymodel7)
r.squaredGLMM(mymodel6)
#r.squaredGLMM(mymodel5)
r.squaredGLMM(mymodel8)
r.squaredGLMM(mymodel3)

#Calculating the % Deviance Explained (~ R2) by each model
100*(deviance(mymodel9) - deviance(mymodel2))/deviance(mymodel9)
100*(deviance(mymodel9) - deviance(mymodel3))/deviance(mymodel9)
100*(deviance(mymodel9) - deviance(mymodel4))/deviance(mymodel9)
100*(deviance(mymodel9) - deviance(mymodel5))/deviance(mymodel9)
100*(deviance(mymodel9) - deviance(mymodel6))/deviance(mymodel9)
100*(deviance(mymodel9) - deviance(mymodel7))/deviance(mymodel9)
100*(deviance(mymodel9) - deviance(mymodel8))/deviance(mymodel9)

## FIGURE
############
boxplot(prod ~ factor(time), ylab = "Productivity") ## productivity decreases with time regardless of richness

#--Models for Decomposition
##Establishing the best random effects
m1 <- lmer(dc~PD2 + time + factor(div) + (1+time+PD2|tank), data=dados)
m2 <- lmer(dc~PD2 + time + factor(div) + (1+time|tank),data=dados)
m3 <- lmer(dc~PD2 + time + factor(div) + (1|tank),data=dados)#best model
m4 <- lmer(dc~PD2 + time + factor(div) + (1+PD2|tank), data=dados)
anova(m1, m2, m3, m4)
(ICtab(m1, m2, m3, m4,  type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE, nobs=504))

##Establishing the best fixed effects model
my2 <- lmer(dc~PD2 + time + factor(div) + (1|tank),data=dados)
my3 <- lmer(dc ~ PD2 + time + (1|tank),data=dados)
my4<-lmer(dc ~ PD2 + factor(div) + (1|tank),data=dados)
my40<-lmer(dc ~ PD2*factor(div) + (1|tank),data=dados)
my5 <- lmer(dc ~ factor(div) + time + (1|tank),data=dados)
my6<-lmer(dc ~ PD2 + (1|tank),data=dados)
my7 <- lmer(dc ~ factor(div) + (1|tank),data=dados)
my8<-lmer(dc ~ time + (1|tank),data=dados)#best model
my9<-lmer(dc ~ 1 + (1|tank),data=dados)#"null" model only with intercept

#Model selection table for the models
(ICtab(my2, my3,my4, my40,my5, my6, my7,my8,my9, 
									 type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE,nobs=504))
r.squaredGLMM(my8)[[2]]
r.squaredGLMM(my3)[[2]]
r.squaredGLMM(my5)[[2]]
r.squaredGLMM(my2)[[2]]
r.squaredGLMM(my6)[[2]]
r.squaredGLMM(my7)[[2]]
r.squaredGLMM(my4)[[2]]
r.squaredGLMM(my40)[[2]]

#--Models for Respiration
##Establishing the best random effects
mo1 <- lmer(resp~PD2 + time + factor(div) + (1+time+PD2|tank), data=dados)
mo2 <- lmer(resp~PD2 + time + factor(div) + (1+time|tank),data=dados)
mo3 <- lmer(resp~PD2 + time + factor(div) + (1|tank),data=dados)#best model
mo4 <- lmer(resp~PD2 + time + factor(div) + (1+PD2|tank), data=dados)

(ICtab(mo1, mo2, mo3, mo4, type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE, nobs=504))

##Establishing the best fixed effects model
mm2 <- lmer(resp~PD2 + time + factor(div) + (1|tank),data=dados)
mm3 <- lmer(resp ~ PD2 + time + (1|tank),data=dados)
mm4<-lmer(resp ~ PD2 + factor(div) + (1|tank),data=dados)
mm14<-lmer(resp ~ PD2 * factor(div) + (1|tank),data=dados)
mm5 <- lmer(resp ~ factor(div) + time + (1|tank),data=dados)
mm6<-lmer(resp ~ PD2 + (1|tank),data=dados)
mm7 <- lmer(resp ~ factor(div) + (1|tank),data=dados)
mm8<-lmer(resp ~ time + (1|tank),data=dados)#best model
mm9<-lmer(resp ~ 1 + (1|tank),data=dados)#"null" model only with intercept

#Model selection table for the models
(ICtab(mm2, mm3,mm4, mm14,mm5, mm6, mm7,mm8,mm9, 
			 type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE,nobs=504))

#Calculating conditional pseudo R2 for each model
r.squaredGLMM(mm8)[[2]]
r.squaredGLMM(mm3)[[2]]
r.squaredGLMM(mm5)[[2]]
r.squaredGLMM(mm2)[[2]]
r.squaredGLMM(mm6)[[2]]
r.squaredGLMM(mm7)[[2]]
r.squaredGLMM(mm4)[[2]]
r.squaredGLMM(mm14)[[2]]