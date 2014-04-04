library(nlme)
library(lmerTest)
library(lme4)
library(bbmle)
library(MuMIn)

##--- Data input and manipulation
dados <- read.csv("data_diogo.csv", header = T, sep = ",")
dados=dados[,-15]
head(dados)
str(dados)


attach(dados)

trata=factor(div)
COMP<-factor(comp)
TIME<-factor(time)
table(trata:COMP)

mode<-lm(prod~trata+trata/COMP+TIME, data=dados)
mod1<-anova(mode)
mode2<-lm(prod~trata+PD+TIME, data=dados)
mod2<-anova(mode2)
pf(mod1[1,3]/mod1[3,3], 2, 18, lower.tail=F)#valor de P para efeito do trataemnto (diversidade)
names(mod2)

modeq<-lm(prod~trata+PD+I(PD^2)+TIME, data=dados)
modq<-anova(modeq)
summary(modeq)
pred.mod<-predict(modeq)
dados$pred<-pred.mod

mq<-lm(prod~PD+I(PD^2), data=dados)
coef(mq)


me<-aggregate(prod, list(trata, PD), "mean")
me$pred<-aggregate(dados$pred, list(trata, PD), "mean")$x
me$pred<-coef(mq)[1]+coef(mq)[2]*me$PD+coef(mq)[3]*me$PD^2

names(me)=c("Diversity", "PD", "Productivity", "pred")

ggplot(me,aes(x=PD, y=Productivity, colour=Diversity))+
	geom_line()

ggplot(me,aes(x=PD, y=Productivity, colour=Diversity))+
	geom_boxplot()

ggplot(me,aes(x=Diversity, y=Productivity))+
	geom_boxplot()

gf<-as.data.frame(cbind(rep(1:2, each=21), rbind(cbind(me$PD, me$Productivity), cbind(me$PD, me$pred))))
names(gf)=c("Group", "PD", "Productivity")

#points(me$PD,me$Productivity, type="l")
#plot(me$PD, me$pred.x, type="l")

ggplot(me, aes(x=PD))+
	geom_line(aes(y=pred), colour="blue")+
	geom_line(aes(y=Productivity), colour="red")
	

#summary(finalModel<-lmer(prod~nest + (1|tank)+(1|time) ,data=dados))
#results <- groupedData(prod~tank|time,outer = ~ trata,dados)#dando erro
#plot(ranef(finalModel))
#rand(finalModel)
#dotchart(ranef(finalModel))

#simulate(finalModel, 999)

#exactLRT 

#getME(finalModel)

##--- models for Productivity
ic<-function(dados){
	mymodel2 <- lmer(prod~PD* time  + (1|factor(tank))+(1|time),data=dados)
	mymodel3 <- lmer(prod ~ PD + time + (1|tank)+(1|time),data=dados)
	mymodel6<-lmer(prod ~ comp + (1|tank)+(1|time),data=dados)
	mymodel8<-lmer(prod ~ time + (1|tank)+(1|time),data=dados)
	mymodel9<-lmer(prod ~ 1 + (1|tank)+(1|time),data=dados)#"null" model only with intercept
	modelos<-list(mymodel2=mymodel2, mymodel6=mymodel6, mymodel8=mymodel8, mymodel3=mymodel3)
	tab1<-do.call(rbind, lapply(modelos, r.squaredGLMM))#Calculating conditional pseudo R2 for each model
	#Model selection table for the models
	tab2<-ICtab(mymodel2, mymodel3, mymodel6, mymodel8,mymodel9, 
							type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE,nobs=84)
	return(list(R2=tab1, AICTab=tab2))
}

lapply(dados1,ic)

##--- models for Decomposition
ic1<-function(dados){
	my2 <- lmer(dc~PD* time  + (1|tank)+(1|time),data=dados)
	my3 <- lmer(dc ~ PD + time + (1|tank)+(1|time),data=dados)
	my6<-lmer(dc ~ comp + (1|tank)+(1|time),data=dados)
	my8<-lmer(dc ~ time + (1|tank)+(1|time),data=dados)
	my9<-lmer(dc ~ 1 + (1|tank)+(1|time),data=dados)#"null" model only with intercept
	modelos<-list(my2=my2, my6=my6, my8=my8, my3=my3)
	tab1<-do.call(rbind, lapply(modelos, r.squaredGLMM))#Calculating conditional pseudo R2 for each model
	#Model selection table for the models
	tab2<-ICtab(my2, my3, my6, my8,my9, 
							type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE,nobs=84)
	return(list(R2=tab1, AICTab=tab2))
}

lapply(dados1,ic1)

##--- models for Respiration
ic2<-function(dados){
	my2 <- lmer(resp~PD* time  + (1|tank)+(1|time),data=dados)
	my3 <- lmer(resp ~ PD + time + (1|tank)+(1|time),data=dados)#--- ver se o tanque entra como fator no comp aleatório 
	my6<-lmer(resp ~ comp + (1|tank)+(1|time),data=dados)#verificar como entrar com o tempo como medida repetida no comp aleatório
	my8<-lmer(resp ~ time + (1|tank)+(1|time),data=dados)#verificar como interpretar interação no caso de medidas repetidas (ex. grafico da prod vs. tempo)
	my9<-lmer(resp ~ 1 + (1|tank)+(1|time),data=dados)#"null" model only with intercept
	modelos<-list(my2=my2, my6=my6, my8=my8, my3=my3)
	tab1<-do.call(rbind, lapply(modelos, r.squaredGLMM))#Calculating conditional pseudo R2 for each model
	#Model selection table for the models
	tab2<-ICtab(my2, my3, my6, my8,my9, 
							type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE,nobs=84)
	return(list(R2=tab1, AICTab=tab2))
}

#incluir interacao de riqueza com PD e outras variaveis de div filog
lapply(dados1,ic2)
