library(lme4)
library(lmerTest)
library(bbmle)
library(MuMIn)

##--- Data input and manipulation
dados <- read.csv("data_diogo.csv", header = T, sep = ";")
dados=dados[,-15]
head(dados)
str(dados)
attach(dados)

dados$PD2<-dados$PD/1000

dados1<-split(dados, dados$div)
str(dados1)
head(dados1$"9")
trata=factor(div)
filo<-factor(PD)
nest<-trata:filo
#summary(lmer(prod~nest + (~tank|time) ,data=dados))
summary(finalModel<-lmer(prod~nest + (1|tank)+(1|time) ,data=dados))
results <- groupedData(prod~tank|time,outer = ~ trata,dados)#dando erro
qqmath(ranef(finalModel))
rand(finalModel)

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
