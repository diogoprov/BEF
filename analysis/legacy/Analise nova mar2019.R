# Mon Mar 25 11:25:45 2019 ------------------------------

#---Script novo com analises Thiago + Pedro

#---Carregando os pacotes
library(bbmle)
library(MuMIn)
library(tidyverse)
library(picante)
library(lme4)

#---Importando os dados
dados <- read.csv("data_diogo.csv", header = T, sep = ",")
head(dados)

#---Working wih dated phylogenies
amy.nature=read.tree("phylog_dated_2002.txt")
amy.nature=drop.tip(amy.nature, c("Vallisneria_americana", "Utricularia_vulgaris"))

#---Manipulando os dados
dados$div <- as.factor(dados$div)
dados$comp <- as.factor(dados$comp)
dados$time <- as.factor(dados$time)

#---Data about treatments, calculating PD
treat=read.table("treat.txt", h=T)

match.phylo.comm(amy.nature, t(treat))

#--Independent on species richness
mpd.nature=mpd(t(treat),cophenetic(amy.nature))#mpd
mntd.nature=mntd(t(treat),cophenetic(amy.nature))#mnnd

d<-dados[order(dados$comp),]#sorting data by composition
head(d)
d$mntd<-rep(mntd.nature,each=12)
d$mpd<-rep(mpd.nature,each=12)#repeat MPD 12 times

#---Explorando os dados


#---Modelo para decomposição
modeL<-lm(dc~trata+trata/COMP+TIME, data=d)
modee1<-anova(modeL)
modee2<-lm(dc~trata+PD+TIME, data=dados)
md2<-anova(modee2)

modeqdec<-lm(dc~trata+PD+I(PD^2)+TIME, data=dados)#model with PD as a quadratic term
anova(modeqdec)#the quadratic term is significant

#---Modelo para respiração
resp.mpd<-lm(resp~trata+mpd+TIME, data=d)
res.resp.lin<-anova(m.resp.lin)
qmpd<-lm(resp~trata+mpd+I(mpd^2)+TIME, data=d)#model with PD as a quadratic term
anova(mdeq)#the quadratic term is significant

#---Modelo para produtividade MPD
prod.pd<-lm(prod~trata+mpd+TIME, data=d)
mq.mpd<-lm(prod~trata+mpd+I(mpd^2)+TIME, data=d)

prod.mpd.misto<-lmer(prod~trata+trata/mpd+(1|TIME), data=d)

#---Multifuncionalidade

