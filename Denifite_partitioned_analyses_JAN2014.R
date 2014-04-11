

##--- Loading required packages

library(bbmle)
library(MuMIn)
library(ggplot2)

##--- Data input and manipulation
dados <- read.csv("data_diogo.csv", header = T, sep = ",")
head(dados)
attach(dados)

trata=factor(div)
COMP<-factor(comp)
TIME<-factor(time)

me<-aggregate(prod, list(trata, PD), "mean")

##--- Linear Model with nested design and repeated measures for Productivity
mode<-lm(prod~trata+trata/COMP+TIME, data=dados)
mod1<-anova(mode)
mode2<-lm(prod~trata+PD+TIME, data=dados)
mod2<-anova(mode2)
pf(mod1[1,3]/mod1[3,3], 2, 18, lower.tail=F)#valor de P para efeito do tratamento (diversidade)
names(mod2)

##-- Quadratic Model with nested design and repeated measures for Productivity
modeq<-lm(prod~trata+PD+I(PD^2)+TIME, data=dados)#model with PD as a quadratic term
anova(modeq)#the quadratic term is significant
summary(modeq)#to get the parameters and estimated values of the model
pred.mod<-predict(modeq)#predicted values for each observation
dados$pred<-pred.mod#binding to the data frame

me$pred<-aggregate(dados$pred, list(trata, PD), "mean")$x

mq<-lm(prod~PD+I(PD^2), data=dados)#to get the values without time, to make the plot 
#of the quadratic model
anova(mq)
summary(mq)
coef(mq)

##-- Plotting

#to make the plot
names(me)=c("Diversity", "PD", "Productivity", "pred")#renaming the columns to make the plot
me$pred<-coef(mq)[1]+coef(mq)[2]*me$PD+coef(mq)[3]*me$PD^2#adding the predicted values 
#for PD by the quadratic model

#width represents PD
ggplot(me,aes(x=PD, y=Productivity, colour=Diversity))+
	geom_line()

ggplot(me,aes(x=PD, y=Productivity, colour=Diversity))+
	geom_boxplot()

#same version but bw, instead
ggplot(me,aes(x=Diversity, y=Productivity))+
	geom_boxplot()

gf<-as.data.frame(cbind(rep(1:2, each=21), rbind(cbind(me$PD, me$Productivity), cbind(me$PD, me$pred))))
names(gf)=c("Group", "PD", "Productivity")

ggplot(gf,aes(x=PD, y=Productivity, group=Group))+
geom_line()

#Final plot showing the variation in PD in red and the quadratic model in blue
ggplot(me, aes(x=PD))+
	geom_line(aes(y=pred), colour="blue")+
	geom_line(aes(y=Productivity), colour="red")

#same version but bw, instead
gf<-as.data.frame(cbind(rep(1:2, each=21), rbind(cbind(me$PD, me$Productivity), cbind(me$PD, me$pred))))
names(gf)=c("Group", "PD", "Productivity")

ggplot(gf,aes(x=PD, y=Productivity, group=Group))+
	geom_point()

##--- Model for Decomposition

dec.plot<-aggregate(dc, list(trata, PD), "mean")
names(dec.plot)=c("Diversity", "PD", "Decomposition")#renaming the columns to make the plot
ggplot(dec.plot,aes(x=PD, y=Decomposition, colour=Diversity))+
	geom_line()

modeL<-lm(dc~trata+trata/COMP+TIME, data=dados)
modee1<-anova(modeL)
modee2<-lm(dc~trata+PD+TIME, data=dados)
md2<-anova(modee2)
pf(modee1[1,3]/modee1[3,3], 2, 18, lower.tail=F)#valor de P para efeito do tratamento (diversidade)
names(md2)

##--- models for Respiration
resp.plot<-aggregate(resp, list(trata, PD), "mean")
names(resp.plot)=c("Diversity", "PD", "Decomposition")#renaming the columns to make the plot
ggplot(resp.plot,aes(x=PD, y=Decomposition, colour=Diversity))+
	geom_line()

m.resp.lin<-lm(resp~trata+PD+TIME, data=dados)
res.resp.lin<-anova(m.resp.lin)
pf(res.resp.lin[1,3]/res.resp.lin[3,3], 2, 18, lower.tail=F)#valor de P para efeito do tratamento (diversidade)

mdeq<-lm(resp~trata+PD+I(PD^2)+TIME, data=dados)#model with PD as a quadratic term
anova(mdeq)#the quadratic term is significant

modq<-lm(resp~PD+I(PD^2), data=dados)#to get the values without time, to make the plot 
#of the quadratic model
anova(modq)
summary(modq)
coef(modq)
resp.pred<-predict(modq)#predicted values for each observation
dados$resp.pred<-resp.pred#binding to the data frame

resp.plot$resp.pred<-aggregate(dados$resp.pred, list(trata, PD), "mean")$x

resp.plot$resp.pred<-coef(modq)[1]+coef(modq)[2]*resp.plot$PD+coef(modq)[3]*resp.plot$PD^2#adding the predicted values 
#for PD by the quadratic model
names(resp.plot)=c("Diversity", "PD", "Respiration", "pred")#renaming the columns to make the plot

#Final plot showing the variation in PD in red and the quadratic model in blue
ggplot(resp.plot, aes(x=PD))+
	geom_line(aes(y=pred), colour="blue")+
	geom_line(aes(y=Respiration), colour="red")

#Just another version with data as points and model as line
ggplot(resp.plot, aes(x=PD)) + 
	geom_line(aes(y=pred), colour="blue", linetype="solid", size=1) + 
	geom_point(aes(y=Respiration), colour="red", size=4, shape=21, fill="white")

##--- Model selection for the other PD variables

#-- Models for Productivity

#-with PD
prod.PD<-lm(prod~trata+PD+TIME, data=dados)
mq.prodPD<-lm(prod~trata+PD+I(PD^2)+TIME, data=dados)

#-With MPD
prod.pd<-lm(prod~trata+mpd+TIME, data=d)
mq.mpd<-lm(prod~trata+mpd+I(mpd^2)+TIME, data=d)

#-With Haed
prod.haed<-lm(prod~trata+haed+TIME, data=d)
mq.haed<-lm(prod~trata+haed+I(haed^2)+TIME, data=d)

#-With Rao
prodrao<-lm(prod~trata+raodiv+TIME, data=d)
mq.rao<-lm(prod~trata+raodiv+I(raodiv^2)+TIME, data=d)

#-With MNTD
prod.mntd<-lm(prod~trata+mntd+TIME, data=d)
mq.mntd<-lm(prod~trata+mntd+I(mntd^2)+TIME, data=d)

#-With Simpson Phylogenetic
prod.SimpPhyl<-lm(prod~trata+SimpsonPhyl+TIME, data=d)
mq.SimpPhyl<-lm(prod~trata+SimpsonPhyl+I(SimpsonPhyl^2)+TIME, data=d)

#-With Shannon
prod.Shannon<-lm(prod~trata+Shannon+TIME, data=d)
mq.Shannon<-lm(prod~trata+Shannon+I(Shannon^2)+TIME, data=d)

#-With Simpson
prod.Simpson<-lm(prod~trata+Simpson+TIME, data=d)
mq.Simpson<-lm(prod~trata+Simpson+I(Simpson^2)+TIME, data=d)

#-With Simpson
prod.Intens_Q<-lm(prod~trata+Intens_Q+TIME, data=d)
mq.Intens_Q<-lm(prod~trata+Intens_Q+I(Intens_Q^2)+TIME, data=d)

#-- Model selection table			
(ICtab(prod.PD,mq.prodPD, prod.pd,mq.mpd,prod.haed,mq.haed,prodrao,mq.rao,prod.mntd,
			 mq.mntd,prod.SimpPhyl, mq.SimpPhyl, prod.Shannon, mq.Shannon, prod.Simpson,
			 mq.Simpson,prod.Intens_Q,mq.Intens_Q,
			 type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE, nobs=504))

#-- Models for Respiration

#-with PD
resp.PD<-lm(resp~trata+PD+TIME, data=dados)
qpd<-lm(resp~trata+PD+I(PD^2)+TIME, data=dados)

#-With MPD
resp.mpd<-lm(resp~trata+mpd+TIME, data=d)
qmpd<-lm(resp~trata+mpd+I(mpd^2)+TIME, data=d)

#-With Haed
resp.haed<-lm(resp~trata+haed+TIME, data=d)
qhaed<-lm(resp~trata+haed+I(haed^2)+TIME, data=d)

#-With Rao
resprao<-lm(resp~trata+raodiv+TIME, data=d)
qrao<-lm(resp~trata+raodiv+I(raodiv^2)+TIME, data=d)

#-With MNTD
resp.mntd<-lm(resp~trata+mntd+TIME, data=d)
qmntd<-lm(resp~trata+mntd+I(mntd^2)+TIME, data=d)

#-With Simpson Phylogenetic
resp.SimpPhyl<-lm(resp~trata+SimpsonPhyl+TIME, data=d)
qSimpPhyl<-lm(resp~trata+SimpsonPhyl+I(SimpsonPhyl^2)+TIME, data=d)

#-With Shannon
resp.Shannon<-lm(resp~trata+Shannon+TIME, data=d)
qShannon<-lm(resp~trata+Shannon+I(Shannon^2)+TIME, data=d)

#-With Simpson
resp.Simpson<-lm(resp~trata+Simpson+TIME, data=d)
qSimpson<-lm(resp~trata+Simpson+I(Simpson^2)+TIME, data=d)

#-With Simpson
resp.Intens_Q<-lm(resp~trata+Intens_Q+TIME, data=d)
qIntens_Q<-lm(resp~trata+Intens_Q+I(Intens_Q^2)+TIME, data=d)

#-- Model selection table			
(ICtab(resp.PD,qpd, qmntd,resp.mntd,qrao,resprao,qhaed,resp.haed,qmpd,resp.mpd,
			resp.SimpPhyl,qSimpPhyl, resp.Shannon, qShannon, resp.Intens_Q, qIntens_Q,
			 type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE, nobs=504))

#-- Models for Decomposition

#-with PD
dec.PD<-lm(dc~trata+PD+TIME, data=dados)
qdec<-lm(dc~trata+PD+I(PD^2)+TIME, data=dados)

#-With MPD
dec.mpd<-lm(dc~trata+mpd+TIME, data=d)
decqmpd<-lm(dc~trata+mpd+I(mpd^2)+TIME, data=d)

#-With Haed
dc.haed<-lm(dc~trata+haed+TIME, data=d)
qdhaed<-lm(dc~trata+haed+I(haed^2)+TIME, data=d)

#-With Rao
dcrao<-lm(dc~trata+raodiv+TIME, data=d)
qdrao<-lm(dc~trata+raodiv+I(raodiv^2)+TIME, data=d)

#-With MNTD
dc.mntd<-lm(dc~trata+mntd+TIME, data=d)
qdmntd<-lm(dc~trata+mntd+I(mntd^2)+TIME, data=d)

#-With Simpson Phylogenetic
dc.SimpPhyl<-lm(dc~trata+SimpsonPhyl+TIME, data=d)
qdSimpPhyl<-lm(dc~trata+SimpsonPhyl+I(SimpsonPhyl^2)+TIME, data=d)

#-With Shannon
dc.Shannon<-lm(dc~trata+Shannon+TIME, data=d)
qdShannon<-lm(dc~trata+Shannon+I(Shannon^2)+TIME, data=d)

#-With Simpson
dc.Simpson<-lm(dc~trata+Simpson+TIME, data=d)
qdSimpson<-lm(dc~trata+Simpson+I(Simpson^2)+TIME, data=d)

#-With Simpson
dc.Intens_Q<-lm(dc~trata+Intens_Q+TIME, data=d)
qdIntens_Q<-lm(dc~trata+Intens_Q+I(Intens_Q^2)+TIME, data=d)

#-- Model selection table			
(ICtab(dec.PD,qdec, qdmntd,dc.mntd,qdrao,dcrao,qdhaed,dc.haed,decqmpd,dec.mpd,
			 dc.SimpPhyl,qdSimpPhyl, dc.Shannon, qdShannon, dc.Intens_Q, qIntens_Q,
			 type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE, nobs=504))

#---Ploting the best model
qdIntens_Q<-lm(dc~Intens_Q+I(Intens_Q^2), data=d)
anova(qdIntens_Q)
summary(qdIntens_Q)
coef(qdIntens_Q)
dec.pred<-predict(qdIntens_Q)#predicted values for each observation
d$dec.pred<-dec.pred#binding to the data frame

dec.plot1<-aggregate(dc, list(trata, d$Intens_Q), "mean")
dec.plot1$dec.pred<-aggregate(d$dec.pred, list(trata, d$Intens_Q), "mean")$x
names(dec.plot1)=c("Diversity", "Intens_Q", "Decompositon", "pred")#renaming the columns to make the plot
dec.plot1$dec.pred<-coef(qdIntens_Q)[1]+coef(qdIntens_Q)[2]*dec.plot1$Intens_Q+coef(qdIntens_Q)[3]*dec.plot1$Intens_Q^2#adding the predicted values 
#for PD by the quadratic model

#Just another version with data as points and model as line
ggplot(dec.plot1, aes(x=Intens_Q)) + 
	geom_line(aes(y=pred), colour="blue", linetype="solid", size=1) + 
	geom_point(aes(y=Decompositon), colour="red", size=4, shape=21, fill="white")