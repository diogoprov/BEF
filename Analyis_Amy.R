#---loading packages
library(spacodiR)
library(foodweb)
library(NetIndices)
library(lme4)
library(picante)
library(ecoPD)
library(apTreeshape)
library(phytools)

#---Data input and manipulation of trees to be dated in Phylocom
phylog.2002=read.tree(text="(((Myriophyllum_verticullatum,Ceratophyllum_demersum),(Elodea_nuttallii,(Potamogeton_zosteriformis,Potamogeton_crispus))),((Lithobates_catesbeianus,Lithobates_clamitans),((Physa_gyrina,Helisoma_trivolis),((Gammarus_americana,Hyalella_azteca),((Gyrinus_sp,Acilius_sp),(Belostoma_flumireum,(Trichocorixa_sp,(Ambrysus_sp,((Notonecta_undulata,Notonecta_sp),Neoplea_striola)))))))));")
phylog.2005=read.tree(text="(((((Sagittaria_rigida,Vallisneria_americana),Elodea_nuttallii),(Potamogeton_crispus,Potamogeton_natans)),(Ceratophyllum_demersum,(Utricularia_vulgaris,Myriophyllum_verticullatum))),((Lithobates_catesbeianus,Lithobates_clamitans),((Helisoma_trivolis,Physa_gyrina),((Crangonyx_richmondensis,Hyallela_azteca),((Coptotomus_sp,(Dineutus_sp,Gyrinus_sp)),(Belostoma_flumireum,((Hesperocorixa_sp,Trichocorixa_sp),(Ambrysus_sp,(Neoplea_striola,(Buenoa_sp,Notonecta_undulata))))))))));")
phyloNodes=makeNodeLabel(phylog.2002)
phyloNodes1=makeNodeLabel(phylog.2005)
write.tree(phyloNodes, file="phylo_2002.txt")
write.tree(phyloNodes1, file="phylo_2005.txt")
plot.phylo(phyloNodes, show.node.label=T)
plot.phylo(phyloNodes1, show.node.label=T)

#---Working wih dated phylogenies
amy.nature=read.tree("phylog_dated_2002.txt")
amy.nature=drop.tip(amy.nature, c("Vallisneria_americana", "Utricularia_vulgaris"))

plot(amy.nature, cex=0.8);axisPhylo()
nodelabels(node=c(127.099976), bg="black", frame="circle", cex=0.3)
c("Node1", "Node2", "Node3","Node4","Node5","Node8","Node10","Node12","Node14","Node15"),
title("Downing & Leibold (2002;2010)", xlab="Mya")
amy.ecology=read.tree("phylog_dated_2005.txt")
plot.phylo(amy.ecology, cex=0.8);axisPhylo()
title("Downing (2005)")

#---Playing with tree shape metrics
treeNature=as.treeshape(amy.nature)
colless.test(treeNature)#the tree is more balanced than predicted by the yule model 
z=ltt(amy.nature, gamma=F)
gammatest(z)#P=0.656
birthdeath(amy.nature)#fit a birth-death model of speciation
yule(amy.nature)#fit a birth model of speciation

#---Data about treatments, calculating PD
treat=read.table("treat.txt", h=T)
#--Independent on species richness
mpd.nature=mpd(t(treat),cophenetic(amy.nature))#mpd
obj=phylo4d(amy.nature, tip.data=treat)
cadotte=haed(obj)#Haed
raodiv=raoD(t(treat), amy.nature)$Dkk
#--Dependent on species richness
PD.nature=pd(t(treat), amy.nature)#Faith's PD
#PDp=pd(t(treat), phy.2002)#Faith's PD with topology only
mntd.nature=mntd(t(treat),cophenetic(amy.nature))#mnnd
SimpsonPhyl=simpson(obj, "phylogenetic")#Simpson Phylogenetic
Shannon=diversity(t(treat), "shannon")#Shannon
Simpson=diversity(t(treat), "simpson")#Simpson
J=sum(cophenetic.phylo(amy.nature))/(colSums(treat)^2)#see_Schweiger et al. 2008 Oecologia

#---Binding phylogenetic diversity variables with data file
d<-dados[order(dados$comp),]#sorting data by composition
head(d)
d$mpd<-rep(mpd.nature,each=12)#repeat MPD 12 times
d$haed<-rep(cadotte,each=12)
d$raodiv<-rep(raodiv,each=12)
d$mntd<-rep(mntd.nature,each=12)
d$SimpsonPhyl<-rep(SimpsonPhyl,each=12)
d$Shannon<-rep(Shannon,each=12)
d$Simpson<-rep(Simpson,each=12)
d$Intens_Q<-rep(J,each=12)

#---Plots variation in PD and MPD
plot(PD.nature[[1]]~PD.nature[[2]], xlab="Species richness", ylab="Faith's PD")
plot(PD.nature[[1]], xlab="Species composition", ylab="Faith's PD")

plot(mpd.nature~PD.nature[[2]], xlab="Species richness", ylab="MPD")
plot(mpd.nature, xlab="Species composition", ylab="MPD")

plot(cadotte~PD.nature[[2]], xlab="Species richness", ylab="Haed")
plot(raodiv~PD.nature[[2]], xlab="Species richness", ylab="Rao Entropy")
plot(mntd.nature~PD.nature[[2]], xlab="Species richness", ylab="MNTD")
plot(SimpsonPhyl~PD.nature[[2]], xlab="Species richness", ylab="Phylogenetic Simpson index")
plot(Shannon~PD.nature[[2]], xlab="Species richness", ylab="Shannon index")
plot(Simpson~PD.nature[[2]], xlab="Species richness", ylab="Simpson index")
plot(J~PD.nature[[2]], xlab="Species richness", ylab="Intensive Quadratic Entropy")

#---Distribution of species assigned to treatments on regional phylogeny
spa=as.spacodi(t(treat))
phy.dotplot(spa, phy=amy.nature, tips.adj=c(0.50,0.55), lab.adj=c(0,1))

#---Food Web indices Nature and FwB papers--> Not used
fw.nature=read.table("foodwebNature.txt", h=T)#including detritus, zoo- and phytoplankton
indices=GenInd(Tij=fw.nature)
write.table(fw.nature, file = "foodwebNature1.csv", append=FALSE, sep=",", col.names=FALSE, row.names=FALSE)
analyse.single(filename="foodwebNature1.csv")
plotweb(col=c("red", "green", "blue", "yellow"), radii=c(5,10,12, 12))

#--Without detritus, zoo- and phytoplankton Nature and FwB papers
fw.naturew=read.table("foodwebNature_w.txt", h=T)#excluding detritus, zoo- and phytoplankton
indices=GenInd(Tij=fw.naturew)
write.table(fw.naturew, file = "foodwebNatureW.csv", append=FALSE, sep=",", col.names=FALSE, row.names=FALSE)
analyse.single(filename="foodwebNatureW.csv")
plotweb(col=c("red", "green", "blue"), radii=c(5,10,12))

#---Ecology 2005
compos=read.table("Ecology2005.txt", h=T)
PD.ecology=pd(t(compos), amy.ecology)#Faith's PD

#---Food Web indices Ecology
fw.ecology=read.table("foodwebEcology_basal.txt", h=T)#including detritus, zoo- and phytoplankton
indices.eco=GenInd(Tij=fw.ecology)
write.table(fw.ecology, file = "foodwebEcology1.csv", append=FALSE, sep=",", col.names=FALSE, row.names=FALSE)
analyse.single(filename="foodwebEcology1.csv")
plotweb(col=c("red", "green", "blue", "yellow", "white"), radii=c(5,10,12, 12, 14))

#--Without detritus, zoo- and phytoplankton Ecology
fw.ecologyw=read.table("foodwebEcology.txt", h=T)#excluding detritus, zoo- and phytoplankton
indices.ecow=GenInd(Tij=fw.ecologyw)
write.table(fw.ecologyw, file = "foodwebecologyW.csv", append=FALSE, sep=",", col.names=FALSE, row.names=FALSE)
analyse.single(filename="foodwebecologyW.csv")
plotweb(col=c("red", "green", "blue", "yellow"), radii=c(5,10,12, 14))

#====data Nandao mandou 6 maio 2013
nandao=read.csv2("dados1.csv")
attach(nandao)
table (div, time) 
(a <- tapply (prod, INDEX = list(factor(div), factor(time)), FUN = mean))
(b <- tapply (prod, INDEX = factor(div), FUN = mean))
m3 <- lmer(prod ~ factor(div) + (1|comp) + (1|treat) + (1|time), data = nandao)
summary(m3)