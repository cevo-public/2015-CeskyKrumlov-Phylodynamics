# Tanja Stadler, Feb 4, 2015

# Maximum likelihood analyses from
# "Insights into the early epidemic spread of Ebola in Sierra Leone provided by viral sequence data"
# T. Stadler, D. Kuhnert, D.A. Rasmussen, L. du Plessis. PLoS Currents: Outbreaks, 2014

path<- "/Users/tstadler/Documents/Data/Uni/Teaching/CeskyKrumlov/R"

library(ape)
library(TreePar)

setwd(paste(path,sep=""))
source("02-Ebola-Functions.R")

############################################
############################################
# read data trees

setwd(paste(path,sep=""))
trees<-read.tree("02-EbolaTrees.trees")

############################################
############################################
# 1) SKYLINE analysis

# estimate death rate:
deathfix<-0
# fix sampling probability to 0.7:
sprobc <- 0.7
# 3 Re intervals (one prior to most ancestral sample, then 50/50 length) 
numbRe<-3  # alternative numbRe<-1

#calculate likelihood for first 10 trees (we did it for all 90)
numbTree<-10 # full analysis with numbTree<-90

likdif<-vector()
estimates<-c()
for (index in 1:numbTree){
test<-trees[[index]]
# drop tips which do not belong to major Sierra Leone clade
droptip<-c("EBOV|KM034550|EM095|SierraLeone_EM|2014-05-25","EBOV|KM034554|G3676|SierraLeone_G|2014-05-27","EBOV|KM034559|G3680|SierraLeone_G|2014-05-28","EBOV|KM034563|G3687|SierraLeone_G|2014-05-28","EBOV|KM034561|G3683|SierraLeone_G|2014-05-28","EBOV|KM034562|G3686|SierraLeone_G|2014-05-28")
test<-drop.tip(test,droptip)
rootheight<-max(getx(test,sersampling=1)[,1])
x<-getx(test,sersampling=1)
times<-x[,1]
ttype<-x[,2]
if (numbRe == 1 ){
out<-optim(c(2,1,0.7),LikShiftsSTTebolaConst,times=times,ttype=ttype,sprobc=sprobc,deathfix=deathfix,cutoff=0)
} else {
out<-optim(c(rep(2,numbRe),1,0.7),LikShiftsSTTebola,times=times,ttype=ttype,sprobc=sprobc,deathfix=deathfix,cutoff=0)}
estimates<-rbind(estimates,parepi(out,sprobc=sprobc,deathfix=deathfix))
}
estimatesraw<-estimates

estimatesmedian<-vector()
estimatesHPD<-vector()
estimatesmean<-vector()
estimatesvar<-vector()
for (i in 1:length(estimates[1,])){
	estimatesmedian<-c(estimatesmedian,median(estimates[,i]))
	estimatesHPD<-cbind(estimatesHPD,HPD(estimates[,i]))
	estimatesmean<-c(estimatesmean,mean(estimates[,i]))
	estimatesvar<-c(estimatesvar,var(estimates[,i]))
}
estimates<-round(rbind(estimatesmedian,estimatesHPD,estimatesmean,estimatesvar,estimatesraw),2)
colnames(estimates)<-c("R1","R2","R3","rateUninf","sampProb")
rownames(estimates)<-c("median","95low","95up","mean","var",1:numbTree)
print(estimates)
         # R1   R2   R3 rateUninf sampProb
# median 1.61 1.18 0.98      4.50      0.7
# 95low  1.53 1.02 0.68      3.26      0.7
# 95up   1.85 2.06 1.72      6.49      0.7
# mean   1.62 1.29 1.04      4.81      0.7
# var    0.02 0.15 0.15      1.54      0.0
# 1      1.79 1.33 0.68      4.29      0.7
# 2      1.85 1.87 1.63      6.97      0.7
# 3      1.53 0.79 1.08      3.26      0.7
# 4      1.59 2.06 0.54      6.49      0.7
# 5      1.54 1.08 1.72      3.60      0.7
# 6      1.61 1.13 0.74      4.54      0.7
# 7      1.61 1.22 1.00      4.95      0.7
# 8      1.65 1.02 0.95      4.45      0.7
# 9      1.70 1.27 0.90      5.75      0.7
# 10     1.37 1.09 1.14      3.74      0.7

# Question 1:
# For tree 1, can you determine if 3 reproductive numbers Re compared to 1 reproductive number Re is more plausible?
# Use functions above and in the Macroevolution tutorial.

############################################
############################################
# 2) EI model analysis
# this analysis takes longish (we did each tree on a seperate node on the cluster)

# here we analyse tree 1 with fixed sampling proportion of .7
i<-1
sprob<-0.7
phylo<-trees[[i]]
droptip<-c("EBOV|KM034550|EM095|SierraLeone_EM|2014-05-25","EBOV|KM034554|G3676|SierraLeone_G|2014-05-27","EBOV|KM034559|G3680|SierraLeone_G|2014-05-28","EBOV|KM034563|G3687|SierraLeone_G|2014-05-28","EBOV|KM034561|G3683|SierraLeone_G|2014-05-28","EBOV|KM034562|G3686|SierraLeone_G|2014-05-28")
phylo<-drop.tip(phylo,droptip)

# initial parameters:
par<-c(1/5,2*1/5,1/5)
# for numerical reasons, we change branch lengths to units of days instead of years:
phylo$edge.length <- phylo$edge.length*365
out<-optim(c(par),LikTypesSTTebolaS,phylo=phylo,sprob=sprob)$par

# transform parameters to R0, incubation time, infectious time, sampling probability: 
epi<-c(out[2]/out[3],1/out[1],1/out[3],sprob)
# epi for tree 1 is 1.338087 2.215416 1.376223 0.700000

# Question 2:
# For tree 1:
# What is the duration of infection when using the model without incubation period and 1 Re?
# When using the model with incubation period: what is the duration of incubation / infectiousness? What is remarkable about the total duration of infection?
# What are the R0 estimates for tree 1 using model with / without incubation period?


############################################
############################################
# 3) superspreader model analysis
# this analysis takes longish (we did each tree on a seperate node on the cluster)

# here we analyse tree 1 with fixed sampling proportion of .7
i<-1
sprob<-0.7
phylo<-trees[[i]]
droptip<-c("EBOV|KM034550|EM095|SierraLeone_EM|2014-05-25","EBOV|KM034554|G3676|SierraLeone_G|2014-05-27","EBOV|KM034559|G3680|SierraLeone_G|2014-05-28","EBOV|KM034563|G3687|SierraLeone_G|2014-05-28","EBOV|KM034561|G3683|SierraLeone_G|2014-05-28","EBOV|KM034562|G3686|SierraLeone_G|2014-05-28")
phylo<-drop.tip(phylo,droptip)
phylo$edge.length <- phylo$edge.length*365

test<-addroot(phylo,(1/24))
x<-getx(test,sersampling=1)
times<-x[,1]
ttype<-x[,2]
# initial values for transmission rate, removal rate, sampling (sampling ignored as fixed to 0.7):
init <- c(2,1,0.7)
out2<-optim(init,LikShiftsSTTebolaConst,times=times,ttype=ttype,sprobc=sprob,deathfix=0,root=0,survival=0)
# initial values for transmission11, transmission12, transmission21, removal
init2<-c(out2$par[1]/2,out2$par[1]/2,out2$par[1]/2,out2$par[2])
# assumption here is: group 1 and group 2 may have different sizes. transmission rate within group 1 is transmission11. transmission rate from group 1 individual into group 2 is transmission12 = transmission11*const with const=sizeGroup2/sizeGroup1. transmission rate within group 2 is transmission22, transmission rate from group 2 individual into group 1 is transmission21 = transmission22/const
# we estimate transmission11, transmission12, transmission21, and we evaluate
# transmission22 = transmission21 * const = transmission21 * transmission12 / transmission11
# see also Stadler & Bonhoeffer, Phil Trans Roy Soc., 2013
outSS<-optim(init2,LikTypesSTTebolaSSpread,phylo=phylo,sprob=sprob)

# for tree 1:
# out2$par transmission 0.3426510 removal 0.2595218 (2.6061846)
# out2$value 318.2714
# outSS$par 0.18973391 0.75521212 0.03998312 0.23648649
# outSS$value 315.7676

# Question 3:
# Is there significance for superspreading at the 5% level? At the 10% level?

# Assuming we know that there is superspreading - from the values in outSS$par, can you determine:
# What is total transmission rate of superspreaders? And normal spreaders?
# How much more do superspreaders transmit compared to normal spreaders?
# What fraction of the population are superspreaders?
# (for answering these questions, make sure you understand what transmission22 = transmission21 * transmission12 / transmission11 means, see above)
