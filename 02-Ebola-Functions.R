# Tanja Stadler, Feb 4, 2015

# functions required to run script 02-Ebola.R

##########################
# summaize parameters

# calculate shortest 95% interval for estimates
HPD  <- function(x)
{
x2<-sort(x)
skip=ceiling(0.05*length(x2))
intlength<-vector()
for (i in 0:skip){
	intlength<-c(intlength,x2[length(x2)-skip+i] - x2[1+i])
	}
t <-which(intlength==min(intlength))[1]
c(x2[t],x2[length(x2)-skip-1+t])
}

# transform parameters returned by the optimization into R0, length of exposed period, length of infectious period, sampling probability
parepi<-function(out,sprobc,deathfix){
	temp<-out$par
	par<-vector()
	sampling<-temp[length(temp)]
	if (sprobc>0) {sampling<-sprobc}
	temp<-temp[-length(temp)]
	death<-temp[length(temp)]
	if (deathfix>0) {death<-deathfix}
	par<-temp[-length(temp)]/death
	par<-c(par,1/death*365,sampling)
	par
}

##########################
# estimate parameters
# all functions with sampling probability 0 prior to most ancestral tip at cutoff=(365*0.0658)
# freq = 0.0001 means that individual an hour before the root was infectious with frequency freq

# estimate the amount of superspreading
LikTypesSTTebolaSSpread<-function(parebola,phylo,sprob=0.7,cutoff=(365*0.0658)){
treesnew<-addroot(phylo,1/24)
out<-LikTypesSTT(parebola,treesnew,fix=rbind(c(4,6,7,8),c(-0.4,-5,0,0),c(0,1,0,0)),sampfrac=c(sprob,sprob),migr=0,unknownStates=TRUE,cutoff=cutoff)
out
}

# Fit BDEI model to tree; estimate all 4 parameters (gamma, lambda, death, sprob)
LikTypesSTTebola<-function(parebola,phylo,freq=0.0001,cutoff=(365*0.0658)){
sprob<-parebola[length(parebola)]
if (sprob>1){out<-10^12} else {
epsi<-10^(-10)
treesnew<-addroot(phylo,1/24)
treesnew$states<-rep(2,length(treesnew$tip.label))
out<-LikTypesSTT(c(parebola[-length(parebola)]),treesnew,fix=rbind(c(1,4,5,7,8),c(epsi,epsi,epsi,epsi,epsi)),sampfrac=c(0,sprob),migr=2,freq=freq,cutoff=cutoff)
}
out
}

# Fit BDEI model to tree; estimate 3 parameters (gamma, lambda, death); fix sprob
LikTypesSTTebolaS<-function(parebola,phylo,sprob,freq=0.0001,cutoff=(365*0.0658)){
epsi<-10^(-10)
treesnew<-addroot(phylo,1/24)
treesnew$states<-rep(2,length(treesnew$tip.label))
out<-LikTypesSTT(c(parebola),treesnew,fix=rbind(c(1,4,5,7,8),c(epsi,epsi,epsi,epsi,epsi)),sampfrac=c(0,sprob),migr=2,freq=freq,cutoff=(365*0.0658))
out
}

# Fit BDEI model to tree; estimate 3 parameters (gamma, lambda, sprob); fix death
LikTypesSTTebolaD<-function(parebola,phylo,death,freq=0.0001,cutoff=(365*0.0658)){
sprob<-parebola[length(parebola)]
if (sprob>1){out<-10^12} else {
epsi<-10^(-10)
treesnew<-addroot(phylo,1/24)
treesnew$states<-rep(2,length(treesnew$tip.label))
out<-LikTypesSTT(c(parebola[-length(parebola)]),treesnew,fix=rbind(c(1,4,5,6,7,8),c(epsi,epsi,epsi,death,epsi,epsi)),sampfrac=c(0,sprob),migr=2,freq=freq,cutoff=cutoff)
#LikTypesSTT(c(1,1,1,1,0.5,0.5,0,0),treesnew,sampfrac=c(0,sprob),migr=2,freq=0.5)
}
out
}

# Fit BDsky model to tree; estimate birth, death, sprob (unless sprobc > 0, then fixed to that value); length(parinit) is number of different birth intervals. First to ancestral sample, then equidistant between first and last sample.
LikShiftsSTTebola <- function(parinit,times,ttype,sprobc,root=1,deathfix=0, cutoff=(365*0.0658)){
sprobctemp<-parinit[(length(parinit))]
if (sprobc==0){sprobc<-sprobctemp}
parinit<-parinit[-(length(parinit))]
birth<-parinit[-(length(parinit))]
death<-parinit[(length(parinit))]
if (deathfix>0) {death<-deathfix}
if (cutoff==0){
shift<-1.001*max(times[which(ttype==0)])} else {shift<-cutoff}
shifttemp<-shift
shiftn<-(length(birth)-1)
sprob<-c(rep(sprobc,shiftn),0)
for (i in (shiftn-1):1){
	shift<-c(i/shiftn*shifttemp,shift)
}	
par<-c(birth,rep(death,length(birth)),shift)
out<-LikShiftsSTT(par,times,ttype,sprob=sprob,root=1)
out
}

# Fit BDsky model to tree; estimate birth, death, sprob (unless sprobc > 0, then fixed to that value); constant birth (ie no change in Re=birth/death).
LikShiftsSTTebolaConst <- function(parinit,times,ttype,sprobc,root=1,deathfix=0,cutoff=(365*0.0658),survival=1){
sprobctemp<-parinit[(length(parinit))]
if (sprobc==0){sprobc<-sprobctemp}
parinit<-parinit[-(length(parinit))]
birth<-parinit[-(length(parinit))]
death<-parinit[(length(parinit))]
if (deathfix>0) {death<-deathfix}
if (cutoff==0){
shift<-1.001*max(times[which(ttype==0)])} else {shift<-cutoff}
shifttemp<-shift
sprob<-c(sprobc,0)
par<-c(c(birth,birth),c(death,death),shift)
out<-LikShiftsSTT(par,times,ttype,sprob=sprob,root=root,survival=survival)
out
}