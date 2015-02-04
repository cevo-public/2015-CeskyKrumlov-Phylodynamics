# Tanja Stadler, Feb 4, 2015

# Solutions to the 3 questions in 02-Ebola.R

###########################
# Question 1

test<-trees[[1]]
# drop tips which do not belong to major Sierra Leone clade
droptip<-c("EBOV|KM034550|EM095|SierraLeone_EM|2014-05-25","EBOV|KM034554|G3676|SierraLeone_G|2014-05-27","EBOV|KM034559|G3680|SierraLeone_G|2014-05-28","EBOV|KM034563|G3687|SierraLeone_G|2014-05-28","EBOV|KM034561|G3683|SierraLeone_G|2014-05-28","EBOV|KM034562|G3686|SierraLeone_G|2014-05-28")
test<-drop.tip(test,droptip)
rootheight<-max(getx(test,sersampling=1)[,1])
x<-getx(test,sersampling=1)
times<-x[,1]
ttype<-x[,2]
out1<-optim(c(2,1,0.7),LikShiftsSTTebolaConst,times=times,ttype=ttype,sprobc=sprobc,deathfix=deathfix,cutoff=0)
out3<-optim(c(rep(2,3),1,0.7),LikShiftsSTTebola,times=times,ttype=ttype,sprobc=sprobc,deathfix=deathfix,cutoff=0)

test<-pchisq(2*(out1$value-out3$value),2)
test
# test is <0.95, ie no significant better fit of three rates compared to 1 rate at 5% level

###########################
# Question 2

# Evaluate commands under 2) EI model analysis

# model without incubation period; total duration of infection:
1/out1$par[2]*365
# incubation period:
epi[2]
# infectious time:
epi[3]
# total duration of infection almost same for model with / without incubation period:
epi[2]+epi[3]

# R0:
out1$par[1]/out1$par[2]
epi[1]

###########################
# Question 3

# evaluate commands under 3) superspreader model analysis

test<-pchisq(2*(out2$value-outSS$value),2)
test
# test is <0.95, ie no significant better fit of superspreading compared to 1 homogeneous mixing at 5% level
# test is >0.9, ie significant better fit of superspreading compared to 1 homogeneous mixing at 10% level

transSS<- outSS$par[1]+outSS$par[2]
transNS<- outSS$par[3]+outSS$par[2]*outSS$par[3]/outSS$par[1]

transSS/transNS
# or
outSS$par[1]/outSS$par[3]
# superspreaders transmit almost 5 times as much

outSS$par[1]/(outSS$par[1]+outSS$par[2])
# about 20% superspreaders
