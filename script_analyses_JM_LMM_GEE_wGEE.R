#### Script Analyses for sex/gender and cognitive decline ####
## MELODEM 
## September 9th, 2020 
## Anaïs Rouanet and Cécile Proust-Lima

rm(list=ls())
setwd("./Path/")

#### 1. Libraries required #####
## please install the packages if needed

library(nlme)
library(JM)
library(lcmm)
library(geeM)
library(weightQuant)
library(splines)


#### 2. Creation of the data set for the analysis #######

# generic names for variables and data
# data = dataframe with the sample
#
# age0C = age at baseline minus 75 in decades
# age0 = age at baseline in years
# age = age in years
# PE = practice effect ie indicator that the visit is the first one
# male = indicator for male (sex/gender)
# EL = binary education level with low=0 and high=1
# time = delay in the study in decades (current time - time of entry of the participant)
# varY = outcome
# ID = Identifier of the participant - should be numeric to avoid any problems with packages
# timeEvent = time to first event between death, dropout or observed at the end of the study
#             indicated in decades
#             a person who dies more than 3 years after the last visit is considered to have dropped out at the last visit
#             (dementia is considered as dropout)
# indEvent = status indicator at the end of the study, which equals 2 if death; 1 if dropout; 0 if still in the cohort at the end of the follow up 
# prevalent = indicator for prevalent dementia at baseline
# dementia = indicator for incident dementia during the follow up
# agedem = age at dementia diagnosis in years
# visit = indicator of visits (1 2 3 ...)
# timeDeath = visit corresponding to time to death
# agedeath = age at death in years (NA if alive)
# ageend = age in years at the end of the study follow-up (NA is dropped out or dead before) i.e at last planned visit (max(visit))
# lastVisitAge = age at the last observed visit in years



#### For instance in paquid: 
# a. loading of the dataset (longitudinal format) without any selection
load("PAQUID_MELODEM.RData")
data <- paq

# specify the name of the dataset for outputs: 
dataset <- "Paquid"
# specify the name of the outcome for outputs:  
outcome <- "Verbal Fluency"

# creation of the variables
data$age0C <- (data$age0 - 75)/10 
data$PE <- as.numeric(data$visit==1) 
data$EL <- as.factor(data$EL)
data$time <- (data$age - data$age0)/10
data$varY <- data$ani_15

theoretical_visits <- unique(data$visit)[order(unique(data$visit))]
# make sure the visit indexes are in consecutive order, otherwise:
data$visit <- sapply(data$visit, function(x) which(theoretical_visits==x))
table(data$visit)

data$timeDeath <- sapply((data$agedeath-data$age0), function(x) ifelse(!is.na(x),
                                                                       ifelse(x<=max(theoretical_visits),min(which(theoretical_visits>=x)),
                                                                              length(theoretical_visits)+1),NA))

# Definition of timeEvent
# initialization
data$indEvent <- -1
data$ageEvent <- -1
# age of event is age of death for those who died (ind = 2)
data$ageEvent[!is.na(data$agedeath)] <- data$agedeath[!is.na(data$agedeath)]
data$indEvent[!is.na(data$agedeath)] <- 2
# age of event is the last visit observed for those not dropped out (ind = 0)
data$ageEvent[!is.na(data$ageend)] <- data$ageend[!is.na(data$ageend)]
data$indEvent[!is.na(data$ageend)] <- 0
# age of event is: min between dementia diagnosis and dropout 
# Careful: participants who die more than 3 years after their last visit are considered as dropped out
# For participants who drop out, add 0.1 to ageEvent to avoid any problem of dropout time concomittant with an observed timepoint
data$inddropout <- ((data$agedeath > data$lastVisitAge+3)| data$dementia==1 | is.na(data$ageend))
data$indEvent[which(data$inddropout==1)] <- 1
data$ageEvent[which(data$inddropout==1)] <- data$lastVisitAge[which(data$inddropout==1)]+0.1
data$ageEvent[which(data$dementia==1)] <- data$agedem[which(data$dementia==1)]

# time to Event from age of event in decades
data$timeEvent <- (data$ageEvent -data$age0)/10
data <- data[order(data$numero, data$age),]


#### 3. Selection of the population #######

# Participants >=65 years old
# with varY observed at baseline
# not demented at baseline
# no missing values for sex/gender, education and age at baseline

# To do so I select the data at baseline in data0 and select the ID to include in ID
data0 <- data[data$time==0,]
ID <- data0$ID[which((data0$age0C >= -1)&(!is.na(data0$varY))&(!(data0$prevalent==1))&(!is.na(data0$EL))&(!is.na(data0$male)))]
dataL <- data[which(data$ID%in%ID),]
# verification that there is no problem: the next computation should give 0
sum(ID - unique(dataL$ID))

# selection of the subjects who were not immediately censored 
# and the non missing observations for varY
dataJM <- dataL[!(dataL$timeEvent==0)&!is.na(dataL$varY)&(dataL$time < dataL$timeEvent),]

## create a  dataset with a unique line per subject
dataS <- unique(dataJM[,c("ID","male","timeEvent","indEvent",
                          "age0C","EL")])


### The sample size is: 
length(dataS$ID)
length(dataJM$ID)
### The number of repeated measures by subject: 
quantile(table(dataJM$ID[!is.na(dataJM$varY)]))


### Description : histogram at all visits and at baseline
hist(dataJM$varY,xlab="outcome",main="varY in STUDY")
hist(dataJM$varY[dataJM$time==0],xlab="outcome",main="varY in STUDY")

### Description: range of the data, mean
summary(dataJM$varY)


# description of gender 
summary(dataS$male)


########## JOINT MODEL WITH A SINGLE TYPE OF EVENT #######


#### Summarize the types of event in a single one
dataS$indSingleEvent <- (!(dataS$indEvent==0))


#### Linear mixed model with LME
DquadLME <- lme( varY ~ time + I(time^2) 
                 + PE + age0C   + male*(time + I(time^2)) 
                 + EL*(time + I(time^2)),  # put the adjustment after that and
                 # if you add some interactions with time, remember to put them also in 
                 # the Dform of the joint model (survival part)
                 random=~ time + I(time^2) | ID, 
                 data = dataJM,method='ML')
summary(DquadLME)


#### Survival model with CoxPH
coxFit <- 
  coxph(formula = Surv(timeEvent, indSingleEvent) ~ male + age0C + EL,
        data = dataS, 
        x = TRUE)
# x : outputs the design matrix associated with covariates
summary(coxFit)


#### Joint Model  with current level
JM_D_level <- 
  jointModel(lmeObject = DquadLME,
             survObject = coxFit,
             timeVar = "time", # time variable in the mixed model
             parameterization = "value", # dependence structure 
             method = "spline-PH-aGH", # baseline risk model
             verbose = TRUE, # trace output
             control = list(GHk = 5))
summary(JM_D_level)


#### Joint Model  with current level and current slope

dForm <- 
  list( fixed = ~ 1 + I(2*(time)) +
          male + male:I(2*(time))+
          EL + EL:I(2*(time)),
        indFixed = c(2:3, 8:11),
        random = ~ 1 + I(2*time),
        indRandom = 2:3)

Binit <- JM_D_level$coefficients
Binit$Dalpha <- 0

JM_D_both <-  jointModel(lmeObject = DquadLME,
                         survObject = coxFit,
                         timeVar = "time", # time variable in the mixed model
                         parameterization = "both", # dependence structure 
                         method = "spline-PH-aGH",# baseline risk model
                         verbose = TRUE, # trace output
                         derivForm = dForm,
                         control = list(GHk = 5),init=Binit)

summary(JM_D_both)
#check Assoct and Assoct.s





############### MIXED MODEL WITH HLME #################


#### Linear mixed model with HLME with quadratic shape of trajectory
DquadHLME <- hlme( varY ~ time + I(time^2) + PE + age0C 
                   + male*(time + I(time^2))
                   + EL*(time + I(time^2)),
                   random=~ time + I(time^2) , subject= "ID", 
                   data = dataJM)
summary(DquadHLME)
WaldMult(DquadHLME,pos=c(8,9))


#### Linear mixed model with HLME with shape of trajectory approximated by splines
nsknots <- as.vector(quantile(dataJM$time[dataJM$time>0],prob=c(0.01,0.33,0.67,0.99)))
dataJM$Z1 <- ns(dataJM$time,knots =nsknots[2:3],Boundary.knots = nsknots[c(1,4)])[,1]
dataJM$Z2 <- ns(dataJM$time,knots =nsknots[2:3],Boundary.knots = nsknots[c(1,4)])[,2]
dataJM$Z3 <- ns(dataJM$time,knots =nsknots[2:3],Boundary.knots = nsknots[c(1,4)])[,3]


DsplHLME <- hlme( varY ~ Z1 + Z2 + Z3 + PE + age0C 
                  + male*(Z1 + Z2 + Z3)
                  + EL*(Z1 + Z2 + Z3),
                  random=~ (Z1 + Z2 + Z3) , subject= "ID", 
                  data = dataJM)

summary(DsplHLME)
DsplHLME$loglik

plot(DsplHLME,which="fit",var.time="time")


############### GEE MODEL WITH GEEM #################


DquadGEE <- geem( varY ~ time + I(time^2) + PE + age0C 
                  + male*(time + I(time^2))
                  + EL*(time + I(time^2)),data=dataJM,
                  id=ID, family=gaussian,corstr="independence",
                  sandwich=T)
DquadGEE
summary(DquadGEE)

DsplGEE <- geem( varY ~ Z1 + Z2 + PE + age0C 
                 + male*(Z1 + Z2)
                 + EL*(Z1 + Z2),data=dataJM,
                 id=ID, family=gaussian,corstr="independence",
                 sandwich=T)
DquadGEE
summary(DsplGEE)

############### weighted GEE MODEL WITH GEEM #################

# Computation of the weights with monotone missing data
# weights_noIMD <- weightsMMD(data=dataJM,Y="varY",X1="male", X2=c("age0C","EL"),
#                             subject="ID", death="timeDeath", time="visit", 
#                             interval.death = 0, name= "weights_noIMD")$data
# weights_noIMD <- weights_noIMD[order(weights_noIMD$ID,weights_noIMD$visit),]

### GEE with no intermittent missing data
# AquadWGEE_noIMD <- geem(varY~ time + I(time^2)  + PE + age0C   
#                         + male*(time + I(time^2))
#                         + EL*(time + I(time^2)), id="ID" ,data = weights_noIMD, family=gaussian, 
#                         corstr="independence", weights= "weights_noIMD",
#                         sandwich=T)



# Computation of the weights with intermittent missing data (IMD)
weights_IMD <- weightsIMD(data=dataJM,Y="varY",X1="male",X2=c("age0C","EL"),
                          subject="ID",
                          death="timeDeath",
                          time="visit",
                          impute=mean(dataJM$varY[dataJM$visit==1],na.rm=T),
                          name="weights_IMD")$data

weights_IMD <- weights_IMD[order(weights_IMD$ID,weights_IMD$visit),]
hist(weights_IMD$weights_IMD,xlab="weights",main='Histogram of weights  (with IMD)')

### GEE with intermittent missing data
DquadwGEE <- geem(varY~ time + I(time^2)  + PE + age0C
                  + male*(time + I(time^2))
                  + EL*(time + I(time^2)),
                  id="ID" ,
                  data = weights_IMD, family=gaussian, 
                  corstr="independence", weights= "weights_IMD",
                  sandwich=T)

############### SAVE THE MODELS #######################

save(DquadGEE, JM_D_both, JM_D_level, DquadHLME, DsplHLME, DquadLME, DquadwGEE, dataJM, dataS,
     file=paste("models_JM_LMM_GEE_wGEE_",dataset,"_",outcome,".RData",sep=""))

# later, you can load those objects (and so not have to rerun them)
# by using for instance:
data_file <- paste("models_JM_LMM_GEE_wGEE_",dataset,"_",outcome,".RData",sep="")
load(data_file)

dataset <- "Paquid"
outcome <- "Verbal Fluency"

summary(DquadGEE)
summary(DquadwGEE)
summary(JM_D_both)
summary(DquadHLME)
############### PREDICTIONS ###########################


#### dataset for new data
datanew <- data.frame(time=seq(0,2,by=0.05),age0=75,EL=0,male=0,PE=0)
datanew$age0C <- (datanew$age0 - 75)/10

#### predictions in JM (with association through level and slope)
datanew$male=0
JMPfemale <- predict(JM_D_both,newdata=datanew)
datanew$male=1
JMPmale <- predict(JM_D_both,newdata=datanew)


### predictions from LMM

datanew$male <-0
Pfemale <- predictY(DquadHLME ,newdata = datanew,draws=T)
datanew$male <-1
Pmale <- predictY(DquadHLME ,newdata = datanew,draws=T)

datanew$male <-0
datanew$Z1 <- ns(datanew$time,knots =nsknots[2:3],Boundary.knots = nsknots[c(1,4)])[,1]
datanew$Z2 <- ns(datanew$time,knots =nsknots[2:3],Boundary.knots = nsknots[c(1,4)])[,2]
datanew$Z3 <- ns(datanew$time,knots =nsknots[2:3],Boundary.knots = nsknots[c(1,4)])[,3]

PSfemale <- predictY(DsplHLME ,newdata = datanew)
datanew$male <-1
PSmale <- predictY(DsplHLME ,newdata = datanew)


####  Predictions from GEE
X <- data.frame(int = rep(1,length(datanew$time)),
                time = datanew$time,
                time2 =(datanew$time^2),
                PE = datanew$PE, 
                age0C = datanew$age0C,
                male = datanew$male,
                EL = datanew$EL)

X$male <- 0
X$maleT <- X$male*X$time
X$maleT2 <- X$male*(X$time)^2
X$ELT <- X$EL*X$time
X$ELT2 <- X$EL*(X$time)^2
predGEEfemale <- as.matrix(X)%*%DquadGEE$beta
X$male <- 1
X$maleT <- X$male*X$time
X$maleT2 <- X$male*(X$time)^2
predGEEmale <- as.matrix(X)%*%DquadGEE$beta

#### Predictions from weighted GEE
X$male <- 0
X$maleT <- X$male*X$time
X$maleT2 <- X$male*(X$time)^2
predwGEEfemale <- as.matrix(X)%*%DquadwGEE$beta
X$male <- 1
X$maleT <- X$male*X$time
X$maleT2 <- X$male*(X$time)^2
predwGEEmale <- as.matrix(X)%*%DquadwGEE$beta



############## GENERAL PLOT FOR PREDICTIONS #####

# Change the name of the file
pdf(file=paste("pred_",dataset,"_",outcome,".pdf",sep=""),height = 5,width=7)

# limits for Y by default
ylim1 <- c(min(dataJM$varY),max(dataJM$varY))
# change the limits for Y axis for better contrast of the curves
#ylim1 <- c(2,8)
# change the limits for Y axis for better contrast of the curves
xlim <- c(min(dataJM$time,na.rm=T),max(dataJM$time,na.rm=T))

### Plot of the predictions for the mixed model (plain) and JM (dashed)
plot(Pfemale,shades=T,xlab="time in decades",ylab=outcome,col="purple3",lwd=2,ylim=ylim1,xlim=xlim,main=paste(outcome," in ",dataset, " by gender (for age0=75,EL=0)",sep=""),bty="n",legend=NULL)
plot(Pmale,shades=T,col="green3",add=T,lwd=2)

lines(JMPfemale ~ datanew$time, type="l",col="purple3",lty=2,lwd=2)
lines(JMPmale ~ datanew$time, col="green3",lty=2,lwd=2)


lines(predGEEfemale ~ datanew$time, type="l",col="purple3",lwd=2,lty=4)
lines(predGEEmale ~ datanew$time, col="green3",lty=4,lwd=2)

lines(predwGEEfemale ~ datanew$time, type="l",col="purple3",lwd=2,lty=3)
lines(predwGEEmale ~ datanew$time, col="green3",lwd=2,lty=3)

legend(legend = c("female","male"),x="topright",bty="n",lty=1,
       col=c("purple3","green3"),lwd=2)

legend(legend = c("GEE","wGEE","LMM","JM"),x="bottomleft",bty="n",lty=c(4,3,1,2),
       col=1,lwd=2)


dev.off()

#### graph comparison between splines and quadratic
pdf(file=paste("pred_",dataset,"_",outcome,"_quad_spl.pdf",sep=""),height = 5,width=7)

### Plot of the predictions for the mixed model (plain) and JM (dashed)
plot(Pfemale,shades=T,xlab="time in decades",ylab=outcome,col="purple3",lwd=2,ylim=ylim1,xlim=xlim,main=paste(outcome," in ",dataset, " by gender (for age0=75,EL=0)",sep=""),bty="n",legend=NULL)
plot(Pmale,shades=T,col="green3",add=T,lwd=2)
plot(PSfemale,shades=T,col="purple3",add=T,lwd=2,lty=5)
plot(PSmale,shades=T,col="green3",add=T,lwd=2,lty=5)

legend(legend = c("female","male"),x="topright",bty="n",lty=1,
       col=c("purple3","green3"),lwd=2)


legend(legend = c("LMM quad","LMM spl"),x="bottomleft",bty="n",lty=c(1,5),
       col=1,lwd=2)


dev.off()



# Change the name of the file
pdf(file=paste("pred_diff_",dataset,"_",outcome,".pdf",sep=""),height = 5,width=7)

# limits for Y by default
ylim1 <- c(-1,1)
# change the limits for Y axis for better contrast of the curves
#ylim1 <- c(2,8)


### Plot of the predictions for the mixed model (plain) and JM (dashed)
plot(Pmale$pred[,1]-Pfemale$pred[,1] ~ datanew$time, type='l',xlim=xlim,xlab="time in decades",ylab=paste("difference in ",outcome),col="blue",lwd=2,ylim=ylim1,main=paste("Difference in ",outcome, " (male-female) \n in ",dataset," (for age0=75,EL=0)",sep=""),bty="n",legend=NULL)

lines(JMPmale-JMPfemale ~ datanew$time, type="l",col="blue",lty=2,lwd=2)


lines(predGEEmale-predGEEfemale ~ datanew$time, type="l",col="blue",lwd=2,lty=4)

lines(predwGEEmale-predwGEEfemale ~ datanew$time, type="l",col="blue",lwd=2,lty=3)

abline(a=0,b=0,col="gray",lwd=2,lty=3)


legend(legend = c("GEE","wGEE","LMM","JM"),x="bottomleft",bty="n",lty=c(4,3,1,2),
       col=1,lwd=2)


dev.off()

### Comparison of Gender estimates across methods ---- 
summary(dataJM$time)
maxtime <- max(dataJM$time)

times_pred  <- c(0, 5, 9, 15)/10
times_pred<-times_pred[which(times_pred<=maxtime)]
diff_gender <- data.frame("time"=rep(times_pred,each=4), "estim"=rep(0,length(times_pred)*4),"se"=rep(0,length(times_pred)*4),
                          "method"=rep(c("LMM","JM","GEE","WGEE"),times=length(times_pred)))

diff_gender$method <- factor(diff_gender$method, levels = c("JM", "LMM", "WGEE", "GEE"))

###computation sd LME
n_param <- 11
LME_V                              <- matrix(0, length(DquadHLME$best), length(DquadHLME$best))
LME_V[upper.tri(LME_V, diag=TRUE)] <- DquadHLME$V
var_LME                            <- diag(LME_V)[c(6,8:9)] 
cov_LME                            <- c(LME_V[6,8], LME_V[6,9], LME_V[8,9])

###computation sd JM
JM_V <- solve(JM_D_both$Hessian)
var_JM                             <- diag(JM_V)[c(6,8:9)] 
cov_JM                            <- c(JM_V[6,8], JM_V[6,9], JM_V[8,9])

for(t in times_pred){
  #Linear mixed model
  diff_gender$estim[diff_gender$time==t & diff_gender$method=='LMM'] <-
    DquadHLME$best[c("male","time:male","I(time^2):male")]%*%c(1,t,t^2)
  
  diff_gender$se[diff_gender$time==t & diff_gender$method=='LMM']    <-
    sqrt(var_LME[1] + var_LME[2]*t^2 + var_LME[3]*(t^2)^2 +
           2*(cov_LME[1]*t + cov_LME[2]*t^2 + cov_LME[3]*t^3))
  
  #Linear joint model
  diff_gender$estim[diff_gender$time==t & diff_gender$method=='JM']  <-
    JM_D_both$coefficients$betas[c("male","time:male","I(time^2):male")]%*%c(1,t,t^2)
  
  diff_gender$se[diff_gender$time==t & diff_gender$method=='JM']     <-
    sqrt(var_JM[1] + var_JM[2]*t^2 + var_JM[3]*(t^2)^2 +
           2*(cov_JM[1]*t + cov_JM[2]*t^2 + cov_JM[3]*t^3))
  
  #GEE
  diff_gender$estim[diff_gender$time==t & diff_gender$method=='GEE'] <-
    DquadGEE$beta[c(6,8:9)]%*%c(1,t,t^2)
  
  diff_gender$se[diff_gender$time==t & diff_gender$method=='GEE']     <-
    sqrt(DquadGEE$var[6,6] + DquadGEE$var[8,8]*t^2 + DquadGEE$var[9,9]*(t^2)^2 +
           2*(DquadGEE$var[6,8]*t + DquadGEE$var[6,9]*t^2 + DquadGEE$var[8,9]*t^3))
  
  #Weighted GEE
  diff_gender$estim[diff_gender$time==t & diff_gender$method=='WGEE']<-
    DquadwGEE$beta[c(6,8:9)]%*%c(1,t,t^2)
  
  diff_gender$se[diff_gender$time==t & diff_gender$method=='WGEE']     <-
    sqrt(DquadwGEE$var[6,6] + DquadwGEE$var[8,8]*t^2 + DquadwGEE$var[9,9]*(t^2)^2 +
           2*(DquadwGEE$var[6,8]*t + DquadwGEE$var[6,9]*t^2 + DquadwGEE$var[8,9]*t^3))
}

# Table of sex/gender estimates for the 4 methods, at baseline, 5 and 9 years
estim <- diff_gender[diff_gender$time%in%c(0, 0.5,0.9)&diff_gender$method%in%c("LMM","JM"),]
estim$pval <- sapply(estim$estim/estim$se, function(x) ifelse(x<0,pnorm(x ,0,1  ,lower.tail = T)*2, pnorm(x ,0,1  ,lower.tail = F)*2 ))
estim


### Cumulative transition intensities for dropout and death ----
#(accounting for competing risks)
library(mstate)

# Definition of the possible transitions between states: Health, dropout, death
# Transition 1: health -> dropout
# Transition 2: health -> death
trans <- trans.illdeath()
colnames(trans) <- c("health", 'Dropout','Death')
rownames(trans) <- c("health", 'Dropout','Death')
trans[2,3]<-NA
trans

# Definition of the death and dropout indicators
table(dataS$indEvent)
dataS$IDeath   <- ifelse(dataS$indEvent ==2,1,0)
dataS$IDropout <- ifelse(dataS$indEvent ==1,1,0)
dataS$Death    <- dataS$timeEvent
dataS$Dropout  <- dataS$timeEvent
head(dataS)

dataAJ <- msprep(data = dataS, trans = trans, time = c(NA, "Dropout", "Death"), 
                 status = c(NA, "IDropout", "IDeath"), keep = c("male"),
                 id="ID")
head(dataAJ)
names(dataAJ)[1] <- "id"


prF <-  subset(dataAJ, male==0)
prM <- subset(dataAJ, male==1)
attr(prF, "trans") <- trans
attr(prM, "trans") <- trans


c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data=prF)
c1 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data=prM)

#estimation of transition hazards and associated (co)variances
msf0 <- msfit(c0, trans=trans)
msf1 <- msfit(c1, trans=trans)

# Estimated cumulative hazard values at all event times 

pt0 <- probtrans(msf0, predt=0)[[1]]
pt1 <- probtrans(msf1, predt=0)[[1]]


### Final plot ----
library(ggplot2)
library(gridExtra)
library(dplyr)
library(patchwork)

traj_plot <- data.frame("Pfemale"=Pfemale$pred[,1], "Pfemale_inf"=Pfemale$pred[,2],"Pfemale_sup"=Pfemale$pred[,3],
                        "times"=datanew$time,
                        "Pmale"=Pmale$pred[,1], "Pmale_inf"=Pmale$pred[,2],"Pmale_sup"=Pmale$pred[,3],
                        "JMPfemale"=JMPfemale, "JMPmale"=JMPmale,
                        "predGEEfemale"=predGEEfemale, "predGEEmale"=predGEEmale,
                        "predwGEEfemale"=predwGEEfemale, "predwGEEmale"=predwGEEmale)

traj_plot_melted_female <- reshape2::melt(traj_plot[, c(4, grep("female", colnames(traj_plot)))], 
                                          id.vars=c("times"), variable.name = "method", value.name = "outcome")
traj_plot_melted_female$method <- forcats::fct_recode(traj_plot_melted_female$method,
                                                      GEE = "predGEEfemale", wGEE = "predwGEEfemale", LMM = "Pfemale", JM = "JMPfemale")

traj_plot_melted_male <- reshape2::melt(traj_plot[, c(4, 5,6,7,9,11,13)], 
                                        id.vars=c("times"), variable.name = "method", value.name = "outcome")
traj_plot_melted_male$method <- forcats::fct_recode(traj_plot_melted_male$method,
                                                    GEE = "predGEEmale", wGEE = "predwGEEmale", LMM = "Pmale", JM = "JMPmale")

traj_plot_melted <- rbind.data.frame(cbind.data.frame(traj_plot_melted_female, sex="female"), 
                                     cbind.data.frame(traj_plot_melted_male, sex="male")
)

IC_index <- c(grep("_inf", traj_plot_melted$method), grep("_sup", traj_plot_melted$method))
traj_plot_melted_pointestimate <- traj_plot_melted[-IC_index, ]
traj_plot_melted_pointestimate$method <- factor(as.character(traj_plot_melted_pointestimate$method),
                                                levels = c("GEE", "wGEE", "LMM", "JM"))

traj_plot_melted_IC <- traj_plot_melted[IC_index, ]
summary(traj_plot_melted_IC)
traj_plot_melted_IC$method <- as.factor(as.character(forcats::fct_recode(traj_plot_melted_IC$method,
                                                                         inf = "Pfemale_inf", 
                                                                         inf = "Pmale_inf", 
                                                                         sup = "Pfemale_sup", 
                                                                         sup = "Pmale_sup")
))

traj_plot_melted_IC <- full_join(
  traj_plot_melted_IC %>% 
    filter(method=="inf") %>% 
    rename(inf = outcome) %>% 
    select(-method),
  traj_plot_melted_IC %>% 
    filter(method=="sup") %>% 
    rename(sup = outcome) %>% 
    select(-method),
  by = c("times", "sex")
)


p_traj <- ggplot(traj_plot_melted_pointestimate, aes(x=times)) +
  geom_ribbon(data = traj_plot_melted_IC,
              aes(ymin = inf, ymax = sup, fill=sex), alpha = 0.2) + 
  geom_line(aes(y=outcome, linetype=method, color=sex), size = 1) +
  scale_color_manual(values = c("purple3", "green3")) + 
  scale_fill_manual(values = c("purple3", "green3")) +
  scale_linetype_manual(values = c("twodash", "dotted","solid","dashed")) + 
  theme_bw() +
  theme(plot.title = element_text(size=15, face = "bold", hjust = 0.5)) +
  ylab(paste(outcome, " in ", dataset, sep="")) + xlab("Time in decades") +
  
  p_traj

p_estimates <- grid.arrange(p1, p2, p3, p4, ncol = 1)

p_estimates2 <- ggplot(diff_gender, aes(estim,  as.factor(method)))+
  xlab(paste("Male vs Female difference in ", outcome, sep='')) +
  geom_point() + 
  geom_segment(limits_low) +
  facet_wrap(~time, ncol = 1, labeller = function(labels){label_both(labels, sep = " = ")}) +
  ylab("Method") +
  theme_bw() 
p_estimates2 <- p_estimates2 + geom_vline(xintercept= 0, linetype="dashed", color="blue", size=1.2) 

data_Johansen_DO <- data.frame("pstate"=c(pt0$pstate2, pt1$pstate2),
                               "pstate_up" = c(pt0$pstate2, pt1$pstate2) + 1.96*c(pt0$se2, pt1$se2),
                               "pstate_low" = c(pt0$pstate2, pt1$pstate2) - 1.96*c(pt0$se2, pt1$se2), 
                               "times" = c(pt0$time, pt1$time),
                               "sex"=c(rep("Female",length(pt0$time)),rep("Male",length(pt1$time))))

data_Johansen_DC <- data.frame("pstate" = c(pt0$pstate3, pt1$pstate3),
                               "pstate_up" = c(pt0$pstate3, pt1$pstate3) + 1.96*c(pt0$se3, pt1$se3),
                               "pstate_low" = c(pt0$pstate3, pt1$pstate3) - 1.96*c(pt0$se3, pt1$se3), 
                               "times" = c(pt0$time, pt1$time),
                               "sex"=c(rep("Female",length(pt0$time)),rep("Male",length(pt1$time))))

p_Johansen_DO <- ggplot(data_Johansen_DO, aes(x=times)) +
  geom_ribbon(aes(ymin = pstate_low, ymax = pstate_up, fill=sex), alpha=0.3) +
  geom_line(aes(y=pstate, color=sex), size = 1)  + ylim(0,1) +
  scale_color_manual(values = c("purple3", "green3")) + 
  scale_fill_manual(values = c("purple3", "green3")) + 
  theme_bw() +
  theme(plot.title = element_text(size=15, face = "bold", hjust = 0.5)) +
  ylab("Dropout probability") + xlab("Time in decades") 


p_Johansen_DC <- ggplot(data_Johansen_DC, aes(x=times)) +
  geom_ribbon(aes(ymin = pstate_low, ymax = pstate_up, fill=sex), alpha=0.3) +
  geom_line(aes(y=pstate, color=sex), size = 1)  + ylim(0,1) +
  scale_color_manual(values = c("purple3", "green3")) + 
  scale_fill_manual(values = c("purple3", "green3")) + 
  theme_bw() +
  theme(plot.title = element_text(size=15, face = "bold", hjust = 0.5)) +
  ylab("Death probability") + xlab("Time in decades") 


(p_final <-  (p_traj + 
                p_estimates2) / (p_Johansen_DO + guides(color="none", fill="none") + p_Johansen_DC + guides(color="none", fill="none")) + 
    patchwork::plot_layout(heights = c(2,1)))

file <- paste("Plot ", outcome, " in ",dataset, ".pdf", sep="")
pdf(file, width=9, height = 7)
print(p_final)
dev.off()
