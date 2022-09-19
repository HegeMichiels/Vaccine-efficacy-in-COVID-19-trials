############
#Library####
############
library(ggplot2)
library(survminer)
library(survival)
library(NMOF)
library(survRM2)

###################
#Import dataset####
###################

#Janssen data ####

data.jnj = read.delim("~/OneDrive - UGent/Per protocol effect/Simulaties Janssen/Janssen data/data.jnj.txt")
data = data.jnj
censor.day = 125 #end of the trial
alpha = 14 #ramp-up period

#Cumulative incidence plot 

p = ggsurvplot(survfit(Surv(time, status) ~ X, data = data.jnj),
               title = "Janssen COVID-19 study",
               risk.table = TRUE,cumevents = TRUE,   break.time.by=7,
               risk.table.height = 0.15, cumevents.height = 0.15,
               risk.table.y.text = FALSE,
               cumevents.y.text = FALSE,
               fontsize = 3,
               risk.table.title = "No. at Risk",
               cumevents.title = "Cumulative No. of Events", 
               xlab = "Days after dose 1", 
               ylab = "Cumulative incidence",xlim=c(0,126),ylim=c(0,0.035),size=1,fun="event")
p$plot = p$plot + theme(plot.title = element_text(hjust = 0.5))
p$plot = p$plot + scale_x_continuous(breaks=c(seq(0,126,7),14))
p$plot = p$plot+ geom_vline(xintercept = 14, linetype="dotted",  color = "pink", size=1)
p$plot = p$plot+ geom_vline(xintercept = 0, linetype="dotted",  color = "grey", size=1)
p$table <- p$table + theme_cleantable()
p$table <- p$table + theme(plot.title = element_text(size = 12))
p$cumevents <- p$cumevents + theme_cleantable()
p$cumevents <- p$cumevents + theme(plot.title = element_text(size = 12))
p

#Pfizer data ####

data.pfizer = read.delim("~/OneDrive - UGent/Per protocol effect/Simulaties Janssen/Janssen data/data.pfizer.txt")
data = data.pfizer
censor.day = 112 #end of the trial
alpha = 28 #ramp-up period

p = ggsurvplot(survfit(Surv(time, status) ~ X, data = data.pfizer),
               title = "Pfizer COVID-19 study",
               risk.table = TRUE,cumevents = TRUE,   break.time.by=7,
               risk.table.height = 0.15, cumevents.height = 0.15,
               risk.table.y.text = FALSE,
               cumevents.y.text = FALSE,
               fontsize = 3,
               risk.table.title = "No. at Risk",
               cumevents.title = "Cumulative No. of Events", 
               xlab = "Days after dose 1", 
               ylab = "Cumulative incidence",xlim=c(0,119),ylim=c(0,0.024),size=1,fun="event")
p$plot = p$plot + theme(plot.title = element_text(hjust = 0.5))
p$plot = p$plot + scale_x_continuous(breaks=c(seq(0,119,7),28))
p$plot = p$plot+ geom_vline(xintercept = 28, linetype="dotted",  color = "pink", size=1)
p$plot = p$plot+ geom_vline(xintercept = 0, linetype="dotted",  color = "grey", size=1)
p$plot = p$plot+ geom_vline(xintercept = 21, linetype="dotted",  color = "grey", size=1)
p$table <- p$table + theme_cleantable()
p$table <- p$table + theme(plot.title = element_text(size = 12))
p$cumevents <- p$cumevents + theme_cleantable()
p$cumevents <- p$cumevents + theme(plot.title = element_text(size = 12))
p

#########################
#ITT vaccine efficacy####
#########################

t = censor.day #day at which the vaccine efficacy will be estimated

#ITT with cumulative incidence as risk measure ####

data = data.jnj #Another dataset (e.g. data.pfizer) can also be used
fit = survfit(Surv(time, status) ~ X, data = data) #Kaplan Meier survival curve
CI.PB= 1-summary(fit,time=c(t), extend = TRUE)$surv[1] #1-Survival probability in placebo arm
CI.V=  1-summary(fit,time=c(t), extend = TRUE)$surv[2] #1-Survival probability in vaccine arm
VE.ITT.CI = 1-CI.V/CI.PB
VE.ITT.CI

#ITT with hazard rate as risk measure ####

data = data.jnj #Another dataset (e.g. data.pfizer) can also be used
data$status = ifelse(data$time<t,yes=data$status,no=0) #administrative censoring after day of estimating the VE
data$time = ifelse(data$time<t,yes=data$time,no=t) #administrative censoring after day of estimating the VE
fit = coxph(Surv(time, status) ~ X, data = data) #Cox proportional hazards model
HR = exp(fit$coef["X"]) #Hazard ratio = exp(beta_X) in Cox model
VE.ITT.HR = 1-HR
VE.ITT.HR

#ITT with incidence rate as risk measure ####

data = data.jnj #Another dataset (e.g. data.pfizer) can also be used
data$status = ifelse(data$time<t,yes=data$status,no=0) #administrative censoring after day of estimating the VE
data$time = ifelse(data$time<t,yes=data$time,no=t) #administrative censoring after day of estimating the VE
data$log_time = log(data$time) #logarithm of the follow-up time will be used as offset in the Poisson model
fit = glm(status~ X, offset=log_time,data=data,family=poisson) #Poisson model
IR = exp(fit$coef["X"]) #Incidence ratio = exp(beta_X) in Poisson model
VE.ITT.IR = 1-IR
VE.ITT.IR

#########################################
#PP (removing cases) vaccine efficacy####
#########################################

#In this PP effect, cases observed during the ramp-up period are removed from the analysis set

t = censor.day #day at which the vaccine efficacy will be estimated
alpha = 14 #End of the ramp-up period (14 days in Janssen trial, 28 days in Pfizer trial)

#PP (removing cases) with cumulative incidence as risk measure ####

data = data.jnj #Another dataset (e.g. data.pfizer) can also be used
data$PP.set = ifelse(data$time<alpha&data$status==1,yes=0,no=1) #indicator for PP set
data = subset(data,PP.set==1) #patients in dataset who were not infected in [0,alpha]
fit = survfit(Surv(time, status) ~ X, data = data) #Kaplan Meier survival curve 
CI.PB = 1-summary(fit,time=c(t), extend = TRUE)$surv[1] #1-Survival probability in placebo arm of the PP dataset
CI.V =  1-summary(fit,time=c(t), extend = TRUE)$surv[2] #1-Survival probability in vaccine arm of the PP dataset
VE.PPremove.CI = 1-CI.V/CI.PB
VE.PPremove.CI

#PP (removing cases) with hazard rate as risk measure ####

data = data.jnj #Another dataset (e.g. data.pfizer) can also be used
data$PP.set = ifelse(data$time<alpha&data$status==1,yes=0,no=1) #indicator for PP set
data = subset(data,PP.set==1) #patients in dataset who were not infected in [0,alpha]
data$status =ifelse(data$time<t,yes=data$status,no=0) #administrative censoring after day of estimating the VE
data$time = ifelse(data$time<t,yes=data$time,no=t) #administrative censoring after day of estimating the VE
fit = coxph(Surv(time, status) ~ X, data = data) #Cox proportional hazards model
HR = exp(fit$coef["X"]) #Hazard ratio = exp(beta_X) in Cox model
VE.PPremove.HR = 1-HR
VE.PPremove.HR

#PP (removing cases) with incidence rate as risk measure ####

data = data.jnj #Another dataset (e.g. data.pfizer) can also be used
data$PP.set = ifelse(data$time<alpha&data$status==1,yes=0,no=1) #indicator for PP set
data = subset(data,PP.set==1) #patients in dataset who were not infected in [0,alpha]
data$status =ifelse(data$time<t,yes=data$status,no=0) #administrative censoring after day of estimating the VE
data$time = ifelse(data$time<t,yes=data$time,no=t) #administrative censoring after day of estimating the VE
data$log_time = log(data$time) #logarithm of the follow-up time will be used as offset in the Poisson model
fit = glm(status~ X, offset=log_time,data=data,family=poisson) #Poisson model
IR = exp(fit$coef["X"]) #Incidence ratio = exp(beta_X) in Poisson model
VE.PPremove.IR = 1-IR
VE.PPremove.IR

##########################################
#PP (censoring cases) vaccine efficacy####
##########################################

#In this PP effect, cases observed during the ramp-up period are censored 

t = censor.day #day at which the vaccine efficacy will be estimated
alpha = 14 #End of the ramp-up period (14 days in Janssen trial, 28 days in Pfizer trial)

#PP (censoring cases) with cumulative incidence as risk measure ####

data = data.jnj #Another dataset (e.g. data.pfizer) can also be used
data$status =ifelse(data$time<t,yes=data$status,no=0) #administrative censoring after day of estimating the VE
data$time = ifelse(data$time<t,yes=data$time,no=t) #administrative censoring after day of estimating the VE
data$status =ifelse(data$time<alpha,yes=0,no=data$status) #censoring if infected before alpha
fit = survfit(Surv(time, status) ~ X, data =data) #Kaplan Meier survival curve 
CI.PB= 1-summary(fit,time=c(t), extend = TRUE)$surv[1] #1-Survival probability in placebo arm of the PP dataset
CI.V =  1-summary(fit,time=c(t), extend = TRUE)$surv[2] #1-Survival probability in vaccine arm of the PP dataset
VE.PPcensor.CI = 1-CI.V/CI.PB
VE.PPcensor.CI

#PP (censoring cases) with hazard rate as risk measure ####

data = data.jnj #Another dataset (e.g. data.pfizer) can also be used
data$status =ifelse(data$time<t,yes=data$status,no=0) #administrative censoring after day of estimating the VE
data$time = ifelse(data$time<t,yes=data$time,no=t) #administrative censoring after day of estimating the VE
data$status =ifelse(data$time<alpha,yes=0,no=data$status) #censoring if infected before alpha
fit = coxph(Surv(time, status) ~ X, data = data) #Cox proportional hazards model
HR = exp(fit$coef["X"]) #Hazard ratio = exp(beta_X) in Cox model
VE.PPcensor.HR = 1-HR
VE.PPcensor.HR

#PP (censoring cases) with incidence rate as risk measure ####

data = data.jnj #Another dataset (e.g. data.pfizer) can also be used
data$status =ifelse(data$time<t,yes=data$status,no=0) #administrative censoring after day of estimating the VE
data$time = ifelse(data$time<t,yes=data$time,no=t) #administrative censoring after day of estimating the VE
data$status =ifelse(data$time<alpha,yes=0,no=data$status) #censoring if infected before alpha
data$log_time = log(data$time) #logarithm of the follow-up time will be used as offset in the Poisson model
fit = glm(status~ X, offset=log_time,data=data,family=poisson) #Poisson model
IR = exp(fit$coef["X"]) #Incidence ratio = exp(beta_X) in Poisson model
VE.PPcensor.IR = 1-IR
VE.PPcensor.IR

###########################################
#VE if ramp-up period can be eliminated####
###########################################

t = censor.day #day at which the vaccine efficacy will be estimated
alpha = 14 #End of the ramp-up period (14 days in Janssen trial, 28 days in Pfizer trial)

#Method 1: Only estimating psi (rho and alpha need to be specified) #### 

data = data.jnj #Another dataset (e.g. data.pfizer) can also be used

#1. Estimate KM survival curves
fit = survfit(Surv(time, status) ~ X, data = data)
summary(fit)

#2 (a) Initial estimate for psi 

#First we define 3 functions: 

#Function that returns the mapped survival probability under vaccination at a time t using the SDM (see paper)
#This probability is estimated using the survival function under placebo and the given parameters rho, psi and alpha
# fit: survival object fitted on the observed data
# rho: specified parameter for the SDM model (Indicates how much weaker the vaccine effect is during the ramp-up time. Should be in the interval [0,1])
# psi: specified parameter for the SDM model (Represents the vaccine effect. Higher values mean higher efficacy. 0 indicates no vaccine effect.)
# alpha: specified parameter for the SDM model (Length of the ramp-up time)
# t: time for which the mapped survival probability needs to be returned 
count.surv = function(fit,rho,psi,alpha,t){
  time1 = t/exp(rho*psi)
  surv1 = ifelse(t<alpha, yes=summary(fit,time=c(time1), extend = TRUE)$surv[1],no=0)
  time2 = (t+alpha*(exp(psi*(1-rho))-1))/exp(psi)
  surv2 = ifelse(t>=alpha,yes=summary(fit,time=c(time2), extend = TRUE)$surv[1],no=0)
  return(surv1+surv2)
}

#Function that returns the squared difference between observed and mapped vaccine survival probability at time t
#function is used to perform a grid search over a range of psi values (rho and alpha need to be specified)
# x: value for psi
# fit: survival object fitted on the observed data
# rho: specified parameter for the SDM model (Indicates how much weaker the vaccine effect is during the ramp-up time. Should be in the interval [0,1])
# alpha: specified parameter for the SDM model (Length of the ramp-up time)
# t: time for which the squared difference between observed vaccine and mapped vaccine survival probability needs to be returned
calculate.diff.1param = function(x,fit,rho,alpha,t){
  psi = as.numeric(x[1])
  #predict infection time under vaccine at time t using the SDM model
  pred = count.surv(fit=fit,rho=rho,psi=psi,alpha=alpha,t=t)
  #observed survival probability under vaccine at time t 
  obs = summary(fit,time=c(t), extend = TRUE)$surv[2] 
  #difference between the 2 probabilities
  diff = pred-obs
  output = diff^2 
  return(output)
}

#Function that returns a start value for psi between 0 and max.psi (rho and alpha need to be specified)
#The mapped and observed survival function are compared at time t (using the squared difference)
#A grid search is performed for 100 psi values in the interval [0,max.psi]
#The value for psi for which the difference between the mapped and observed survival function is minimized, is returned
# data: observed dataset
# fit: survival object fitted on the observed data
# rho: specified parameter for the SDM model (Indicates how much weaker the vaccine effect is during the ramp-up time. Should be in the interval [0,1])
# alpha: specified parameter for the SDM model (Length of the ramp-up time)
# t: time for which the squared difference between observed and mapped vaccine survival probability needs to be calculated
# max.psi: maximum possible psi value
initial.1param = function(data,fit,rho,alpha,t,max.psi){
  res <- gridSearch(fun=calculate.diff.1param,fit=fit,rho=rho,alpha=alpha,t=t, lower = 0, upper = max.psi, npar = 1, n = 100)
  initial.psi = res$minlevels
  return(initial.psi)
}  

#Now we can estimate an initial value for psi by minimizing the squared difference in
#survival probabilities at the end of the trial 
#rho, alpha and a maximum value for psi need to be specified:
rho = 0
alpha = 14
max.psi = 10
t = censor.day
initial.psi = initial.1param(data,fit,rho,alpha,t=t,max.psi=max.psi)
initial.psi

#2 (b) Update of this initial psi estimate 

#First we define 3 functions: 

#Function that simulates data under vaccine according to the SDM with specified rho,psi and alpha values
#This simulated dataset can be used to calculate the restricted mean survival time
# n: number of patients 
# fit: survival object fitted on the observed data
# rho: specified parameter for the SDM model (Indicates how much weaker the vaccine effect is during the ramp-up time. Should be in the interval [0,1])
# psi: specified parameter for the SDM model (Represents the vaccine effect. Higher values mean higher efficacy. 0 indicates no vaccine effect.)
# alpha: specified parameter for the SDM model (Length of the ramp-up time)
# censor.day: last visit of the trial
counterfactual.data = function(n,fit,rho,psi,alpha,censor.day){
  tdom<-seq(1,censor.day,by=1)
  #failure times under vaccine
  failtimes.V = c()
  u = runif(n) #draws from uniform distribution
  Surv =  c()  #survival function under vaccine 
  for(t in tdom){
    time1 = t/exp(rho*psi)
    surv1 = ifelse(t<alpha, yes=summary(fit,time=c(time1), extend = TRUE)$surv[1],no=0)
    time2 = (t+alpha*(exp(psi*(1-rho))-1))/exp(psi)
    surv2 = ifelse(t>=alpha,yes=summary(fit,time=c(time2), extend = TRUE)$surv[1],no=0)
    Surv[t] = surv1+surv2
  }
  for(i in 1:n){
    colsums = colSums(outer(Surv, u[i], `>`))
    failtimes.V[i] = ifelse(colsums==0, yes=0, no = tdom[colsums])  #Failure times are randomly drawn from the mapped survival function
  }
  data= data.frame(patient=1:n,X=0,time=failtimes.V)
  data$status =ifelse(data$time<censor.day,yes=1,no=0) #Patients are administratively censored after last visit of the trial
  return(data)
}

#Function that returns the squared difference in restricted mean survival time till time L
#function is used by optim
# x: value for psi
# fit: survival object fitted on the observed data
# data: observed dataset
# rho: specified parameter for the SDM model (Indicates how much weaker the vaccine effect is during the ramp-up time. Should be in the interval [0,1])
# alpha: specified parameter for the SDM model (Length of the ramp-up time)
# L: restricted mean survival time is estimated from baseline till time L
calculate.diff.RMST.1param = function(x,fit,data,rho,alpha,L){
  psi = as.numeric(x[1])
  data.count = counterfactual.data(n=20000,fit=fit,rho=rho,psi=psi,alpha=alpha,censor.day=censor.day) #simulated data under vaccine according to the SDM with specified rho,psi and alpha values
  data.count$X = 0 
  data.long = rbind(subset(data,X==1),data.count)
  obj = rmst2(time=data.long$time,status=data.long$status,arm=data.long$X,tau=L) #Restricted mean survival times for the observed and mapped vaccine arm
  diff = obj$RMST.arm1$result["RMST","Est."]-obj$RMST.arm0$result["RMST","Est."] #Difference in restricted mean survival time
  output = diff^2 + ifelse(psi<0,yes=1000,no=0) #penalty if psi is negative
  return(output)
}


#Function that estimates value for psi (rho and alpha need to be specified)
#Survival curves are compared using restricted mean survival time at time L
# data: observed dataset
# fit: survival object fitted on the observed data
# initial.values = c(psi.initial) initial value for psi
# rho: specified parameter for the SDM model (Indicates how much weaker the vaccine effect is during the ramp-up time. Should be in the interval [0,1])
# alpha: specified parameter for the SDM model (Length of the ramp-up time)
# L: restricted mean survival time is estimated from baseline till time L
estimate.RMST.1param = function(data,fit,initial.values,rho,alpha,L){
  psi = optim(par =initial.values,fn = calculate.diff.RMST.1param,fit=fit,data=data,rho=rho,alpha=alpha,L=L)$par
  return(psi)
}


#Now we can update the initial value for psi by minimizing the squared difference in
#in restricted mean survival times over de entire duration of the trial
psi = estimate.RMST.1param(data,fit,initial.values=c(initial.psi),rho,alpha,L=censor.day)
parameters = c(rho,psi,alpha)
names(parameters) = c("rho","psi","alpha")
parameters

#To check how good the SDM with the obtained parameters approximates the observed survival curve, a Kaplan Meier curve can be plotted
#Therefore, we define a function that plots the original CI curves and the predicted vaccine curve for certain rho,psi and alpha
# fit: survival object fitted on the observed data
# data: observed dataset
# rho: specified parameter for the SDM model (Indicates how much weaker the vaccine effect is during the ramp-up time. Should be in the interval [0,1])
# psi: specified parameter for the SDM model (Represents the vaccine effect. Higher values mean higher efficacy. 0 indicates no vaccine effect.)
# alpha: specified parameter for the SDM model (Length of the ramp-up time)
# censor.day: last visit of the trial
plot.SDM = function(fit,data,rho,psi,alpha,censor.day){
  fit = survfit(Surv(time, status) ~ X, data = data)
  data.count = counterfactual.data(n=100000,fit=fit,rho=rho,psi=psi,alpha=alpha,censor.day=censor.day)
  data.copy = data
  data.copy$X = ifelse(data.copy$X==0,yes="Placebo",no="Vaccine")
  data.count$X = "Modelled vaccine curve"
  data.long <<- rbind(data.copy,data.count)
  data.long$X <- factor(data.long$X, levels=c("Placebo", "Vaccine", "Modelled vaccine curve"))  
  p = ggsurvplot(survfit(Surv(time, status) ~ X, data = data.long), title = " ",
                 xlab = "Days", 
                 ylab = "Cumulative incidence",size=1,fun="event",palette = c("red", "blue","pink"))$plot 
  p = p+ geom_vline(xintercept = alpha, linetype="dotted",  color = "red", size=1)
  print(p)
}

rho = parameters["rho"]
psi = parameters["psi"]
alpha = parameters["alpha"]
plot.SDM(fit,data,rho,psi,alpha,censor.day)

#Finally, the vaccine efficacy 'if the ramp-up time can be eliminated' can be calculated by setting the ramp-up period to 0 days:

rho = parameters["rho"]
psi = parameters["psi"]
CI.V = 1-count.surv(fit,rho,psi,alpha=0,t) #1-Survival probability in vaccine arm if there would be no ramp-up time (alpha=0)
CI.PB = 1-summary(fit,time=c(t),extend=TRUE)$surv[1] #1-Survival probability in placebo arm 
VE.hypothetical.1param.CI = 1-(CI.V)/(CI.PB)
VE.hypothetical.1param.CI

#Method 2: Estimating rho and psi (alpha needs to be specified) #### 

data = data.jnj #Another dataset (e.g. data.pfizer) can also be used

#1. Estimate KM survival curves
fit = survfit(Surv(time, status) ~ X, data = data)
summary(fit)

#2 (a) Initial estimates for rho and psi

#First we define 2 functions: 

#Function that returns the sum of squared differences between observed and mapped vaccine survival probability at the visits in t.vector
#function is used to perform a grid search over a range of psi and rho values (alpha needs to be specified)
# x: vector with values for rho, psi
# fit: survival object fitted on the observed data
# alpha: specified parameter for the SDM model (Length of the ramp-up time)
# t.vector: vector with times for which the squared difference between observed vaccine and mapped vaccine survival probability needs to be returned
calculate.diff.2param = function(x,fit,alpha,t.vector){
  rho = as.numeric(x[1])
  psi = as.numeric(x[2])
  t1 = t.vector[1]
  t2 = t.vector[2]
  #predict infection time under vaccine at times t1 and t2 using the SDM model
  pred1 = count.surv(fit=fit,rho=rho,psi=psi,alpha=alpha,t=t1)
  pred2 = count.surv(fit=fit,rho=rho,psi=psi,alpha=alpha,t=t2)
  #observed survival probabilities under vaccine at times t1 and t2
  obs1 = summary(fit,time=c(t1), extend = TRUE)$surv[2] 
  obs2 = summary(fit,time=c(t2), extend = TRUE)$surv[2]
  #difference between the 2 probabilities
  diff1 = pred1-obs1
  diff2 = pred2-obs2
  output = diff1^2 + diff2^2 #sum of squared differences
  return(output)
}

#Function that returns a start value for rho between 0 and max.rho and a start value for psi between 0 and max.psi (alpha needs to be specified)
#The mapped and observed survival function are compared at the times in t.vector (using the sum of squared differences)
#A grid search is performed for 50 rho values in the interval [0,max.rho] and 50 psi values in the interval [0,max.psi]
#The values for rho and psi for which the sum of squared differences between the mapped and observed survival function is minimized, are returned
# data: observed dataset
# fit: survival object fitted on the observed data
# rho: specified parameter for the SDM model (Indicates how much weaker the vaccine effect is during the ramp-up time. Should be in the interval [0,1])
# alpha: specified parameter for the SDM model (Length of the ramp-up time)
# t.vector: times for which the squared difference between observed and mapped vaccine survival probability needs to be calculated
# max.rho: maximum possible rho value
# max.psi: maximum possible psi value
initial.2param = function(data,fit,alpha,t.vector,max.rho,max.psi){
  res <- gridSearch(fun=calculate.diff.2param,fit=fit,alpha=alpha,t.vector=t.vector, lower = c(0,0), upper = c(max.rho,max.psi), npar = 2, n = 50)
  initial.rho = res$minlevels[1]
  initial.psi = res$minlevels[2]
  output = c(initial.rho,initial.psi)
  names(output) = c("rho","psi")
  return(output)
}  

#Now we can estimate a initial values for rho and psi by minimizing the sum of squared differences in
#survival probabilities at the end of the ramp-up period and the end of the trial 
#alpha and a maximum value for rrho and psi need to be specified:
alpha = 14
max.rho = 1
max.psi = 10
t.vector = c(alpha,censor.day)
initial = initial.2param(data, fit,alpha,t.vector=t.vector,max.rho=max.rho,max.psi=max.psi)
initial.rho = initial["rho"]
initial.rho
initial.psi = initial["psi"]
initial.psi

#2 (b) Update of these initial rho and psi estimates 

#First we define 2 functions: 

#Function that returns the sum of squared differences in restricted mean survival time till the first and second time in L.vector
#function is used by optim
# x: vector with values for rho and psi
# fit: survival object fitted on the observed data
# data: observed dataset
# alpha: specified parameter for the SDM model (Length of the ramp-up time)
# L.vector: vector with two visits, restricted mean survival time is estimated from baseline till time these visits
calculate.diff.RMST.2param = function(x,fit,data,alpha,L.vector){
  rho = as.numeric(x[1])
  psi = as.numeric(x[2])
  data.count = counterfactual.data(n=20000,fit=fit,rho=rho,psi=psi,alpha=alpha,censor.day=censor.day) #simulated data under vaccine according to the SDM with specified rho,psi and alpha values
  data.count$X = 0
  data.long = rbind(subset(data,X==1),data.count) 
  obj = rmst2(time=data.long$time,status=data.long$status,arm=data.long$X,tau=L.vector[1]) #Restricted mean survival times for the observed and mapped vaccine arm from baseline till first visit in L.vector
  diff1 = obj$RMST.arm1$result["RMST","Est."]-obj$RMST.arm0$result["RMST","Est."] #Difference in restricted mean survival time
  obj = rmst2(time=data.long$time,status=data.long$status,arm=data.long$X,tau=L.vector[2]) #Restricted mean survival times for the observed and mapped vaccine arm from baseline till second visit in L.vector
  diff2 = obj$RMST.arm1$result["RMST","Est."]-obj$RMST.arm0$result["RMST","Est."] #Difference in restricted mean survival time
  output =diff1^2 + diff2^2 +ifelse(psi<0,yes=1000,no=0) +ifelse(rho<0,yes=1000,no=0)+ifelse(rho>1,yes=1000,no=0) #sum of squared differences and penalties if psi is negative or rho outside [0,1]
  return(output)
}


#Function that estimates values for rho and psi (alpha needs to be specified)
#Survival curves are compared using restricted mean survival time till the first and second time in L.vector
# data: observed dataset
# fit: survival object fitted on the observed data
# initial.values: vector with initial values for rho and psi
# alpha: specified parameter for the SDM model (Length of the ramp-up time)
# L.vector: vector with two visits, restricted mean survival time is estimated from baseline till time these visits
estimate.RMST.2param = function(data,fit,initial.values,alpha,L.vector){
  param = optim(par =initial.values,fn = calculate.diff.RMST.2param,fit=fit,data=data,alpha=alpha,L=L.vector)$par
  param = c(param,alpha)
  names(param) = c("rho","psi","alpha")
  return(param)
}


#Now we can update the initial values for rho and psi by minimizing the sum of squared differences
#in restricted mean survival times over halve and the entire duration of the trial
L.vector = c(censor.day/2,censor.day)
param =  estimate.RMST.2param(data,fit,initial.values=c(initial.rho,initial.psi),alpha = alpha,L.vector=L.vector)
param = c(param)    
names(param) = c("rho","psi","alpha")


#To check how good the SDM with the obtained parameters approximates the observed survival curve, a Kaplan-Meier curve can be plotted
rho = parameters["rho"]
psi = parameters["psi"]
alpha = parameters["alpha"]
plot.SDM(fit,data,rho,psi,alpha,censor.day)

#Finally, the vaccine efficacy 'if the ramp-up time can be eliminated' can be calculated by setting the ramp-up period to 0 days:
CI.V = 1-count.surv(fit,rho,psi,alpha=0,t) #1-Survival probability in vaccine arm if there would be no ramp-up time (alpha=0)
CI.PB = 1-summary(fit,time=c(t),extend=TRUE)$surv[1] #1-Survival probability in placebo arm 
VE.hypothetical.2param.CI = 1-(CI.V)/(CI.PB)
VE.hypothetical.2param.CI

#Method 3: Estimating rho, psi and alpha #### 

data = data.jnj #Another dataset (e.g. data.pfizer) can also be used

#1. Estimate KM survival curves
fit = survfit(Surv(time, status) ~ X, data = data)
summary(fit)

#2 (a) Initial estimates for rho, psi and alpha

#First we define 2 functions: 

#Function that returns the sum of squared differences between observed and mapped vaccine survival probability at the visits in t.vector
#function is used to perform a grid search over a range of psi, rho and alpha values 
# x: vector with values for rho, psi and alpha
# fit: survival object fitted on the observed data
# t.vector: vector with times for which the squared difference between observed vaccine and mapped vaccine survival probability needs to be returned
calculate.diff.3param = function(x,fit,t.vector){
  rho = as.numeric(x[1])
  psi = as.numeric(x[2])
  alpha = as.numeric(x[3])
  t1 = t.vector[1]
  t2 = t.vector[2]
  t3 = t.vector[3]
  #predict infection time under vaccine at times t1, t2 and t3 using the SDM model
  pred1 = count.surv(fit=fit,rho=rho,psi=psi,alpha=alpha,t=t1)
  pred2 = count.surv(fit=fit,rho=rho,psi=psi,alpha=alpha,t=t2)
  pred3 = count.surv(fit=fit,rho=rho,psi=psi,alpha=alpha,t=t3)
  #observed survival probabilities under vaccine at times t1, t2 and t3
  obs1 = summary(fit,time=c(t1), extend = TRUE)$surv[2] 
  obs2 = summary(fit,time=c(t2), extend = TRUE)$surv[2]
  obs3 = summary(fit,time=c(t3), extend = TRUE)$surv[2]
  #difference between the 2 probabilities
  diff1 = pred1-obs1
  diff2 = pred2-obs2
  diff3 = pred3-obs3
  output = diff1^2 + diff2^2 +diff3^2  #sum of squared differences
  return(output)
}


#Function that returns a start value for rho between 0 and max.rho, a start value for psi between 0 and max.psi and a start value for alpha between 0 and max.alpha
#The mapped and observed survival function are compared at the times in t.vector (using the sum of squared differences)
#A grid search is performed for 20 rho values in the interval [0,max.rho], 20 psi values in the interval [0,max.psi] and 20 alpha values in the interval [0,max.alpha]
#The values for rho, psi and alpha for which the sum of squared differences between the mapped and observed survival function is minimized, is returned
# data: observed dataset
# fit: survival object fitted on the observed data
# t.vector: vector with times for which the squared difference between observed and mapped vaccine survival probability needs to be calculated
# max.rho: maximum possible rho value
# max.psi: maximum possible psi value
# max.alpha: maximum possible alpha value
initial.3param = function(data,fit,t.vector,max.rho,max.psi,max.alpha){
  res <- gridSearch(fun=calculate.diff.3param,fit=fit,t.vector=t.vector, lower = c(0,0,0), upper = c(max.rho,max.psi,max.alpha), npar = 3, n = 20)
  initial.rho = res$minlevels[1]
  initial.psi = res$minlevels[2]
  initial.alpha = res$minlevels[3]
  output = c(initial.rho,initial.psi,initial.alpha)
  names(output) = c("rho","psi","alpha")
  return(output)
}

#Now we can estimate a initial values for rho, psi and alpha by minimizing the sum of squared differences in
#survival probabilities on a fourth, halfway and at the end of the trial 
# a maximum value for rho, psi and alpha need to be specified:
max.rho = 1
max.psi = 10
max.alpha = censor.day/4
t.vector = c(censor.day/4,censor.day/2,censor.day)
initial = initial.3param(data, fit,t.vector=t.vector,max.rho=max.rho,max.psi=max.psi,max.alpha = max.alpha)
initial.rho = initial["rho"]
initial.rho
initial.psi = initial["psi"]
initial.psi
initial.alpha = initial["alpha"]
initial.alpha

#2 (b) Update of these initial rho, psi and alpha estimates 

#First we define 2 functions: 

#Function that returns the sum of squared differences in restricted mean survival time till the first, second and third time in L.vector
#function is used by optim
# x: vector with values for rho, psi and alpha
# fit: survival object fitted on the observed data
# data: observed dataset
# L.vector: vector with three visits, restricted mean survival time is estimated from baseline till  these visits
calculate.diff.RMST.3param = function(x,fit,data,L.vector){
  rho = as.numeric(x[1])
  psi = as.numeric(x[2])
  alpha = as.numeric(x[3])
  data.count = counterfactual.data(n=20000,fit=fit,rho=rho,psi=psi,alpha=alpha,censor.day=censor.day) #simulated data under vaccine according to the SDM with specified rho,psi and alpha values
  data.count$X = 0
  data.long = rbind(subset(data,X==1),data.count)
  obj = rmst2(time=data.long$time,status=data.long$status,arm=data.long$X,tau=L.vector[1]) #Restricted mean survival times for the observed and mapped vaccine arm from baseline till first visit in L.vector
  diff1 = obj$RMST.arm1$result["RMST","Est."]-obj$RMST.arm0$result["RMST","Est."] #Difference in restricted mean survival time
  obj = rmst2(time=data.long$time,status=data.long$status,arm=data.long$X,tau=L.vector[2]) #Restricted mean survival times for the observed and mapped vaccine arm from baseline till second visit in L.vector
  diff2 = obj$RMST.arm1$result["RMST","Est."]-obj$RMST.arm0$result["RMST","Est."] #Difference in restricted mean survival time
  obj = rmst2(time=data.long$time,status=data.long$status,arm=data.long$X,tau=L.vector[3]) #Restricted mean survival times for the observed and mapped vaccine arm from baseline till third visit in L.vector
  diff3 = obj$RMST.arm1$result["RMST","Est."]-obj$RMST.arm0$result["RMST","Est."] #Difference in restricted mean survival time
  output =diff1^2 + diff2^2 + diff3^2 +ifelse(psi<0,yes=1000,no=0) +ifelse(rho<0,yes=1000,no=0)+ifelse(rho>1,yes=1000,no=0) +ifelse(alpha<0,yes=1000,no=0)+ifelse(alpha>censor.day/2,yes=1000,no=0) #sum of squared differences and penalties if psi is negative or rho outside [0,1] or alpha outside [0,censor.day/2] 
  return(output)
}


#Function that estimates values for rho, psi and alpha 
#Survival curves are compared using restricted mean survival times in L.vector
# data: observed dataset
# fit: survival object fitted on the observed data
# initial.values: vector with initial values for rho, psi and alpha
# L.vector: vector with three visits, restricted mean survival time is estimated from baseline till time these visits
estimate.RMST.3param = function(data,fit,initial.values,L.vector){
  param = optim(par =initial.values,fn = calculate.diff.RMST.3param,fit=fit,data=data,L=L.vector)$par
  param = c(param)
  names(param) = c("rho","psi","alpha")
  return(param)
}


#Now we can update the initial values for rho and psi by minimizing the sum of squared differences
#in restricted mean survival times over halve and the entire duration of the trial
L.vector = c(censor.day/4,censor.day/2,censor.day)
param = estimate.RMST.3param(data,fit,initial.values=c(initial.rho,initial.psi,initial.alpha),L.vector=L.vector)
param = c(param)    
names(param) = c("rho","psi","alpha")


#To check how good the SDM with the obtained parameters approximates the observed survival curve, a Kaplan-Meier curve can be plotted
rho = parameters["rho"]
psi = parameters["psi"]
alpha = parameters["alpha"]
plot.SDM(fit,data,rho,psi,alpha,censor.day)

#Finally, the vaccine efficacy 'if the ramp-up time can be eliminated' can be calculated by setting the ramp-up period to 0 days:
CI.V = 1-count.surv(fit,rho,psi,alpha=0,t) #1-Survival probability in vaccine arm if there would be no ramp-up time (alpha=0)
CI.PB = 1-summary(fit,time=c(t),extend=TRUE)$surv[1] #1-Survival probability in placebo arm 
VE.hypothetical.3param.CI = 1-(CI.V)/(CI.PB)
VE.hypothetical.3param.CI
