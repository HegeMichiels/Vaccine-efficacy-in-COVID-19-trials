library(ggplot2)
library(survival)
library(survminer)

set.seed(1) #Seed to make the code reproducible 

#Number of patients at risk, cases and censored patients #####

#Number of patients at risk 
at.risk.V = c(21314,21230,21054,20481,19314,18377,17702,17186,15464,14038,12169,9591,6403,3374,1463,398,0)  #Patients at risk (per 7 days) in the vaccine arm
at.risk.PB = c(21258,21170,20970,20366,19209,18218,17578,17025,15290,13876,11994,9471,6294,3301,1449,398,0)  #Patients at risk (per 7 days) in the placebo arm

#Cumulative number of cases 
cum.cases.V = c(0,21,37,39,41,42,42,43,44,47,48,48,49,49,50,50,50) #Cumulative number of cases (per 7 days) in the vaccine arm
cum.cases.PB = c(0,25,55,73,97,123,143,166,192,212,235,249,257,267,274,275,275)  #Cumulative number of cases (per 7 days) in the placebo arm

#Number of cases 
cases.V = c()
cases.PB = c()
cases.V[1] = 0
cases.PB[1] = 0
for(i in 2:length(cum.cases.PB)){
  cases.V[i] = cum.cases.V[i]-cum.cases.V[i-1]
  cases.PB[i] = cum.cases.PB[i]-cum.cases.PB[i-1]
}
cases.V #Cases (per 7 days) in the vaccine arm
cases.PB #Cases (per 7 days) in the placebo arm

#Number of censored patients
censored.V = c()
censored.PB = c()
censored.V[1] = 0
censored.PB[1] = 0
for(i in 2:length(cum.cases.PB)){
  censored.V[i] = at.risk.V[i-1] - at.risk.V[i] - cases.V[i]
  censored.PB[i] = at.risk.PB[i-1] - at.risk.PB[i] - cases.PB[i]
}
censored.V  #Censored patients (per 7 days) in the vaccine arm
censored.PB #Censored patients (per 7 days) in the placebo arm

#Construct dataset ####

visits = seq(0,112,by= 7)  #Study visits in the Pfizer trial
censor.day = 112 #End of the trial
alpha = 28 #Ramp-up period

#Placebo data
PB.patients = 1:at.risk.PB[1]
PB.patients.at.risk = PB.patients #Patients at risk in the placebo arm at baseline
PB.X = rep(0,times = at.risk.PB[1])
PB.status = rep(0,times = at.risk.PB[1]) #at the start, nobody is censored
PB.time = rep(112,times = at.risk.PB[1]) #everybody's follow-up time is set at the end of the trial
for(i in 2:length(visits)){ #for-loop over all visits
  t.start = visits[i-1]+1 #start of this time interval
  t.end = visits[i] #end of this time interval
  #Cases in this interval
  cases = cases.PB[i] #number of cases in this interval
  cases.patients = sample(PB.patients.at.risk,size=cases) #patients who are cases in this interval are randomly chosen
  cases.time = sample(t.start:t.end,size=cases,replace=TRUE) #time of infection is chosen randomly in this interval for every patient
  PB.time[cases.patients] = cases.time
  PB.status[cases.patients] = 1 #for all cases, the status is set to 1
  PB.patients.at.risk = setdiff(PB.patients.at.risk, cases.patients) #cases are removed from the at risk set
  #Censored patients in this interval
  censored = censored.PB[i] #number of censored patients in this interval
  censored.patients = sample(PB.patients.at.risk,size=censored) #patients who are censored in this interval are randomly chosen
  censored.time = sample(t.start:t.end,size=censored,replace=TRUE) #time of censoring is chosen randomly in this interval for every patient
  PB.time[censored.patients] = censored.time
  PB.status[censored.patients] = 0
  PB.patients.at.risk = setdiff(PB.patients.at.risk, censored.patients) #censored patients are removed from the at risk set
}
PB.data = data.frame(PB.patients,PB.X,PB.time,PB.status,stringsAsFactors = FALSE)
colnames(PB.data) = c("patient","X","time","status")
head(PB.data)

#Vaccine data
V.patients = 1:at.risk.V[1]
V.patients.at.risk = V.patients #Patients at risk in the vaccine arm at baseline
V.X = rep(1,times = at.risk.V[1])
V.status = rep(0,times = at.risk.V[1]) #at the start, nobody is censored
V.time = rep(112,times = at.risk.V[1]) #everybody's follow-up time is set at the end of the trial
for(i in 2:length(visits)){
  t.start = visits[i-1]+1 #start of this time interval
  t.end = visits[i] #end of this time interval
  #Cases in this interval
  cases = cases.V[i] #number of cases in this interval
  cases.patients = sample(V.patients.at.risk,size=cases) #patients who are cases in this interval are randomly chosen
  cases.time = sample(t.start:t.end,size=cases,replace=TRUE) #time of infection is chosen randomly in this interval for every patient
  V.time[cases.patients] = cases.time
  V.status[cases.patients] = 1
  V.patients.at.risk = setdiff(V.patients.at.risk, cases.patients) #cases are removed from the at risk set
  #Censored patients in this interval
  censored = censored.V[i] #number of censored patients in this interval
  censored.patients = sample(V.patients.at.risk,size=censored) #patients who are censored in this interval are randomly chosen
  censored.time = sample(t.start:t.end,size=censored,replace=TRUE) #time of censoring is chosen randomly in this interval for every patient
  V.time[censored.patients] = censored.time
  V.status[censored.patients] = 0
  V.patients.at.risk = setdiff(V.patients.at.risk, censored.patients) #censored patients are removed from the at risk set
}
V.data = data.frame(V.patients,V.X,V.time,V.status,stringsAsFactors = FALSE)
colnames(V.data) = c("patient","X","time","status")
head(V.data)

data = rbind(PB.data,V.data)
data.pfizer = data
head(data)

#Save the dataset
setwd("~/OneDrive - UGent/Per protocol effect/Simulaties Janssen/Janssen data")
write.table(data.pfizer, file = "data.pfizer.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)


#Cumulative incidence plot #####

data.plot = data
data.plot$X = ifelse(data.plot$X==0,yes="Placebo",no="Vaccine")

p = ggsurvplot(survfit(Surv(time, status) ~ X, data = data.plot),
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
p$plot = p$plot + scale_x_continuous(breaks=c(seq(0,119,7),alpha))
p$plot = p$plot+ geom_vline(xintercept = alpha, linetype="dotted",  color = "pink", size=1)
p$plot = p$plot+ geom_vline(xintercept = 0, linetype="dotted",  color = "grey", size=1)
p$plot = p$plot+ geom_vline(xintercept = 21, linetype="dotted",  color = "grey", size=1)
p$table <- p$table + theme_cleantable()
p$table <- p$table + theme(plot.title = element_text(size = 12))
p$cumevents <- p$cumevents + theme_cleantable()
p$cumevents <- p$cumevents + theme(plot.title = element_text(size = 12))
p

