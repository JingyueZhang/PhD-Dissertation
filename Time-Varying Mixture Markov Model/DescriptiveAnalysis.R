rm(list=ls()) # clear memory
library("markovchain")
library(TraMineR) # This library is designed for analyzing social-science sequence (i.e. life course sequence). The library is required for creating sequence object
library(seqHMM) # This library is used for Mixture MM and MM
library(ggplot2)
library(plyr)
##################################################################################################
# This script is used for descriptive analysis of the elderly population travel behavior analysis
# The data for this analysis is from five waves of NHTS (1990,1995,2001,2009,2017)
##################################################################################################

##################################################################################################
######################################  READ DATA ################################################
##################################################################################################
tempdat=read.csv(file="C:/Users/jiz13007/OneDrive - University of Connecticut/Pattern Recognition/NHTS sequence data/data with social demographic/transportation_project_full_data_with_weight.csv",header=T,sep=",")
dim(tempdat)
table(tempdat$R_AGE)
sum(tempdat$R_AGE>=65)/dim(tempdat)[1]

olddat=tempdat[tempdat$R_AGE>=65,]
dim(olddat)
write.csv(olddat, file = "C:/Users/jiz13007/OneDrive - University of Connecticut/Pattern Recognition/R script/Time varying MM/data/elderly data.csv")

rm(tempdat)
table(olddat$WAVE)
round(table(olddat$WAVE)*100/dim(olddat)[1],2)

olddat=read.csv(file="C:/Users/jiz13007/OneDrive - University of Connecticut/Pattern Recognition/R script/Time varying MM/data/elderly data.csv",header=T,sep=",")
dim(olddat)
###################################################################################################
#######################################  Activity-Travel Pattern ##################################
###################################################################################################
tseq<- seqdecomp(olddat,which( colnames(olddat)=="AGGSEQUENCE" ))
## convert to multiple rows 
tseq2<-seqdecomp(tseq, sep = "")

## label of states
data.lab <- c("Auto", "Public Transit", "Non-motorized",
              "Other Mode", "Home", "Mandatory", "Maintenance", "Discretionary", "Initial Home")
## value of states
data.alphab <-c("A","B","C","D","E","F","G","H","T")

myseq<- seqdef(tseq2, alphabet = data.alphab,
               labels = data.lab, xtstep = 9)

rm(tseq)
rm(tseq2)

#### Aggregate to five categories
varlist<-names(myseq)
myseq2=myseq
for (i in varlist){
  myseq2[[i]] = as.factor(revalue(myseq[[i]], c("T"="E","A"="A","B"="A","C"="A","D"="A",
                                                "E"="E","F"="F","G"="G","H"="H")))}

data.lab2 <- c("Travel", "Home", "Mandatory", "Maintenance", "Discretionary")
## value of states
data.alphab2 <-c("A","E","F","G","H")
aggseq<-seqdef(myseq2, alphabet = data.alphab2,
               labels = data.lab2, xtstep = 5)

### Plot
windows()
seqdplot(aggseq, border = NA, main = 'Activity-Travel Pattern of the Elderly')


###################################################################################################
#######################################  Transition over Time    ##################################
###################################################################################################

TTRAN=seqtrate(aggseq, sel.states = c("A","E","F","G","H"), time.varying = TRUE, weighted = FALSE,
               lag = 1, with.missing = FALSE, count = TRUE)

trancnt<-do.call(rbind, Map(data.frame,
                            ##A
                            AtoA=TTRAN[1,1,], AtoE=TTRAN[1,2,], AtoF=TTRAN[1,3,], AtoG=TTRAN[1,4,], AtoH=TTRAN[1,5,],
                            ##E
                            EtoA=TTRAN[2,1,], EtoE=TTRAN[2,2,], EtoF=TTRAN[2,3,], EtoG=TTRAN[2,4,], EtoH=TTRAN[2,5,], 
                            ##F
                            FtoA=TTRAN[3,1,], FtoE=TTRAN[3,2,], FtoF=TTRAN[3,3,], FtoG=TTRAN[3,4,], FtoH=TTRAN[3,5,], 
                            ##G
                            GtoA=TTRAN[4,1,], GtoE=TTRAN[4,2,], GtoF=TTRAN[4,3,], GtoG=TTRAN[4,4,], GtoH=TTRAN[4,5,], 
                            ##H
                            HtoA=TTRAN[5,1,], HtoE=TTRAN[5,2,], HtoF=TTRAN[5,3,], HtoG=TTRAN[5,4,], HtoH=TTRAN[5,5,]
                            ))

trancnt1<-transform(trancnt, sum=rowSums(trancnt))
table(trancnt1$sum)

## plot transition frequency over time
# Use the polygon function to draw the graph
time<-seq(1,95, by=1)

### CREATE PLOTS

varlist<-names(trancnt)

for (n in seq(1,length(varlist),by=3)){
  windows()
  par(mfrow=c(3,1))
  for (i in varlist[n:(n+2)]){
    plot( time , trancnt[[i]], col=rgb(0.2,0.1,0.5,0.9) , type="o" , lwd=3 , xlab="Time" , ylab=paste("OBS: ", i, sep=" ") , pch=20)
  }}


###################################################################################################
#######################################  State Count over Time  ###################################
###################################################################################################
n<-(1:96)
dataname1=sprintf("T%d",n)
colnames(aggseq)=dataname1
statename=c("A","E","F","G","H")

statecnt<-c()
for (stename in statename){
  Xcnt<-c()
  for (i in dataname1){
    cnt<-length(which(aggseq[[i]]==stename))
    names(cnt)<- i
    Xcnt<-c(Xcnt,cnt)}
  statecnt<-cbind(statecnt, Xcnt)
  colnames(statecnt)[ncol(statecnt)] <- paste(stename,"count",sep="_")}
View(statecnt)
table(rowSums(statecnt))

## plot transition frequency over time
time1<-seq(1,96, by=1)

varlist1<-colnames(statecnt)

windows()
par(mfrow=c(3,1))
for (i in 1:3){
  plot( time1 , statecnt[,i], col=rgb(0.2,0.1,0.5,0.9) , type="o" , lwd=3 , ylim = c(min(statecnt[,i]),max(statecnt[,i])),xlab="Time" , ylab=paste("OBS: ", varlist1[i], sep=" ") , pch=20)}

windows()
par(mfrow=c(2,1))
for (i in 4:5){
  plot( time1 , statecnt[,i], col=rgb(0.2,0.1,0.5,0.9) , type="o" , lwd=3 , ylim = c(min(statecnt[,i]),max(statecnt[,i])),xlab="Time" , ylab=paste("OBS: ", varlist1[i], sep=" ") , pch=20)}


########################################################################################################
############################################## Grouping ################################################
########################################################################################################
table(olddat$WORKER)
olddat$R_WORKER=ifelse((olddat$WORKER == 1),"worker", "non-worker")
table(olddat$R_WORKER)
windows()
seqdplot(aggseq, group = olddat$R_WORKER, border = NA)

sum(olddat$R_AGE <= 75)
sum(olddat$R_AGE > 75)
olddat$R_AGE2=ifelse((olddat$R_AGE <= 75),"younger", "older")
table(olddat$R_AGE2)
windows()
seqdplot(aggseq, group = olddat$R_AGE2, border = NA)
