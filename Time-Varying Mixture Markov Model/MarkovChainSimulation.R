rm(list=ls()) # clear memory
library("markovchain")
library(TraMineR) # This library is designed for analyzing social-science sequence (i.e. life course sequence). The library is required for creating sequence object
library(seqHMM) # This library is used for Mixture MM and MM
library(gtools)
library('cluster')
library(plyr)
library(HiddenMarkov)
library("optimization")
## This script is used for a simulation study analysis for time-varying mixture Markov model
## The model estimation process applied an EM algorithm
###################################################################################################
############################ Simulate Time-Varying Markov Chain ###################################
###################################################################################################
#Step 1: Define basic inputs
#m:number of state
m<-3
#Sname: name of the state
Sname<-c("A","B","C")
#n:number of subjects
n<-250
#T: length of the sequence
lenseq<-100

#x: exogeneous variables
x_mid<-c(rep(0,30),rep(1,30), rep(0,40)) # mid-day indicator
x_end<-c(rep(0,60),rep(1,40)) # evening time indicator
exog<-do.call(rbind, Map(data.frame,
                      xmid=x_mid,
                      xend=x_end))
table(exog)
dim(exog)
#p: number of exogenous variable
numexo<-dim(exog)[2]

#numcomp: number of components
numcomp<-2

######################################################################################################
#############################         VALUES FOR PARAMETERS       ############################
######################################################################################################
## STEP1: CREATE COLUMN NAMES OF TRANSITION/INITIAL MATRIX, EXOGENEOUS VARIABLES ETC.

# column names of transition probabilities
TRANLIST<-c()
for (lter in Sname){
  for (lter2 in Sname){
    tranname<-paste0(lter,"to",lter2)
    TRANLIST<-c(TRANLIST,tranname)}}

numpair<-length(TRANLIST)


#Column names of initial probabilities
INITLIST<-c()
for (j in 1:m){
  initname<-paste(Sname[j],"0",sep = "")
  INITLIST<-c(INITLIST,initname)
}

#column names of exogeneous variables
XLIST<-names(x)

#column names of random components
COMPLIST<-c()
for (numc in 1:numcomp){
  compname<-paste("Component",numc,sep = "_")
  COMPLIST<-c(COMPLIST,compname)
}


# STEP2 : SET INITIAL VALUES FOR PARAMETERS

# Set initial values for coefficients of exogeneous variables (gamma)
initgamma_function<-function(initgamma,num_state,num_pair,num_exogeneous){
  
  gamma<-matrix(c((rep(initgamma,num_pair))),
                nrow = num_pair, ncol = num_exogeneous,byrow = TRUE,
                dimnames = list(TRANLIST,XLIST))
  
  for (indx in 1:num_state){
    same_indx=(indx-1)*num_state+indx
    gamma[same_indx,]=0
  }
  
  return(gamma)}

gamma<-initgamma_function(c(2,3),m,numpair,numexo)
dfgamma<-as.data.frame(gamma)

## SET CONSTRAINS OF TRANSITION (if some transitions are not possible, make it always be 0)
# dfgamma$CONSTRAIN=c(rep(c(rep(0,8),1),8),rep(0,9))
# table(dfgamma$CONSTRAIN)
# 
# dfgamma[(dfgamma$CONSTRAIN == 1),1:numexo]<-NA
#dfgamma$CONSTRAIN=0 # No constrain since it is 15 minutes interval any transition can happen
# constrain_list<-c("EtoF","EtoG","EtoH",
#                   "FtoE","FtoG","FtoH",
#                   "GtoE","GtoF","GtoH",
#                   "HtoE","HtoG","HtoF")
# dfgamma[constrain_list, "CONSTRAIN"]=1
# dfgamma[(dfgamma$CONSTRAIN == 1), 1:numexo]<-NA
# dfgamma

# Set initial values for component sepecific intercepts (c)
c_g1<-c(0,-2,-3,-4,0,-5,-6,-7,0)
c_g2<-c(0,1,2,2,0,1,1,2,0)
c<-matrix(cbind(c_g1,c_g2),
          nrow=numpair, ncol=numcomp,byrow=FALSE,
          dimnames= list(TRANLIST, COMPLIST))

### last state is constraint to be 0
c0<-matrix(c(rep(4,(m-1)),0,rep(4,(m-1)),0), 
           nrow=m, ncol=numcomp, byrow = FALSE,dimnames = list(INITLIST,COMPLIST))

# Set initial values for component probability (pi)
pi<-c(1/2, 1/2)


######################################################################################################
#############################   CALCULATE INITIAL/TRANSITION PROBABILITIES       #####################
######################################################################################################
#Calculate initial probabilities
delta_values<-function(c0_maxtrix,num_state,num_component){
  delta<-matrix(-99, nrow = num_state, ncol = num_component,
                dimnames = list(Sname,COMPLIST))
  for (g in 1:dim(c0_maxtrix)[2]){# calculate initial probablities for each component
    for (j in 1:dim(c0_maxtrix)[1]){
      cj0g<-c0_maxtrix[j,g]
      deltaj0g<-(exp(cj0g)/(1+(sum(exp(c0_maxtrix[1:(num_state-1),g])))))
      delta[j,g]<-deltaj0g}}
  return (delta)}

delta<-delta_values(c0,m,numcomp)
sum(delta[,2])

round(delta,digits=3)


######### transition probability is a funtion of gamma and c
#### Calculate transition probability at each time point
# if initial tranisiton probability value is 0, then do not change, keep as 0

# CREATE COLUMN NAMES FOR TRANSITON PROBABILITIES

#qtjkg=exp(x_t*gamma_jk+c_jk)/1+sum(h=1:m; h is not j){exp(x_t*gamma_jh+c_jh)}
q_function<-function(gamma_matrix,c_matrix,exog_matrix,num_state,num_component,seqlength){
  TGLIST<-c()
  for (g in 1:num_component){
    for(t in 1:seqlength){
      varname<-paste("T",t,"Comp",g, sep="_")
      #print(t)
      TGLIST<-c(TGLIST,varname)}}
  
  num_pair<-num_state*num_state # number of transition pairs
  p<-dim(exog_matrix)[2]   # number of exogenous variables
  
  q<-matrix(-99, nrow = num_pair, ncol = (seqlength*num_component),byrow=TRUE,
            dimnames = list(TRANLIST,TGLIST)) #create a matrix to keep transition probability
  
  
  for (t in 1:seqlength){# for each time point
    exog_t<-as.vector(exog_matrix[t,]) # a vector (1*p) at time point t
    for (g in 1: num_component){ # for each component
      for(j in 1:num_state){# for each state
        strt_indx<-(j-1)*num_state+1
        end_indx<-j*num_state
        
        c_jg<-c_matrix[strt_indx:end_indx,g] #1 by numstate component specific intercept 
        gamma_j<-as.matrix(gamma_matrix[strt_indx:end_indx,1:p]) # numstate by (p) matrix, the last column is constrain
        temp_gamma<-rowSums(gamma_j%*%diag(exog_t)) #1 by 5
        constrain_indx<-which(is.na(temp_gamma))
        ## within the state j
        for (k in 1:num_state){ # loop through each pair
          if(k %in% constrain_indx){
            q_jkg=0
            q[strt_indx+k-1,(t+(g-1)*seqlength)]=q_jkg}
          else{
            c_jkg<-c_jg[k] # single number
            gamma_jk<-temp_gamma[k] # 1 single value
            q_jkg<-exp(gamma_jk+c_jkg)/(1+sum(exp(temp_gamma[-c(j,constrain_indx)]+c_jg[-c(j,constrain_indx)])))
            q[strt_indx+k-1,(t+(g-1)*seqlength)]=q_jkg}
        } # close of k loop
      }# close of j loop
    } # close of g loop
  }# close of t loop
  
  return (q)
}  # end of q function

q<-q_function(dfgamma,c,exog,m,numcomp,lenseq)
dim(q)

View(round(q,digits=4))

#############################################################
####  POPULATION BASED SIMULATION APPROACH ##################
###############################################################

## Initial Value
set.seed(1234)
initsim1<-sample(Sname, 250, replace = TRUE, prob = delta[,1])
table(initsim1)
set.seed(1234)
initsim2<-sample(Sname, 250, replace = TRUE, prob = delta[,2])
table(initsim2)


## Function for Markov Chain Generation
DataGen<-function(cnt,len, deltamatrix,qmatrix, statename,Rseed){
  ## Initial value
  set.seed(Rseed)
  oldword <-sample(statename,cnt,replace = TRUE, prob = deltamatrix)
  
  dataholder=matrix(-999, nrow=cnt, ncol = (len-1),dimnames = list(seq(1,cnt),seq(2,len)))
  dataholder=cbind(oldword,dataholder)
  for (t in 2:len){ # loop through each time point
    #for(t in 2:5){
    #for(i in 1:1){
    for (i in 1:m){ # loop through each state
      key = statename[i]
      tempProb = qmatrix[(((i-1)*m+1):(i*m)),t] # take the transition probability of that particular state at time t
      #print(tempProb)
      tempN=sum(oldword==key)
      #print (tempN)
      set.seed(Rseed)
      tempNew<-sample(Sname,tempN,replace=TRUE, prob = tempProb)
      #print(table(tempNew))
      ## append data based on key
      dataholder[dataholder[,(t-1)]==key,t]=tempNew}
    oldword = dataholder[,t] # update oldword
    #print (table(oldword))
  }
  
  return (dataholder)}

sim1<-DataGen(250,lenseq,delta[,1],q[,1:100],Sname,(1234))
dim(sim1)

round(table(sim1[,1])/250, digits = 4)
round(table(sim1[,2])/250, digits = 4)
round(table(sim1[,62])/250, digits = 4)
round(table(sim1[,82])/250, digits = 4)

sim2<-DataGen(250,lenseq,delta[,2],q[,101:200],Sname,(1234))

########################################################################
#########  EXACT FREQUENCIES BASED SIMULATION ##########################
########################################################################

## Function for Markov Chain Generation
FreqGen<-function(cnt,len, deltamatrix,qmatrix, statename,Rseed){
  round_preserve_sum <- function(x) {
    y <- floor(x)
    indices <- tail(order(x-y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    return(y)
  }
  
  listGen <- function(statename, freqcnt,seed){
    listholder<-c()
    for (i in 1:m){
      key=statename[i]
      keyN=freqcnt[i]
      templist<-rep(key, keyN)
      listholder<-c(listholder,templist)
    }
    set.seed(seed)
    ranlist<-sample(listholder)
    return(ranlist)
  }
  
  ## set initial value
  freqN=round_preserve_sum(deltamatrix*cnt)
  oldword<-listGen(statename,freqN,Rseed )
  table(oldword)
  
  dataholder=matrix(-999, nrow=cnt, ncol = (len-1),dimnames = list(seq(1,cnt),seq(2,len)))
  dataholder=cbind(oldword,dataholder)
  for (t in 2:len){ # loop through each time point
    #for(t in 2:2){
    #for(i in 1:3){
    for (i in 1:m){ # loop through each state
      key = statename[i]
      tempProb = qmatrix[(((i-1)*m+1):(i*m)),t] # take the transition probability of that particular state at time t
      #print(tempProb)
      tempN=sum(oldword==key)
      #print (tempN)
      freqN2=round_preserve_sum(tempProb*tempN)
      #print (freqN2)
      tempNew=listGen(statename, freqN2,Rseed)
      #print (table(tempNew))
      ## append data based on key
      dataholder[dataholder[,(t-1)]==key,t]=tempNew}
    oldword = dataholder[,t] # update oldword
    #print (table(oldword))
  }
  
  return (dataholder)}


Fsim1<-FreqGen(250,lenseq,delta[,1],q[,1:100],Sname,(1234))
dim(Fsim1)

round(table(Fsim1[,1])/250, digits = 4)
round(table(Fsim1[,2])/250, digits = 4)
round(table(Fsim1[,62])/250, digits = 4)
round(table(Fsim1[,82])/250, digits = 4)

Fsim2<-FreqGen(250,lenseq,delta[,2],q[,101:200],Sname,(1234))

###################################################################################################
####################################### Plot Pattern ##################################
###################################################################################################

data.lab2 <- c("A", "B", "C")
## value of states
data.alphab2 <-c("A","B","C")
seqSim1<-seqdef(sim1, alphabet = data.alphab2,
                labels = data.lab2, xtstep = 3)

FseqSim1<-seqdef(Fsim1, alphabet = data.alphab2,
                 labels = data.lab2, xtstep = 3)

### Plot
windows()
seqdplot(seqSim1, border = NA, main = 'Pattern of Group 1 (Sample Based)')
windows()
seqdplot(FseqSim1, border = NA, main = 'Pattern of Group 1 (Frequency Based)')


seqSim2<-seqdef(sim2, alphabet = data.alphab2,
                labels = data.lab2, xtstep = 3)

FseqSim2<-seqdef(Fsim2, alphabet = data.alphab2,
                 labels = data.lab2, xtstep = 3)

### Plot
windows()
seqdplot(seqSim2, border = NA, main = 'Pattern of Group 2 (Sample Based)')
windows()
seqdplot(FseqSim2, border = NA, main = 'Pattern of Group 2 (Frequency Based)')

###################################################################################################
#######################################  Transition over Time    ##################################
###################################################################################################
# Sample based 
seqSim1_strt<-seqdef(sim1[,1:30], alphabet = data.alphab2,
                     labels = data.lab2, xtstep = 3)
sim1_transtrt=seqtrate(seqSim1_strt)

seqSim1_mid<-seqdef(sim1[,31:60], alphabet = data.alphab2,
                    labels = data.lab2, xtstep = 3)
sim1_tranmid=seqtrate(seqSim1_mid)

seqSim1_end<-seqdef(sim1[,61:100], alphabet = data.alphab2,
                    labels = data.lab2, xtstep = 3)
sim1_tranend=seqtrate(seqSim1_end)

round(sim1_transtrt,digits = 4)
round(sim1_tranmid,digits = 4)
round(sim1_tranend,digits = 4)

# Frequency based 
FseqSim1_strt<-seqdef(Fsim1[,1:30], alphabet = data.alphab2,
                      labels = data.lab2, xtstep = 3)
Fsim1_transtrt=seqtrate(FseqSim1_strt)

FseqSim1_mid<-seqdef(Fsim1[,31:60], alphabet = data.alphab2,
                     labels = data.lab2, xtstep = 3)
Fsim1_tranmid=seqtrate(FseqSim1_mid)

FseqSim1_end<-seqdef(Fsim1[,61:100], alphabet = data.alphab2,
                     labels = data.lab2, xtstep = 3)
Fsim1_tranend=seqtrate(FseqSim1_end)

round(Fsim1_transtrt,digits =4)
round(Fsim1_tranmid,digits = 4)
round(Fsim1_tranend,digits =4)



## Group 2
# Sample based 
seqSim2_strt<-seqdef(sim2[,1:30], alphabet = data.alphab2,
                     labels = data.lab2, xtstep = 3)
sim2_transtrt=seqtrate(seqSim2_strt)

seqSim2_mid<-seqdef(sim2[,31:60], alphabet = data.alphab2,
                    labels = data.lab2, xtstep = 3)
sim2_tranmid=seqtrate(seqSim2_mid)

seqSim2_end<-seqdef(sim2[,61:100], alphabet = data.alphab2,
                    labels = data.lab2, xtstep = 3)
sim2_tranend=seqtrate(seqSim2_end)

round(sim2_transtrt,digits = 4)
round(sim2_tranmid,digits = 4)
round(sim2_tranend,digits = 4)

# Frequency based 
FseqSim2_strt<-seqdef(Fsim2[,1:30], alphabet = data.alphab2,
                      labels = data.lab2, xtstep = 3)
Fsim2_transtrt=seqtrate(FseqSim2_strt)

FseqSim2_mid<-seqdef(Fsim2[,31:60], alphabet = data.alphab2,
                     labels = data.lab2, xtstep = 3)
Fsim2_tranmid=seqtrate(FseqSim2_mid)

FseqSim2_end<-seqdef(Fsim2[,61:100], alphabet = data.alphab2,
                     labels = data.lab2, xtstep = 3)
Fsim2_tranend=seqtrate(FseqSim2_end)

round(Fsim2_transtrt,digits = 4)
round(Fsim2_tranmid,digits = 4)
round(Fsim2_tranend,digits = 4)


#####################################################################################################
#################################  OUTPUT DATA ######################################################
#####################################################################################################
saveRDS(sim1,file="C:/Users/jiz13007/OneDrive - University of Connecticut/Pattern Recognition/R script/Time varying MM/Simulation study/sim1_new.rds")
saveRDS(Fsim1,file="C:/Users/jiz13007/OneDrive - University of Connecticut/Pattern Recognition/R script/Time varying MM/Simulation study/Fsim1_new.rds")
saveRDS(sim2,file="C:/Users/jiz13007/OneDrive - University of Connecticut/Pattern Recognition/R script/Time varying MM/Simulation study/sim2_new.rds")
saveRDS(Fsim2,file="C:/Users/jiz13007/OneDrive - University of Connecticut/Pattern Recognition/R script/Time varying MM/Simulation study/Fsim2_new.rds")


# SAVE PARAMETERS

simpara<-structure(list(gamma=dfgamma, c=c, c0=c0, delta=delta, q=q))
saveRDS(Fsim2,file="C:/Users/jiz13007/OneDrive - University of Connecticut/Pattern Recognition/R script/Time varying MM/Simulation study/Simulation Parameters.rds")

