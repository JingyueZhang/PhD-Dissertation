## This script is used for estimating time-varying mixture Markov model
## The EM algorithm used by this script is developed by Jingyue
rm(list=ls()) # clear memory
library("markovchain")
library(TraMineR) # This library is designed for analyzing social-science sequence (i.e. life course sequence). The library is required for creating sequence object
library(seqHMM) # This library is used for Mixture MM and MM
library(gtools)
library('cluster')
library(plyr)
library(HiddenMarkov)
library("optimization")
library("nloptr")
library(Rsolnp)
### Read simulated data
### Check Script for generating Markov chains based on time-varying transition probabilites "MarkovChainSimulation.R"
Fsim1<-readRDS(file="C:/Users/jiz13007/OneDrive - University of Connecticut/Pattern Recognition/R script/Time varying MM/Simulation study/Fsim1_new.rds")
Fsim2<-readRDS(file="C:/Users/jiz13007/OneDrive - University of Connecticut/Pattern Recognition/R script/Time varying MM/Simulation study/Fsim2_new.rds")

Fsim1<-as.data.frame(Fsim1)
Fsim2<-as.data.frame(Fsim2)

Fsim1$group=1
Fsim2$group=2
## Combine two groups of data
aggseq=rbind(Fsim1, Fsim2)
dim(aggseq)


######################################################################################################
###########################                Basic input            ####################################
######################################################################################################
# sequential data 
myseq=aggseq[sample(nrow(aggseq)),]
colnames(myseq)[1:(ncol(myseq)-1)]=seq(0,99)
View(myseq)
rm(aggseq)
#m:number of state
m<-3
#Sname: name of the state
Sname<-c("A","B","C")
#n:number of subjects
n<-dim(myseq)[1]
#numcomp: number of components
numcomp<-2
#T: length of the sequence
lenseq<-100

#x: exogeneous variables
x_mid<-c(rep(0,30),rep(1,30),rep(0,40))
x_end<-c(rep(0,60),rep(1,40)) # for test
exog<-do.call(rbind, Map(data.frame,
                         XMID=x_mid,
                         XEND=x_end))

dim(exog)
#p: number of exogenous variable
numexo<-dim(exog)[2]
epsilon=0.00001
######################################################################################################
#############################         INITIAL VALUES FOR PARAMETERS       ############################
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
XLIST<-names(exog)

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

gamma<-initgamma_function(c(4,8),m,numpair,numexo)
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
c_value<-function(initc,num_state,num_pair,num_component){
  c<-matrix(c(rep(initc,num_pair)),
            nrow = num_pair, ncol = num_component,byrow = TRUE,
            dimnames = list(TRANLIST,COMPLIST))
  
  for (indx in 1:num_state){
    same_indx=(indx-1)*num_state+indx
    c[same_indx,]=0}
  
  return(c)}

c<-c_value(c(-1,-2),m,numpair,numcomp)

### last state is constraint to be 0
c0<-matrix(c(rep(1,(m-1)),0,rep(1,(m-1)),0), 
           nrow=m, ncol=numcomp, byrow = FALSE,dimnames = list(INITLIST,COMPLIST))

# Set initial values for component probability (pi)
pi<-c(0.2, 0.8)

initgamma1=dfgamma
initC1=c
initc01=c0
initPi1=pi

######################################################################################################
#################################   CALCULATE INITIAL PROBABILITIES     ##############################
######################################################################################################
#Calculate initial probabilities
delta_values<-function(c0_maxtrix,num_state,num_component){
  delta<-matrix(-99, nrow = num_state, ncol = num_component, dimnames = list(INITLIST, COMPLIST))
  for (g in 1:dim(c0_maxtrix)[2]){# calculate initial probablities for each component
    for (j in 1:dim(c0_maxtrix)[1]){
      cj0g<-c0_maxtrix[j,g]
      deltaj0g<-(exp(cj0g)/(1+(sum(exp(c0_maxtrix[1:(num_state-1),g])))))
      delta[j,g]<-deltaj0g}}
  return (delta)}

######### transition probability is a funtion of gamma and c
#### Calculate transition probability at each time point
# if initial tranisiton probability value is 0, then do not change, keep as 0

######################################################################################################
#############################   CALCULATE TRANSITION PROBABILITIES       #############################
######################################################################################################
#qtjkg=exp(x_t*gamma_jk+c_jk)/1+sum(h=1:m; h is not j){exp(x_t*gamma_jh+c_jh)}
q_function<-function(gamma_matrix,c_matrix,exog_matrix,num_state,num_component,seqlength){
  # CREATE COLUMN NAMES FOR TRANSITON PROBABILITIES
  TGLIST<-c()
  for (g in 1:num_component){
    for(t in 1:(seqlength-1)){ # 1 to T
      varname<-paste("T",t,"Comp",g, sep="_")
      #print(t)
      TGLIST<-c(TGLIST,varname)}}
  
  num_pair<-num_state*num_state # number of transition pairs
  p<-dim(exog_matrix)[2]   # number of exogenous variables
  
  q<-matrix(-99, nrow = num_pair, ncol = ((seqlength-1)*num_component),byrow=TRUE,
            dimnames = list(TRANLIST,TGLIST)) #create a matrix to keep transition probability
  
  
  for (t in (1:(seqlength-1))){# no need to calcuate q_0
    exog_t<-as.vector(exog_matrix[(t+1),]) # a vector (1*p) at time point t
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
            q[strt_indx+k-1,(t+(g-1)*(seqlength-1))]=q_jkg}
          else{
            c_jkg<-c_jg[k] # single number
            gamma_jk<-temp_gamma[k] # 1 single value
            q_jkg<-exp(gamma_jk+c_jkg)/(1+sum(exp(temp_gamma[-c(j,constrain_indx)]+c_jg[-c(j,constrain_indx)])))
            q[strt_indx+k-1,(t+(g-1)*(seqlength-1))]=q_jkg}
        } # close of k loop
      }# close of j loop
    } # close of g loop
  }# close of t loop
  
  return (q)
}  # end of q function

######################################################################################################
#####################################   LOG JOINT PROBABILITIES       ################################
######################################################################################################
log_joint<-function(seqdata,delta_matrix, q_matrix, num_component,seqlen,num_state){
  state_name<-Sname
  num_pair<-num_state*num_state
  ### STEP 1: Calculate the first forward probabilities at time 0
  logjoint_0<-c()
  # calculate initial forward probabilities: t=0
  for (g in 1: num_component){ # loop through each component
    logjoint_g0<-c() # for each component, reset matrix
    for (j in 1:length(state_name)){ # loop through each state
      init_s<-state_name[j]
      delta_jg<-delta_matrix[j,g]
      if(any(delta_jg==0)){
        delta_jg=delta_jg+epsilon
      }
      logdelta_jg<-ifelse(seqdata[,1] == init_s, log(delta_jg), 0)
      logjoint_g0<-cbind(logjoint_g0,logdelta_jg)
      colnames(logjoint_g0)[ncol(logjoint_g0)] <- paste(init_s,0,"C",g,sep="_")
    } # out of j loop
    
    # sum over j to get final joint log probability
    temp_sum0<-rowSums(logjoint_g0)
    logjoint_0<-cbind(logjoint_0, temp_sum0)
    colnames(logjoint_0)[ncol(logjoint_0)]<-paste("T_0","C",g, sep="_")
  }
  
  logjoint_matrix<-c()
  logjoint_matrix<-cbind(logjoint_matrix, logjoint_0)
  
  ###############################################  
  for(t in 1:(seqlen-1)){ # from t=1 to T
    #for(t in 1:3){
    logjoint_t<-c() # for each t point reset matrix
    for(comp in 1:num_component){
      logjoint_g<-c() # for each component reset matrix
      q_tg<-q_matrix[,(comp-1)*(ncol(q_matrix)/2)+t]
      #print(colnames(q_matrix)[(comp-1)*(ncol(q_matrix)/2)+t])
      
      if(any(q_tg==0)){
        q_tg=q_tg+epsilon # avoid log(0)
      }
      
      previous_jointg<-logjoint_matrix[,(ncol(logjoint_matrix)-num_component+comp)]
      #print(colnames(logjoint_matrix)[(ncol(logjoint_matrix)-num_component+comp)])
      for(j in 1:num_state){
        current_j=state_name[j]
        q_tjg<-q_tg[((j-1)*num_state+1):(j*num_state)]
        for(k in 1:num_state){
          next_k<-state_name[k]
          q_tjkg<-q_tjg[k]
          logprob_tjkg<-ifelse((seqdata[,t]==current_j)&(seqdata[,t+1]==next_k),log(q_tjkg),0)
          logjoint_g<-cbind(logjoint_g,logprob_tjkg)
          colnames(logjoint_g)[ncol(logjoint_g)]<-paste(current_j,next_k,"T",t,"C",comp,sep="_")
        }
      } # end of k loop
      temp_sum<-rowSums(logjoint_g)+previous_jointg
      logjoint_t<-cbind(logjoint_t,(temp_sum))
      colnames(logjoint_t)[ncol(logjoint_t)]<-paste("T",t,"C",comp,sep="_")
    }
    logjoint_matrix<-cbind(logjoint_matrix,(logjoint_t))
  }
  final_prob<-logjoint_matrix[,(ncol(logjoint_matrix)-num_component+1):(ncol(logjoint_matrix))]
  
  return(structure(list(log_prob=logjoint_matrix, logjoint=final_prob)))
}



##########################################################################################################################
###########################################  Calculate posterior probabilities ###########################################
##########################################################################################################################
poster_probs<-function(seqdata,joint_prob, pi_values, lenseq, state_names, num_state,num_component){
  
  ## Apply log-sum-exp trick to aviod underflow
  ## eta = pi_g*exp(log(p_ig)-max_g(p_ig))/ sum{over g}[exp(log(p_ig)-max_g(p_ig))]
  maxprob=apply(joint_prob, 1, max) 
  #View(maxprob)
  logdiff<-c()
  for (comp in 1:num_component){
    temp_log<-joint_prob[,comp]
    diff_log<-temp_log-maxprob
    logdiff<-cbind(logdiff,diff_log)
    colnames(logdiff)[ncol(logdiff)]<-paste("log_diff","C",comp, sep="_")
  }
  
  #View(logdiff)
  expdiff<-exp(logdiff)
  for (comp in 1:num_component){
    expdiff[logdiff[,comp]<=-10,comp]=0
  }
  
  #View(expdiff)
  # Calculate posterior probabilities
  exp_matrix<-c()
  for (comp in 1:num_component){
    temp_pi<-pi_values[comp]
    temp_exp<-temp_pi*expdiff[,comp]
    exp_matrix<-cbind(exp_matrix,temp_exp)
    colnames(exp_matrix)[ncol(exp_matrix)]<-paste("exp","C",comp, sep = "_")
  }
  #View(exp_matrix)
  
  eta<-c()
  for (comp in 1:num_component){
    temp_eta<-exp_matrix[,comp]/rowSums(exp_matrix)
    eta<-cbind(eta,temp_eta)
    colnames(eta)[ncol(eta)]<-paste("eta","C",comp, sep = "_")
  }
  #View(eta)
  
  return(eta)
}


##########################################################################################################
############################################### M STEP ###################################################
##########################################################################################################

# Maximuize three terms separately
# try both numerical solution and exact solution based on Lagrangian function

#################################################################################################
#################################### Update pi ##################################################
#################################################################################################
pi_lagrangian<-function(eta_matrix,num_component){
  updated_pi<-c() 
  for (comp in 1: num_component){
    temp_pi<-sum(eta_matrix[,comp])/dim(eta_matrix)[1]
    updated_pi<-c(updated_pi,temp_pi)
  }
  return(updated_pi)
}


pi_opt<-function(pi_para){
  
  if(any(pi_para==0)){
    pi_para=pi_para+epsilon
  }
  
  sum_value<-c()
  for(comp in 1:numcomp){
    temp_pi<-pi_para[comp]
    temp_value<-sum(log(temp_pi)*exp_eta[,comp])
    sum_value<-c(sum_value,temp_value)
  }
  return(-sum(sum_value))
}

pi_equal<-function(pi_para){
  constr=sum(pi_para)
}



# updated_pi.O<-solnp(c(0.1,0.9), #starting values (random - obviously need to be positive and sum to 1)
#                     pi_opt, #function to optimise
#                     eqfun=pi_equal, #equality function
#                     eqB=1,   #the equality constraint
#                     LB=c(0,0), #lower bound for parameters i.e. greater than zero
#                     UB=c(1,1))

# updated_pi.O$par
# updated_pi.O$value
# Both solutions give the same values of pi



#################################################################################################
#################################### Update delta ###############################################
#################################################################################################
delta_lagrangian<-function(eta_matrix, state_names, num_state,num_component){
  temp_matrix<-c()
  for(comp in 1:num_component){
    temp_eta<-eta_matrix[,comp]
    for(j in 1:num_state){
      init_j<-state_names[j]
      temp_v<-ifelse(myseq[,1]==init_j,1,0)
      temp_value<-temp_eta*temp_v
      temp_matrix<-cbind(temp_matrix,temp_value)
      colnames(temp_matrix)[ncol(temp_matrix)]<-paste(init_j,"C",comp, sep="_")
    }}
  
  delta_new<-c()
  for(comp in 1:num_component){
    temp_delta<-c()
    temp_matrix2<-temp_matrix[,((comp-1)*num_state+1):(comp*num_state)]
    for (j in 1:num_state){
      init_j<-state_names[j]
      temp_value2=sum(temp_matrix2[,j])/sum(eta_matrix[,comp])
      temp_delta<-c(temp_delta,temp_value2)
    }
    delta_new<-cbind(delta_new,temp_delta)
    colnames(delta_new)[ncol(delta_new)]<-paste("component",comp,sep="_")
  }
  
  rownames(delta_new)=state_names
  return(delta_new)
}



delta_opt<-function(delta_para){
  delta_matrix<-matrix(delta_para,nrow=m, ncol=numcomp)
  sum_matrix2<-c()
  for (comp in 1:numcomp){
    sum_matrix<-c()
    temp_delta<-delta_matrix[,comp]
    if(any(temp_delta==0)){
      temp_delta=temp_delta+epsilon
    }
    temp_eta<-exp_eta[,comp]
    for(j in 1:m){
      temp_delta_j=temp_delta[j]
      init_j<-Sname[j]
      temp_value<-ifelse(myseq[,1]==init_j,log(temp_delta_j)*temp_eta,0)
      sum_matrix<-cbind(sum_matrix,temp_value)
      colnames(sum_matrix)[ncol(sum_matrix)]<-paste(init_j,"C",comp, sep="_")
    }
    sum_matrix2<-cbind(sum_matrix2,rowSums(sum_matrix))
    colnames(sum_matrix2)[ncol(sum_matrix2)]<-paste("component",comp, sep="_")
  }
  return(-sum(sum_matrix2))
}

delta_equal<-function(delta_para){
  
  delta_matrix<-matrix(delta_para,nrow=m, ncol=numcomp)
  constr=colSums(delta_matrix)
  #return(constr)
}

# delta_opt(as.vector(delta))
# delta_equal(as.vector(delta))
# LB_val=rep(0, m*numcomp)
# UB_val=rep(1, m*numcomp)
# eqB_val=rep(1,numcomp)
# 
# updated_delta.O<-solnp(as.vector(delta), #starting values (random - obviously need to be positive and sum to 1)
#                        fun=delta_opt, #function to optimise
#                        eqfun=delta_equal, #equality function
#                        eqB=eqB_val,   #the equality constraint
#                        LB=LB_val, #lower bound for parameters i.e. greater than zero
#                        UB=UB_val)

# delta_opt(as.vector(updated_delta.L))
# delta_opt(updated_delta.O$pars)
## Both Lagrangian method and opt algorithm give the same result


#################################################################################################
#################################### Update gamma and c #########################################
#################################################################################################

#opt_fun = sum over{N;G;T=1:T;j=1:M; k=1:M}{eta_ig*mu_itjkg,log(q_itjkg)}
#q_itjkg=exp(gamma_jkg*exog_ijkg+c_g)/sum{over G}{exp(gamma_jkg*exog_ijkg+c_g)})
opt_fun<-function(para){
  consRow<-c()
  for (i in 1:m){
    conval=(i-1)*m+i
    consRow<-c(consRow,conval)
  }
  para = as.matrix(para, (numpair-m),(numexo+numcomp))
  newpara <- matrix(0,nrow=numpair,ncol=(numexo+numcomp))  
  newpara[-consRow,] <- para
  expGamma=newpara[,(1:(dim(newpara)[2]-numcomp))]
  expC=newpara[,((dim(newpara)[2]-numcomp+1):dim(newpara)[2])]
  # transition probability is a function of gamma and c
  expQ<-q_function(expGamma,expC,exog,m,numcomp,lenseq)
  # sum over n, T, G
  strtindx=1
  endindx=strtindx+numpair-1
  
  sum_matrix<-c()
  for (t in 1:(lenseq-1)){
    for (comp in 1:numcomp){
      temp_q = expQ[,(((comp-1)*(ncol(expQ)/2))+t)] # 1 * numpair vector
      temp_eta = exp_eta[,comp] # N * 1 vector
      
      if(any(temp_q==0)){
        temp_q=temp_q+epsilon
      }
      
      sum_matrix2<-c()
      for (j in 1:m){
        q_jtg=temp_q[((j-1)*m+1):(j*m)]
        current_j<-Sname[j]
        for(k in 1:m){
          next_k<-Sname[k]
          q_jktg=q_jtg[k]
          temp_mu<-ifelse((myseq[,t]==current_j)&(myseq[,t+1]==next_k),1,0)
          temp_val<-temp_eta*temp_mu
          temp_val2<-temp_val*log(q_jktg)
          sum_matrix2<-cbind(sum_matrix2,temp_val2)
          colnames(sum_matrix2)[ncol(sum_matrix2)]<-paste(current_j,next_k,"T",t,"C",comp, sep="_")
        }
      }
      sum_matrix<-cbind(sum_matrix,rowSums(sum_matrix2))
      colnames(sum_matrix)[ncol(sum_matrix)]<-paste("T",t,"C",comp, sep="_")}}
  #print(exp_eta)
  return(sum(sum_matrix))
}

#####################################################################################
#!!!!!!!!!!!!!!!!!!!!!!!!  REMOVE CONSTRAINED PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
######################################################################################
# consRow<-c()
# for (i in 1:m){
#   conval=(i-1)*m+i
#   consRow<-c(consRow,conval)
# }
# initpara=cbind(dfgamma[,1:numexo],c)
# initpara=initpara[-(consRow),]
# dim(initpara)
# 
# testopt=opt_fun(initpara)
# testopt
# 
# exppar<-optim(par=as.matrix(initpara),fn=opt_fun,control=list(fnscale=-1))
# exppar$value
# exppar$par
# 
# LB=matrix(-20, nrow=nrow(initpara), ncol=ncol(initpara))
# UB=matrix(20, nrow=nrow(initpara), ncol=ncol(initpara))
# exppar_limit<-optim(par=as.matrix(initpara),
#                     fn=opt_fun,control=list(fnscale=-1),
#                     lower=LB,upper=UB, method = "L-BFGS-B")
# para=exppar_limit$par
# initpara
# exppar_limit$par
# exppar_limit$value
# 
# opt_fun(exppar$par)


#################################################################################################
###################################### Q Value Calculation ######################################
#################################################################################################
maxFun<-function(pi_new,delta_new,optval){
  
  # part 1: posterior part
  obj1<-(-pi_opt(pi_new))
  # part 2: initial probabilities part
  obj2<-(-delta_opt(as.vector(delta_new)))
  
  # part 3: generated from optimal function
  
  objval=obj1+obj2+optval
  return(objval)
}

# objval=maxFun(updated_pi.L,updated_delta.L,exppar$value)
# objval


###################################################################################################
################################### RANDOM START POINT ############################################
###################################################################################################
multipoint<-function(initgamma, initC, initC0, sd=sd_value,
                     change_initC=TRUE,change_initC0=FALSE, seed=seed_value){
  set.seed(seed)
  restart_gamma<-c()
  nelement<-(nrow(initgamma)*ncol(initgamma))
  means_list<-as.vector(as.matrix(initgamma))
  
  for(i in 1:nelement){
    temp_mean<-means_list[i]
    temp_value<-rnorm(1,mean=temp_mean,sd=sd)
    restart_gamma<-c(restart_gamma,temp_value)
  }
  
  dfregamma<-matrix(restart_gamma,nrow=nrow(initgamma), 
                    ncol=ncol(initgamma))
  colnames(dfregamma)<-colnames(initgamma)
  rownames(dfregamma)<-rownames(initgamma)
  
  if(change_initC==TRUE){
    restart_C<-c()
    nelement_C<-(nrow(initC)*ncol(initC))
    means_list_C<-as.vector(as.matrix(initC))
    
    for(i in 1:nelement_C){
      temp_mean_C<-means_list_C[i]
      temp_value_C<-rnorm(1,mean=temp_mean_C,sd=sd)
      restart_C<-c(restart_C,temp_value_C)
    }
    
    dfreC<-matrix(restart_C,nrow=nrow(initC), 
                  ncol=ncol(initC))
    colnames(dfreC)=colnames(initC)
    rownames(dfreC)=rownames(initC)
  }else{dfreC=initC}
  
  if(change_initC0==TRUE){
    restart_C0<-c()
    nelement_C0<-(nrow(initC0)*ncol(initC0))
    means_list_C0<-as.vector(as.matrix(initC0))
    
    for(i in 1:nelement_C0){
      temp_mean_C0<-means_list_C0[i]
      temp_value_C0<-rnorm(1,mean=temp_mean_C0,sd=sd)
      restart_C0<-c(restart_C0,temp_value_C0)
    }
    
    dfreC0<-matrix(restart_C0,nrow=nrow(initC0), 
                   ncol=ncol(initC0))
    colnames(dfreC0)=colnames(initC0)
    rownames(dfreC0)=rownames(initC0)
  }
  
  else{dfreC0=initC0}
  
  return(structure(list(dfregamma=dfregamma,dfreC=dfreC,dfreC0=dfreC0)))
}


#test_restart<-multipoint(initgamma1, initC1, initc01, sd=2, seed = 1234)


####################################################################################################
####################################  SINGLE EM FUNCTION ###########################################
####################################################################################################

em_single<-function(initgamma, initC, initC0,initPi, restart_iter,thred){
  
  objvalue_old=-1000000000000
  restart_status=FALSE
  
  obj_list<-c()
  identifier_list<-c()
  opt_list<-c()
  
  ## Calculate initial parameters for kicking start the estimation
  
  delta<-delta_values(initC0,m,numcomp)
  
  q<-q_function(initgamma,initC,exog,m,numcomp,lenseq)
  
  pi<-initPi
  dfgamma=initgamma
  
  c=initC
  c0=initC0
  
  consRow<-c()
  for (i in 1:m){
    conval=(i-1)*m+i
    consRow<-c(consRow,conval)
  }
  
  iter<-0
  while(restart_status == FALSE){
    ########################### Start of E Step #################
    logprob<-log_joint(myseq,delta, q, numcomp,lenseq,m)          
    print (sum(is.na(logprob$logjoint)))
    
    exp_eta<<-poster_probs(myseq,logprob$logjoint, pi, lenseq, Sname, m,numcomp)
    print (sum(is.na(exp_eta)))
    #print (exp_eta)
    #print (paste(iter, "exp_eta"))
    ########################  Start of M step ######################
    
    pi_next<-pi_lagrangian(exp_eta,numcomp)
    delta_next<-delta_lagrangian(exp_eta, Sname, m,numcomp)
    initpara=cbind(dfgamma[,1:numexo],c)
    initpara=initpara[-(consRow),]
    
    exppar<-optim(par=as.matrix(initpara),fn=opt_fun,control=list(fnscale=-1))
    objvalue=maxFun(pi_next,delta_next,exppar$value)
    ################ Check change of loglikelihood ######################
    identifier=objvalue-objvalue_old
    obj_list<-c(obj_list,objvalue)
    identifier_list<-c(identifier_list,identifier)
    opt_list<-c(opt_list,exppar$value)
    
    print(obj_list)
    print(identifier_list)
    print(opt_list)
    
    ###################### Update parameters ########################
    objvalue_old=objvalue
    temppara = as.matrix(exppar$par, (numpair-m),(numexo+numcomp))
    temppara1 <- matrix(0,nrow=numpair,ncol=(numexo+numcomp))  
    temppara1[-consRow,] <- temppara
    colnames(temppara1)=colnames(cbind(dfgamma[,1:numexo],c))
    rownames(temppara1)=rownames(cbind(dfgamma[,1:numexo],c))
    #update gamma and c
    dfgamma=temppara1[,(1:(dim(temppara1)[2]-numcomp))]
    c=temppara1[,((dim(temppara1)[2]-numcomp+1):dim(temppara1)[2])]
    # update q
    q<-q_function(dfgamma,c,exog,m,numcomp,lenseq)
    #dim(q)
    #View(q)
    # update delta
    delta=delta_next
    # update pi
    pi=pi_next
    
    print(dfgamma)
    print(c)
    print(pi)
    
    #####################################################################
    iter=iter+1
    # stop the loop if reach to the threshold or maximum number of iteration
    if(iter==restart_iter){
      restart_status <- TRUE
    }
    if (identifier<=thred){
      restart_status <- TRUE
    }
  }
  
  
  ######################## Return parameters #######################
  return (structure(list(objval = obj_list, gredient = identifier_list, opt=opt_list, iter = iter, 
                         gamma = dfgamma, c=c, pi=pi,delta=delta,final_objval=objvalue)))
}

#sink('C:/Users/jiz13007/OneDrive - University of Connecticut/Pattern Recognition/R script/Time varying MM/Simulation study/testoutput.txt')
#test_em<-em_single(initgamma1,initC1,initc01,initPi1,3,0.00001)
#sink()

####################################################################################################
####################################  EM FUNCTION #####################################################
####################################################################################################
em<-function(thre,res_thre,maxiter, initgamma, initC, initC0,initPi,num_restart, maxiter_res,sd_value,seed_value){
  objvalue_old=-1000000000000
  finish=FALSE
  
  obj_list<-c()
  identifier_list<-c()
  opt_list<-c()
  
  restart_cnt<-1
  
  
  opt_result<-em_single(initgamma,initC,initC0,initPi,maxiter_res,res_thre)
  max_obj<-opt_result$final_objval
  
  opt_initpara<-structure(list(gamma=initgamma, c=initC, c0=initC0,pi=initPi))
  opt_para<-structure(list(gamma=opt_result$gamma,c=opt_result$c,pi=opt_result$pi,delta=opt_result$delta ))
  
  initgamma_res=initgamma
  initC_res=initC
  initC0_res=initC0
  initPi_res=initPi
  
  
  multi_obj=max_obj
  while (restart_cnt<=num_restart){
    # generate new start points
    new_points<-multipoint(initgamma_res, initC_res, initC0_res, sd=sd_value,seed=seed_value)
    # update parameters
    initgamma_res<-new_points$dfregamma 
    initC_res<-new_points$dfreC
    initC0_res<-new_points$dfreC0
    # em algorithm
    current_result<-em_single(initgamma_res,initC_res,initC0_res,initPi_res,maxiter_res,res_thre)
    
    
    current_obj<-current_result$final_objval
    current_initpara<-structure(list(gamma=initgamma_res, c=initC_res, c0=initC0_res,pi=initPi_res))
    current_para<-structure(list(gamma=current_result$gamma,c=current_result$c,
                                 pi=current_result$pi,delta=current_result$delta ))
    
    
    multi_obj<-c(multi_obj,current_obj)
    
    if(current_obj>max_obj){ # if new start points result in greater loglikelood
      max_obj=current_obj
      opt_result=current_result
      opt_initpara=current_initpara
      opt_para=current_para
    }
    restart_cnt=restart_cnt+1
  }#after all the random start, gives us the opt random starts
  
  ### USE THE OPT START POINT TO CONTINUE ESTIMATION
  
  ## Calculate initial parameters for kicking start the estimation
  
  delta<-opt_para$delta
  
  q<-q_function(opt_para$gamma,opt_para$c,exog,m,numcomp,lenseq)
  
  pi<-opt_para$pi
  
  dfgamma=opt_para$gamma
  
  c=opt_para$c
  
  consRow<-c()
  for (i in 1:m){
    conval=(i-1)*m+i
    consRow<-c(consRow,conval)
  }
  
  iter<-0
  while(finish == FALSE){
    
    ########################### Start of E Step #################
    logprob<-log_joint(myseq,delta, q, numcomp,lenseq,m)          
    print (sum(is.na(logprob$logjoint)))
    
    exp_eta<<-poster_probs(myseq,logprob$logjoint, pi, lenseq, Sname, m,numcomp)
    print (sum(is.na(exp_eta)))
    
    ########################  Start of M step ######################
    
    pi_next<-pi_lagrangian(exp_eta,numcomp)
    delta_next<-delta_lagrangian(exp_eta, Sname, m,numcomp)
    
    initpara=cbind(dfgamma[,1:numexo],c)
    initpara=initpara[-(consRow),]
    
    exppar<-optim(par=as.matrix(initpara),fn=opt_fun,control=list(fnscale=-1))
    objvalue=maxFun(pi_next,delta_next,exppar$value)
    
    ################ Check change of loglikelihood ######################
    identifier=objvalue-objvalue_old
    obj_list<-c(obj_list,objvalue)
    identifier_list<-c(identifier_list,identifier)
    opt_list<-c(opt_list,exppar$value)
    
    print(obj_list)
    print(identifier_list)
    print(opt_list)
    
    ###################### Update parameters ########################
    objvalue_old=objvalue
    temppara = as.matrix(exppar$par, (numpair-m),(numexo+numcomp))
    temppara1 <- matrix(0,nrow=numpair,ncol=(numexo+numcomp))  
    temppara1[-consRow,] <- temppara
    colnames(temppara1)=colnames(cbind(dfgamma[,1:numexo],c))
    rownames(temppara1)=rownames(cbind(dfgamma[,1:numexo],c))
    #update gamma and c
    dfgamma=temppara1[,(1:(dim(temppara1)[2]-numcomp))]
    c=temppara1[,((dim(temppara1)[2]-numcomp+1):dim(temppara1)[2])]
    # update q
    q<-q_function(dfgamma,c,exog,m,numcomp,lenseq)
    #dim(q)
    #View(q)
    # update delta
    delta=delta_next
    # update pi
    pi=pi_next
    
    print(dfgamma)
    print(c)
    print(pi)
    
    #####################################################################
    iter=iter+1
    # stop the loop if reach to the threshold or maximum number of iteration
    
    if (identifier<=thre){
      finish <- TRUE
    }
    if(iter==maxiter){
      finish <- TRUE
    }
  }
  
  ################ Final Parameters ##############################
  logprob_final<-log_joint(myseq,delta, q, numcomp,lenseq,m)          
  print (sum(is.na(logprob_final$logjoint)))
  print ("Done 1")
  exp_eta_final<-poster_probs(myseq,logprob_final$logjoint, pi, lenseq, Sname, m,numcomp)
  print (sum(is.na(exp_eta_final)))
  print ("Done 2")
  cluster_final<-colnames(logprob_final$logjoint)[apply(logprob_final$logjoint,1,which.max)]
  print ("Done 3")
  
  ######################## Return parameters #######################
  return (structure(list(objval = obj_list, gredient = identifier_list, opt=opt_list, iter = iter, 
                         gamma = dfgamma, c=c, pi=pi,delta=delta,q=q, logprob_final=logprob_final, 
                         eta=exp_eta_final, cluster=cluster_final, 
                         opt_initpara=opt_initpara, opt_strtpara=opt_para,
                         strt_objval=multi_obj)))
}

myseq1=myseq
myseq=myseq1[,1:100]
############################ Model Estimation ######################
start_time <- Sys.time()
output<-em(0.000001,0.000001,30,initgamma1, initC1,initc01,initPi1,10,4,3,1234)

end_time <- Sys.time()
runtime<-end_time-start_time
print(runtime)
View(output)

saveRDS(output, file = "C:/Users/jiz13007/OneDrive - University of Connecticut/Pattern Recognition/R script/Time varying MM/Simulation study/EstimationResult_new.rds")

