## This script is used for fitting Mixture Homogeneous Markov Model
## 4 MMM will be fitted for each subsample

rm(list=ls()) # clear memory
library(TraMineR) # This library is designed for analyzing social-science sequence (i.e. life course sequence). The library is required for creating sequence object
library(seqHMM) # This library is used for Mixture MM and MM
library(plyr)

for (samnum in c(93,98,100)){
  ################################################################################################################
  #########################                   LOAD SUBSAMPLE DATA       ##########################################
  ################################################################################################################
  datafname<-paste("C:/Users/jiz13007/Documents/Pattern Recognition/NHTS sequence data/data with social demographic/subdata/subsample_",samnum,".csv")
  subdat=read.csv(file=datafname,header=T,sep=",")
  #subdat<-subdat1[1:300,]
  
  dim(subdat)
  
  
  ################################################################################################################
  #####################                    CREATE SEQUENCE OBJECT       ##########################################
  ################################################################################################################
  dataseq<- seqdecomp(subdat,which( colnames(subdat)=="TXTSEQUENCE" ))
  #head(dataseq)
  ## convert to multiple rows 
  dataseq2<-seqdecomp(dataseq, sep = "")
  
  ## label of states
  data.lab <- c("Auto", "Public Transit", "Non-motorized",
                "Other Mode", "Home", "Mandatory", "Maintenance", "Discretionary", "Initial Home")
  ## value of states
  data.alphab <-c("A","B","C","D","E","F","G","H","T")
  
  txtseq<- seqdef(dataseq2, alphabet = data.alphab,
                  labels = data.lab, xtstep = 9)
  
  ################################################################################################################
  #############################                   AGGREGATE STATES      ##########################################
  ################################################################################################################
  varlist<-names(txtseq)
  txtseq2=txtseq
  for (i in varlist){
    txtseq2[[i]] = as.factor(revalue(txtseq[[i]], c("T"="E","A"="A","B"="A","C"="A","D"="A",
                                                    "E"="E","F"="F","G"="G","H"="H")))}
  
  data.lab2 <- c("Auto", "Home", "Mandatory", "Maintenance", "Discretionary")
  ## value of states
  data.alphab2 <-c("A","E","F","G","H")
  aggseq<-seqdef(txtseq2, alphabet = data.alphab2,
                 labels = data.lab2, xtstep = 5)
  
  
  # windows()
  # seqdplot(txtseq, border = NA, main = 'Simulated Sequence State Distribution-Nine States')
  # 
  # windows()
  # seqdplot(aggseq, border = NA, main = 'Observed Sequence State Distribution-Nine States')
  
  ################################################################################################################
  ###################                    FIT MIXTURE MARKOV MODEL       ##########################################
  ################################################################################################################
  
  set.seed(1234)
  inimx<-c(0.08,0.89,	0.01,	0.01,	0.01)
  tramx<-seqtrate(aggseq)
  #tramx<-matrix(c(0.9445,	0,	0,	0,	0.0192,	0.0112,	0.0197,	0.0054,	0,
  #0,	0.9788,	0,	0,	0.0082,	0.0067,	0.0048,	0.0015,	0,
  #0,	0,	0.9349,	0,	0.0224,	0.0101,	0.0167,	0.0159,	0,
  #0,	0,	0,	0.9787,	0.0062,	0.0062,	0.0055,	0.0034,	0,
  #0.0001,	0.0001,	0.0001,	0.0001,	0.9996,	0,	0,	0,	0,
  #0.0001,	0.0001,	0.0001,	0.0001,	0,	0.9996,	0,	0,	0,
  #0.0001,	0.0001,	0.0001,	0.0001,	0,	0,	0.9996,	0,	0,
  #0.0001,	0.0001,	0.0001,	0.0001,	0,	0,	0,	0.9996,	0,
  #0.0001,	0.0001,	0.0001,	0.0001,	0,	0,	0,	0,	0.9996), nrow=9, ncol=9, byrow=T)
  
  
  table_BIC=c()   ### keep BIC value
  table_mmmfit<-list() ### keep model fitting result
  strttime_mmm<- proc.time()
  ### LOOP for model component equals to 2,3,4,5##
  for (i in 2:5){
    # setup model
    tryCatch({
      mmm_data<- build_mmm(observations = aggseq,n_clusters = i,
                           transition_probs=tramx, initial_probs=inimx,data=aggseq)
      
      # fit model
      mmm_fit_em <- fit_model(mmm_data, local_step = TRUE,
                              control_em = list(restart = list(times = 50)),threads = 6)
      
      # Obtain BIC values
      BIC=BIC(mmm_fit_em$model)
      table_BIC<-c(table_BIC,BIC)
      #table_BIC<-c(BIC(sample_1_cluster_2$model), BIC(sample_1_cluster_3$model), BIC(sample_1_cluster_4$model))
      
      # rename the model result
      table_mmmfit[[i]]<-mmm_fit_em},
      error=function(e){})
  }
  table_BIC
  
  ################################################################################################################
  #################################              CONVERT MODEL FIT      ##########################################
  ################################################################################################################
  #table_mmmfit<-c()
  #table_mmmfit<-c(table_mmmfit,cluster_2)
  #table_mmmfit<-list(table_mmmfit,cluster_3)
  #table_BIC<-c()
  #table_BIC<-c(table_BIC,BIC(cluster_2$model))
  #table_BIC<-c(table_BIC,BIC(cluster_3$model))
  # Minimum value of BIC
  compid<-which.min(table_BIC)
  # final model result
  mmmfit<-table_mmmfit[[(compid+1)]]
  
  ## attach cluster number back to sample data
  clu=summary(mmmfit$model)$most_probable_cluster
  subdat[[paste0("CLUSTER",sep="")]] <- clu
  subdat$sample_name=paste("sample",samnum,sep="_")
  
  # obtain transition probabilities 
  tran_mmm=mmmfit$model$transition_probs
  # convert transition probabbilities to dataframe (i.e., first row is the
  # transition probabilities of cluster 1)
  df_tran_mmm <- data.frame(matrix(unlist(tran_mmm), nrow=(compid+1), byrow=T))
  # change variable name
  varnew<-c("AtoA",	"EtoA",	"FtoA",	"GtoA",	"HtoA",	
             "AtoE",  "EtoE",	"FtoE",	"GtoE",	"HtoE",
             "AtoF",	"EtoF",	"FtoF",	"GtoF",	"HtoF",
             "AtoG",	"EtoG",	"FtoG",	"GtoG",	"HtoG",
             "AtoH",	"EtoH",	"FtoH",	"GtoH",	"HtoH")
  
  colnames(df_tran_mmm) <- varnew
  df_tran_mmm$num_cluster=(compid+1)
  # number of persons within each cluster
  df_tran_mmm$cntper_clu=-99
  for(comp in 1:(compid+1)){
    cntclu<-as.vector(table(clu))
    print (cntclu[comp])
    df_tran_mmm[comp,"cntper_clu"]=cntclu[comp]}
  df_tran_mmm$BIC=table_BIC[compid]
  
  # sample name
  df_tran_mmm$sample_name=paste("sample",samnum,sep="_")
  #df_tran_mmm$sample_name=paste("sample",1,sep="_")
  
  # obtain initial probabilities
  init_mmm=mmmfit$model$initial_probs
  df_init_mmm <- data.frame(matrix(unlist(init_mmm), nrow=(compid+1), byrow=T))
  init_varlist<-c("InitA","InitE","InitF","InitG","InitH")
  colnames(df_init_mmm) <- init_varlist
  df_init_mmm$sample_name=paste("sample",samnum,sep="_")
  #df_init_mmm$sample_name=paste("sample",1,sep="_")
  
  ################################################################################################################
  #################################                   WRITE DATA        ##########################################
  ################################################################################################################
  #1.write sample data with clustering result
  #2.write transition probability dataframe
  #3.write initial probability dataframe
  #4.write BIC table
  #5.write model fitting result of the final model
  #6.write model fitting result of all models
  path="C:/Users/jiz13007/Documents/Pattern Recognition/Result MMM/Aggregate/"
  
  #1.write sample data with clustering result
  samname<-paste("agg_sample_",samnum,".csv",sep="")
  filepath1<-paste(path,"subdata result/",samname, sep = "")
  write.csv(subdat,file=filepath1,row.names=T)
  
  #2.write transition probability dataframe
  tranname<-paste("agg_sample_",samnum,"_tran",".csv", sep = "")
  filepath2<-paste(path,"transition result/",tranname, sep = "")
  write.csv(df_tran_mmm,file=filepath2,row.names=T)
  
  #3.write initial probability dataframe
  initname<-paste("agg_sample_",samnum,"_init",".csv", sep = "")
  filepath3<-paste(path,"initial result/",initname, sep = "")
  write.csv(df_init_mmm,file=filepath3,row.names=T)
  
  #4.write BIC table
  BICname<-paste("agg_sample_",samnum,"_BIC",".csv", sep = "")
  filepath4<-paste(path,"all model fit result/",BICname, sep = "")
  write.csv(table_BIC,file=filepath4,row.names=T)
  
  #5.write model fitting result of the final model
  fmodelname<-paste("agg_sample_",samnum,"_final_model", sep = "")
  filepath5<-paste(path,"final model fit result/",fmodelname,".rds",sep = "")
  saveRDS(mmmfit,file=filepath5)
  
  #6.write model fitting result of all models
  allmodelname<-paste("agg_sample_",samnum,"_all_model", sep = "")
  filepath6<-paste(path,"all model fit result/",allmodelname,".rds",sep = "")
  saveRDS(table_mmmfit,file=filepath6)
  
  ################################################################################################################
  ####################             PRINT RUN TIME AND SAMPLE NUMBER     ##########################################
  ################################################################################################################
  print(samnum)
  print (proc.time()-strttime_mmm)
  rm(list=ls()) # clear memory
}
################################################################################################################
######################                    END OF SCRIPT     ##########################################
################################################################################################################
