rm(list=ls()) # clear memory
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(seqHMM)
library(TraMineR)
library(ggplot2)

# load final combined clustering and sequence data
################################################################################################################
####################################          LOAD Data         ################################################
################################################################################################################

myMMM<-read.csv(file="C:/Users/jiz13007/Documents/Pattern Recognition/Book Chapter/Data/Result/Final Result/myMMM_all.csv")
dim(myMMM)

table(myMMM$sample_name)
myMMM$CLUSTER<-paste0("Cluster"," ", as.character(myMMM$X))
table(myMMM$CLUSTER)

# features
varlist<-c("AtoA",	"EtoA",	"FtoA",	"GtoA",	"HtoA",
           "AtoE",	"EtoE",
           "AtoF",	"FtoF",
           "AtoG",	"GtoG",
           "AtoH",  "HtoH")

d <- dist(myMMM[,varlist], method = "euclidean")

# # compute divisive hierarchical clustering
# hc4 <- diana(myMMM[,varlist])
# # Divise coefficient; amount of clustering structure found
# hc4$dc
# # plot dendrogram
# windows()
# pltree(hc4, cex = 0.6, hang = -1, main = "Dendrogram of diana")
# methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

# function to compute coefficient
ac <- function(x) {
  agnes(myMMM[,varlist], method = x)$ac
}

map_dbl(m, ac)

hc3 <- agnes(myMMM[,varlist], method = "ward")
windows()
pltree(hc3, cex = 0.6, hang = -1, main = "Dendrogram of agnes")

subgrp3<-cutree(as.hclust(hc3), k = 3)
# Number of members in each cluster
table(subgrp3)

subgrp4<-cutree(as.hclust(hc3), k = 4)
# Number of members in each cluster
table(subgrp4)

subgrp5<-cutree(as.hclust(hc3), k = 5)
# Number of members in each cluster
table(subgrp5)

subgrp<-ifelse(subgrp4==1,1,
               ifelse(subgrp4==2,2,
                      ifelse(subgrp4==3,3,
                             ifelse(subgrp4==4,2,0))))
subgrp2<-ifelse(subgrp4==1,1,
               ifelse(subgrp4==2,2,
                      ifelse(subgrp4==3,3,
                             ifelse(subgrp4==4,4,0))))
table(subgrp)
myMMM$ward_clu3<-subgrp
myMMM$ward_clu4<-subgrp2


#### Load sequence data
perclust=readRDS(file="C:/Users/jiz13007/Documents/Pattern Recognition/Book Chapter/Data/Result/Final Result/perclust_all.rds")
table(perclust$num_cluster)

table(perclust$sample_name)
perclust1<-perclust[perclust$ward_clu3==1,]
dim(perclust1)
perclust2<-perclust[perclust$ward_clu3==2,]
dim(perclust2)
perclust3<-perclust[perclust$ward_clu3==3,]
dim(perclust3)


rm(perclust)


saveRDS(perclust1,file="C:/Users/jiz13007/Documents/Pattern Recognition/Book Chapter/Data/Result/Final Result/threeclu_clu1.rds")
saveRDS(perclust2,file="C:/Users/jiz13007/Documents/Pattern Recognition/Book Chapter/Data/Result/Final Result/threeclu_clu2.rds")
saveRDS(perclust3,file="C:/Users/jiz13007/Documents/Pattern Recognition/Book Chapter/Data/Result/Final Result/threeclu_clu3.rds")
################################ Creat Plots  ######################################################
rm(list=ls()) # clear memory
perclust<-readRDS(file="C:/Users/jiz13007/Documents/Pattern Recognition/Book Chapter/Data/Result/Final Result/threeclu_clu3.rds")
dim(perclust)
read.sequence <- as.character(perclust[,which( colnames(perclust)=="TXTSEQUENCE" )])
str.bin <- strsplit(read.sequence, NULL)
length(str.bin)
rm(perclust)

ts.series <- sapply(str.bin, FUN = function(x){
  y = rep(0,length(x))
  y[x=="A"] = "A"
  y[x=="B"] = "A"
  y[x=="C"] = "A"
  y[x=="D"] = "A"
  y[x=="E"] = "E"
  y[x=="F"] = "F"
  y[x=="G"] = "G"
  y[x=="H"] = "H"
  y[x=="T"] = "E"
  return(y)
})

dim(ts.series)

labeluse <- c("Travel",
              "Home",
              "Mandatory",
              "Maintenance",
              "Discretionary")

ts.density <- apply(ts.series, 1, FUN = function(x){
  table( factor(x, levels = c("A","E","F","G","H"),
                labels = labeluse) )})/dim(ts.series)[2]

############################################ DENSITY STACK DATA ##############################
df_density=as.data.frame(ts.density)
stack_density<-stack(df_density)
dim(stack_density)
stack_density$TYPE=rep((c("Travel",
                          "Home",
                          "Mandatory",
                          "Maintenance",
                          "Discretionary")),1440)

stack_density$time= rep(1:1440, each=5)

windows()
ggplot(stack_density, aes(x=time, y=values, fill=TYPE)) + 
  geom_area(alpha=0.6 , size=1,colour="black")+ggtitle("Cluster 1")

### Three Clusters
saveRDS(stack_density, file ="C:/Users/jiz13007/Documents/Pattern Recognition/Book Chapter/Data/Result/Final Result/stacked_density_3clusters_clu3.rds" )
saveRDS(df_density, file ="C:/Users/jiz13007/Documents/Pattern Recognition/Book Chapter/Data/Result/Final Result/density_3clusters_clu3.rds" )


#################################################################################################################
test= perclust[perclust$sample_name=="sample_24",]
test1=test[test$CLUSTER== "Cluster 2",]

dim(test1)
table(test$CLUSTER)
read.sequence <- as.character(test1[,which( colnames(test1)=="TXTSEQUENCE" )])
str.bin <- strsplit(read.sequence, NULL)
length(str.bin)

ts.series <- sapply(str.bin, FUN = function(x){
  y = rep(0,length(x))
  y[x=="A"] = "A"
  y[x=="B"] = "A"
  y[x=="C"] = "A"
  y[x=="D"] = "A"
  y[x=="E"] = "E"
  y[x=="F"] = "F"
  y[x=="G"] = "G"
  y[x=="H"] = "H"
  y[x=="T"] = "E"
  return(y)
})

dim(ts.series)

labeluse <- c("Travel",
              "Home",
              "Mandatory",
              "Maintenance",
              "Discretionary")

ts.density <- apply(ts.series, 1, FUN = function(x){
  table( factor(x, levels = c("A","E","F","G","H"),
                labels = labeluse) )})/dim(ts.series)[2]


df_density=as.data.frame(ts.density)
stack_density<-stack(df_density)
dim(stack_density)
stack_density$TYPE=rep((c("Travel",
                          "Home",
                          "Mandatory",
                          "Maintenance",
                          "Discretionary")),1440)

stack_density$time= rep(1:1440, each=5)

windows()
ggplot(stack_density, aes(x=time, y=values, fill=TYPE)) + 
  geom_area(alpha=0.6 , size=1,colour="black")+ggtitle("1 node")
