######### COUNT EVENTS FOR EACH CLUSTER ################

rm(list=ls()) # clear memory
library(ggplot2)
library(RColorBrewer)
# library(seqHMM)
# library(TraMineR)
perclust<-readRDS(file="C:/Users/jiz13007/Documents/Pattern Recognition/Result MMM/Aggregate/Final Result/threeclu_clu3.rds")
dim(perclust)

################################ Creat Plots (three states)  ##################################################
read.sequence <- as.character(perclust[,which( colnames(perclust)=="TXTSEQUENCE" )])
str.bin <- strsplit(read.sequence, NULL)
str.bin1<-str.bin[1:52645]
str.bin1<-str.bin[52646:105289]
str.bin1<-str.bin[105290:157933]
str.bin1<-str.bin[157934:length(str.bin)]
sum(length(str.bin1),length(str.bin2),length(str.bin3),length(str.bin4))
rm(perclust)
rm(read.sequence)
rm(str.bin)

ts.series1 <- sapply(str.bin1, FUN = function(x){
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

dim(ts.series1)

labeluse <- c("Travel",
              "Home",
              "Mandatory",
              "Maintenance",
              "Discretionary")

ts.count1 <- apply(ts.series1, 1, FUN = function(x){
  table( factor(x, levels = c("A","E","F","G","H"),
                labels = labeluse) )})

sum(ts.count1[,1440])


saveRDS(ts.count1,file="C:/Users/jiz13007/Documents/Pattern Recognition/Result MMM/Aggregate/Final Result/Count_Event_Cluster3_part4.rds")

ts.count1<-readRDS(file="C:/Users/jiz13007/Documents/Pattern Recognition/Result MMM/Aggregate/Final Result/Count_Event_Cluster3_part1.rds")
ts.count2<-readRDS(file="C:/Users/jiz13007/Documents/Pattern Recognition/Result MMM/Aggregate/Final Result/Count_Event_Cluster3_part2.rds")
ts.count3<-readRDS(file="C:/Users/jiz13007/Documents/Pattern Recognition/Result MMM/Aggregate/Final Result/Count_Event_Cluster3_part3.rds")
ts.count4<-readRDS(file="C:/Users/jiz13007/Documents/Pattern Recognition/Result MMM/Aggregate/Final Result/Count_Event_Cluster3_part4.rds")
ts.count=ts.count1+ts.count2+ts.count3+ts.count4

sum(ts.count[,1440])
saveRDS(ts.count1,file="C:/Users/jiz13007/Documents/Pattern Recognition/Result MMM/Aggregate/Final Result/Count_Event_Cluster3.rds")


##### Calculate percentage of time spent for each event
Rsumcnt<-rowSums(ts.count)
Rsumper<-Rsumcnt/sum(Rsumcnt)
print (paste(sum(Rsumcnt), 210579*1440, sep=":"))
#print (paste(sum(Rsumcnt), dim(perclust)[1]*1440, sep=":"))
round(Rsumper, digits = 4)

rm(list=ls()) # clear memory
