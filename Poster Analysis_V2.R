## THIS SCRIPT IS FOR POST ANALYSIS OF TRANSPORTATION PROJECT - TRAVEL BEHAVIOR ANALYSIS

# rm(list=ls()) # clear memory
library(plyr)
require(ggplot2)
library(survey)
library("MASS")
library(descr)
library("RColorBrewer")
require(gridExtra)
library(cowplot)
# load final clustering and sequence data
perclust<-readRDS(file="C:/Users/jiz13007/Documents/Pattern Recognition/Result MMM/Aggregate/Final Result/perclust_all.rds")
dim(perclust)

######## CHANGE THIS TO YOUR LABEL ################################
perclust$final_3clust<-ifelse(perclust$ward_clu3==1,"c1-night discretionary",
                              ifelse(perclust$ward_clu3==2,"c2-home and work",
                                     ifelse(perclust$ward_clu3==3,"c3-in home",
                                            0)))

################################################################################################################
##################################      CREATE GENERATION VARIABLE    ##########################################
################################################################################################################
perclust$GI <- ifelse((perclust$WAVE == 1)&(perclust$R_AGE >= 66)&(perclust$R_AGE <= 89),1, 
                      ifelse((perclust$WAVE == 2)&(perclust$R_AGE >= 71)&(perclust$R_AGE <= 94),1,
                             ifelse((perclust$WAVE == 3)&(perclust$R_AGE >= 77)&(perclust$R_AGE <= 100),1,
                                    ifelse((perclust$WAVE == 4)&(perclust$R_AGE >= 85)&(perclust$R_AGE <= 108),1,
                                           ifelse((perclust$WAVE == 5)&(perclust$R_AGE >= 93)&(perclust$R_AGE <= 116),1,0)))))
table(perclust$GI)
table(subset(perclust$R_AGE, perclust$GI==1))
table(subset(perclust$R_AGE, (perclust$GI==1)&(perclust$WAVE==1)))


perclust$SG <- ifelse((perclust$WAVE == 1)&(perclust$R_AGE >= 47)&(perclust$R_AGE <= 65),1, 
                      ifelse((perclust$WAVE == 2)&(perclust$R_AGE >= 52)&(perclust$R_AGE <= 70),1,
                             ifelse((perclust$WAVE == 3)&(perclust$R_AGE >= 58)&(perclust$R_AGE <= 76),1,
                                    ifelse((perclust$WAVE == 4)&(perclust$R_AGE >= 66)&(perclust$R_AGE <= 84),1,
                                           ifelse((perclust$WAVE == 5)&(perclust$R_AGE >= 74)&(perclust$R_AGE <= 92),1,0)))))

table(perclust$SG)
table(subset(perclust$R_AGE, perclust$SG==1))
table(subset(perclust$R_AGE, (perclust$SG==1)&(perclust$WAVE==4)))


perclust$BB <- ifelse((perclust$WAVE == 1)&(perclust$R_AGE >= 26)&(perclust$R_AGE <= 46),1, 
                      ifelse((perclust$WAVE == 2)&(perclust$R_AGE >= 31)&(perclust$R_AGE <= 51),1,
                             ifelse((perclust$WAVE == 3)&(perclust$R_AGE >= 37)&(perclust$R_AGE <= 57),1,
                                    ifelse((perclust$WAVE == 4)&(perclust$R_AGE >= 45)&(perclust$R_AGE <= 65),1,
                                           ifelse((perclust$WAVE == 5)&(perclust$R_AGE >= 53)&(perclust$R_AGE <= 73),1,0)))))

table(perclust$BB)
table(subset(perclust$R_AGE, perclust$BB==1))
table(subset(perclust$R_AGE, (perclust$BB==1)&(perclust$WAVE==1)))


perclust$GX <- ifelse((perclust$WAVE == 1)&(perclust$R_AGE >= 9)&(perclust$R_AGE <= 25),1, 
                      ifelse((perclust$WAVE == 2)&(perclust$R_AGE >= 14)&(perclust$R_AGE <= 30),1,
                             ifelse((perclust$WAVE == 3)&(perclust$R_AGE >= 20)&(perclust$R_AGE <= 36),1,
                                    ifelse((perclust$WAVE == 4)&(perclust$R_AGE >= 28)&(perclust$R_AGE <= 44),1,
                                           ifelse((perclust$WAVE == 5)&(perclust$R_AGE >= 36)&(perclust$R_AGE <= 52),1,0)))))

table(perclust$GX)
table(subset(perclust$R_AGE, perclust$GX==1))
table(subset(perclust$R_AGE, (perclust$GX==1)&(perclust$WAVE==1)))



perclust$MN <- ifelse((perclust$WAVE == 1)&(perclust$R_AGE >= 0)&(perclust$R_AGE <= 8),1, 
                      ifelse((perclust$WAVE == 2)&(perclust$R_AGE >= 0)&(perclust$R_AGE <= 13),1,
                             ifelse((perclust$WAVE == 3)&(perclust$R_AGE >= 1)&(perclust$R_AGE <= 19),1,
                                    ifelse((perclust$WAVE == 4)&(perclust$R_AGE >= 9)&(perclust$R_AGE <= 27),1,
                                           ifelse((perclust$WAVE == 5)&(perclust$R_AGE >= 17)&(perclust$R_AGE <= 35),1,0)))))

table(perclust$MN)
table(subset(perclust$R_AGE, perclust$MN==1))
table(subset(perclust$R_AGE, (perclust$MN==1)&(perclust$WAVE==3)))

sum(perclust$GI,perclust$SG, perclust$BB, perclust$GX, perclust$MN)
for(gene1 in list("GI","SG","BB","GX","MN")){
  for(gene2 in list("GI","SG","BB","GX","MN")){
    if(gene1!=gene2){
      print (paste(gene1,gene2, sep = ","))
      print (sum((perclust[[gene1]]==1)&(perclust[[gene2]]==1)))
    }
  }
}

perclust$GENERATION<-ifelse(perclust$GI==1,1,
                            ifelse(perclust$SG==1,2,
                                   ifelse(perclust$BB==1,3,
                                          ifelse(perclust$GX==1,4,
                                                 ifelse(perclust$MN==1,5,0)))))

table(perclust$GENERATION)
table(subset(perclust$GENERATION, perclust$GI==1))

################################################################################################################
##########################################     GENERATION OVERALL    ##########################################
################################################################################################################
numclu=3
cluname<-c("c1-in home", "c2-night discretionary","c3-home and work")
wavename<-c("1990","1995","2001","2009","2017")
perclust$R_WTPERFIN=perclust$WTPERFIN ### WTPERFIN is the weight
perclust$R_WTPERFIN[is.na(perclust$R_WTPERFIN)]=1 # Replace the NAN value of weight to 1
sum(is.na(perclust$R_WTPERFIN))
sum(is.na(perclust$WTPERFIN))
table(perclust$ward_clu3)


### THIS PART IS THE OVERALL GENERATION DISTRIBUTION FOR EACH CLUSTER
### Each bar represents the proportion of generation fall into a particular cluster
### i.e., 80% of Baby Boomers fall into cluster 1

wavewt<-xtabs(perclust$R_WTPERFIN ~ perclust$GENERATION + perclust$final_3clust) 

wavewt_matrix<-matrix(NA, nrow =dim(wavewt)[1] , ncol =dim(wavewt)[2] , 
                      dimnames = list(c("GI","Silence Generation","Baby Boomer","Generation X","Millennial")
                                      ,cluname))

for (i in 1:dim(wavewt_matrix)[1]){
  wavewt_matrix[i,]=wavewt[i,]/margin.table(wavewt, 1)[i]
}
wavewt_matrix
sum(wavewt_matrix[1,])




barplot(wavewt_matrix,ylim = c(0,1),
        xlab="cluster", col=brewer.pal(n = dim(wavewt_matrix)[1], name = "Set3"),
        ylab = "proportion",
        legend.text  = c("GI","Silence Generation","Baby Boomer","Generation X","Millennial"),
        args.legend = list(x = "top"), beside=TRUE)




###################################################################################################################
#######################################   GENERATION CHANGE OVER TIME      ############################################
###################################################################################################################
## THIS PART IS TO SHOW THE CHANGE OF PROPORTION OF PERSONS FALL INTO EACH CLUSTER OVER LAST THREE DECADES
### FOR EACH GENERATION
### EACH PLOT REPRESENTS THE CHANGE OF EACH GENERATION
wave1<-perclust[perclust$WAVE==1,]
wave2<-perclust[perclust$WAVE==2,]
wave3<-perclust[perclust$WAVE==3,]
wave4<-perclust[perclust$WAVE==4,]
wave5<-perclust[perclust$WAVE==5,]

Rgenwt1<-xtabs(wave1$R_WTPERFIN ~ wave1$GENERATION + wave1$final_3clust)  # weighted cross table

Rgenwt1
Rgenwt1<-rbind(Rgenwt1, rep(0,dim(Rgenwt1)[2]))
Rgenwt1

Rgen_matrix1<-matrix(NA, nrow = dim(Rgenwt1)[1], ncol = dim(Rgenwt1)[2], 
                     dimnames = list(c("GI","Silence Generation","Baby Boomer","Generation X","Millennial"),
                                     cluname))

for (i in 1:dim(Rgenwt1)[1]){
  Rgen_matrix1[i,]=Rgenwt1[i,]/margin.table(Rgenwt1, 1)[i]
}
Rgen_matrix1
sum(Rgen_matrix1[1,])

Rgenwt2<-xtabs(wave2$R_WTPERFIN ~ wave2$GENERATION + wave2$final_3clust)  # weighted cross table
Rgenwt2
Rgenwt2<-rbind(Rgenwt2,rep(0,dim(Rgenwt2)[2]))

Rgen_matrix2<-matrix(NA, nrow = dim(Rgenwt2)[1], ncol = dim(Rgenwt2)[2], 
                     dimnames = list(c("GI","Silence Generation","Baby Boomer","Generation X","Millennial"),
                                     cluname))

for (i in 1:dim(Rgenwt2)[1]){
  Rgen_matrix2[i,]=Rgenwt2[i,]/margin.table(Rgenwt2, 1)[i]
}
Rgen_matrix2
sum(Rgen_matrix2[1,])

Rgenwt3<-xtabs(wave3$R_WTPERFIN ~ wave3$GENERATION + wave3$final_3clust)  # weighted cross table
Rgenwt3

Rgen_matrix3<-matrix(NA, nrow = dim(Rgenwt3)[1], ncol = dim(Rgenwt3)[2], 
                     dimnames = list(c("GI","Silence Generation","Baby Boomer","Generation X","Millennial"),
                                     cluname))

for (i in 1:dim(Rgenwt3)[1]){
  Rgen_matrix3[i,]=Rgenwt3[i,]/margin.table(Rgenwt3, 1)[i]
}
Rgen_matrix3
sum(Rgen_matrix3[1,])

Rgenwt4<-xtabs(wave4$R_WTPERFIN ~ wave4$GENERATION + wave4$final_3clust)  # weighted cross table
Rgenwt4

Rgen_matrix4<-matrix(NA, nrow = dim(Rgenwt4)[1], ncol = dim(Rgenwt4)[2] 
                     ,dimnames = list(c("GI","Silence Generation","Baby Boomer","Generation X","Millennial"),
                                      cluname))


for (i in 1:dim(Rgenwt4)[1]){
  Rgen_matrix4[i,]=Rgenwt4[i,]/margin.table(Rgenwt4, 1)[i]
}
Rgen_matrix4
sum(Rgen_matrix4[1,])

Rgenwt5<-xtabs(wave5$R_WTPERFIN ~ wave5$GENERATION + wave5$final_3clust)  # weighted cross table
Rgenwt5
Rgenwt5<-rbind(rep(0,dim(Rgenwt5)[2]), Rgenwt5)

Rgen_matrix5<-matrix(NA, nrow = dim(Rgenwt5)[1], ncol = dim(Rgenwt5)[2]
                     ,dimnames = list(c("GI","Silence Generation","Baby Boomer","Generation X","Millennial"),
                                      cluname))


for (i in 1:dim(Rgenwt5)[1]){
  Rgen_matrix5[i,]=Rgenwt5[i,]/margin.table(Rgenwt5, 1)[i]
}
Rgen_matrix5
sum(Rgen_matrix5[2,])

####################################              Create Plots    #############################################

#### generation 1 - GI#####

ComGI<-rbind(Rgen_matrix1[1,],Rgen_matrix2[1,],Rgen_matrix3[1,],Rgen_matrix4[1,],Rgen_matrix5[1,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp1<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ ylim(0, 1)

###### 2 - Silence Generation ################
ComGI<-rbind(Rgen_matrix1[2,],Rgen_matrix2[2,],Rgen_matrix3[2,],Rgen_matrix4[2,],Rgen_matrix5[2,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp2<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ ylim(0, 1)

###### 3 - Baby Boomers ################
ComGI<-rbind(Rgen_matrix1[3,],Rgen_matrix2[3,],Rgen_matrix3[3,],Rgen_matrix4[3,],Rgen_matrix5[3,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp3<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ ylim(0, 1)


###### 4 - Generation X ################
ComGI<-rbind(Rgen_matrix1[4,],Rgen_matrix2[4,],Rgen_matrix3[4,],Rgen_matrix4[4,],Rgen_matrix5[4,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp4<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ ylim(0, 1)

###### 5 - Millennials ################
ComGI<-rbind(Rgen_matrix1[5,],Rgen_matrix2[5,],Rgen_matrix3[5,],Rgen_matrix4[5,],Rgen_matrix5[5,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp5<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ ylim(0, 1)


#### COMBINE ALL PLOTS #######################
Rlegend <- get_legend(Rp1)
Rprow <- plot_grid( Rp1 + theme(legend.position="none"),
                    Rp2 + theme(legend.position="none"),
                    Rp3 + theme(legend.position="none"),
                    Rp4 + theme(legend.position="none"),
                    Rp5 + theme(legend.position="none"),
                    Rlegend,
                    align = 'vh',
                    labels = c("         GI",
                               "Silent Generation",
                               "Baby Boomers",
                               "Generation X",
                               "   Millennials"),
                    label_size = 14,
                    hjust = -1,
                    vjust = 1,
                    nrow = 3,
                    ncol = 2,
                    scale = 0.9
)
Rprow


################################################################################################################
######################################     WORKER    ##########################################
################################################################################################################
table(perclust$WORKER)
perclust$WORKER2=perclust$WORKER
perclust$WORKER2[perclust$WORKER==-9]=NA
sum(is.na(perclust$WORKER2))

wavewt<-xtabs(perclust$R_WTPERFIN ~ perclust$WORKER2 + perclust$final_3clust) 

wavewt_matrix<-matrix(NA, nrow =dim(wavewt)[1] , ncol =dim(wavewt)[2] , 
                      dimnames = list(c("Worker","Non-worker")
                                      ,cluname))

for (i in 1:2){
  wavewt_matrix[i,]=wavewt[i,]/margin.table(wavewt, 1)[i]
}
wavewt_matrix
sum(wavewt_matrix[1,])

### THE BAR CHART FOR OVERALL WORKER DISTRIBUTION BY CLUSTER

barplot(wavewt_matrix,ylim = c(0,1),
        #main="Percent of Age by Cluster (Unweighted)",
        xlab="cluster", col=c("red","blue"),
        ylab = "proportion",
        legend.text  = c("Worker","Non-worker"),
        args.legend = list(x = "topleft"), beside=TRUE)

wave1<-perclust[perclust$WAVE==1,]
wave2<-perclust[perclust$WAVE==2,]
wave3<-perclust[perclust$WAVE==3,]
wave4<-perclust[perclust$WAVE==4,]
wave5<-perclust[perclust$WAVE==5,]


###################################################################################################
##################################       WORKER CHANGE OVER TIME           ###############################
###################################################################################################
cluname<-c("c1-in home", "c2-night discretionary","c3-home and work")
numclu<-3

Ragewt1<-xtabs(wave1$R_WTPERFIN ~ wave1$WORKER2 + wave1$final_3clust)  # weighted cross table

Ragewt1

Rage_matrix1<-matrix(NA, nrow = dim(Ragewt1)[1], ncol = dim(Ragewt1)[2], 
                     dimnames = list( c("Worker","Non-worker"),
                                      cluname))

for (i in 1:dim(Ragewt1)[1]){
  Rage_matrix1[i ,]=Ragewt1[i ,]/margin.table(Ragewt1, 1)[i]
}
Rage_matrix1
sum(Rage_matrix1[1,])

Ragewt2<-xtabs(wave2$R_WTPERFIN ~ wave2$WORKER2 + wave2$final_3clust)  # weighted cross table
Ragewt2

Rage_matrix2<-matrix(NA, nrow = dim(Ragewt2)[1], ncol = dim(Ragewt2)[2], 
                     dimnames = list( c("Worker","Non-worker"),
                                      cluname))

for (i in 1:dim(Ragewt2)[1]){
  Rage_matrix2[i ,]=Ragewt2[i ,]/margin.table(Ragewt2, 1)[i]
}
Rage_matrix2
sum(Rage_matrix2[1,])

Ragewt3<-xtabs(wave3$R_WTPERFIN ~ wave3$WORKER2 + wave3$final_3clust)  # weighted cross table
Ragewt3

Rage_matrix3<-matrix(NA, nrow = dim(Ragewt3)[1], ncol = dim(Ragewt3)[2], 
                     dimnames = list( c("Worker","Non-worker"),
                                      cluname))

for (i in 1:dim(Ragewt3)[1]){
  Rage_matrix3[i ,]=Ragewt3[i ,]/margin.table(Ragewt3, 1)[i]
}
Rage_matrix3
sum(Rage_matrix3[1,])

Ragewt4<-xtabs(wave4$R_WTPERFIN ~ wave4$WORKER2+ wave4$final_3clust)  # weighted cross table


Rage_matrix4<-matrix(NA, nrow = dim(Ragewt4)[1], ncol = dim(Ragewt4)[2] 
                     ,dimnames = list( c("Worker","Non-worker"),
                                       cluname))


for (i in 1:dim(Ragewt4)[1]){
  Rage_matrix4[i ,]=Ragewt4[i ,]/margin.table(Ragewt4, 1)[i]
}
Rage_matrix4
sum(Rage_matrix4[1,])

Ragewt5<-xtabs(wave5$R_WTPERFIN ~ wave5$WORKER2 + wave5$final_3clust)  # weighted cross table
Ragewt5


Rage_matrix5<-matrix(NA, nrow = dim(Ragewt5)[1], ncol = dim(Ragewt5)[2]
                     ,dimnames = list( c("Worker","Non-worker"),
                                       cluname))


for (i in 1:dim(Ragewt5)[1]){
  Rage_matrix5[i ,]=Ragewt5[i ,]/margin.table(Ragewt5, 1)[i]
}
Rage_matrix5
sum(Rage_matrix5[1,])

####### 1 - WORKER ###############
ComGI<-rbind(Rage_matrix1[1,],Rage_matrix2[1,],Rage_matrix3[1,],Rage_matrix4[1,],Rage_matrix5[1,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp1<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ylim(0,1)

###### 2 - NON-WORKER ################
ComGI<-rbind(Rage_matrix1[2,],Rage_matrix2[2,],Rage_matrix3[2,],Rage_matrix4[2,],Rage_matrix5[2,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp2<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ylim(0,1)

######## COMBINED PLOT #############
Rlegend <- get_legend(Rp1)
Rprow <- plot_grid( Rp1 + theme(legend.position="none"),
                    Rp2 + theme(legend.position="none"),
                    Rlegend,
                    rel_widths = c(1,1,0.3),
                    align = 'vh',
                    labels = c("Worker","Non-workers"),
                    label_size = 14,
                    hjust = -1,
                    vjust = 1,
                    nrow = 1,
                    ncol = 3,
                    scale = 0.9
)
Rprow

################################################################################################################
######################################     GENDER    ##########################################
################################################################################################################
table(perclust$R_SEX)
perclust$R_SEX2=perclust$R_SEX
perclust$R_SEX2[(perclust$R_SEX==-7)|(perclust$R_SEX==-8)]=NA
sum(is.na(perclust$R_SEX2))

wavewt<-xtabs(perclust$R_WTPERFIN ~ perclust$R_SEX2 + perclust$final_3clust) 

wavewt_matrix<-matrix(NA, nrow =dim(wavewt)[1] , ncol =dim(wavewt)[2] , 
                      dimnames = list(c("Male","Female")
                                      ,cluname))

for (i in 1:2){
  wavewt_matrix[i,]=wavewt[i,]/margin.table(wavewt, 1)[i]
}
wavewt_matrix
sum(wavewt_matrix[1,])


## OVERALL GENDER DISTRIBUTION BY CLUSTER

barplot(wavewt_matrix,ylim = c(0,1),
        #main="Percent of Age by Cluster (Unweighted)",
        xlab="cluster", col=c("red","blue"),
        ylab="proportion",
        legend.text  = c("Male","Female"),
        args.legend = list(x = "topleft"), beside=TRUE)

###################################################################################################
##################################      GENDER CHANGE OVER TIME           #########################
###################################################################################################
wave1<-perclust[perclust$WAVE==1,]
wave2<-perclust[perclust$WAVE==2,]
wave3<-perclust[perclust$WAVE==3,]
wave4<-perclust[perclust$WAVE==4,]
wave5<-perclust[perclust$WAVE==5,]

numclu<-3

Ragewt1<-xtabs(wave1$R_WTPERFIN ~ wave1$R_SEX2 + wave1$final_3clust)  # weighted cross table

Ragewt1

Rage_matrix1<-matrix(NA, nrow = dim(Ragewt1)[1], ncol = dim(Ragewt1)[2], 
                     dimnames = list( c("Male","Female"),
                                      cluname))

for (i in 1:dim(Ragewt1)[1]){
  Rage_matrix1[i,]=Ragewt1[i,]/margin.table(Ragewt1, 1)[i]
}
Rage_matrix1
sum(Rage_matrix1[1,])

Ragewt2<-xtabs(wave2$R_WTPERFIN ~ wave2$R_SEX2 + wave2$final_3clust)  # weighted cross table
Ragewt2

Rage_matrix2<-matrix(NA, nrow = dim(Ragewt2)[1], ncol = dim(Ragewt2)[2], 
                     dimnames = list( c("Male","Female"),
                                      cluname))

for (i in 1:dim(Ragewt2)[1]){
  Rage_matrix2[i,]=Ragewt2[i,]/margin.table(Ragewt2, 1)[i]
}
Rage_matrix2
sum(Rage_matrix2[1,])

Ragewt3<-xtabs(wave3$R_WTPERFIN ~ wave3$R_SEX2 + wave3$final_3clust)  # weighted cross table
Ragewt3

Rage_matrix3<-matrix(NA, nrow = dim(Ragewt3)[1], ncol = dim(Ragewt3)[2], 
                     dimnames = list( c("Male","Female"),
                                      cluname))

for (i in 1:dim(Ragewt3)[1]){
  Rage_matrix3[i,]=Ragewt3[i,]/margin.table(Ragewt3, 1)[i]
}
Rage_matrix3
sum(Rage_matrix3[1,])

Ragewt4<-xtabs(wave4$R_WTPERFIN ~ wave4$R_SEX2+ wave4$final_3clust)  # weighted cross table


Rage_matrix4<-matrix(NA, nrow = dim(Ragewt4)[1], ncol = dim(Ragewt4)[2] 
                     ,dimnames = list( c("Male","Female"),
                                       cluname))


for (i in 1:dim(Ragewt4)[1]){
  Rage_matrix4[i,]=Ragewt4[i,]/margin.table(Ragewt4, 1)[i]
}
Rage_matrix4
sum(Rage_matrix4[1,])

Ragewt5<-xtabs(wave5$R_WTPERFIN ~ wave5$R_SEX2 + wave5$final_3clust)  # weighted cross table
Ragewt5


Rage_matrix5<-matrix(NA, nrow = dim(Ragewt5)[1], ncol = dim(Ragewt5)[2]
                     ,dimnames = list( c("Male","Female"),
                                       cluname))


for (i in 1:dim(Ragewt5)[1]){
  Rage_matrix5[i,]=Ragewt5[i,]/margin.table(Ragewt5, 1)[i]
}
Rage_matrix5
sum(Rage_matrix5[1,])


####### 1 -MALE ###############
ComGI<-rbind(Rage_matrix1[1,],Rage_matrix2[1,],Rage_matrix3[1,],Rage_matrix4[1,],Rage_matrix5[1,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp1<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ylim(0,1)

###### 2 -FEMALE ################
ComGI<-rbind(Rage_matrix1[2,],Rage_matrix2[2,],Rage_matrix3[2,],Rage_matrix4[2,],Rage_matrix5[2,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp2<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ylim(0,1)


##### COMBINED ###########
Rlegend <- get_legend(Rp1)
Rprow <- plot_grid( Rp1 + theme(legend.position="none"),
                    Rp2 + theme(legend.position="none"),
                    Rlegend,
                    rel_widths = c(1,1,0.3),
                    align = 'vh',
                    labels = c("Male","Female"),
                    label_size = 14,
                    hjust = -1,
                    vjust = 1,
                    nrow = 1,
                    ncol = 3,
                    scale = 0.9
)

Rprow

################################################################################################################
######################################    INCOME VARIABLE    ##########################################
################################################################################################################


#####CREATE AGGREGATE INCOME LEVEL
perclust$INCL25K<-ifelse(((perclust$HHFAMINC>0)&(perclust$HHFAMINC<=5)&(perclust$WAVE!=5)),1,
                         ifelse(((perclust$HHFAMINC>0)&(perclust$HHFAMINC<=3)&(perclust$WAVE==5)),1,0
                         ))
table(perclust$INCL25K)
sum(((perclust$HHFAMINC>0)&(perclust$HHFAMINC<=5)&(perclust$WAVE!=5))|
      ((perclust$HHFAMINC>0)&(perclust$HHFAMINC<=3)&(perclust$WAVE==5)))

perclust$INC25T50K<-ifelse((perclust$HHFAMINC>=6)&(perclust$HHFAMINC<=10)&(perclust$WAVE!=5),1,
                           ifelse((((perclust$HHFAMINC==4)|(perclust$HHFAMINC==5))&(perclust$WAVE==5)),1,
                                  0))
table(perclust$INC25T50K)
sum(((perclust$HHFAMINC>=6)&(perclust$HHFAMINC<=10)&(perclust$WAVE!=5))|
      (((perclust$HHFAMINC==4)|(perclust$HHFAMINC==5))&(perclust$WAVE==5)))

perclust$INC50T75K<-ifelse((perclust$HHFAMINC>=11)&(perclust$HHFAMINC<=15)&(perclust$WAVE!=5),1,
                           ifelse(((perclust$HHFAMINC==6)&(perclust$WAVE==5)),1,0))
table(perclust$INC50T75K)
sum(((perclust$HHFAMINC>=11)&(perclust$HHFAMINC<=15)&(perclust$WAVE!=5))|
      ((perclust$HHFAMINC==6)&(perclust$WAVE==5)))

perclust$INC75KT100K<-ifelse((perclust$HHFAMINC>=16)&(perclust$HHFAMINC<=17)&(perclust$WAVE!=5),1,
                             ifelse(((perclust$HHFAMINC==7)&(perclust$WAVE==5)),1,0))

table(perclust$INC75KT100K)
sum(((perclust$HHFAMINC>=16)&(perclust$HHFAMINC<=17)&(perclust$WAVE!=5))|
      ((perclust$HHFAMINC==7)&(perclust$WAVE==5)))

perclust$INCG100K<-ifelse((perclust$HHFAMINC>=18)&(perclust$WAVE!=5),1,
                          ifelse(((perclust$HHFAMINC>=8)&(perclust$WAVE==5)),1,0))
table(perclust$INCG100K)
sum(((perclust$HHFAMINC>=18)&(perclust$WAVE!=5))|
      ((perclust$HHFAMINC>=8)&(perclust$WAVE==5)))

perclust$CATEINC<-ifelse((perclust$INCL25K==1),1,
                         ifelse((perclust$INC25T50K==1),2,
                                ifelse((perclust$INC50T75K==1),3,
                                       ifelse((perclust$INCG100K==1),5,4
                                       ))))

perclust$CATEINC[perclust$HHFAMINC<0]=NA
table(perclust$CATEINC)
sum(is.na(perclust$CATEINC))
table(perclust$HHFAMINC)
perclust$CATEINC<-factor(perclust$CATEINC)
table(subset(perclust$CATEINC, perclust$WAVE==5))
## END OF CREATE INCOME VARIABLE



wavewt<-xtabs(perclust$R_WTPERFIN ~ perclust$CATEINC + perclust$final_3clust) 

wavewt_matrix<-matrix(NA, nrow =dim(wavewt)[1] , ncol =dim(wavewt)[2] , 
                      dimnames = list(c("     25K-", 
                                        "25K-55K",
                                        "55K-75K",
                                        "75K-100K",
                                        "100K+")
                                      ,cluname))

for (i in 1:dim(wavewt)[1]){
  wavewt_matrix[i,]=wavewt[i,]/margin.table(wavewt, 1)[i]
}
wavewt_matrix
sum(wavewt_matrix[1,])

## OVERALL INCOME DISTRIBUTION BY CLUSTER

barplot(wavewt_matrix,ylim = c(0,1),
        #main="Percent of Age by Cluster (Unweighted)",
        xlab="cluster", col=brewer.pal(n = dim(wavewt_matrix)[1], name = "Set3"),
        ylab="proportion",
        legend.text  = c("25K-", "25K-55K","55K-75K","75K-100K","100K+"),
        args.legend = list(x = "top"), beside=TRUE)
################################################################################################################
######################################    INCOME CHANGE OVER TIME   ##########################################
################################################################################################################
table(perclust$CATEINC)
wave1<-perclust[perclust$WAVE==1,]
wave2<-perclust[perclust$WAVE==2,]
wave3<-perclust[perclust$WAVE==3,]
wave4<-perclust[perclust$WAVE==4,]
wave5<-perclust[perclust$WAVE==5,]


table(wave5$HHFAMINC)
Ragewt1<-xtabs(wave1$R_WTPERFIN ~ wave1$CATEINC + wave1$final_3clust)  # weighted cross table

Ragewt1

Rage_matrix1<-matrix(NA, nrow = dim(Ragewt1)[1], ncol = dim(Ragewt1)[2], 
                     dimnames = list( c("25K-", "25K-55K","55K-75K","75K-100K","100K+"),
                                      cluname))

for (i in 1:dim(Ragewt1)[1]){
  Rage_matrix1[i,]=Ragewt1[i,]/margin.table(Ragewt1, 1)[i]
}
Rage_matrix1
sum(Rage_matrix1[1,])

Ragewt2<-xtabs(wave2$R_WTPERFIN ~ wave2$CATEINC + wave2$final_3clust)  # weighted cross table
Ragewt2

Rage_matrix2<-matrix(NA, nrow = dim(Ragewt2)[1], ncol = dim(Ragewt2)[2], 
                     dimnames = list( c("25K-", "25K-55K","55K-75K","75K-100K","100K+"),
                                      cluname))

for (i in 1:dim(Ragewt2)[1]){
  Rage_matrix2[i,]=Ragewt2[i,]/margin.table(Ragewt2, 1)[i]
}
Rage_matrix2
sum(Rage_matrix2[1,])

Ragewt3<-xtabs(wave3$R_WTPERFIN ~ wave3$CATEINC + wave3$final_3clust)  # weighted cross table
Ragewt3

Rage_matrix3<-matrix(NA, nrow = dim(Ragewt3)[1], ncol = dim(Ragewt3)[2], 
                     dimnames = list( c("25K-", "25K-55K","55K-75K","75K-100K","100K+"),
                                      cluname))

for (i in 1:dim(Ragewt3)[1]){
  Rage_matrix3[i,]=Ragewt3[i,]/margin.table(Ragewt3, 1)[i]
}
Rage_matrix3
sum(Rage_matrix3[1,])

Ragewt4<-xtabs(wave4$R_WTPERFIN ~ wave4$CATEINC+ wave4$final_3clust)  # weighted cross table


Rage_matrix4<-matrix(NA, nrow = dim(Ragewt4)[1], ncol = dim(Ragewt4)[2] 
                     ,dimnames = list( c("25K-", "25K-55K","55K-75K","75K-100K","100K+"),
                                       cluname))


for (i in 1:dim(Ragewt4)[1]){
  Rage_matrix4[i,]=Ragewt4[i,]/margin.table(Ragewt4, 1)[i]
}
Rage_matrix4
sum(Rage_matrix4[1,])



Ragewt5<-xtabs(wave5$R_WTPERFIN ~ wave5$CATEINC + wave5$final_3clust)  # weighted cross table
Ragewt5


Rage_matrix5<-matrix(NA, nrow = dim(Ragewt5)[1], ncol = dim(Ragewt5)[2]
                     ,dimnames = list( c("25K-", "25K-55K","55K-75K","75K-100K","100K+"),
                                       cluname))


for (i in 1:dim(Ragewt5)[1]){
  Rage_matrix5[i,]=Ragewt5[i,]/margin.table(Ragewt5, 1)[i]
}
Rage_matrix5
sum(Rage_matrix5[1,])


####### 1 - LESS THEN 25K###############
ComGI<-rbind(Rage_matrix1[1,],Rage_matrix2[1,],Rage_matrix3[1,],Rage_matrix4[1,],Rage_matrix5[1,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp1<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ylim(0,1)

###### 2- 25 TO 50 K ################
ComGI<-rbind(Rage_matrix1[2,],Rage_matrix2[2,],Rage_matrix3[2,],Rage_matrix4[2,],Rage_matrix5[2,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp2<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ylim(0,1)

###### 3 - 50 TO 75K ################
ComGI<-rbind(Rage_matrix1[3,],Rage_matrix2[3,],Rage_matrix3[3,],Rage_matrix4[3,],Rage_matrix5[3,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp3<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ylim(0,1)


###### 4 -75 TO 100 K ################
ComGI<-rbind(Rage_matrix1[4,],Rage_matrix2[4,],Rage_matrix3[4,],Rage_matrix4[4,],Rage_matrix5[4,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp4<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ylim(0,1)

###### 5 100K AND ABOVE ################
ComGI<-rbind(Rage_matrix1[5,],Rage_matrix2[5,],Rage_matrix3[5,],Rage_matrix4[5,],Rage_matrix5[5,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp5<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ylim(0,1)
Rlegend <- get_legend(Rp1)

##### COMBINED PLOT
Rprow <- plot_grid( Rp1 + theme(legend.position="none"),
                    Rp2 + theme(legend.position="none"),
                    Rp3 + theme(legend.position="none"),
                    Rp4 + theme(legend.position="none"),
                    Rp5 + theme(legend.position="none"),
                    Rlegend,
                    #rel_widths = c(1,1,0.3),
                    align = 'vh',
                    labels = c("    25K-", "25K-55K","55K-75K","75K-100K","100K+"),
                    label_size = 14,
                    hjust = -1,
                    vjust = 1,
                    nrow = 3,
                    ncol = 2,
                    scale = 0.9
)

Rprow


################################################################################################################
#############################                         AGE             ###########################################
################################################################################################################

perclust$A18T25=ifelse((perclust$R_AGE>=18)&(perclust$R_AGE<=25),1,0)
table(perclust$A18T25)
sum((perclust$R_AGE>=18)&(perclust$R_AGE<=25))
perclust$A26T35=ifelse((perclust$R_AGE>25)&(perclust$R_AGE<=35),1,0)
table(perclust$A26T35)
sum((perclust$R_AGE>25)&(perclust$R_AGE<=35))
perclust$A36T45=ifelse((perclust$R_AGE>35)&(perclust$R_AGE<=45),1,0)
table(perclust$A36T45)
sum((perclust$R_AGE>35)&(perclust$R_AGE<=45))
perclust$A46T55=ifelse((perclust$R_AGE>45)&(perclust$R_AGE<=55),1,0)
table(perclust$A46T55)
sum((perclust$R_AGE>45)&(perclust$R_AGE<=55))
perclust$A56T65=ifelse((perclust$R_AGE>55)&(perclust$R_AGE<=65),1,0)
table(perclust$A56T65)
sum((perclust$R_AGE>55)&(perclust$R_AGE<=65))
perclust$AG65=ifelse((perclust$R_AGE>65),1,0)
table(perclust$AG65)
sum(perclust$R_AGE>65)
sum((perclust$A18T25==1)|(perclust$A26T35==1)|(perclust$A36T45==1)|(perclust$A46T55==1)
    |(perclust$A56T65==1)|(perclust$AG65==1))
perclust$CATEAGE<-ifelse(perclust$A18T25==1,1,
                         ifelse(perclust$A26T35==1,2,
                                ifelse(perclust$A36T45==1,3,
                                       ifelse(perclust$A46T55==1,4,
                                              ifelse(perclust$A56T65==1,5,
                                                     ifelse(perclust$AG65==1,6,0))))))
table(perclust$CATEAGE)
perclust$CATEAGE<-factor(perclust$CATEAGE)


wavewt<-xtabs(perclust$R_WTPERFIN ~ perclust$CATEAGE + perclust$final_3clust) 

wavewt_matrix<-matrix(NA, nrow =dim(wavewt)[1] , ncol =dim(wavewt)[2] , 
                      dimnames = list(c("18to25","26to35","36to45","46to55","56to65","65+")
                                      ,cluname))

for (i in 1:6){
  wavewt_matrix[i,]=wavewt[i,]/margin.table(wavewt, 1)[i]
}
wavewt_matrix
sum(wavewt_matrix[1,])

### OVERALL AGE DISTRIBUTION BY CLUSTER
# windows()
barplot(wavewt_matrix,ylim = c(0,1),
        #main="Percent of Age by Cluster (Unweighted)",
        xlab="cluster", col=brewer.pal(n = dim(wavewt_matrix)[1], name = "Set3"),
        ylab = "proportion",
        legend.text  = c("18to25","26to35","36to45","46to55","56to65","65+"),
        args.legend = list(x = "top"), beside=TRUE)


################################################################################################################
######################################     AGE CHANGE OVER TIME    ##########################################
################################################################################################################
wave1<-perclust[perclust$WAVE==1,]
wave2<-perclust[perclust$WAVE==2,]
wave3<-perclust[perclust$WAVE==3,]
wave4<-perclust[perclust$WAVE==4,]
wave5<-perclust[perclust$WAVE==5,]


Ragewt1<-xtabs(wave1$R_WTPERFIN ~ wave1$CATEAGE + wave1$final_3clust)  # weighted cross table

Ragewt1

Rage_matrix1<-matrix(NA, nrow = dim(Ragewt1)[1], ncol = dim(Ragewt1)[2], 
                     dimnames = list(c("18to25","26to35","36to45","46to55","56to65","65+"),
                                     cluname))

for (i in 1:dim(Ragewt1)[1]){
  Rage_matrix1[i,]=Ragewt1[i,]/margin.table(Ragewt1, 1)[i]
}
Rage_matrix1
sum(Rage_matrix1[1,])

Ragewt2<-xtabs(wave2$R_WTPERFIN ~ wave2$CATEAGE + wave2$final_3clust)  # weighted cross table
Ragewt2

Rage_matrix2<-matrix(NA, nrow = dim(Ragewt2)[1], ncol = dim(Ragewt2)[2], 
                     dimnames = list(c("18to25","26to35","36to45","46to55","56to65","65+"),
                                     cluname))

for (i in 1:dim(Ragewt2)[1]){
  Rage_matrix2[i,]=Ragewt2[i,]/margin.table(Ragewt2, 1)[i]
}
Rage_matrix2
sum(Rage_matrix2[1,])

Ragewt3<-xtabs(wave3$R_WTPERFIN ~ wave3$CATEAGE + wave3$final_3clust)  # weighted cross table
Ragewt3

Rage_matrix3<-matrix(NA, nrow = dim(Ragewt3)[1], ncol = dim(Ragewt3)[2], 
                     dimnames = list(c("18to25","26to35","36to45","46to55","56to65","65+"),
                                     cluname))

for (i in 1:dim(Ragewt3)[1]){
  Rage_matrix3[i,]=Ragewt3[i,]/margin.table(Ragewt3, 1)[i]
}
Rage_matrix3
sum(Rage_matrix3[1,])

Ragewt4<-xtabs(wave4$R_WTPERFIN ~ wave4$CATEAGE+ wave4$final_3clust)  # weighted cross table


Rage_matrix4<-matrix(NA, nrow = dim(Ragewt4)[1], ncol = dim(Ragewt4)[2] 
                     ,dimnames = list(c("18to25","26to35","36to45","46to55","56to65","65+"),
                                      cluname))


for (i in 1:dim(Ragewt4)[1]){
  Rage_matrix4[i,]=Ragewt4[i,]/margin.table(Ragewt4, 1)[i]
}
Rage_matrix4
sum(Rage_matrix4[1,])

Ragewt5<-xtabs(wave5$R_WTPERFIN ~ wave5$CATEAGE + wave5$final_3clust)  # weighted cross table
Ragewt5


Rage_matrix5<-matrix(NA, nrow = dim(Ragewt5)[1], ncol = dim(Ragewt5)[2]
                     ,dimnames = list(c("18to25","26to35","36to45","46to55","56to65","65+"),
                                      cluname))


for (i in 1:dim(Ragewt5)[1]){
  Rage_matrix5[i,]=Ragewt5[i,]/margin.table(Ragewt5, 1)[i]
}
Rage_matrix5
sum(Rage_matrix5[1,])

####### 1 - 18 TO 25###############
ComGI<-rbind(Rage_matrix1[1,],Rage_matrix2[1,],Rage_matrix3[1,],Rage_matrix4[1,],Rage_matrix5[1,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp1<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ylim(0,1)

###### 2 - 26 TO 35 ################
ComGI<-rbind(Rage_matrix1[2,],Rage_matrix2[2,],Rage_matrix3[2,],Rage_matrix4[2,],Rage_matrix5[2,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp2<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ylim(0,1)

###### 3 - 36 TO 45 ################
ComGI<-rbind(Rage_matrix1[3,],Rage_matrix2[3,],Rage_matrix3[3,],Rage_matrix4[3,],Rage_matrix5[3,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp3<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ylim(0,1)


###### 4 -46 TO 55 ################
ComGI<-rbind(Rage_matrix1[4,],Rage_matrix2[4,],Rage_matrix3[4,],Rage_matrix4[4,],Rage_matrix5[4,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp4<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ylim(0,1)

###### 5 56 TO 65 ################
ComGI<-rbind(Rage_matrix1[5,],Rage_matrix2[5,],Rage_matrix3[5,],Rage_matrix4[5,],Rage_matrix5[5,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp5<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ylim(0,1)
Rlegend <- get_legend(Rp1)

###### 6 - OVER 65 YEARS OLD################
ComGI<-rbind(Rage_matrix1[6,],Rage_matrix2[6,],Rage_matrix3[6,],Rage_matrix4[6,],Rage_matrix5[6,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp6<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ylim(0,1)
Rlegend <- get_legend(Rp1)

### COMBINED PLOT
Rprow <- plot_grid( Rp1 + theme(legend.position="none"),
                    Rp2 + theme(legend.position="none"),
                    Rp3 + theme(legend.position="none"),
                    Rp4 + theme(legend.position="none"),
                    Rp5 + theme(legend.position="none"),
                    Rp6 + theme(legend.position="none"),
                    #Rlegend,
                    #rel_widths = c(1,1,0.3),
                    align = 'vh',
                    labels = c("18to25","26to35","36to45","46to55","56to65","65+"),
                    label_size = 14,
                    hjust = -1,
                    vjust = 1,
                    nrow = 3,
                    ncol = 2,
                    scale = 0.9
)

Rprow


################################################################################################################
######################################    HHSIZE       #####################################
################################################################################################################
table(perclust$HHSIZE)
perclust$HHSIZE1=ifelse(perclust$HHSIZE==0,0,
                        ifelse(perclust$HHSIZE==1,1,
                               ifelse(perclust$HHSIZE==2,2,
                                      ifelse(perclust$HHSIZE==3,3,
                                             ifelse(perclust$HHSIZE>3,4,-99)))))
table(perclust$HHSIZE1)


wt1<-xtabs(perclust$R_WTPERFIN ~ perclust$HHSIZE1 + perclust$final_3clust)
wt_matrix<-matrix(NA, nrow = dim(wt1)[1], ncol = dim(wt1)[2], 
                  dimnames = list( c("1","2","3","4+"),
                                   cluname))

for (i in 1:dim(wt1)[1]){
  wt_matrix[i,]<-wt1[i,]/margin.table(wt1,1)[i]
}
wt_matrix
sum(wt_matrix[1,])


barplot(wt_matrix,
        #main="Percent of Age by Cluster (Unweighted)",
        xlab="cluster", col=brewer.pal(n = dim(wt_matrix)[1], name = "Set3"),
        ylab = "proportion",
        legend.text  = c("1","2","3","4+"),
        args.legend = list(x = "top"), beside=TRUE)


###################################################################################################
####################################  HHSIZE CHANGE OVER TIME   ###################################
###################################################################################################
wave1<-perclust[perclust$WAVE==1,]
wave2<-perclust[perclust$WAVE==2,]
wave3<-perclust[perclust$WAVE==3,]
wave4<-perclust[perclust$WAVE==4,]
wave5<-perclust[perclust$WAVE==5,]


Ragewt1<-xtabs(wave1$R_WTPERFIN ~ wave1$HHSIZE1 + wave1$final_3clust)  # weighted cross table

Ragewt1

Rage_matrix1<-matrix(NA, nrow = dim(Ragewt1)[1], ncol = dim(Ragewt1)[2], 
                     dimnames = list( c("1","2","3","4+"),
                                      cluname))

for (i in 1:dim(Ragewt1)[1]){
  Rage_matrix1[i,]=Ragewt1[i,]/margin.table(Ragewt1, 1)[i]
}
Rage_matrix1
sum(Rage_matrix1[1,])

Ragewt2<-xtabs(wave2$R_WTPERFIN ~ wave2$HHSIZE1 + wave2$final_3clust)  # weighted cross table
Ragewt2

Rage_matrix2<-matrix(NA, nrow = dim(Ragewt2)[1], ncol = dim(Ragewt2)[2], 
                     dimnames = list( c("1","2","3","4+"),
                                      cluname))

for (i in 1:dim(Ragewt2)[1]){
  Rage_matrix2[i,]=Ragewt2[i,]/margin.table(Ragewt2, 1)[i]
}
Rage_matrix2
sum(Rage_matrix2[1,])

Ragewt3<-xtabs(wave3$R_WTPERFIN ~ wave3$HHSIZE1 + wave3$final_3clust)  # weighted cross table
Ragewt3

Rage_matrix3<-matrix(NA, nrow = dim(Ragewt3)[1], ncol = dim(Ragewt3)[2], 
                     dimnames = list( c("1","2","3","4+"),
                                      cluname))

for (i in 1:dim(Ragewt3)[1]){
  Rage_matrix3[i,]=Ragewt3[i,]/margin.table(Ragewt3, 1)[i]
}
Rage_matrix3
sum(Rage_matrix3[1,])

Ragewt4<-xtabs(wave4$R_WTPERFIN ~ wave4$HHSIZE1+ wave4$final_3clust)  # weighted cross table


Rage_matrix4<-matrix(NA, nrow = dim(Ragewt4)[1], ncol = dim(Ragewt4)[2] 
                     ,dimnames = list( c("1","2","3","4+"),
                                       cluname))


for (i in 1:dim(Ragewt4)[1]){
  Rage_matrix4[i,]=Ragewt4[i,]/margin.table(Ragewt4, 1)[i]
}
Rage_matrix4
sum(Rage_matrix4[1,])



Ragewt5<-xtabs(wave5$R_WTPERFIN ~ wave5$HHSIZE1 + wave5$final_3clust)  # weighted cross table
Ragewt5


Rage_matrix5<-matrix(NA, nrow = dim(Ragewt5)[1], ncol = dim(Ragewt5)[2]
                     ,dimnames = list( c("1","2","3","4+"),
                                       cluname))


for (i in 1:dim(Ragewt5)[1]){
  Rage_matrix5[i,]=Ragewt5[i,]/margin.table(Ragewt5, 1)[i]
}
Rage_matrix5
sum(Rage_matrix5[1,])

 
####### 1 - ONE PERSON HOUSEHOLD ###############
ComGI<-rbind(Rage_matrix1[1,],Rage_matrix2[1,],Rage_matrix3[1,],Rage_matrix4[1,],Rage_matrix5[1,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp1<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ylim(0,1)

###### 2 - 2 PERSON HOUSEHOLD ################
ComGI<-rbind(Rage_matrix1[2,],Rage_matrix2[2,],Rage_matrix3[2,],Rage_matrix4[2,],Rage_matrix5[2,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp2<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ylim(0,1)

###### 3 - 3 PERSON HOUSEHOLD ################
ComGI<-rbind(Rage_matrix1[3,],Rage_matrix2[3,],Rage_matrix3[3,],Rage_matrix4[3,],Rage_matrix5[3,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp3<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ylim(0,1)


###### 4 - FOUR AND MORE PERSON HOUSEHOLD ################
ComGI<-rbind(Rage_matrix1[4,],Rage_matrix2[4,],Rage_matrix3[4,],Rage_matrix4[4,],Rage_matrix5[4,])
dfComGI=stack(as.data.frame(ComGI))
dfComGI$year<-rep(c("1990","1995","2001","2009","2017"),
                  dim(ComGI)[2])
names(dfComGI)<-c("proportion","cluster","year")
Rp4<-ggplot(data=dfComGI, aes(x=year, y=proportion, group=cluster)) +geom_line(aes(color=cluster))+
  geom_point(aes(color=cluster))+guides(fill=guide_legend(title="Cluster"))+ylim(0,1)


##### COMBINED PLOT
Rprow <- plot_grid( Rp1 + theme(legend.position="none"),
                    Rp2 + theme(legend.position="none"),
                    Rp3 + theme(legend.position="none"),
                    Rp4 + theme(legend.position="none"),
                    Rp5 + theme(legend.position="none"),
                    Rp6 + theme(legend.position="none"),
                    #Rlegend,
                    #rel_widths = c(1,1,0.3),
                    align = 'vh',
                    labels = c("1","2","3","4+"),
                    label_size = 14,
                    hjust = -1,
                    vjust = 1,
                    nrow = 2,
                    ncol = 2,
                    scale = 0.9
)

Rprow
