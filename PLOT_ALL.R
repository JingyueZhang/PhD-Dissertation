rm(list=ls()) # clear memory
library(ggplot2)
library(RColorBrewer)
require(gridExtra)
library(cowplot)

stack_density1=readRDS(
  file="C:/Users/jiz13007/Documents/Pattern Recognition/2018 TRB paper/Data and Script/Aggregate/Final Result/stacked_density_3clusters_clu1.rds")
stack_density2=readRDS(
  file="C:/Users/jiz13007/Documents/Pattern Recognition/2018 TRB paper/Data and Script/Aggregate/Final Result/stacked_density_3clusters_clu1.rds")
stack_density3=readRDS(
  file="C:/Users/jiz13007/Documents/Pattern Recognition/2018 TRB paper/Data and Script/Aggregate/Final Result/stacked_density_3clusters_clu1.rds")

perclust=readRDS(file="C:/Users/jiz13007/Documents/Pattern Recognition/2018 TRB paper/Data and Script/Aggregate/Final Result/perclust_all.rds")
n1<-sum(perclust$ward_clu3==1)
n2<-sum(perclust$ward_clu3==2)
n3<-sum(perclust$ward_clu3==3)
# stack_density1$TYPE1<-ifelse(stack_density1$TYPE=="Travel",5,
#                              ifelse(stack_density1$TYPE=="Home",1,
#                                     ifelse(stack_density1$TYPE=="Mandatory",2,
#                                            ifelse(stack_density1$TYPE=="Maintenance",3,
#                                                   ifelse(stack_density1$TYPE=="Discretionary",4,0)))))
# 
# table(stack_density1$TYPE)
# table(stack_density1$TYPE1)

#Rstack_density1<-stack_density1[with(stack_density1, order(stack_density1$time, stack_density1$TYPE1)), ]

# Recode time from 4 am to 4 am of the next day
dim(stack_density1)
time_list<-c()
i=1
hr=4
minu=1
while (i<=1440){
  tim=paste(hr,":",if(minu<10){"0"},minu, sep="")
  time_list<-c(time_list, rep(tim,5))
  minu=minu+1
  if(minu>59){
    hr=hr+1
    minu=0
  }
  if(hr>23){
    hr=0
  }
  
  i=i+1
}

stack_density1$time_new=time_list
stack_density2$time_new=time_list
stack_density3$time_new=time_list

windows()
plot1

plot1<-ggplot(stack_density1, aes(x=time, y=(values*100), fill=TYPE)) + 
  geom_area(alpha=0.6 , size=1,colour="black") +
  ylab("Percentage")+ xlab("Time")
plot2<-ggplot(stack_density2, aes(x=time, y=(values*100), fill=TYPE)) + 
  geom_area(alpha=0.6 , size=1,colour="black")+
  ylab("Percentage")+ xlab("Time")
plot3<-ggplot(stack_density3, aes(x=time, y=(values*100), fill=TYPE)) + 
  geom_area(alpha=0.6 , size=1,colour="black")+
  ylab("Percentage")+ xlab("Time")

Rlegend <- get_legend(plot1 + theme(legend.position="bottom")+ theme(legend.text=element_text(size=12)))

Rprow<-plot_grid( plot1 + theme(legend.position="none"),
                  plot2 + theme(legend.position="none"),
                  plot3 + theme(legend.position="none"),
                  align = 'vh',
                  labels = c(paste0("Cluster 1 ","(N=",n1,")"),
                             paste0("Cluster 2 ","(N=",n2,")"),
                             paste0("Cluster 3 ","(N=",n3,")")),
                  label_size = 14,
                  hjust = -1,
                  vjust = 1,
                  nrow = 1,
                  ncol = 3,
                  scale = 1
)
windows()
Rprow


windows()
Rprow1 <- plot_grid( Rprow,
                     Rlegend,
                     align = 'vh',
                     rel_heights = c(1, .2),
                     nrow = 2,
                     ncol = 1,
                     scale = 1
)
Rprow1