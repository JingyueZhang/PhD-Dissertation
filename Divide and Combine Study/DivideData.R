### Read in data from csv file
mydat=read.csv(file="C:/Users/jiz13007/Documents/Pattern Recognition/Book Chapter/Data/elderly data.csv",header=T,sep=",")
dim(mydat)
names(mydat)
### Randomly resample data
set.seed(1234)
mydat <- mydat[sample(nrow(mydat)),]
for (i in 1:50){
  subdat<-mydat[((i-1)*1273+1):(i*1273),]
  print(dim(subdat))
  filename<-paste("C:/Users/jiz13007/Documents/Pattern Recognition/Book Chapter/Data/Split Data/subsample_",i,".csv")
  write.csv(subdat, file = filename)}

filename<-paste("C:/Users/jiz13007/Documents/Pattern Recognition/Book Chapter/Data/Split Data/subsample_",50,".csv")
subdata50=read.csv(file=filename,header=T,sep=",")
dim(subdata50)
names(subdata50)

subdata50_new<-subdata50[!is.na(subdata50$UNIPERID),]
dim(subdata50_new)

write.csv(subdata50_new, file = filename)

for (i in 1:50){
  print (paste(((i-1)*1273+1), ":",(i*1273), sep = ""))
}
