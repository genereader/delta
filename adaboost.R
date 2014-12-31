library(ada)
tdata <- read.table("trainData.txt")
pdata <- read.table("predictData.txt")
nc <- dim(tdata)[2]
colnames(pdata) <- colnames(tdata[,2:nc])
adamodel <- ada(x=tdata[,2:nc],y=tdata[,1],iter=100)
adapred.fwd <- predict(adamodel, newdata=pdata,type="probs")
pdata[,seq(4,nc,4)-1] <- -pdata[,seq(4,nc,4)-1]
adapred.bkwd <- predict(adamodel, newdata=pdata,type="probs")
prob.mat <- cbind(adapred.fwd[,2],adapred.bkwd[,2])
write.table(apply(prob.mat,1,min),"pred",quote=F,row.names=F,col.names=F)
