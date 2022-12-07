
library(randomForest)    

hurdata<-read.table("HUR06_17.csv",header=TRUE,sep=",")    


set.seed(5678)
hur_rf<-randomForest(average ~ X+Y+DEM+CSIRO,data=hurdata, ntree=100, important=TRUE, proximity=TRUE)
hur_rf
varImpPlot(hur_rf)

hur_rf
importance(hur_rf)

hur_pre<-predict(hur_rf,hurdata)
hur_rf_result<-cbind(hurdata$average,hur_rf$predicted)
write.table(hur_rf_result,"hur_rf_result1.csv",sep=",")


prdata<-read.table("PR06_17.csv",header=TRUE,sep=",")   

set.seed(5678)
pr_rf<-randomForest(average~X+Y+DEM+bcc+CSIRO,data=prdata, ntree=100, important=TRUE, proximity=TRUE)
varImpPlot(pr_rf)

pr_rf
importance(pr_rf)

pr_pre<-predict(pr_rf,prdata)
pr_rf_result<-cbind(prdata$average,pr_rf$predicted)
write.table(pr_rf_result,"pr_rf_result1.csv",sep=",")


psdata<-read.table("PS06_17.csv",header=TRUE,sep=",")  

set.seed(5678)
ps_rf<-randomForest(average~X+DEM+bcc+CSIRO,data=psdata, ntree=100, important=TRUE, proximity=TRUE)

ps_rf
importance(ps_rf)


ps_pre<-predict(ps_rf,psdata)
ps_rf_result<-cbind(psdata$average,ps_rf$predicted)
write.table(ps_rf_result,"ps_rf_result1.csv",sep=",")


tsdata<-read.table("TS06_17.csv",header=TRUE,sep=",")  

set.seed(8765)
ts_rf1<-randomForest(average~Y+bcc+CSIRO+GFDL,data=tsdata, ntree=100, important=TRUE, proximity=TRUE)
ts_rf1
importance(ts_rf1)

ts_pre<-predict(ts_rf1,tsdata)
ts_rf_result<-cbind(tsdata$average,ts_rf1$predicted)
write.table(ts_rf_result,"ts_rf_result1.csv",sep=",")