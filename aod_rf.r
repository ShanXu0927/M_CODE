library(randomForest) 

traindata<-read.table("EM05_20aod_select.csv",header=TRUE,sep=",")  

set.seed(1234)    
rf_ntree1<-randomForest(aod ~ bc_a+bc_b+ch4_a+ch4_b+co_a+co_b+nh3_a+nh3_b+nmvoc_a+nmvoc_b+no_a+no_b+oc_a+oc_b+so2_a+so2_b+lon+lat,data=traindata, ntree=200, important=TRUE, proximity=TRUE)
plot(rf_ntree1)

mytraindata<-cbind(traindata$bc_a,traindata$bc_b,traindata$ch4_a,traindata$ch4_b,traindata$co_a,traindata$co_b,traindata$nh3_a,traindata$nh3_b,traindata$nmvoc_a,traindata$nmvoc_b,traindata$no_a,traindata$no_b,traindata$oc_a,traindata$oc_b,traindata$so2_a,traindata$so2_b,traindata$lon,traindata$lat)      #trainx,cbind根据列进行合并，即叠加所有列
result=rfcv(mytraindata,traindata$aod,cv.flod=10)

with(result,plot(n.var,error.cv,log="x",type="o",lwd=2))

set.seed(4321)
aod_rf<-randomForest(aod ~ bc_a+bc_b+ch4_a+ch4_b+co_a+co_b+nh3_a+nh3_b+nmvoc_a+nmvoc_b+no_a+no_b+oc_a+oc_b+so2_a+so2_b+lon+lat,data=traindata, ntree=100, important=TRUE, proximity=TRUE)

aod_rf
importance(aod_rf)

set.seed(4321)
aod_rf1<-randomForest(aod ~ bc_a+ch4_a+co_a+nh3_a+nmvoc_a+no_a+oc_a+so2_a+lon+lat,data=traindata, ntree=100, important=TRUE, proximity=TRUE)

aod_rf1
importance(aod_rf1)

aod_pre<-predict(aod_rf1,traindata)
aod_rf_result<-cbind(traindata$aod,aod_pre)
write.table(aod_rf_result,"aod05_20_rf_result.csv",sep=",")

aoddata<-read.table("EM05_20aod_select.csv",header=TRUE,sep=",")

# cross validation
library("caret")  
set.seed(4321) 
folds<-createFolds(y=aoddata$ID,k=10) 
aod_cv_mai<-{}
aod_cv_pre<-{}
for(i in 1:10){  
  traindata<-aoddata[-folds[[i]],]  
  testdata<-aoddata[folds[[i]],]
  aodcv_rf <- randomForest(aod ~ bc_a+ch4_a+co_a+nh3_a+nmvoc_a+no_a+oc_a+so2_a+lon+lat, data=traindata, ntree=100,importance=TRUE, proximity=TRUE)
  #rf <- randomForest(Species ~ ., data=training, ntree=100, proximity=TRUE) 
  test_pre<-predict(aodcv_rf,testdata)  #pre
  aod_cv_mai<-c(aod_cv_mai,testdata$aod)
  aod_cv_pre<-c(aod_cv_pre,test_pre)

}  
aod_cv_result<-cbind(aod_cv_mai,aod_cv_pre)
write.table(aod_cv_result,"aod05_20_cv_result.csv",sep=",")