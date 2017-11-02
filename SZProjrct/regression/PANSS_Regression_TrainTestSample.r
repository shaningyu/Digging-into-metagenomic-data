workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

#data preprocessing
motu <- read.table("../../SZData/mOTU.0.05.profile",header = 1,row.names = 1)
state <- read.table("../../SZData/state171.txt",header = 1,row.names = 1)
phe <- read.table("../../SZData/panss.txt",header = 1,row.names = 1)
phe <- na.omit(phe)

colnames(motu) <- gsub("\\.","-",colnames(motu))
motu <- motu[,pmatch(rownames(phe),colnames(motu))]
state <- state[pmatch(rownames(phe),rownames(state)),]
case_motu <- motu[,which(state$state2==1)]

set.seed(1)
train_ca <- sample(1:90,67)
train_case <- case_motu[,train_ca]
test_ca <- setdiff(1:90,train_ca)
test_case <- case_motu[,test_ca]

train.x <- train_case
test.x <- test_case
train.y <- phe[pmatch(colnames(train.x),rownames(phe)),]
test.y <- phe[pmatch(colnames(test.x),rownames(phe)),]

train.x <- t(train.x)
test.x <- t(test.x)
library(glmnet)

grid=10^ seq (10,-2, length =1000)
pseR2 <- matrix(0,ncol(phe),1)
lasso.pred <- matrix(0,23,ncol(phe))

##by changing alpha, the pseR2 is better when it's 0.15. 
for(i in 1:ncol(phe)){
  set.seed(1)
  cv.out <- cv.glmnet(train.x,train.y[,i],family = "gaussian", alpha = 0.15,lambda = grid,nfolds = 5)
  bestlam <- cv.out$lambda.min
  lasso.pred[,i] <- predict(cv.out$glmnet.fit,s=bestlam,newx = test.x,type="response")
  pseR2[i,1] <- 1-mean((lasso.pred[,i]-test.y[,i])^2)/mean((test.y[,i]-mean(test.y[,i]))^2)
  rownames(pseR2) <- colnames(phe)
  
##the value 0.15 is a little bit small, but it's fine for check the prediction.
    if(pseR2[i,1]>0.15){ 
    lasso.coef <- predict(cv.out$glmnet.fit,type='coefficients',s=bestlam)[1:ncol(test.x),]
    lasso.f <- lasso.coef[lasso.coef!=0]
    name <- paste(colnames(test.y)[i],"lasso_regress",sep = "_")
    write.table(lasso.f,name,quote=F,row.names = T,sep = "\t")
  }
}



