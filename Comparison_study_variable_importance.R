########################################
######## Script for Able ###############
########################################

# 5.1 Comparison Study - Variable importance
# without CV

# drug is given as argument to R script

args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

######################################################

.libPaths("~/Rlibs")

require(nnet)
require(randomForest)
require(gbm)
require(caret)
require(pROC)
require(rpart)
require(caretEnsemble)
require(adabag)
require(glmnet)
require(RSNNS)

######################################################

# load customized functions

source( "caret_boosted_nnet.R" )
source( "caret_boosted_logicreg.R" )


# methods parallel with caretEnsemble

load( "Subsample.RData" ) # index list for my_control

set.seed(17)

my_control <- trainControl(
  method="boot632",
  number=200,
  repeats=1,
  savePredictions="final",
  classProbs=TRUE,
  index=index, # rows which are used for training iterations
  summaryFunction=twoClassSummary,
  seed=NULL
)



# get data
load( paste("IC50_data_", drugname, ".RData", sep="") ) # name "data"
  
  # only for drugs with at least 4 differences in response, only drugs with 40 celllines
  if(length(which(data[,1]==1))>3 && length(which(data[,1]==0))>3 && length(data[,1])==40){
    
    data[,1]= as.factor(ifelse(data[,1]==1, "j", "n")) # nnet needs factor with two classes (not 0 and 1)
    
    # formula
    n <- colnames(data)
      
      # spilt in training and test set
      
      train = data
      levels(train[,1])=c("j", "n")
      train[,1] = relevel(train[,1], "j")
      
      # set name for output
      name = paste("result_", drugname, sep="")
      
      # do parallel caret model calculation with ROC
      set.seed(17)
      assign(name, caretList(x=train[,-1], y=train[,1],
                             trControl=my_control,
                             metric="ROC",
                             methodList=c("rf"),# bagged tree/random forest
                             tuneList=list(
                               # boosted tree
                               caretModelSpec(method="gbm", tuneGrid=data.frame(.n.trees=c(1000,2000,3000), 
                                                                                .shrinkage=rep(c(0.1,0.01,0.001),3), 
                                                                                .n.minobsinnode=5, 
                                                                                .interaction.depth=6)),
                               # tree
                               caretModelSpec(method="rpart"),
                               # neural net
                               caretModelSpec(method="nnet", tuneGrid=data.frame(.decay=sort(rep(c( 5e-3,  5e-5,  5e-7),5)), 
                                                                                 .size=seq(2, 10, by=2)), 
                                              trace=F),
                               # bagged neural net
                               caretModelSpec(method="avNNet", tuneGrid=data.frame(.decay=sort(rep(c(5e-3,  5e-5,  5e-7),5)), 
                                                                                   .size=seq(2, 10, by=2), .bag=T), 
                                              trace=F),
                               # neural net 3 hidden layers
                               caretModelSpec(method="mlpWeightDecayML", tuneGrid=data.frame(.layer1=seq(10,100, by=10), 
                                                                                             .layer2=seq(5,50, by=5), 
                                                                                             .layer3=seq(2,20, by=2), 
                                                                                             .decay=sort(rep(c(5e-3, 5e-5, 5e-7),10)))),
                               # logic regression
                               caretModelSpec(method="logreg"),
                               # bagged logic regression
                               caretModelSpec(method="logicBag", tuneGrid=data.frame(.ntrees=1,
                                                                                     .nleaves=c(5, 10, 15),
                                                                                     .glm.if.1tree=T)),
                               # linear model --> elastic net
                               caretModelSpec(method="glmnet"),
                               # boosted tree --> Breiman method
                               caretModelSpec(method='AdaBoost.M1', 
                                              tuneGrid=data.frame(.coeflearn="Breiman", 
                                                                  .maxdepth=seq(1,9, by=2),
                                                                  .mfinal=50))
                             ),
                             continue_on_fail=T
      )
      )
      
      # boosted logic regression
      
      logreg.grid <- expand.grid(.treesize =10,
                                 .ntrees = 1,
                                 .B=c(8,10,12,14)
      )
      
      set.seed(17)
      my_boostLogreg = test_caret_logreg = train(y=train[,1], x=train[,-1], 
                                                 method=boostLogreg, 
                                                 tuneGrid = logreg.grid, 
                                                 metric="ROC",
                                                 trControl = my_control)
      
      # boosted neural network
      
      nnet.grid <- expand.grid(.decay =c(5e-3, 5e-5, 5e-7),
                               .size = 10,
                               .B=c(10,20,30,50,70)
      )
      
      set.seed(17)
      my_boostNNet = train(y=train[,1], x=train[,-1], 
                           method=boostNNet, 
                           tuneGrid = nnet.grid, 
                           metric="ROC",
                           trace=F, 
                           trControl = my_control)
      
      # add to result list
      
      result = get(name)
      result[["boostLogreg"]] = my_boostLogreg
      result[["boostNNet"]] = my_boostNNet
      
      # save results
      
      filepath="Results_200BS_all"
      dir.create(filepath)
      save(result, file=paste(filepath,"/CCLE_melanoma_result_", drugname, ".RData", sep=""))
      
    }