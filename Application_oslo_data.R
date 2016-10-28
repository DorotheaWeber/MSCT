########################################
######## Script for Able ###############
########################################

# 5.2 Application to Oslo data

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
require(glmnet)
require(adabag)
require(RSNNS)

######################################################

# load customized functions

source( "caret_boosted_nnet.R" )
source( "caret_boosted_logicreg.R" )


# methods parallel with caretEnsemble

load( "Subsample2.RData" )

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
load( paste("train_", drugname, ".RData", sep="") ) # name "data"

  # formula
  n <- colnames(data)
    
    excl = c(which(colnames(data) == "TTN_MUT"), which(colnames(data) == "OBSCN_MUT"))
    train = data[,-excl]
    levels(train[,1])=c("j", "n")
    train[,1] = relevel(train[,1], "j")
    
    # set name for output
    name = paste("result_kinome_", drugname, sep="")
    
    # do parallel caret model calculation with ROC
    set.seed(17)
    assign(name, caretList(x=train[,-1], y=train[,1],
                           trControl=my_control,
                           metric="ROC",
                           methodList=c("rf"),# bagged tree
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
                             # neural net several hidden layers
                             caretModelSpec(method="mlpWeightDecayML", tuneGrid=data.frame(.layer1=seq(10,100, by=10), 
                                                                                           .layer2=seq(5,50, by=5), 
                                                                                           .layer3=seq(2,20, by=2), 
                                                                                           .decay=sort(rep(c(5e-3, 5e-5, 5e-7),10)))),
                             # logic regression
                             caretModelSpec(method="logreg"),
                             # bagged logic regression
                             caretModelSpec(method="logicBag", tuneGrid=data.frame(.ntrees=1,
                                                                                   .nleaves=c(5, 10, 15) )),
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
    
    result = get(name)
    result[["boostLogreg"]] = my_boostLogreg
    result[["boostNNet"]] = my_boostNNet
    
    filepath=paste("Results_kinome_",drugname, sep="")
    dir.create(filepath)
    save(result, file=paste(filepath, "/CCLE_melanoma_result_kinome_excl_", drugname, ".RData", sep=""))