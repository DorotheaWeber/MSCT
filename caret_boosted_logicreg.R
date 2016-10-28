################################################
######## caret function boosted logreg #########
################################################


load("caret_logreg.R")

boostLogreg = caret.logreg

boostLogreg$label = "Boosted Logic Regression"

boostLogreg$type = "Classification"    

boostLogreg$parameters = data.frame(parameter=c("treesize", "ntrees", "B"), class=rep("numeric",3), label=c("Number of Leaves", "Number of Trees", "Boosting Steps"))

boostLogreg$grid = function (x, y, len = NULL, search = "grid") {
  expand.grid(ntrees = (1:3) + 1, treesize = 2^(1 + (1:len)), B=15)
}
  

boostLogreg$fit = function (x, y, wts, param, lev, last, classProbs, ...) 
{
  isReg <- is.numeric(y)
  if (is.factor(y)) {
     y <- ifelse(y == levels(y)[1], 1, 0)
  }
  
  dist = rep(1/length(y), length(y))
  
  model = NULL
  f_boost = matrix(NA, nrow=param$B, ncol=length(y))
  
  for(i in 1:param$B){
    
    model[[i]] = LogicReg:::logreg(resp=y, bin=x, wgt=dist, type=3, select=1, 
                                   ntrees = param$ntrees, tree.control = logreg.tree.control(treesize = param$treesize), ...)
    prob = LogicReg:::predict.logreg(model[[i]])
    
    f_boost[i,] = 0.5 * log((prob)/(1-prob))
    
    dist = dist * exp(-ifelse(y==1, 1, -1) * f_boost[i,])
    dist = dist/sum(dist)
    
  }
  
  F_all = apply(f_boost, 2, sum)
  
  class_prob.tmp = ifelse(is.na(exp(F_all)/(1 + exp(F_all)))==T,
                          1-(1/(1 + exp(F_all))), exp(F_all)/(1 + exp(F_all)))
  
  probs = cbind(class_prob.tmp, 1-class_prob.tmp)
  colnames(probs) = c("j", "n")
  out = list(models=model, prob=probs)
  
  out
  
}

boostLogreg$predict = function (modelFit, newdata, submodels = NULL) 
{
  require(LogicReg)
  B = length(modelFit$models)
  
  f_boost = matrix(NA, nrow=B, ncol=dim(newdata)[1])
  
  for(i in 1:B){
    class(modelFit$models[[i]]) <- "logreg"
    prob = LogicReg:::predict.logreg(modelFit$models[[i]], newbin=newdata)
    f_boost[i,] = 0.5 * log(prob/(1-prob))
  }
  
  F_all = apply(f_boost, 2, sum)
  
  out = ifelse(F_all>=0, "j", "n")
  
  out
}

boostLogreg$prob = function (modelFit, newdata, submodels = NULL) 
{
  require(LogicReg)
  
  B = length(modelFit$models)
  
  f_boost = matrix(NA, nrow=B, ncol=dim(newdata)[1])
  
  for(i in 1:B){
    class(modelFit$models[[i]]) <- "logreg"
    prob = LogicReg:::predict.logreg(modelFit$models[[i]], newbin=newdata)
    f_boost[i,] = 0.5 * log(prob/(1-prob))
  }
  
  F_all = apply(f_boost, 2, sum)
  
  class_prob.tmp = ifelse(is.na(exp(F_all)/(1 + exp(F_all)))==T,
                         1-(1/(1 + exp(F_all))), exp(F_all)/(1 + exp(F_all)))
  
  class_prob = cbind(class_prob.tmp, 1-class_prob.tmp)
  
  colnames(class_prob) = c("j", "n")
  
  out = class_prob
  out
}