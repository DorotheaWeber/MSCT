##############################################
######## caret function boosted nnet #########
##############################################


load("caret_nnet.R")

boostNNet = caret.nnet

boostNNet$label = "Boosted Neural Network"

boostNNet$type = "Classification"    

boostNNet$parameters = data.frame(parameter=c("size", "decay", "B"), class=rep("numeric",3), label=c("Hidden Units", "Weight Decay", "Boosting Steps"))

boostNNet$fit = function (x, y, wts, param, lev, last, classProbs, ...) 
{
  library(nnet)
  dat <- if (is.data.frame(x)) 
    x
  else as.data.frame(x)
  dat$.outcome <- y
    
    dist = rep(1/length(y), length(y))
  
    f_boost = matrix(NA, nrow=param$B, ncol=length(y))
    model = NULL
    
    for(i in 1:param$B){
      
      model[[i]] = nnet:::nnet(.outcome ~ ., data = dat, weights=dist, size=param$size, decay=param$decay)
      prob = predict(model[[i]], type="raw") # probability for "no mutation"
      
      f_boost[i,] = 0.5 * log((1-prob)/(prob))
      
      dist = dist * exp(-ifelse(dat$.outcome=="j", 1, -1) * f_boost[i,])
      dist = dist/sum(dist)
      
    }
    
    F_all = apply(f_boost, 2, sum)
    
    class_prob.tmp = ifelse(is.na(exp(F_all)/(1 + exp(F_all)))==T,
                            1-(1/(1 + exp(F_all))), exp(F_all)/(1 + exp(F_all)))
    
    probs = cbind(class_prob.tmp, 1-class_prob.tmp)
    colnames(probs) = c("j", "n")
    out=list(models=model, prob = probs)
  out
}

boostNNet$predict = function (modelFit, newdata, submodels = NULL) 
{
  require(nnet)
  B = length(modelFit$models)
  
  f_boost = matrix(NA, nrow=B, ncol=dim(newdata)[1])
  
  for(i in 1:B){
    class(modelFit$models[[i]]) <- "nnet"
    prob = nnet:::predict.nnet(modelFit$models[[i]], newdata, type="raw")
    f_boost[i,] = 0.5 * log((1-prob)/prob)
  }
  
  F_all = apply(f_boost, 2, sum)
  
  out = ifelse(F_all>=0, "j", "n")
  out
}

boostNNet$prob = function (modelFit, newdata, submodels = NULL) 
{
  require(nnet)
  
  B = length(modelFit$models)
  
  f_boost = matrix(NA, nrow=B, ncol=dim(newdata)[1])
  
  for(i in 1:B){
    class(modelFit$models[[i]]) <- "nnet"
    prob = nnet:::predict.nnet(modelFit$models[[i]], newdata, type="raw")
    f_boost[i,] = 0.5 * log((1-prob)/prob)
  }
  
  F_all = apply(f_boost, 2, sum)
  
  class_prob.tmp = ifelse(is.na(exp(F_all)/(1 + exp(F_all)))==T,
                          1-(1/(1 + exp(F_all))), exp(F_all)/(1 + exp(F_all)))
  
  class_prob = cbind(class_prob.tmp, 1-class_prob.tmp)
  
  colnames(class_prob) = c("n", "j")
  out=class_prob
  out
}


boostNNet$tags = c("Neural Network", "Ensemble Model", "Boosting", "L2 Regularization", "Accepts Case Weights")