############################################################
########## neuralnet with .362 bootstrap ###################
############################################################

boot.neuralnet = function(f, data, excl, constwts, B, hiddenK, priorK){ 
  
  fit.all.P = neuralnet(f, data, hidden=hiddenK, exclude=excl, constant.weights=constwts, 
                        linear.output=F, rep=1, err.fct = "sse")
  
  #err.P = fit.all.P$result.matrix["error",]
  err.P = mean((fit.all.P$net.result[[1]]-data[,1])^2)
  
  e632.P = 0
  
  for(i in 1:B){
    # bootstrap dataset with replacment
    
    bootS = sample(1:dim(data)[1], dim(data)[1], replace = T)
    
    data.boot = data[bootS,]
    
    dataT.boot = data[-unique(bootS),]
    
    fit.boot = neuralnet(f, data.boot, hidden=hiddenK, exclude=excl, constant.weights=constwts, 
                         linear.output=F, rep=1, err.fct="sse")
    
    pred = compute(fit.boot, dataT.boot[,-1])$net.result
    
    sse = mean((pred-dataT.boot[,1])^2)
    
    e632.P = e632.P + 0.632*sse
  }
  
  e632.P = e632.P/B
  
  e632.final.P = e632.P + 0.368 * err.P
  
  
  # without prior info for comparison
  fit.all = neuralnet(f, data, hidden=(hiddenK-priorK),
                      linear.output=F, rep=1, err.fct = "sse")
  
  #err = fit.all$result.matrix["error",]
  err = mean((fit.all$net.result[[1]]-data[,1])^2)
  
  e632 = 0
  
  for(i in 1:B){
    # bootstrap dataset with replacment
    
    bootS = sample(1:dim(data)[1], dim(data)[1], replace = T)
    
    data.boot = data[bootS,]
    
    dataT.boot = data[-unique(bootS),]
    
    fit.boot = neuralnet(f, data.boot, hidden=(hiddenK-priorK),
                         linear.output=F, rep=1, err.fct="sse")
    
    pred = compute(fit.boot, dataT.boot[,-1])$net.result
    
    sse = mean((pred-dataT.boot[,1])^2)
    
    e632 = e632 + 0.632*sse
  }
  
  e632 = e632/B
  
  e632.final = e632 + 0.368 * err
  
  out = list(e632.prior = e632.final.P, e632 = e632.final)
  out
}