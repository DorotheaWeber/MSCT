#######################################
######### prior info neural net #######
#######################################

# 5.3 Pathway models

require(neuralnet)
require(devtools)

# load genes of pathways --> prior info

load("Prior_info_pathways.RData")

load("targetpathwaysCCLE.RData")

# load function for variable importance for nnet

source_url('https://gist.githubusercontent.com/fawda123/6206737/raw/d6f365c283a8cae23fb20892dc223bc5764d50c7/gar_fun.r')

# load boot_neuralnet function

source("boot_neuralnet.r")

#####################################################################################################################

# create data file

compounds = c(2,3,4,8,9,14,16, 17,18, 20, 21, 22,24)

load(paste("IC50_data_", ND[1],".RData", sep=""))
data = ifelse(data=="0", 0, 1)
data = as.data.frame(data)
colnames(data)[which(colnames(data)=="BRAF.MC_MUT")]="BRAF_MUT"

data_all_drugs = data
colnames(data_all_drugs)[1]="AAG"

for(i in compounds){
  
  load(paste("IC50_data_", ND[i],".RData", sep=""))
  data = ifelse(data=="0", 0, 1)
  data = as.data.frame(data)
  
  data_all_drugs = cbind(data[,1], data_all_drugs)
  
  colnames(data_all_drugs)[1]=ND[i]
}


###################################################

# univariate

BrierScore_univ = matrix(0, ncol=14, nrow=10)

for(i in 1:14){

for (j in 1:10){
  
  training = data_all_drugs[-CV_Splits[[i]][[j]], c(i,15:103)]
  test = data_all_drugs[CV_Splits[[i]][[j]], c(i,15:103)]
  
  n=colnames(data_all_drugs)
  
  f <- as.formula(paste(n[i], "~", paste(n[!n %in% ND], collapse = " + ")))
  
  Model_uni = neuralnet(f, training, hidden=80,
                             linear.output=F, rep=200)
  prediction = compute(Model_uni, test[,-1])$net.result
  
  assign(paste("BS",j, sep=""), 
         apply((prediction-test[,1])^2, 2, mean) )
}


BrierScore_univ[,i] = c(BS1 , BS2 ,BS3 , BS4 , BS5 ,BS6 , BS7 , BS8 ,BS9 , BS10)

}

BrierScore_univ_mean = apply(BrierScore_univ, 2, mean)
names(BrierScore_univ_mean) = colnames(data_all_drugs)[1:14]
colnames(BrierScore_univ) = colnames(data_all_drugs)[1:14]

write.table(BrierScore_univ_mean, "BrierScore_univ_mean.txt", sep=";")
write.table(BrierScore_univ, "BrierScore_univ.txt", sep=";")


###################################################

# formula  for neural net

n=colnames(data_all_drugs)

f <- as.formula(paste(paste(n[n %in% ND], collapse = " + "), "~", paste(n[!n %in% ND], collapse = " + ")))

Model_multiple = neuralnet(f, data_all_drugs, hidden=80, exclude=excl, constant.weights=constwts,
                           linear.output=F, rep=200)

save(Model_multiple, file="Model_multiple.RData")


for (j in 1:10){
  
  training = data_all_drugs[-CV_Splits[["PF2341066"]][[j]],]
  test = data_all_drugs[CV_Splits[["PF2341066"]][[j]],]
  
  Model_multiple = neuralnet(f, training, hidden=80,
                                  linear.output=F, rep=200)
  prediction = compute(Model_multiple, test[,15:103])$net.result
  
  assign(paste("BS",j, sep=""), 
         apply((prediction-test[,1:14])^2, 2, mean) )
}

BrierScore_multiple = matrix(c(BS1 , BS2 ,BS3 , BS4 , BS5 ,BS6 , BS7 , BS8 ,BS9 , BS10), ncol=14)
colnames(BrierScore_multiple) = colnames(data_all_drugs)[1:14]
BrierScore_multiple_mean = (BS1 + BS2 + BS3 + BS4 + BS5 + BS6 + BS7 + BS8 + BS9 + BS10)/10
names(BrierScore_multiple_mean) = colnames(data_all_drugs)[1:14]

write.table(BrierScore_multiple_mean, "BrierScore_multiple_mean.txt", sep=";")
write.table(BrierScore_multiple, "BrierScore_multiple.txt", sep=";")


##################################################################

# include pathway info

# exclude matrix

genepathways.ana = list()
predictors = n[!n %in% ND]

for(l in 1:334){
  if(length(which(predictors %in% paste(pathids.genes[pathids.unique[l]][[1]]$hgnc_symbol, "_MUT", sep="")))>3){
    genepathways.ana[[names(pathids.genes[pathids.unique[l]])]] = predictors[which(predictors %in% paste(pathids.genes[pathids.unique[l]][[1]]$hgnc_symbol, "_MUT", sep=""))]
    }
}

genepathways.unique = unique(genepathways.ana)

n.total = length(genepathways.unique[[1]])

excl = matrix(NA, ncol=3, nrow=89-n.total)

# layer
excl[,1] = 1
# which exclude
excl[,2] = seq(1, 89, by=1)[-which(predictors %in% genepathways.unique[[1]])]
# output
excl[,3] = rep(1, (89-n.total))


for(l in 2:22){
  n.total = length(which(predictors %in% genepathways.unique[[l]]))
  
  tmp = matrix(NA, ncol=3, nrow=89-n.total)
  
  # layer
  tmp[,1] = 1
  # which exclude
  tmp[,2] = seq(1, 89, by=1)[-which(predictors %in% genepathways.unique[[l]])]
  # output
  tmp[,3] = rep(l, (89-n.total))
  
  excl=rbind(excl, tmp)
}


# constant weights

constwts = rep(0,dim(excl)[1])

# neuralnet

#Boot_multiple = boot.neuralnet(f, data_all_drugs, hidden = 80, excl, constwts, 10)

Model_multiple = neuralnet(f, data_all_drugs, hidden=80, exclude=excl, constant.weights=constwts,
                              linear.output=F, rep=200)

save(priorModel_multiple, file="priorModel_multiple.RData")

for (j in 1:10){
  
  training = data_all_drugs[-CV_Splits[["PF2341066"]][[j]],]
  test = data_all_drugs[CV_Splits[["PF2341066"]][[j]],]
  
  Model_multiple = neuralnet(f, training, hidden=80,
                                  linear.output=F, rep=200)
  prediction = compute(Model_multiple, test[,15:103])$net.result
  
  assign(paste("BS",j, sep=""), 
         apply((prediction-test[,1:14])^2, 2, mean) )
}

BrierScore_multiple = matrix(c(BS1 , BS2 ,BS3 , BS4 , BS5 ,BS6 , BS7 , BS8 ,BS9 , BS10), ncol=14)
colnames(BrierScore_multiple) = colnames(data_all_drugs)[1:14]
BrierScore_multiple_mean = (BS1 + BS2 + BS3 + BS4 + BS5 + BS6 + BS7 + BS8 + BS9 + BS10)/10
names(BrierScore_multiple_mean) = colnames(data_all_drugs)[1:14]

write.table(BrierScore_multiple_mean, "BrierScore_multiple_mean_pathway.txt", sep=";")
write.table(BrierScore_multiple, "BrierScore_multiple_pathway.txt", sep=";")

#######################################################

BrierScore_univ_prior = matrix(0, ncol=14, nrow=10)

for(i in 1:14){
  
  for (j in 1:10){
    
    training = data_all_drugs[-CV_Splits[[i]][[j]], c(i,15:103)]
    test = data_all_drugs[CV_Splits[[i]][[j]], c(i,15:103)]
    
    n=colnames(data_all_drugs)
    
    f <- as.formula(paste(n[i], "~", paste(n[!n %in% ND], collapse = " + ")))
    
    priorModel_uni = neuralnet(f, training, hidden=80,exclude=excl, constant.weights=constwts,
                          linear.output=F, rep=200)
    prediction = compute(priorModel_uni, test[,-1])$net.result
    
    assign(paste("BS",j, sep=""), 
           apply((prediction-test[,1])^2, 2, mean) )
  }
  
  
  BrierScore_univ_prior[,i] = c(BS1 , BS2 ,BS3 , BS4 , BS5 ,BS6 , BS7 , BS8 ,BS9 , BS10)
  
}

BrierScore_univ_mean_prior = apply(BrierScore_univ_prior, 2, mean)
names(BrierScore_univ_mean_prior) = colnames(data_all_drugs)[1:14]
colnames(BrierScore_univ_prior) = colnames(data_all_drugs)[1:14]

write.table(BrierScore_univ_mean_prior, "BrierScore_univ_mean_prior.txt", sep=";")
write.table(BrierScore_univ_prior, "BrierScore_univ_prior.txt", sep=";")

