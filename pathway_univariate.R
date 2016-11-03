#######################################
######### prior info nerual net #######
#######################################

# 5.2 Pathway models

require(neuralnet)
require(devtools)

# load genes of pathways --> prior info

load("Prior_info_pathways.RData")

# load function for variable importance for nnet

source_url('https://gist.githubusercontent.com/fawda123/6206737/raw/d6f365c283a8cae23fb20892dc223bc5764d50c7/gar_fun.r')

# load boot_neuralnet function

source("boot_neuralnet.r")

#####################################################################################################################

# AAG
# targets: HSP90

load("IC50_data_AAG.RData")

data = ifelse(data=="0", 0, 1)
data = as.data.frame(data)
colnames(data)[which(colnames(data)=="BRAF.MC_MUT")]="BRAF_MUT"

# formula
n <- colnames(data)
predictors=n[2:90]
f <- as.formula(paste("IC50 ~", paste(n[!n %in% "IC50"], collapse = " + ")))

# exclude matrix

HSP90.gene = unique(ifelse(prior.info$HSP90=="MYH11_MUT","GHR_MUT",prior.info$HSP90))

excl = matrix(NA, ncol=3, nrow=(89-length(HSP90.gene)))

# layer
excl[,1] = 1
# which exclude
excl[,2] = seq(1, 89, by=1)[-which(predictors %in% prior.info$HSP90)]
# output
excl[,3] = 1

# constant weights

constwts = rep(0,55)

Boot_AAG = boot.neuralnet(f, data, excl, constwts, 200)
priorModel_AAG = neuralnet(f, data, hidden=8, exclude=excl, constant.weights=constwts,
                           linear.output=F, rep=200)
for(i in 1:200){
  priorModel_AAG$weights[[i]][[1]]=ifelse(is.na(priorModel_AAG$weights[[i]][[1]]), 0, priorModel_AAG$weights[[i]][[1]])
}

vI_prior_neuralnet_AAG = gar.fun('IC50',priorModel_AAG, bar.plot = F)$rel.imp

Model_AAG = neuralnet(f, data, hidden=8, linear.output=F, rep=200)

vI_neuralnet_AAG = gar.fun('IC50',Model_AAG, bar.plot = F)$rel.imp

knot_weights_prior = 0

for(i in 1:200){
  knot_weights_prior = knot_weights_prior + priorModel_AAG$weights[[i]][[2]]
}

knot_weights = 0

for(i in 1:200){
  knot_weights = knot_weights + Model_AAG$weights[[i]][[2]]
}

prior.result.AAG = list(BrierScore=Boot_AAG, VarImp = list(prior=vI_prior_neuralnet_AAG, all=vI_neuralnet_AAG),
                        KnotWeight = list(all=knot_weights, prior=knot_weights_prior),
                        targets = predictors[which(predictors %in% prior.info$HSP90)],
                        result.weights = list(prior = priorModel_AAG$weights, all =  Model_AAG$weights))

BS_AAG = rep(0,10)
for (j in 1:10){
  not = CV_Splits[[1]][[j]]
  training = data[-which(rownames(data) %in% not), ]
  test = data[not, ]

  priorModel_AAG = neuralnet(f, training, hidden=8, exclude=excl, constant.weights=constwts,
                                                       linear.output=F, rep=200)
  
  prediction = compute(priorModel_AAG, test[,-1])$net.result
  
  BS_AAG[j] = mean((prediction-test[,1])^2)
}

BS_AAG_mean = mean(BS_AAG)


BS_AAG_no = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_AAG_no[j] = mean((prediction-test[,1])^2)
}

BS_AAG_mean_no = mean(BS_AAG_no)

#####################################################################################################################

# PD_0325901
# targets: MEK

load("IC50_data_PD_0325901.RData")

data = ifelse(data=="0", 0, 1)
data = as.data.frame(data)
colnames(data)[which(colnames(data)=="BRAF.MC_MUT")]="BRAF_MUT"

# formula
n <- colnames(data)
predictors=n[2:90]
f <- as.formula(paste("IC50 ~", paste(n[!n %in% "IC50"], collapse = " + ")))

# exclude matrix

MEK.gene = unique(ifelse(prior.info$MEK=="MYH11_MUT","GHR_MUT",prior.info$MEK))

excl = matrix(NA, ncol=3, nrow=(89-length(MEK.gene)))

# layer
excl[,1] = 1
# which exclude
excl[,2] = seq(1, 89, by=1)[-which(predictors %in% prior.info$MEK)]
# output
excl[,3] = 1

# constant weights

constwts = rep(0,64)

Boot_PD_0325901 = boot.neuralnet(f, data, excl, constwts, 200)
priorModel_PD_0325901 = neuralnet(f, data, hidden=8, exclude=excl, constant.weights=constwts,
                           linear.output=F, rep=200)
for(i in 1:200){
  priorModel_PD_0325901$weights[[i]][[1]]=ifelse(is.na(priorModel_PD_0325901$weights[[i]][[1]]), 
                                                 0, priorModel_PD_0325901$weights[[i]][[1]])
}

vI_prior_neuralnet_PD_0325901 = gar.fun('IC50',priorModel_PD_0325901, bar.plot = F)$rel.imp

Model_PD_0325901 = neuralnet(f, data, hidden=8, linear.output=F, rep=200)

vI_neuralnet_PD_0325901 = gar.fun('IC50', Model_PD_0325901, bar.plot = F)$rel.imp

knot_weights_prior = 0

for(i in 1:200){
  knot_weights_prior = knot_weights_prior + priorModel_PD_0325901$weights[[i]][[2]]
}

knot_weights = 0

for(i in 1:200){
  knot_weights = knot_weights + Model_PD_0325901$weights[[i]][[2]]
}

prior.result.PD_0325901 = list(BrierScore=Boot_PD_0325901, 
                               VarImp=list(prior=vI_prior_neuralnet_PD_0325901, all=vI_neuralnet_PD_0325901), 
                               KnotWeight = list(all=knot_weights, prior=knot_weights_prior),
                               targets = predictors[which(predictors %in% prior.info$MEK)],
                               result.weights = list(prior = priorModel_PD_0325901$weights,
                                                     all =  Model_PD_0325901$weights))

BS_PD_0325901 = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8, exclude=excl, constant.weights=constwts,
                             linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_PD_0325901[j] = mean((prediction-test[,1])^2)
}

BS_PD_0325901_mean = mean(BS_PD_0325901)


BS_PD_0325901_no = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_PD_0325901_no[j] = mean((prediction-test[,1])^2)
}

BS_PD_0325901_mean_no = mean(BS_PD_0325901_no)

#####################################################################################################################

# Sorafenib
# targets: BRAF, PDGFRB, KIT, RAF, FLT4, KDR, VEGFR1

# prior.info$KIT = prior.info$PDGFRB
# prior.info$VEGFR1 = prior.info$FLT4

load("IC50_data_Sorafenib.RData")

data = ifelse(data=="0", 0, 1)
data = as.data.frame(data)
colnames(data)[which(colnames(data)=="BRAF.MC_MUT")]="BRAF_MUT"

# formula
n <- colnames(data)
predictors=n[2:90]
f <- as.formula(paste("IC50 ~", paste(n[!n %in% "IC50"], collapse = " + ")))

# exclude matrix

BRAF.gene = unique(ifelse(prior.info$BRAF=="MYH11_MUT","GHR_MUT",prior.info$BRAF))
PDGFRB.gene = unique(ifelse(prior.info$PDGFRB=="MYH11_MUT","GHR_MUT",prior.info$PDGFRB))
RAF.gene = unique(ifelse(prior.info$RAF=="MYH11_MUT","GHR_MUT",prior.info$RAF))
KDR.gene = unique(ifelse(prior.info$KDR=="MYH11_MUT","GHR_MUT",prior.info$KDR))
VEGFR1.gene = unique(ifelse(prior.info$VEGFR1=="MYH11_MUT","GHR_MUT",prior.info$VEGFR1))

n.total = ((89-length(BRAF.gene))+(89-length(PDGFRB.gene))+(89-length(RAF.gene))+(89-length(KDR.gene))+
             (89-length(VEGFR1.gene)))

excl = matrix(NA, ncol=3, nrow=n.total)

# layer
excl[,1] = 1
# which exclude
excl[,2] = c(seq(1, 89, by=1)[-which(predictors %in% prior.info$BRAF)],
             seq(1, 89, by=1)[-which(predictors %in% prior.info$PDGFRB)],
             seq(1, 89, by=1)[-which(predictors %in% prior.info$RAF)],
             seq(1, 89, by=1)[-which(predictors %in% prior.info$KDR)],
             seq(1, 89, by=1)[-which(predictors %in% prior.info$VEGFR1)])
# output
excl[,3] = c(rep(1, 89-length(BRAF.gene)), rep(2, 89-length(PDGFRB.gene)),
             rep(4, 89-length(RAF.gene)), rep(6, 89-length(KDR.gene)), rep(7, 89-length(VEGFR1.gene)))


prior.Sorafenib = unique(c(prior.info$BRAF, prior.info$PDGFRB, prior.info$RAF,
                           prior.info$KDR, prior.info$VEGFR1))
prior.Sorafenib = prior.Sorafenib[-which(prior.Sorafenib=="MYH11_MUT")]

excl2 = matrix(NA, ncol=3, nrow=(89-length(prior.Sorafenib)))

excl2[,1] = 1
# which exclude
excl2[,2] = c(seq(1, 89, by=1)[-which(predictors %in% prior.Sorafenib)])
# output
excl2[,3] = c(rep(1, (89-length(prior.Sorafenib))))

# constant weights

constwts = rep(0,n.total)

constwts2 = rep(0,length(prior.Sorafenib))

# every pathway one hidden knot

Boot_Sorafenib = boot.neuralnet(f, data, excl, constwts, 200)

priorModel_Sorafenib = neuralnet(f, data, hidden=8, exclude=excl, constant.weights=constwts,
                               linear.output=F, rep=200)
for(i in 1:200){
  priorModel_Sorafenib$weights[[i]][[1]]=ifelse(is.na(priorModel_Sorafenib$weights[[i]][[1]]), 
                                              0, priorModel_Sorafenib$weights[[i]][[1]])
}

vI_prior_neuralnet_Sorafenib = gar.fun('IC50',priorModel_Sorafenib, bar.plot = F)$rel.imp

Model_Sorafenib = neuralnet(f, data, hidden=8, linear.output=F, rep=200)

vI_neuralnet_Sorafenib = gar.fun('IC50', Model_Sorafenib, bar.plot = F)$rel.imp

# all pathways in one hidden knot

Boot_Sorafenib_all = boot.neuralnet(f, data, excl2, constwts2, 200)

priorModel_Sorafenib_all = neuralnet(f, data, hidden=8, exclude=excl2, constant.weights=constwts2,
                                   linear.output=F, rep=200)
for(i in 1:200){
  priorModel_Sorafenib_all$weights[[i]][[1]]=ifelse(is.na(priorModel_Sorafenib_all$weights[[i]][[1]]), 
                                                  0, priorModel_Sorafenib_all$weights[[i]][[1]])
}

vI_prior_neuralnet_Sorafenib_all = gar.fun('IC50',priorModel_Sorafenib_all, bar.plot = F)$rel.imp

knot_weights_prior = 0

for(i in 1:200){
  knot_weights_prior = knot_weights_prior + priorModel_Sorafenib$weights[[i]][[2]]
}

knot_weights = 0

for(i in 1:200){
  knot_weights = knot_weights + Model_Sorafenib$weights[[i]][[2]]
}

knot_weights_prior2 = 0

for(i in 1:200){
  knot_weights_prior2 = knot_weights_prior2 + priorModel_Sorafenib_all$weights[[i]][[2]]
}


prior.result.Sorafenib = list(BrierScore=Boot_Sorafenib, 
                               VarImp=list(prior=vI_prior_neuralnet_Sorafenib, all=vI_neuralnet_Sorafenib), 
                               KnotWeight = list(all=knot_weights, prior=knot_weights_prior),
                              targets = list(BRAF = predictors[which(predictors %in% prior.info$BRAF)],
                                             PDGFRB = predictors[which(predictors %in% prior.info$PDGFRB)],
                                             KIT = predictors[which(predictors %in% prior.info$KIT)],
                                             RAF = predictors[which(predictors %in% prior.info$RAF)],
                                             FLT4 = predictors[which(predictors %in% prior.info$FLT4)],
                                             KDR = predictors[which(predictors %in% prior.info$KDR)],
                                             VEGFR1 = predictors[which(predictors %in% prior.info$VEGFR1)]),
                              result.weights = list(prior = priorModel_Sorafenib$weights,
                                                    all =  Model_Sorafenib$weights))

prior.result.Sorafenib.2 = list(BrierScore=Boot_Sorafenib_all, 
                              VarImp=list(prior=vI_prior_neuralnet_Sorafenib_all, all=vI_neuralnet_Sorafenib), 
                              KnotWeight = list(all=knot_weights, prior=knot_weights_prior2),
                              targets = prior.Sorafenib,
                              result.weights = list(prior = priorModel_Sorafenib_all$weights,
                                                    all =  Model_Sorafenib$weights))

BS_Sorafenib = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8, exclude=excl, constant.weights=constwts,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_Sorafenib[j] = mean((prediction-test[,1])^2)
}

BS_Sorafenib_mean = mean(BS_Sorafenib)


BS_Sorafenib_2 = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8, exclude=excl2, constant.weights=constwts2,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_Sorafenib_2[j] = mean((prediction-test[,1])^2)
}

BS_Sorafenib_mean2 = mean(BS_Sorafenib_2)

BS_Sorafenib_no = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_Sorafenib_no[j] = mean((prediction-test[,1])^2)
}

BS_Sorafenib_mean_no = mean(BS_Sorafenib_no)


#####################################################################################################################

# AZD0530
# targets: SRC, ABL, EGFR

load("IC50_data_AZD0530.RData")

data = ifelse(data=="0", 0, 1)
data = as.data.frame(data)
colnames(data)[which(colnames(data)=="BRAF.MC_MUT")]="BRAF_MUT"

# formula
n <- colnames(data)
predictors=n[2:90]
f <- as.formula(paste("IC50 ~", paste(n[!n %in% "IC50"], collapse = " + ")))

# exclude matrix

SRC.gene = unique(ifelse(prior.info$SRC=="MYH11_MUT","GHR_MUT",prior.info$SRC))
ABL.gene = unique(ifelse(prior.info$ABL=="MYH11_MUT","GHR_MUT",prior.info$ABL))
EGFR.gene = unique(ifelse(prior.info$EGFR=="MYH11_MUT","GHR_MUT",prior.info$EGFR))

n.total = ((89-length(SRC.gene))+(89-length(ABL.gene))+(89-length(EGFR.gene)))

excl = matrix(NA, ncol=3, nrow=n.total)

# layer
excl[,1] = 1
# which exclude
excl[,2] = c(seq(1, 89, by=1)[-which(predictors %in% prior.info$SRC)],
             seq(1, 89, by=1)[-which(predictors %in% prior.info$ABL)],
             seq(1, 89, by=1)[-which(predictors %in% prior.info$EGFR)])
# output
excl[,3] = c(rep(1, (89-length(SRC.gene))), rep(2, (89-length(ABL.gene))),
             rep(3, (89-length(EGFR.gene))))

# one prior hidden knot

prior.AZD0530 = unique(c(prior.info$SRC, prior.info$ABL, prior.info$EGFR))
prior.AZD0530 = prior.AZD0530[-which(prior.AZD0530=="MYH11_MUT")]

excl2 = matrix(NA, ncol=3, nrow=(89-length(prior.AZD0530)))

excl2[,1] = 1
# which exclude
excl2[,2] = c(seq(1, 89, by=1)[-which(predictors %in% prior.AZD0530)])
# output
excl2[,3] = c(rep(1, (89-length(prior.AZD0530))))

# constant weights

constwts = rep(0,n.total)

constwts2 = rep(0,length(prior.AZD0530))

# every pathway one hidden knot

Boot_AZD0530 = boot.neuralnet(f, data, excl, constwts, 200)

priorModel_AZD0530 = neuralnet(f, data, hidden=8, exclude=excl, constant.weights=constwts,
                                  linear.output=F, rep=200)
for(i in 1:200){
  priorModel_AZD0530$weights[[i]][[1]]=ifelse(is.na(priorModel_AZD0530$weights[[i]][[1]]), 
                                                 0, priorModel_AZD0530$weights[[i]][[1]])
}

vI_prior_neuralnet_AZD0530 = gar.fun('IC50',priorModel_AZD0530, bar.plot = F)$rel.imp

Model_AZD0530 = neuralnet(f, data, hidden=8, linear.output=F, rep=200)

vI_neuralnet_AZD0530 = gar.fun('IC50', Model_AZD0530, bar.plot = F)$rel.imp

# all pathways in one hidden knot

Boot_AZD0530_all = boot.neuralnet(f, data, excl2, constwts2, 200)

priorModel_AZD0530_all = neuralnet(f, data, hidden=8, exclude=excl2, constant.weights=constwts2,
                               linear.output=F, rep=200)
for(i in 1:200){
  priorModel_AZD0530_all$weights[[i]][[1]]=ifelse(is.na(priorModel_AZD0530_all$weights[[i]][[1]]), 
                                              0, priorModel_AZD0530_all$weights[[i]][[1]])
}

vI_prior_neuralnet_AZD0530_all = gar.fun('IC50',priorModel_AZD0530_all, bar.plot = F)$rel.imp


knot_weights_prior = 0

for(i in 1:200){
  knot_weights_prior = knot_weights_prior + priorModel_AZD0530$weights[[i]][[2]]
}

knot_weights = 0

for(i in 1:200){
  knot_weights = knot_weights + Model_AZD0530$weights[[i]][[2]]
}

knot_weights_prior2 = 0

for(i in 1:200){
  knot_weights_prior2 = knot_weights_prior2 + priorModel_AZD0530_all$weights[[i]][[2]]
}


prior.result.AZD0530 = list(BrierScore=Boot_AZD0530, 
                              VarImp=list(prior=vI_prior_neuralnet_AZD0530, all=vI_neuralnet_AZD0530), 
                              KnotWeight = list(all=knot_weights, prior=knot_weights_prior),
                            targets = list(SRC = predictors[which(predictors %in% prior.info$SRC)],
                                           ABL = predictors[which(predictors %in% prior.info$ABL)],
                                           EGFR = predictors[which(predictors %in% prior.info$EGFR)]),
                            result.weights = list(prior = priorModel_AZD0530$weights,
                                                  all =  Model_AZD0530$weights))

prior.result.AZD0530.2 = list(BrierScore=Boot_AZD0530_all, 
                               VarImp=list(prior=vI_prior_neuralnet_AZD0530_all, all=vI_neuralnet_AZD0530), 
                               KnotWeight = list(all=knot_weights, prior=knot_weights_prior2),
                              targets = prior.AZD0530,
                              result.weights = list(prior = priorModel_AZD0530_all$weights,
                                                    all =  Model_AZD0530$weights))



BS_AZD0530 = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8, exclude=excl, constant.weights=constwts,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_AZD0530[j] = mean((prediction-test[,1])^2)
}

BS_AZD0530_mean = mean(BS_AZD0530)

BS_AZD0530_2 = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8, exclude=excl2, constant.weights=constwts2,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_AZD0530_2[j] = mean((prediction-test[,1])^2)
}

BS_AZD0530_mean2 = mean(BS_AZD0530_2)

BS_AZD0530_no = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_AZD0530_no[j] = mean((prediction-test[,1])^2)
}

BS_AZD0530_mean_no = mean(BS_AZD0530_no)

#####################################################################################################################

# TKI258
# targets: EGFR, KDR, VEGFR1, PDGFRB, FGFR1

load("IC50_data_TKI258.RData")

data = ifelse(data=="0", 0, 1)
data = as.data.frame(data)
colnames(data)[which(colnames(data)=="BRAF.MC_MUT")]="BRAF_MUT"

# formula
n <- colnames(data)
predictors=n[2:90]
f <- as.formula(paste("IC50 ~", paste(n[!n %in% "IC50"], collapse = " + ")))

# exclude matrix

EGFR.gene = unique(ifelse(prior.info$EGFR=="MYH11_MUT","GHR_MUT",prior.info$EGFR))
PDGFRB.gene = unique(ifelse(prior.info$PDGFRB=="MYH11_MUT","GHR_MUT",prior.info$PDGFRB))
FGFR1.gene = unique(ifelse(prior.info$FGFR1=="MYH11_MUT","GHR_MUT",prior.info$FGFR1))
KDR.gene = unique(ifelse(prior.info$KDR=="MYH11_MUT","GHR_MUT",prior.info$KDR))
VEGFR1.gene = unique(ifelse(prior.info$VEGFR1=="MYH11_MUT","GHR_MUT",prior.info$VEGFR1))

n.total = ((89-length(EGFR.gene))+(89-length(PDGFRB.gene))+
             (89-length(FGFR1.gene))+(89-length(KDR.gene))+
             (89-length(VEGFR1.gene)))

excl = matrix(NA, ncol=3, nrow=n.total)

# layer
excl[,1] = 1
# which exclude
excl[,2] = c(seq(1, 89, by=1)[-which(predictors %in% prior.info$EGFR)],
             seq(1, 89, by=1)[-which(predictors %in% prior.info$PDGFRB)],
             seq(1, 89, by=1)[-which(predictors %in% prior.info$FGFR1)],
             seq(1, 89, by=1)[-which(predictors %in% prior.info$KDR)],
             seq(1, 89, by=1)[-which(predictors %in% prior.info$VEGFR1)])
# output
excl[,3] = c(rep(1, (89-length(EGFR.gene))), rep(2, (89-length(PDGFRB.gene))),
             rep(3, (89-length(FGFR1.gene))), rep(4, (89-length(KDR.gene))),
             rep(5, (89-length(VEGFR1.gene))))

# constant weights

constwts = rep(0,n.total)

# one prior hidden knot

prior.TKI258 = unique(c(prior.info$EGFR, prior.info$PDGFRB, prior.info$FGFR1, prior.info$KDR, prior.info$VEGFR1))
prior.TKI258 = prior.TKI258[-which(prior.TKI258=="MYH11_MUT")]

excl2 = matrix(NA, ncol=3, nrow=(89-length(prior.TKI258)))

excl2[,1] = 1
# which exclude
excl2[,2] = c(seq(1, 89, by=1)[-which(predictors %in% prior.TKI258)])
# output
excl2[,3] = c(rep(1, (89-length(prior.TKI258))))


Boot_TKI258 = boot.neuralnet(f, data, excl, constwts, 200)

priorModel_TKI258 = neuralnet(f, data, hidden=8, exclude=excl, constant.weights=constwts,
                               linear.output=F, rep=200)
for(i in 1:200){
  priorModel_TKI258$weights[[i]][[1]]=ifelse(is.na(priorModel_TKI258$weights[[i]][[1]]), 
                                              0, priorModel_TKI258$weights[[i]][[1]])
}

vI_prior_neuralnet_TKI258 = gar.fun('IC50',priorModel_TKI258, bar.plot = F)$rel.imp

Model_TKI258 = neuralnet(f, data, hidden=8, linear.output=F, rep=200)

vI_neuralnet_TKI258 = gar.fun('IC50', Model_TKI258, bar.plot = F)$rel.imp

# all pathways in one knot

Boot_TKI258_all = boot.neuralnet(f, data, excl2, constwts2, 200)

priorModel_TKI258_all = neuralnet(f, data, hidden=8, exclude=excl2, constant.weights=constwts2,
                              linear.output=F, rep=200)
for(i in 1:200){
  priorModel_TKI258_all$weights[[i]][[1]]=ifelse(is.na(priorModel_TKI258_all$weights[[i]][[1]]), 
                                             0, priorModel_TKI258_all$weights[[i]][[1]])
}

vI_prior_neuralnet_TKI258_all = gar.fun('IC50',priorModel_TKI258_all, bar.plot = F)$rel.imp


knot_weights_prior = 0

for(i in 1:200){
  knot_weights_prior = knot_weights_prior + priorModel_TKI258$weights[[i]][[2]]
}

knot_weights = 0

for(i in 1:200){
  knot_weights = knot_weights + Model_TKI258$weights[[i]][[2]]
}

knot_weights_prior2 = 0

for(i in 1:200){
  knot_weights_prior2 = knot_weights_prior2 + priorModel_TKI258_all$weights[[i]][[2]]
}


prior.result.TKI258 = list(BrierScore=Boot_TKI258, 
                            VarImp=list(prior=vI_prior_neuralnet_TKI258, all=vI_neuralnet_TKI258), 
                            KnotWeight = list(all=knot_weights, prior=knot_weights_prior),
                           targets = list(KDR = predictors[which(predictors %in% prior.info$KDR)],
                                          VEGFR1 = predictors[which(predictors %in% prior.info$VEGFR1)],
                                          EGFR = predictors[which(predictors %in% prior.info$EGFR)],
                                          PDGFRB = predictors[which(predictors %in% prior.info$PDGFRB)],
                                          FGFR1 = predictors[which(predictors %in% prior.info$FGFR1)]),
                           result.weights = list(prior = priorModel_TKI258$weights,
                                                 all =  Model_TKI258$weights))

prior.result.TKI258.2 = list(BrierScore=Boot_TKI258_all, 
                              VarImp=list(prior=vI_prior_neuralnet_TKI258_all, all=vI_neuralnet_TKI258), 
                              KnotWeight = list(all=knot_weights, prior=knot_weights_prior2),
                             targets = prior.TKI258,
                             result.weights = list(prior = priorModel_TKI258$weights,
                                                   all =  Model_TKI258$weights))



BS_TKI258 = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8, exclude=excl, constant.weights=constwts,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_TKI258[j] = mean((prediction-test[,1])^2)
}

BS_TKI258_mean = mean(BS_TKI258)


BS_TKI258_2 = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8, exclude=excl2, constant.weights=constwts2,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_TKI258_2[j] = mean((prediction-test[,1])^2)
}

BS_TKI258_mean2 = mean(BS_TKI258_2)

BS_TKI258_no = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_TKI258_no[j] = mean((prediction-test[,1])^2)
}

BS_TKI258_mean_no = mean(BS_TKI258_no)

#####################################################################################################################

# Lapatinib
# targets: EGFR

load("IC50_data_Lapatinib.RData")

data = ifelse(data=="0", 0, 1)
data = as.data.frame(data)
colnames(data)[which(colnames(data)=="BRAF.MC_MUT")]="BRAF_MUT"

# formula
n <- colnames(data)
predictors=n[2:90]
f <- as.formula(paste("IC50 ~", paste(n[!n %in% "IC50"], collapse = " + ")))

# exclude matrix

EGFR.gene = prior.info$EGFR[-which(prior.info$EGFR=="MYH11_MUT")]

n.total = 89-length(EGFR.gene)

excl = matrix(NA, ncol=3, nrow=n.total)

# layer
excl[,1] = 1
# which exclude
excl[,2] = seq(1, 89, by=1)[-which(predictors %in% prior.info$EGFR)]
# output
excl[,3] = 1



# constant weights

constwts = rep(0,n.total)

Boot_Lapatinib = boot.neuralnet(f, data, excl, constwts, 200)

priorModel_Lapatinib = neuralnet(f, data, hidden=8, exclude=excl, constant.weights=constwts,
                              linear.output=F, rep=200)
for(i in 1:200){
  priorModel_Lapatinib$weights[[i]][[1]]=ifelse(is.na(priorModel_Lapatinib$weights[[i]][[1]]), 
                                             0, priorModel_Lapatinib$weights[[i]][[1]])
}

vI_prior_neuralnet_Lapatinib = gar.fun('IC50',priorModel_Lapatinib, bar.plot = F)$rel.imp

Model_Lapatinib = neuralnet(f, data, hidden=8, linear.output=F, rep=200)

vI_neuralnet_Lapatinib = gar.fun('IC50', Model_Lapatinib, bar.plot = F)$rel.imp


knot_weights_prior = 0

for(i in 1:200){
  knot_weights_prior = knot_weights_prior + priorModel_Lapatinib$weights[[i]][[2]]
}

knot_weights = 0

for(i in 1:200){
  knot_weights = knot_weights + Model_Lapatinib$weights[[i]][[2]]
}

prior.result.Lapatinib = list(BrierScore=Boot_Lapatinib, 
                           VarImp=list(prior=vI_prior_neuralnet_Lapatinib, all=vI_neuralnet_Lapatinib), 
                           KnotWeight = list(all=knot_weights, prior=knot_weights_prior),
                           targets = predictors[which(predictors %in% prior.info$EGFR)],
                           result.weights = list(prior = priorModel_Lapatinib$weights, all =  Model_Lapatinib$weights))



BS_Lapatinib = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8, exclude=excl, constant.weights=constwts,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_Lapatinib[j] = mean((prediction-test[,1])^2)
}

BS_Lapatinib_mean = mean(BS_Lapatinib)


BS_Lapatinib_no = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_Lapatinib_no[j] = mean((prediction-test[,1])^2)
}

BS_Lapatinib_mean_no = mean(BS_Lapatinib_no)

#####################################################################################################################

# PHA-665752
# targets: CMET

load("IC50_data_PHA_665752.RData")

data = ifelse(data=="0", 0, 1)
data = as.data.frame(data)
colnames(data)[which(colnames(data)=="BRAF.MC_MUT")]="BRAF_MUT"

# formula
n <- colnames(data)
predictors=n[2:90]
f <- as.formula(paste("IC50 ~", paste(n[!n %in% "IC50"], collapse = " + ")))

# exclude matrix

CMET.gene = prior.info$CMET[-which(prior.info$CMET=="MYH11_MUT")]

n.total = 89-length(CMET.gene)

excl = matrix(NA, ncol=3, nrow=n.total)

# layer
excl[,1] = 1
# which exclude
excl[,2] = seq(1, 89, by=1)[-which(predictors %in% prior.info$CMET)]
# output
excl[,3] = 1



# constant weights

constwts = rep(0,n.total)

Boot_PHA_665752 = boot.neuralnet(f, data, excl, constwts, 200)

priorModel_PHA_665752 = neuralnet(f, data, hidden=8, exclude=excl, constant.weights=constwts,
                                 linear.output=F, rep=200)
for(i in 1:200){
  priorModel_PHA_665752$weights[[i]][[1]]=ifelse(is.na(priorModel_PHA_665752$weights[[i]][[1]]), 
                                                0, priorModel_PHA_665752$weights[[i]][[1]])
}

vI_prior_neuralnet_PHA_665752 = gar.fun('IC50',priorModel_PHA_665752, bar.plot = F)$rel.imp

Model_PHA_665752 = neuralnet(f, data, hidden=8, linear.output=F, rep=200)

vI_neuralnet_PHA_665752 = gar.fun('IC50', Model_PHA_665752, bar.plot = F)$rel.imp


knot_weights_prior = 0

for(i in 1:200){
  knot_weights_prior = knot_weights_prior + priorModel_PHA_665752$weights[[i]][[2]]
}

knot_weights = 0

for(i in 1:200){
  knot_weights = knot_weights + Model_PHA_665752$weights[[i]][[2]]
}

prior.result.PHA_665752 = list(BrierScore=Boot_PHA_665752, 
                              VarImp=list(prior=vI_prior_neuralnet_PHA_665752, all=vI_neuralnet_PHA_665752), 
                              KnotWeight = list(all=knot_weights, prior=knot_weights_prior),
                              targets = predictors[which(predictors %in% prior.info$CMET)],
                              result.weights = list(prior = priorModel_PHA_665752$weights, 
                                                    all =  Model_PHA_665752$weights))


BS_PHA_665752 = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8, exclude=excl, constant.weights=constwts,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_PHA_665752[j] = mean((prediction-test[,1])^2)
}

BS_PHA_665752_mean = mean(BS_PHA_665752)


BS_PHA_665752_no = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_PHA_665752_no[j] = mean((prediction-test[,1])^2)
}

BS_PHA_665752_mean_no = mean(BS_PHA_665752_no)

#####################################################################################################################

# PF-2341066
# targets: CMET

load("IC50_data_PF_2341066.RData")

data = ifelse(data=="0", 0, 1)
data = as.data.frame(data)
colnames(data)[which(colnames(data)=="BRAF.MC_MUT")]="BRAF_MUT"

# formula
n <- colnames(data)
predictors=n[2:90]
f <- as.formula(paste("IC50 ~", paste(n[!n %in% "IC50"], collapse = " + ")))

# exclude matrix

n.total = 89-length(CMET.gene)

excl = matrix(NA, ncol=3, nrow=n.total)

# layer
excl[,1] = 1
# which exclude
excl[,2] = seq(1, 89, by=1)[-which(predictors %in% prior.info$CMET)]
# output
excl[,3] = 1



# constant weights

constwts = rep(0,n.total)

Boot_PF_2341066 = boot.neuralnet(f, data, excl, constwts, 200)

priorModel_PF_2341066 = neuralnet(f, data, hidden=8, exclude=excl, constant.weights=constwts,
                                  linear.output=F, rep=200)
for(i in 1:200){
  priorModel_PF_2341066$weights[[i]][[1]]=ifelse(is.na(priorModel_PF_2341066$weights[[i]][[1]]), 
                                                 0, priorModel_PF_2341066$weights[[i]][[1]])
}

vI_prior_neuralnet_PF_2341066 = gar.fun('IC50',priorModel_PF_2341066, bar.plot = F)$rel.imp

Model_PF_2341066 = neuralnet(f, data, hidden=8, linear.output=F, rep=200)

vI_neuralnet_PF_2341066 = gar.fun('IC50', Model_PF_2341066, bar.plot = F)$rel.imp


knot_weights_prior = 0

for(i in 1:200){
  knot_weights_prior = knot_weights_prior + priorModel_PF_2341066$weights[[i]][[2]]
}

knot_weights = 0

for(i in 1:200){
  knot_weights = knot_weights + Model_PF_2341066$weights[[i]][[2]]
}

prior.result.PF_2341066 = list(BrierScore=Boot_PF_2341066, 
                               VarImp=list(prior=vI_prior_neuralnet_PF_2341066, all=vI_neuralnet_PF_2341066), 
                               KnotWeight = list(all=knot_weights, prior=knot_weights_prior),
                               targets = predictors[which(predictors %in% prior.info$CMET)],
                               result.weights = list(prior = priorModel_PF_2341066$weights, 
                                                     all =  Model_PF_2341066$weights))


BS_PF_2341066 = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8, exclude=excl, constant.weights=constwts,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_PF_2341066[j] = mean((prediction-test[,1])^2)
}

BS_PF_2341066_mean = mean(BS_PF_2341066)


BS_PF_2341066_no = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_PF_2341066_no[j] = mean((prediction-test[,1])^2)
}

BS_PF_2341066_mean_no = mean(BS_PF_2341066_no)

#####################################################################################################################

# PLX4720
# targets: RAF

load("IC50_data_PLX4720.RData")

data = ifelse(data=="0", 0, 1)
data = as.data.frame(data)
colnames(data)[which(colnames(data)=="BRAF.MC_MUT")]="BRAF_MUT"

# formula
n <- colnames(data)
predictors=n[2:90]
f <- as.formula(paste("IC50 ~", paste(n[!n %in% "IC50"], collapse = " + ")))

# exclude matrix

n.total = 89-length(RAF.gene)

excl = matrix(NA, ncol=3, nrow=n.total)

# layer
excl[,1] = 1
# which exclude
excl[,2] = seq(1, 89, by=1)[-which(predictors %in% prior.info$RAF)]
# output
excl[,3] = 1

# constant weights

constwts = rep(0,n.total)

Boot_PLX4720 = boot.neuralnet(f, data, excl, constwts, 200)

priorModel_PLX4720 = neuralnet(f, data, hidden=8, exclude=excl, constant.weights=constwts,
                                  linear.output=F, rep=200)
for(i in 1:200){
  priorModel_PLX4720$weights[[i]][[1]]=ifelse(is.na(priorModel_PLX4720$weights[[i]][[1]]), 
                                                 0, priorModel_PLX4720$weights[[i]][[1]])
}

vI_prior_neuralnet_PLX4720 = gar.fun('IC50',priorModel_PLX4720, bar.plot = F)$rel.imp

Model_PLX4720 = neuralnet(f, data, hidden=8, linear.output=F, rep=200)

vI_neuralnet_PLX4720 = gar.fun('IC50', Model_PLX4720, bar.plot = F)$rel.imp


knot_weights_prior = 0

for(i in 1:200){
  knot_weights_prior = knot_weights_prior + priorModel_PLX4720$weights[[i]][[2]]
}

knot_weights = 0

for(i in 1:200){
  knot_weights = knot_weights + Model_PLX4720$weights[[i]][[2]]
}

prior.result.PLX4720 = list(BrierScore=Boot_PLX4720, 
                            VarImp=list(prior=vI_prior_neuralnet_PLX4720, all=vI_neuralnet_PLX4720), 
                            KnotWeight = list(all=knot_weights, prior=knot_weights_prior),
                            targets = predictors[which(predictors %in% prior.info$RAF)],
                            result.weights = list(prior = priorModel_PLX4720$weights, all =  Model_PLX4720$weights))



BS_PLX4720 = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8, exclude=excl, constant.weights=constwts,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_PLX4720[j] = mean((prediction-test[,1])^2)
}

BS_PLX4720_mean = mean(BS_PLX4720)


BS_PLX4720_no = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_PLX4720_no[j] = mean((prediction-test[,1])^2)
}

BS_PLX4720_mean_no = mean(BS_PLX4720_no)


#####################################################################################################################

# AZD6244
# targets: MEK

load("IC50_data_AZD6244.RData")

data = ifelse(data=="0", 0, 1)
data = as.data.frame(data)
colnames(data)[which(colnames(data)=="BRAF.MC_MUT")]="BRAF_MUT"

# formula
n <- colnames(data)
predictors=n[2:90]
f <- as.formula(paste("IC50 ~", paste(n[!n %in% "IC50"], collapse = " + ")))

# exclude matrix

n.total = 89-length(MEK.gene)
excl = matrix(NA, ncol=3, nrow=n.total)

# layer
excl[,1] = 1
# which exclude
excl[,2] = seq(1, 89, by=1)[-which(predictors %in% prior.info$MEK)]
# output
excl[,3] = 1

# constant weights

constwts = rep(0,n.total)

Boot_AZD6244 = boot.neuralnet(f, data, excl, constwts, 200)

priorModel_AZD6244 = neuralnet(f, data, hidden=8, exclude=excl, constant.weights=constwts,
                               linear.output=F, rep=200)
for(i in 1:200){
  priorModel_AZD6244$weights[[i]][[1]]=ifelse(is.na(priorModel_AZD6244$weights[[i]][[1]]), 
                                              0, priorModel_AZD6244$weights[[i]][[1]])
}

vI_prior_neuralnet_AZD6244 = gar.fun('IC50',priorModel_AZD6244, bar.plot = F)$rel.imp

Model_AZD6244 = neuralnet(f, data, hidden=8, linear.output=F, rep=200)

vI_neuralnet_AZD6244 = gar.fun('IC50', Model_AZD6244, bar.plot = F)$rel.imp

knot_weights_prior = 0

for(i in 1:200){
  knot_weights_prior = knot_weights_prior + priorModel_AZD6244$weights[[i]][[2]]
}

knot_weights = 0

for(i in 1:200){
  knot_weights = knot_weights + Model_AZD6244$weights[[i]][[2]]
}

prior.result.AZD6244 = list(BrierScore=Boot_AZD6244, 
                            VarImp=list(prior=vI_prior_neuralnet_AZD6244, all=vI_neuralnet_AZD6244), 
                            KnotWeight = list(all=knot_weights, prior=knot_weights_prior),
                            targets = predictors[which(predictors %in% prior.info$MEK)],
                            result.weights = list(prior = priorModel_AZD6244$weights, all =  Model_AZD6244$weights))


BS_AZD6244 = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8, exclude=excl, constant.weights=constwts,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_PLX4720[j] = mean((prediction-test[,1])^2)
}

BS_AZD6244_mean = mean(BS_AZD6244)


BS_AZD6244_no = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_AZD6244_no[j] = mean((prediction-test[,1])^2)
}

BS_AZD6244_mean_no = mean(BS_AZD6244_no)

#####################################################################################################################

# ZD6474
# targets: VEGFR1

load("IC50_data_ZD_6474.RData")

data = ifelse(data=="0", 0, 1)
data = as.data.frame(data)
colnames(data)[which(colnames(data)=="BRAF.MC_MUT")]="BRAF_MUT"

# formula
n <- colnames(data)
predictors=n[2:90]
f <- as.formula(paste("IC50 ~", paste(n[!n %in% "IC50"], collapse = " + ")))

# exclude matrix

n.total = 89-length(VEGFR1.gene)
excl = matrix(NA, ncol=3, nrow=n.total)

# layer
excl[,1] = 1
# which exclude
excl[,2] = seq(1, 89, by=1)[-which(predictors %in% prior.info$VEGFR1)]
# output
excl[,3] = 1

# constant weights

constwts = rep(0,n.total)

Boot_ZD6474 = boot.neuralnet(f, data, excl, constwts, 200)

priorModel_ZD6474 = neuralnet(f, data, hidden=8, exclude=excl, constant.weights=constwts,
                               linear.output=F, rep=200)
for(i in 1:200){
  priorModel_ZD6474$weights[[i]][[1]]=ifelse(is.na(priorModel_ZD6474$weights[[i]][[1]]), 
                                              0, priorModel_ZD6474$weights[[i]][[1]])
}

vI_prior_neuralnet_ZD6474 = gar.fun('IC50',priorModel_ZD6474, bar.plot = F)$rel.imp

Model_ZD6474 = neuralnet(f, data, hidden=8, linear.output=F, rep=200)

vI_neuralnet_ZD6474 = gar.fun('IC50', Model_ZD6474, bar.plot = F)$rel.imp

knot_weights_prior = 0

for(i in 1:200){
  knot_weights_prior = knot_weights_prior + priorModel_ZD6474$weights[[i]][[2]]
}

knot_weights = 0

for(i in 1:200){
  knot_weights = knot_weights + Model_ZD6474$weights[[i]][[2]]
}

prior.result.ZD6474 = list(BrierScore=Boot_ZD6474, 
                            VarImp=list(prior=vI_prior_neuralnet_ZD6474, all=vI_neuralnet_ZD6474), 
                            KnotWeight = list(all=knot_weights, prior=knot_weights_prior),
                           targets = predictors[which(predictors %in% prior.info$VEGFR1)],
                           result.weights = list(prior = priorModel_ZD6474$weights, all =  Model_ZD6474$weights))


BS_ZD6474 = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8, exclude=excl, constant.weights=constwts,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_ZD6474[j] = mean((prediction-test[,1])^2)
}

BS_ZD6474_mean = mean(BS_ZD6474)


BS_ZD6474_no = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_ZD6474_no[j] = mean((prediction-test[,1])^2)
}

BS_ZD6474_mean_no = mean(BS_ZD6474_no)

#####################################################################################################################

# AEW541
# targets: IGF1R

load("IC50_data_AEW541.RData")

data = ifelse(data=="0", 0, 1)
data = as.data.frame(data)
colnames(data)[which(colnames(data)=="BRAF.MC_MUT")]="BRAF_MUT"

# formula
n <- colnames(data)
predictors=n[2:90]
f <- as.formula(paste("IC50 ~", paste(n[!n %in% "IC50"], collapse = " + ")))

# exclude matrix

IGF1R.gene = prior.info$IGF1R[-which(prior.info$IGF1R=="MYH11_MUT")]

n.total = 89-length(IGF1R.gene)
excl = matrix(NA, ncol=3, nrow=n.total)

# layer
excl[,1] = 1
# which exclude
excl[,2] = seq(1, 89, by=1)[-which(predictors %in% prior.info$IGF1R)]
# output
excl[,3] = 1

# constant weights

constwts = rep(0,n.total)

Boot_AEW541 = boot.neuralnet(f, data, excl, constwts, 200)

priorModel_AEW541 = neuralnet(f, data, hidden=8, exclude=excl, constant.weights=constwts,
                              linear.output=F, rep=200)
for(i in 1:200){
  priorModel_AEW541$weights[[i]][[1]]=ifelse(is.na(priorModel_AEW541$weights[[i]][[1]]), 
                                             0, priorModel_AEW541$weights[[i]][[1]])
}

vI_prior_neuralnet_AEW541 = gar.fun('IC50',priorModel_AEW541, bar.plot = F)$rel.imp

Model_AEW541 = neuralnet(f, data, hidden=8, linear.output=F, rep=200)

vI_neuralnet_AEW541 = gar.fun('IC50', Model_AEW541, bar.plot = F)$rel.imp


knot_weights_prior = 0

for(i in 1:200){
  knot_weights_prior = knot_weights_prior + priorModel_AEW541$weights[[i]][[2]]
}

knot_weights = 0

for(i in 1:200){
  knot_weights = knot_weights + Model_AEW541$weights[[i]][[2]]
}

prior.result.AEW541 = list(BrierScore=Boot_AEW541, 
                           VarImp=list(prior=vI_prior_neuralnet_AEW541, all=vI_neuralnet_AEW541), 
                           KnotWeight = list(all=knot_weights, prior=knot_weights_prior),
                           targets = predictors[which(predictors %in% prior.info$IGF1R)],
                           result.weights = list(prior = priorModel_AEW541$weights, all =  Model_AEW541$weights))


BS_AEW541 = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8, exclude=excl, constant.weights=constwts,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_AEW541[j] = mean((prediction-test[,1])^2)
}

BS_AEW541_mean = mean(BS_AEW541)


BS_AEW541_no = rep(0,10)
for (j in 1:10){
  not = c(((j-1)*4+1):(j*4))
  training = data[-not, ]
  test = data[not, ]
  
  priorModel = neuralnet(f, training, hidden=8,
                         linear.output=F, rep=200)
  
  prediction = compute(priorModel, test[,-1])$net.result
  
  BS_AEW541_no[j] = mean((prediction-test[,1])^2)
}

BS_AEW541_mean_no = mean(BS_AEW541_no)

#####################################################################################################################

# save results of prior analysis

results_prior = list(AEW541 = prior.result.AEW541,
                     AAG = prior.result.AAG,
                     Lapatinib = prior.result.Lapatinib,
                     PHA_665752 = prior.result.PHA_665752,
                     PF_2341066 = prior.result.PF_2341066,
                     AZD0530 = prior.result.AZD0530,
                     AZD0530_1knot = prior.result.AZD0530.2,
                     Sorafenib = prior.result.Sorafenib,
                     Sorafenib_1knot = prior.result.Sorafenib.2,
                     TKI258 = prior.result.TKI258,
                     TKI258_1knot = prior.result.TKI258.2,
                     PLX4720 = prior.result.PLX4720,
                     PD_0325901 = prior.result.PD_0325901,
                     AZD6244 = prior.result.AZD6244,
                     ZD6474 = prior.result.ZD6474
                     )

# save results

save(results_prior, file="Pathway_prior_results.RData")

############################################################################

for(i in 1:15){
  print(names(results_prior)[i])
  print(results_prior[[i]]$BrierScore)
}

# Plots




############################################################################################################

##### weight of target genes

dif.weight.AEW541 = cbind(prior = results_prior$AEW541$VarImp$prior[results_prior$AEW541$targets],
                          all= results_prior$AEW541$VarImp$all[results_prior$AEW541$targets])

dif.weight.AAG = cbind(prior=results_prior$AAG$VarImp$prior[results_prior$AAG$targets],
                          all=results_prior$AAG$VarImp$all[results_prior$AAG$targets])

dif.weight.Lapatinib = cbind(prior=results_prior$Lapatinib$VarImp$prior[results_prior$Lapatinib$targets],
                             all=results_prior$Lapatinib$VarImp$all[results_prior$Lapatinib$targets])

dif.weight.PD_0325901 = cbind(prior=results_prior$PD_0325901$VarImp$prior[results_prior$PD_0325901$targets],
                              all=results_prior$PD_0325901$VarImp$all[results_prior$PD_0325901$targets])

dif.weight.PHA_665752 = cbind(prior=results_prior$PHA_665752$VarImp$prior[results_prior$PHA_665752$targets],
                              all=results_prior$PHA_665752$VarImp$all[results_prior$PHA_665752$targets])

dif.weight.PF_2341066 = cbind(prior=results_prior$PF_2341066$VarImp$prior[results_prior$PF_2341066$targets],
                              all=results_prior$PF_2341066$VarImp$all[results_prior$PF_2341066$targets])

dif.weight.AZD0530 = cbind(prior=results_prior$AZD0530_1knot$VarImp$prior[results_prior$AZD0530_1knot$targets],
                           all=results_prior$AZD0530_1knot$VarImp$all[results_prior$AZD0530_1knot$targets])

dif.weight.Sorafenib = cbind(prior=results_prior$Sorafenib_1knot$VarImp$prior[results_prior$Sorafenib_1knot$targets],
                             all=results_prior$Sorafenib_1knot$VarImp$all[results_prior$Sorafenib_1knot$targets])

dif.weight.TKI258 = cbind(prior=results_prior$TKI258_1knot$VarImp$prior[results_prior$TKI258_1knot$targets],
                          all=results_prior$TKI258_1knot$VarImp$all[results_prior$TKI258_1knot$targets])

dif.weight.PLX4720 = cbind(prior=results_prior$PLX4720$VarImp$prior[results_prior$PLX4720$targets],
                           all=results_prior$PLX4720$VarImp$all[results_prior$PLX4720$targets])

dif.weight.AZD6244 = cbind(prior=results_prior$AZD6244$VarImp$prior[results_prior$AZD6244$targets],
                          all=results_prior$AZD6244$VarImp$all[results_prior$AZD6244$targets])

dif.weight.ZD6474 = cbind(prior=results_prior$ZD6474$VarImp$prior[results_prior$ZD6474$targets],
                           all=results_prior$ZD6474$VarImp$all[results_prior$ZD6474$targets])

target.weights = list(AEW541 = dif.weight.AEW541,
                      AAG = dif.weight.AAG,
                      Lapatinib = dif.weight.Lapatinib,
                      PHA_665752 = dif.weight.PHA_665752,
                      PF_2341066 = dif.weight.PF_2341066,
                      AZD0530 = dif.weight.AZD0530,
                      Sorafenib = dif.weight.Sorafenib,
                      TKI258 = dif.weight.TKI258,
                      PLX4720 = dif.weight.PLX4720,
                      PD_0325901 = dif.weight.PD_0325901,
                      AZD6244 = dif.weight.AZD6244,
                      ZD6474 = dif.weight.ZD6474
)

save(target.weights, file="Pathway_target_weights.RData")

BS_all = c(BS_AAG_mean, BS_AAG_mean_no,
           BS_TKI258_mean, BS_TKI258_mean2, BS_TKI258_mean_no, 
           BS_Sorafenib_mean, BS_Sorafenib_mean2, BS_Sorafenib_mean_no,
           BS_PD_0325901_mean, BS_PD_0325901_mean_no,
           BS_AZD0530_mean, BS_AZD0530_mean2, BS_AZD0530_mean_no,
           BS_ZD6474_mean, BS_ZD6474_mean_no,
           BS_PLX4720_mean, BS_PLX4720_mean_no,
           BS_PHA_665752_mean, BS_PHA_665752_mean_no,
           BS_Lapatinib_mean, BS_Lapatinib_mean_no,
           BS_AEW541_mean, BS_AEW541_mean_no,
           BS_AZD6244_mean, BS_AZD6244_mean_no,
           BS_PF_2341066_mean, BS_PF_2341066_mean_no)
names(BS_all) = c("BS_AAG_mean", "BS_AAG_mean_no",
                  "BS_TKI258_mean", "BS_TKI258_mean2", "BS_TKI258_mean_no", 
                  "BS_Sorafenib_mean", "BS_Sorafenib_mean2", "BS_Sorafenib_mean_no",
                  "BS_PD_0325901_mean", "BS_PD_0325901_mean_no",
                  "BS_AZD0530_mean", "BS_AZD0530_mean2", "BS_AZD0530_mean_no",
                  "BS_ZD6474_mean", "BS_ZD6474_mean_no",
                  "BS_PLX4720_mean", "BS_PLX4720_mean_no",
                  "BS_PHA_665752_mean", "BS_PHA_665752_mean_no",
                  "BS_Lapatinib_mean", "BS_Lapatinib_mean_no",
                  "BS_AEW541_mean", "BS_AEW541_mean_no",
                  "BS_AZD6244_mean", "BS_AZD6244_mean_no",
                  "BS_PF_2341066_mean", "BS_PF_2341066_mean_no")
save(BS_all, file="Pathway_target_BS.RData")

load("Pathway_target_weights.RData")