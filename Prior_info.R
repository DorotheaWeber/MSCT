#######################################
######### genes same response #########
#######################################

CCLE_melanoma_MUT <- read.table("CCLE_melanoma_MUT.csv", header=T, quote="\"", as.is=T)
rownames(CCLE_melanoma_MUT) = CCLE_melanoma_MUT[,1]

load(paste("IC50_", "data_AAG", ".RData", sep=""))

n <- colnames(data)[-1]
cellL = CCLE_melanoma_drug[rownames(data),1]

same.response = matrix(NA,nrow=dim(CCLE_melanoma_MUT)[1], ncol=89)
colnames(same.response)=n
genes.CCLE = rownames(CCLE_melanoma_MUT)

for(j in 1:89){
  for(i in 1:dim(CCLE_melanoma_MUT)[1]){
    if(all(CCLE_melanoma_MUT[i, cellL]==CCLE_melanoma_MUT[n[j], cellL])){
      same.response[i,j] = rownames(CCLE_melanoma_MUT)[i]
    }
  }
}

tmp=rep(NA, 89)

for(i in 1:89){
  if(length(which(is.na(same.response[,i])==F))>1){
    tmp[i] = T
  }
}

which(tmp==T) # 31, 48

same.response[which(!is.na(same.response[,31])),31]

same.response[which(!is.na(same.response[,48])),48]

CCLE_melanoma_MUT[same.response[which(!is.na(same.response[,31])),31],cellL]

# Genes with same response: GHR, MYH11

genes.analysis = c(n, "MYH11_MUT")
genes.analysis[which(genes.analysis=="BRAF.MC_MUT")]="BRAF_MUT"


########################################
####### prior info pathways ############
########################################

load("targetpathwaysCCLE.RData")

for(j in 1:length(targets)){
  
  tmp = rep(0, length(pathids.genes))
  
  for(i in 1:length(pathids.genes)){
    if(any(pathids.genes[[i]]$hgnc_symbol==targets[j])){
      tmp[i]=1
    }
  }
  
  pathids.genes.tmp = which(tmp==1)
  
  pathids.tmp = names(pathids.names)[pathids.genes.tmp]
  
  all.genes.pathway.tmp = pathids.genes[[pathids.tmp[1]]]$hgnc_symbol
  
  if(length(pathids.genes.tmp)>1){
      for(i in 2:length(pathids.genes.tmp)){
        all.genes.pathway.tmp = c(all.genes.pathway.tmp, pathids.genes[[pathids.tmp[i]]]$hgnc_symbol)
    }
  }
  
  genes.pathway.tmp = unique(all.genes.pathway.tmp)
  
  genes.pathway.tmp = paste(genes.pathway.tmp, "_MUT", sep="")
  
  name = paste("prior.info.", targets[j], sep="")
  assign(name, genes.analysis[which(genes.analysis %in% genes.pathway.tmp)])
  
}

# merge different genes to target

target.uni = unique(names(targets))

prior.info.HSP90 = unique(c(prior.info.HSP90AA1, prior.info.HSP90AB1, prior.info.HSP90B1))
prior.info.HDAC = unique(c(prior.info.HDAC1, prior.info.HDAC2, prior.info.HDAC3,
                           prior.info.HDAC4, prior.info.HDAC5, prior.info.HDAC6,
                           prior.info.HDAC7, prior.info.HDAC8, prior.info.HDAC9,
                           prior.info.HDAC10, prior.info.HDAC11))
prior.info.CDK46 = unique(c(prior.info.CDK4, prior.info.CDK6))
prior.info.gammaS = unique(c(prior.info.PSEN1, prior.info.PSEN2, 
                             prior.info.APH1A, prior.info.APH1B,
                             prior.info.NCSTN))

prior.info = list(gammaS=prior.info.gammaS, CDK46 = prior.info.CDK46,
                  HDAC = prior.info.HDAC, HSP90=prior.info.HSP90,
                  ABL = prior.info.ABL1, BRAF = prior.info.BRAF,
                  EGFR = prior.info.EGFR, KIT = prior.info.KIT,
                  ERBB2 = prior.info.ERBB2, FGFR1 = prior.info.FGFR1,
                  VEGFR1 = prior.info.FLT1, FLT4 = prior.info.FLT4,
                  IGF1R = prior.info.IGF1R, KDR = prior.info.KDR,
                  MEK = prior.info.MAP2K1, MDM2 =prior.info.MDM2,
                  CMET = prior.info.MET, PDGFRB = prior.info.PDGFRB,
                  RAF = prior.info.RAF1, XIAP = prior.info.XIAP,
                  SRC = prior.info.SRC, TUBB = prior.info.TUBB1
                  )

save(prior.info, file="Prior_info_pathways.RData")
