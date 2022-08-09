##this script is to try SingleCellNet
##only binary data works
library(singleCellNet)

##load files
args <- commandArgs(T)
message (' - Rscript SingleCellNet.R input_training_meta input_training_exp input_testing_exp output_dir')
stTM_real <- as.character(args[1])
expTMraw_real <- as.character(args[2])
expQuery_real <- as.character(args[3])
input_output_dir <- as.character(args[4])

##read files
message (' - we will load files with count matrix')
##read training
stTM_real <- read.csv(stTM_real,row.names = 1)
expTMraw_real <- read.csv(expTMraw_real,row.names=1)
##read testing
expQuery_real <- read.csv(expQuery_real,row.names = 1)
##add a new col to show the name
stTM_real$cell <- row.names(stTM_real)
stTM_real$cell <- gsub(':','.',stTM_real$cell)
stTM_real$cell <- gsub('-','.',stTM_real$cell)

expTMraw_real <- as.matrix((expTMraw_real>0)+0)
expQuery_real <- as.matrix((expQuery_real>0)+0)

##training using default setting
system.time(class_info<-scn_train(stTrain = stTM_real, expTrain = expTMraw_real, nTopGenes = 10, nRand = 70, nTrees = 1000, nTopGenePairs = 25, dLevel = "cell_type", colName_samp = "cell"))

##prediction
system.time(crPBMC <- scn_predict(class_info[['cnProc']], expQuery_real, nrand = 0))

write.csv(crPBMC,paste0(input_output_dir,'/','opt_prediction.csv'),quote = F)













