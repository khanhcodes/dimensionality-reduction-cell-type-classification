## sample ##


##udpating 080122 only conduct projection

##this script is to do the embeddings of cells and features between training and independent testing dataset 

##updating 032522 this script is to check the scPred

##1. we need to obtain the reference embedding first
#1) we first use all the train set (90%) the other 10% is the independent testing dataset to do the training
#2) now we already have train set (PCs by cells) and train cell meta

##redo the processing
##1) all data (gene by cells) divided into independent and training
##2) transfer the train set to be Seurat object without any normalization
##3) use the Seurat to do the reducation to obtain the gene and cell embeddings respectively
##4) use the cell embeddings to have a training
##5) use the gene embedings to make a new embedings for the independent testing datasets.
##6) do the prediction

##load library
library("Seurat")
library("magrittr")
library('matrixStats')
library('Matrix')
library('harmony')

args = commandArgs(T)

input_training_cell_by_PC_fl <- as.character(args[1])
##provide a file that has been conducted dimension reduction, which is provided by KN

input_training_gene_by_PC_fl <- as.character(args[2])
##provide a file that has been conducted dimenstion reduction, which is provided by KN

input_training_object_fl <- as.character(args[3])
##gene by cell
##'/scratch/hy17471/rice_altlas_scATAC_seq_042021/08_training_classifiers_031822/02_SPmarker_032322/working_dir_topCTver/Step1_prepare_data_dir/Step1_3_split_data_dir/output_dir/step1_split_train_indetest_opt_dir/opt_exp_train.csv'
##we need this one that is required for the projection
##this object contains the cell by pc assays that have been normalized and scaled, which is provided by HD

input_testing_rds_fl <- as.character(args[4])
##it is an object file (gene by cell)
##we need to make sure this testing file has been normalized but do not need to be scaled since it would be scaled in the following script based on the mean and sd of the training dataset, which should be done by HD
##'/scratch/hy17471/rice_altlas_scATAC_seq_042021/08_training_classifiers_031822/02_SPmarker_032322/working_dir_topCTver/Step1_prepare_data_dir/Step1_3_split_data_dir/output_dir/step1_split_train_indetest_opt_dir/opt_exp_test.csv'

input_output_dir <- as.character(args[5])

#####################
##set some parameters
max.iter.harmony = 20
seed=123

###########
##read data
message ('- load files')
cellEmbeddings_ref <- t(read.csv(input_training_cell_by_PC_fl,row.names = 1))
##read reference dimension reduction file

reference_obj <- readRDS(input_training_object_fl)
##this obj is normalized and scaled

loadings_ref <- as.matrix(read.csv(input_training_gene_by_PC_fl,row.names = 1))
##read gene by pc matrix file

testing_dt <- readRDS(input_testing_rds_fl)
##read the testing matrix file

##############
##process data
message (' - process to obtain testing dataset')

##correspond feature name
message (' - correspond features')
rownames(testing_dt) <- gsub('_','-',rownames(testing_dt))
shared_ID_test_ref <- intersect(rownames(testing_dt),rownames(loadings_ref))

##data process
message (' - get mean and std')
assay_ref <- GetAssay(reference_obj)
means_ref <- rowMeans(assay_ref)

assay_mtx_ref <- as.matrix(GetAssay(reference_obj[,])[,]) 
std_ref <- rowSds(assay_mtx_ref)

names(std_ref) <- rownames(assay_ref) -> names(std_ref)

means_ref <- means_ref[shared_ID_test_ref]
std_ref <- std_ref[shared_ID_test_ref]

loadings_ref <- loadings_ref[shared_ID_test_ref, ]

testing_data <- as.matrix(testing_dt[shared_ID_test_ref,])

message (' - scale testing data')
testing_data <- Matrix::t(testing_data)
scaled_testing_data <- scale(testing_data, means_ref, std_ref) ##row is cell and col is gene

##keep in mind the loadings_ref should be matrix
message (' - project to reference')
test_embeddings <- scaled_testing_data %*% as.matrix(loadings_ref)

##create a dataset of reference and testing
message (' - create dataset for reference and testing')
dataset <- factor(c(rep("reference", nrow(cellEmbeddings_ref)), rep("test", nrow(test_embeddings))), 
                  levels = c("reference", "test"))
rownames(cellEmbeddings_ref) <- paste0("ref_", rownames(cellEmbeddings_ref))
rownames(test_embeddings) <- paste0("test_", rownames(test_embeddings))
eigenspace <- as.data.frame(rbind(cellEmbeddings_ref, test_embeddings))
#> head(eigenspace)
#PC_1        PC_2       PC_3       PC_4
#ref_CCACGGACATGCTGGC-1   3.5462928  4.24607540  0.7342539  0.4024824
#ref_GATTCAGGTCTCCACT-1   5.1154612 -1.93469009 -2.5112122 -0.1026473
#ref_TGTATTCTCTATGTGG-1 -32.3288611 -0.03741952 -4.1121667  7.3990149
#ref_GTCTCGTGTTCACGGC-1  -0.9768477 -6.43608814 14.1769009  0.2691508
#ref_TATCAGGCACCAACCG-1   0.7108476  8.61825733  0.8995970 -0.1301161
#ref_TTTGCGCAGTAAGTAC-1   5.6391219 -2.81575179 -2.1088155 -0.2297144

meta_data <- data.frame(rownames(eigenspace), dataset = dataset)
#> head(meta_data)
#rownames.eigenspace.   dataset
#1 ref_CCACGGACATGCTGGC-1 reference
#2 ref_GATTCAGGTCTCCACT-1 reference
#3 ref_TGTATTCTCTATGTGG-1 reference
#4 ref_GTCTCGTGTTCACGGC-1 reference
#5 ref_TATCAGGCACCAACCG-1 reference
#6 ref_TTTGCGCAGTAAGTAC-1 reference

##do the harmony embeddings for the testing dataset
set.seed(seed)
message (' - do the harmony')
harmony_embeddings <- HarmonyMatrix(eigenspace, 
                                    meta_data, 
                                    'dataset', 
                                    do_pca = FALSE, 
                                    reference_values = "reference",
                                    max.iter.harmony = max.iter.harmony)

test_embeddings_aligned <- harmony_embeddings[dataset == "test", , drop = FALSE]

##save the training and testing embeddings
##remove the 'ref_' and 'test_'
test_embeddings_aligned_t <- t(test_embeddings_aligned)
colnames(test_embeddings_aligned_t) <- gsub('test_','',colnames(test_embeddings_aligned_t))
write.csv(test_embeddings_aligned_t,paste0(input_output_dir,'/opt_indep_testing_embeddings.csv'),quote = F)















