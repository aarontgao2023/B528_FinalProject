# setup
library(tidyverse)
library(NMF)
library(Seurat)
library(Matrix)
data_process= function(x){
  # removed cells less than 2,500 genes were detected
  x = x[colSums(x>0)>2500]
  # remove cells that below threshold
  threshold = mean(colSums(x>0))-2*sd(colSums(x>0))
  x = x[colSums(x>0)>threshold]
  # excluded genes with Ea<4
  tpm = (2^x-1)*10
  Ea=log2(rowMeans(replace(tpm, tpm == 0, NA), na.rm = TRUE)+1)
  Ea = as.data.frame(Ea)
  row_index=which(Ea$Ea>4)
  x = x[row_index,]
  # centering
  x=sweep(x,1,rowMeans(x))
}

# data import
# 1287 male
d='/Users/tianchuangao/Desktop/kidney/1287'
genefile <- paste0(d,"/features.tsv.gz")
barcodefile <- paste0(d,"/barcodes.tsv.gz")
matrixfile <- paste0(d,"/matrix.mtx.gz")   
gn <- read.table(genefile,as.is=T)
cn <- readLines(barcodefile)
matrix_1287 <- as.matrix(readMM(matrixfile))
rownames(matrix_1287) = gn$V2
colnames(matrix_1287) = cn
matrix_1287_after = data_process(as.data.frame(matrix_1287))
matrix_1287_after[matrix_1287_after<0] = 0.000000001
save(matrix_1287_after, file = '/Users/tianchuangao/Desktop/kidney/1287.RData')

# 0606 female
d='/Users/tianchuangao/Desktop/kidney/0606/1'
genefile <- paste0(d,"/features.tsv.gz")
barcodefile <- paste0(d,"/barcodes.tsv.gz")
matrixfile <- paste0(d,"/matrix.mtx.gz")   
gn <- read.table(genefile,as.is=T)
cn <- readLines(barcodefile)
matrix_0606_1 <- as.matrix(readMM(matrixfile))
rownames(matrix_0606_1) = gn$V2
colnames(matrix_0606_1) = cn
matrix_0606_1_after = data_process(as.data.frame(matrix_0606_1))
matrix_0606_1_after[matrix_0606_1_after<0] = 0.000000001
save(matrix_0606_1_after, file = '/Users/tianchuangao/Desktop/kidney/0606_1.RData')

d='/Users/tianchuangao/Desktop/kidney/0606/2'
genefile <- paste0(d,"/features.tsv.gz")
barcodefile <- paste0(d,"/barcodes.tsv.gz")
matrixfile <- paste0(d,"/matrix.mtx.gz")   
gn <- read.table(genefile,as.is=T)
cn <- readLines(barcodefile)
matrix_0606_2 <- as.matrix(readMM(matrixfile))
rownames(matrix_0606_2) = gn$V2
colnames(matrix_0606_2) = cn
matrix_0606_2_after = data_process(as.data.frame(matrix_0606_2))
matrix_0606_2_after[matrix_0606_2_after<0] = 0.000000001
save(matrix_0606_2_after, file = '/Users/tianchuangao/Desktop/kidney/0606_2.RData')

d='/Users/tianchuangao/Desktop/kidney/0606/3'
genefile <- paste0(d,"/features.tsv.gz")
barcodefile <- paste0(d,"/barcodes.tsv.gz")
matrixfile <- paste0(d,"/matrix.mtx.gz")   
gn <- read.table(genefile,as.is=T)
cn <- readLines(barcodefile)
matrix_0606_3 <- as.matrix(readMM(matrixfile))
rownames(matrix_0606_3) = gn$V2
colnames(matrix_0606_3) = cn
matrix_0606_3_after = data_process(as.data.frame(matrix_0606_3))
matrix_0606_3_after[matrix_0606_3_after<0] = 0.000000001
save(matrix_0606_3_after, file = '/Users/tianchuangao/Desktop/kidney/0606_3.RData')


matrix_0601 = cbind(matrix_0606_1,matrix_0606_2)
# nmf
load('1287.RData')
estimate_1 = nmf(matrix_1287_after, rank=2:8, method="brunet", nrun=10, seed=123)
result = nmf(matrix_1287_after,rank=3)
fs_my = extractFeatures(result,30L)
fs_my = lapply(fs_my,function(x)rownames(result)[x])
fs_my

load('0606_1.RData')
estimate_0606_1 = nmf(matrix_0606_1_after, rank=2:6, method="brunet", nrun=10, seed=123)
plot(estimate_0606_1)
result_0606_1 = nmf(matrix_0606_1_after,rank=3)
fs_0606_1 = extractFeatures(result_0606_1,30L)
fs_0606_1 = lapply(fs_0606_1,function(x)rownames(result_0606_1)[x])
fs_0606_1


