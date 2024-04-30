# Install repotools (need to be done only once)
source('http://renozao.github.io/repotools/install.R')

library(repotools)
# latest development version
install.pkgs('NMF')


# ======================================================
# Library
# ======================================================
library(tidyverse)
library(NMF)
library(Seurat)
library(Matrix)
library(topGO)
library(reshape2)
library(pheatmap)
library(RColorBrewer)

# ======================================================
# Function to perform quality control
# ======================================================
data_process= function(x){
  x = as.data.frame(x)
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
# ======================================================
# Function to select rank
# ======================================================
rank_selection = function(matrix,a,b){
  estimate = nmf(matrix, rank=a:b, method="brunet", nrun=10, seed=123)
  diff = c(diff(estimate$measures$cophenetic),NA)
  coph = data.frame(estimate$measures$rank, diff)
  coph %>%
    filter(diff < 0) %>%
    head(n=1) %>%
    pull(estimate.measures.rank)
}

