
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
library(scales)

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
rank_selection = function(matrix,a,b,nrun = 10, seed = 123){
  estimate = nmf(matrix, rank=a:b, method="brunet", nrun, seed)
  diff_coph = c(diff(estimate$measures$cophenetic),NA)
  coph = data.frame(estimate$measures$rank, diff_coph)
  coph = coph %>%
    filter(diff_coph < 0) %>%
    head(n=1) %>%
    pull(estimate.measures.rank)
  diff_rss = c(diff(diff(estimate$measures$rss)))
  diff_diff_rss = c(NA,NA,diff(diff_rss),NA)
  rss = data.frame(estimate$measures$rank, diff_diff_rss)
  rss = rss %>%
    filter(diff_diff_rss < 0) %>%
    head(n=1) %>%
    pull(estimate.measures.rank)
  disper = data.frame(estimate$measures$rank, estimate$measures$dispersion)%>%
    arrange(desc(estimate$measures$dispersion))%>%
    head(n=1) %>%
    pull(estimate.measures.rank)
  equal_value <- function(a, b, c) {
    if (a == b) {
      return(a)
    } else if (a == c) {
      return(a)
    } else if (b == c) {
      return(b)
    } else {
      return(NULL)
    }
    p = equal_value(coph,disper,rss)
    return(p)
  }
}

# ======================================================
# Function to Perform NMF and get gene metaprograms
# ======================================================
nmf_process = function(matrix,rank,method){
  result = nmf(matrix, rank = rank,method = method)
  fs = extractFeatures(result,30L)
  fs = lapply(fs,function(x)rownames(result)[x])
  fs
}

# ======================================================
# Function to perform GO Enrich Analysis
# ======================================================
GOenrich = function(allgene,diffgene){
  sep = ':.*'
  gl <- sub(sep, '', diffgene)
  back <- sub(sep, '', allgene)
  geneList <- factor(as.integer(back %in% gl))
  names(geneList) <- back
  back <- sub(sep, '', allgene)
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                geneSel = function(a) {a},
                annot = annFUN.org,
                mapping = "org.Hs.eg.db",
                ID = "Symbol")
  resultFisher <- topGO::runTest(GOdata, algorithm = "classic", statistic = "fisher")
  sigres <- topGO::GenTable(
    GOdata,
    classicFisher = resultFisher,
    topNodes = length(resultFisher@score),
    orderBy = "classicFisher",
    numChar = 1000
  )
  sigres$classicFisher[sigres$classicFisher == "< 1e-30"] <- 0
  sigres <- sigres[sigres$Annotated >= 10,]
  sigres$FDR <- stats::p.adjust(sigres$classicFisher, method = "fdr")
  ptcount <- 0
  fc <- ((sigres[, "Significant"] + ptcount) / (sum(GOdata@allScores[GOdata@feasible] == 1) + ptcount)) / ((sigres[, "Annotated"] + ptcount) / (sum(GOdata@feasible) +ptcount))
  sigres <- data.frame(sigres, FC = fc)
  sigres <- sigres[order(sigres$FDR,-sigres$FC),]
}

# ======================================================
# Function to plot GO Enrichment Heatmap
# ======================================================
plotGOEnrich <-
  function(goRes,
           n = 5,
           sortByFDR = TRUE,
           fdr.cutoff = 0.05,
           fc.cutoff = 2) {
    d <- lapply(names(goRes), function(i) {
      cbind(Cluster = i, goRes[[i]])
    })
    d <- do.call(rbind, d)
    if (sortByFDR) {
      d <- d[order(d$FDR, -d$FC), ]
    } else {
      d <- d[order(-d$FC, d$FDR), ]
    }
    cd <- do.call(rbind, sapply(sort(unique(d$Cluster)), function(i) {
      tmp <- d[d$Cluster == i, ]
      tmp <-
        tmp[tmp$FDR < fdr.cutoff & tmp$FC > fc.cutoff, , drop = FALSE]
      if (nrow(tmp) > 0) {
        tmp[seq_len(min(n, nrow(tmp))), ]
      } else {
        NULL
      }
    }, simplify = FALSE))
    ut <- unique(cd$Term)
    d <- d[d$Term %in% ut, c('Cluster', 'Term', 'FDR', 'FC')]
    
    d <- d[d$FDR < fdr.cutoff & d$FC > fc.cutoff, , drop = FALSE]
    
    dmat <- dcast(d, Term ~ Cluster)
    rownames(dmat) <- dmat[, 1]
    dmat <- as.matrix(dmat[, -1, drop = FALSE])
    dmat <- is.na(dmat)
    
    v <- sapply(seq_len(nrow(dmat)), function(i)
      which.min(dmat[i, ]))
    names(v) <- rownames(dmat)
    
    pd <- melt(dmat)
    colnames(pd) <- c('Term', 'Cluster', 'enc')
    pd$Cluster = as.character(pd$Cluster)
    pd$Term <-
      factor(as.character(pd$Term), levels = names(v[rev(order(v))]))
    pd$enc <- ifelse(pd$enc, 'Non-significant', 'Significant')
    p <-
      ggplot(pd, aes(x = pd[,2], y = pd[,1], fill = pd[,3])) + geom_tile() + theme_classic() + scale_fill_manual(values =
                                                                                                                   c('darkblue', 'orange')) +
      theme(legend.position = 'right',
            text = element_text(size = 12)) +
      scale_y_discrete(position = "right") +
      xlab('Cluster') + ylab('GO Terms')
    return(p)
    
  }
