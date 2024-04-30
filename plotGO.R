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
plotGOEnrich(goRes = res, n = 5, fdr.cutoff = .1, fc.cutoff = 2)


