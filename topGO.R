library(topGO)
data(geneList)
sampleGOdata <- new("topGOdata",
                      description = "Simple session", ontology = "BP",
                      allGenes = geneList, geneSel = topDiffGenes,
                      nodeSize = 10,
                      annot = annFUN.db, affyLib = affyLib)

genels = as.vector(fs_my[[1]])
geneList
head(geneList)
