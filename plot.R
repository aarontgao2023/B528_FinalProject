# ======================================================
# Cell-Cell Correlation
# ======================================================

matrix = read.table("GSM3905417_MUV41.txt")
matrix = data_process(matrix)
cm = cor(matrix)
set.seed(12345)
hclu <- cutree(hclust(as.dist(1-cm),method = "ward.D"), k = 3)
cellorder = names(sort(hclu))
cellorder <- unlist(cellorder)
max_2 = max(cm[cm != max(cm)])
cm[cm==1]=max_2
cm = rescale(cm, to= c(-0.1,0.2))
pheatmap(cm[cellorder,cellorder], 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_rownames=F,
         show_colnames=F,
         fontsize = 15,
         height = 3,
         width = 4.5,
         filename = "muv41_corr.png",
         color = colorRampPalette(c("dark blue","blue","green","yellow","red","dark red"))(200))

# ======================================================
# PCA Scatter Plot
# ======================================================

library(splines)
pca_lmFilter <- function(genebycellmat, topnum = 1e3){
  d = genebycellmat
  m = rowMeans(d)
  v = apply(d,1,sd)
  resid = resid(lm(v~bs(m)))
  d = d[names(sort(resid,decreasing = T))[1:topnum], ]
  prcomp(t(d),scale. = T)
}
pcres = pca_lmFilter(matrix, topnum = 1000)
pc = as.data.frame(pcres$x[,1:30])
plotdata <- data.frame(PC1 = pc$PC1, PC2 = pc$PC2, cluster = hclu)
plotdata$cluster = as.factor(plotdata$cluster)
ggplot() + geom_point(data = plotdata, aes(x = PC1, y = PC2, color = cluster),size=0.2)+  theme_classic() +  
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 7),
        axis.ticks = element_line(colour="black", linetype = 1, size= 0.25),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6,  vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        plot.title = element_text(face = "plain", size = 8, hjust = 0.5),
        plot.margin = unit(c(3, 3, 4, 4), "pt"),
        legend.background = element_blank())+xlab("Principal component 1") +ylab("Principal component 2")
ggsave("pca_scatter_muv41.png",height=2, width = 3)


# ======================================================
# PCA Heatmap
# ======================================================
plotpc = apply(pc,2,function(x){rescale(x, to = c(-30,30))})
pheatmap(pc[cellorder,1:5],
         show_rownames=F,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         fontsize=15,
         angle_col = 45,
         height = 3,
         width = 4.5,
         filename = "pcaheatmap_muv41.png")

# ======================================================
# NMF Heatmap
# ======================================================
matrix[matrix<0]=0.000001
res = nmf(matrix, rank = 3)
coefmap(res,tracks = NA,filename= "nmfheatmap.png",hclustfun = "ward")
basismap(res,filename= "nmfheatmap1.png",hclustfun = "ward", tracks = NA)
fs = extractFeatures(res,30L)
fs = lapply(fs,function(x)rownames(res)[x])

k1 <- fs[[1]]
mat1 <- matrix[rownames(matrix) %in% k1, ]
mean1 = as.data.frame(colMeans(mat1))
k2 <- fs[[2]]
mat2 <- matrix[rownames(matrix) %in% k2, ] 
mean2 = as.data.frame(colMeans(mat2))
k3 <- fs[[3]]
mat3 <- matrix[rownames(matrix) %in% k3, ] 
mean3 = as.data.frame(colMeans(mat3))
nmfplot = data.frame(mean1,mean2,mean3)
nmfplot = apply(nmfplot,2,function(x){rescale(x, to = c(-4,4))})
colnames(nmfplot) = c("NMF Meta1","NMF Meta2","NMF Meta3")
pheatmap(nmfplot[cellorder,],
         show_rownames=F,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         fontsize=15,
         angle_col = 0,
         height = 3,
         width = 4.5,
         color = colorRampPalette(c("blue","white","red"))(500),
         filename = "nmfheatmap_muv41.png")


# ======================================================
# barplot
# ======================================================
barplot_mouse = readxl::read_xlsx("/Users/tianchuangao/Desktop/Book1.xlsx", sheet = "Sheet3")

ggplot(data = barplot_mouse,
       mapping = aes(x = reorder(Organ,Number), y = Number,fill=Organ)) + 
  geom_col() +
  guides(fill = "none")+ 
  labs(x="Organ",y = "Number of samples") +
  coord_flip()+
  theme_classic()+
  theme(text = element_text(size = 15))
ggsave("barplot_mouse.png",height = 3,width = 4.5)


barplot_human = readxl::read_xlsx("/Users/tianchuangao/Desktop/Book1.xlsx", sheet = "Sheet1")
ggplot(data = barplot_human,
       mapping = aes(x = reorder(Organ,Number), y = Number,fill=Organ)) + 
  geom_col() +
  guides(fill = "none")+ 
  labs(x="Organ",y = "Number of samples") +
  coord_flip()+
  theme_classic()+
  theme(text = element_text(size = 15))
ggsave("barplot_human.png",height = 5,width = 5)

