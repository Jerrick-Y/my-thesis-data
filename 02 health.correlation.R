rm(list = ls())
####获得当前路径#### 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
####correlation####
library(linkET)
library(ggplot2)
library(RColorBrewer)
library(corrplot)
library(scales)
muscle <- read.delim2('02 corr.phy.all.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
muscle <- as.data.frame(lapply(muscle, as.numeric)) #[,-c(7,12,15,17,19)]
#mtcars
col <- colorRampPalette(c('blue','white','red'))(200)
#plot(1:10,rep(1,10),col = col,pch = 16,cex = 2)
p <- correlate(muscle, type = "lower", method = "pearson") %>% 
  qcorrplot() +
  geom_square(data = filter_func(type = "lower")) +
  geom_mark(data = filter_func(type = "upper"), sep = '\n', size = 3) +
  geom_diag_label(size = 4, colour = "black", angle = 0) +
  scale_fill_gradientn(colours = col, 
                       values = rescale(c(-0.35, 0, 1)),
                       name = NULL,
                       limits = c(-0.35, 1),
                       breaks = c(-0.3, 0, 0.5, 1),
                       labels = c(-0.3, 0, 0.5, 1)) +
  theme_void()
p
ggsave("02.correlation5.pdf", plot = p, width = 8.8, height = 8.8, limitsize = FALSE)

library(PerformanceAnalytics)
qPCR <- read.delim2('02.health.all.txt')
qPCR <- as.data.frame(lapply(qPCR, as.numeric))
########
pdf('02.correlation.all.spearman.2.pdf', width = 13, height = 13)
chart.Correlation(qPCR, histogram = TRUE, pch = 15, method = c( "spearman"))
dev.off()
cor(qPCR,method = c( "spearman"))

library(tidyr)
library(dplyr)
library(factoextra)
PCA.qPCR <- read.delim2('02.health.all.txt')
PCA.qPCR <- as.data.frame(lapply(PCA.qPCR, as.numeric))
#PCA.qPCR <- apply(PCA.qPCR, 2, log)
pca_qPCR <- prcomp(PCA.qPCR, scale = TRUE)
fviz_eig(pca_qPCR)
fviz_contrib(pca_qPCR, choice = "var", axes = 1)
qPCR.PC1 <- pca_qPCR$x %>% 
  as.data.frame() %>% 
  select(PC1) %>%
  rename(qPCR = PC1) %>%
  #rownames_to_column("FishID") %>%
  uncount(1)
write.table(qPCR.PC1,file = "03.PC1.2.txt",sep = "\t", quote = FALSE, row.names = F)

PCA.phy <- read.delim2('PCA.physio.txt')
PCA.phy <- as.data.frame(lapply(PCA.phy[,-(1:3)], as.numeric))
#PCA.qPCR <- apply(PCA.qPCR, 2, log)
pca_phy <- prcomp(PCA.phy, scale = TRUE)
fviz_eig(pca_phy)
fviz_contrib(pca_phy, choice = "var", axes = 1)
phy.PC1 <- pca_phy$x %>% 
  as.data.frame() %>% 
  select(PC1,PC2) %>%
  rename(phy = PC1) %>%
  rename(phy2 = PC2) %>%
  #rownames_to_column("FishID") %>%
  uncount(1)
write.table(phy.PC1,file = "phy.PC.txt",sep = "\t", quote = FALSE, row.names = F)

