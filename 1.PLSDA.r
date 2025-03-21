rm(list = ls())
####path#### 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#path <- getwd()
#path
####packages####
library(plyr, warn.conflicts = F)
library(ggplot2)
library(ggrepel)
library(ggforce)
#library(ape)
#library(vegan)
#library(GUniFrac)
#library(phyloseq)
library(mixOmics)
library(dplyr)
library(ropls)
####data####
asv.all <- read.delim('../1.ASV_table/tax.new/ASV_table_gut_CFe.tax.new.txt', row.names = 1)
asv.all$taxonomy <- NULL
group <- read.delim('../1.ASV_table/tax.new/ASV_table_group.txt')
group <- group[group$type == "Intestine" & 
                 group$Diets == "C" & 
                 group$Mineral == "Fe", ]
rownames(group) <- group$Sample
# Extract the sample names from the filtered group data frame
selected_samples <- group$Sample
# Select only the columns in the ASV table that exist in selected_samples
asv <- asv.all[, colnames(asv.all) %in% selected_samples, drop = FALSE]
asv <- asv[, match(selected_samples, colnames(asv))]
# Transpose the ASV table
asv <- data.frame(t(asv))
####PLS-DA####
plsda <- opls(asv, group$level, predI = 5, orthoI = 0)
# 提取 R2X 和 R2Y
R2X <- round(plsda@summaryDF$`R2X(cum)`[nrow(plsda@summaryDF)], 2)
R2Y <- round(plsda@summaryDF$`R2Y(cum)`[nrow(plsda@summaryDF)], 2)
Q2 <- round(plsda@summaryDF$`Q2(cum)`[nrow(plsda@summaryDF)], 2)
R2X
R2Y
Q2
pQ2 <- round(plsda@summaryDF$pQ2[nrow(plsda@summaryDF)], 2)
pR2Y <- round(plsda@summaryDF$pR2Y[nrow(plsda@summaryDF)], 2)
pQ2
pR2Y
plsda.breast <- plsda(asv, group$level, ncomp = 5, scale = TRUE)
P1 <- paste0("P1: (", sprintf("%.2f", round(plsda.breast$prop_expl_var$X[1], 4) * 100), "%)")
P1
P2 <- paste0("P2: (", sprintf("%.2f", round(plsda.breast$prop_expl_var$X[2], 4) * 100), "%)")
P2
P3 <- paste0("P3: (", sprintf("%.2f", round(plsda.breast$prop_expl_var$X[3], 4) * 100), "%)")
P3
P4 <- paste0("P4: (", sprintf("%.2f", round(plsda.breast$prop_expl_var$X[4], 4) * 100), "%)")
P4
####plot####
site <- as.data.frame(plsda.breast$variates)[1:4]
sample_site <- site
colnames(sample_site) <- c("PCoA1", "PCoA2", "PCoA3", "PCoA4")
sample_site$level <- group$level
sample_site[, 1:4] <- lapply(sample_site[, 1:4], as.numeric)
####group.order####
#sample_site$type <- factor(group$type, levels = unique(group$type))
sample_site$level <- factor(sample_site$level, levels = unique(group$level))
#sample_site$copper <- factor(group$copper, levels = unique(group$copper))
####midpoint.PC1 PC2####
find_hull.1 <- function(sample_site) sample_site[chull(sample_site$PCoA1, sample_site$PCoA2),]
hulls.1 <- ddply(sample_site, "level", find_hull.1)
PCoA1_mean <- tapply(sample_site$PCoA1,sample_site$level, mean)
PCoA2_mean <- tapply(sample_site$PCoA2,sample_site$level, mean)
mean_point12 <- rbind(PCoA1_mean,PCoA2_mean)
mean_point12 <- data.frame(mean_point12)
t1 <- t(mean_point12)
t2 <- as.data.frame(t1)
t2$time <- table(sample_site$level)
####midpoint.PC1 PC3####
find_hull.2 <- function(sample_site) sample_site[chull(sample_site$PCoA1, sample_site$PCoA2),]
hulls.2 <- ddply(sample_site, "level", find_hull.2)
PCoA3_mean <- tapply(sample_site$PCoA3,sample_site$level, mean)
mean_point13 <- rbind(PCoA1_mean,PCoA3_mean)
mean_point13 <- data.frame(mean_point13)
t3 <- t(mean_point13)
t4 <- as.data.frame(t3)
t4$time <- table(sample_site$level)
####midpoint.PC2 PC3####
find_hull.3 <- function(sample_site) sample_site[chull(sample_site$PCoA1, sample_site$PCoA2),]
hulls.3 <- ddply(sample_site, "level", find_hull.3)
mean_point23 <- rbind(PCoA2_mean,PCoA3_mean)
mean_point23 <- data.frame(mean_point23)
t5 <- t(mean_point23)
t6 <- as.data.frame(t5)
t6$time <- table(sample_site$level)
####midpoint.PC1 PC4####
find_hull.4 <- function(sample_site) sample_site[chull(sample_site$PCoA1, sample_site$PCoA4),]
hulls.4 <- ddply(sample_site, "level", find_hull.4)
PCoA4_mean <- tapply(sample_site$PCoA4,sample_site$level, mean)
mean_point14 <- rbind(PCoA1_mean,PCoA4_mean)
mean_point14 <- data.frame(mean_point14)
t7 <- t(mean_point14)
t8 <- as.data.frame(t7)
t8$time <- table(sample_site$level)
color3 <- c('#9B2227', '#808000', '#006073','#fa4572', '#874efd',
            '#EE9B00', '#CC6602', '#3CB371','#1e90ff', '#90fc27')
color4 <- c('#9B2227', '#808000', '#006073','#fa4572', '#874efd',
            '#EE9B00', '#CC6602', '#3CB371','#1e90ff', '#90fc27',
            '#FFC0CB', '#0ddFFF', '#FF6347', '#FFD700', '#DA70D6',
            '#C9CBFF', '#6495ED', '#FFA07A', '#20B2AA', '#FFE4B5')
#c('#FFA500', '#32CD32', '#BA55D3', '#40E0D0', '#FF69B4', '#FFDAB9', '#00BFFF', '#FF4500', '#7B68EE', '#FFE4B5')
#c('#FF8C00', '#3CB371', '#8A2BE2', '#FF00FF', '#00FFFF', '#FF6384', '#4BC0C0', '#C9CBFF', '#FFCD56', '#A0D468')
color2 <- c('#9B2227', '#808000', '#006073','#fa4572', '#874efd',
            '#EE9B00', '#CC6602', '#1e90ff')
color <- c( '#FFC0CB', '#0ddFFF', '#FF6347', '#FFD700', '#DA70D6',
            '#C9CBFF', '#6495ED', '#FFA07A', '#20B2AA', '#FFE4B5')
colors <- rev(color)
color2s <- rev(color2)
color3s <- rev(color3)
color4s <- rev(color4)
####PC1 PC2####
sample_site$sample_name <- rownames(sample_site)
p1 <- ggplot(sample_site, aes(PCoA1, PCoA2, color = level)) +
  geom_vline(xintercept = 0, color = 'black', linewidth = 0.5, linetype = 2) + 
  geom_hline(yintercept = 0, color = 'black', linewidth = 0.5, linetype = 2) +  
  geom_point(aes(color = level, shape = level, size = level)) +
  geom_segment(aes(x = PCoA1, xend = rep(t2[, 1], times = t2[, 3]),
                   y = PCoA2, yend = rep(t2[, 2], times = t2[, 3])),
               linetype = 2, linewidth = 0.5, show.legend = F) +
  scale_shape_manual(values = c(#5,8,11,
    17,15,16,18,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16)) + 
  scale_color_manual(values = color3) +
  scale_size_manual(values = c(2,2,2,2.5,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3)) +
  stat_ellipse(aes(color = level), linetype = 1, linewidth = 0.50,
               level = 0.45, show.legend = F, type = "norm") +   
  labs(x = P1, y = P2) +
  labs(title = "PLS-DA") +
  ggprism::theme_prism()  +
# geom_text_repel(aes(label = sample_name), size = 2, box.padding = 0.3, point.padding = 0.3,
#                 segment.color = "grey", segment.size = 0.5, max.overlaps = 20) + 
# geom_text(data = data.frame(x = 10, y = 10, label = "R2 = 0.216, P = 0.053"),
#           mapping = aes(x = x, y = y, label = label), size = 2.8, fontface = 2, inherit.aes = FALSE) +
  theme(legend.text = element_text(face = "bold", size = 10), 
        plot.title = element_text(size = 12),
        axis.line = element_blank(), 
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 11),
        axis.ticks = element_line(linewidth = 0.5),
        axis.ticks.length = unit(0.1, "cm"),
        panel.background = element_rect(colour = "black", linetype = "solid", linewidth = 0.50))
p1
ggsave("1.PLS-DA.12.pdf", p1, width = 110, height = 85, units = 'mm')
####PC1 PC3####
p2 <- ggplot(sample_site, aes(PCoA1, PCoA3, color = level)) +
  geom_vline(xintercept = 0, color = 'black', linewidth = 0.5, linetype = 2) + 
  geom_hline(yintercept = 0, color = 'black', linewidth = 0.5, linetype = 2) +  
  geom_point(aes(color = level, shape = level, size = level)) +
  geom_segment(aes(x = PCoA1, xend = rep(t4[, 1], times = t4[, 3]),
                   y = PCoA3, yend = rep(t4[, 2], times = t4[, 3])),
               linetype = 2, linewidth = 0.5, show.legend = F) +
  scale_shape_manual(values = c(#5,8,11,
    17,15,16,18,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16)) + 
  scale_color_manual(values = color3) +
  scale_size_manual(values = c(2,2,2,2.5,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3)) +
  stat_ellipse(aes(color = level), linetype = 1, linewidth = 0.50,
               level = 0.4, show.legend = F, type = "norm") +   
  labs(x = P1, y = P3) +
  labs(title = "PLS-DA") +
  ggprism::theme_prism() +
  # geom_text_repel(aes(label = sample_name), size = 2, box.padding = 0.3, point.padding = 0.3,
  #                 segment.color = "grey", segment.size = 0.5, max.overlaps = 20) + 
  # geom_text(data = data.frame(x = 0, y = 0.25, label = "R2 = 0.501, P < 0.001"),
  #           mapping = aes(x = x, y = y, label = label), size = 2.8, fontface = 2, inherit.aes = FALSE) +
  theme(legend.text = element_text(face = "bold", size = 10), 
        plot.title = element_text(size = 12),
        axis.line = element_blank(), 
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 11),
        axis.ticks = element_line(linewidth = 0.5),
        axis.ticks.length = unit(0.1, "cm"),
        panel.background = element_rect(colour = "black", linetype = "solid", linewidth = 0.50))
p2
ggsave("1.PLS-DA.13.pdf", p2, width = 110, height = 85, units = 'mm')
####PC2 PC3####
p3 <- ggplot(sample_site, aes(PCoA2, PCoA3, color = level)) +
  geom_vline(xintercept = 0, color = 'black', linewidth = 0.5, linetype = 2) + 
  geom_hline(yintercept = 0, color = 'black', linewidth = 0.5, linetype = 2) +  
  geom_point(aes(color = level, shape = level, size = level)) +
  geom_segment(aes(x = PCoA2, xend = rep(t6[, 1], times = t6[, 3]),
                   y = PCoA3, yend = rep(t6[, 2], times = t6[, 3])),
               linetype = 2, linewidth = 0.5, show.legend = F) +
  scale_shape_manual(values = c(#5,8,11,
    17,15,16,18,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16)) + 
  scale_color_manual(values = color3) +
  scale_size_manual(values = c(2,2,2,2.5,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3)) +
  stat_ellipse(aes(color = level), linetype = 1, linewidth = 0.50,
               level = 0.2, show.legend = F, type = "norm") +   
  labs(x = "PC2 (12.19%)", y = "PC3 (8.03%)") +
  labs(title = "PLS-DA") +
  ggprism::theme_prism() +
  # geom_text_repel(aes(label = sample_name), size = 2, box.padding = 0.3, point.padding = 0.3,
  #                 segment.color = "grey", segment.size = 0.5, max.overlaps = 20) + 
  # geom_text(data = data.frame(x = 0, y = 0.25, label = "R2 = 0.501, P < 0.001"),
  #           mapping = aes(x = x, y = y, label = label), size = 2.8, fontface = 2, inherit.aes = FALSE) +
  theme(legend.text = element_text(face = "bold", size = 10), 
        plot.title = element_text(size = 12),
        axis.line = element_blank(), 
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 11),
        axis.ticks = element_line(linewidth = 0.5),
        axis.ticks.length = unit(0.1, "cm"),
        panel.background = element_rect(colour = "black", linetype = "solid", linewidth = 0.50))
p3
#ggsave("1.PC2_PC3_weighted.pdf", p3, width = 110, height = 85, units = 'mm')
####PC1 PC4####
p4 <- ggplot(sample_site, aes(PCoA1, PCoA4, color = level)) +
  geom_vline(xintercept = 0, color = 'black', linewidth = 0.5, linetype = 2) + 
  geom_hline(yintercept = 0, color = 'black', linewidth = 0.5, linetype = 2) +  
  geom_point(aes(color = level, shape = level, size = level)) +
  geom_segment(aes(x = PCoA1, xend = rep(t8[, 1], times = t8[, 3]),
                   y = PCoA4, yend = rep(t8[, 2], times = t8[, 3])),
               linetype = 2, linewidth = 0.5, show.legend = F) +
  scale_shape_manual(values = c(#5,8,11,
    17,15,16,18,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16)) + 
  scale_color_manual(values = color3) +
  scale_size_manual(values = c(2,2,2,2.5,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3)) +
  stat_ellipse(aes(color = level), linetype = 1, linewidth = 0.50,
               level = 0.2, show.legend = F, type = "norm") +    
  labs(x = "PC1 (48.13%)", y = "PC4 (7.25%)") +
  labs(title = "PCoA") +
  ggprism::theme_prism() +
  # geom_text_repel(aes(label = sample_name), size = 2, box.padding = 0.3, point.padding = 0.3,
  #                 segment.color = "grey", segment.size = 0.5, max.overlaps = 20) + 
  # geom_text(data = data.frame(x = 0, y = 0.25, label = "R2 = 0.501, P < 0.001"),
  #           mapping = aes(x = x, y = y, label = label), size = 2.8, fontface = 2, inherit.aes = FALSE) +
  theme(legend.text = element_text(face = "bold", size = 10), 
        plot.title = element_text(size = 12),
        axis.line = element_blank(), 
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 11),
        axis.ticks = element_line(linewidth = 0.5),
        axis.ticks.length = unit(0.1, "cm"),
        panel.background = element_rect(colour = "black", linetype = "solid", linewidth = 0.50))
p4
#ggsave("1.PC1_PC4_weighted.pdf", p4, width = 110, height = 85, units = 'mm')
sessionInfo()