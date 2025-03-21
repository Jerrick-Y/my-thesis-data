#rm(list = ls())
####path#### 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#path <- getwd()
#path
####packages####
library(dplyr)
library(ape)
library(phyloseq)
#detach("package:phyloseq", unload = TRUE)
####data####
asv.all <- read.delim('../1.ASV_table/tax.new/ASV_table_gut_CFe.tax.new.txt', row.names = 1)
asv.all$taxonomy <- NULL
group <- read.delim('../1.ASV_table/tax.new/ASV_table_group.txt')
group <- group[group$type == "Intestine" & 
                 group$Diets == "C" & 
                 group$Mineral == "Fe", ]
# Extract the sample names from the filtered group data frame
selected_samples <- group$Sample
# Select only the columns in the ASV table that exist in selected_samples
asv <- asv.all[, colnames(asv.all) %in% selected_samples, drop = FALSE]
asv <- asv[, match(selected_samples, colnames(asv))]
asv <- data.frame(t(asv))
asv$Sample <- rownames(asv)
group$type <- as.factor(group$level)
merged_data <- inner_join(group, asv, by = "Sample")
#ASV_sum <- merged_data %>%
#  mutate(row_total = rowSums(across(starts_with("ASV")))) %>%
#  mutate(across(starts_with("ASV"), ~ .x / row_total * 1000000)) %>%
#  mutate(scaled_row_total = rowSums(across(starts_with("ASV"))))
ASV.sum <- merged_data %>%
  group_by(factor(type, level = unique(merged_data$type))) %>%
  summarise(across(starts_with("ASV"), sum)) %>%
  rename(type = 'factor(type, level = unique(merged_data$type))')
ASV.sum <- data.frame(ASV.sum)
rownames(ASV.sum) <- ASV.sum$type
ASV.sum <- ASV.sum[, -1]
tree <- read.tree("../1.ASV_table/fasta/ASV_table_gut_CFe_tree.rooted.nwk")
colnames(asv) <- tree$tip.label
df <- phyloseq(otu_table(ASV.sum, taxa_are_rows = F), phy_tree(tree)) #phyloseq
dist.mat <- phyloseq::UniFrac(df, weighted = T) #weighted
uniFrac.dis <- as.dist(dist.mat)
upgma <- hclust(uniFrac.dis, method = "average")
upgma_phylo <- as.phylo(upgma)
write.tree(upgma_phylo, file = "weighted.upgma.group.nwk")
weighted.tree <- read.tree("weighted.upgma.group.nwk")
detach("package:phyloseq", unload = TRUE)
library(ggplot2)
library(ggtree)
library(treeio)
group$type <- factor(group$type, level = unique(group$type))
groupInfo <- split(group$type, group$type)
weighted.tree <- groupOTU(weighted.tree, groupInfo, group_name = 'group1')
color3 <- c('#9B2227', '#808000', '#006073','#fa4572', '#874efd',
            '#EE9B00', '#CC6602', '#3CB371','#1e90ff', '#90fc27')
color4 <- c('#9B2227', '#808000', '#006073','#fa4572', '#874efd',
            '#EE9B00', '#CC6602', '#3CB371','#1e90ff', '#90fc27',
            '#FFC0CB', '#0ddFFF', '#FF6347', '#FFD700', '#DA70D6',
            '#C9CBFF', '#6495ED', '#FFA07A', '#20B2AA', '#FFE4B5','#000000')
color5 <- c('#9B2227', '#808000', '#006073','#fa4572', '#874efd',
            '#EE9B00', '#CC6602', '#3CB371','#1e90ff', '#90fc27',
            '#FFC0CB', '#0ddFFF',
            '#9B2227', '#808000', '#006073','#fa4572', '#874efd',
            '#EE9B00', '#CC6602', '#3CB371','#1e90ff', '#90fc27',
            '#FFC0CB', '#0ddFFF', '#000000', '#000000')
#c('#FFA500', '#32CD32', '#BA55D3', '#40E0D0', '#FF69B4', '#FFDAB9', '#00BFFF', '#FF4500', '#7B68EE', '#FFE4B5')
#c('#FF8C00', '#3CB371', '#8A2BE2', '#FF00FF', '#00FFFF', '#FF6384', '#4BC0C0', '#C9CBFF', '#FFCD56', '#A0D468')
color2 <- c('#9B2227', '#808000', '#006073','#fa4572', '#874efd',
            '#EE9B00', '#CC6602', '#1e90ff')
color <- c('#9B2227', '#808000')
colors <- rev(color)
color2s <- rev(color2)
color3s <- rev(color3)
color4s <- rev(color4)
weighted.tree$edge.length <- round(weighted.tree$edge.length, 4)
p.weighted <- ggtree(weighted.tree, layout = "rectangular", aes(color = group1)) +
  ggtitle("UPGMA") + 
  geom_tiplab(size = 2.8, fontface = 2, hjust = -.1) + 
  scale_color_manual(values = color3, breaks = unique(group$type)) +
  geom_text(aes(label= branch.length, x =  branch), vjust = 1.7, 
            size = 2.5, fontface = 2) + 
  geom_nodelab(aes(subset =!isTip, label = node),
               hjust = -0.5, color =  '#000000', size = 2) +
  #geom_tiplab(aes(label = node), size = 2, color = '#000000', offset = 0, hjust = 0) + 
  #geom_hilight(node = 235, fill = "#808000",alpha = 0.2) +
  ggprism::theme_prism()  +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "none",
        legend.text = element_text(face = "bold", size = 10), 
        legend.background = element_blank(),
        axis.line = element_blank(), 
        panel.background = element_blank(),
        axis.ticks = element_line(linetype = "blank"),
        axis.text = element_blank()) +
  labs(colour = NULL) 
p.weighted
p <- flip(p.weighted, 3, 4)
#p <- flip(p, 6, 34)
#p <- flip(p, 6, 34)
#p
p.1 <- p +
  #geom_cladelabel(node = 229, label = "Intestine", barsize = 1, color = "#3CB371", offset = 2.1, offset.text = 0.7) +
  #geom_cladelabel(node = 153,label = "Water", barsize = 1, color = "#FF6347", offset = 2.1, offset.text = 0.7) +
  #geom_cladelabel(node = 155,label = "WS", barsize = 1, color = "#874efd", offset = 2.1, offset.text = 0.7) +
  geom_hilight(node = 41, fill = "#3CB371", alpha = 0.2) +
  geom_hilight(node = 47, fill = "#C9CBFF", alpha = 0.2) +
  #geom_hilight(node = 41, fill = "#808000", alpha = 0.2) +
  geom_hilight(node = 28, fill = "#FF6347", alpha = 0.2) +
  geom_hilight(node = 29, fill = "#0ddFFF", alpha = 0.2) +
  geom_hilight(node = 1, fill = "#874efd", alpha = 0.2) 
# geom_strip(87, 70, barsize = 1,color = "#3CB371",alpha = 0.3, hjust = 0.5, label = "Intestine",  angle = 0, offset = 3, offset.text = 1,fontsize = 2.8) +
# geom_strip(87, 75, barsize = 1,color = "#FF6347",alpha = 0.3, hjust = 0.5, vjust = 0.5, label = "Water",  angle = -50, offset = 3, offset.text = 1,fontsize = 2.8) +
# geom_strip(65, 70, barsize = 1,color = "#874efd",alpha = 0.3, hjust = 0.5, vjust = 0.5, label = "WS", angle = 40, offset = 3, offset.text = 1,fontsize = 2.8) 
p.1
ggsave("3.weighted.upgma.group.pdf", p,  width = 60, height = 60, units = 'mm')
####unweighted####
un.dist.mat <- phyloseq::UniFrac(df, weighted = F)
un.uniFrac.dis <- as.dist(un.dist.mat)
un.upgma <- hclust(un.uniFrac.dis, method = "average")
un.upgma_phylo <- as.phylo(un.upgma)
write.tree(un.upgma_phylo, file = "unweighted.upgma.group.nwk")
unweighted.tree <- read.tree("unweighted.upgma.group.nwk")
unweighted.tree <- groupOTU(unweighted.tree, groupInfo, group_name = 'group1')
#rectangular
unweighted.tree$edge.length <- round(unweighted.tree$edge.length, 4)
p.unweighted <- ggtree(unweighted.tree, layout = "rectangular", aes(color = group1)) +
  ggtitle("UPGMA") + 
  geom_tiplab(size = 2.8, fontface = 2, hjust = -.1) + 
  scale_color_manual(values = color5, breaks = unique(group$type)) +
  geom_text(aes(label= branch.length, x =  branch), vjust = 1.7, 
            size = 2.5, fontface = 2) + 
  geom_nodelab(aes(subset =!isTip, label = node),
               hjust = -0.5, color =  '#000000', size = 2) +
  #geom_tiplab(aes(label = node), size = 2, color = '#000000', offset = 0, hjust = 0) + 
  #geom_hilight(node = 235, fill = "#808000",alpha = 0.2) +
  ggprism::theme_prism()  +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "none",
        legend.text = element_text(face = "bold", size = 10), 
        legend.background = element_blank(),
        axis.line = element_blank(), 
        panel.background = element_blank(),
        axis.ticks = element_line(linetype = "blank"),
        axis.text = element_blank()) +
  labs(colour = NULL) 
p.unweighted 
p.p <- flip(p.unweighted, 3, 4) 
p.p <- flip(p.p, 8, 35) 
p.p <- flip(p.p, 7, 34) 
p.p <- flip(p.p, 17, 43) 
p.p <- flip(p.p, 44, 47) 
p.p
p.p.1 <- p.p +
  geom_hilight(node = 36, fill = "#3CB371", alpha = 0.2) +
  geom_hilight(node = 32, fill = "#C9CBFF", alpha = 0.2) +
  #geom_hilight(node = 41, fill = "#808000", alpha = 0.2) +
  geom_hilight(node = 44, fill = "#FF6347", alpha = 0.2) +
  geom_hilight(node = 47, fill = "#0ddFFF", alpha = 0.2) +
  geom_hilight(node = 2, fill = "#874efd", alpha = 0.2) 
#  geom_strip(120, 140, barsize = 1,color = "#3CB371",alpha = 0.3, hjust = 0.5, label = "Intestine",  angle = 0, offset = 3, offset.text = 1,fontsize = 2.8) +
#  geom_strip(8, 1, barsize = 1,color = "#FF6347",alpha = 0.3, hjust = 0.5, label = "Water",  angle = 0, offset = 3, offset.text = 1,fontsize = 2.8) +
#  geom_strip(2, 7, barsize = 1,color = "#874efd",alpha = 0.3, hjust = 0.5, vjust = 0.5, label = "WS", angle = 90, offset = 3, offset.text = 1,fontsize = 2.8) 
p.p.1
ggsave("3.unweighted.upgma.group.pdf", p.p, width = 60, height = 60, units = 'mm')
####clean###
#rm(list = ls())
#"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), 
#"mcquitty" (= WPGMA), 
#"median" (= WPGMC) 
#"centroid" (= UPGMC).
sessionInfo()
