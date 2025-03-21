rm(list = ls())
####path#### 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
####package####
library(tidyverse)
library(vegan)
library(ggprism)
library(ggplot2)
library(reshape2)
library(ggalluvial)
library(dplyr)
library(conflicted)
#load(".RData")
####table####
abasv <- read.delim('../1.ASV_table/tax.new/ASV_table_gut_CFe.tax.new.txt', row.names = 1)
group <- read.delim('../1.ASV_table/tax.new/ASV_table_group.txt')
group <- group[group$type == "Intestine" & 
                 group$Diets == "C" & 
                 group$Mineral == "Fe", ]
# Extract the sample names from the filtered group data frame
selected_samples <- group$Sample
# Select only the columns in the ASV table that exist in selected_samples
asv <- abasv[, colnames(abasv) %in% selected_samples, drop = FALSE]
asv <- asv[, match(selected_samples, colnames(asv))]
tax <- read.delim('../1.ASV_table/tax/ASV_table_gut_CFe.tax.txt', row.names = 1)
reasv <- cbind(asv[, c(1:NCOL(asv))])
reasv <- t(reasv)
reasv <- reasv/rowSums(reasv)
reasv <- data.frame(t(reasv))
dat <- merge(x = reasv[, c(1:(NCOL(reasv)))], y = tax, by = "row.names")
dat <- dplyr::rename(dat, OTUID = Row.names)
#group <- read.delim('../ASV_table_group_PC.txt')
group1 <- group
row.names(group1) <- group1$Sample
group1 <- group1[, c("level"), drop = FALSE]
####color####
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
color <- c('#9B2227', '#808000')
colors <- rev(color)
color2s <- rev(color2)
color3s <- rev(color3)
color4s <- rev(color4)
####phylum####
#class(dat$C1_1)
phylum <- aggregate(dat[, c(2:(NCOL(dat) - 7))], by = list(dat$Phylum), FUN = sum)
row.names(phylum) <- phylum$Group.1
phylum <- dplyr::select(phylum,-Group.1)
order <- sort(rowSums(phylum[,1:ncol(phylum)]),
              index.return = TRUE, decreasing = T)
c.phylum <- phylum[order$ix,]
####phylum.Top10####
d.phylum <- rbind(colSums(c.phylum[11:as.numeric(length(rownames(c.phylum))),]),c.phylum[10:1,])
rownames(d.phylum)[1] <- "Others"
d.phylum <- as.data.frame(t(d.phylum))
d.phylum$Unassigned <- NULL
d.phylum$Others <- NULL
merged_table <- merge(d.phylum, group1, by = "row.names", sort = FALSE)
rownames(merged_table) <- merged_table[, 1]
merged_table <- merged_table[, -1]
d.phylum.group <- merged_table %>%
  group_by(factor(level, levels = unique(merged_table$level))) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  rename(level = 'factor(level, levels = unique(merged_table$level))')
d.phylum.group <- as.data.frame(d.phylum.group)
rownames(d.phylum.group) <- d.phylum.group[, 1]
d.phylum.group <- d.phylum.group[, -1]
write.table(t(d.phylum.group), "1.phylum.group.txt", sep = "\t", col.names = NA, quote = F)
write.table(t(d.phylum), "1.phylum.sample.txt", sep = "\t", col.names = NA, quote = F)
####phylum.melt.long####
data.type.phylum.group <- melt(t(t(d.phylum.group)),id.vars = 'phylum')
type.phylum.group <- ggplot(data.type.phylum.group,aes(x = Var1, y = value, alluvium = Var2, stratum = Var2)) +
  geom_alluvium(aes(color = Var2),alpha = 1,width = 0.6, linetype = 1, lwd = 0.5, fill = 'transparent') + 
  geom_bar(aes(x = Var1, y = value, fill = Var2), position = "stack", stat = "identity", width = 0.6) +
  scale_color_manual(values = color3s) +
  scale_fill_manual(values = color3s) +
  labs(x = "", y = "Relative abundance") +
  labs(title = "Phylum") +
  ggprism::theme_prism() +
  theme(legend.text = element_text(face = "bold", size = 9),
        legend.key.size = unit(0.7, "lines"),
        legend.key.width = unit(0.45, "lines"),
        plot.title = element_text(size = 11),
        axis.ticks.length = unit(0.1, "cm"),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 9,
                                   color = "black", face = "bold"),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.5)) +
  scale_y_continuous(expand = c(0,0)) 
type.phylum.group
ggsave("1.phylum.group.pdf", type.phylum.group, width = 100, height = 75, units = 'mm')
#genus####
genus <- aggregate(dat[, c(2:(NCOL(dat) - 7))], by = list(dat$Genus), FUN = sum)
row.names(genus) <- genus$Group.1
genus <- dplyr::select(genus,-Group.1)
order <- sort(rowSums(genus[,1:ncol(genus)]),
              index.return = TRUE, decreasing = T)
c.genus <- genus[order$ix,]
write.table(c.genus, "genus.all.sample.txt", sep = "\t", col.names = T, quote = F)
####genus.Top10####
d.genus <- rbind(colSums(c.genus[11:as.numeric(length(rownames(c.genus))),]),c.genus[10:1,])
rownames(d.genus)[1] <- "Others"
d.genus <- as.data.frame(t(d.genus))
d.genus$Unassigned <- NULL
d.genus$Others <- NULL
merged_table <- merge(d.genus, group1, by = "row.names", sort = FALSE)
rownames(merged_table) <- merged_table[, 1]
merged_table <- merged_table[, -1]
d.genus.group <- merged_table %>%
  group_by(factor(level, levels = unique(merged_table$level))) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  rename(level = 'factor(level, levels = unique(merged_table$level))')
d.genus.group <- as.data.frame(d.genus.group)
rownames(d.genus.group) <- d.genus.group[, 1]
d.genus.group <- d.genus.group[, -1]
write.table(t(d.genus.group), "1.genus.group.txt", sep = "\t", col.names = NA, quote = F)
write.table(t(d.genus), "1.genus.sample.txt", sep = "\t", col.names = NA, quote = F)
####genus.melt.long####
data.type.genus.group <- melt(t(t(d.genus.group)),id.vars = 'genus')
data.type.genus.group$Var2 <- gsub("\\.\\.", "-", data.type.genus.group$Var2)
data.type.genus.group$Var2 <- factor(data.type.genus.group$Var2, levels = unique(data.type.genus.group$Var2))
type.genus.group <- ggplot(data.type.genus.group,aes(x = Var1, y = value, alluvium = Var2, stratum = Var2)) +
  geom_alluvium(aes(color = Var2),alpha = 1,width = 0.6, linetype = 1, lwd = 0.5, fill = 'transparent') + 
  geom_bar(aes(x = Var1, y = value, fill = Var2), position = "stack", stat = "identity", width = 0.6) +
  scale_color_manual(values = color3s) +
  scale_fill_manual(values = color3s) +
  labs(x = "", y = "Relative abundance") +
  labs(title = "Genus") +
  ggprism::theme_prism() +
  theme(legend.text = element_text(face = "bold.italic", size = 9),
        legend.key.size = unit(0.7, "lines"),
        legend.key.width = unit(0.45, "lines"),
        plot.title = element_text(size = 11),
        axis.ticks.length = unit(0.1, "cm"),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 9,
                                   color = "black", face = "bold"),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.5)) +
  scale_y_continuous(expand = c(0,0)) 
type.genus.group
ggsave("2.genus.group.pdf", type.genus.group, width = 90, height = 75, units = 'mm')
#species####
species <- aggregate(dat[, c(2:(NCOL(dat) - 7))], by = list(dat$Species), FUN = sum)
row.names(species) <- species$Group.1
species <- dplyr::select(species,-Group.1)
order <- sort(rowSums(species[,1:ncol(species)]),
              index.return = TRUE, decreasing = T)
c.species <- species[order$ix,]
write.table(c.species, "species.all.sample.txt", sep = "\t", col.names = T, quote = F)
####species.Top10####
d.species <- rbind(colSums(c.species[11:as.numeric(length(rownames(c.species))),]),c.species[10:1,])
rownames(d.species)[1] <- "Others"
d.species <- as.data.frame(t(d.species))
d.species$Unassigned <- NULL
d.species$Others <- NULL
merged_table <- merge(d.species, group1, by = "row.names", sort = FALSE)
rownames(merged_table) <- merged_table[, 1]
merged_table <- merged_table[, -1]
d.species.group <- merged_table %>%
  group_by(factor(level, levels = unique(merged_table$level))) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  rename(level = 'factor(level, levels = unique(merged_table$level))')
d.species.group <- as.data.frame(d.species.group)
rownames(d.species.group) <- d.species.group[, 1]
d.species.group <- d.species.group[, -1]
write.table(t(d.species.group), "1.species.group.txt", sep = "\t", col.names = NA, quote = F)
write.table(t(d.species), "1.species.sample.txt", sep = "\t", col.names = NA, quote = F)
####species.melt.long####
data.type.species.group <- melt(t(t(d.species.group)),id.vars = 'species')
process_string <- function(str) {
  str <- as.character(str)
  if (grepl("__", str) && grepl("_", str)) {
    result <- str} else if (grepl("_", str)) {
      result <- gsub("(?<!_)_(?!_)", " ", str, perl = TRUE)} else {
        result <- str}
  result <- as.character(result)
  return(result)}
data.type.species.group$Var2 <- vapply(data.type.species.group$Var2, process_string, FUN.VALUE = character(1), USE.NAMES = FALSE)
data.type.species.group$Var2 <- factor(data.type.species.group$Var2, levels = unique(data.type.species.group$Var2))
type.species.group <- ggplot(data.type.species.group,aes(x = Var1, y = value, alluvium = Var2, stratum = Var2)) +
  geom_alluvium(aes(color = Var2),alpha = 1,width = 0.6, linetype = 1, lwd = 0.5, fill = 'transparent') + 
  geom_bar(aes(x = Var1, y = value, fill = Var2), position = "stack", stat = "identity", width = 0.6) +
  scale_color_manual(values = color3s) +
  scale_fill_manual(values = color3s) +
  labs(x = "", y = "Relative abundance") +
  labs(title = "Species") +
  ggprism::theme_prism() +
  theme(legend.text = element_text(face = "bold.italic", size = 9),
        legend.key.size = unit(0.7, "lines"),
        legend.key.width = unit(0.45, "lines"),
        plot.title = element_text(size = 11),
        axis.ticks.length = unit(0.1, "cm"),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 9,
                                   color = "black", face = "bold"),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.5)) +
  scale_y_continuous(expand = c(0,0)) 
type.species.group
ggsave("3.species.group.pdf", type.species.group, width = 105, height = 75, units = 'mm')
#species.20####
species.20 <- aggregate(dat[, c(2:(NCOL(dat) - 7))], by = list(dat$Species), FUN = sum)
row.names(species.20) <- species.20$Group.1
species.20 <- dplyr::select(species.20,-Group.1)
order <- sort(rowSums(species.20[,1:ncol(species.20)]),
              index.return = TRUE, decreasing = T)
c.species.20 <- species.20[order$ix,]
#write.table(c.species.20, "species.20.all.sample.txt", sep = "\t", col.names = T, quote = F)
####species.20.Top10####
d.species.20 <- rbind(colSums(c.species.20[21:as.numeric(length(rownames(c.species.20))),]),c.species.20[20:1,])
rownames(d.species.20)[1] <- "Others"
d.species.20 <- as.data.frame(t(d.species.20))
d.species.20$Unassigned <- NULL
d.species.20$Others <- NULL
merged_table <- merge(d.species.20, group1, by = "row.names", sort = FALSE)
rownames(merged_table) <- merged_table[, 1]
merged_table <- merged_table[, -1]
d.species.20.group <- merged_table %>%
  group_by(factor(level, levels = unique(merged_table$level))) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  rename(level = 'factor(level, levels = unique(merged_table$level))')
d.species.20.group <- as.data.frame(d.species.20.group)
rownames(d.species.20.group) <- d.species.20.group[, 1]
d.species.20.group <- d.species.20.group[, -1]
write.table(t(d.species.20.group), "1.species.20.group.txt", sep = "\t", col.names = NA, quote = F)
write.table(t(d.species.20), "1.species.20.sample.txt", sep = "\t", col.names = NA, quote = F)
####species.20.melt.long####
data.type.species.20.group <- melt(t(t(d.species.20.group)),id.vars = 'species.20')
process_string <- function(str) {
  str <- as.character(str)
  if (grepl("__", str) && grepl("_", str)) {
    result <- str} else if (grepl("_", str)) {
      result <- gsub("(?<!_)_(?!_)", " ", str, perl = TRUE)} else {
        result <- str}
  result <- as.character(result)
  return(result)}
data.type.species.20.group$Var2 <- vapply(data.type.species.20.group$Var2, process_string, FUN.VALUE = character(1), USE.NAMES = FALSE)
data.type.species.20.group$Var2 <- gsub("\\.\\.", "-", data.type.species.20.group$Var2)
data.type.species.20.group$Var2 <- gsub("_UCG-014", " UCG-014", data.type.species.20.group$Var2)
data.type.species.20.group$Var2 <- factor(data.type.species.20.group$Var2, levels = unique(data.type.species.20.group$Var2))
type.species.20.group <- ggplot(data.type.species.20.group,aes(x = Var1, y = value, alluvium = Var2, stratum = Var2)) +
  geom_alluvium(aes(color = Var2),alpha = 1,width = 0.6, linetype = 1, lwd = 0.5, fill = 'transparent') + 
  geom_bar(aes(x = Var1, y = value, fill = Var2), position = "stack", stat = "identity", width = 0.6) +
  scale_color_manual(values = color4s) +
  scale_fill_manual(values = color4s) +
  labs(x = "", y = "Relative abundance") +
  labs(title = "Species") +
  ggprism::theme_prism() +
  theme(legend.text = element_text(face = "bold.italic", size = 9),
        legend.key.size = unit(0.7, "lines"),
        legend.key.width = unit(0.45, "lines"),
        plot.title = element_text(size = 11),
        axis.ticks.length = unit(0.1, "cm"),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 9,
                                   color = "black", face = "bold"),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.5)) +
  scale_y_continuous(expand = c(0,0)) 
type.species.20.group
ggsave("4.species.20.group.pdf", type.species.20.group, width = 105, height = 75, units = 'mm')


####clean####
#rm(list = ls())
sessionInfo()