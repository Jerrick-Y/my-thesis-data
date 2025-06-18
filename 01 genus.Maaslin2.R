rm(list = ls())
####获得当前路径#### 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(11111)
library(Maaslin2)
library(dplyr)
library(tidyverse)
abasv <- read.delim2('../02_core.file/ASV.new.txt', row.names = 1)
tax <- abasv %>% select(-contains("_"))
tax <- tax %>% separate(taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") 
tax$Phylum <- gsub('(p__)','',tax$Phylum)
tax$Genus <- gsub('(g__)','',tax$Genus)
dat <- merge(x = abasv[, 1:12], y = tax, by = "row.names")
dat <- dplyr::rename(dat, OTUID = Row.names)
genus <- aggregate(dat[, 2:13], by = list(dat$Genus), FUN = sum)
row.names(genus) <- genus$Group.1
genus <- dplyr::select(genus,-Group.1)
input_data <- genus
input_data <- data.frame(lapply(input_data, as.numeric), row.names = rownames(input_data))
class(input_data$C1_1)
#colSums(input_data)
#rownames(input_data) <- rownames(genus)
#rowSums(input_data != 0)
input_metadata <- read.delim2('../02_core.file/ASV.group.txt', row.names = 1)
input_metadata$level <- as.character(input_metadata$level)
fit_data <- Maaslin2(input_data, input_metadata, 'C.vs2',  min_prevalence = 0,
                     min_abundance = 0, #max_significance = 0.05,
                     transform = "LOG", normalization  = "TSS",
                     fixed_effects = c('level'), #random_effects = c('tank'),
                     reference = "level,Con", standardize = FALSE)
fit_data.YYC <- Maaslin2(input_data, input_metadata, 'YYC.vs2', min_prevalence = 0,
                       min_abundance = 0, #max_significance = 0.05,
                       transform = "LOG", normalization  = "TSS",
                       fixed_effects = c('level'), #random_effects = c('tank'),
                       reference = "level,YYC", standardize = FALSE)
fit_data.YKN <- Maaslin2(input_data, input_metadata, 'YKN.vs2', min_prevalence = 0,
                         min_abundance = 0, #max_significance = 0.05,
                         transform = "LOG", normalization  = "TSS",
                         fixed_effects = c('level'), #random_effects = c('tank'),
                         reference = "level,YKN", standardize = FALSE)
fit_data <- read.delim('./C.vs2/all_results.tsv')
#fit_data$feature <- gsub("\\.L", "-L", fit_data$feature)
fit_data.YYC <- read.delim('./YYC.vs2/all_results.tsv')
#fit_data.YYC$feature <- gsub("\\.L", "-L", fit_data.YYC$feature)
fit_data.YKN <- read.delim('./YKN.vs2/all_results.tsv')
#fit_data.YKN$feature <- gsub("\\.L", "-L", fit_data.YKN$feature)
library(dplyr)
genus.27 <- fit_data
genus.27$ref <- c('con')
#genus.27.filter <- genus.27[genus.27$feature %in% rownames(feature.66), ]
genus.27.filter <- genus.27 %>% mutate(qval = ifelse(qval > 0.25, 1, qval))
genus.27.YYC <- fit_data.YYC
genus.27.YYC$ref <- c('YYC')
#genus.27.GS.filter <- genus.27.GS[genus.27.GS$feature %in% rownames(feature.66), ]
genus.27.YYC.filter <- genus.27.YYC %>% mutate(qval = ifelse(qval > 0.25, 1, qval))
genus.27.YYC.filter <- genus.27.YYC.filter %>% filter(value != "Con")

genus.27.YKN <- fit_data.YKN
genus.27.YKN$ref <- c('YKN')
#genus.27.YKN.filter <- genus.27.YKN[genus.27.YKN$feature %in% rownames(feature.66), ]
genus.27.YKN.filter <- genus.27.YKN %>% mutate(qval = ifelse(qval > 0.25, 1, qval))
genus.27.YKN.filter <- genus.27.YKN.filter %>% filter(value != "WC") %>% filter(value != "WS")
filter.0.25 <- rbind(genus.27.filter, genus.27.YYC.filter)
# 首先，找出feature列中出现3次的值
feature_counts <- filter.0.25 %>%
  count(feature) %>%
  filter(n == 3) %>%
  pull(feature)
feature_counts
# 然后，删除满足条件的行
filter.0.25.3 <- filter.0.25 %>%
  group_by(feature) %>%
  filter(!feature %in% feature_counts | sum(qval) != 3) %>%
  ungroup()
phylum.genus <- tax[, c('Phylum','Genus')]
row.names(phylum.genus) <- NULL
phylum.genus <- phylum.genus %>%
  distinct(Genus, .keep_all = TRUE)
filter.0.25.3.1 <- filter.0.25.3 %>%
  left_join(phylum.genus %>% select(Genus, Phylum), by = c("feature" = "Genus"))
write.table(filter.0.25.3.1, file = 'filter.0.25.txt', quote = F, sep = '\t', row.names = F)