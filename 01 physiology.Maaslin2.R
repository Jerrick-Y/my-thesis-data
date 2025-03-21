rm(list = ls())
####获得当前路径#### 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(11111)
library(Maaslin2)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(forcats)
ab.genus <- read.delim('ab.CSO2O3.txt', row.names = 1)
ab.genus <- select(ab.genus, -matches("O2.", ignore.case = TRUE))
phylums <- ab.genus[,-(1:54)]
re.genus <- ab.genus[, 1:54]
row_names <- row.names(re.genus)
input_data <- re.genus
input_data <- apply(input_data, 2, as.numeric)
colSums(input_data)
rownames(input_data) <- rownames(re.genus)
rowSums(input_data != 0)
input_metadata <- read.delim2('CSO.group.txt', row.names = 1)
input_metadata$level <- as.character(input_metadata$level)
input_metadata$type <- as.character(input_metadata$type)
#input_metadata$tank <- as.character(input_metadata$tank)
input_metadata$Group <- as.character(input_metadata$Group)

#input_metadata$copper2 <- as.character(input_metadata$copper2)
#input_metadata$copper3 <- as.character(input_metadata$copper3)
input_metadata$Physiology <- as.numeric(input_metadata$Physiology)
#input_metadata$OMn <- as.numeric(input_metadata$OMn)
#O_cleaned <- input_metadata[complete.cases(input_metadata$O), ]
#OM_cleaned <- input_metadata[complete.cases(input_metadata$OMn), ]
#O1_cleaned <- input_metadata[complete.cases(input_metadata$O1.vs), ]
#O2_cleaned <- input_metadata[complete.cases(input_metadata$O2.vs), ]
#O3_cleaned <- input_metadata[complete.cases(input_metadata$O3.vs), ]
#SO1_cleaned <- input_metadata[complete.cases(input_metadata$SO1), ]
#SO2_cleaned <- input_metadata[complete.cases(input_metadata$SO2), ]

fit_data <- Maaslin2(input_data, input_metadata, 'physio',  min_prevalence = 0,
                     min_abundance = 0, #max_significance = 0.05,
                     transform = "LOG", normalization  = "TSS",
                     fixed_effects = c('Physiology'), #random_effects = c('tank'),
                     standardize = FALSE)
physio.all <- fit_data$results %>% rename_with(~ gsub("zero", "0", .)) #%>% filter(qval <= 0.25) 
physio.all$feature <- gsub("\\.", "_", physio.all$feature)
feature.list <- read.delim('feature.lis.109.txt')
physio.109 <-  physio.all[physio.all$feature %in% feature.list$feature, ]
physio.109$qval[physio.109$qval > 0.25] <- 1
physio.109$score <- -log(physio.109$qval) * sign(physio.109$coef)
physio.46 <- physio.109[physio.109$feature %in% Masslin2.45$feature, ]
physio.60 <- physio.109[physio.109$feature %in% Masslin2.60$feature, ]
physio.46.bind <- physio.46 %>% rename_with(~ gsub("name", "Ref", .)) %>% select(-N, -N.not.0)
physio.46.bind <- physio.46.bind[, names(Masslin2.45)]
Masslin2.45.all <- rbind(Masslin2.45, physio.46.bind)
Masslin2.45.all$metadata <- gsub("Physiology", "D", Masslin2.45.all$metadata)
Masslin2.45.all$feature <- factor(Masslin2.45.all$feature, levels = feature.46$genus)


physio.60.bind <- physio.60 %>% rename_with(~ gsub("name", "Ref", .)) %>% select(-N, -N.not.0)
physio.60.bind <- physio.60.bind[, names(Masslin2.60[,-10])]
physio.60.bind <- physio.60.bind %>%
  left_join(Masslin2.60 %>% select(feature, phylum), by = c("feature" = "feature"))
Masslin2.60.all <- rbind(Masslin2.60, physio.60.bind)


Masslin2.60.all <- arrange(Masslin2.60.all, feature)
Masslin2.60.all <- arrange(Masslin2.60.all, phylum)
Masslin2.60.all$feature <- factor(Masslin2.60.all$feature, unique(Masslin2.60$feature))
# 按照因子的顺序排序
Masslin2.45.all <- Masslin2.45.all[order(Masslin2.45.all$feature), ]
Masslin2.45.all$feature <- factor(Masslin2.45.all$feature)
Masslin2.45.all$feature <- fct_rev(Masslin2.45.all$feature)
#Maaslin2.plot
# 创建所有可能的 metadata 和 feature 组合
all_combinations <- expand.grid(
  metadata = unique(Masslin2.45.all$metadata),
  feature = unique(Masslin2.45.all$feature)
)

# 将所有组合与现有数据合并
Masslin2.45_full <- all_combinations %>%
  left_join(Masslin2.45.all, by = c("metadata", "feature"))

# 替换 NA 值为灰色
Masslin2.45_full$score <- ifelse(is.na(Masslin2.45_full$score), 0, Masslin2.45_full$score)

# 将分数列转换为因子，以便正确处理颜色填充
Masslin2.45_full$score <- as.numeric(Masslin2.45_full$score)

counts <- Masslin2.45_full %>%
  mutate(sign = ifelse(score > 0, "+", ifelse(score < 0, "-", "0"))) %>%
  group_by(metadata, sign) %>%
  summarise(count = n()) %>%
  spread(sign, count, fill = 0)

# 组合 "+" 和 "-" 数量成字符
counts <- counts %>%
  mutate(label = paste0("+", `+`, "\n-", `-`))
feature.p <- ggplot(Masslin2.45_full, aes(x = metadata, y = feature, fill = score)) +
  geom_tile(color = "gray96") +  # 设置瓦片颜色为填充值
  scale_fill_gradient2(low = "blue", mid = "grey90", high = "red", 
                       name = NULL, limits = c(-3.7, 4.2),
                       breaks = c(-3.6, -2,  0, 2, 4),
                       labels = c(-3.6, -2, 0, 2, 4),
                       guide = guide_colorbar(barheight = 12)) +
  theme_void() +  
  theme(axis.text.x = NULL,
        axis.text.y = element_text(size = 10, colour = "black", face = "bold.italic", hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2)) + 
  theme(legend.text = element_text(size = 10, face = "bold"), 
        legend.position = "left", 
        legend.justification = "left") +
  geom_text(data = subset(Masslin2.45_full, score > 0), aes(label = "+"), color = "black", size = 4, fontface = "bold") +
  geom_text(data = subset(Masslin2.45_full, score < 0 & score != 0), aes(label = "-"), color = "black", size = 3, fontface = "bold")
feature.p
ggsave("feature.p.pdf", feature.p, width = 4.2, height = 7.1)

library(patchwork)
library(gridExtra)
library(grid)
# 设置每个图表的宽度
widths <- unit(c(3.6, 1.7, 4.1), "null")
# 使用 grid.arrange 并设置宽度和高度
p.1 <- grid.arrange(
     arrangeGrob(feature.p, nrow = 1), 
     arrangeGrob(line.p, nrow = 1), 
    arrangeGrob(function.p, nrow = 1), 
     ncol = 3, 
    widths = widths)
ggsave("all.1.phy.2.pdf", p.1, width = 11.2, height = 8)

#Maaslin2.60.plot
# 创建所有可能的 metadata 和 feature 组合
Masslin2.60.all$feature <- factor(Masslin2.60.all$feature)
Masslin2.60.all$feature <- fct_rev(Masslin2.60.all$feature)

all_combinations.60 <- expand.grid(
  metadata = unique(Masslin2.60.all$metadata),
  feature = unique(Masslin2.60.all$feature)
)

# 将所有组合与现有数据合并
Masslin2.60_full <- all_combinations.60 %>%
  left_join(Masslin2.60.all, by = c("metadata", "feature"))

# 替换 NA 值为灰色
Masslin2.60_full$score <- ifelse(is.na(Masslin2.60_full$score), 0, Masslin2.60_full$score)

# 将分数列转换为因子，以便正确处理颜色填充
Masslin2.60_full$score <- as.numeric(Masslin2.60_full$score)

counts.60 <- Masslin2.60_full %>%
  mutate(sign = ifelse(score > 0, "+", ifelse(score < 0, "-", "0"))) %>%
  group_by(metadata, sign) %>%
  summarise(count = n()) %>%
  spread(sign, count, fill = 0)

# 组合 "+" 和 "-" 数量成字符
counts.60 <- counts.60 %>%
  mutate(label = paste0("+", `+`, "\n-", `-`))
Masslin2.60_full$feature <- factor(Masslin2.60_full$feature)

feature.p.60 <- ggplot(Masslin2.60_full, aes(x = metadata, y = feature, fill = score)) +
  geom_tile(color = "gray96") +  # 设置瓦片颜色为填充值
  scale_fill_gradient2(low = "blue", mid = "grey90", high = "red", 
                       name = NULL, limits = c(-3.9, 4.2),
                       breaks = c(-3.8, -2,  0, 2, 4),
                       labels = c(-3.8, -2,  0, 2, 4),
                       guide = guide_colorbar(barheight = 12)) +
  theme_void() +  
  theme(axis.text.x = NULL,
        axis.text.y = element_text(size = 9, colour = "black", face = "bold.italic", hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2)) + 
  theme(legend.text = element_text(size = 9, face = "bold"), 
        legend.position = "left", 
        legend.justification = "left") +
  geom_text(data = subset(Masslin2.60_full, score > 0), aes(label = "+"), color = "black", size = 3.5, fontface = "bold") +
  geom_text(data = subset(Masslin2.60_full, score < 0 & score != 0), aes(label = "-"), color = "black", size = 3.5, fontface = "bold")
feature.p.60
ggsave("feature.p.60.4.pdf", feature.p.60, width = 5.1, height = 9)

