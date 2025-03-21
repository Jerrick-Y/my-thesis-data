rm(list = ls())
####获得当前路径#### 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(11111)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(forcats)
#read genus table#
otu.tss <- read.delim('CSO3.TSS.txt')
#read  maaslin2#
Maaslin2.all <- read.delim('Masslin2.all.txt')
#read list of maaslin2#
Maaslin2.list <- read.delim('feature.lis.109.txt')
#read list of function#
feature.46 <- read.delim('feature.list.46.txt')
#read function#
function.46 <- read.delim('function.list.46.txt')
#long shape
function.46.long <- function.46 %>%
  melt(id.vars = c(1:3), 
       measure.vars = c(4:ncol(function.46)), 
       variable.name = "Sample", 
       value.name = "functions") %>%
       filter(!apply(. == "", 1, any)) %>%
        select(-Sample)
#function.label
function.label <- read.delim('function.label.txt')
#add label to function.46
# 使用 left_join 将 function.label 的标记添加到 function.46.long 中
function.label$functions <- gsub(" ", "_", trimws(function.label$functions))
function.46.long$functions <- gsub(" ", "_", trimws(function.46.long$functions))
function.46.long <- function.46.long %>% 
                    left_join(function.label, by = "functions", suffix = c("", ".label"))
# 创建一个因子变量，以 feature.46$genus 的顺序排序
function.46.long$genus <- factor(function.46.long$genus, levels = feature.46$genus)
# 按照因子的顺序排序
function.46.long <- function.46.long[order(function.46.long$genus), ]
write.table(function.46.long, file = 'function.46.long.txt', quote = F, sep = '\t', row.names = F)

#filter genus table#
samples.105 <- otu.tss[otu.tss$genus %in% Maaslin2.list$feature, ]
samples.45 <- samples.105[samples.105$genus %in% feature.46$genus, ]
samples.60 <- samples.105[!samples.105$genus %in% feature.46$genus, ]
Masslin2.45 <- Maaslin2.all[Maaslin2.all$feature %in% samples.45$genus, ]
Masslin2.60 <- Maaslin2.all[Maaslin2.all$feature %in% samples.60$genus, ]
Masslin2.60 <- Masslin2.60 %>%
               left_join(samples.105 %>% select(genus, phylum), by = c("feature" = "genus")) %>%
                arrange(feature, phylum)
Masslin2.60 <- arrange(Masslin2.60, feature)
Masslin2.60 <- arrange(Masslin2.60, phylum)
Masslin2.60$feature <- factor(Masslin2.60$feature, unique(Masslin2.60$feature))
#calculate score#
Masslin2.45$score <- -log(Masslin2.45$qval) * sign(Masslin2.45$coef)
Masslin2.60$score <- -log(Masslin2.60$qval) * sign(Masslin2.60$coef)
Masslin2.45$feature <- factor(Masslin2.45$feature, levels = feature.46$genus)
# 按照因子的顺序排序
Masslin2.45 <- Masslin2.45[order(Masslin2.45$feature), ]
Masslin2.45$feature <- factor(Masslin2.45$feature)
Masslin2.45$feature <- fct_rev(Masslin2.45$feature)
#Maaslin2.plot
# 创建所有可能的 metadata 和 feature 组合
all_combinations <- expand.grid(
  metadata = unique(Masslin2.45$metadata),
  feature = unique(Masslin2.45$feature)
)

# 将所有组合与现有数据合并
Masslin2.45_full <- all_combinations %>%
  left_join(Masslin2.45, by = c("metadata", "feature"))

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
                       name = NULL, limits = c(-3.7, 3.8),
                       breaks = c(-3.6, -1.5,  0, 1.5, 3.7),
                       labels = c(-3.6, -1.5, 0, 1.5, 3.7),
                       guide = guide_colorbar(barheight = 12)) +
  theme_void() +  
  theme(axis.text.x = NULL,
        axis.text.y = element_text(size = 10, colour = "black", face = "bold.italic", hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2)) + 
  theme(legend.text = element_text(size = 10, face = "bold"), 
        legend.position = "left", 
        legend.justification = "left") +
  geom_text(data = subset(Masslin2.45_full, score > 0), aes(label = "+"), color = "black", size = 3, fontface = "bold") +
  geom_text(data = subset(Masslin2.45_full, score < 0 & score != 0), aes(label = "-"), color = "black", size = 3, fontface = "bold")
feature.p
ggsave("feature.p.pdf", feature.p, width = 4.2, height = 7.1)
# function数据框是
function.36 <- function.46.long[,c(1,3,4,5)]
# 创建一个新的数据框，用于绘制线段
#line

function.list <- read.delim('functions.36.txt')
function.list$value <- rep(1)
function.list$x <- rep('D')
function.list$functions <- factor(function.list$functions, levels = rev(unique(function.list$functions)))
function.p <- ggplot(function.list, aes(x = x, y = functions, fill = labels)) +
  geom_tile(color = "gray96") +  # 设置瓦片颜色为填充值
  theme_void() +  
  scale_y_discrete(position = "right") +  # 将 y 轴标签移到右侧
  theme(axis.text.x = NULL,
        axis.text.y.right = element_text(size = 10, colour = "black", face = "bold", hjust = 0), # 右侧y轴文本
        legend.text = element_text(size = 10, face = "bold"), 
        legend.position = "right", 
        legend.justification = "right")
function.p
ggsave("function.pp.pdf", function.p, width = 4.5, height = 6)
y1 <- ggplot_build(feature.p)$data[[1]][,c('x','y')]
axis_left <- data.frame(features = Masslin2.45_full$feature, y1 = y1$y)
y2 <- ggplot_build(function.p)$data[[1]][,c('x','y')]
axis_right <- data.frame(functions = function.list$functions,y2 = y2$y)
df_line <- function.46.long[,c(1,3,4,5)]
df_line <- merge(df_line, axis_right, by.x = 'functions', by.y = 'functions')
df_line <- merge(df_line, axis_left, by.x = 'genus', by.y = 'features')
df_line$x1 <- rep(2)
df_line$x2 <- rep(4)
#df_line$y2 <- df_line$y2 + 4
df_line <- df_line %>%
  mutate(y2 = case_when(
    y2 > 24 ~ y2 + 9,
    y2 >= 16 & y2 <= 24 ~ y2 + 6,
    y2 >= 11 & y2 <= 15 ~ y2 + 3,
    y2 <= 10 ~ y2))
color30 <- c('#acc2d9',  '#56ae57','#b2996e', '#894585', 
             '#d4ffff', '#388004', '#efb435', "#00bFcF", '#ffd8b1', '#1f6357',
             '#acc2d9',  '#56ae57','#b2996e', '#894585', 
             '#d4ffff', '#388004', '#efb435', "#00bFcF", '#ffd8b1', '#1f6357',
             '#acc2d9',  '#56ae57','#b2996e', '#894585', 
             '#d4ffff', '#388004', '#efb435', "#00bFcF", '#ffd8b1', '#1f6357',
             '#acc2d9',  '#56ae57','#b2996e', '#894585', 
             '#d4ffff', '#388004', '#efb435', "#00bFcF", '#ffd8b1', '#1f6357',
             '#acc2d9',  '#56ae57','#b2996e', '#894585', 
             '#d4ffff', '#388004', '#efb435', "#00bFcF", '#ffd8b1', '#1f6357',
             '#3778bf', '#ff0789', '#a9a9a9', '#430541',
             '#ffb2d0', '#ad900d', '#f6688e', '#850e04', '#aaffc3', '#f97306',
             '#76fda8',  '#41fdfe', '#0c1793', '#a50055', '#ad03de',
             '#aeff6e',  '#fffd01', '#0165fc', '#f97306')
df_line <- df_line[order(df_line$y2, decreasing = TRUE), ]

# 重新定义labels的因子水平
df_line$labels.1 <- factor(df_line$functions, levels = unique(df_line$functions))
line.p <- ggplot(df_line) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, color = labels.1), linewidth = 0.5, show.legend = FALSE) +
  geom_point(aes(x = x1, y = y1, fill = phylum), size = 3, shape = 21, show.legend = FALSE) +
  geom_point(aes(x = x2, y = y2, fill = labels), size = 3, shape = 21, show.legend = FALSE) +
  theme_void() + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_color_manual(values = color30) +  # 使用自定义的颜色
  theme(legend.text = element_text(size = 10, face = "bold"), 
        legend.position = "right")
line.p  
ggsave("line.p.pdf", line.p, width = 1, height = 6.8)

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
ggsave("all.1.pdf", p.1, width = 10.7, height = 8)

#Maaslin2.60.plot
# 创建所有可能的 metadata 和 feature 组合
Masslin2.60$feature <- factor(Masslin2.60$feature)
Masslin2.60$feature <- fct_rev(Masslin2.60$feature)

all_combinations.60 <- expand.grid(
  metadata = unique(Masslin2.60$metadata),
  feature = unique(Masslin2.60$feature)
)

# 将所有组合与现有数据合并
Masslin2.60_full <- all_combinations.60 %>%
  left_join(Masslin2.60, by = c("metadata", "feature"))

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
                       name = NULL, limits = c(-3.9, 4.0),
                       breaks = c(-3.8, -2,  0, 2, 3.9),
                       labels = c(-3.8, -2,  0, 2, 3.9),
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
ggsave("feature.p.60.pdf", feature.p.60, width = 4.2, height = 7.1)
y1.60 <- ggplot_build(feature.p.60)$data[[1]][,c('x','y')]
axis_left.60 <- data.frame(features = Masslin2.60_full$feature, y1 = y1.60$y)
df_line.60 <- Masslin2.60_full[,c(2,10)]
df_line.60 <- merge(df_line.60, axis_left.60, by.x = 'feature', by.y = 'features')
df_line.60$x1 <- rep(2)
df_line.60 <- df_line.60 %>% filter(complete.cases(df_line.60))
df_line.60 <- df_line.60 %>% distinct()
color31 <- c('#acc2d9',  '#56ae57','#b2996e', '#894585', '#efb435', "#00bFcF", '#f97306', '#1f6357',
             '#3778bf',  '#430541',
             '#ffb2d0', '#ad900d', '#f6688e', '#850e04', '#aaffc3', '#f97306',
             '#76fda8',  '#41fdfe', '#0c1793', '#a50055', '#ad03de',
             '#aeff6e',  '#fffd01', '#0165fc', '#f97306')

line.p.60 <- ggplot(df_line.60) +
  geom_point(aes(x = x1, y = y1, fill = phylum, color = phylum), size = 3, shape = 21, show.legend = T) +
  theme_void() + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_color_manual(values = color31) + 
  scale_fill_manual(values = color31) +  # 使用自定义的颜色
  theme(legend.text = element_text(size = 10, face = "bold"), 
        legend.position = "right")
line.p.60  
library(patchwork)
library(gridExtra)
library(grid)
# 设置每个图表的宽度
widths <- unit(c(5, 2), "null")
# 使用 grid.arrange 并设置宽度和高度
p.2 <- grid.arrange(
  arrangeGrob(feature.p.60, nrow = 1), 
  arrangeGrob(line.p.60, nrow = 1), 
  ncol = 2, 
  widths = widths)
ggsave("all.60.pdf", p.2, width = 6.5, height = 9)
####boxplot####
data.box <- select(samples.105, -matches("F|WR", ignore.case = F))
data.box$letter <- apply(data.box[, 7:ncol(data.box)], 1, function(x) sum(x != 0))
data.box <- data.box %>% mutate(letter = paste0("N.not.0: ", letter, "/54"))
data.box$max <- apply(data.box[, c(7:(ncol(data.box)-1))], 1, max, na.rm = TRUE)
data.box.level <- read.delim('CSO.box.group.txt')
#data
box.data <- data.frame(t(data.box[,-c(1:5,ncol(data.box),(ncol(data.box)-1))]))
# 将第一行设置为列名
colnames(box.data) <- as.character(box.data[1, ])
# 移除现在已经作为列名的第一行
box.data <- box.data[-1, ]
box.data <- box.data %>% mutate(sample = rownames(.)) %>% select(sample, everything())
box.data.letter <- data.box[,c('genus','letter','max')]
box.data.letter$level <- c('S30')  
#box.data <- box.data[,-c(2:4,110:112)]
box.data <- left_join(data.box.level, box.data, by = 'sample')
box.data <-  box.data %>% 
  mutate(type = factor(type, levels = unique(box.data$type))) %>%
  mutate(level = factor(level, levels = unique(box.data$level))) %>%
  mutate(copper = factor(copper, levels = unique(box.data$copper)))
box.data.a <- box.data %>% gather(key = adv, value = value, 5:ncol(box.data))  %>%
  mutate(value = as.numeric(value)) 
box.data.a <- box.data.a %>% 
  mutate(type = factor(type, levels = unique(box.data.a$type))) %>%
  mutate(level = factor(level, levels = unique(box.data.a$level))) %>%
  mutate(copper = factor(copper, levels = unique(box.data.a$copper)))
box.data.a$value1 <- box.data.a$value*100
#class(box.data.a$value)
box.data.a.stats <- box.data.a %>%
  group_by(adv, copper) %>%
  summarise(N = length(value1),
            Mean = mean(value1),
            se = sd(value1) / sqrt(n()),
            Max = max(value1),
            ci_lower = Mean - qt(0.975, n() - 1) * (sd(value1) / sqrt(n())),
            ci_upper = Mean + qt(0.975, n() - 1) * (sd(value1) / sqrt(n())))
box.data.letter$max2 <- box.data.letter$max*100
#box.data.a <- box.data.a %>% mutate(size = ifelse(value1 == 0, 0.2, 0.8))
p <- ggplot(box.data.a, aes(x = copper, y = value1)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(aes(colour = type), size = 1.6, shape = 16, position = position_jitter(0.2)) +
  facet_wrap(~adv, scales = "free", ncol = 5) +
  geom_text(data = box.data.letter, aes(x = level, y = max2 * 1.0, label = letter), size = 3) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.2)), 
                     labels = scales::label_number(accuracy = 0.01)) +
  labs(y = "Relative abundance (%)") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold.italic"),legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-5, -5, -5, -5))
#p
ggsave("105.box.2.pdf",p, width = 9.5, height = 37, dpi = 3000, limitsize = FALSE)
#dev.off()
####core.feature####
core.list <- read.delim('core.28.list.txt')
data.core <- otu.tss[otu.tss$genus %in% core.list$genus, ]
data.core <- select(data.core, -matches("F|WR", ignore.case = F))
data.core$letter <- apply(data.core[, 7:ncol(data.core)], 1, function(x) sum(x != 0))
data.core <- data.core %>% mutate(letter = paste0("N.not.0: ", letter, "/54"))
data.core$max <- apply(data.core[, c(7:(ncol(data.core)-1))], 1, max, na.rm = TRUE)
data.core.level <- read.delim('CSO.box.group.txt')
#data
core.data <- data.frame(t(data.core[,-c(1:5,ncol(data.core),(ncol(data.core)-1))]))
# 将第一行设置为列名
colnames(core.data) <- as.character(core.data[1, ])
# 移除现在已经作为列名的第一行
core.data <- core.data[-1, ]
core.data <- core.data %>% mutate(sample = rownames(.)) %>% select(sample, everything())
core.data.letter <- data.core[,c('genus','letter','max')]
core.data.letter$level <- c('S30')  
core.data.letter$adv <- core.data.letter$genus
#box.data <- box.data[,-c(2:4,110:112)]
core.data <- left_join(data.core.level, core.data, by = 'sample')
core.data <-  core.data %>% 
  mutate(type = factor(type, levels = unique(core.data$type))) %>%
  mutate(level = factor(level, levels = unique(core.data$level))) %>%
  mutate(copper = factor(copper, levels = unique(core.data$copper)))
core.data.a <- core.data %>% gather(key = adv, value = value, 5:ncol(core.data))  %>%
  mutate(value = as.numeric(value)) #
core.data.a <- core.data.a %>% 
  mutate(type = factor(type, levels = unique(core.data.a$type))) %>%
  mutate(level = factor(level, levels = unique(core.data.a$level))) %>%
  mutate(copper = factor(copper, levels = unique(core.data.a$copper)))
core.data.a$value1 <- core.data.a$value*100
#class(core.data.a$value)
core.data.a.stats <- core.data.a %>%
  group_by(adv, copper) %>%
  summarise(N = length(value1),
            Mean = mean(value1),
            se = sd(value1) / sqrt(n()),
            Max = max(value1),
            ci_lower = Mean - qt(0.975, n() - 1) * (sd(value1) / sqrt(n())),
            ci_upper = Mean + qt(0.975, n() - 1) * (sd(value1) / sqrt(n())))
core.data.letter$max2 <- core.data.letter$max*100
#core.data.a <- core.data.a %>% mutate(size = ifelse(value1 == 0, 0.2, 0.8))
p.core <- ggplot(core.data.a, aes(x = copper, y = value1)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(aes(colour = type), size = 1.6, shape = 16, position = position_jitter(0.2)) +
  facet_wrap(~adv, scales = "free", ncol = 5) +
  geom_text(data = core.data.letter, aes(x = level, y = max2 * 1.0, label = letter), size = 3) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.2)), 
                     labels = scales::label_number(accuracy = 0.01)) +
  labs(y = "Relative abundance (%)") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold.italic"),legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-5, -5, -5, -5))
#p.core
ggsave("core.28.box.pdf", p.core, width = 9.5, height = 11, dpi = 3000, limitsize = FALSE)
