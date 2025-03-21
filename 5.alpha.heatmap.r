rm(list = ls())
####获得当前路径#### 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(11111)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(ggh4x)
data.alpha.all <- read.delim('alpha.heat.txt')
data.label.all <- read.delim('alpha.label.txt')
data.alpha <- data.alpha.all[,1:6]
data.label <- data.label.all[,1:6]

data.alpha.long <- data.alpha %>%
  gather(key = adv, value = value, 4:ncol(data.alpha))  %>%
  mutate(alpha = factor(alpha)) %>%
  mutate(adv = factor(adv))
data.label.long <- data.label %>%
  gather(key = adv, value = value, 4:ncol(data.alpha))  %>%
  mutate(alpha = factor(alpha)) %>%
  mutate(adv = factor(adv))
data.alpha.long$label <- data.label.long$value
means <- data.alpha.long %>%
  group_by(alpha) %>%
  summarise(mean_value = mean(value))
sds <- data.alpha.long %>%
  group_by(alpha) %>%
  summarise(sd_value = sd(value))

# 将均值和标准差合并到原始数据框中
data.alpha.long <- left_join(data.alpha.long, means, by = "alpha")
data.alpha.long <- left_join(data.alpha.long, sds, by = "alpha")

# 计算 z-score
data.alpha.long <- data.alpha.long %>%
  mutate(z_score = (value - mean_value) / sd_value)

data.alpha.long$adv <- factor(data.alpha.long$adv, levels = unique(data.alpha.long$adv))
data.alpha.long$type <- factor(data.alpha.long$type, levels = unique(data.alpha.long$type))
data.alpha.long$alpha <- factor(data.alpha.long$alpha, levels = unique(data.alpha.long$alpha))
data.p <- ggplot(data.alpha.long, aes(x = adv, y = type, fill = z_score)) +
  theme_void() + 
  geom_tile(color = "gray96") +
  scale_fill_gradient2(low = "blue", mid = "grey90", high = "red", midpoint = 0,
                       name = NULL, limits = c(-1.5, 1.81),
                       breaks = c(-1.4, -1.0, -0.5, 0, 0.5, 1, 1.5),
                       labels = c(-1.4, -1.0, -0.5, 0, 0.5, 1, 1.5),
                       guide = guide_colorbar(barheight = 12)) +
  labs(x = NULL, y = NULL) +  
  facet_wrap(~alpha, ncol = 1, strip.position = "right") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 11, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 10, colour = "black", face = "bold", hjust = 1),
        legend.text = element_text(size = 10, face = "bold"), 
        legend.position = "right",
        legend.justification = "left",
        strip.text.y = element_text(angle = -90, hjust = 0.5,size = 10, 
                                    colour = "black", face = "bold", margin = margin(r = 10)),
        strip.background = element_rect(fill = "gray90", color = NA)) +
  geom_text(data = subset(data.alpha.long), aes(label = label), color = "black", size = 4, fontface = "bold") 
data.p
ggsave("5.1.alpha.heat.pdf", data.p, width = 2.7, height = 4.5)