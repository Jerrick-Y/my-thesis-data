rm(list = ls())
####获得当前路径#### 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
path <- getwd()
path
####加载包####
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
#genus####
genus.all <- read.delim('core.heat.txt')
#genus.list <- read.delim('2.genus.list.txt')
# 筛选出 genus.all 中存在于 genus.list$genus 中的行
#genus.filtered <- genus.all[genus.all$genus %in% genus.list$genus, ]
#genus.filtered <- genus.filtered[order(rowSums(genus.filtered[,7:NCOL(genus.filtered)]), decreasing = F), ]
# 找到 genus.all 中不存在于 genus.list$genus 中的行，并将它们合并为 "Others"
#genus.others <- genus.all[!genus.all$genus %in% genus.list$genus, ]
#others_row <- c("Others", colSums(genus.others[, -(1:6)], na.rm = TRUE))  # 创建 "Others" 行并计算列的总和
#genus.others <- rbind(genus.others[,-(1:5)], others_row)  
# 合并 "Others" 行和其他行的数据
#genus.heat <- rbind(genus.filtered[,6:NCOL(genus.filtered)], genus.others[NROW(genus.others),])
#genus.heat[,-1] <- sapply(genus.heat[, -1], as.numeric)
#colSums(genus.heat[,-1])
#sapply(genus.heat, class)
row.names(genus.all) <- genus.all$genus
genus.heat <- genus.all[,5:NCOL(genus.all)]
#heatmap.genus####
genus.heat.long <- genus.heat %>%
                    gather(key = adv, value = value, 2:ncol(genus.heat))  %>%
                    mutate(genus = factor(genus)) %>%
                    mutate(adv = factor(adv))
genus.heat.long$adv <- factor(genus.heat.long$adv, levels = unique(genus.heat.long$adv))
genus.heat.long$genus <- factor(genus.heat.long$genus, levels = unique(genus.heat.long$genus))
data.genus <- ggplot(genus.heat.long, aes(x = adv, y = genus, fill = value)) +
  theme_void() + 
  geom_tile(color = "gray96") +
  labs(x = NULL, y = NULL) +  
  scale_fill_gradientn(colors = c("white", "blue", "red"),
                       values = rescale(c(0, 15, 57)),
                       name = NULL, limits = c(0, 57),
                       #breaks = c(0, 0.015, 0.03, 0.25, 0.55),
                       #labels = c(0, 0.015, 0.03, 0.25, 0.55),
                       guide = guide_colorbar(barheight = 8, barwidth = 1,
                                              title.position = "top",
                                              title.hjust = 0.5,
                                              label.position = "right")) +
  #theme_minimal() +
   theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold", angle = 90,hjust = 1),
        axis.text.y = element_text(size = 9, colour = "black", face = "bold.italic", hjust = 1),
        legend.text = element_text(size = 10, face = "bold"), 
        legend.position = "right",
        legend.justification = "left",
        strip.text.y = element_text(angle = -90, hjust = 0.5,size = 9, 
                                    colour = "black", face = "bold", margin = margin(r = 10)),
        strip.background = element_rect(fill = "gray90", color = NA)) +
  geom_text(data = genus.heat.long, 
            aes(label = ifelse(value == 0, "-", 
                               ifelse(value > 10, sprintf("%.1f", value ), 
                                      sprintf("%.2f", value))), 
                color = ifelse(value > 10, "white", "black")), 
            size = 2.5) +
  scale_color_identity() 
  data.genus
ggsave("9.core.heat.pdf", data.genus, width = 5.5, height = 7.5)

####清空####
#rm(list = ls())
