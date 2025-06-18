rm(list = ls())
####获得当前路径#### 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#path <- getwd()
#path
####加载包####
#library(here) # [CRAN] A replacement for 'file.path()', locating the files relative to the project root
library(conflicted)
library(tidyverse)
library(ggplot2)
#library(cowplot)
library(reshape2)
re.genus <- read.delim('../04_top10/genus.all.sample.txt', row.names = 1)
feature.20 <- read.delim('box.letter.txt')
feature.letter <- read.delim('box.letter.txt')
feature.letter <- feature.letter %>% mutate(N.not.0 = paste0("N.not.0: ", N.not.0))
feature.20 <- re.genus[rownames(re.genus) %in% feature.20$adv, ]
feature.20 <- feature.20 %>% mutate(genus = rownames(.)) %>% select(genus, everything())
row.names(feature.20) <- NULL
feature.list <- feature.letter[feature.letter$level == "YYC", ]
feature.list$adv <- factor(feature.list$adv, levels = unique(feature.list$adv), ordered = FALSE)
feature.20 <- feature.20[match( feature.list$adv, feature.20$genus), ]
row.names(feature.20) <- feature.20$genus
write.table(feature.20, 'genus.20.txt', sep = '\t', row.names = F,  quote = F)
write.table(feature.letter, 'genus.letter.txt', sep = '\t', row.names = F,  quote = F)
data.20 <- feature.20
data.20 <- data.20[, -1]
level <- read.delim2('../02_core.file/ASV.group.txt', row.names = 1)
level$level <- factor(level$level,levels = unique(level$level), ordered = TRUE )
#level$copper <- factor(level$copper,levels = unique(level$copper), ordered = TRUE )
level$sample <- rownames(level)
####导入数据####
#转长数据
data.20 <- data.frame(t(data.20))
data.20$sample <- rownames(data.20)
data.20$type <- level$level
#data$copper <- level$copper
data.20$type <- factor(data.20$type,levels = unique(data.20$type), ordered = TRUE )
#data$copper <- factor(data$copper,levels = unique(data$copper), ordered = TRUE )
alpha.a <- melt(data.20, variable.name = "adv", id.vars = names(data.20)[21:ncol(data.20)], measure.vars = 1:20)
alpha.a$type <- factor(alpha.a$type,levels = unique(alpha.a$type), ordered = F )
#alpha.a$copper <- factor(alpha.a$copper,levels = unique(alpha.a$copper), ordered = TRUE )
alpha.a$value <- alpha.a$value *100
#alpha.a$adv <- gsub("\\.L", "-L", alpha.a$adv)
alpha.a$value.1 <- alpha.a$value
alpha.a <- alpha.a %>% 
  mutate(value.1 = ifelse(value == 0, "", as.character(value.1)))
feature.letter <- feature.letter %>% mutate(Phylum = ifelse(level == "YYC", Phylum, ""),
                                          N.not.0 = ifelse(level == "YYC", N.not.0, ""))
alpha.a.stats <- alpha.a %>%
  group_by(adv, type) %>%
  summarise(N = length(value),
            Mean = mean(value),
            #se = sd(value) / sqrt(n()),
            Max = max(value),
            ci_lower = Mean - qt(0.975, n() - 1) * (sd(value) / sqrt(n())),
            ci_upper = Mean + qt(0.975, n() - 1) * (sd(value) / sqrt(n())))
alpha.a.stats$leyel <- alpha.a.stats$type
alpha.letter <-  alpha.a.stats %>%  group_by(adv) %>%  summarise(Max = max(Max))
feature.letter <- merge(feature.letter, alpha.letter[ , c('adv','Max')], by = 'adv')
class(alpha.a$type)
class(alpha.a$adv)
alpha.a$adv <- as.factor(alpha.a$adv)
p <- ggplot(alpha.a, aes(x = type, y = value)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(aes(colour = type), size = 1.6, shape = 16, position = position_jitter(0.2)) +
  geom_text(data = feature.letter, aes(x = level, y = Max * 1.0, label = letter), size = 3) +
  geom_text(data = feature.letter, aes(x = level, y = Max * 1.1, label = Phylum), size = 3) +
  geom_text(data = feature.letter, aes(x = level, y = Max * 1.2, label = N.not.0), size = 3) +
  facet_wrap(~adv, scales = "free", ncol = 4) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.2)), 
                     labels = scales::label_number(accuracy = 0.001)) +
  labs(y = "Relative abundance (%)") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold.italic"),legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-5, -5, -5, -5))
p
ggsave("20.box.pdf",p, width = 8, height = 9.5)
####clean####
rm(list = ls())
