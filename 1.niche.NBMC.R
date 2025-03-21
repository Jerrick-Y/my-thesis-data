####生态位宽度####
library(spaa)
library(EcolUtils)
library(vegan)
library(dplyr)
genus.rela <- read.delim('feature.species.2.txt', row.names = 1)
tax <- read.delim('../../1.ASV_table/tax/ASV_table.tax.txt', row.names = 1)
unique_tax <- tax[!duplicated(tax$Species), ]
asv.taxonomy <- unique_tax
write.table(asv.taxonomy, 'tax.species.merge.txt', quote = F, sep = '\t', row.names = F)
row.names(asv.taxonomy) <- asv.taxonomy$Species
group.1 <- read.delim('group.2.txt')
group.1 <- group.1[#(group$type == "Water" | group$type == "WS") & 
  #group.1$Diets == "C" & group.1$Mineral == "Fe" & 
  group.1$type2 != '' , ]
group.1 <- group.1 %>% 
  mutate(Group = case_when(
    row_number() <= 12 | (row_number() >= 73 & row_number() <= 84) ~ paste0('C', Group),
    TRUE ~ Group
  ))
row.names(group.1) <- group.1$Sample
selected_samples <- group.1$Sample
group.1$Group <- factor(group.1$Group, levels = unique(group.1$Group))
genus.rela <- genus.rela[, colnames(genus.rela) %in% selected_samples, drop = FALSE]
genus.rela <- genus.rela[, match(selected_samples, colnames(genus.rela))]
dune <- t(genus.rela[,1:NCOL(genus.rela)])
niche_width <- niche.width(dune, method = 'levins')
write.table(t(niche_width), 'niche_width.mean.txt', sep = '\t', row.names = T,col.names = F, quote = FALSE)
boxplot(unlist(niche_width), ylab = 'niche breadth index')
#set.seed(123)
dune10 <- dune * 10000000
dune10 <- round(dune10, 0)
dune.group <- aggregate(dune10, by = list(group.1$Group), FUN = sum)
rownames(dune.group) <- dune.group$Group.1
dune.group$Group.1 <- NULL
RS <- t(dune.group)
rowSumRS <- rowSums(RS)
#rowSumRS
a <- RS / rowSumRS
a2 <- a * a
rowSuma2 <- rowSums(a2)
niche <- 1 / rowSuma2
niche.A <- (niche - 1) / (NCOL(RS) - 1)
class(niche.A)
niche.A_df <- as.data.frame(niche.A)
write.xlsx(niche.A_df, "./species.niche.A.xlsx",sep = "", rowNames = TRUE)
otu <- t(dune10)
otu[otu > 0] <- 1
niche <- read.xlsx('species.niche.A.xlsx')
if (ncol(niche) > 0) {
  rownames(niche) <- niche[, 1]
  niche <- niche[, -1, drop = FALSE]
}
niche <- as.matrix(niche)
niche[is.na(niche)] <- 0
niche.matrix <- otu * niche[,1]
Bcom <- colSums(niche.matrix)/colSums(otu)
class(Bcom)
Bcom_df <- as.data.frame(Bcom)
write.xlsx(Bcom_df, "./Bcom.xlsx",sep = "", rowNames = TRUE)
spec_gen.all <- spec.gen(dune10, niche.width.method = 'levins',
                         perm.method = 'quasiswap', n = 100, 
                         probs = c(0.025, 0.975))
merge.table.all <- merge(spec_gen.all, genus.rela, by = 'row.names')
# 提取 genus.rela 的行名
genus_rela_rownames <- rownames(genus.rela)
str(merge.table.all$Row.names)
# 按照 genus.rela 的行名顺序对 merge.table.all 进行排序
# 由于 merge 函数使用 'row.names' 合并时，行名存储在名为 'Row.names' 的列中
merge.table.all <- merge.table.all[match(genus_rela_rownames, merge.table.all$Row.names), ]
rownames(merge.table.all) <- NULL
merge.table.all <- inner_join(merge.table.all, asv.taxonomy, by = c("Row.names" = "Species"))
merge.table.all <- merge.table.all %>% rename_with(~ "genus.2", matches("Row.names"))
write.xlsx(merge.table.all, "./NBMC.all.xlsx",sep = "", rowNames = FALSE)
