#rm(list = ls())
####获得当前路径#### 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#path <- getwd()
#path
####加载包####
set.seed(1515)
library(phyloseq)
library(ggClusterNet)
library(tidyverse)
library(Biostrings)
library(ggrepel)
library(tidyfst)
library(patchwork)
library(igraph)
library(ggsci)
library("pulsar")
library(ggrepel)
library(openxlsx)
library(stringr)
abasv <- read.delim('feature.species.2.txt',row.names = 1)
group <- read.delim('group.2.txt')
group <- group[#(group$type == "Water" | group$type == "WS") & 
                 group$Diets == "C" & group$Mineral == "Fe" & 
                 group$Group != 'WS' , ]
#group$level <- paste(group$Diets, group$level, sep = "")
row.names(group) <- group$Sample
selected_samples <- group$Sample
# Select only the columns in the ASV table that exist in selected_samples
first_column <- abasv[, 1, drop = FALSE]
asv <- abasv[, -1, drop = FALSE]
asv <- asv[, colnames(asv) %in% selected_samples, drop = FALSE]
asv <- asv[, match(selected_samples, colnames(asv))]
asv <- cbind(first_column, asv)
tax <- read.delim('../../1.ASV_table/tax/ASV_table.tax.txt', row.names = 1)
unique_tax <- tax[!duplicated(tax$Species), ]
asv.taxonomy <- unique_tax
asv.taxonomy <- asv.taxonomy[rownames(abasv), ]
write.table(asv.taxonomy, 'tax.species.merge.txt', quote = F, sep = '\t', row.names = F)
row.names(asv.taxonomy) <- asv.taxonomy$Species
asv$RowNames <- rownames(asv)
re.ASV <- inner_join(asv, asv.taxonomy, by = c("genus.2" = "Species"))
rownames(re.ASV) <- re.ASV$RowNames
re.ASV$RowNames <- NULL
ab.ASV <- re.ASV
#re.ASV.filter####
re.ASV$Genus <- NULL
re.ASV$Family <- NULL
re.ASV$Order <- NULL
re.ASV$Class <- NULL
re.ASV$Kingdom <- NULL
# Calculate the maximum value for each row from column 2 to the second last column
max_values <- apply(re.ASV[, 2:(ncol(re.ASV) - 1)], 1, max)

# Count the number of non - zero values for each row from column 2 to the second last column
non_zero_counts <- rowSums(re.ASV[, 2:(ncol(re.ASV) - 1)] != 0)

# Create a logical vector to keep track of rows that meet the previous conditions:
# 1. At least one non - zero value in columns 2 to the second last column
# 2. The maximum value in columns 2 to the second last column is greater than or equal to 0.0001
# 3. The number of non - zero values in columns 2 to the second last column is at least 10
keep_rows <- non_zero_counts > 0 & max_values >= 0.0001 & non_zero_counts >= 6

# Initialize a logical vector to record whether each row has at least two non - zero values 
# and mean >= 0.0001 in any group of 6 columns
has_two_or_more_nonzeros_and_mean <- rep(FALSE, nrow(re.ASV))

# Loop through each group of 6 columns from column 2 to the second last column
for(i in seq(2, ncol(re.ASV) - 1, by = 6)){
  # Determine the end column index of the current group, ensuring it does not exceed the second last column
  end_col <- min(i + 5, ncol(re.ASV) - 1)
  
  # Check if each row has at least three non - zero values in the current group of 6 columns
  non_zero_counts_in_group <- rowSums(re.ASV[, i:end_col] != 0)
  two_or_more_nonzeros_in_group <- non_zero_counts_in_group >= 3
  
  # Check if the number of non - zero values in the group is exactly 1 for any row
  single_non_zero_rows <- non_zero_counts_in_group == 1
  
  # If a row has exactly one non - zero value in the group, set all values in the group to 0
  if(any(single_non_zero_rows)){
    re.ASV[single_non_zero_rows, i:end_col] <- 0
  }
  
  # Calculate the mean value for each row in the current group of 6 columns
  group_means <- rowMeans(re.ASV[, i:end_col])
  
  # Check if the mean of each row in the current group is greater than or equal to 0.0001
  mean_condition <- group_means >= 0.0001
  
  # Combine the non - zero count condition and the mean condition
  combined_condition <- two_or_more_nonzeros_in_group & mean_condition
  
  # Update the has_two_or_more_nonzeros_and_mean vector. 
  # Use the logical OR operator to combine, so if a row meets the condition in any group, it will be marked as TRUE
  has_two_or_more_nonzeros_and_mean <- has_two_or_more_nonzeros_and_mean | combined_condition
}
# Combine the previous conditions and the new condition to update the keep_rows vector
keep_rows <- keep_rows & has_two_or_more_nonzeros_and_mean

# Subset the data frame based on the final logical vector
re.ASV <- re.ASV[keep_rows, ]
re.ASV.1 <- ab.ASV[rownames(re.ASV), colnames(re.ASV)]
# Print the processed data frame
#selected_cols <- re.ASV[, c(2:c(NCOL(re.ASV)-1))]
#colSums(selected_cols == 0)
#colSums(selected_cols != 0)
otutable <- re.ASV[, c(2:c(NCOL(re.ASV)-1))]
taxonomy <- re.ASV[, c(1, which(colnames(re.ASV) == "Phylum"))]
# Define the desired order for Phylum column in taxonomy
desired_phylum_order <- unique(taxonomy$Phylum)
# Convert the Phylum column to a factor with the specified order
taxonomy$Phylum <- factor(taxonomy$Phylum, levels = desired_phylum_order)
# Define the desired order for the Group column in group data frame
# Use the unique values in the Group column as the order
desired_group_order <- unique(group$Group)
# Convert the Group column in group data frame to a factor with the specified order
group$Group <- factor(group$Group, levels = desired_group_order)
# Ensure the row order of otutable matches the order in taxonomy
otutable <- otutable[match(taxonomy[, 1], rownames(otutable)), ]
pss <- phyloseq(sample_data(group),
               otu_table(as.matrix(otutable), taxa_are_rows=TRUE),
              tax_table(as.matrix(taxonomy)))
write.table(re.ASV, 'CFe.net.txt', quote = F, sep = '\t', row.names = F)

####path####
otupath = "./"
netpath = paste(otupath,"/network.CFe/",sep = "")
dir.create(netpath)
####网络分析主函数####
tab.r = network.pip(
    ps = pss,
    N = nrow(re.ASV),
    # ra = 0.05,
    big = TRUE,
    select_layout = FALSE,
    layout_net = "model_maptree2",
    r.threshold = 0.8,
    p.threshold = 0.05,
    maxnode = 2,
    method = "spearman",
    label = TRUE,
    lab = "elements",
    group = "Group",
    fill = "Phylum",
    size = "igraph.degree",
    zipi = TRUE,
    ram.net = TRUE,
    clu_method = "cluster_fast_greedy",
    step = 50,
    R = 10,
    ncpus = 8)
saveRDS(tab.r,paste0(netpath,"network.pip.sparcc.rds"))
#readRDS####
tab.r = readRDS(paste0(netpath,"network.pip.sparcc.rds"))
dat = tab.r[[2]]
cortab = dat$net.cor.matrix$cortab
# 大型相关矩阵跑出来不容易，建议保存，方便各种网络性质的计算
saveRDS(cortab,paste0(netpath,"cor.matrix.all.group.rds"))
cor = readRDS(paste0(netpath,"cor.matrix.all.group.rds"))
#-提取全部图片的存储对象
#pre-plot####
plot = tab.r[[1]]
# 提取网络图可视化结果
p0 = plot[[1]]
ggsave(paste0(netpath,"plot.network.2.pdf"),p0,width = 20,height = 9)
#ggsave(paste0(netpath,"plot.network2.pdf"),p0,width = 30,height = 28)
#zipi展示
plot[[2]]
#与随机网络的比对
plot[[3]]
#Net parameters####
#网络属性计算-丰富的网络属性，16个
i = 1
id = names(cor)
for (i in 1:length(id)) {
  igraph= cor[[id[i]]] %>% make_igraph()
  dat = net_properties.4(igraph,n.hub = F)
  head(dat,n = 16)
  colnames(dat) = id[i]
  if (i == 1) {
    dat2 = dat
  } else{
    dat2 = cbind(dat2,dat)}}
#head(dat2)
FileName <- paste(netpath, "net.network.attribute.data.xlsx", sep = "")
write.xlsx(dat2, FileName, rowNames = TRUE)
class(dat2)
#计算单个样本网络属性用于和其他指标关联
for (i in 1:length(id)) {
  pst = pss %>% subset_samples.wt("Group",id[i]) %>% remove.zero()
  dat.f = netproperties.sample(pst = pst,cor = cor[[id[i]]])
  # head(dat.f)
  if (i == 1) {
    dat.f2 = dat.f
  } else{
    dat.f2 = rbind(dat.f2,dat.f) }}
#head(dat.f2)
FileName <- paste(netpath,"net.network.attribute.data.sample.xlsx", sep = "")
write.xlsx(dat.f2,FileName, rowNames = TRUE)
map= sample_data(pss)
map$ID = row.names(map)
map = map %>% as_tibble()
dat3 = dat.f2 %>% rownames_to_column("ID") %>% inner_join(map,by = "ID")
FileName <- paste(netpath,"net.network.attribute.data.sample.add.group.info.xlsx", sep = "")
write.xlsx(dat3,FileName, rowNames = TRUE)
#节点属性计算
for (i in 1:length(id)) {
  igraph= cor[[id[i]]] %>% make_igraph()
  nodepro = node_properties(igraph) %>% as.data.frame()
  nodepro$Group = id[i]
  head(nodepro)
  colnames(nodepro) = paste0(colnames(nodepro),".",id[i])
  nodepro = nodepro %>%
    as.data.frame() %>%
    rownames_to_column("genus")
  # head(dat.f)
  if (i == 1) {
    nodepro2 = nodepro
  } else{
   nodepro2 = nodepro2 %>% full_join(nodepro,by = "genus")}}
#head(nodepro2)
FileName <- paste(netpath,"net.node.attribute.data.sample.xlsx", sep = "")
write.xlsx(nodepro2,FileName,rowNames = TRUE)
#可定制网络输出####
#####group 顺序####
desired_order <- c("Fe0", "SFe90", "OFe45" ,"OFe90")
#tab.r####
replace_rules <- c("nodes" = "N")
tab.r[[2]][["net.cor.matrix"]][["node"]][["label"]] <-  str_replace_all(tab.r[[2]][["net.cor.matrix"]][["node"]][["label"]], replace_rules)
tab.r[[2]][["net.cor.matrix"]][["edge"]][["label"]] <-  str_replace_all(tab.r[[2]][["net.cor.matrix"]][["edge"]][["label"]], replace_rules)
# 提取需要处理的 label 向量
label_vector <- tab.r[[2]][["net.cor.matrix"]][["node"]][["label"]]
# 遍历每个 label
for (i in seq_along(label_vector)) {
  label <- label_vector[i]
  # 提取 label 开头的分组信息
  group <- strsplit(label, ": ")[[1]][1]
  # 检查分组信息是否在 dat2 的列名中
  if (group %in% colnames(dat2)) {
    # 获取对应列的 L 值
    l_value <- dat2["num.edges(L)", group]
    # 将 links: 替换为 L: <对应的值>
    label <- gsub("links: ", paste0("L: ", l_value), label)
  }
  # 更新当前 label
  label_vector[i] <- label
}

# 将处理后的结果放回原来的位置
tab.r[[2]][["net.cor.matrix"]][["node"]][["label"]] <- label_vector
# 提取需要处理的 label 向量
label_vector <- tab.r[[2]][["net.cor.matrix"]][["edge"]][["label"]]
# 遍历每个 label
for (i in seq_along(label_vector)) {
  label <- label_vector[i]
  # 提取 label 开头的分组信息
  group <- strsplit(label, ": ")[[1]][1]
  # 检查分组信息是否在 dat2 的列名中
  if (group %in% colnames(dat2)) {
    # 获取对应列的 L 值
    l_value <- dat2["num.edges(L)", group]
    # 将 links: 替换为 L: <对应的值>
    label <- gsub("links: ", paste0("L: ", l_value), label)
  }
  # 更新当前 label
  label_vector[i] <- label
}
# 将处理后的结果放回原来的位置
tab.r[[2]][["net.cor.matrix"]][["edge"]][["label"]] <- label_vector
tab.r[[2]][["net.cor.matrix"]][["node"]][["label"]] <- factor(tab.r[[2]][["net.cor.matrix"]][["node"]][["label"]], levels = unique(tab.r[[2]][["net.cor.matrix"]][["node"]][["label"]]))
tab.r[[2]][["net.cor.matrix"]][["edge"]][["group"]] <- factor(tab.r[[2]][["net.cor.matrix"]][["edge"]][["group"]], levels = desired_order)
tab.r[[2]][["net.cor.matrix"]][["edge"]] <- tab.r[[2]][["net.cor.matrix"]][["edge"]][order(tab.r[[2]][["net.cor.matrix"]][["edge"]][["group"]]),]
rownames(tab.r[[2]][["net.cor.matrix"]][["edge"]]) <- NULL
tab.r[[2]][["net.cor.matrix"]][["edge"]][["label"]] <- factor(tab.r[[2]][["net.cor.matrix"]][["edge"]][["label"]], levels = unique(tab.r[[2]][["net.cor.matrix"]][["edge"]][["label"]]))
dat = tab.r[[2]]
dat$net.cor.matrix$node$Group <- factor(dat$net.cor.matrix$node$Group, levels = desired_order)
dat$net.cor.matrix$edge$Group <- factor(dat$net.cor.matrix$edge$Group, levels = desired_order)
dat$net.cor.matrix$node$group <- factor(dat$net.cor.matrix$node$group, levels = desired_order)
dat$net.cor.matrix$edge$group <- factor(dat$net.cor.matrix$edge$group, levels = desired_order)
node = dat$net.cor.matrix$node
node$sign <- nich.label$sign[match(node$genus.2, nich.label$genus.2)]
capitalize_first <- function(s) {
  if (nchar(s) > 0) {
    paste0(toupper(substr(s, 1, 1)), tolower(substr(s, 2, nchar(s))))
  } else {
    s
  }
}
node$sign <- sapply(node$sign, capitalize_first)
node$sign[node$sign == "Non significant"] <- "Neutral"
edge = dat$net.cor.matrix$edge
#node$Group <- factor(node$Group, levels = desired_order)
#edge$Group <- factor(edge$Group, levels = desired_order)
#node$group <- factor(node$group, levels = desired_order)
#edge$group <- factor(edge$group, levels = desired_order)
#node2$group <- factor(node2$group, levels = desired_order)
#node2$Group <- factor(node2$Group, levels = desired_order)
node2  = add.id.facet(node,"Group")
node2$Num <- as.numeric(row.names(node2))
#node2$label <- str_replace_all(node2$label, replace_rules)
node2$Group <- factor(node2$Group, levels = desired_order)
#node$label <- factor(node$label, levels = unique(node$label))
node2$label <- factor(node2$label, levels = unique(node2$label))
#head(edge)
#head(node)
#head(node2)
#关键hubs####
nich.label <- read.xlsx('NBMC.all.xlsx',rowNames = F)
nich.label$Species <- nich.label$genus.2
ab.ASV$sign <- nich.label$sign[match(ab.ASV$genus.2, nich.label$genus.2)]
#Fe0####
ps.C <- pss
ps.C@otu_table <- ps.C@otu_table[,c(1:6,25:30)]
ps.C@sam_data <- ps.C@sam_data[c(1:6,25:30),]
result <- corMicro(ps = ps.C,
                   N = 50,
                   method.scale = "TMM",
                   r.threshold=0.8,
                   p.threshold=0.05,
                   method = "spearman")
cor1 = result[[1]]
result4 = nodeEdge(cor = cor1)
#提取edge文件
edge.1 = result4[[1]]
#--提取节点文件
node.1 = result4[[2]]
dim(edge.1)
#edge.1$weight
igraph  = igraph::graph_from_data_frame(edge.1, directed = FALSE, vertices = node.1)
hub.C = hub_score(igraph)$vector %>%
  sort(decreasing = TRUE) %>%
  head(20) %>%
  as.data.frame()
colnames(hub.C) = "hub_sca"
ggplot(hub.C) +
  geom_bar(aes(x = hub_sca,y = reorder(row.names(hub.C),hub_sca)),
           stat = "identity",fill = "#4DAF4A")
common_rows <- intersect(row.names(hub.C), row.names(re.ASV.1))
merged_data <- data.frame(hub.C, re.ASV.1[common_rows, ])
columns_to_add <- c("Class", "Order", "Family", "Genus",  "Species", "sign")
for (col in columns_to_add) {
  merged_data[[col]] <- nich.label[[col]][match(merged_data$genus.2, nich.label$genus.2)]
}
merged_data$type <- c('Fe0')
rownames(merged_data) <- paste0('Fe0.', rownames(merged_data))
write.table(merged_data, file = 'Fe0.20.0.8.txt', quote = F, sep = '\t', row.names = F)
#SFe90####
ps.S <- pss
ps.S@otu_table <- ps.S@otu_table[, c(7:12,31:36)]
ps.S@sam_data <- ps.S@sam_data[c(7:12,31:36), ]
result.S <- corMicro(ps = ps.S,
                    N = 50,
                    method.scale = "TMM",
                    r.threshold=0.8,
                    p.threshold=0.05,
                    method = "spearman")
cor2 = result.S[[1]]
result5 = nodeEdge(cor = cor2)
#提取edge文件
edge.2 = result5[[1]]
#--提取节点文件
node.2 = result5[[2]]
igraph.2  = igraph::graph_from_data_frame(edge.2, directed = FALSE, vertices = node.2)
hub.S = hub_score(igraph.2)$vector %>%
  sort(decreasing = TRUE) %>%
  head(22) %>%
  as.data.frame()
colnames(hub.S) = "hub_sca"
ggplot(hub.S) +
  geom_bar(aes(x = hub_sca,y = reorder(row.names(hub.S),hub_sca)),
           stat = "identity",fill = "#4DAF4A")
common_rows2 <- intersect(row.names(hub.S), row.names(re.ASV.1))
merged_data2 <- data.frame(hub.S, re.ASV.1[common_rows2, ])
for (col in columns_to_add) {
  merged_data2[[col]] <- nich.label[[col]][match(merged_data2$genus.2, nich.label$genus.2)]
}
merged_data2$type <- c('SFe90')
rownames(merged_data2) <- paste0('SFe90.', rownames(merged_data2))
write.table(merged_data2, file = 'SFe90.22.0.8.txt', quote = F, sep = '\t', row.names = F)
#OFe45####
ps.O2 <- pss
ps.O2@otu_table <- ps.O2@otu_table[, c(13:18,37:42)]
ps.O2@sam_data <- ps.O2@sam_data[c(13:18,37:42), ]
result.O2 <- corMicro(ps = ps.O2,
                    N = 50,
                    method.scale = "TMM",
                    r.threshold=0.8,
                    p.threshold=0.05,
                    method = "spearman")
cor2 = result.O2[[1]]
result5 = nodeEdge(cor = cor2)
#提取edge文件
edge.2 = result5[[1]]
#--提取节点文件
node.2 = result5[[2]]
igraph.3  = igraph::graph_from_data_frame(edge.2, directed = FALSE, vertices = node.2)
hub.O2 = hub_score(igraph.3)$vector %>%
  sort(decreasing = TRUE) %>%
  head(22) %>%
  as.data.frame()
colnames(hub.O2) = "hub_sca"
ggplot(hub.O2) +
  geom_bar(aes(x = hub_sca,y = reorder(row.names(hub.O2),hub_sca)),
           stat = "identity",fill = "#4DAF4A")
common_rows2 <- intersect(row.names(hub.O2), row.names(re.ASV.1))
merged_data3 <- data.frame(hub.O2, re.ASV.1[common_rows2, ])
for (col in columns_to_add) {
  merged_data3[[col]] <- nich.label[[col]][match(merged_data3$genus.2, nich.label$genus.2)]
}
merged_data3$type <- c('OFe45')
rownames(merged_data3) <- paste0('OFe45.', rownames(merged_data3))
write.table(merged_data3, file = 'OFe45.22.0.8.txt', quote = F, sep = '\t', row.names = F)
#OFe90####
ps.O3 <- pss
ps.O3@otu_table <- ps.O3@otu_table[, c(19:24,43:48)]
ps.O3@sam_data <- ps.O3@sam_data[c(19:24,43:48), ]
result.O3 <- corMicro(ps = ps.O3,
                    N = 50,
                    method.scale = "TMM",
                    r.threshold=0.8,
                    p.threshold=0.05,
                    method = "spearman",
                    ncpus = 8)
cor2 = result.O3[[1]]
result5 = nodeEdge(cor = cor2)
#提取edge文件
edge.2 = result5[[1]]
#--提取节点文件
node.2 = result5[[2]]
igraph.O3  = igraph::graph_from_data_frame(edge.2, directed = FALSE, vertices = node.2)
hub.O3 = hub_score(igraph.O3)$vector %>%
  sort(decreasing = TRUE) %>%
  head(2) %>%
  as.data.frame()
colnames(hub.O3) = "hub_sca"
ggplot(hub.O3) +
  geom_bar(aes(x = hub_sca,y = reorder(row.names(hub.O3),hub_sca)),
           stat = "identity",fill = "#4DAF4A")
common_rows2 <- intersect(row.names(hub.O3), row.names(re.ASV.1))
merged_data4 <- data.frame(hub.O3, re.ASV.1[common_rows2, ])
for (col in columns_to_add) {
  merged_data4[[col]] <- nich.label[[col]][match(merged_data4$genus.2, nich.label$genus.2)]
}
merged_data4$type <- c('OFe90')
rownames(merged_data4) <- paste0('OFe90.', rownames(merged_data4))
write.table(merged_data4, file = 'OFe90.2.0.8.txt', quote = F, sep = '\t', row.names = F)
combined_data <- rbind(merged_data, merged_data2, merged_data3, merged_data4)
write.xlsx(combined_data, file = 'Hub.0.6.xlsx', quote = F, sep = '\t', rowNames = T)
#color####
color30 <- c('#d32c1f', '#acc2d9', "#00bFcF", '#56ae57', '#a8ff04','#b2996e', '#894585', 
             '#d4ffff', '#fcfc81', '#388004', '#efb435', '#ff08e8', 
              '#1f6357',  "#EEE8AA", '#ffd8b1','#ff0789', '#3778bf', '#430541','#a9a9a9',
             '#850e04', '#aaffc3', '#f97306', '#ffb2d0', '#ad900d', '#f6688e',
             '#76fda8',  '#41fdfe', '#0165fc', '#0c1793', '#a50055', '#ad03de','#aeff6e', 
             '#fffd01', '#0165fc')

####Network####
#hub10 <- read.delim('Hubstop.txt')
#merged_data2 <- merge(node2, hub10, by.x=c("elements", "Group"), by.y=c("genus", "type"), all=TRUE)
#merged_data2 <- merged_data2[order(merged_data2$Num),]
#row.names(merged_data2) <- NULL
#merged_data2$genus <- merged_data2$elements
#merged_data2$genus[is.na(merged_data2$label2)] <- ""
#huball <- merged_data2
#huball[is.na(huball)] <- ""
#huball.clean <- huball %>%
#  filter(!(rowSums(is.na(huball) | huball == "") > 0))
#write.table(huball, file = 'hubs.marker.txt', quote = F, sep = '\t', row.names = F)
#write.table(huball.clean, file = 'hubs.clean.txt', quote = F, sep = '\t', row.names = F)
#get_size <- function(x) {
#  ifelse(x < 50, 
#         scales::rescale(x, to = c(0.8, 2), from = c(0, 50)),
#         ifelse(x < 100, 
#                scales::rescale(x, to = c(2, 3.5), from = c(50, 100)),
#                ifelse(x < 150, 
#                       scales::rescale(x, to = c(3.5, 4.5), from = c(100, 150)),
#                       scales::rescale(x, to = c(4.5, 5), from = c(150, max(node$igraph.degree)))
#                )
#         )
#  )
#}
#node$size <- get_size(node$igraph.degree)
p <- ggplot() + 
  theme_bw(base_size = 20) +
  geom_segment(data = edge, aes(x = X1, y = Y1, xend = X2, yend = Y2,color = cor),
               linewidth = 0.3,alpha = 1) +
  geom_point(data = node, aes(X1, X2, fill = Phylum, size = igraph.degree),
             pch = 21, color = "gray40") +
  facet_wrap(.~ label, scales="free",nrow = 2) +
 # geom_text_repel(data = merged_data2, aes(X1, X2,label = label2),pch = 21) +
  #geom_text(data = merged_data2, aes(X1, X2,label = label2), size = 2) +
  scale_colour_manual(values = c("#6D98B5","#D48852")) +
  scale_fill_manual(values = color30) +
  scale_size(range = c(0.5, 4.5)) +
  #scale_size_identity() +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        legend.background = element_rect(colour = NA),
        panel.background = element_rect(fill =  NA),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())
p
ggsave(paste0(netpath,"plot.network.pdf"),p,width = 15.8, height = 8)
color3 <-  c("#0000dc", "#dddddd", "#dc0000")
p.nich <- ggplot() + 
  theme_bw(base_size = 20) +
  geom_segment(data = edge, aes(x = X1, y = Y1, xend = X2, yend = Y2,color = cor),
               linewidth = 0.3,alpha = 1) +
  geom_point(data = node, aes(X1, X2, fill = sign, size = igraph.degree),
             pch = 21, color = "gray40") +
  facet_wrap(.~ label, scales="free",nrow = 2) +
  # geom_text_repel(data = merged_data2, aes(X1, X2,label = label2),pch = 21) +
  #geom_text(data = merged_data2, aes(X1, X2,label = label2), size = 2) +
  scale_colour_manual(values = c("#6D98B5","#D48852")) +
  scale_fill_manual(values = color3) +
  scale_size(range = c(0.5, 4.5)) +
  #scale_size_identity() +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        legend.background = element_rect(colour = NA),
        panel.background = element_rect(fill =  NA),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())
p.nich
ggsave(paste0(netpath,"plot.network.nich.pdf"), p.nich, width = 10.5, height = 8)
write.xlsx(node,'node.xlsx',rowNames = F)
#-8.6 zipi可视化-定制####
dat.z = dat$zipi.data
#head(dat.z)
x1<- c(0, 0.62,0,0.62)
x2<- c( 0.62,1,0.62,1)
y1<- c(-Inf,2.5,2.5,-Inf)
y2 <- c(2.5,Inf,Inf,2.5)
lab <- c("peripheral",'Network hubs','Module hubs','Connectors')
roles.colors <- c("#E6E6FA","#DCDCDC","#F5FFFA", "#FAEBD7")
tab = data.frame(x1 = x1,y1 = y1,x2 = x2,y2 = y2,lab = lab)
tem = dat.z$group %>% unique() %>% length()
for ( i in 1:tem) {
  if (i == 1) {
    tab2 = tab
  } else{
    tab2 = rbind(tab2,tab)}}
dat.z$group <- factor(dat.z$group, levels = desired_order)
p.zipi <- ggplot() +
  geom_rect(data=tab2, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill = lab))+
  guides(fill=guide_legend(title="Topological roles")) +
  scale_fill_manual(values = roles.colors)+
  geom_point(data=dat.z,aes(x=p, y=z,color=module), size = 1.2) + 
  theme_bw()+
  guides(color= 'none') +
  ggrepel::geom_text_repel(data = dat.z, aes(x = p, y = z, color = module,label=label),
                           size = 1.2,, max.overlaps = 20, force = 1, force_pull = 0.1)+
  facet_wrap(.~group, scale='free', nrow = 2) +
  #facet_grid(.~ group, scale='free', ncol = 2) +
  theme(strip.background = element_rect(fill = "white"))+
  xlab("Participation Coefficient") +
  ylab("Within-module connectivity z-score")
p.zipi
ggsave(paste0(netpath,"plot.zipi.pdf"),p.zipi, width = 155,height = 115, units = 'mm')
dat.z.1 <- dat.z
dat.z.1 <- dat.z.1[dat.z.1$label != '',]
dat.z.1$rownames <- row.names(dat.z.1)
dat.z.1.1 <- inner_join(dat.z.1, re.ASV.1, by = c("label" = "genus.2"))
row.names(dat.z.1.1) <- dat.z.1.1$rownames
dat.z.1.1 <- dat.z.1.1[row.names(dat.z.1),]
columns_to_add <- c("Class", "Order", "Family", "Genus",  "Species", "sign")
for (col in columns_to_add) {
  dat.z.1.1[[col]] <- nich.label[[col]][match(dat.z.1.1$label, nich.label$genus.2)]
}
#write.table(dat.z.1.1, file = './network.CFe/plot.zipi.data.txt', quote = F, sep = '\t', row.names = T)
write.xlsx(dat.z.1.1, file = './network.CFe/plot.zipi.data.1.xlsx', quote = F, sep = '\t', rowNames = T)
# 8.7 随机网络，幂率分布#####
dat.r = dat$random.net.data
dat.r$g <- factor(dat.r$g, levels = desired_order)
p3 <- ggplot(dat.r) +
  geom_point(aes(x = ID,y = network,
                 group =group,fill = group),pch = 21,size = 1.5) +
  geom_smooth(aes(x = ID,y = network,group =group,color = group))+
  #facet_grid(.~g, scales = "free") +
  facet_wrap(.~g, scales = 'free', nrow = 2) +
  theme_bw() + 
  theme(plot.margin=unit(c(0,0,0,0), "cm"))
p3
ggsave(paste0(netpath,"plot.幂律分布.pdf"),p3, width = 155,height = 110, units = 'mm')
#多网络比对-网络显著性####
dat = module.compare.net.pip(
  ps = NULL,
  corg = cor,
  degree = TRUE,
  zipi = FALSE,
  r.threshold= 0.8,
  p.threshold=0.05,
  method = "spearman",
  padj = F,
  n = 8)
res = dat[[1]]
#head(res)
FileName <- paste(netpath,"net.compare.diff.sig.xlsx", sep = "")
write.xlsx(res,FileName,rowNames = TRUE)
#网络稳定性-模块比对####
res1 = module.compare.m(
  ps = pss,
  Top = 250,
  corg = cor,
  zipi = FALSE,
  zoom = 0.2,
  degree = TRUE,
  r.threshold= 0.8,
  p.threshold=0.05,
  method = "spearman",
  padj = F,
  n = 8)

#不同分组使用一个圆圈展示，圆圈内一个点代表一个模块，相连接的模块代表了相似的模块。
p.res1 = res1[[1]]
p.res1
ggsave(paste0(netpath,"plot.modules.pdf"),p.res1, width = 5,height = 4)
#--提取模块的OTU，分组等的对应信息
dat1 = res1[[2]]
#head(dat1)
#模块相似度结果表格
dat2 = res1[[3]]
#head(dat2)
dat2$m1 = dat2$module1 %>% strsplit("model") %>%
  sapply(`[`, 1)
dat2$m2 = dat2$module2 %>% strsplit("model") %>%
  sapply(`[`, 1)
dat2$cross = paste(dat2$m1, "\nVs\n", dat2$m2, sep = "")
# head(dat2)
dat2 = dat2 %>% filter(module1 != "none")
p2 = ggplot(dat2) + 
  geom_bar(aes(x = cross,fill = cross), width = 0.6) +
  labs(x = NULL,  y = "Numbers of similar modules")+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 15)) +
  geom_text(data = data.frame(x = 1, y = 9.5, label = "10"),
            mapping = aes(x = x, y = y, label = label), color = 'white',
           size = 5, lineheight = 1.45, fontface = 2, inherit.aes = FALSE) +
  geom_text(data = data.frame(x = 2, y = 13.5, label = "14"),
            mapping = aes(x = x, y = y, label = label), color = 'white',
            size = 5, lineheight = 1.45, fontface = 2, inherit.aes = FALSE) +
  geom_text(data = data.frame(x = 3, y = 3.5, label = "4"),
            mapping = aes(x = x, y = y, label = label), color = 'white',
            size = 5, lineheight = 1.45, fontface = 2, inherit.aes = FALSE) +
  geom_text(data = data.frame(x = 4, y = 2.5, label = "3"),
            mapping = aes(x = x, y = y, label = label), color = 'white',
            size = 5, lineheight = 1.45, fontface = 2, inherit.aes = FALSE) +
  geom_text(data = data.frame(x = 5, y = 5.5, label = "6"),
            mapping = aes(x = x, y = y, label = label), color = 'white',
            size = 5, lineheight = 1.45, fontface = 2, inherit.aes = FALSE) +
  geom_text(data = data.frame(x = 6, y = 0.5, label = "1"),
            mapping = aes(x = x, y = y, label = label), color = 'white',
            size = 5, lineheight = 1.45, fontface = 2, inherit.aes = FALSE) +
 ggprism::theme_prism() +
 theme(axis.text.x = element_text(angle = 0), legend.position = "none") 
p2
ggsave(paste0(netpath,"plot.similar modules.pdf"),p2, width = 4,height = 4)
#网络稳定性-鲁棒性####
res2= Robustness.Targeted.removal(ps = pss,
                                  corg = cor,
                                  degree = TRUE,
                                  Top = nrow(re.asv),
                                  r.threshold= 0.8,
                                  p.threshold=0.05,
                                  method = "spearman",
                                  zipi = FALSE)
library(gridExtra)
res2[[2]][["Group"]] <- factor(res2[[2]][["Group"]], levels = desired_order)
res2[[3]][["data"]][["Group"]] <- factor(res2[[3]][["data"]][["Group"]], levels = desired_order)
res2[[4]][["data"]][["Group"]] <- factor(res2[[4]][["data"]][["Group"]], levels = desired_order)
p3 = res2[[3]]
p3.2 <- res2[[4]]  +  theme_bw() 
p3 <- p3 +  theme_bw() 
#p3
p3.3 <- grid.arrange(p3, p3.2, ncol = 2)
#提取数据
dat4 = res2[[2]]
dir.create("./Robustness_Targeted_removal/")
path = paste(netpath,"/Robustness_Targeted_removal/",sep = "")
fs::dir_create(path)
write.xlsx(dat4, paste(path,"targeted_removal_network.xlsx",sep = ""), rowNames = TRUE)
ggsave(paste(path,"targeted_removal_network.pdf",sep = ""), p3.3,width = 8.5,height = 3.5)
#去除随机节点
res3 = Robustness.Random.removal(ps = pss,
                                 corg = cortab,
                                 r.threshold= 0.8,
                                 Top = 0.5)
res3[[2]][["Group"]] <- factor(res3[[2]][["Group"]], levels = desired_order)
res3[[3]][["data"]][["Group"]] <- factor(res3[[3]][["data"]][["Group"]], levels = desired_order)
res3[[4]][["data"]][["Group"]] <- factor(res3[[4]][["data"]][["Group"]], levels = desired_order)

p4 = res3[[3]] +  theme_bw()
p4.2 <- res3[[4]]  +  theme_bw() 
p4.3 <- grid.arrange(p4, p4.2, ncol = 2)
#提取数据
dat5 = res3[[2]]
# head(dat5)
path = paste(netpath,"/Robustness_Random_removal/",sep = "")
fs::dir_create(path)
write.xlsx(dat5,paste(path,"random_removal_network.xlsx",sep = "") , rowNames = TRUE)
ggsave(paste(path,"random_removal_network.pdf",sep = ""),  p4.3, width = 8.5, height = 3.5)
#网络稳定性-网络抗毁性#####
#ps <- pss
res6 = natural.con.microp (
  ps = pss,
  Top = 350,
  r.threshold= 0.8,
  p.threshold=0.05,
  method = "spearman",
  norm = F,
  end = 200,# 小于网络包含的节点数量
  start = 0,
  con.method = "pulsar"
)
p7 = res6[[1]]
#p7[["data"]] <- p7[["data"]][c(1:150,201:350,401:550,601:750), ]
p7[["data"]][["Group"]] <- factor(p7[["data"]][["Group"]], levels = desired_order)
p7 <- p7 + geom_point(alpha = 0.3, size = 1.5) +
 # scale_x_discrete(limits=c("C","S","SO","O1", "O2", "O3", "O4", "Feeds", "WR")) +
  #geom_smooth(method = 'loess', formula = y ~ x) +
  theme(axis.title = element_text(size = 13, face = "bold", color = 'black'), 
        axis.text = element_text(size = 10, face = "bold", color = 'black'), 
        plot.title = element_text(face = "bold"),
    legend.text = element_text(size = 10, face = "bold"), 
    legend.key = element_rect(fill = NA),
    legend.background = element_rect(fill = NA)) +
  labs(x = "Num of remove nodes", y = "Natural connectivity", colour = NULL)
p7
res6[[2]][["data"]][["Group"]] <- factor(res6[[2]][["data"]][["Group"]], levels = desired_order)
dat8  = res6[[2]]
path = paste(netpath,"/Natural_connectivity/",sep = "")
fs::dir_create(path)
write.xlsx(dat8,paste(path,"/Natural_connectivity.xlsx",sep = ""), rowNames = TRUE)
ggsave(paste(path,"/Natural_connectivity.pdf",sep = ""),  p7,width = 4.6, height = 3.75)
