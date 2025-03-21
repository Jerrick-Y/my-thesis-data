rm(list = ls())
####path#### 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#path <- getwd()
#path
####packages####
library(vegan)
library(picante)
####data####
asv <- read.delim('../1.ASV_table/ASV_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
asv$taxonomy <- NULL
asv <- t(asv)
tree <- read.tree("../1.ASV_table/fasta/ASV_table_tree.rooted.nwk")
tree$tip.label <- gsub("'", "", tree$tip.label)
tree_tips <- tree$tip.label
tree_tips
alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  ASVs <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Shan.divers <- exp(1)^Shannon
  Pielou <- Shannon / log(ASVs, base)
  Simpson <- diversity(x, index = 'simpson')
  sim.divers <- 1/(1 - Simpson)
  sim.even <- 1/(ASVs * (1 - Simpson))
  inv.sim <- simpson.unb(x, inverse = TRUE)
  fish.a <- fisher.alpha(x, MARGIN = 1)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  result <- data.frame( ASVs, Chao1, ACE, Shannon, Shan.divers, Pielou, 
                        Simpson, sim.divers, inv.sim, sim.even,
                        fish.a, goods_coverage)
  if (!is.null(tree)) {
    Faith_PD <- pd(x, tree, include.root = TRUE)[1]
    names(Faith_PD) <- 'Faith_PD'
    result <- cbind(result, Faith_PD)
  }
  result
}
alpha_all <- alpha(asv, tree)
alpha.2 <- cbind(samples = rownames(alpha_all), alpha_all)
alpha.3 <- alpha.2[, c(1,2,6,9,14,4,13)]
#rownames(alpha_all) <- NULL
write.table(alpha.2, file = "alpha.txt",sep = "\t", quote = FALSE, row.names = F)
write.table(alpha.3, file = "alpha.1.txt",sep = "\t", quote = FALSE, row.names = F)
####clean####
rm(list = ls())