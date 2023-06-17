early_TPO_scp1 <- tpo[, which(ridges$state=="early" & ridges$class=="TPO" & ridges$cluster=="5, nonmyelinating SCPs")]@assays$RNA@data
late_TPO_scp1 <- tpo[, which(ridges$state=="late" & ridges$class=="TPO" & ridges$cluster=="6, myelinating SCPs")]@assays$RNA@data
early_TPO_scp2 <- tpo[, which(ridges$state=="early" & ridges$class=="TPO" & ridges$cluster=="5, nonmyelinating SCPs")]@assays$RNA@data
late_TPO_scp2 <- tpo[, which(ridges$state=="late" & ridges$class=="TPO" & ridges$cluster=="6, myelinating SCPs")]@assays$RNA@data

profile1 <- data.frame(genes = rownames(early_TPO_scp1))
profile2 <- data.frame(genes = rownames(early_TPO_scp2))
rownames(profile1) <- profile1$genes
rownames(profile2) <- profile2$genes
profile1$early_mean <- apply(early_TPO_scp1, 1, mean)
profile1$late_mean <- apply(late_TPO_scp1, 1, mean)
profile2$early_mean <- apply(early_TPO_scp2, 1, mean)
profile2$late_mean <- apply(late_TPO_scp2, 1, mean)

profile1 <- profile1[-which(profile1$early_mean==0 | profile1$late_mean==0),]
profile2 <- profile1[-which(profile2$early_mean==0 | profile2$late_mean==0),]
profile1$FC <- profile1$late_mean / profile1$early_mean
profile2$FC <- profile2$late_mean / profile2$early_mean
profile1$logFC <- log2(profile1$FC)
profile2$logFC <- log2(profile2$FC)
profile1 <- profile1[order(profile1$logFC, decreasing = TRUE),]
profile2 <- profile2[order(profile2$logFC, decreasing = TRUE),]
profile1$pvalue <- apply(cbind(early_TPO_scp1[profile1$genes,], late_TPO_scp1[profile1$genes,]), 1,
                         function(x){
                             return(wilcox.test(x[1:17], x[18:73])[[3]])
                         })
profile2$pvalue <- apply(cbind(early_TPO_scp2[profile2$genes,], late_TPO_scp2[profile2$genes,]), 1,
                         function(x){
                             return(wilcox.test(x[1:17], x[18:73])[[3]])
                         })
profile1$fdr <- p.adjust(profile1$pvalue, method = "BH")
profile2$fdr <- p.adjust(profile2$pvalue, method = "BH")
write.table(profile1$genes[which(profile1$fdr<0.05&abs(profile1$logFC)>1)], file = "scp1_process_gene.txt", sep = "\n", row.names = F, col.names = F, quote = F)
write.table(profile2$genes[which(profile2$fdr<0.05&abs(profile2$logFC)>1)], file = "scp2_process_gene.txt", sep = "\n", row.names = F, col.names = F, quote = F)


library(clusterProfiler)
#data <- data.frame(SYMBOL = profile1$genes, logFC = profile1$logFC)
geneList1 <- profile1$logFC
names(geneList1) <- profile1$genes
geneList2 <- profile2$logFC
names(geneList2) <- profile2$genes

all <- read.gmt("mh.all.v2023.1.Mm.symbols.gmt")
g1 <- GSEA(geneList1, TERM2GENE = all)
g2 <- GSEA(geneList2, TERM2GENE = all)

dotplot(g1, x = "GeneRatio")
dotplot(g2, x = "GeneRatio")
