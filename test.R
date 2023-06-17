library('Seurat')
library('dplyr')
library('monocle')
library('scales')
load(file = "./Myelofibrosis scRNA-seq Supplementary Data /Robjects&Markdown/tpo.umap.seurat2.rds")
tpo <- UpdateSeuratObject(tpo)
pdf("umap.pdf")
DimPlot(tpo, reduction = "umap")
dev.off()
coun <- c()
for(i in 1:13834){
    coun <- c(coun, sum(tpo[i]@assays$RNA@counts!=0))
}
library(pheatmap)
showGenes <- c("Cd34", "Ly6a", "Ly6c1", "Sparcl1", "Eng", "Pdgfra",
               "Lepr", "Cxcl12", "Kitl", "Pdgfrb", "Sp7", "Bglap",
               "Bglap2", "Alpl", "Mog", "Mal", "Sox10", "Mobp")
avg_exp <- AverageExpression(tpo, verbose=T) #verbose是日志显示选项
avg_exp <- avg_exp$RNA
all_avg_exp <- avg_exp[showGenes,]
pdf("cluster_pheatmap.pdf")
pheatmap(all_avg_exp, cluster_rows=F, cluster_cols=T, scale="row", fontsize_col = 15, angle_col = 315)
dev.off()
MSC_showGeans <- c("Mgp","Egr1","Col3a1","Srgn","Ccl19","Tnc",
                   "Wif1","Col8a1","Mmp14","Igf1r","Gapdh",
                   "Aldoa","Neat1","Malat1","Cxcl10","Cxcl9",
                   "Irgm1","Isg15")
DoHeatmap(tpo[,which(tpo@meta.data$Clusters %in% c(1, 2, 3, 4))],
          features = MSC_showGeans,
          group.by = "anno", slot = "data", assay = 'RNA')
MSC_avg_exp <- avg_exp[MSC_showGeans, c(1, 2, 3, 4)]
pdf("MSC_pheatmap.pdf")
pheatmap(MSC_avg_exp, cluster_rows=F, cluster_cols=F, scale="row", fontsize_col = 15, angle_col = 315)
dev.off()
SCP_showGenes <- c("Klk6","S100b","Mdm2","Anxa5","Cdkn1c",
                   "Opalin","Ptgds","Il33","Cd9")
SCP_avg_exp <- avg_exp[SCP_showGenes, c(5, 6)]
pdf("SCP_pheatmap.pdf")
pheatmap(SCP_avg_exp, cluster_rows=F, cluster_cols=F, scale="row", fontsize_col = 15, angle_col = 315)
dev.off()

seurat_data <- readRDS("seurat_data.rds")
library(SingleR)
library(celldex)
nhd.se <- NovershternHematopoieticData()
nhd.se

expr <- GetAssayData(seurat_data, slot="data")
seurat_nhd <- SingleR(test = expr, ref = nhd.se,
                      labels = nhd.se$label.main)
ind <- which(!seurat_data$seurat_clusters %in% c(0, 2, 6, 7, 9, 11, 12, 13))
table(seurat_nhd$labels, seurat_data$seurat_clusters)

expr2 <- GetAssayData(tpo, slot="data")
seurat_nhd2 <- SingleR(test = expr2, ref = nhd.se,
                      labels = nhd.se$label.main)
ind <- which(!seurat_data$seurat_clusters %in% c(0, 2, 6, 7, 9, 11, 12, 13))
table(seurat_nhd2$labels, tpo@meta.data$Clusters)

library(GSVA)
geneList <- read.table("../data142240/down-target.txt")
geneList <- geneList[, 1]
geneList <- geneList[-1]

#安装homologene这个R包
install.packages('homologene')
#加载homologene这个R包
library(homologene)
#这里以小鼠的三个基因为例
#更多基因方法是一样的
#使用homologene函进行转换
#@genelist是要转换的基因列表
#@inTax是输入的基因列表所属的物种号，9606是人
#@outTax是要转换成的物种号，10090是小鼠
geneList <- homologene(geneList, inTax = 9606, outTax = 10090)[, 2]
geneList <- list("score" = c(NA, geneList))
g <- gsva(expr, geneList, method="ssgsea", kcdf="Gaussian", abs.ranking=T)
g2 <- gsva(expr2, geneList, method="ssgsea", kcdf="Gaussian", abs.ranking=T)
library(ggplot2)
library(ggridges)
ridges <- data.frame(score = g2[1,], cluster = tpo@active.ident, class = tpo$orig.ident, state = tpo@meta.data$state)
ggplot(ridges[which(ridges$state=="early"),], aes(x = score, y = cluster, fill = class)) +
    geom_density_ridges(alpha = 0.5) 

ggplot(ridges[which(ridges$state=="late"),], aes(x = score, y = cluster, fill = class)) +
    geom_density_ridges(alpha = 0.5)

pdf("ridges.pdf")
ggplot(ridges[which(ridges$class=="TPO"),], aes(x = score, y = cluster, fill = state)) +
    geom_density_ridges(alpha = 0.5)+
    theme(axis.text = element_text(size = 20))
dev.off()

ggplot(ridges[which(ridges$class=="EV"),], aes(x = score, y = cluster, fill = state)) +
    geom_density_ridges(alpha = 0.5)

early_test <- c()
for(cluster in levels(tpo@active.ident)){
    early_test <- c(early_test,
              wilcox.test(g2[which(tpo@active.ident==cluster & tpo@meta.data$orig.ident=="EV" & tpo@meta.data$state=="early")],
                          g2[which(tpo@active.ident==cluster & tpo@meta.data$orig.ident=="TPO" & tpo@meta.data$state=="early")])[[3]])
}
late_test <- c()
for(cluster in levels(tpo@active.ident)){
    late_test <- c(late_test,
                    wilcox.test(g2[which(tpo@active.ident==cluster & tpo@meta.data$orig.ident=="EV" & tpo@meta.data$state=="late")],
                                g2[which(tpo@active.ident==cluster & tpo@meta.data$orig.ident=="TPO" & tpo@meta.data$state=="late")])[[3]])
}
tpo_test <- c()
for(cluster in levels(tpo@active.ident)){
    tpo_test <- c(tpo_test,
                   wilcox.test(g2[which(tpo@active.ident==cluster & tpo@meta.data$orig.ident=="TPO" & tpo@meta.data$state=="early")],
                               g2[which(tpo@active.ident==cluster & tpo@meta.data$orig.ident=="TPO" & tpo@meta.data$state=="late")])[[3]])
}
EV_test <- c()
for(cluster in levels(tpo@active.ident)){
    EV_test <- c(EV_test,
                  wilcox.test(g2[which(tpo@active.ident==cluster & tpo@meta.data$orig.ident=="EV" & tpo@meta.data$state=="early")],
                              g2[which(tpo@active.ident==cluster & tpo@meta.data$orig.ident=="EV" & tpo@meta.data$state=="late")])[[3]])
}
names(early_test) <- levels(tpo@active.ident)
names(late_test) <- levels(tpo@active.ident)
names(tpo_test) <- levels(tpo@active.ident)
names(EV_test) <- levels(tpo@active.ident)
early_test[which(early_test < 0.05)]
late_test[which(late_test < 0.05)]
tpo_test[which(tpo_test < 0.05)]
EV_test[which(EV_test < 0.05)]

ridges$profile <- paste0(ridges$class, "_", ridges$state)
ridges$profile <- factor(ridges$profile,
                            levels = c("EV_early", "EV_late",
                                       "TPO_early", "TPO_late"))

cluster_count <- table(ridges[which(ridges$class=="TPO"), "state"], ridges[which(ridges$class=="TPO"), "cluster"])
prop1 <- prop.table(cluster_count, 1)
prop1 <- data.frame(prop1)
prop1$Freq <- prop1$Freq*100
colnames(prop1)[3] <- "percentage"
prop1$label <- paste0(sprintf("%.1f", prop1$percentage), "%")
#看细胞占比变化
library(ggalluvial)
prop1$Var1 <- factor(rep(c("pre-MF", "MF"), time = 8), levels = c("pre-MF", "MF"))
pdf("cell_proporation.pdf")
ggplot(prop1, aes(x = Var1, y = percentage, fill = Var2, alluvium = Var2, stratum = Var2))+
    geom_alluvium(aes(fill = Var2), alpha = .5,width = 0.5)+
    geom_stratum(aes(fill = Var2),width = 0.5)+
    scale_y_continuous(labels = scales::percent_format(scale = 1))+
    geom_text(aes(label=label), vjust=1.2, hjust = 0.5, size=3, position = "stack", color="black")+
    labs(x = "state", fill = "cluster")+
    theme(axis.text.x = element_text(size = 15), legend.text = element_text(size = 15))
dev.off()

fisher_test <- list()
for(i in 1:ncol(cluster_count)){
    fisher_test[[colnames(cluster_count)[i]]] <- fisher.test(matrix(c(cluster_count[1, i], sum(cluster_count[1, -i]),
                                             cluster_count[2, i], sum(cluster_count[2, -i])), nrow = 2))[[1]]
}


ridges$sample <- tpo@meta.data$group

ridges$state <- factor(ridges$state, levels = c("early", "late"))
pplist1 <- list()
for(cluster_ in levels(ridges$cluster)){
    clu_data <- ridges[which(ridges$cluster == cluster_ & ridges$class == "TPO"), c(1, 4)]
    p <- ggplot(clu_data, aes(x = state, y = score, fill = state))+
        geom_boxplot()+
        labs(title = cluster_)+
        theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 15), legend.text = element_text(size = 15))
    compare_means(score ~ state, data = clu_data, size = 2)
    my_comparisons <- list(c("early", "late"))
    pplist1[[cluster_]] <- p + stat_compare_means(comparison = my_comparisons)
}

pdf("SCP_boxplot.pdf")
plot_grid(pplist1[[5]], pplist1[[6]])
dev.off()

plot_grid(pplist1[[1]], pplist1[[2]])
plot_grid(pplist1[[3]], pplist1[[4]])
plot_grid(pplist1[[5]], pplist1[[6]])
plot_grid(pplist1[[7]], pplist1[[8]])


#看一下细胞占比的差异情况
sample <- c("TPO_early1", "TPO_early2", "TPO_late1", "TPO_late2")
group <- rep(c("TPO_early", "TPO_late"), each = 2)
samples <- data.frame(sample, group)
rownames(samples) <- samples$sample
prop <- prop.table(table(ridges$sample, ridges$cluster)[5:8,], 1)
prop <- data.frame(prop)
library(reshape2)
#长数据转为宽数据
prop <- dcast(prop, Var1~Var2, value.var = "Freq")
rownames(prop) <- prop[, 1]
prop <- prop[, -1]
prop$sample <- samples[rownames(prop), "sample"]
prop$group <- factor(samples[rownames(prop), "group"], levels = c("TPO_early", "TPO_late"))

pplist2 <- list()
groups <- levels(ridges$cluster)
library(dplyr)
library(ggpubr)
library(cowplot)

for(group_ in groups){
    cellprop_ <- prop %>% select(one_of(c("sample", "group", group_)))
    colnames(cellprop_) <- c("sample", "group", "percent")
    cellprop_$percent <- as.numeric(cellprop_$percent)
    cellprop_ <- cellprop_ %>%
                 group_by(group) %>%
                 mutate(upper =  quantile(percent, 0.75), 
                        lower = quantile(percent, 0.25),
                        mean = mean(percent),
                        median = median(percent))#上下分位数
    pp1 <- ggplot(cellprop_, aes(x = group, y = percent))+
        geom_jitter(shape = 21, aes(fill = group), width = 0.25)+
        stat_summary(fun = mean, geom = "point", color = "grey60")+
        theme_cowplot()+
        labs(title = group_, y = "percentage")+
        geom_errorbar(aes(ymin = lower, ymax = upper), col = "grey60", width = 1)
    labely <- max(cellprop_$percent)
    compare_means(percent ~ group, data = cellprop_)
    my_comparisons <- list(c("TPO_early", "TPO_late"))
    pp1 <- pp1 + stat_compare_means(comparisons = my_comparisons, size = 3, method = "t.test")
    pplist2[[group_]] <- pp1
}
plot_grid(pplist2[["1, adipogenic MSC"]],
          pplist2[["2, osteogenic MSC"]],
          pplist2[["3, transition MSC"]],
          pplist2[["4, interferon high MSCs"]],
          pplist2[["5, nonmyelinating SCPs"]],
          pplist2[["6, myelinating SCPs"]],
          pplist2[["7, OLC"]],
          pplist2[["8, endothelial cells"]])

plot_grid(pplist2[["5, nonmyelinating SCPs"]])
plot_grid(pplist2[["6, myelinating SCPs"]])
for(cluster_ in names(pplist2)){
    p <- pplist2[[cluster_]]
    print(p)
}

#cellchat
library(Seurat)
library(dplyr)
library(SeuratData)
library(patchwork) #最强大的拼图包
library(ggplot2)
library(CellChat)
library(ggalluvial)
library(svglite)
options(stringsAsFactors = F) #输入数据不自动转换成因子（防止数据格式错误）

cellchatDB <- CellChatDB.mouse
# 查看数据库具体信息
CellChatDB$interaction[1:4,1:4]
head(cellchatDB$cofactor)
head(cellchatDB$complex)
head(cellchatDB$geneInfo)

data.input <- tpo@assays$RNA@data
cell.use <- list()
for(profile_ in levels(ridges$profile)){
    cell.use[[profile_]] <- rownames(ridges)[ridges$profile == profile_]
}

#预处理用于细胞通信分析的表达数据
cellchat <- list()
for(profile_ in names(cell.use)){
    cellchat[[profile_]] <- createCellChat(object = data.input[, cell.use[[profile_]]],
                                           meta = ridges[cell.use[[profile_]],],
                                           group.by = "cluster")
    #设置cluster为默认细胞标识
    cellchat[[profile_]] <- setIdent(cellchat[[profile_]], ident.use = "cluster")
    cellchat[[profile_]]@DB <- cellchatDB
    # subset the expression data of signaling genes for saving computation cost
    cellchat[[profile_]] <- subsetData(cellchat[[profile_]])
    cellchat[[profile_]] <- identifyOverExpressedGenes(cellchat[[profile_]])
    cellchat[[profile_]] <- identifyOverExpressedInteractions(cellchat[[profile_]])
    cellchat[[profile_]] <- projectData(cellchat[[profile_]], PPI.mouse)
    #细胞通信网络的推断
    #计算通信概率并推断cellchat网络
    cellchat[[profile_]] <- computeCommunProb(cellchat[[profile_]], raw.use = TRUE)
    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    #cellchat[[profile_]] <- filterCommunication(cellchat[[profile_]], min.cells = 10)
    #在信号通路级别推断细胞-细胞通信
    cellchat[[profile_]] <- computeCommunProbPathway(cellchat[[profile_]])
    #计算整合的细胞通信网络
    cellchat[[profile_]] <- aggregateNet(cellchat[[profile_]])
}

#可视化
#使用圆图显示任意两个细胞组之间的相互作用次数或总交互强度
inds <- c(3, 4)
for(ind in inds){
    groupSize <- as.numeric(table(cellchat[[ind]]@idents))
    #netVisual_circle(cellchat[[ind]]@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    netVisual_circle(cellchat[[ind]]@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
}


#检查每个细胞组发送的信号
mat <- cellchat[[profile_]]@net$weight
par(mfrow = c(2, 4), xpd=TRUE)
for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#使用层次结构图、圆图或和弦图可视化每个信号通路
##层次结构图
pathways.show <- c("CXCL")
vertex.receiver = seq(1,4)
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver)
##圈图
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
##和弦图
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
##热图
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

# 计算每个配体受体对整体信号通路的贡献，并可视化由单个配体受体对调节的细胞通信
netAnalysis_contribution(cellchat, signaling = pathways.show)
##可视化由单个配体受体对调节的细胞-细胞通信
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
### show one ligand-receptor pair
LR.show <- pairLR.CXCL[1,] 
## Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
## Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

