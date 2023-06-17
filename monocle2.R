library(monocle)
expr_matrix <- as(as.matrix(tpo@assays$RNA@counts), "sparseMatrix")
p_data <- tpo@meta.data
f_data <- data.frame(gene_short_name = rownames(tpo), row.names = rownames(tpo))
pd <- new("AnnotatedDataFrame", p_data)
fd <- new("AnnotatedDataFrame", f_data)
#构建CellDataSet对象
cds <- newCellDataSet(
    expr_matrix,
    phenoData = pd,
    featureData = fd,
    lowerDetectionLimit = 0.5,
    expressionFamily = negbinomial.size()
)
#估计size factor和离散度
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
#过滤低质量的细胞
#seurat对象过滤过了
#过滤基因
cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- rownames(subset(fData(cds), num_cells_expressed >= 10))
#轨迹定义基因选择及可视化和构建轨迹
#选择定义过程的基因
#使用monocle选择的高变基因
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 &
                         dispersion_empirical >= 1*dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)
#使用dpFeature选择排序基因
diff <- differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~anno", cores = 1)
head(diff)
deg <- subset(diff, qval < 0.01)
deg <- deg[order(deg$qval, decreasing = F),]
head(deg)
#轨迹构建基因可视化
ordergene <- rownames(deg)[1:2000]
cds <- setOrderingFilter(cds, ordergene)
plot_ordering_genes(cds)

#降维
cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")
#拟时间轴轨迹构建和在拟时间内排列细胞
#修改R包源码,把第35行的if语句删掉
trace('project2MST', edit = T, where = asNamespace("monocle"))
cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by = "Pseudotime", show_backbone = TRUE)
plot_cell_trajectory(cds, color_by = "anno", show_backbone = TRUE)

