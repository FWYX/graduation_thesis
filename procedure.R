data <- read.table("./GSE142240-GPL16384_series_matrix.txt",
                   comment.char="!", header=T)
library(GEOquery)
GPL <- getGEO("GPL16384", destdir=".")
ref <- Table(GPL)
ref <- ref[which(ref$`Species Scientific Name`=="Homo sapiens"),]
data <- data[data$ID_REF%in%ref$ID,]
rownames(data) <- data$ID_REF
data <- data[,-1]
weak_effect <- 1:7
strong_effect <- 8:14
group_list <- factor(rep(c("weak", "strong"), each = 7), levels=c("weak", "strong"))
colnames(data) <- paste(group_list, 1:ncol(data), sep="")
library(tidyr)
frame <- gather(data=data, key="sample_name", value="expr")
frame$group <- as.factor(rep(c("weak", "strong"), each = 7*nrow(data)))
#批次效应
library(ggplot2)
#boxplot图
ggplot(data=frame, aes(x=sample_name, y=expr, color=group))+
    geom_boxplot()+
    theme(axis.text.x=element_text(angle=30))
#热图
library(pheatmap)
pheatmap(data, cluster_row = FALSE, show_rownames = FALSE,
         scale = "row")

#片间标准化
BiocManager::install("preprocessCore", configure.args="--disable-threading")
library(preprocessCore)
quan_norm <- normalize.quantiles(as.matrix(data))
colnames(quan_norm) <- colnames(data)
rownames(quan_norm) <- rownames(data)
write.csv(quan_norm, "quan_norm.csv")
#boxplot图
library(tidyr)
frame2 <- gather(data=data.frame(quan_norm), key="sample_name", value="expr")
frame2$group <- as.factor(rep(c("weak", "strong"), each = 7*nrow(data)))
library(ggplot2)
ggplot(data=frame2, aes(x=sample_name, y=expr, color=group))+
    geom_boxplot()+
    theme(axis.text.x=element_text(angle=30))
#热图
library(pheatmap)
pheatmap(quan_norm, cluster_row = FALSE, show_rownames = FALSE,
         scale = "row")
#差异分析
res <- data.frame(row.names=rownames(quan_norm))
res$mean_weak <- apply(as.matrix(quan_norm[, 1:7]), 1, mean)
res$mean_strong <- apply(as.matrix(quan_norm[, 8:14]), 1, mean)
res$log2_fold <- apply(res[, 1:2], 1,
                       function(x){return(log2(x[1]/x[2]))})
res$p.value <- apply(quan_norm, 1,
                    function(x){return(t.test(x[1:7], x[8:14])[[3]])})
res$change <- apply(res, 1,
                    function(x){
                        if((as.numeric(x[4])>0.05)|(abs(as.numeric(x[3]))<0.5))
                            return("stable")
                        else{
                            if(as.numeric(x[3])>=0.5)
                                return("up")
                            if(as.numeric(x[3])<=(-0.5))
                                return("down")
                        }
                    })
res$change <- factor(res$change, levels=c("stable", "down", "up"))
ggplot(res, aes(x=log2_fold, y=-log10(p.value), color=change))+
    geom_point(alpha=0.4)+
    scale_colour_manual(values=c("grey", "blue", "red"))+
    geom_vline(xintercept=c(-0.5, 0.5), lty=2, col="black", lwd=0.8)+
    geom_hline(yintercept=-log10(0.05), lty=2, col="black", lwd=0.8)+
    labs(x="log2(fold change)", y="-log10(p-value)")

diff <- res[which(res$change!="stable"),]

#绘制差异基因热图
library(pheatmap)
type <- group_list
names(type) <- colnames(ori_data)
type <- as.data.frame(type)
dif_data <- quan_norm[rownames(diff),]
pheatmap(dif_data, annotation=type, show_rownames=F, cluster_cols=F)
p <- pheatmap(dif_data,
         annotation=type,
         cluster_cols = T,
         color = colorRampPalette(c("green", "white", "red"))(50),
         show_colnames = T,
         show_rownames = F,
         cluster_rows = T,
         cutree_rows = 3,
         scale="row",  #矫正
         #border_color ="NA",
         fontsize = 12,
         fontsize_row=10,
         fontsize_col=10)
#提取三个miRNA集
row_cluster <- cutree(p$tree_row, k = 3)
miRNA_set <- list()
for(i in 1:3){
    miRNA_set[[i]] <- names(row_cluster)[which(row_cluster==i)]
}

#miRNA集热图
pheatmap(dif_data[miRNA_set[[1]],], annotation=type,
         show_rownames=F, cluster_cols=F, scale="row",
         color = colorRampPalette(c("green", "white", "red"))(50))
pheatmap(dif_data[miRNA_set[[2]],], annotation=type,
         show_rownames=F, cluster_cols=F, scale="row",
         color = colorRampPalette(c("green", "white", "red"))(50))
pheatmap(dif_data[miRNA_set[[3]],], annotation=type,
         show_rownames=F, cluster_cols=F, scale="row",
         color = colorRampPalette(c("green", "white", "red"))(50))
names(miRNA_set) <- c("up2", "down", "up1")
#svm-rfe
library(e1071)
library(caret)
library(kernlab)
set.seed(123)
Profile=rfe(x=t(dif_data),
            y=as.numeric(group_list),
            sizes = c(2,4,6,8, seq(10,40,by=3)),
            rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
            methods="svmRadial")
par(las=1)
x = Profile$results$Variables
y = Profile$results$RMSE
plot(x, y, xlab="Variables", ylab="RMSE (Cross-Validation)", col="darkgreen")
lines(x, y, col="darkgreen")
wmin <- which.min(y)
wmin.x=x[wmin]
wmin.y=y[wmin]
points(wmin.x, wmin.y, col="blue", pch=16)
text(wmin.x, wmin.y, paste0('N=',wmin.x), pos=2, col=2)
svm_rfe_var <- Profile$optVariables

#LASSO
library(glmnet)
library(foreign)
f1 <- glmnet(x=t(dif_data), y=as.numeric(group_list), family="binomial", alpha=1, nlambda=100)
cvfit <- cv.glmnet(x=t(dif_data), y=as.numeric(group_list), family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
plot(cvfit)
coef=coef(f1, s = cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lasso_var=lassoGene[-1]

#WGCNA

write.table(ref[which(ref$ID%in%miRNA_set[[1]]), c(1, 10, 9)], "miRNA-set1.txt", sep="\t", row.names=FALSE)
write.table(ref[which(ref$ID%in%miRNA_set[[2]]), c(1, 10, 9)], "miRNA-set2.txt", sep="\t", row.names=FALSE)
write.table(ref[which(ref$ID%in%miRNA_set[[3]]), c(1, 10, 9)], "miRNA-set3.txt", sep="\t", row.names=FALSE)


#limma
library(limma)
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(quan_norm)
contrast.matrix <- makeContrasts("weak-strong", levels=design)
##step1
fit <- lmFit(quan_norm, design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)  
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)
res <- nrDEG
res$change <- apply(res, 1,
                    function(x){
                        if((as.numeric(x[5])>0.1)|(abs(as.numeric(x[1]))<0.75))
                            return("stable")
                        else{
                            if(as.numeric(x[3])>=0.75)
                                return("up")
                            if(as.numeric(x[3])<=(-0.75))
                                return("down")
                        }
                    })
res$change <- factor(res$change, levels=c("stable", "down", "up"))
ggplot(res, aes(x=logFC, y=-log10(adj.P.Val), color=change))+
    geom_point(alpha=0.4)+
    scale_colour_manual(values=c("grey", "blue", "red"))+
    geom_vline(xintercept=c(-0.75, 0.75), lty=2, col="black", lwd=0.8)+
    geom_hline(yintercept=-log10(0.1), lty=2, col="black", lwd=0.8)+
    labs(x="log2(fold change)", y="-log10(p-adj.value)")

diff <- nrDEG[which(res$change!="stable"),]

#绘制差异基因热图
library(pheatmap)
dif_data <- quan_norm[rownames(diff),]
colnames(dif_data) <- c(paste0("resis", 1:7), paste0("sens", 1:7))
type <- factor(c(rep("resistance", time = 7), rep("sensitive", time = 7)))
names(type) <- colnames(dif_data)
type <- as.data.frame(type)
pheatmap(dif_data, annotation=type, show_rownames=F, cluster_cols=F)
p <- pheatmap(dif_data,
              annotation=type,
              cluster_cols = T,
              color = colorRampPalette(c("green", "white", "red"))(50),
              show_colnames = T,
              show_rownames = F,
              cluster_rows = T,
              cutree_rows = 3,
              scale="row",  #矫正
              #border_color ="NA",
              fontsize = 12,
              fontsize_row=10,
              fontsize_col=10)
#提取三个miRNA集
row_cluster <- cutree(p$tree_row, k = 3)
miRNA_set <- list()
for(i in 1:3){
    miRNA_set[[i]] <- names(row_cluster)[which(row_cluster==i)]
}

#miRNA集热图
pheatmap(dif_data[miRNA_set[[1]],], annotation=type,
         show_rownames=F, cluster_cols=F, scale="row",
         color = colorRampPalette(c("green", "white", "red"))(50))
pheatmap(dif_data[miRNA_set[[2]],], annotation=type,
         show_rownames=F, cluster_cols=F, scale="row",
         color = colorRampPalette(c("green", "white", "red"))(50))
pheatmap(dif_data[miRNA_set[[3]],], annotation=type,
         show_rownames=F, cluster_cols=F, scale="row",
         color = colorRampPalette(c("green", "white", "red"))(50))
names(miRNA_set) <- c("down", "up1", "up2")
