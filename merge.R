filenames <- c("GSE107844.txt", "GSE153434_all.counts.txt", "GSE52093.txt", "GSE147026_mRNA-ADVSCK.All.anno.txt")
data <- list()
for(i in 1:4){
  print(i)
  data[[i]] <- read.table(filenames[i], header=T, sep="\t", quote="")
}

#保留data[[2]]的有用信息
data[[2]] <- data[[2]][,c(1, 3:22)]
#数据过滤
for(i in 1:2){
  data[[i]] <- data[[i]][apply(data[[i]][,-1], 1, function(x)return(sum(x)>0)),]
}
normal_sample <- list()
normal_sample[[1]] <- c(4:6)
normal_sample[[2]] <- c(1:10)
case_sample <- list()
case_sample[[1]] <- c(1:3)
case_sample[[2]] <- c(11:20)
#共有基因
togene<-intersect(data[[1]][,1], data[[2]][,1])

for(i in 1:2){
  print(i)
  #重复基因取均值
  dat <- matrix(nrow=length(togene), ncol=(ncol(data[[i]])-1))
  colnames(dat) <- colnames(data[[i]])[-1]
  rownames(dat) <- togene
  for(gene in togene){
    dat[gene,] <- apply(data[[i]][data[[i]][,1]==gene, -1], 2, mean)
  }
  
  if(i == 1){
    normal_merge <- dat[togene,normal_sample[[i]]]
    case_merge <- dat[togene,case_sample[[i]]]
  }
  else{
    normal_merge <- cbind(normal_merge, dat[togene,normal_sample[[i]]])
    case_merge <- cbind(case_merge, dat[togene,case_sample[[i]]])
  }
}

colnames(normal_merge) <- paste0("normal", 1:ncol(normal_merge))
colnames(case_merge) <- paste0("case", 1:ncol(case_merge))
data_merge <- as.data.frame(cbind(normal_merge, case_merge))
write.csv(data_merge, "data_merge.csv")
#去批次效应
library(sva)
group_list <- read.table("grouplist.txt", header=T, row.names = 1)
modcombat = model.matrix(~1, data = group_list)
batch = group_list$batch
combat_edata = ComBat(dat=data_merge, batch=batch, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
write.csv(combat_edata, "combat_edata.csv")

library(preprocessCore)
quan_norm <- normalize.quantiles(as.matrix(data_merge))
colnames(quan_norm) <- colnames(data_merge)
rownames(quan_norm) <- rownames(data_merge)
write.csv(quan_norm, "quan_norm.csv")
#找差异基因
merge_data <- read.csv("quan_norm.csv", header=T, row.names=1)
general_data <- data.frame(gene=rownames(merge_data))
general_data$normal_mean <- apply(merge_data, 1, function(x){return(mean(x[1:13]))})
general_data$case_mean <- apply(merge_data, 1, function(x){return(mean(x[14:26]))})
general_data$FC <- apply(general_data, 1, function(x){return(as.numeric(x[3])/as.numeric(x[2]))})
general_data$pValue <- apply(merge_data, 1,
                             function(x){return(t.test(x[1:13], x[14:26])[[3]])})
general_data$fdr <- p.adjust(general_data$pValue, method="BH")
general_data$logFC <- log2(general_data$FC)
deg <- general_data$gene[(general_data$FC>2|general_data$FC<0.5)&general_data$pValue<0.05]

deg_data <- merge_data[deg,]
#绘制差异基因热图
library(pheatmap)
type <- c(rep("N", 13), rep("C", 13))
names(type) <- colnames(deg_data)
type <- as.data.frame(type)
pheatmap(deg_data, annotation=type, show_rownames=F, cluster_cols=F)
pheatmap(deg_data,
         annotation=type,
         cluster_cols = T,
         color = colorRampPalette(c("green", "white", "red"))(50),
         show_colnames = T,
         show_rownames = F,
         cluster_rows = T,
         #cutree_rows = 2,
         scale="row",  #矫正
         #border_color ="NA",
         fontsize = 12,
         fontsize_row=10,
         fontsize_col=10)
#GSEA
library(clusterProfiler)
library(limma)
library(org.Hs.eg.db)
library(enrichplot)
gmt <- read.gmt("KEGG(1).gmt")
logFC <- general_data$logFC
names(logFC) <- general_data$gene
logFC <- sort(logFC, decreasing = T)
#genes <- names(logFC)
#genes <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#logFC <- logFC[genes$SYMBOL]
#names(logFC) <- genes$ENTREZID
kk <- GSEA(logFC, TERM2GENE = gmt, pvalueCutoff = 1)
kkTab <- as.data.frame(kk)
kkTab <- kkTab[kkTab$pvalue<0.05,]

#显示高富集的通路
termNum=5
kkUP <- kkTab[kkTab$NES>0,]
#kkUP <- kkUP[sort(kkUP$enrichmentScore, index.return=T)$ix,]
if(nrow(kkUP) >= termNum){
  showTerm <- row.names(kkUP)[1:termNum]
  gseaplot <- gseaplot2(kk, showTerm, base_size = 8, title="Enriched in Case")
  #pdf(file = "GSEA.highRisk1.pdf", width=7, height=5.5)
  print(gseaplot)
  #dev.off()
}else{
  showTerm <- row.names(kkUP)
  gseaplot <- gseaplot2(kk, showTerm, base_size = 8, title="Enriched in Case")
  #pdf(file = "GSEA.highRisk1.pdf", width=7, height=5.5)
  print(gseaplot)
  #dev.off()
}

kkDOWN <- kkTab[kkTab$NES<0,]
#kkDOWN <- kkDOWN[sort(kkDOWN$enrichmentScore, index.return=T)$ix,]
if(nrow(kkDOWN) >= termNum){
  showTerm <- row.names(kkDOWN)[1:termNum]
  gseaplot <- gseaplot2(kk, showTerm, base_size = 8, title="Reduced in Case")
  #pdf(file = "GSEA.highRisk2.pdf", width=7, height=5.5)
  print(gseaplot)
  #dev.off()
}else{
  showTerm <- row.names(kkDOWN)
  gseaplot <- gseaplot2(kk, showTerm, base_size = 8, title="Enriched in Case")
  #pdf(file = "GSEA.highRisk2.pdf", width=7, height=5.5)
  print(gseaplot)
  #dev.off()
}
#取cycle通路的基因
genes <- gmt[gmt$term=="KEGG_CELL_CYCLE", 2]
data <- as.matrix(read.csv("quan_norm.csv", header=T, row.names=1))
cycle_data <- data[genes,]

#构建随机森林
#找差异基因
library(randomForest)
library(ROCR)
library(genefilter)
library(Hmisc)
cyc_general_data <- data.frame(gene=rownames(cycle_data))
cyc_general_data$normal_mean <- apply(cycle_data, 1, function(x){return(mean(x[1:13]))})
cyc_general_data$case_mean <- apply(cycle_data, 1, function(x){return(mean(x[14:26]))})
cyc_general_data$FC <- apply(cyc_general_data, 1, function(x){return(as.numeric(x[3])/as.numeric(x[2]))})
cyc_general_data$pValue <- apply(cycle_data, 1,
                             function(x){return(t.test(x[1:13], x[14:26])[[3]])})
cyc_general_data$fdr <- p.adjust(cyc_general_data$pValue, method="BH")
cyc_general_data$logFC <- log2(cyc_general_data$FC)
deg <- cyc_general_data$gene[cyc_general_data$FC>2&cyc_general_data$pValue<0.05]
write.csv(cyc_general_data[general_data$gene%in%deg,], "cyc_deg_outlook.csv")
cyc_deg_data <- cycle_data[deg,]
scale_cyc_deg_data <- t(apply(cyc_deg_data, 1, scale))
colnames(scale_cyc_deg_data) <- colnames(cycle_data)
write.csv(cyc_general_data, "cycle_general_data.csv")
pheatmap(cyc_deg_data,
         annotation=type,
         cluster_cols = F,
         color = colorRampPalette(c("green", "white", "red"))(50),
         show_colnames = T,
         show_rownames = T,
         cluster_rows = T,
         #cutree_rows = 2,
         scale="row",  #矫正
         #border_color ="NA",
         fontsize = 12,
         fontsize_row=10,
         fontsize_col=10)

data <- as.data.frame(t(cyc_deg_data))
data$class <- as.factor(c(rep("N", 13), rep("C", 13)))
set.seed(123)
rf <- randomForest(class~., data=data, importance=T, proximity=T, ntree=500)
plot(rf)
important_f <- importance(rf)
varImpPlot(x=rf,sort=TRUE,n.var=nrow(rf$importance), main="Random Forest Variance Importance")
#选最重要的前25%个基因
rf_var <- names(sort(important_f[,3], decreasing=T)[1:floor(22/4)])
#svm-ref
library(e1071)
library(caret)
library(kernlab)
data <- cyc_deg_data
class <- as.factor(c(rep("N", 13), rep("C", 13)))
set.seed(123)
Profile=rfe(x=t(data),
            y=as.numeric(class),
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
#lasso
library(glmnet)
library(foreign)
data <- cyc_deg_data
class <- as.factor(c(rep("N", 13), rep("C", 13)))
f1 <- glmnet(x=t(data), y=as.numeric(class), family="binomial", alpha=1, nlambda=100)
cvfit <- cv.glmnet(x=t(data), y=as.numeric(class), family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
plot(cvfit)
coef=coef(f1, s = cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lasso_var=lassoGene[-1]

#遗传算法
library(mcga)
library(e1071)
library(caret)
data <- t(cyc_deg_data)
myfun <- function(select, data1=data){
  if(sum(select) == 0){
    return(0)
  }
  else{
    data2 <- data1[,select*(1:22)]
    data2 <- as.data.frame(data2)
    data2$class <- as.factor(c(rep("N", 13), rep("C", 13)))
    svm_model <- svm(class~., data=data1, type="C", kernel="radial")
    svm_pred <- predict(svm_model, newdata=data2[, -ncol(data)])
    acc <- confusionMatrix(data1$class, svm_pred)$overall[1]
    return((acc*10000)/sum(select+10000))
  }
}

#k倍交叉证实的准确率作为适应度函数
myfun <- function(select, data1=data){
  #print(select)
  if(sum(select) == 0){
    return(0)
  }
  else{
    data2 <- data1[,select*(1:22)]
    data2 <- as.data.frame(data2)
    if(sum(select)==1){
      colnames(data2) <- colnames(data1)[select*(1:22)]
    }
    data2$class <- as.factor(c(rep("N", 13), rep("C", 13)))
    #设置折叠数
    k <- 3
    #用列表储存分层折叠后的数据集
    dataSet <- list()
    #分层折叠时以fold_index为参考索引
    fold_index <- matrix(0, nrow=2, ncol=(k+1))
    base_num <- 13%/%k
    mod_num <- 13%%k
    fold_index[c(1,2), 1] <- 1
    for(i in 1:2){
      for(j in 2:(k+1)){
        if(j <= mod_num){
          fold_index[i, j] <- fold_index[i, j-1] + base_num + 1
        }
        else{
          fold_index[i, j] <- fold_index[i, j-1] + base_num
        }
      }
    }
    #分层折叠
    #将混乱的索引按类储存进列表
    index <- list()
    classes <- c("N", "C")
    for(i in 1:2){
      index[[i]] <- sample(which(data2$class==classes[i]),
                           length(which(data2$class==classes[i])),
                           replace=FALSE)
    }
    #将分层折叠的数据存入dataSet
    for(i in 1:k){
      dataSet[[i]] <- data2[1,]
      for(j in 1:2){
        dataSet[[i]] <- rbind(dataSet[[i]],
                              data2[index[[j]][fold_index[j,i]:(fold_index[j,i+1]-1)],])
      }
      dataSet[[i]] <- dataSet[[i]][-1, ]
    }
    #两个做训练集，一个做测试集
    acc <- c()
    for(i in 1:3){
      test_set <- dataSet[[i]]
      test_class <- test_set$class
      test_set <- as.data.frame(test_set[,-ncol(test_set)])
      if(sum(select)==1){
        colnames(test_set) <- colnames(data1)[select*(1:22)]
      }
      left <- dataSet[-i]
      train_set <- rbind(left[[1]], left[[2]])
      svm_model <- svm(class~., data=train_set, type="C", kernel="radial")
      svm_pred <- predict(svm_model, newdata=test_set)
      t_acc <- confusionMatrix(test_class, svm_pred)$overall[1]
      acc <- c(acc, t_acc)
    }
    return(mean(acc)*100/(sum(select)+100))
  }
}

m <- ga(fitness = myfun,#适应度函数
        popSize = 100,#每一代100个个体
        type = "binary",#数据类型
        names = colnames(data),
        crossover = gabin_uCrossover,
        nBits = 22,#变量的个数
        run = 25,#迭代50次，如果最大值仍没有改进则停止
        maxiter = 50,#最多迭代多少代
        pmutation = 0.5,#变异的概率，这里是25%
        monitor = plot,#每次迭代时绘制结果
        seed = 54321
        )
result <- summary(m)
result$solution
ga_var <- list()
for(i in 1:nrow(result$solution)){
  ga_var[[i]] <- colnames(result$solution)[result$solution[i,]==1]
}
#哑铃图
library(UpSetR)
library(ggplot2)
library(grid)
library(plyr)
listInput <- list(rf_opt=rf_var, svm_rfe_opt=svm_rfe_var,
                  lasso_opt=lasso_var,
                  ga_opt1=ga_var[[1]],
                  ga_opt2=ga_var[[2]],
                  ga_opt3=ga_var[[3]],
                  ga_opt4=ga_var[[4]],
                  ga_opt5=ga_var[[5]])
upset(fromList(listInput), order.by="freq")
feature_sum <- data.frame(gene=Reduce(union, listInput))
features_count <- c()
for(gene in feature_sum$gene){
  a <- 0
  for(j in 1:length(listInput)){
    if(gene%in%listInput[[j]]){
      a <- a + 1
    }
  }
  features_count <- c(features_count, a)
}
feature_sum$counts <- features_count
feature_sum <- feature_sum[sort(feature_sum$counts, decreasing=T, index.return=T)$ix,]
feature_sum <- feature_sum[feature_sum$counts>3,]
rownames(feature_sum) <- 1:nrow(feature_sum)
feature_sum$gene <- factor(feature_sum$gene, levels=feature_sum$gene)
ggplot(feature_sum) +
  aes(x = gene, fill = gene, weight = counts) +
  geom_bar() +
  scale_fill_hue(direction = 1) +
  theme_minimal()
