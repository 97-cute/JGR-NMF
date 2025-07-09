# 导入Seurat包
library(Seurat)

# 查看当前工作目录
getwd()

####4.3.2


# 加载必要的库
library(NMF)        # 非负矩阵分解
library(ggplot2)    # 可视化
library(cluster)     # 聚类分析
library(igraph)      # 图形操作
library(dplyr) 
library(flexclust)  # 如果kmeansFunc或后续分析需要此库
library(fpc)        # 用于计算ARI



# 1. 读取数据
# 假设基因表达矩阵和坐标矩阵的路径为 "gene_expression_matrix.csv" 和 "coordinate_matrix.csv"
# 指定要读取的文件所在位置和文件名称
h5_file <- "E:/文献及数据/数据/乳腺癌/ST/Visium_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
# 读取h5格式的文件（使用Read10X_h5函数读取h5格式的单细胞数据文件）
seurat_data <- Read10X_h5(file = h5_file)
expr_matrix<-as.matrix(seurat_data)
TC <- read.csv("E:/文献及数据/数据/乳腺癌/ST/metadata.csv", header = TRUE)
true_labels <- TC[,2]
TC <- read.csv(metadata_file, header = TRUE)
coord_matrix <- as.matrix(read.csv("E:/文献及数据/数据/乳腺癌/ST/spatial/tissue_positions_list.csv", header = F))
coord_matrix<-coord_matrix[,c(1,5:6)]
# 将第一列设置为行名，并删除该列
rownames(coord_matrix) <- coord_matrix[, 1]
coord_matrix <- coord_matrix[, -1]  # 删除第一列
coord_matrix<-coord_matrix[colnames(expr_matrix),]
colnames(coord_matrix) <- c("x", "y")
valid_indices <- !is.na(true_labels)



# 将其转换为数据框以便进行操作
coord_df <- as.data.frame(coord_matrix)
coord_df$x <- as.numeric(coord_df$x)
coord_df$y <- as.numeric(coord_df$y)
str(coord_df)
#####对位置进行处理####
coord_df$x <- coord_df$x - min(coord_df$x) 
coord_df$y <- coord_df$y - min(coord_df$y) 
scaleFactor <- max(coord_df$x, coord_df$y)
coord_df$x <- coord_df$x / scaleFactor
coord_df$y <- coord_df$y / scaleFactor
coord_df <- as.matrix(coord_df[, c("x", "y")])
# 设定零值阈值（50% 为零的基因将被删除）
zero_threshold <- 80
expr_matrix <- expr_matrix[rowMeans(expr_matrix == 0) < zero_threshold, ]
# 检查 expr_matrix 的行名
# 设置行名为基因名，假设你有100个基因
rownames(expr_matrix) <- paste0("Gene", 1:nrow(expr_matrix))
rawData <- list("dataset1" = expr_matrix)
library(rliger)
ifnb_liger <- createLiger(rawData, remove.missing = F, 
                          take.gene.union = T, verbose = F)
class(ifnb_liger)
# 归一化数据
ifnb_liger <- rliger::normalize(ifnb_liger)

class(ifnb_liger)
# 选择变量基因,高可变基因
# ifnb_liger <- selectGenes(ifnb_liger,num.genes = 2000)
ifnb_liger <- selectGenes(ifnb_liger)
# 运行 scaleNotCenter
ifnb_liger <- scaleNotCenter(ifnb_liger)
numK <- 50  # 潜在维度数

#ifnb_liger <-runIntegration(ifnb_liger, k = numK, verbose = FALSE)
Y<-as.matrix(ifnb_liger@datasets[["dataset1"]]@scaleData)

# 初始化 B 矩阵
set.seed(123)
n_rows <- nrow(Y)
n_cols <- ncol(Y)


B <- matrix(runif(n_rows * numK ), nrow=n_rows, ncol=numK )
P <- matrix(runif(n_cols *numK), nrow=n_cols, ncol=numK )
library(RANN)
library(Matrix)
createA <- function(locationList) {
  nSlices = length(locationList)
  AList = list()
  for (islice in 1:nSlices) {
    location = as.data.frame(locationList[[islice]])
    norm_cords = location[, c("x", "y")]
    rownames(norm_cords) <- rownames(location)
    ineibor = 50
    near_data = nn2(norm_cords[, 1:2], k = ineibor)
    neibors = near_data$nn.idx
    neibors = neibors[, -1]
    Nmat = Matrix(0, nrow = nrow(neibors), ncol = nrow(neibors), sparse = TRUE)
    for (icol in 1:ncol(neibors)) {
      edges = data.frame(i = 1:nrow(neibors), j = neibors[, icol])
      adjacency = sparseMatrix(i = as.integer(edges$i),
                               j = as.integer(edges$j), x = 1, dims = rep(nrow(neibors), 2), use.last.ij = TRUE)
      Nmat = Nmat + adjacency
    }
    Nmat = Nmat * t(Nmat)
    rownames(Nmat) = colnames(Nmat) = rownames(norm_cords)
    AList[[islice]] = Nmat
  }
  return(AList)
}
# 通知图构建开始
cat(paste0("### 图构建开始! ...\n"))
colnames(coord_df) <- c("x", "y")
coord_df <- as.data.frame(coord_df)
locationList <- list(coord_df)
AList <- createA(locationList)  # 创建邻接矩阵 A
A<-as.matrix(AList[[1]])
D <- diag(rowSums( A))
L <- D-A


################################自适应knn构建###################

h5_file <- "E:/文献及数据/数据/乳腺癌/ST/Visium_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
# 读取h5格式的文件（使用Read10X_h5函数读取h5格式的单细胞数据文件）
data <- Read10X_h5(file = h5_file)
data<-as.matrix(data)


#data<-data[,grep("epi3",colnames(data))]
# process the data using Seurat
library(Seurat)
library(magrittr)
library(dplyr)
# obj<-CreateSeuratObject(counts=data,min.cells=3,min.features=200)
obj<-CreateSeuratObject(counts=data,min.cells=3,min.features=0)
obj <- NormalizeData(obj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims=1:10)
#构建自适应K近邻图
library(aKNNO)
obj<-FindNeighbors_aKNNO(obj,verbose = T)
# 基于优化自适应k近邻图的聚类
obj<-FindClusters(obj,graph.name="aKNN_O",verbose=T)
# visualization
DimPlot(obj,label=T,group.by="aKNN_O_res.0.8")+ggtitle("aKNNO")+NoLegend() -> p_aKNNO
#DotPlot(obj,features=c("Krt19","Krt8","Top2a","Mki67","Muc2","Clca1","Lgr5","Ascl2","Fabp6","Apoa1","Gip","Chga","Chgb","Adh6a","Hck","Lrmp","Defa21","Defa22"))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) -> p_marker
#p_aKNNO+p_marker
p_aKNNO
#提取细胞与细胞分类
a<-obj@meta.data
aKNN_graph <- obj@graphs$aKNN_O
A1 <- as.matrix(aKNN_graph)######构建拉普拉斯的矩阵
class(A1)
#########################根据自适应图构建拉普拉斯图########################
#  A1是邻接矩阵为密集矩阵，类似"dgCMatrix"，是稀疏矩阵。"matrix"或"array"，则是密集矩阵

# 1. 计算度矩阵 D
degree_vector <- rowSums(A1)  # 计算每个节点的度
D1 <- diag(degree_vector)  # 创建度矩阵 D
# 2. 计算拉普拉斯矩阵 L
laplacian_matrix <- D1 - A1


rm(list = setdiff(ls(), c("Y","A","D","A1","D1","B","P","true_labels")))

# 加载必要的库
library(NMF)        # 非负矩阵分解
library(ggplot2)    # 可视化
library(cluster)     # 聚类分析
library(igraph)      # 图形操作
library(dplyr) 
library(flexclust)  # 如果kmeansFunc或后续分析需要此库
library(fpc)        # 用于计算ARI


# 加载必要的库
library(proxy)
k_clusters <- 8   # 聚类数量


kmeansFunc<-function (data, k) 
{
  set.seed(12345678)
  if (nrow(data) < 3e+05) {
    numStart = 100
  } else {
    numStart = 1
  }
  cl <- suppressWarnings(try(kmeans(data, k, nstart = numStart, 
                                    iter.max = 100), silent = TRUE))
  
  # 返回聚类标签和聚类质心
  return(list(kmeans = cl$cluster, centers = cl$centers))
}




GNMF<-function(B, P, max_iter,epsilon,a1,a2){
  for (iter in 1:max_iter) {
    # 保存前一次的W和H
    B_old <- B
    P_old <- P
    # 更新H矩阵
    B <- (Y%*%P)*B/(B%*%t(P)%*%P)  
    B[B < 0] <- 0  # 保证非负性
    # 更新W矩阵
    P <- (t(Y)%*%B+a1*(A%*%P)+a2*(A1%*%P))*P/(P%*%t(B)%*%B+a1*(D%*%P)+a2*(D1%*%P))
    P[P < 0] <- 0  # 保证非负性
    # 计算参数变化
    B_change <-max(abs(B - B_old))
    P_change <- max(abs(P - P_old))
    # 收敛条件检查
    if (P_change < epsilon) {
      #cat("收敛于迭代", iter, "次", "\n")
      break
    }
    # 打印当前损失值（可选）
    # if (iter %% 50 == 0) {
    #   cat("迭代", iter, "次", "\n")
    # }
  }
  return(P = P)
}

library(flexclust)
library(fpc) 

start_time <- Sys.time()


# # 定义 a1 和 a2 的值范围
# a1_values <- seq(0.1, 1, length.out = 10)
# a2_values <- seq(0.1, 1, length.out = 10)
# 
# # 嵌套循环遍历 a1 和 a2 的值
# for (a1 in a1_values) {
#   for (a2 in a2_values) {
#     # 调用 GNMF 函数进行计算
#     result <- GNMF(B, P, max_iter = 500, epsilon = 0.05, a1 = a1, a2 = a2)
#     # 使用 kmeansFunc 对结果进行聚类
#     s1 <- kmeansFunc(result, k_clusters)
#     # 提取聚类结果
#     J1 <- s1$kmeans  # 假设这是正确的提取方式
#     # 计算 Rand Index
#     ri <- randIndex(table(J1, true_labels), correct = TRUE)
#     
#     # 输出 a1, a2 和对应的 Rand Index
#     cat("a1 =", a1, ", a2 =", a2, ", Rand Index =", ri, "\n")
#   }
# }
# 定义 a1 和 a2 的值范围

# a1_values <-seq(0.1, 1, length.out = 10)
a1_values <- c(0.2)
a2_values <- c(0)
# a2_values <- seq(0, 1, length.out = 11)
# 创建一个空的数据框，用于存储结果
results_df <- data.frame(a1 = numeric(0), a2 = numeric(0), RandIndex = numeric(0))
start_time <- Sys.time()
for (a1 in a1_values) {
  for (a2 in a2_values) {
    # 调用 GNMF 函数进行计算
    result <- GNMF(B, P, max_iter = 500, epsilon = 0.05, a1 = a1, a2 = a2)
    
    # 使用 kmeansFunc 对结果进行聚类
    s1 <- kmeansFunc(result, k_clusters)
    
    # 提取聚类结果
    J1 <- s1$kmeans  # 假设这是正确的提取方式
    
    # 计算 Rand Index
    ri <- randIndex(table(J1, true_labels), correct = TRUE)
    
    # 将 a1, a2 和 Rand Index 的值添加到数据框中
    results_df <- rbind(results_df, data.frame(a1 = a1, a2 = a2, RandIndex = ri))
    
    # 输出 a1, a2 和对应的 Rand Index
    cat("a1 =", a1, ", a2 =", a2, ", Rand Index =", ri, "\n")
  }
}
end_time <- Sys.time()
# 安装必要的包（如果还没有安装的话）
# devtools::install_local("E:/文献及数据/第二篇/代码/BayesSpace/aricode-master.zip")
# install.packages("aricode")
# # 安装并加载cluster包（如果需要使用其他聚类评估方法）
# install.packages("cluster")



# 
# 
# # 加载包
# library(aricode)
# library(cluster)
# library(fpc) 
# # 计算归一化互信息 (NMI)
# nmi_value <- NMI(as.factor(true_labels), as.factor(J1))
# print(nmi_value)
# ARI <- randIndex(table(J1, true_labels), correct = TRUE)
# print(ARI)
# 
# write.csv(J1,"E:/文献及数据/第二篇/结果/乳腺癌圆/提出_cluster.csv")
# 
# 
# 
write.csv(J1,"E:/文献及数据/第二篇/结果/乳腺癌圆/提出_消融cluster2.csv")



