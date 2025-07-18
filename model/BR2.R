# Load Seurat package
library(Seurat)

# Check current working directory
getwd()

####4.3.2

# Load necessary libraries
library(NMF)        # Non-negative Matrix Factorization
library(ggplot2)    # Visualization
library(cluster)    # Cluster analysis
library(igraph)     # Graph operations
library(dplyr) 
library(flexclust)  # Required for kmeansFunc or subsequent analysis
library(fpc)        # For calculating ARI

# 1. Read data
# Assuming gene expression matrix and coordinate matrix paths are "gene_expression_matrix.csv" and "coordinate_matrix.csv"
# Specify file locations and names
h5_file <- "./Visium_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
# Read h5 format file (using Read10X_h5 function to read single-cell data in h5 format)
seurat_data <- Read10X_h5(file = h5_file)
expr_matrix <- as.matrix(seurat_data)
TC <- read.csv("./metadata.csv", header = TRUE)
true_labels <- TC[,2]
TC <- read.csv(metadata_file, header = TRUE)
coord_matrix <- as.matrix(read.csv("./spatial/tissue_positions_list.csv", header = F))
coord_matrix <- coord_matrix[,c(1,5:6)]
# Set first column as row names and remove that column
rownames(coord_matrix) <- coord_matrix[, 1]
coord_matrix <- coord_matrix[, -1]  # Remove first column
coord_matrix <- coord_matrix[colnames(expr_matrix),]
colnames(coord_matrix) <- c("x", "y")
valid_indices <- !is.na(true_labels)

# Convert to dataframe for operations
coord_df <- as.data.frame(coord_matrix)
coord_dfx<−as.numeric(coorddfx <- as.numeric(coord_dfx)
coord_dfy<−as.numeric(coorddfy <- as.numeric(coord_dfy)
str(coord_df)
##### Process coordinates ####
coord_dfx<−coorddfx <- coord_dfx - min(coord_df$x) 
coord_dfy<−coorddfy <- coord_dfy - min(coord_df$y) 
scaleFactor <- max(coord_dfx,coorddfx, coord_dfy)
coord_dfx<−coorddfx <- coord_dfx / scaleFactor
coord_dfy<−coorddfy <- coord_dfy / scaleFactor
coord_df <- as.matrix(coord_df[, c("x", "y")])
# Set zero-value threshold (genes with >80% zeros will be removed)
zero_threshold <- 80
expr_matrix <- expr_matrix[rowMeans(expr_matrix == 0) < zero_threshold, ]
# Check row names of expr_matrix
# Set row names as gene names, assuming 100 genes
rownames(expr_matrix) <- paste0("Gene", 1:nrow(expr_matrix))
rawData <- list("dataset1" = expr_matrix)
library(rliger)
ifnb_liger <- createLiger(rawData, remove.missing = F, 
                          take.gene.union = T, verbose = F)
class(ifnb_liger)
# Normalize data
ifnb_liger <- rliger::normalize(ifnb_liger)

class(ifnb_liger)
# Select variable genes, highly variable genes
# ifnb_liger <- selectGenes(ifnb_liger,num.genes = 2000)
ifnb_liger <- selectGenes(ifnb_liger)
# Run scaleNotCenter
ifnb_liger <- scaleNotCenter(ifnb_liger)
numK <- 50  # Number of latent dimensions

#ifnb_liger <-runIntegration(ifnb_liger, k = numK, verbose = FALSE)
Y <- as.matrix(ifnb_liger@datasets[["dataset1"]]@scaleData)

# Initialize B matrix
set.seed(123)
n_rows <- nrow(Y)
n_cols <- ncol(Y)

B <- matrix(runif(n_rows * numK), nrow=n_rows, ncol=numK)
P <- matrix(runif(n_cols *numK), nrow=n_cols, ncol=numK)
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
# Notify graph construction start
cat(paste0("### Graph construction started! ...\n"))
colnames(coord_df) <- c("x", "y")
coord_df <- as.data.frame(coord_df)
locationList <- list(coord_df)
AList <- createA(locationList)  # Create adjacency matrix A
A <- as.matrix(AList[[1]])
D <- diag(rowSums(A))
L <- D-A

################################ Adaptive knn construction ###################

h5_file <- "E:/文献及数据/数据/乳腺癌/ST/Visium_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
# Read h5 format file (using Read10X_h5 function to read single-cell data in h5 format)
data <- Read10X_h5(file = h5_file)
data <- as.matrix(data)

# Process data using Seurat
library(Seurat)
library(magrittr)
library(dplyr)
obj <- CreateSeuratObject(counts=data, min.cells=3, min.features=200)
obj <- NormalizeData(obj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims=1:10)
# Construct adaptive K-nearest neighbor graph
library(aKNNO)
obj <- FindNeighbors_aKNNO(obj, verbose = T)
# Clustering based on optimized adaptive k-nearest neighbor graph
obj <- FindClusters(obj, graph.name="aKNN_O", verbose=T)
# Visualization
DimPlot(obj, label=T, group.by="aKNN_O_res.0.8")+ggtitle("aKNNO")+NoLegend() -> p_aKNNO
p_aKNNO
# Extract cell classifications
a <- obj@meta.data
aKNN_graph <- obj@graphs$aKNN_O
A1 <- as.matrix(aKNN_graph) ###### Construct Laplacian matrix
class(A1)
######################### Construct Laplacian graph based on adaptive graph ########################
# A1 is adjacency matrix as dense matrix, similar to "dgCMatrix" which is sparse matrix. "matrix" or "array" are dense matrices.

# 1. Calculate degree matrix D
degree_vector <- rowSums(A1)  # Calculate degree for each node
D1 <- diag(degree_vector)  # Create degree matrix D
# 2. Calculate Laplacian matrix L
laplacian_matrix <- D1 - A1

rm(list = setdiff(ls(), c("Y","A","D","A1","D1","B","P","true_labels")))

# Load necessary libraries
library(NMF)        # Non-negative Matrix Factorization
library(ggplot2)    # Visualization
library(cluster)    # Cluster analysis
library(igraph)     # Graph operations
library(dplyr) 
library(flexclust)  # Required for kmeansFunc or subsequent analysis
library(fpc)        # For calculating ARI

# Load necessary libraries
library(proxy)
k_clusters <- 8   # Number of clusters

kmeansFunc <- function(data, k) {
  set.seed(12345678)
  if (nrow(data) < 3e+05) {
    numStart = 100
  } else {
    numStart = 1
  }
  cl <- suppressWarnings(try(kmeans(data, k, nstart = numStart, 
                                    iter.max = 100), silent = TRUE))
  
  # Return cluster labels and centroids
  return(list(kmeans = clcluster,centers=clcluster, centers = clcenters))
}

GNMF <- function(B, P, max_iter, epsilon, a1, a2){
  for (iter in 1:max_iter) {
    # Save previous W and H
    B_old <- B
    P_old <- P
    # Update H matrix
    B <- (Y%*%P)*B/(B%*%t(P)%*%P)  
    B[B < 0] <- 0  # Ensure non-negativity
    # Update W matrix
    P <- (t(Y)%*%B + a1*(A%*%P) + a2*(A1%*%P))*P/(P%*%t(B)%*%B + a1*(D%*%P) + a2*(D1%*%P))
    P[P < 0] <- 0  # Ensure non-negativity
    # Calculate parameter changes
    B_change <- max(abs(B - B_old))
    P_change <- max(abs(P - P_old))
    # Check convergence condition
    if (P_change < epsilon) {
      break
    }
  }
  return(P = P)
}

library(flexclust)
library(fpc) 

start_time <- Sys.time()

# Define a1 and a2 value ranges
a1_values <- c(0.2)
a2_values <- c(0)
# Create empty dataframe to store results
results_df <- data.frame(a1 = numeric(0), a2 = numeric(0), RandIndex = numeric(0))
start_time <- Sys.time()
for (a1 in a1_values) {
  for (a2 in a2_values) {
    # Call GNMF function
    result <- GNMF(B, P, max_iter = 500, epsilon = 0.05, a1 = a1, a2 = a2)
    
    # Use kmeansFunc for clustering
    s1 <- kmeansFunc(result, k_clusters)
    
    # Extract clustering results
    J1 <- s1$kmeans  # Assuming this is the correct extraction method
    
    # Calculate Rand Index
    ri <- randIndex(table(J1, true_labels), correct = TRUE)
    
    # Add a1, a2 and Rand Index values to dataframe
    results_df <- rbind(results_df, data.frame(a1 = a1, a2 = a2, RandIndex = ri))
    
    # Output a1, a2 and corresponding Rand Index
    cat("a1 =", a1, ", a2 =", a2, ", Rand Index =", ri, "\n")
  }
}
end_time <- Sys.time()

write.csv(J1,"./cluster2.csv")
