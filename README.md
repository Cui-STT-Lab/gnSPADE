# Reference-free deconvolution incorporating gene network structures in spatial transcriptomics
gnSPADE is a reference-free SPAtial DEconvolution method incorporating gene network structures for spatial transcriptomics. gnSPADE requires the gene count matrix and gene network structures derived from existing databases such as the PPI network, KEGG pathway and Reactome. It builds upon the latent Dirichlet allocation (LDA) model to recover cell-type compositions and transcriptional profiles within each spatial location, without relying on external single-cell reference information.
![gnSPADE_overview](https://github.com/user-attachments/assets/19c719ed-3380-4b8a-b936-eaea8b84d046)

## Installation

You can install the development version of gwSPADE from [GitHub](https://github.com/Cui-STT-Lab/gwSPADE) with:

```r
# install.packages("devtools")
devtools::install_github("Cui-STT-Lab/gnSPADE")
```

#### **Run gnSPADE with Example Data**
## Run gnSPADE with Example Data

The function WLDA takes the spatial transcriptomics data matrix `corpus` (spots x genes) as inputs.

The `corpus` needs to be non-transformed counts (or round the data to the nearest integers).
```r
n_topics <- k <- 4
n_docs <- 1000
num_nodes = n_vocab <- V <- 100
blocks <- list(1:10, 21:30, 41:50, 61:70)
gene_names <- paste0("gene", 1:100)

library(MASS)
n <- 10   # Dimension of the covariance matrix (number of variables)
rho <- 0.8  # Lag-1 autocorrelation coefficient
# Generate the AR(1) covariance matrix
ar1_cov_matrix <- matrix(0, n, n)  # Initialize an n x n matrix
for (i in 1:n) {
  for (j in 1:n) {
    ar1_cov_matrix[i, j] <- rho^abs(i - j)
  }
}

create_cov_matrix2 <- function(block_indices, size, matrix) {
  cov_matrix <- diag(size)  # Start with an identity matrix
  cov_matrix[block_indices, block_indices] <- matrix
  return(cov_matrix)
}
cov_matrices <- lapply(blocks, create_cov_matrix2, size = num_nodes, ar1_cov_matrix)

cor_matrix = matrix(0,100,100)
for (i in seq_along(cov_matrices)) {
  cor_matrix = cor_matrix+cov_matrices[[i]]
}
diag(cor_matrix)=0
rownames(cor_matrix) <- gene_names
colnames(cor_matrix) <- gene_names

correlation_list <- lapply(1:nrow(cor_matrix), function(i) {
  # Find the indices of genes with correlation > 0.5 (including itself)
  correlated_indices <- which(cor_matrix[i, ] > 0.5)
  # Get the gene names corresponding to these indices
  correlated_genes <- gene_names[correlated_indices]
  return(correlated_genes)
})
# Assign names to the list
names(correlation_list) <- gene_names
diag(cor_matrix)=1

set.seed(2012)
a1 = exp(mvrnorm(1, mu = rep(2,100), Sigma = cov_matrices[[1]]))
a2 = exp(mvrnorm(1, mu = rep(2,100), Sigma = cov_matrices[[2]]))
a3 = exp(mvrnorm(1, mu = rep(2,100), Sigma = cov_matrices[[3]]))
a4 = exp(mvrnorm(1, mu = rep(2,100), Sigma = cov_matrices[[4]]))

b1 = gtools::rdirichlet(1, a1)
b2 = gtools::rdirichlet(1, a2)
b3 = gtools::rdirichlet(1, a3)
b4 = gtools::rdirichlet(1, a4)
beta = rbind(b1,b2,b3,b4)
beta.true = beta
rowSums(beta.true)
rownames(beta.true) = paste0('Cell', seq(nrow(beta.true)))
colnames(beta.true) = paste0('gene', seq(ncol(beta.true)))
theta.true <- gtools::rdirichlet(n_docs, rep(1/n_topics, n_topics))
colnames(theta.true) = rownames(beta.true)
rownames(theta.true) = paste0('Doc', seq(nrow(theta.true)))
docs <- matrix(0, nrow = n_docs, ncol = n_vocab)
ksai <- 1000 # average words per doc
for (i in 1:n_docs) {
  # draw topics for each word
  tops <- rmultinom(1, rpois(1, ksai), theta.true[i, ]) # draw words
  for (j in 1:n_topics) {
    docs[i, ] <- docs[i, ] + rmultinom(1, tops[j], beta.true[j, ])
  }
}
rownames(docs) = paste0('Doc', seq(nrow(docs)))
colnames(docs) = paste0('gene', seq(ncol(docs)))

corpus = docs

#devtools::install_github("Cui-STT-Lab/gnSPADE")
library(gnSPADE)
gnModel = WRLDA(docs, k=4, lambda = 1, neiPath = correlation_list)
corPlot_matchall(gnModel$phi, beta.true, gnModel$theta, theta.true)
```
<img width="647" alt="image" src="https://github.com/user-attachments/assets/5281958c-347a-4ca2-b44b-e622dd0f4632" />


