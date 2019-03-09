# Read scaled data
X_std <- read.csv("synthetic_data.csv", header=TRUE)

# Perform PCA
pca <- prcomp(X_std, scale = FALSE)

# Loading matrix/ eigenvectors
pca_loading <- pca$rotation


# Compare SPCA and prcomp
spca_nonr <- spca(X_std,K=10,type="predictor",sparse="penalty",trace=TRUE,para=rep(0, times = 10)) 

# Loading matrix/ eigenvectors
spca_nonr_loading <- spca_nonr$loadings

difference = max(abs(pca_loading - spca_nonr_loading))

# Perform proper sparse PCA on the data
spca_2 <- spca(X_std,K=10,type="predictor",sparse="varnum",trace=TRUE,para=c(4,4,2, rep(0, times = 7))) 
spca_2_loading <- spca_2$loadings
spca_2_pev <- spca_2$pev