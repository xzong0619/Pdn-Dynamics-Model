# Following the example of 
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

library("FactoMineR")
library("factoextra")

# Get the data
data(decathlon2)
head(decathlon2)

decathlon2.active <- decathlon2[1:23, 1:10]
head(decathlon2.active[, 1:6], 4)

# Perform PCA
res.pca <- PCA(decathlon2.active, graph = FALSE)

print(res.pca)

# Extra eigenvalues
eig.val <- get_eigenvalue(res.pca)
eig.val

# Make Scree Plot
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

# Extra the active variables
var <- get_pca_var(res.pca)
var

# Coordinates, I think these are loadings
head(var$coord)
# Plot coordinates of variables
fviz_pca_var(res.pca, col.var = "black")


# Cos2: quality on the factore map
head(var$cos2)
library("corrplot")
corrplot(var$cos2, is.corr=FALSE)

# Contributions to the principal components
head(var$contrib)

# Total cos2 of variables on Dim.1 and Dim.2
fviz_cos2(res.pca, choice = "var", axes = 1:2)

# Color by cos2 values: quality on the factor map
fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
