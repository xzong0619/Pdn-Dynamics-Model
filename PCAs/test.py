import numpy as np
from sklearn.decomposition import SparsePCA
from sklearn.decomposition import PCA 
from sklearn.preprocessing import StandardScaler 
'''
Create the dataset
'''
# set the random seed

np.random.seed(seed=0)

# Number of data points
n = 100

# Number of parameters
p = 10

V1 = np.random.randn(n) * np.sqrt(290)
V2 = np.random.randn(n) * np.sqrt(300)
eps = np.random.randn(n) 
V3 = -0.3 * V1 + 0.925 * V2 + eps

X = np.zeros((n, p))
for i in [0,1,2,3]:
    X[:,i] = V1 + np.random.randn(n) 
for i in [4,5,6,7]:
    X[:,i] = V2 + np.random.randn(n) 
for i in [8,9]:
    X[:,i] = V3 + np.random.randn(n) 

# Normalize the data
X_std = StandardScaler().fit_transform(X)
    
#%% PCA provide by Sklearn
pca = PCA()    
Xpc = pca.fit_transform(X_std) 
eig_vals = pca.explained_variance_ #eigenvalues 
eig_vecs = pca.components_  # eigenvector, PC loadings
var_exp = pca.explained_variance_ratio_ #explained variance ratio
cum_var_exp = np.cumsum(var_exp) #cumulative variance ratio

PC_loadings = []
# Create a set of bars at each position
for i, pci in enumerate(eig_vecs):
    
    vals = np.abs(pci)/np.sum(np.abs(pci))
    PC_loadings.append(vals)
    
PC_loadings = np.array(PC_loadings)

#%% Use PCA SVD descomposition
u, s, vh = np.linalg.svd(X_std, full_matrices=True)
smat = np.zeros((n, p))
smat[:p, :p] = np.diag(s)
z = np.dot(u, smat) # principle component matrix, = Xpc
v = np.transpose(vh)

flag1 = np.allclose(X_std, np.dot(z, vh)) # check if two matrices are the same
print(flag1)

flag2 = np.allclose(np.abs(z), np.abs(Xpc)) # interesting that some signs are flipped 
print(flag2)

flag3 = np.allclose(np.abs(eig_vecs), np.abs(vh)) # the row of vh correspond to loading of the PCs
print(flag3)

#%% Use SPCA with all zero parameters
# SPCA seems completely off 
SPCA= SparsePCA(alpha = 0, ridge_alpha  = 0,  n_components=p, random_state=0)
SPCA.fit(X_std) 
Xspca = SPCA.transform(X_std)
spca_c = SPCA.components_ 
flag4 = np.allclose(Xpc, Xspca)
print(flag4)

flag5 = np.allclose(eig_vecs, spca_c)
print(flag5)