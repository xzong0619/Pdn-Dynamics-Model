## Use to plot figure 2 in Zou SPCA paper

## library and data
library(elasticnet)
data(pitprops)

# Lambda1 used in the paper
lambda1_paper <- c(0.06, 0.16, 0.1, 0.5, 0.5, 0.5) 



## Plot PC function,
## takes in the index of PC and the range of lambda1 used in plotting 
plot_PCi <- function(iPC, lambda_range) {
  PCn <- paste("PC",iPC, sep = '')
  pc_pev <- c()
  
  # lambda1 used to plot for this PC
  lambda1 <- seq(0, lambda_range, length = 20)
  
  # Set the lambdas for previous PCs
  if (iPC >= 1){
    # previous lambdas, Set all other PC pently as 0, no lasso regularization
    pre_lambdas = rep(0,times = iPC-1)
    # previous lambdas, Set all other PC pently as the previously determined values
    pre_lambdas = lambda1_paper[1:iPC-1]
  } else {
    pre_lambdas = c()
  }
  
  # Loop through each lambda for plotting
  for (i in 1:length(lambda1)){
    
    pc <- spca(pitprops, K=iPC, type="Gram", sparse="penalty", trace=TRUE, 
               para = c(pre_lambdas, lambda1[i])) 
    
    pc_pev[i] <- pc$pev[iPC]
  }
  
  # save the plot to png
  png(paste(PCn, ".png", sep = '') ,width=5,height=5,units="in",res=900)
  plot(lambda1, pc_pev, main = PCn, ylab = 'PEV',  col = "red")
  lines(lambda1, pc_pev, type ='l', lty=2)
  lines(lambda1_paper[iPC] * rep(1, times = length((lambda1))), pc_pev, type ='l', lty=2, col = 'blue')
  dev.off()
}

# PC 1
plot_PCi(1,3.5)

# PC 2
plot_PCi(2,3)

# PC 3
plot_PCi(3,2)

# PC 4
plot_PCi(4,2)

# PC 5
plot_PCi(5,1.1)

# PC 6
plot_PCi(6,1.1)