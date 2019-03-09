## The example given in Zou elasticnet package

## library and data
library(elasticnet)
data(pitprops)

out1<-spca(pitprops,K=6,type="Gram",sparse="penalty",trace=TRUE,para=c(0.06,0.16,0.1,0.5,0.5,0.5)) 
## print the object out1 out1
out1

out2<-spca(pitprops,K=6,type="Gram",sparse="varnum",trace=TRUE,para=c(7,4,4,1,1,1)) 
out2 
## to see the contents of out2 
names(out2) 
## to get the loadings 
out2$loadings
print.spca(out2)

# Compare SPCA and prcomp, not possible, cannnot use cor matrix instead of data matrix
out3<-spca(pitprops,K=6,type="Gram",sparse="penalty",trace=TRUE,para=rep(0, times = 13)) 
