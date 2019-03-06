library(elasticnet)
data(pitprops)
out1<-spca(pitprops,K=6,type="Gram",sparse="penalty",trace=TRUE,para=c(0.06,0.16,0.1,0.5,0.5,0.5))
## print the object out1
out1
out2<-spca(pitprops,K=6,type="Gram",sparse="varnum",trace=TRUE,para=c(7,4,4,1,1,1))
out2
## to see the contents of out2
names(out2)
## to get the loadings
out2$loadings

## to get the figure 2 PEV for PC1
out3<-spca(pitprops,K=1,type="Gram",sparse="penalty",trace=TRUE,para=c(3.4))
out3$pev


## to get the figure 2 PEV for PC2
out4<-spca(pitprops,K=4,type="Gram",sparse="penalty",trace=TRUE,para=c(0.06, 0.16, 0.1, 0.5))
out4$pev



