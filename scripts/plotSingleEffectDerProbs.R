library('viridis')
input <- get(load('smallEffectInsensitivityDerProbs.Rdata'))

## extract parameter table
pt <- do.call(rbind,sapply(input,function(X)X[1]))
gammas <- unique(pt$gamma)

## extract derived state probs
dp <- sapply(input,function(X)X[[2]][3,])
split.dps <- list() ## gah I don't remember basic R functions so here's a for loop
for(i in seq_along(gammas)){
    split.dps[[i]] <- dp[,pt$gamma==gammas[i]]
}



pdf('../figures/singleEffectDerProbs.pdf',width=8,height=8)
plot(NA,xlim=c(0,1),ylim=c(0,1),bty='n',type='l')
cols <- viridis(length(split.dps))
for(j in 1:5){
    f <- as.numeric(row.names(split.dps[[j]]))
    points(f,rowMeans(split.dps[[j]]),col=cols[j],pch=20,cex=0.3)
    g <- gammas[j]
    exp.dps <- (exp(g)-exp(g*f))/(exp(g)-1)
    lines(f,exp.dps,col=cols[j],lwd=2,lty=2)
}
dev.off()
