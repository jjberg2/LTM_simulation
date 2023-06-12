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


plot(NA,xlim=c(0,1),ylim=c(0,1),bty='n',type='l')

cols <- viridis(length(split.dps))
for(j in seq_along(split.dps)){
    freqs <- as.numeric(row.names(split.dps[[j]]))
    lines(freqs,rowMeans(split.dps[[j]]))
    g <- gammas[j]
    exp.dps <- exp(g
    lines(freqs
}
