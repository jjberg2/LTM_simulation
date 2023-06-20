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



pdf('figures/singleEffectDerProbs1.pdf',width=8,height=8)
draw.n <- length(split.dps)
plot(
    NA,
    xlim=c(0,1),
    ylim=c(0,1),
    bty='n',
    type='l',
    xlab='Frequency of Risk Increasing Allele',
    ylab='P(risk allele is derived | frequency of risk allele)'
)
cols <- viridis(length(split.dps))
for(j in 1:draw.n){
    f <- as.numeric(row.names(split.dps[[j]]))
    unambiguous <- f!=0.5
    points(f[unambiguous],rowMeans(split.dps[[j]])[unambiguous],col=cols[j],pch=20,cex=0.3)
    g <- gammas[j]
    exp.dps <- (exp(g)-exp(g*f))/(exp(g)-1)
    lines(f,exp.dps,col=cols[j],lwd=2,lty=2)
}
legend(
    'bottomleft',
    legend=small.gamma[1:draw.n],
    pch=20,
    col=cols
)
text(
    x=0.2,
    y=0.4,
    labels=expression((e^y-e^(y*x))/(e^y-e))
)
dev.off()






pdf('figures/singleEffectDerProbs2.pdf',width=8,height=8)
draw.n <- 4
plot(
    NA,
    xlim=c(0,1),
    ylim=c(0,1),
    bty='n',
    type='l',
    xlab='Frequency of Risk Increasing Allele',
    ylab='P(risk allele is derived | frequency of risk allele)'
)
cols <- viridis(length(split.dps))
for(j in 1:draw.n){
    f <- as.numeric(row.names(split.dps[[j]]))
    unambiguous <- f!=0.5
    points(f[unambiguous],rowMeans(split.dps[[j]])[unambiguous],col=cols[j],pch=20,cex=0.3)
    g <- gammas[j]
    exp.dps <- (exp(g)-exp(g*f))/(exp(g)-1)
    lines(f,exp.dps,col=cols[j],lwd=2,lty=2)
}
legend(
    'bottomleft',
    legend=small.gamma[1:draw.n],
    pch=20,
    col=cols
)
text(
    x=0.2,
    y=0.4,
    labels=expression((e^y-e^(y*x))/(e^y-e))
)
dev.off()









tf <- F
if(tf){
    draw.n <- length(split.dps)
} else {
    draw.n <- 4
}
