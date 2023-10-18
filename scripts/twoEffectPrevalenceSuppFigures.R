rm(list=ls())
library('MetBrewer')
source('scripts/solveTwoEffect.R')
recover.flag <- FALSE


## main single effect params
L <- 1e7
Ne <- 1e3
u <- 1e-8
C <- 1
my.bt <- 0.8
h2 <- 1/2

## ratio of large effect observed var to small effect
var.ratio <- c(1/2,1,2)

as <- 1 
my.als <- exp(seq(log(4),log(200),length.out=200))
theta <- 4*Ne*u


solns <- list()
tmp.output <- list()
output <- list()
for( j in seq_along(var.ratio)){
    last.Ll <- numeric()
    last.Ll[1] <- 90000
    last.gs <- numeric()
    last.gs[1] <- (L - last.Ll[1])/L
    last.bs <- numeric()
    last.bs[1] <- 0.7
    for( i in seq_along(my.als)){
        solns[[i]] <- solveTwoEffect(
            bs=last.bs[i],
            bt=my.bt,
            Ne=Ne,
            as=as,
            al=my.als[i],
            L=L,
            gs=last.gs[i],
            last.tstar=NULL,
            h2=h2,
            u=u,
            C=C,
            LL.soln=TRUE,
            var.ratio=var.ratio[j],
            equalize.observed.vars=TRUE
        )
        tmp.output[[i]] <- makeOutput(solns[[i]]) 
        last.gs[i+1]  <- tmp.output[[i]]['gs']
        last.bs[i+1] <- tmp.output[[i]]['bs']
        if(i %% 10 == 0) print(i)
    }
    output[[j]] <- as.data.frame(do.call(rbind,tmp.output))
}





norm.ft <- 1/(4*Ne*C)*log((1+my.bt/(1-my.bt)))
norm.Vg <- 8*Ne*L*u*my.bt / log((1+my.bt/(1-my.bt)))
norm.Vt <- norm.Vg/h2
norm.astd <- 1 / sqrt ( norm.Vt )
norm.dens <- 2*L*u*my.bt*norm.astd / (h2*C) 
norm.prev <- 1-pnorm(dnorminv(norm.dens))




my.cols = met.brewer('Hiroshige',3)
my.als <- c(1,output[[1]][,'al'])
tmp.prev <- sapply(output,function(X) X$prev)
my.prevs <- rbind(rep(norm.prev,ncol(tmp.prev)),tmp.prev)
plot(
    NA,
    xlim = range(my.als),
    ylim = c(0,max(my.prevs)),
    type='l'
)
for(i in 1:ncol(my.prevs)){
    lines(
        my.als,
        my.prevs[,i],
        col = my.cols[i],
        lty = 1,
        lwd = 2
    )
}
