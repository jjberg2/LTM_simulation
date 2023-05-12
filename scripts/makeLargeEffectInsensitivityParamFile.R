## length of average gene in bp: 67kpb

## total mutation rate to LOFs: 0.0256
## per base mutation rate: 1.3e-8
## approx number of LOF mut opportunities: 0.0256/1e-6 ~ 25000

## setwd('~/Documents/academics/liability-model/LTM_simulation/')
library('viridis')
source('scripts/solveSingleEffect.R')
source('scripts/largeEffectInsensitivityParams.R')




anchor.out <- list()
anchor.pi <- numeric()
anchor.L <- numeric()
for(j in seq_along(my.thr)){
    if(j==1) {
        init.pi <- 1
        init.L <- 500
    } else{
        init.pi <- anchor.pi[j-1]
        init.L <- anchor.L[j-1]
    }
    anchor.out[[j]] <- SolvePoissonGivenThrAndPrev(
        myTheta=my.theta,
        thr=my.thr[j],
        fitCost=anchor.cost,
        Ne=my.Ne,
        h2=h2,
        prev=target.prev,
        alphal=1,
        init.pi=init.pi,
        init.L=init.L,
        verbose.output=TRUE
    )
    anchor.pi[j] <- anchor.out[[j]]$this.pi
    anchor.L[j] <- anchor.out[[j]]$target
}



k <- 0
sim.out <- list()
for(j in seq_along(my.thr)){
    sim.out[[j]] <- list()
    for(i in seq_along(sim.costs)){
        sim.out[[j]][[i]] <- SolvePoissonGivenThr(myTheta=my.theta,thr=my.thr[j],fitCost=sim.costs[i],Ne=my.Ne,L=anchor.L[j],h2=h2,alphal=1,verbose.output=TRUE)
        ## print(paste(j,i,sep=':'))
        k <- k+1
    }
}
print(k)

my.table <- expand.grid(sim.costs,my.thr,my.theta,h2,my.Ne)
colnames(my.table) <- c('cost','thr','theta','h2','Ne')
env.sd <- numeric()
my.mean <- numeric()
pois.prev <- numeric()
deltaR <- numeric()
target.size <- numeric()
for(j in seq_along(sim.out)){
    for(i in seq_along(sim.out[[j]])){
        my.out <- sim.out[[j]][[i]]
        if(length(my.out)==1){
            my.mean <- append(my.mean,NA)
            env.sd <- append(env.sd,NA)
            pois.prev <- append(pois.prev,NA)
            deltaR <- append(deltaR,NA)
            target.size <- append(target.size,NA)
        } else{
            my.mean <- append(my.mean,my.out$mean)
            env.sd <- append(env.sd,my.out$env.sd)
            pois.prev <- append(pois.prev,my.out$prev)
            deltaR <- append(deltaR,my.out$this.pi)
            target.size <- append(target.size,my.out$target)
        }
    }
}
my.rho <- my.table[,'thr']/(2*target.size)
my.b <- 1-2*my.rho
my.gamma <- log((1+my.b)/(1-my.b))
new.table <- cbind(my.table,'target.size'=target.size,'mean'=my.mean,'rho'=my.rho,'b'=my.b,'gamma'=my.gamma,'env.sd'=env.sd,'pois.prev'=pois.prev,'deltaR'=deltaR)
new.table <- new.table[!is.na(new.table$mean),]

write.table(new.table,file='largeEffectInsensitivityParamTable.txt')
