## length of average gene in bp: 67kpb

## total mutation rate to LOFs: 0.0256
## per base mutation rate: 1.3e-8
## approx number of LOF mut opportunities: 0.0256/1e-6 ~ 25000






## setwd('~/Documents/academics/liability-model/LTM_simulation/')
library('viridis')
source('scripts/solveSingleEffect.R')
source('scripts/largeEffectInsensitivityParams.R')


sim.out <- list()
for(j in seq_along(my.thr)){
    sim.out[[j]] <- list()
    for(i in seq_along(sim.costs)){
        sim.out[[j]][[i]] <- SolvePoissonGivenThr(myTheta=my.theta,thr=my.thr[j],fitCost=sim.costs[i],Ne=my.Ne,L=L,h2=h2,alphal=1,verbose.output=TRUE)
        print(paste(j,i,sep=':'))
    }
}

my.table <- expand.grid(sim.costs,my.thr,L,my.theta,h2,my.Ne)
colnames(my.table) <- c('cost','thr','target.size','theta','h2','Ne')
my.rho <- my.table[,'thr']/(2*my.table[,'target.size'])
my.b <- 1-2*my.rho
my.gamma <- log((1+my.b)/(1-my.b))
env.sd <- numeric()
my.mean <- numeric()
pois.prev <- numeric()
deltaR <- numeric()
for(j in seq_along(sim.out)){
    for(i in seq_along(sim.out[[j]])){
        my.out <- sim.out[[j]][[i]]
        if(length(my.out)==1){
            my.mean <- append(my.mean,NA)
            env.sd <- append(env.sd,NA)
            pois.prev <- append(pois.prev,NA)
            deltaR <- append(deltaR,NA)
        } else{
            my.mean <- append(my.mean,my.out$mean)
            env.sd <- append(env.sd,my.out$env.sd)
            pois.prev <- append(pois.prev,my.out$prev)
            deltaR <- append(deltaR,my.out$this.pi)
        }
    }
}

new.table <- cbind(my.table,'mean'=my.mean,'rho'=my.rho,'b'=my.b,'gamma'=my.gamma,'env.sd'=env.sd,'pois.prev'=pois.prev,'deltaR'=deltaR)
new.table <- new.table[!is.na(new.table$mean),]

write.table(new.table,file='largeEffectInsensitivityParamTable.txt')
