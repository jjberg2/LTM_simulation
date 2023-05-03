# setwd('~/Documents/academics/liability-model/LTM_simulation/')
library('viridis')
source('scripts/solveSingleEffect.R')

my.thr <- c(1,2,5,10,15,20,50)
my.Ne <- 5000
L <- 10000
h2=0.5
mu <- 1e-6
my.theta <- 4*my.Ne*mu



## sim.costs <- seq(0.01,1,by=0.01)
sim.costs <- c(0.25,0.5,0.75,0.95)
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
for(j in seq_along(sim.out)){
    for(i in seq_along(sim.out[[j]])){
        my.out <- sim.out[[j]][[i]]
        if(length(my.out)==1){
            my.mean <- append(my.mean,NA)
            env.sd <- append(env.sd,NA)
            pois.prev <- append(pois.prev,NA)
        } else{
            my.mean <- append(my.mean,my.out$mean)
            env.sd <- append(env.sd,my.out$env.sd)
            pois.prev <- append(pois.prev,my.out$prev)
        }
    }
}
new.table <- cbind(my.table,'rho'=my.rho,'b'=my.b,'gamma'=my.gamma,'env.sd'=env.sd,'pois.prev'=pois.prev)

write.table(new.table,file='largeEffectInsensitivityParamTable.txt')
