library(viridis)



small <- get(load('smallEffectInsensitivityResultsTable.Rdata'))
large <- get(load('largeEffectInsensitivityResultsTable.Rdata'))





small.thr <- unique(small$thr)
small$thetaL <- small$theta*small$target.size
small$sim.genVarPerSite <- small$sim.genVar/small$thetaL
small.cols <- seq_along(small.thr)




large.thr <- unique(large$thr)
large$thetaL <- large$theta*large$target.size
large$sim.genVarPerSite <- large$sim.genVar/large$thetaL
large.cols <- viridis(length(large.thr))

large$pois.Var  <- 
large$sim.stdGenVarPerSite  <- NA
for(i in 1:length(large.thr)){
    these <- large$thr == large.thr[i]
    large$sim.stdGenVarPerSite[these]  <- large$sim.genVarPerSite[these]/tail(large$sim.genVarPerSite[these],1)
}






plot(
    NA,
    xlim=c(0,1),
    ylim=c(0,max(small$sim.genVarPerSite))
)
for(i in 1:length(small.thr)){
    these <- small$thr == small.thr[i]
    points(
        small$cost[these],
        small$sim.genVarPerSite[these],
        col=small.cols[i],
        pch=20
    )
}



plot(
    NA,
    xlim=c(0,1),
    ylim=c(0,max(large$sim.genVarPerSite))
)
for(i in 1:length(large.thr)){
    these <- large$thr == large.thr[i]
    points(
        large$cost[these],
        large$sim.genVarPerSite[these],
        col=large.cols[i],
        pch=20
    )
}













op <- par(mfrow=c(1,2))

plot(
    NA,
    xlim=c(0,1),
    ylim=c(0,max(large$sim.stdGenVarPerSite)*1.05),
    xlab='Fitness Cost',
    ylab='Relative Per Site Average Heterozygosity'
)
abline(h = 1, lty = 3 )
for(i in 1:length(large.thr)){
    these <- large$thr == large.thr[i]
    points(
        large$cost[these],
        large$sim.stdGenVarPerSite[these],
        col=large.cols[i],
        pch=20
    )
    lines(
        large$cost[these],
        large$sim.stdGenVarPerSite[these],
        col=large.cols[i]
    )
}



plot(
    NA,
    xlim=c(0,1),
    ylim=c(min(large$pois.prev)*0.95,max(large$pois.prev)*1.05),
    xlab='Fitness Cost',
    ylab='Prevalence',
    log='y'
)
for(i in 1:length(large.thr)){
    these <- large$thr == large.thr[i]
    points(
        large$cost[these],
        large$sim.prev[these],
        col=large.cols[i],
        pch=20
    )
    lines(
        large$cost[these],
        large$sim.prev[these],
        col=large.cols[i]
    )
    lines(
        large$cost[these],
        large$pois.prev[these],
        col=large.cols[i],
        lty=2
    )
}
