library(viridis)



small <- get(load('smallEffectInsensitivityResultsTable.Rdata'))
large <- get(load('largeEffectInsensitivityResultsTable.Rdata'))

large.cost <- sort(unique(large$cost))
large.thr <- unique(large$thr)
large$thetaL <- large$theta*large$target.size
large$sim.genVarPerSite <- large$sim.genVar/large$thetaL
large.cols <- viridis(length(large.thr))




## this is sort of wrong but close
large$sim.addRiskVar <- large$sim.genVarPerSite*large$sim.deltaR^2



large$pois.Var  <- large$mean
large$sim.stdGenVarPerSite  <- NA
large$sim.scaledDeltaR  <- NA
large$sim.stdAddRiskVar <- NA
for(i in 1:length(large.thr)){
    these <- large$thr == large.thr[i]
    large$sim.stdGenVarPerSite[these]  <- large$sim.genVarPerSite[these]/tail(large$sim.genVarPerSite[these],1)
    large$sim.scaledDeltaR[these]  <- large$sim.deltaR[these]/tail(large$sim.deltaR[these],1)
    ## sort of wrong
    large$sim.stdAddRiskVar[these]  <- large$sim.addRiskVar[these]/tail(large$sim.addRiskVar[these],1)
}
large$sim.stdEffect <- 1/sqrt(large$sim.genVar)
## large$sim.mean/large$sim.nSeg
odds1 <- (large$sim.prev + large$sim.deltaR)/(1-large$sim.prev - large$sim.deltaR)
odds0 <- (large$sim.prev)/(1-large$sim.prev)
large$sim.OR <- odds1/odds0
large$sim.logOR  <- log(large$sim.OR)









op <- par(mfrow=c(1,3))

## plot 1
plot(
    NA,
    xlim=c(0,1),
    ylim=c(0,max(large$sim.stdGenVarPerSite)*1.05),
    xlab='Fitness Cost',
    ylab='Relative Per Site Average Heterozygosity'
)
abline(
    h = 1,
    lty = 2,
    lwd = 2
)
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
lines(
    x=large.cost,
    y=1/large.cost,
    lty=3,
    lwd=2
)
legend(
    'topright',
    legend=large.thr,
    col=large.cols,
    pch=20
)

## plot 2

plot(
    NA,
    xlim=c(0,1),
    ylim=c(min(large$sim.scaledDeltaR^2)*0.95,max(large$sim.scaledDeltaR^2)*1.05),
    xlab='Fitness Cost',
    ylab='Additive Risk Var',
    log='y'
)
for(i in 1:length(large.thr)){
    these <- large$thr == large.thr[i]
    points(
        large$cost[these],
        large$sim.scaledDeltaR[these]^2,
        col=large.cols[i],
        pch=20
    )
    lines(
        large$cost[these],
        large$sim.scaledDeltaR[these]^2,
        col=large.cols[i]
    )
    lines(
        large$cost[these],
        large$sim.scaledDeltaR[these]^2,
        col=large.cols[i],
        lty=2
    )
}
abline(
    h=1,
    lty=3,
    lwd=2
)
lines(
    x=large.cost,
    y=1/large.cost,
    lty=2,
    lwd=2
)



## plot 3
plot(
    NA,
    xlim=c(0,1),
    ylim=c(min(large$sim.stdAddRiskVar)*0.95,max(large$sim.stdAddRiskVar)*1.05),
    xlab='Fitness Cost',
    ylab='Scaled Risk Effect',
    log='y'
)
for(i in 1:length(large.thr)){
    these <- large$thr == large.thr[i]
    points(
        large$cost[these],
        large$sim.stdAddRiskVar[these],
        col=large.cols[i],
        pch=20
    )
    lines(
        large$cost[these],
        large$sim.stdAddRiskVar[these],
        col=large.cols[i]
    )
}
abline(
    h=1,
    lty=3,
    lwd=2
)
lines(
    x=large.cost,
    y=1/large.cost,
    lty=2,
    lwd=2
)




















## plot 4
plot(
    NA,
    xlim=c(0,1),
    ylim=c(min(large$sim.deltaR)*0.95,max(large$sim.deltaR)*1.05),
    xlab='Fitness Cost',
    ylab='Risk Effect',
    log='y'
)
for(i in 1:length(large.thr)){
    these <- large$thr == large.thr[i]
    points(
        large$cost[these],
        large$sim.deltaR[these],
        col=large.cols[i],
        pch=20
    )
    lines(
        large$cost[these],
        large$sim.deltaR[these],
        col=large.cols[i]
    )
    lines(
        large$cost[these],
        large$deltaR[these],
        col=large.cols[i],
        lty=2
    )
}






## plot 5
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

## plot 6
plot(
    NA,
    xlim=c(0,1),
    ylim=c(0,5),
    xlab='Fitness Cost',
    ylab='Odds Ratio'
)
for(i in 1:length(large.thr)){
    these <- large$thr == large.thr[i]
    points(
        large$cost[these],
        large$sim.OR[these],
        col=large.cols[i],
        pch=20
    )
    lines(
        large$cost[these],
        large$sim.OR[these],
        col=large.cols[i]
    )
    lines(
        large$cost[these],
        large$sim.OR[these],
        col=large.cols[i],
        lty=2
    )
}











### small effects ###







small.thr <- unique(small$thr)
small$thetaL <- small$theta*small$target.size
small$sim.genVarPerSite <- small$sim.genVar/small$thetaL
small.cols <- viridis(length(small.thr))


small$sim.stdGenVarPerSite  <- NA
for(i in 1:length(small.thr)){
    these <- small$thr == small.thr[i]
    small$sim.stdGenVarPerSite[these]  <- small$sim.genVarPerSite[these]/tail(small$sim.genVarPerSite[these],1)
}





op <- par(mfrow=c(1,2))
plot(
    NA,
    xlim=c(0,1),
    ylim=c(0,max(small$sim.stdGenVarPerSite)*1.05),
    xlab='Fitness Cost',
    ylab='Relative Per Site Average Heterozygosity'
)
abline(h = 1, lty = 3 )
for(i in 1:length(small.thr)){
    these <- small$thr == small.thr[i]
    points(
        small$cost[these],
        small$sim.stdGenVarPerSite[these],
        col=small.cols[i],
        pch=20
    )
    lines(
        small$cost[these],
        small$sim.stdGenVarPerSite[these],
        col=small.cols[i]
    )
}

plot(
    NA,
    xlim=c(0,1),
    ylim=c(min(small$sim.prev)*0.95,max(small$sim.prev)*1.05),
    xlab='Fitness Cost',
    ylab='Prevalence',
    log='y'
)
for(i in 1:length(small.thr)){
    these <- small$thr == small.thr[i]
    points(
        small$cost[these],
        small$sim.prev[these],
        col=small.cols[i],
        pch=20
    )
    lines(
        small$cost[these],
        small$sim.prev[these],
        col=small.cols[i]
    )
    lines(
        small$cost[these],
        small$pois.prev[these],
        col=small.cols[i],
        lty=2
    )
}


