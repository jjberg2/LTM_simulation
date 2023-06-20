library(viridis)



input <- get(load('smallEffectInsensitivityResultsTable.Rdata'))
small  <- input[input$sim.prev < 1/2,]
small.thr <- unique(small$thr)
small.gamma <- unique(small$gamma)
small.b <- unique(small$b)
small$thetaL <- small$theta*small$target.size
## small$sim.genVarPerSite <- small$sim.genVar/small$thetaL

small$sim.addRiskHet <- small$sim.siteVar*small$sim.deltaR

small$sim.stdGenVarPerSite  <- NA
small$sim.scaledDeltaR  <- NA
small$sim.stdAddRiskHet <- NA
for(i in 1:length(small.thr)){
    these <- small$thr == small.thr[i] & small$sim.prev < 1/2
    small$sim.stdGenVarPerSite[these]  <- small$sim.siteVar[these]/tail(small$sim.siteVar[these],1)
    small$sim.scaledDeltaR[these]  <- small$sim.deltaR[these]/tail(small$sim.deltaR[these],1)
    small$sim.stdAddRiskHet[these]  <- small$sim.addRiskHet[these]/tail(small$sim.addRiskHet[these],1)
}
small.cols <- seq_along(small.thr)







small$genVarPerSite <- (2*small$b*small$theta*small$target.size/small$gamma)/small$target.size
small$sim.genVarPerSite <- small$sim.genVar/small$target.size

plot(
    NA,
    xlim=c(0,1),
    ylim=c(0,0.0021),
    bty='n',
    xlab='Fitness Cost',
    ylab='Genetic Variance Per Site'
)
for(i in 1:length(small.b)){
    these <- small$b==small.b[i] & small$sim.prev < 1/2
    lines(
        small[these,'cost'],
        small[these,'genVarPerSite'],
        col=small.cols[i]
    )
    points(
        small[these,'cost'],
        small[these,'sim.genVarPerSite'],
        pch=20,
        col=small.cols[i]
    )
}
legend(
    'bottomleft',
    legend=small.gamma,
    pch=20,
    col=small.cols
)


smallAvgs <- data.frame(b= small.b)
for(i in 1:length(small.b)){
    these <- small$b==small.b[i] & small$sim.prev < 1/2
    smallAvgs$genVar[i] <- mean(small$genVar[these])
    smallAvgs$sim.genVar[i] <- mean(small$sim.genVar[these])    
}


















large <- get(load('largeEffectInsensitivityResultsTableBackup.Rdata'))
large.cost <- sort(unique(large$cost))
large.thr <- unique(large$thr)
large$thetaL <- large$theta*large$target.size
large$sim.genVarPerSite <- large$sim.genVar/large$thetaL
large.cols <- viridis(length(large.thr))


## this is sort of wrong but close
large$sim.addRiskHet <- large$sim.genVarPerSite*large$sim.deltaR
large$pois.Var  <- large$mean
large$sim.stdGenVarPerSite  <- NA
large$sim.scaledDeltaR  <- NA
large$sim.stdAddRiskHet <- NA
for(i in 1:length(large.thr)){
    these <- large$thr == large.thr[i]
    large$sim.stdGenVarPerSite[these]  <- large$sim.genVarPerSite[these]/tail(large$sim.genVarPerSite[these],1)
    large$sim.scaledDeltaR[these]  <- large$sim.deltaR[these]/tail(large$sim.deltaR[these],1)
    ## sort of wrong
    large$sim.stdAddRiskHet[these]  <- large$sim.addRiskHet[these]/tail(large$sim.addRiskHet[these],1)
}
large$sim.stdEffect <- 1/sqrt(large$sim.genVar)
## large$sim.mean/large$sim.nSeg
odds1 <- (large$sim.prev + large$sim.deltaR)/(1-large$sim.prev - large$sim.deltaR)
odds0 <- (large$sim.prev)/(1-large$sim.prev)
large$sim.OR <- odds1/odds0
large$sim.logOR  <- log(large$sim.OR)








pdf('figures/costInsensitivityFigure.pdf',width=20,height=12)

op <- par(mfrow=c(1,3))
## plot 1
plot(
    NA,
    xlim=c(0,1),
    ylim=c(min(large$sim.stdGenVarPerSite)*0.95,11),
    xlab='Fitness Cost',
    ylab='Relative Per Site Average Heterozygosity',
    yaxt='n'
)
axis(
    side=2,
    at=seq(1,12)
)
abline(
    h = 1,
    lty = 2,
    lwd = 2
)
lines(
    x=large.cost,
    y=1/large.cost,
    lty=3,
    lwd=2
)
legend(
    'topright',
    legend=round(large.thr,0),
    col=large.cols,
    pch=20,
    bty='n'
)
legend(
    x=0.95,
    y=10,
    legend=small.b,
    col=small.cols,
    pch=20,
    bty='n'
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
for(i in 1:length(small.thr)){
    these.small <- small$thr == small.thr[i]
    points(
        small$cost[these.small],
        small$sim.stdGenVarPerSite[these.small],
        col=small.cols[i],
        pch=20
    )
    lines(
        small$cost[these.small],
        small$sim.stdGenVarPerSite[these.small],
        col=small.cols[i]
    )
}


## plot 2

plot(
    NA,
    xlim=c(0,1),
    ylim=c(min(large$sim.scaledDeltaR)*0.95,max(large$sim.scaledDeltaR)*1.05),
    xlab='Fitness Cost',
    ylab='Scaled Risk Effect',
    log='y'
)
for(i in 1:length(large.thr)){
    these <- large$thr == large.thr[i]
    points(
        large$cost[these],
        large$sim.scaledDeltaR[these],
        col=large.cols[i],
        pch=20
    )
    lines(
        large$cost[these],
        large$sim.scaledDeltaR[these],
        col=large.cols[i]
    )
   
}
for(i in 1:length(small.thr)){
    these.small <- small$thr == small.thr[i]
    points(
        small$cost[these.small],
        small$sim.scaledDeltaR[these.small],
        col=small.cols[i],
        pch=20
    )
    lines(
        small$cost[these.small],
        small$sim.scaledDeltaR[these.small],
        col=small.cols[i]
    )
}
abline(
    h=1,
    lty=3,
    lwd=2
)
fine.costs <- seq(0.001,1,length.out=1000)
lines(
    x=fine.costs,
    y=1/fine.costs,
    lty=2,
    lwd=2
)



## plot 3
plot(
    NA,
    xlim=c(0,1),
    ylim=c(min(large$sim.stdAddRiskHet)*0.95,max(large$sim.stdAddRiskHet)*1.05),
    xlab='Fitness Cost',
    ylab='Additive Risk Var',
    log='y'
)
## for(i in 1:length(large.thr)){
##     these <- large$thr == large.thr[i]
##     points(
##         large$cost[these],
##         large$sim.stdAddRiskHet[these],
##         col=large.cols[i],
##         pch=20
##     )
##     lines(
##         large$cost[these],
##         large$sim.stdAddRiskHet[these],
##         col=large.cols[i]
##     )
## }
for(i in 1:length(small.thr)){
    these.small <- small$thr == small.thr[i]
    points(
        small$cost[these.small],
        small$sim.stdAddRiskHet[these.small],
        col=small.cols[i],
        pch=20
    )
    lines(
        small$cost[these.small],
        small$sim.stdAddRiskHet[these.small],
        col=small.cols[i]
    )
}
abline(
    h=1,
    lty=3,
    lwd=2
)
lines(
    x=fine.costs,
    y=1/fine.costs,
    lty=2,
    lwd=3
)

dev.off()



















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

legend(
    'topright',
    legend=large.thr,
    col=large.cols,
    pch=20
)



## plot 6
plot(
    NA,
    xlim=c(0,1),
    ylim=c(0,50),
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


