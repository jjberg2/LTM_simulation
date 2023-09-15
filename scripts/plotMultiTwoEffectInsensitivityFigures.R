library(viridis)
library(wesanderson)




input <- get(load('twoEffectCostInsensitivityN1000ResultsTable.Rdata'))
d  <- input[input$sim.prev < 1/2 & input$sim.fixedLarge == 0,]



d.tmp <- split(d,round(d$bt,2))
als <- sapply(d.tmp,function(X) unique(X$al))
ds <- lapply(d.tmp,function(X) split(X,X$al))
bt <- as.numeric(names(d.tmp))

## compute relative per site heterozygosities
compute_rels <- function(SIMS){
    for( i in 1:length(SIMS)){
        lethal <- SIMS[[i]]$C==1
        SIMS[[i]]$rel.sim.siteVarSmall <- SIMS[[i]]$sim.siteVarSmall / SIMS[[i]]$sim.siteVarSmall[lethal]
        SIMS[[i]]$rel.sim.siteVarLarge <- SIMS[[i]]$sim.siteVarLarge / SIMS[[i]]$sim.siteVarLarge[lethal]
        SIMS[[i]]$rel.sim.derFreqLarge <- SIMS[[i]]$sim.derFreqLarge / SIMS[[i]]$sim.derFreqLarge[lethal]
        ## compute relative prevalence
        SIMS[[i]]$rel.sim.prev  <- SIMS[[i]]$sim.prev / SIMS[[i]]$sim.prev[lethal]
    }
    SIMS
}
ds <- lapply(ds,compute_rels)
ds.clumped <- do.call(rbind,ds)

min.C <- min(sapply(ds.clumped, function(X) min(X$C)))
max.relprev <- max(sapply(ds.clumped, function(X) max(X$rel.sim.prev)))
max.deltaRLarge <- max(sapply(ds.clumped, function(X) max(X$sim.deltaRLarge)))
max.deltaRSmall <- max(sapply(ds.clumped, function(X) max(X$sim.deltaRSmall)))
max.h2s <- max(sapply(ds.clumped, function(X) max(X$sim.h2s)))
max.h2l <- max(sapply(ds.clumped, function(X) max(X$sim.h2l)))


my.axis.cex = 1.2
my.pt.cex = 1.3
my.seq = seq(0,1,length.out = 1000)
my.cols <- wes_palette('Zissou1',5)
for(j in 1:length(ds)){
    pdf(paste('figures/twoEffectFigures/costInsensitivity/costInsensitivityFigure_bt=' , sprintf('%.2f',bt[j]), '.pdf',sep=''),width=20,height=12)
    par(mfrow = c(1,2))
    plot(
        NA,
        bty = 'n',
        xlim = c(0,1),
        ylim = c(0,1/min.C),
        xlab = '',
        ylab = '',
        cex.lab = 1.4
    )
    mtext(
      side = 1 ,
      text = 'Fitness Cost of Disease',
      line = 2.7,
      cex = 1.4
    )
    mtext(
      side = 2 ,
      text = 'Relative per site contribution to variance',
      line = 2.7,
      cex = 1.4
    )
    text(
        x = 0.84,
        y = 17.2,
        labels = 'Type of site',
        cex = 1.2
    )
    legend(
        x = 0.8,
        y = 17,
        legend = c('Small effect', 'Large effect'),
        pch = c(19,17),
        bty = 'n',
        cex = 1.2
    )
    text(
      x = 0.87,
      y = 13.2,
      labels = 'Large effect size',
      cex = 1.2
    )
    legend(
      x = 0.8,
      y = 13 ,
      col = my.cols,
      pch = 17,
      legend = round(als[,j],2),
      bty = 'n',
      cex = 1.2
    )
    lines(
        x = my.seq,
        y = 1/my.seq,
        lty = 3 ,
        lwd = 1
    )
    abline(
        h = 1 ,
        lty = 2 ,
        lwd = 1
    )
    for(i in 1:length(ds[[j]])){
        points(
            x = ds[[j]][[i]]$C,
            y = ds[[j]][[i]]$rel.sim.siteVarSmall,
            pch = 19,
            col = my.cols[i] ,
            cex = my.pt.cex
        )
        points(
            x = ds[[j]][[i]]$C,
            y = ds[[j]][[i]]$rel.sim.siteVarLarge,
            pch = 17,
            col = my.cols[i] ,
            cex = my.pt.cex
        )
    }
    plot(
        NA,
        ylim = c(0,max.relprev),
        bty = 'n',
        xlim = c(0,1),
        pch = 19,
        xlab = '',
        ylab = '' ,
        cex.lab = 1.4
    )
    mtext(
      side = 1 ,
      text = 'Fitness Cost of Disease',
      line = 2.7,
      cex = 1.4
    )
    mtext(
      side = 2 ,
      text = 'Relative Prevalence',
      line = 2.7,
      cex = 1.4
    )
    for(i in 1:length(ds[[j]])){
        points(
            x = ds[[j]][[i]]$C,
            y = ds[[j]][[i]]$rel.sim.prev,
            pch = 19,
            col = my.cols[i],
            cex = my.pt.cex
        )
    }
    lines(
        x = my.seq,
        y = 1/my.seq,
        lty = 3,
        lwd = 1
    )
    abline(
        h = 1 ,
        lty = 2 ,
        lwd = 1
    )
    dev.off()
}



for(j in 1:length(ds)){
    pdf(paste('figures/twoEffectFigures/riskEffect/riskEffectFigure_bt=', sprintf('%.2f',bt[j]), '.pdf',sep=''),width=20,height=12)
    par(mfrow=c(1,2))
    plot(
        NA,
        ylim = c(0,0.01),
        bty = 'n',
        xlim = c(0,1),
        pch = 19,
        xlab = '',
        ylab = ''
    )
    mtext(
      side = 3 ,
      text = 'Small effect loci',
      line = 1,
      cex = 1.4
    )
    mtext(
      side = 1 ,
      text = 'Fitness Cost of Disease',
      line = 2.7,
      cex = 1.4
    )
    mtext(
      side = 2 ,
      text = 'Effect Size on Risk Scale',
      line = 2.7,
      cex = 1.4
    )
    for(i in 1:length(ds[[j]])){
        points(
            x = ds[[j]][[i]]$C,
            y = ds[[j]][[i]]$sim.deltaRSmall,
            pch = 19,
            col = my.cols[i],
            cex = my.pt.cex
        )
    }
    plot(
        NA,
        ylim = c(0,1),
        bty = 'n',
        xlim = c(0,1),
        pch = 19,
        xlab = '',
        ylab = ''
    )
    mtext(
      side = 3 ,
      text = 'Large effect loci',
      line = 1,
      cex = 1.4
    )
    mtext(
      side = 1 ,
      text = 'Fitness Cost of Disease',
      line = 2.7,
      cex = 1.4
    )
    mtext(
      side = 2 ,
      text = 'Effect Size on Risk Scale',
      line = 2.7,
      cex = 1.4
    )
    text(
      x = 0.93 , 
      y = 1.02 ,
      labels = 'Relative size of' ,
      cex = 1.2
    )
    text(
      x = 0.94 , 
      y = 0.99 ,
      labels = 'large effect sites' ,
      cex = 1.2
    )
    legend(
        x = 0.87,
        y = 0.99,
        col = my.cols,
        pch = 19,
        legend = round(als[,j],2),
        bty = 'n',
        cex = 1.2
    )
    for(i in 1:length(ds[[j]])){
        points(
            x = ds[[j]][[i]]$C,
            y = ds[[j]][[i]]$sim.deltaRLarge,
            pch = 19,
            col = my.cols[i] ,
            cex = my.pt.cex
        )
    }
    dev.off()
}







for(j in 1:length(ds)){
    pdf(paste('figures/twoEffectFigures/selCoef/selCoefFigure_bt=', sprintf('%.2f',bt[j]), '.pdf',sep=''),width=20,height=12)
    par(mfrow=c(1,2))
    plot(
        NA,
        ylim = c(0,0.002),
        bty = 'n',
        xlim = c(0,1),
        pch = 19,
        xlab = '',
        ylab = '' ,
        cex.axis = my.axis.cex
    )
    mtext(
      side = 3 ,
      text = 'Small effect loci',
      line = 1,
      cex = 1.4
    )
    mtext(
      side = 1 ,
      text = 'Fitness Cost of Disease',
      line = 2.7,
      cex = 1.4
    )
    mtext(
      side = 2 ,
      text = 'Product of risk effect and fitness cost  (selection coefficient)',
      line = 2.7,
      cex = 1.4
    )
    for(i in 1:length(ds[[j]])){
        points(
            x = ds[[j]][[i]]$C,
            y = ds[[j]][[i]]$sim.deltaRSmall*ds[[j]][[i]]$C,
            pch = 19,
            col = my.cols[i] ,
            cex = my.pt.cex
        )
    }
    plot(
        NA,
        ylim = c(0,1),
        bty = 'n',
        xlim = c(0,1),
        pch = 19,
        xlab = '',
        ylab = '',
        cex.axis = my.axis.cex
    )
    mtext(
      side = 3 ,
      text = 'Large effect loci',
      line = 1,
      cex = 1.4
    )
    mtext(
      side = 1 ,
      text = 'Fitness Cost of Disease',
      line = 2.7,
      cex = 1.4
    )
    mtext(
      side = 2 ,
      text = 'Product of risk effect and fitness cost (selection coefficient)',
      line = 2.7,
      cex = 1.4
    )
    abline(a=0,b=1,lty=2)
    text(
      x = 0.21 , 
      y = 1.02 ,
      labels = 'Relative size of' ,
      cex = 1.2
    )
    text(
      x = 0.22 , 
      y = 0.99 ,
      labels = 'large effect sites' ,
      cex = 1.2
    )
    legend(
      x = 0.13,
      y = 0.99,
      col = my.cols,
      pch = 19,
      legend = round(als[,j],2),
      bty = 'n',
      cex = 1.2
    )
    for(i in 1:length(ds[[j]])){
        points(
            x = ds[[j]][[i]]$C,
            y = ds[[j]][[i]]$sim.deltaRLarge*ds[[j]][[i]]$C,
            pch = 19,
            col = my.cols[i] ,
            cex = my.pt.cex
        )
    }
    dev.off()
}











for(j in 1:length(ds)){
    pdf(paste('figures/twoEffectFigures/h2Contrib/h2ContribFigure_bt=', sprintf('%.2f',bt[j]), '.pdf', sep = ''),width=20,height=12)
    par(mfrow=c(1,2))
    plot(
        NA,
        ylim = c(0,max.h2s),
        bty = 'n',
        xlim = c(0,1),
        pch = 19,
        xlab = '',
        ylab = '',
        cex.axis = my.axis.cex
    )
    mtext(
      side = 3 ,
      text = 'Small effect loci',
      line = 1,
      cex = 1.4
    )
    mtext(
      side = 1 ,
      text = 'Fitness cost of disease',
      line = 2.7,
      cex = 1.4
    )
    mtext(
      side = 2 ,
      text = 'Contribution to heritability',
      line = 2.7,
      cex = 1.4
    )
    for(i in 1:length(ds[[j]])){
        points(
            x = ds[[j]][[i]]$C,
            y = ds[[j]][[i]]$sim.h2s,
            pch = 19,
            col = my.cols[i] ,
            cex = my.pt.cex
        )
    }
    plot(
        NA,
        ylim = c(0,max.h2l),
        bty = 'n',
        xlim = c(0,1),
        pch = 19,
        xlab = '',
        ylab = ''
    )
    mtext(
      side = 3 ,
      text = 'Large effect loci',
      line = 1,
      cex = 1.4
    )
    mtext(
      side = 1 ,
      text = 'Fitness cost of disease',
      line = 2.7,
      cex = 1.4
    )
    mtext(
      side = 2 ,
      text = 'Contribution to heritability',
      line = 2.7,
      cex = 1.4
    )
    text(
      x = 0.93 , 
      y = 0.057 ,
      labels = 'Relative size of' ,
      cex = 1.2
    )
    text(
      x = 0.94 , 
      y = 0.0555 ,
      labels = 'large effect sites' ,
      cex = 1.2
    )
    legend(
      x = 0.87,
      y = 0.055,
      col = my.cols,
      pch = 19,
      legend = round(als[,j],2),
      bty = 'n',
      cex = 1.2
    )
    for(i in 1:length(ds[[j]])){
        points(
            x = ds[[j]][[i]]$C,
            y = ds[[j]][[i]]$sim.h2l,
            pch = 19,
            col = my.cols[i],
            cex = my.pt.cex
        )
    }
    dev.off()
}





for(j in 1:length(ds)){
  my.fixed = unlist(sapply(ds[[j]], function(X) X$sim.fixedSmall))
  my.means = unlist(sapply(ds[[j]],function(X) X$sim.mean))
  y.lower = min(my.means,my.fixed) - (min(my.means,my.fixed) %% 100)
  y.upper = max(my.means,my.fixed) + (100-(max(my.means,my.fixed) %% 100)) + 100
  pdf(paste('figures/twoEffectFigures/means/meansFigure_bt=', sprintf('%.2f',bt[j]), '.pdf', sep = ''),width=20,height=12)
  ## par(mfrow=c(1,2))
  plot(
    NA,
    ylim = c(y.lower,y.upper),
    bty = 'n',
    xlim = c(0,1),
    pch = 19,
    xlab = '',
    ylab = '',
    cex.axis = my.axis.cex
  )
  mtext(
    side = 1 ,
    text = 'Fitness cost of disease',
    line = 2.7,
    cex = 1.4
  )
  mtext(
    side = 2 ,
    text = 'Liability',
    line = 2.7,
    cex = 1.4
  )
  ##for(i in 1:length(ds[[j]])){
  i=3
  points(
    x = ds[[j]][[i]]$C,
    y = ds[[j]][[i]]$sim.fixedSmall,
    pch = 17,
    col = my.cols[i] ,
    cex = my.pt.cex
  )
  points(
    x = ds[[j]][[i]]$C,
    y = ds[[j]][[i]]$sim.mean,
    pch = 19,
    col = my.cols[i] ,
    cex = my.pt.cex
  )
  abline(h = ds[[j]][[i]]$thr[1],lty=2)
  legend(
   'topright',
    col = c(rep(my.cols[i],2),'black'),
    pch = c(17,19,-1),
    lty = c(0,0,2),
    legend = c('Fixed','Mean','Threshold'),
    bty = 'n',
    cex = 1.2
  )
  dev.off()
}






if(FALSE){







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


}
