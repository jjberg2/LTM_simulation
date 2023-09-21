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
        SIMS[[i]]$rel.sim.riskFreqLarge <- SIMS[[i]]$sim.riskFreqLarge / SIMS[[i]]$sim.riskFreqLarge[lethal]
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




for(j in 1:length(ds)){
    y.lower = 0
    y.upper = 1/min(ds[[1]][[1]]$C)
    pdf(paste('figures/twoEffectFigures/riskFreqLarge/riskFreqLarge_bt=', sprintf('%.2f',bt[j]), '.pdf', sep = ''),width=20,height=12)
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
        text = 'Relative Mean Frequency of Large Effect Alleles',
        line = 2.7,
        cex = 1.4
    )
    for(i in 1:length(ds[[j]])){
        points(
            x = ds[[j]][[i]]$C,
            y = ds[[j]][[i]]$rel.sim.riskFreqLarge,
            pch = 17,
            col = my.cols[i] ,
            cex = my.pt.cex
        )
    }
    lines(1:9999/10000,10000/(1:9999),lty = 2)
    legend(
        'topright',
        col = my.cols,
        pch = 19,
        legend = round(als[,j],2),
        bty = 'n',
        cex = 1.2
    )
    dev.off()
}


















