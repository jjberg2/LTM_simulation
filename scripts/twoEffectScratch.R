re.sim <- TRUE
if(re.sim) source('scripts/twoEffectPrevalenceSolutions.R')
solns <- get(load(file='solutions/twoEffectPrevalenceSolutions_solutions.Robj'))
output <- get(load(file='solutions/twoEffectPrevalenceSolutions_output.Robj'))


my.prevs <- list()
my.multi.norm.prevs <- list()
my.deltass <- list()
my.deltals <- list()
my.ys <- list()
my.std.as <- list()
my.std.al <- list()
my.pgl <- list()
my.tstar <- list()
my.multi.norm.tstar <- list()
for (l in seq_along(my.bt)) {
  my.prevs[[l]] <- list()
  my.multi.norm.prevs[[l]] <- list()
  my.deltass[[l]] <- list()
  my.deltals[[l]] <- list()
  my.multi.norm.tstar[[l]] <- list()
  my.std.as[[l]] <- list()
  my.ys[[l]] <- list()
  my.std.al[[l]] <- list()
  my.pgl[[l]] <- list()
  my.tstar[[l]] <- list()
  for (k in seq_along(h2)) {
    my.prevs[[l]][[k]] <- sapply(output[[l]][[k]], function(X)
      X$prev)
    my.multi.norm.prevs[[l]][[k]] <-
      sapply(output[[l]][[k]], function(X)
        X$multi.norm.prev)
    my.std.as[[l]][[k]] <- sapply(output[[l]][[k]], function(X)
      X$std.as)
    my.ys[[l]][[k]] <- sapply(output[[l]][[k]], function(X)
      X$ys)
    my.std.al[[l]][[k]] <- sapply(output[[l]][[k]], function(X)
      X$std.al)
    my.deltals[[l]][[k]] <- sapply(output[[l]][[k]], function(X)
      X$deltal)
    my.deltass[[l]][[k]] <- sapply(output[[l]][[k]], function(X)
      X$deltas)
    my.pgl[[l]][[k]] <- sapply(output[[l]][[k]], function(X)
      X$pgl)
    my.tstar[[l]][[k]] <- sapply(output[[l]][[k]], function(X)
      X$tstar)
    my.multi.norm.tstar[[l]][[k]] <-
      sapply(output[[l]][[k]], function(X)
        X$multi.norm.tstar)
  }
}





my.cols = met.brewer('Isfahan2', 2)
my.pch = c(21, 22)
xmax <- 0.04
ymax <- 0.02
plot(
  NA,
  bty = 'n',
  ylim = c(0,ymax),
  xlim = c(0,xmax)
)
for(ll in seq_along(output)) {
  for (kk in seq_along(output[[ll]])) {
    std.fts <- sapply(output[[ll]][[kk]], function(X)
      X$std.ft)
    prevs <- sapply(output[[ll]][[kk]], function(X)
      X$prev)
    matplot(
      std.fts,
      prevs,
      col = my.cols,
      lty = kk,
      lwd = 2,
      type = 'l',
      add = TRUE
    )
    ft.range <- seq(1e-9,xmax,length.out=1000)
  }
}
lines(
  x = ft.range,
  y = 1-pnorm(dnorminv(ft.range))
)



these.bts <- sapply(output,function(X) unique(round(X[[1]][[1]]$bt,2)))
these.var.ratios <- sapply(output[[1]][[1]],function(X) unique(round(X$var.ratio,2)) )
these.h2 <- sapply(output[[1]],function(X) unique(round(X[[1]]$h2,2)))
layout(matrix(c(1,2),ncol=2),width=c(5,1))
par(mar=c(5,4,4,1)+0.1)
my.cols = met.brewer('Isfahan2', 2)
my.pch = c(21, 22)
xmax <- 1
ymin <- 0
ymax <- 1
plot(
  NA,
  bty = 'n',
  ylim = c(ymin,ymax),
  xlim = c(0,xmax)
)

for(ll in seq_along(output)) {
  for (kk in seq_along(output[[ll]])) {
    
    deltals <- sapply(output[[ll]][[kk]], function(X)
      X$deltal)
    h2s.est <- sapply(output[[ll]][[kk]], function(X)
      X$h2s.est)
    h2s <- sapply(output[[ll]][[kk]], function(X)
      X$h2s)
    matplot(
      deltals,
      h2s.est,
      col = my.cols,
      lty = kk,
      lwd = 2,
      type = 'l',
      add = TRUE
    )
    matplot(
      deltals,
      h2s,
      col = my.cols,
      lty = 3,
      lwd = 1/2,
      type = 'l',
      add = TRUE
    )
  }
}
par(mar=c(0,0,0,0))
plot.new()
legend( 'top', bty='n', lty=rep(1,2), leg=these.var.ratios, cex=1.4, col=my.cols, lwd=2 )
legend( 'center', bty='n', lty=c(1,2), leg=these.h2, cex=1.4, lwd=2 )




these.bts <- sapply(output,function(X) unique(round(X[[1]][[1]]$bt,2)))
these.var.ratios <- sapply(output[[1]][[1]],function(X) unique(round(X$var.ratio,2)) )
these.h2 <- sapply(output[[1]],function(X) unique(round(X[[1]]$h2,2)))
layout(matrix(c(1,2),ncol=2),width=c(5,1))
par(mar=c(5,4,4,1)+0.1)
my.cols = met.brewer('Isfahan2', length(these.var.ratios))
my.pch = c(21, 22)
xmax <- 1
ymin <- 0
ymax <- 2
plot(
  NA,
  bty = 'n',
  ylim = c(ymin,ymax),
  xlim = c(0,xmax)
)
for(ll in seq_along(output)) {
  if(ll>1) break
  for (kk in seq_along(output[[ll]])) {
    deltals <- sapply(output[[ll]][[kk]], function(X)
      X$deltal)
    h2s.est <- sapply(output[[ll]][[kk]], function(X)
      X$h2s.est)
    h2s <- sapply(output[[ll]][[kk]], function(X)
      X$h2s)
    h2all.est <- sapply(output[[ll]][[kk]], function(X)
      X$h2all.est)
    h2 <- sapply(output[[ll]][[kk]], function(X)
      X$h2)
    matplot(
      deltals,
      h2s.est/h2s,
      col = my.cols,
      lty = kk,
      lwd = 2,
      type = 'l',
      add = TRUE
    )
    matplot(
      deltals,
      h2all.est/h2,
      col = my.cols,
      lty = kk,
      lwd = 2,
      type = 'l',
      add = TRUE
    )
    # matplot(
    #   deltals,
    #   h2s,
    #   col = my.cols,
    #   lty = 3,
    #   lwd = 1/2,
    #   type = 'l',
    #   add = TRUE
    # )
  }
}
par(mar=c(0,0,0,0))
plot.new()
legend( 'top', bty='n', lty=rep(1,2), leg=these.var.ratios, cex=1.4, col=my.cols, lwd=2 )
legend( 'center', bty='n', lty=c(1,2), leg=these.h2, cex=1.4, lwd=2 )






pdf(
  'figures/suppFigures/twoEffectStdas.pdf',
  width = 12,
  height = 6
)
par(mfrow = c(1, 2))
for (ll in seq_along(my.bt)) {
  
  norm.ft <-
    1 / (4 * Ne * cost) * log((1 + my.bt[ll] / (1 - my.bt[ll])))
  norm.Vg <-
    8 * Ne * L * u * my.bt[ll] / log((1 + my.bt[ll] / (1 - my.bt[ll])))
  norm.Vt <- norm.Vg / h2
  norm.astd <- 1 / sqrt (norm.Vt)
  
  plot(
    NA,
    xlim = range(c(0, my.als)),
    ylim = c(0, max(unlist(my.std.as[[ll]]))),
    ylab = "Standardized effect (small effect loci)",
    xlab = "Ratio of large to small effect size",
    main = paste('b_T = ', my.bt[ll], sep = '')
  )
  
  # if (ll == 1) {
  #   legend(
  #     'bottomright',
  #     col = my.cols,
  #     lty = 1 ,
  #     lwd = 2,
  #     legend = var.ratio,
  #     bty = 'n'
  #   )
  #   text(x = 195,
  #        y = 0.00047,
  #        labels = 'Variance')
  #   text(x = 195,
  #        y = 0.00035,
  #        labels = 'Ratio')
  #   legend(
  #     x = 120,
  #     y = 0.0003,
  #     lty = 1:2 ,
  #     lwd = 2,
  #     legend = round(h2, 2),
  #     bty = 'n'
  #   )
  #   legend(
  #     x = 114,
  #     y = 0.0003,
  #     pch = my.pch, 
  #     legend = c('',''),
  #     bty = 'n'
  #   )
  #   text(x = 135,
  #        y = 0.00032,
  #        labels = 'Heritability')
  # }
  for (kk in seq_along(h2)) {
    for (jj in seq_along(var.ratio)) {
      plot.als <- output[[ll]][[kk]][[jj]][, 'al']
      lines(
        plot.als,
        my.std.as[[ll]][[kk]][, jj],
        col = my.cols[jj],
        lty = kk,
        lwd = 2
      )
      points(x = 1,
             y = norm.astd[[kk]],
             pch = my.pch[[kk]])
    }
  }
}
dev.off()


pdf(
  'figures/suppFigures/twoEffectStdal.pdf',
  width = 12,
  height = 6
)
par(mfrow = c(1, 2))
for (ll in seq_along(my.bt)) {
  
  plot(
    NA,
    xlim = range(c(0, my.als)),
    ylim = c(0, max(unlist(my.std.al[[ll]]))),
    ylab = "Standardized effect (large effect loci)",
    xlab = "Ratio of large to small effect size",
    main = paste('b_T = ', my.bt[ll], sep = '')
  )
  for (kk in seq_along(h2)) {
    for (jj in seq_along(var.ratio)) {
      plot.als <- output[[ll]][[kk]][[jj]][, 'al']
      lines(
        plot.als,
        my.std.al[[ll]][[kk]][, jj],
        col = my.cols[jj],
        lty = kk,
        lwd = 2
      )
    }
  }
}
dev.off()



pdf(
  'figures/suppFigures/twoEffectDeltass.pdf',
  width = 12,
  height = 6
)
par(mfrow = c(1, 2))
for (ll in seq_along(my.bt)) {
  plot(
    NA,
    xlim = range(c(0, my.als)),
    ylim = c(0, max(unlist(my.deltass[[ll]]))),
    ylab = "Risk Effect (small effect loci)",
    xlab = "Ratio of large to small effect size",
    main = paste('b_T = ', my.bt[ll], sep = '')
  )
  for (kk in seq_along(h2)) {
    for (jj in seq_along(var.ratio)) {
      plot.als <- output[[ll]][[kk]][[jj]][, 'al']
      lines(
        plot.als,
        my.deltass[[ll]][[kk]][, jj],
        col = my.cols[jj],
        lty = kk,
        lwd = 2
      )
    }
  }
}
dev.off()



pdf(
  'figures/suppFigures/twoEffectDeltals.pdf',
  width = 12,
  height = 6
)
par(mfrow = c(1, 2))
for (ll in seq_along(my.bt)) {
  
  plot(
    NA,
    xlim = range(c(0, my.als)),
    ylim = c(0, max(unlist(my.deltals[[ll]]))),
    ylab = "Standardized effect (large effect loci)",
    xlab = "Ratio of large to small effect size",
    main = paste('b_T = ', my.bt[ll], sep = '')
  )
  for (kk in seq_along(h2)) {
    for (jj in seq_along(var.ratio)) {
      plot.als <- output[[ll]][[kk]][[jj]][, 'al']
      lines(
        plot.als,
        my.deltals[[ll]][[kk]][, jj],
        col = my.cols[jj],
        lty = kk,
        lwd = 2
      )
    }
  }
}
dev.off()


pdf(
  'figures/suppFigures/twoEffectys.pdf',
  width = 12,
  height = 6
)
par(mfrow = c(1, 2))
for (ll in seq_along(my.bt)) {
  plot(
    NA,
    xlim = range(c(0, my.als)),
    ylim = c(0, max(unlist(my.ys[[ll]]))),
    ylab = "Scaled selection coefficient (small effect loci)",
    xlab = "Ratio of large to small effect size",
    main = paste('b_T = ', my.bt[ll], sep = '')
  )
  for (kk in seq_along(h2)) {
    for (jj in seq_along(var.ratio)) {
      plot.als <- output[[ll]][[kk]][[jj]][, 'al']
      lines(
        plot.als,
        my.ys[[ll]][[kk]][, jj],
        col = my.cols[jj],
        lty = kk,
        lwd = 2
      )
    }
  }
}
dev.off()


pdf(
  'figures/suppFigures/twoEffectpgl.pdf',
  width = 12,
  height = 6
)
par(mfrow = c(1, 2))
for (ll in seq_along(my.bt)) {
  plot(
    NA,
    xlim = range(c(0, my.als)),
    ylim = c(0, max(unlist(my.pgl[[ll]]))),
    ylab = "Scaled selection coefficient (small effect loci)",
    xlab = "Ratio of large to small effect size",
    main = paste('b_T = ', my.bt[ll], sep = '')
  )
  if (ll == 1) {
    legend(
      'bottomright',
      col = my.cols,
      lty = 1 ,
      lwd = 2,
      legend = var.ratio,
      bty = 'n'
    )
    text(x = 195,
         y = 0.00047,
         labels = 'Variance')
    text(x = 195,
         y = 0.00035,
         labels = 'Ratio')
  }
  for (kk in seq_along(h2)) {
    for (jj in seq_along(var.ratio)) {
      plot.als <- output[[ll]][[kk]][[jj]][, 'al']
      lines(
        plot.als,
        my.pgl[[ll]][[kk]][, jj],
        col = my.cols[jj],
        lty = kk,
        lwd = 2
      )
    }
  }
}
dev.off()


pdf(
  'figures/suppFigures/twoEffectTStar.pdf',
  width = 12,
  height = 6
)
par(mfrow = c(1, 2))
for (ll in seq_along(my.bt)) {
  plot(
    NA,
    xlim = range(c(0, my.als)),
    ylim = c(0, max(unlist(my.tstar[[ll]]))),
    ylab = "Thresold distance in standardized units",
    xlab = "Ratio of large to small effect size",
    main = paste('b_T = ', my.bt[ll], sep = '')
  )
  if (ll == 1) {
    legend(
      x = 335,
      y = 0.6,
      col = my.cols,
      lty = 1 ,
      lwd = 2,
      legend = var.ratio,
      bty = 'n'
    )
    text(x = 360,
         y = 0.62,
         labels = 'Variance')
    text(x = 360,
         y = 0.6,
         labels = 'Ratio')
    legend(
      x = 335,
      y = 0.48,
      lty = 1:2 ,
      lwd = 2,
      legend = round(h2, 2),
      bty = 'n'
    )
    legend(
      x = 335,
      y = 0.48,
      pch = my.pch, 
      legend = c('',''),
      bty = 'n'
    )
    text(x = 360,
         y = 0.48,
         labels = 'Heritability')
  }
  for (kk in seq_along(h2)) {
    for (jj in seq_along(var.ratio)) {
      plot.als <- output[[ll]][[kk]][[jj]][, 'al']
      lines(
        plot.als,
        my.tstar[[ll]][[kk]][, jj],
        col = my.cols[jj],
        lty = kk,
        lwd = 2
      )
      lines(
        plot.als,
        my.multi.norm.tstar[[ll]][[kk]][, jj],
        col = my.cols[jj],
        lty = 3,
        lwd = 2/3
      )
    }
  }
}
dev.off()


ymax <- 4.5
xmax <- 0.04
plot(
  NA,
  bty = 'n',
  ylim = c(0,ymax),
  xlim = c(0,xmax),
  ylab = 'Standardized Threshold Distance',
  xlab = 'Standardized Threshold Density'
)
for(ll in seq_along(output)) {
  for (kk in seq_along(output[[ll]])) {
    std.fts <- sapply(output[[ll]][[kk]], function(X)
      X$std.ft)
    tstars <- sapply(output[[ll]][[kk]], function(X)
      X$tstar)
    matplot(
      std.fts,
      tstars,
      col = my.cols,
      lty = 1,
      lwd = 2,
      type = 'l',
      add = TRUE
    )
    ft.range <- seq(1e-9,xmax,length.out=1000)
  }
}
lines(
  x = ft.range,
  y = dnorminv(ft.range)
)


my.cols = met.brewer('Isfahan2', 3)
ymax <- 0.02
xmax <- 0.04
plot(
  NA,
  bty = 'n',
  ylim = c(0,ymax),
  xlim = c(0,xmax),
  ylab = 'Prevalence',
  xlab = 'Standardized Threshold Density'
)
for(ll in seq_along(output)) {
  for (kk in seq_along(output[[ll]])) {
    std.fts <- sapply(output[[ll]][[kk]], function(X)
      X$std.ft)
    prevs <- sapply(output[[ll]][[kk]], function(X)
      X$prev)
    matplot(
      std.fts,
      prevs,
      col = my.cols,
      lty = kk,
      lwd = 2,
      type = 'l',
      add = TRUE
    )
    ft.range <- seq(1e-9,xmax,length.out=1000)
  }
}
lines(
  x = ft.range,
  y = 1-pnorm(dnorminv(ft.range))
)



ymax <- 8
xmax <- 30
plot(
  NA,
  bty = 'n',
  ylim = c(0,ymax),
  xlim = c(0,xmax),
  ylab = 'Standardized Threshold Distance',
  xlab = 'Scaled large selection coefficient'
)
for(ll in seq_along(output)) {
  for (kk in seq_along(output[[ll]])) {
    deltals <- sapply(output[[ll]][[kk]], function(X)
      X$deltal)
    tstars <- sapply(output[[ll]][[kk]], function(X)
      X$tstar)
    matplot(
      deltals*Ne*cost,
      tstars,
      col = my.cols,
      lty = 1,
      lwd = 2,
      type = 'l',
      add = TRUE
    )
    ft.range <- seq(1e-9,xmax,length.out=1000)
  }
}


ymax <- 1
xmax <- 0.2
plot(
  NA,
  bty = 'n',
  ylim = c(0,ymax),
  xlim = c(0,xmax),
  ylab = 'Proportion of liability variance',
  xlab = 'PAR of large effect alleles'
)
for(ll in seq_along(output)) {
  for (kk in seq_along(output[[ll]])) {
    deltals <- sapply(output[[ll]][[kk]], function(X)
      X$deltal)
    pgls <- sapply(output[[ll]][[kk]], function(X)
      X$pgl)
    matplot(
      deltals,
      pgls,
      col = my.cols,
      lty = kk,
      lwd = 2,
      type = 'l',
      add = TRUE
    )
    ft.range <- seq(1e-9,xmax,length.out=1000)
  }
}


ymax <- 0.007 
xmax <- 0.02
plot(
  NA,
  bty = 'n',
  ylim = c(0,ymax),
  xlim = c(0,xmax),
  ylab = 'Standardized Threshold Distance',
  xlab = 'Large effect size coefficient'
)
for(ll in seq_along(output)) {
  for (kk in seq_along(output[[ll]])) {
    deltals <- sapply(output[[ll]][[kk]], function(X)
      X$deltal)
    prevs <- sapply(output[[ll]][[kk]], function(X)
      X$prev)
    matplot(
      deltals,
      prevs,
      col = my.cols,
      lty = kk,
      lwd = 2,
      type = 'l',
      add = TRUE
    )
    ft.range <- seq(1e-9,xmax,length.out=1000)
  }
}



ymax <- 30
xmax <- 0.02
plot(
  NA,
  bty = 'n',
  ylim = c(0,ymax),
  xlim = c(0,xmax),
  ylab = 'Mean number of large effect alleles per individual',
  xlab = 'PAR of large effect sites'
)
for(ll in seq_along(output)) {
  for (kk in seq_along(output[[ll]])) {
    deltals <- sapply(output[[ll]][[kk]], function(X)
      X$deltal)
    mean.nls <- sapply(output[[ll]][[kk]], function(X)
      X$mean.nl)
    matplot(
      deltals,
      mean.nls,
      col = my.cols,
      lty = kk,
      lwd = 2,
      type = 'l',
      add = TRUE
    )
  }
}

ymax <- 0.05
xmax <- 0.02
plot(
  NA,
  bty = 'n',
  ylim = c(0,ymax),
  xlim = c(0,xmax),
  ylab = 'Mean number of new large effect mutations per individual',
  xlab = 'PAR of large effect sites'
)
for(ll in seq_along(output)) {
  for (kk in seq_along(output[[ll]])) {
    deltals <- sapply(output[[ll]][[kk]], function(X)
      X$deltal)
    bigU_large <- sapply(output[[ll]][[kk]], function(X)
      X$bigU_large)
    matplot(
      deltals,
      bigU_large,
      col = my.cols,
      lty = kk,
      lwd = 2,
      type = 'l',
      add = TRUE
    )
  }
}