sim.new <- T
if (sim.new) {
  rm(list = ls())
  library('MetBrewer')
  source('scripts/solveTwoEffect.R')
  recover.flag <- F
  
  
  
  ## main single effect params
  L <- 1e7
  Ne <- 3e3
  u <- 1e-8
  cost <- 1 / 2 ##c(1,1/2,1,1/2)
  my.bt <- seq(0.3, 0.9, length.out = 4)  ##c(0.6,0.6,0.8,0.8)
  h2 <- c(1 / 6, 5 / 6)
  
  ## ratio of large effect observed var to small effect
  var.ratio <- c(1 / 2, 1, 2)
  
  as <- 1
  my.als <- exp(seq(log(10), log(400), length.out = 60))
  theta <- 4 * Ne * u
  
  
  tmp.output <- list()
  nulls <- list()
  output <- list()
  my.prevs <- list()
  my.multi.norm.prevs <- list()
  my.deltass <- list()
  my.deltals <- list()
  solns <- list()
  my.ys <- list()
  my.std.as <- list()
  my.std.al <- list()
  my.pgl <- list()
  my.tstar <- list()
  my.multi.norm.tstar <- list()
  for (l in seq_along(my.bt)) {
    output[[l]] <- list()
    my.prevs[[l]] <- list()
    my.multi.norm.prevs[[l]] <- list()
    my.deltass[[l]] <- list()
    my.deltals[[l]] <- list()
    solns[[l]] <- list()
    my.std.as[[l]] <- list()
    my.ys[[l]] <- list()
    my.std.al[[l]] <- list()
    my.pgl[[l]] <- list()
    my.tstar[[l]] <- list()
    my.multi.norm.tstar[[l]] <- list()
    for (k in seq_along(h2)) {
      output[[l]][[k]] <- list()
      solns[[l]][[k]] <- list()
      for (j in seq_along(var.ratio)) {
        last.gs <- numeric()
        last.gs[1] <- 0.99
        last.bs <- numeric()
        last.bs[1] <- my.bt[l] - 0.05
        solns[[l]][[k]][[j]] <- list()
        for (i in seq_along(my.als)) {
          ## if(l==1 && k==3 && j==1 && i==1) stop('Stop')
          solns[[l]][[k]][[j]][[i]] <- solveTwoEffect(
            bs = last.bs[i],
            bt = my.bt[l],
            Ne = Ne,
            as = as,
            al = my.als[i],
            L = L,
            gs = last.gs[i],
            last.tstar = NULL,
            h2 = h2[k],
            u = u,
            C = cost,
            LL.soln = TRUE,
            var.ratio = var.ratio[j],
            equalize.observed.vars = TRUE
          )
          # if (length(solns[[i]]) == 1) {
          #   nulls <- append(nulls, list(c(l, k, j, i)))
          #   last.gs[i + 1] <- last.gs[i]
          #   last.bs[i + 1] <- last.bs[i]
          #   stop('stop')
          #   next
          # } else{
          # }
          tmp.output[[i]] <- makeOutput(solns[[l]][[k]][[j]][[i]])
          last.gs[i + 1]  <- tmp.output[[i]]['gs']
          last.bs[i + 1] <- tmp.output[[i]]['bs']
          if (i %% 10 == 0)
            print(paste(l, k, j, i, sep = ':'))
        }
        output[[l]][[k]][[j]] <-
          as.data.frame(do.call(rbind, tmp.output))
      }
      my.prevs[[l]][[k]] <- sapply(output[[l]][[k]], function(X)
        X$prev)
      my.multi.norm.prevs[[l]][[k]] <- sapply(output[[l]][[k]], function(X)
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
      my.multi.norm.tstar[[l]][[k]] <- sapply(output[[l]][[k]], function(X)
        X$multi.norm.tstar)
    }
  }
}
  
save(solns,file='solutions/twoEffectPrevalenceSolutions_solutions.Robj')
save(output,file='solutions/twoEffectPrevalenceSolutions_output.Robj')





my.cols = met.brewer('Isfahan2', 3)
my.pch = c(21, 22)

pdf(
  'figures/suppFigures/twoEffectPrevalence.pdf',
  width = 12,
  height = 12
)
par(mfrow = c(2, 2))
for (ll in seq_along(my.bt)) {
  norm.ft <-
    1 / (4 * Ne * cost) * log((1 + my.bt[ll] / (1 - my.bt[ll])))
  norm.Vg <-
    8 * Ne * L * u * my.bt[ll] / log((1 + my.bt[ll] / (1 - my.bt[ll])))
  norm.Vt <- norm.Vg / h2
  norm.astd <- 1 / sqrt (norm.Vt)
  norm.dens <- 2 * L * u * my.bt[ll] * norm.astd / (h2 * cost)
  norm.prev <- 1 - pnorm(dnorminv(norm.dens))
  
  
  plot(
    NA,
    xlim = range(c(0, my.als)),
    ylim = c(0, max(unlist(my.prevs[[ll]]))),
    ylab = "Disease Prevalence",
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
    text(x = 385,
         y = 0.00047,
         labels = 'Variance')
    text(x = 385,
         y = 0.00035,
         labels = 'Ratio')
    legend(
      x = 260,
      y = 0.0003,
      lty = 1:2 ,
      lwd = 2,
      legend = round(h2, 2),
      bty = 'n'
    )
    legend(
      x = 244,
      y = 0.0003,
      pch = my.pch, 
      legend = c('',''),
      bty = 'n'
    )
    text(x = 285,
         y = 0.00032,
         labels = 'Heritability')
  }
  for (kk in seq_along(h2)) {
    for (jj in seq_along(var.ratio)) {
      plot.als <- output[[ll]][[kk]][[jj]][, 'al']
      lines(
        plot.als,
        my.prevs[[ll]][[kk]][, jj],
        col = my.cols[jj],
        lty = kk,
        lwd = 2
      )
      lines(
        plot.als,
        my.multi.norm.prevs[[ll]][[kk]][, jj],
        col = my.cols[jj],
        lty = 3,
        lwd = 2/3
      )
      points(x = 1,
             y = norm.prev[[kk]],
             pch = my.pch[[kk]])
    }
  }
}
dev.off()




pdf(
  'figures/suppFigures/twoEffectStdas.pdf',
  width = 12,
  height = 12
)
par(mfrow = c(2, 2))
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
  height = 12
)
par(mfrow = c(2, 2))
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
  height = 12
)
par(mfrow = c(2, 2))
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
  height = 12
)
par(mfrow = c(2, 2))
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
  'figures/suppFigures/twoEffectDeltals.pdf',
  width = 12,
  height = 12
)
par(mfrow = c(2, 2))
for (ll in seq_along(my.bt)) {
  
  plot(
    NA,
    xlim = range(c(0, my.als)),
    ylim = c(0, max(unlist(my.deltals[[ll]]))),
    ylab = "Risk Effect (large effect loci)",
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
  height = 12
)
par(mfrow = c(2, 2))
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
  height = 12
)
par(mfrow = c(2, 2))
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
  'figures/suppFigures/twoEffectpgl.pdf',
  width = 12,
  height = 12
)
par(mfrow = c(2, 2))
for (ll in seq_along(my.bt)) {
  plot(
    NA,
    xlim = range(c(0, my.als)),
    ylim = c(0, max(unlist(my.pgl[[ll]]))),
    ylab = "Proportion of genetic variance in liability explained by large effect sites",
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
  height = 12
)
par(mfrow = c(2, 2))
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


xmax <- 0.04
plot(
  NA,
  bty = 'n',
  ylim = c(0,0.02),
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
  y = 1-pnorm(dnorminv(ft.range))
)




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




