rm(list=ls())
library('MetBrewer')
source('scripts/solveTwoEffect.R')
re.sim <- F
if(re.sim) source('scripts/twoEffectPrevalenceSolutions.R')
solns <- get(load(file='solutions/twoEffectPrevalenceSolutions_solutions.Robj'))
output <- get(load(file='solutions/twoEffectPrevalenceSolutions_output.Robj'))


my.bt <- sapply(output,function(X) unique(round(X[[1]][[1]]$bt,2)))
h2 <- sapply(output[[1]],function(X) unique(round(X[[1]]$h2,3)))
var.ratio <- sapply(output[[1]][[1]],function(X) unique(X$var.ratio))
Ne <- unique(output[[1]][[1]][[1]]$Ne)
cost <- unique(output[[1]][[1]][[1]]$C)
L <- unique(round(output[[1]][[1]][[1]]$Ls + output[[1]][[1]][[1]]$Ll,0))
u <- unique(output[[1]][[1]][[1]]$u)
my.als <- unique(output[[1]][[1]][[1]]$al)

my.deltals <- list()
my.prevs <- list()
my.multi.norm.prevs <- list()
my.std.as <- list()
my.std.fts <- list()
for (l in seq_along(my.bt)) {
  my.multi.norm.prevs[[l]] <- list()
  my.prevs[[l]] <- list()
  my.deltals[[l]] <- list()
  my.std.as[[l]] <- list()
  my.std.fts[[l]] <- list()
  for (k in seq_along(h2)) {
    my.deltals[[l]][[k]] <- sapply(output[[l]][[k]], function(X)
      X$deltal)
    my.prevs[[l]][[k]] <- sapply(output[[l]][[k]], function(X)
      X$prev)
    my.multi.norm.prevs[[l]][[k]] <- sapply(output[[l]][[k]], function(X)
      X$multi.norm.prev)
    my.std.as[[l]][[k]] <- sapply(output[[l]][[k]], function(X)
      X$std.as)
    my.std.fts[[l]][[k]] <- sapply(output[[l]][[k]], function(X)
      X$std.ft)
  }
}

my.pch = c(21, 22)
my.cols = met.brewer('Isfahan2', length(var.ratio))
pdf(
  'figures/suppFigures/twoEffectPrevalence.pdf',
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
  norm.dens <- 2 * L * u * my.bt[ll] * norm.astd / (h2 * cost)
  norm.prev <- 1 - pnorm(dnorminv(norm.dens))
  plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, max(unlist(my.prevs[[ll]]))),
    ylab = "Disease prevalence",
    xlab = "PAR of large effect alleles",
    main = paste('b_T = ', my.bt[ll], sep = '')
  )
  
  if (ll == 1) {
    legend(
      'bottomright',
      col = my.cols,
      lty = 1 ,
      lwd = 2,
      legend = round(var.ratio,2),
      bty = 'n'
    )
    text(x = 0.94,
         y = 0.00062,
         labels = 'Variance')
    text(x = 0.94,
         y = 0.0005,
         labels = 'Ratio')
    legend(
      x = 0.55,
      y = 0.0005,
      lty = 1:2 ,
      lwd = 2,
      legend = round(h2, 2),
      bty = 'n'
    )
    legend(
      x = 0.5,
      y = 0.0005,
      pch = my.pch, 
      legend = c('',''),
      bty = 'n'
    )
    text(x = 0.62,
         y = 0.0005,
         labels = 'Heritability')
  }
  for (kk in seq_along(h2)) {
    for (jj in seq_along(var.ratio)) {
      ##plot.als <- output[[ll]][[kk]][[jj]][, 'al']
      plot.these <- 1:min(which(my.deltals[[ll]][[kk]][, jj] > (1-max(my.prevs[[ll]][[kk]][, jj])-0.01)))
      lines(
        my.deltals[[ll]][[kk]][plot.these, jj],
        my.prevs[[ll]][[kk]][plot.these, jj],
        col = my.cols[jj],
        lty = kk,
        lwd = 2
      )
      lines(
        my.deltals[[ll]][[kk]][plot.these, jj],
        my.multi.norm.prevs[[ll]][[kk]][plot.these, jj],
        col = my.cols[jj],
        lty = 3,
        lwd = 2/3
      )
      points(x = 0,
             y = norm.prev[[kk]],
             pch = my.pch[[kk]])
    }
  }
}
dev.off()





my.pch = c(21, 22)
my.cols = met.brewer('Isfahan2', length(var.ratio))
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
  ## norm.dens <- 2 * L * u * my.bt[ll] * norm.astd / (h2 * cost)
  ## norm.prev <- 1 - pnorm(dnorminv(norm.dens))
  plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, max(unlist(my.std.as[[ll]]))),
    ylab = "Standadrized liability effect of small effect sites",
    xlab = "PAR of large effect sites",
    main = paste('b_T = ', my.bt[ll], sep = '')
  )
  
  if (ll == 1) {
    legend(
      'bottomright',
      col = my.cols,
      lty = 1 ,
      lwd = 2,
      legend = round(var.ratio,2),
      bty = 'n'
    )
    text(x = 0.94,
         y = 0.0055,
         labels = 'Variance')
    text(x = 0.94,
         y = 0.0043,
         labels = 'Ratio')
    legend(
      x = 0.55,
      y = 0.004,
      lty = 1:2 ,
      lwd = 2,
      legend = round(h2, 2),
      bty = 'n'
    )
    legend(
      x = 0.5,
      y = 0.004,
      pch = my.pch, 
      legend = c('',''),
      bty = 'n'
    )
    text(x = 0.62,
         y = 0.004,
         labels = 'Heritability')
  }
  for (kk in seq_along(h2)) {
    for (jj in seq_along(var.ratio)) {
      ##plot.als <- output[[ll]][[kk]][[jj]][, 'al']
      plot.these <- 1:min(which(my.deltals[[ll]][[kk]][, jj] > (1-max(my.prevs[[ll]][[kk]][, jj])-0.01)))
      lines(
        my.deltals[[ll]][[kk]][plot.these, jj],
        my.std.as[[ll]][[kk]][plot.these, jj],
        col = my.cols[jj],
        lty = kk,
        lwd = 2
      )
      ## these.deltals <- my.deltals[[ll]][[kk]][, jj]
      abline(
        h=norm.astd[jj],
        col = 'black',
        lty = 3,
        lwd = 2
      )
    }
  }
}
dev.off()



my.pch = c(21, 22)
my.cols = met.brewer('Isfahan2', length(var.ratio))
pdf(
  'figures/suppFigures/twoEffectStdFts.pdf',
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
  norm.dens <- 2 * L * u * my.bt[ll] * norm.astd / (h2 * cost)
  norm.prev <- 1 - pnorm(dnorminv(norm.dens))
  plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, max(unlist(my.std.fts[[ll]]))),
    ylab = "Disease prevalence",
    xlab = "PAR of large effect alleles",
    main = paste('b_T = ', my.bt[ll], sep = '')
  )
  
  if (ll == 1) {
    legend(
      'bottomright',
      col = my.cols,
      lty = 1 ,
      lwd = 2,
      legend = round(var.ratio,2),
      bty = 'n'
    )
    text(x = 0.94,
         y = 0.00062,
         labels = 'Variance')
    text(x = 0.94,
         y = 0.0005,
         labels = 'Ratio')
    legend(
      x = 0.55,
      y = 0.0005,
      lty = 1:2 ,
      lwd = 2,
      legend = round(h2, 2),
      bty = 'n'
    )
    legend(
      x = 0.5,
      y = 0.0005,
      pch = my.pch, 
      legend = c('',''),
      bty = 'n'
    )
    text(x = 0.62,
         y = 0.0005,
         labels = 'Heritability')
  }
  for (kk in seq_along(h2)) {
    for (jj in seq_along(var.ratio)) {
      ##plot.als <- output[[ll]][[kk]][[jj]][, 'al']
      plot.these <- 1:min(which(my.deltals[[ll]][[kk]][, jj] > (1-max(my.prevs[[ll]][[kk]][, jj])-0.01)))
      lines(
        my.deltals[[ll]][[kk]][plot.these, jj],
        my.std.fts[[ll]][[kk]][plot.these, jj],
        col = my.cols[jj],
        lty = kk,
        lwd = 2
      )
      these.deltals <- my.deltals[[ll]][[kk]][plot.these, jj]
      abline(
        h=norm.dens[jj],
        col = 'black',
        lty = 3,
        lwd = 2
      )
      # points(x = 0,
      #        y = norm.prev[[kk]],
      #        pch = my.pch[[kk]])
    }
  }
}
dev.off()


