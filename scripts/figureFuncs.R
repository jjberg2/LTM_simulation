makeEmptyFrame <-
  function(h2 = 0.5,
           my.mean = 0,
           xlim = c(-2, 3),
           n.points = 1000,
           y.max = NULL,
           y.max.factor = 1.1,
           return.y.max = FALSE,
           make.plot = TRUE,
           plot.x.axis = TRUE) {
    my.x <- seq(from = xlim[1],
                to = xlim[2],
                length.out = n.points)
    plot.x <- c(xlim[1], my.x, xlim[2])
    plot.y <- c(0, dnorm(my.x, my.mean, sd = sqrt(h2)), 0)
    if (is.null(y.max))
      y.max <- max(plot.y) * y.max.factor
    if (make.plot) {
      plot(
        NA,
        xlim = xlim,
        ylim = c(0, y.max),
        ylab = '',
        xlab = '',
        bty = 'n',
        xaxt = 'n',
        yaxt = 'n'
      )
      axis(2, labels = FALSE)
      if (plot.x.axis)
        axis(1)
    }
    if (return.y.max)
      y.max
    ## polygon(
    ##     x=plot.x,
    ##     y=plot.y,
    ##     col=adjustcolor('forestgreen',alpha.f=0.05),
    ##     border=NA
    ## )
    ## polygon(
    ##     x=plot.x,
    ##     y=plot.y,
    ##     col=adjustcolor('forestgreen',alpha.f=0.4),
    ##     border=NA,
    ##     density=120,
    ##     angle=315
    ## )
  }

makeGenLi <-
  function(h2 = 0.6,
           my.mean = 0,
           xlim = c(-2, 3),
           n.points = 1000,
           y.max = NULL,
           y.max.factor = 1.1,
           plot.x.axis = TRUE,
           my.color = 'forestgreen',
           axis.cex = 1.6,
           w.thr = FALSE,
           prev = 0.01,
           shade = TRUE,
           plot.y.axis = TRUE) {
    my.x <- seq(from = xlim[1],
                to = xlim[2],
                length.out = n.points)
    plot.x <- c(xlim[1], my.x, xlim[2])
    plot.y <- c(0, dnorm(my.x, my.mean, sd = sqrt(h2)), 0)
    if (is.null(y.max))
      y.max <- max(plot.y) * y.max.factor
    plot(
      NA,
      xlim = xlim,
      ylim = c(0, y.max),
      ylab = '',
      xlab = '',
      bty = 'n',
      xaxt = 'n',
      yaxt = 'n'
    )
    if (plot.y.axis)
      axis(2, labels = FALSE)
    if (plot.x.axis)
      axis(1, cex.axis = axis.cex, labels = FALSE)
    polygon(
      x = plot.x,
      y = plot.y,
      col = adjustcolor(my.color, alpha.f = 0.05),
      border = NA
    )
    if (shade) {
      polygon(
        x = plot.x,
        y = plot.y,
        col = adjustcolor(my.color, alpha.f = 0.4),
        border = NA,
        density = 120,
        angle = 315
      )
    }
    if (w.thr) {
      thr <- qnorm(1 - prev)
      abline(v = thr, lty = 1, lwd = 2)
    }
    return(y.max)
  }
addEnvRiskDemo <-
  function(h2 = 0.5,
           my.g,
           prev = 0.01,
           xlim = c(-2, 3),
           y.max = 1) {
    #recover()
    env.x <- seq(xlim[1], xlim[2], length.out = 1000)
    env.y <-
      dnorm(env.x, mean = my.g, sd = sqrt(1 - h2)) * y.max / dnorm(my.g, my.g, sd =
                                                                     sqrt(1 - h2))
    lines(env.x, env.y, col = 'grey', lwd = 2)
    
    thr <- qnorm(1 - prev)
    dis.x <- env.x[env.x > thr]
    dis.y <- env.y[env.x > thr]
    plot.dis.x <- c(dis.x, xlim[2], thr)
    plot.dis.y <- c(dis.y, 0, 0)
    
    polygon(
      x = plot.dis.x,
      y = plot.dis.y,
      col = 'darkgrey',
      border = NA
    )
    
  }
addRiskCurve <-
  function(h2 = 0.5,
           my.mean = 0,
           phen.var = 1,
           t.pos = NULL,
           prev,
           xlim = c(-3, 4),
           y.max,
           fit.cost = 0.5,
           my.col = 'darkorchid',
           fit.surface = FALSE) {
    if (is.null(t.pos)) {
      t.pos <- qnorm(1 - prev, mean = my.mean, sd = sqrt(phen.var))
    }
    
    my.x <- seq(xlim[1], xlim[2], length.out = 1000)
    my.risk <- pnorm(my.x, t.pos, sd = sqrt(1 - h2)) * y.max
    if (fit.surface) {
      lines(my.x,
            y.max - fit.cost * my.risk,
            lwd = 2,
            col = my.col)
    } else {
      lines(my.x, my.risk, lwd = 2, col = my.col)
    }
    
  }
makeEnvRiskShading <-
  function(env.var = 0.5,
           my.mean = 0,
           xlim = c(-3, 4),
           n.points = 1000,
           y.max = NULL,
           y.max.factor = 1.1,
           plot.x.axis = TRUE,
           prev = 0.01,
           my.color = 'deeppink3') {
    my.x <- seq(from = xlim[1],
                to = xlim[2],
                length.out = n.points)
    
    if (is.null(t.pos))
      t.pos <- qnorm(1 - prev, mean = my.mean, sd = sqrt(phen.var))
    
    plot.x <- c(xlim[1], my.x, xlim[2])
    plot.y <- c(0, dnorm(my.x, my.mean, sd = sqrt(env.var)), 0)
    
    ##healthy.x <- seq(from=xlim[1],to=t.pos,length.out=n.points)
    x.upper.bound <-
      max(t.pos, min(xlim[2], qnorm(1 - 1e-4, my.mean, sqrt(env.var))))
    disease.x <-
      seq(from = t.pos,
          to = x.upper.bound,
          length.out = n.points)
    
    ##plot.healthy.x <- c(xlim[1],healthy.x,t.pos)
    plot.disease.x <- c(t.pos, disease.x, x.upper.bound)
    
    #plot.healthy.y <- c(0,dnorm(healthy.x,my.mean,sd=sqrt(env.var)),0)
    plot.disease.y <-
      c(0, dnorm(disease.x, my.mean, sd = sqrt(env.var)), 0)
    
    
    if (is.null(y.max))
      y.max <- max(plot.y) * y.max.factor
    plot(
      NA,
      xlim = xlim,
      ylim = c(0, y.max),
      ylab = '',
      xlab = '',
      bty = 'n',
      xaxt = 'n',
      yaxt = 'n'
    )
    axis(2, labels = FALSE)
    if (plot.x.axis)
      axis(1, labels = FALSE)
    polygon(
      x = plot.x,
      y = plot.y,
      col = adjustcolor(my.color, alpha.f = 0.05),
      border = NA
    )
    polygon(
      x = plot.disease.x,
      y = plot.disease.y,
      col = adjustcolor(my.color, alpha.f = 0.4),
      border = NA,
      density = 120,
      angle = 315
    )
    return(y.max)
  }
makePhenLiNoThr <-
  function(phen.var = 1,
           my.mean = 0,
           xlim = c(-3, 4),
           n.points = 1000,
           y.max = NULL,
           y.max.factor = 1.1,
           prev = 0.01,
           plot.x.axis = TRUE,
           axis.cex = 1,
           t.pos = NULL,
           shade = c('both', 'healthy', 'disease')) {
    if (is.null(t.pos))
      t.pos <- qnorm(1 - prev, mean = my.mean, sd = sqrt(phen.var))
    
    tmp.x <- seq(from = xlim[1],
                 to = xlim[2],
                 length.out = n.points)
    plot.healthy.y <- c(0, dnorm(tmp.x, my.mean, sd = sqrt(phen.var)), 0)
    plot.healthy.x <- c(xlim[1], tmp.x, xlim[2])
    
    if (is.null(y.max))
      y.max <- max(plot.healthy.y) * y.max.factor
    plot(
      NA,
      xlim = xlim,
      ylim = c(0, y.max),
      ylab = '',
      xlab = '',
      bty = 'n',
      xaxt = 'n',
      yaxt = 'n'
    )
    axis(2, cex.axis = axis.cex, labels = FALSE)
    if (plot.x.axis)
      axis(1, cex.axis = axis.cex, labels = FALSE)
    
    polygon(
      x = plot.healthy.x,
      y = plot.healthy.y,
      col = adjustcolor('blue', alpha.f = 0.05),
      border = NA
    )
    polygon(
      x = plot.healthy.x,
      y = plot.healthy.y,
      col = adjustcolor('blue', alpha.f = 0.4),
      border = NA,
      density = 120,
      angle = 315
    )
  }
makePhenLi <-
  function(phen.var = 1,
           my.mean = 0,
           xlim = c(-3, 4),
           n.points = 1000,
           y.max = NULL,
           y.max.factor = 1.1,
           prev = 0.01,
           plot.x.axis = TRUE,
           axis.cex = 1,
           t.pos = NULL,
           shade = c('both', 'healthy', 'disease')) {
    if (is.null(t.pos))
      t.pos <- qnorm(1 - prev, mean = my.mean, sd = sqrt(phen.var))
    
    healthy.x <- seq(from = xlim[1],
                     to = t.pos,
                     length.out = n.points)
    disease.x <- seq(from = t.pos,
                     to = xlim[2],
                     length.out = n.points)
    
    plot.healthy.x <- c(xlim[1], healthy.x, t.pos)
    plot.disease.x <- c(t.pos, disease.x, xlim[2])
    
    plot.healthy.y <-
      c(0, dnorm(healthy.x, my.mean, sd = sqrt(phen.var)), 0)
    plot.disease.y <-
      c(0, dnorm(disease.x, my.mean, sd = sqrt(phen.var)), 0)
    
    if (is.null(y.max))
      y.max <- max(c(plot.healthy.y, plot.disease.y)) * y.max.factor
    plot(
      NA,
      xlim = xlim,
      ylim = c(0, y.max),
      ylab = '',
      xlab = '',
      bty = 'n',
      xaxt = 'n',
      yaxt = 'n'
    )
    axis(2, cex.axis = axis.cex, labels = FALSE)
    if (plot.x.axis)
      axis(1, cex.axis = axis.cex, labels = FALSE)
    
    polygon(
      x = plot.healthy.x,
      y = plot.healthy.y,
      col = adjustcolor('blue', alpha.f = 0.05),
      border = NA
    )
    if ('healthy' %in% shade | 'both' %in% shade) {
      polygon(
        x = plot.healthy.x,
        y = plot.healthy.y,
        col = adjustcolor('blue', alpha.f = 0.4),
        border = NA,
        density = 120,
        angle = 315
      )
    }
    polygon(
      x = plot.disease.x,
      y = plot.disease.y,
      col = adjustcolor('red', alpha.f = 0.05),
      border = NA
    )
    if ('disease' %in% shade | 'both' %in% shade) {
      polygon(
        x = plot.disease.x,
        y = plot.disease.y,
        col = adjustcolor('red', alpha.f = 0.4),
        border = NA,
        density = 120,
        angle = 315
      )
    }
    abline(v = t.pos, lwd = 2)
    
  }

makePhenWEffect <-
  function(phen.var = 1,
           prev = 0.01,
           risk.effect = 1 / (40000 * 0.5),
           alpha = NULL,
           t.pos = NULL,
           my.mean = 0,
           xlim = c(-3, 4),
           n.points = 1000,
           y.max = NULL,
           y.max.factor = 1.1,
           return.stuff = FALSE,
           make.plot = TRUE,
           bty = 'n',
           yaxt = 'n',
           xaxt = NULL,
           draw.thresholds = TRUE,
           lwd = 3) {
    ## recover()
    
    if (is.null(t.pos))
      t.pos <- qnorm(1 - prev, mean = my.mean, sd = sqrt(phen.var))
    
    if (is.null(alpha)) {
      effect.pos <-
        qnorm(1 - (prev + risk.effect),
              mean = my.mean,
              sd = sqrt(phen.var))
      alpha <- t.pos - effect.pos
    } else {
      effect.pos <- t.pos - alpha
    }
    
    healthy.x <- seq(from = xlim[1],
                     to = effect.pos,
                     length.out = n.points)
    effect.x <- seq(from = effect.pos,
                    to = t.pos,
                    length.out = n.points)
    disease.x <- seq(from = t.pos,
                     to = xlim[2],
                     length.out = n.points)
    
    plot.healthy.x <- c(xlim[1], healthy.x, effect.pos)
    plot.effect.x <- c(effect.pos, effect.x, t.pos)
    plot.disease.x <- c(t.pos, disease.x, xlim[2])
    
    plot.healthy.y <-
      c(0, dnorm(healthy.x, my.mean, sd = sqrt(phen.var)), 0)
    plot.effect.y <-
      c(0, dnorm(effect.x, my.mean, sd = sqrt(phen.var)), 0)
    plot.disease.y <-
      c(0, dnorm(disease.x, my.mean, sd = sqrt(phen.var)), 0)
    
    
    if (is.null(y.max))
      y.max <- max(plot.healthy.y) * y.max.factor
    
    if (make.plot) {
      plot(
        NA,
        xlim = xlim,
        ylim = c(0, y.max),
        ylab = '',
        xlab = '',
        bty = bty,
        yaxt = yaxt,
        xaxt = xaxt
      )
      polygon(
        x = plot.healthy.x,
        y = plot.healthy.y,
        col = adjustcolor('blue', alpha.f = 0.05),
        border = NA
      )
      polygon(
        x = plot.healthy.x,
        y = plot.healthy.y,
        col = adjustcolor('blue', alpha.f = 0.4),
        border = NA,
        density = 120,
        angle = 315
      )
      
      polygon(
        x = plot.effect.x,
        y = plot.effect.y,
        col = rgb(1, 0, 1, 0.05),
        border = NA
      )
      polygon(
        x = plot.effect.x,
        y = plot.effect.y,
        col = rgb(1, 0, 1, 0.4),
        border = NA,
        density = 120,
        angle = 0
      )
      
      polygon(
        x = plot.disease.x,
        y = plot.disease.y,
        col = rgb(1, 0, 0, 0.05),
        border = NA
      )
      polygon(
        x = plot.disease.x,
        y = plot.disease.y,
        col = rgb(1, 0, 0, 0.4),
        border = NA,
        density = 120,
        angle = 90
      )
      if (draw.thresholds) {
        abline(v = t.pos, lwd = lwd)
        abline(v = effect.pos, lty = 2, lwd = lwd)
      }
    }
    
    if (return.stuff)
      list(
        alpha = alpha,
        y.max = y.max,
        plot.effect.x = plot.effect.x,
        plot.effect.y = plot.effect.y
      )
    
  }

makePhenWDiploEffect <-
  function(phen.var = 1,
           prev = 0.01,
           t.pos = NULL,
           risk.effect = 1 / (40000 * 0.5),
           alpha = NULL,
           my.mean = 0,
           xlim = c(-3, 4),
           n.points = 1000,
           y.max = NULL,
           y.max.factor = 1.1,
           make.plot = TRUE) {
    ## recover()
    
    if (is.null(t.pos))
      t.pos <- qnorm(1 - prev, mean = my.mean, sd = sqrt(phen.var))
    
    het.pos <- t.pos - alpha
    hom.pos <- t.pos - 2 * alpha
    
    
    healthy.x <- seq(from = xlim[1],
                     to = hom.pos,
                     length.out = n.points)
    hom.x <- seq(from = hom.pos,
                 to = het.pos,
                 length.out = n.points)
    het.x <- seq(from = het.pos,
                 to = t.pos,
                 length.out = n.points)
    disease.x <- seq(from = t.pos,
                     to = xlim[2],
                     length.out = n.points)
    
    plot.healthy.x <- c(xlim[1], healthy.x, hom.pos)
    plot.hom.x <- c(hom.pos, hom.x, het.pos)
    plot.het.x <- c(het.pos, het.x, t.pos)
    plot.disease.x <- c(t.pos, disease.x, xlim[2])
    
    plot.healthy.y <-
      c(0, dnorm(healthy.x, my.mean, sd = sqrt(phen.var)), 0)
    plot.hom.y <- c(0, dnorm(hom.x, my.mean, sd = sqrt(phen.var)), 0)
    plot.het.y <- c(0, dnorm(het.x, my.mean, sd = sqrt(phen.var)), 0)
    plot.disease.y <-
      c(0, dnorm(disease.x, my.mean, sd = sqrt(phen.var)), 0)
    
    
    if (is.null(y.max))
      y.max <- max(plot.healthy.y) * y.max.factor
    
    if (make.plot) {
      plot(
        NA,
        xlim = xlim,
        ylim = c(0, y.max),
        ylab = '',
        xlab = '',
        bty = 'n'
      )
      polygon(
        x = plot.healthy.x,
        y = plot.healthy.y,
        col = adjustcolor('blue', alpha.f = 0.05),
        border = NA
      )
      polygon(
        x = plot.healthy.x,
        y = plot.healthy.y,
        col = adjustcolor('blue', alpha.f = 0.4),
        border = NA,
        density = 120,
        angle = 315
      )
      
      polygon(
        x = plot.hom.x,
        y = plot.hom.y,
        col = rgb(0.4, 0, 0.8, 0.05),
        border = NA
      )
      polygon(
        x = plot.hom.x,
        y = plot.hom.y,
        col = rgb(0.4, 0, 0.8, 0.4),
        border = NA,
        density = 120,
        angle = 0
      )
      
      polygon(
        x = plot.het.x,
        y = plot.het.y,
        col = rgb(0.8, 0, 0.4, 0.05),
        border = NA
      )
      polygon(
        x = plot.het.x,
        y = plot.het.y,
        col = rgb(0.8, 0, 0.4, 0.4),
        border = NA,
        density = 120,
        angle = 90
      )
      
      polygon(
        x = plot.disease.x,
        y = plot.disease.y,
        col = rgb(1, 0, 0, 0.05),
        border = NA
      )
      polygon(
        x = plot.disease.x,
        y = plot.disease.y,
        col = rgb(1, 0, 0, 0.4),
        border = NA,
        density = 120,
        angle = 45
      )
      abline(v = t.pos)
      abline(v = het.pos, lty = 2)
      abline(v = hom.pos, lty = 3)
      
    }
    
  }

makePhenWHighlightHet <-
  function(phen.var = 1,
           prev = 0.01,
           t.pos = NULL,
           risk.effect = 1 / (40000 * 0.5),
           alpha = NULL,
           my.mean = 0,
           xlim = c(-3, 4),
           n.points = 1000,
           y.max = NULL,
           y.max.factor = 1.1,
           make.plot = TRUE) {
    ## recover()
    
    if (is.null(t.pos))
      t.pos <- qnorm(1 - prev, mean = my.mean, sd = sqrt(phen.var))
    
    het.pos <- t.pos - alpha
    hom.pos <- t.pos - 2 * alpha
    
    
    healthy.x <- seq(from = xlim[1],
                     to = hom.pos,
                     length.out = n.points)
    hom.x <- seq(from = hom.pos,
                 to = het.pos,
                 length.out = n.points)
    het.x <- seq(from = het.pos,
                 to = t.pos,
                 length.out = n.points)
    disease.x <- seq(from = t.pos,
                     to = xlim[2],
                     length.out = n.points)
    
    plot.healthy.x <- c(xlim[1], healthy.x, hom.pos)
    plot.hom.x <- c(hom.pos, hom.x, het.pos)
    plot.het.x <- c(het.pos, het.x, t.pos)
    plot.disease.x <- c(t.pos, disease.x, xlim[2])
    
    plot.healthy.y <-
      c(0, dnorm(healthy.x, my.mean, sd = sqrt(phen.var)), 0)
    plot.hom.y <- c(0, dnorm(hom.x, my.mean, sd = sqrt(phen.var)), 0)
    plot.het.y <- c(0, dnorm(het.x, my.mean, sd = sqrt(phen.var)), 0)
    plot.disease.y <-
      c(0, dnorm(disease.x, my.mean, sd = sqrt(phen.var)), 0)
    
    
    if (is.null(y.max))
      y.max <- max(plot.healthy.y) * y.max.factor
    
    if (make.plot) {
      plot(
        NA,
        xlim = xlim,
        ylim = c(0, y.max),
        ylab = '',
        xlab = '',
        bty = 'n'
      )
      polygon(
        x = plot.healthy.x,
        y = plot.healthy.y,
        col = adjustcolor('blue', alpha.f = 0.05),
        border = NA
      )
      ## polygon(
      ##     x=plot.healthy.x,
      ##     y=plot.healthy.y,
      ##     col=adjustcolor('blue',alpha.f=0.4),
      ##     border=NA,
      ##     density=120,
      ##     angle=315
      ## )
      
      polygon(
        x = plot.hom.x,
        y = plot.hom.y,
        col = rgb(0.4, 0, 0.8, 0.05),
        border = NA
      )
      ## polygon(
      ##     x=plot.hom.x,
      ##     y=plot.hom.y,
      ##     col=rgb(0.4,0,0.8,0.4),
      ##     border=NA,
      ##     density=120,
      ##     angle=0
      ## )
      
      polygon(
        x = plot.het.x,
        y = plot.het.y,
        col = rgb(0.8, 0, 0.4, 0.05),
        border = NA
      )
      polygon(
        x = plot.het.x,
        y = plot.het.y,
        col = rgb(0.8, 0, 0.4, 0.4),
        border = NA,
        density = 120,
        angle = 90
      )
      
      polygon(
        x = plot.disease.x,
        y = plot.disease.y,
        col = rgb(1, 0, 0, 0.05),
        border = NA
      )
      ## polygon(
      ##     x=plot.disease.x,
      ##     y=plot.disease.y,
      ##     col=rgb(1,0,0,0.4),
      ##     border=NA,
      ##     density=120,
      ##     angle=45
      ## )
      abline(v = t.pos)
      abline(v = het.pos, lty = 2)
      abline(v = hom.pos, lty = 3)
      
    }
    
  }


if (FALSE) {
  pdf('figures/emptyFrame.pdf',
      height = 7,
      width = 7)
  y.max <- makeEmptyFrame(xlim = c(-3, 4), return.y.max = TRUE)
  dev.off()
  
  pdf('figures/simpleGenLi.pdf',
      height = 7,
      width = 7)
  makeGenLi(xlim = c(-3, 4))
  dev.off()
  
  pdf('figures/simplePhenLi.pdf',
      height = 7,
      width = 7)
  makePhenLi(xlim = c(-3, 4), y.max = y.max)
  dev.off()
  
  pdf('figures/simplePhenLiWithEffect.pdf',
      height = 7,
      width = 7)
  alpha <-
    makePhenWEffect(
      xlim = c(-3, 4),
      y.max = y.max,
      risk.effect = 0.01,
      return.alpha = TRUE
    )
  dev.off()
  
  t.pos <- qnorm(0.99)
  alpha <- GetAlpha(0.01, 0.01)
  pdf(
    'figures/simplePhenLiWithEffectEnvChange.pdf',
    height = 7,
    width = 7
  )
  alpha <-
    makePhenWEffect(
      xlim = c(-3, 4),
      t.pos = t.pos,
      alpha = alpha,
      my.mean = 0.2
    )
  dev.off()
  
  t.pos <- qnorm(0.99)
  alpha <- 0.01
  pdf(
    'figures/simplePhenLiWithEffectShiftLeft.pdf',
    height = 7,
    width = 7
  )
  makePhenWEffect(
    xlim = c(2, 3),
    y.max = y.max,
    prev = 0.1,
    alpha = alpha,
    t.pos = t.pos
  )
  dev.off()
  
  pdf(
    'figures/simplePhenLiWithEffectShiftRight.pdf',
    height = 7,
    width = 7
  )
  makePhenWEffect(xlim = c(-3, 4),
                  y.max = y.max,
                  prev = 0.01)
  dev.off()
  
  
  alpha <- GetAlpha(0.001, 0.01)
  pdf('figures/zooomedInBothEffectsx0.pdf',
      height = 7,
      width = 7)
  makePhenWDiploEffect(
    xlim = c(2.1, 2.5),
    prev = 0.01,
    alpha = alpha,
    y.max = 0.05
  )
  dev.off()
  
  
  t.pos <- qnorm(0.99)
  pdf('figures/zooomedInBothEffectsx05.pdf',
      height = 7,
      width = 7)
  makePhenWDiploEffect(
    xlim = c(2.1, 2.5),
    my.mean = -alpha,
    t.pos = t.pos,
    alpha = alpha,
    y.max = 0.05
  )
  dev.off()
  
  
  pdf('figures/zooomedInBothEffectsx1.pdf',
      height = 7,
      width = 7)
  makePhenWDiploEffect(
    xlim = c(2.1, 2.5),
    my.mean = -2 * alpha,
    t.pos = t.pos,
    alpha = alpha,
    y.max = 0.05
  )
  dev.off()
  
  
  alpha <- GetAlpha(0.001, 0.01)
  pdf('figures/zooomedInHighlightHet.pdf',
      height = 7,
      width = 7)
  makePhenWHighlightHet(
    xlim = c(2.1, 2.5),
    prev = 0.01,
    alpha = alpha,
    y.max = 0.05
  )
  dev.off()
  
  
  
  
  ### robustness
  alpha <- GetAlpha(0.001, 0.01)
  t.pos <- qnorm(0.99)
  pdf(
    'figures/simplePhenLiWithEffectZoomed.pdf',
    height = 7,
    width = 7
  )
  makePhenWEffect(
    xlim = c(2.1, 2.5),
    t.pos = t.pos,
    y.max = 0.05,
    alpha = alpha
  )
  dev.off()
  
  
  alpha <- GetAlpha(0.001, 0.01)
  pdf(
    'figures/simplePhenLiWithEffectZoomedLowerRisk.pdf',
    height = 7,
    width = 7
  )
  makePhenWEffect(
    xlim = c(2.1, 2.5),
    t.pos = t.pos,
    my.mean = t.pos - qnorm(0.995),
    y.max = 0.05,
    alpha = alpha
  )
  dev.off()
  
  
  
  
  
  source('scripts/freqSpecFuncs.R')
  
  derivedBias <- function(x, gamma) {
    1 / (1 + exp(2 * gamma * (x - 1)))
  }
  
  
  N <- 2e4
  my.gamma <- exp(seq(from = log(1e-3), log(1e3), length.out = 1e4))
  a.few.x <- c(1 / N, 0.1, 0.25, 0.5, 0.75, 0.9, 1 - 1 / N)
  my.biases <- lapply(a.few.x, function(X)
    derivedBias(X, my.gamma))
  my.cols <-
    c('black',
      'red',
      'green',
      'blue',
      'orange',
      'purple',
      'brown',
      'gold')
  
  pdf('figures/derivedBiases.pdf',
      width = 9.22,
      height = 5.86)
  plot(
    NA,
    xlim = c(1e-3, 1e3),
    ylim = c(0.5, 1),
    bty = 'n',
    log = 'x',
    xaxt = 'n',
    xlab = '',
    ylab = ''
  )
  polygon(
    x = c(1e-3, 1e-3, 1e-1, 1e-1),
    y = c(0.5, 1, 1, 0.5),
    col = 'grey',
    border = NA
  )
  polygon(
    x = c(1e1, 1e1, 1e3, 1e3),
    y = c(0.5, 1, 1, 0.5),
    col = 'grey',
    border = NA
  )
  at.these <- 10 ^ seq(-3, 3)
  my.digits <- c(3, 2, 1, 0, 0, 0, 0)
  axis(1, at = at.these, labels = FALSE)
  for (j in 1:length(at.these)) {
    axis(1,
         at = at.these[j],
         labels = formatC(at.these[j], format = 'f', digits = my.digits[j]))
  }
  for (i in 1:length(my.biases)) {
    lines(
      x = my.gamma,
      y = my.biases[[i]],
      col = my.cols[i],
      lwd = 2
    )
  }
  legend(
    x = 1e-3,
    y = 1,
    legend = c('De novo',
               0.1,
               0.25,
               0.5,
               0.75,
               0.9,
               'About to Fix'),
    col = my.cols,
    lwd = 2,
    bty = 'n'
  )
  dev.off()
  
  
  
  N <- 2e4
  my.x <- seq(1 / N, 1 - 1 / N, 1 / N)
  common.x <- seq(0.01, 0.99, 1 / N)
  a.few.gamma <- c(0.01, 0.1, 0.5, 1, 2, 3, 10, 100)
  
  
  
  bias.by.gamma <-
    lapply(a.few.gamma, function(GAMMA)
      derivedBias(my.x, GAMMA))
  
  
  pdf('figures/biasByFreq.pdf',
      width = 10,
      height = 6)
  plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, 1),
    bty = 'n',
    xlab = '',
    ylab = ''
  )
  for (i in 1:length(my.biases)) {
    lines(
      x = my.x,
      y = bias.by.gamma[[i]],
      col = my.cols[i],
      lwd = 2
    )
  }
  legend(
    'bottomright',
    legend = a.few.gamma,
    lwd = 2,
    col = my.cols,
    bty = 'n'
  )
  dev.off()
  
  
  
  
  
  full.dist <- list()
  for (i in 1:length(a.few.gamma)) {
    this.gamma <- a.few.gamma[i]
    approx <- ifelse(this.gamma > 500, TRUE, FALSE)
    denom <-
      integrate(
        function(X)
          VarSpec(X, this.gamma, theta = 1, approx = approx),
        lower = 1 / N,
        upper = 1 - 1 / N
      )$value
    nums <- sapply(my.x,
                   function(Y)
                     integrate(
                       function(X)
                         VarSpec(X, this.gamma, theta = 1, approx = approx),
                       lower = 1 / N,
                       upper = Y
                     )$value)
    full.dist[[i]] <- nums / denom
  }
  
  pdf('figures/cumulativeVariance.pdf',
      width = 10,
      height = 6)
  plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, 1),
    bty = 'n',
    xlab = '',
    ylab = ''
  )
  for (i in 1:length(full.dist)) {
    lines(x = my.x,
          full.dist[[i]],
          col = my.cols[i],
          lwd = 2)
  }
  legend(
    'bottomright',
    legend = a.few.gamma,
    lwd = 2,
    col = my.cols,
    bty = 'n'
  )
  dev.off()
  
  
  
  
  fitCost <- 0.5
  var.per.site <- list()
  for (i in 1:length(a.few.gamma)) {
    this.gamma <- a.few.gamma[i]
    approx <- ifelse(this.gamma > 500, TRUE, FALSE)
    var.per.site[[i]] <- sapply(my.x,
                                function(X)
                                  VarSpec(X, this.gamma, theta = 1, approx = approx) * (this.gamma /
                                                                                          (4 * N * fitCost)) ^ 2)
  }
  
  plot(
    NA,
    xlim = c(0, 1),
    ylim = c(1e-15, max(unlist(var.per.site)) * 1.1),
    bty = 'n',
    xlab = '',
    ylab = '',
    log = 'y'
  )
  for (i in 1:length(full.dist)) {
    lines(x = my.x,
          var.per.site[[i]],
          col = my.cols[i],
          lwd = 2)
  }
  
  
  legend(
    'bottomright',
    legend = a.few.gamma,
    lwd = 2,
    col = my.cols,
    bty = 'n'
  )
  
  
  
  
  
  
  PowerDensity <- function(x, Nsamp, Qsamp, R, pii, gamma, theta, N) {
    power.thr <- qchisq(5e-8, 1, lower.tail = FALSE)
    non.c <- ChiSqNonCentral(
      N = Nsamp,
      Q = Qsamp,
      R = R,
      x = x,
      pii = pii
    )
    power.x <- pchisq(power.thr, 1, ncp = non.c, lower.tail = FALSE)
    x.dens <- SiteSpec(x, gamma, theta, epsilon = 1 / N)
    x.dens * power.x
  }
  
  GetPower <-
    function(Nsamp,
             Qsamp,
             R,
             pii,
             gamma,
             theta,
             N,
             lower = NULL,
             upper = NULL,
             return.upper = FALSE) {
      test.upper <- PowerDensity(upper,
                                 N.samp,
                                 Q.samp,
                                 prev,
                                 this.gamma / (4 * N * fitCost),
                                 this.gamma,
                                 4e-4,
                                 N)
      if (is.null(lower))
        lower <- 1 / N
      if (is.null(upper))
        upper <- 1 - 1 / N
      if (is.na(test.upper)) {
        while (is.na(test.upper)) {
          test.upper <- PowerDensity(upper,
                                     N.samp,
                                     Q.samp,
                                     prev,
                                     this.gamma / (4 * N * fitCost),
                                     this.gamma,
                                     4e-4,
                                     N)
          upper <- upper - 1 / N
        }
      }
      my.power <- integrate(function(X) {
        PowerDensity(X,
                     N.samp,
                     Q.samp,
                     prev,
                     this.gamma / (4 * N * fitCost),
                     this.gamma,
                     4e-4,
                     N)
      },
      lower = lower,
      upper = upper)$value
      if (return.upper)
        c(my.power, upper)
      else
        my.power
    }
  
  
  
  
  
  my.gamma <- exp(seq(from = log(1e-3), log(1e3), length.out = 1e3))
  fitCost <- 0.05
  seq.power <- numeric()
  geno.power <- numeric()
  large.s.var <- numeric()
  N.samp <- 200000
  Q.samp <- 0.4
  prev <- 0.01
  n.gamma <- length(my.gamma)
  pb <- txtProgressBar(min = 1, max = n.gamma)
  seq.lower <- 1 / N
  current.seq.upper <- 1 - 1 / N
  geno.lower <- 0.01
  current.geno.upper <- 0.99
  for (i in 1:n.gamma) {
    this.gamma <- my.gamma[i]
    pii <- this.gamma / (2 * N * fitCost)
    geno.out <-
      GetPower(
        N.samp,
        Q.samp,
        prev,
        pii,
        this.gamma,
        4e-4,
        N,
        lower = 0.05,
        upper = current.geno.upper,
        return.upper = TRUE
      )
    geno.power[i] <- geno.out[1]
    current.geno.upper <- geno.out[2]
    seq.out <-
      GetPower(
        N.samp,
        Q.samp,
        prev,
        pii,
        this.gamma,
        4e-4,
        N,
        lower = 1 / N,
        upper = current.seq.upper,
        return.upper = TRUE
      )
    seq.power[i] <- seq.out[1]
    current.seq.upper <- seq.out[2]
    setTxtProgressBar(pb, i)
  }
  
  
  
  
  pdf('figures/powerAsFunctionOfGamma.pdf',
      width = 10,
      height = 6)
  plot(
    NA,
    xlim = c(10 ^ -3, 10 ^ 3),
    ylim = c(0, 1),
    xlab = '',
    ylab = '',
    bty = 'n',
    log = 'x',
    xaxt = 'n'
  )
  axis(
    1,
    at = 10 ^ seq(-3, 3),
    labels = c('0.001', '0.01', '0.1', '1', '10', '100', '1000')
  )
  lines(my.gamma,
        geno.power,
        lwd = 2,
        col = 'red')
  lines(my.gamma,
        seq.power,
        lwd = 2)
  legend(
    'topleft',
    legend = c('sequencing', 'genotyping'),
    lwd = c(2, 2),
    col = c('black', 'red'),
    bty = 'n'
  )
  dev.off()
  
  
  
  
  power.density <- list()
  for (i in 1:length(a.few.gamma)) {
    this.gamma <- a.few.gamma[i]
    pii <- this.gamma / (2 * N * fitCost)
    power.density[[i]] <-
      PowerDensity(my.x, N.samp, Q.samp, prev, pii, this.gamma, 4e-4, N)
  }
  
  
  pdf('figures/Gamma1HitDensity.pdf',
      width = 10,
      height = 6)
  plot(
    my.x,
    power.density[[4]],
    lwd = 2,
    type = 'l',
    xlab = '',
    ylab = '',
    bty = 'n'
  )
  dev.off()
  
  
  pdf('figures/Gamma2HitDensity.pdf',
      width = 10,
      height = 6)
  plot(
    my.x,
    power.density[[5]],
    lwd = 2,
    type = 'l',
    xlab = '',
    ylab = '',
    bty = 'n'
  )
  dev.off()
  
  
  pdf('figures/Gamma3HitDensity.pdf',
      width = 10,
      height = 6)
  plot(
    my.x,
    power.density[[6]],
    lwd = 2,
    type = 'l',
    xlab = '',
    ylab = '',
    bty = 'n'
  )
  dev.off()
  
  
  pdf('figures/Gamma10HitDensity.pdf',
      width = 10,
      height = 6)
  plot(
    my.x,
    power.density[[7]],
    lwd = 2,
    type = 'l',
    xlab = '',
    ylab = '',
    bty = 'n'
  )
  lines(
    my.x,
    ifelse(my.x > 0.01 & my.x < 0.99, power.density[[7]], 0),
    lty = 2,
    lwd = 3,
    col = 'red'
  )
  dev.off()
  
  pdf('figures/Gamma100HitDensity.pdf',
      width = 10,
      height = 6)
  plot(
    my.x,
    power.density[[8]],
    lwd = 2,
    type = 'l',
    xlab = '',
    ylab = '',
    bty = 'n'
  )
  lines(
    my.x,
    ifelse(my.x > 0.01 & my.x < 0.99, power.density[[7]], 0),
    lty = 2,
    lwd = 3,
    col = 'red'
  )
  dev.off()
  
  
  
  
  
  library(RColorBrewer)
  my.cols <- c(rev(brewer.pal(5, 'YlOrRd')), brewer.pal(5, 'YlGnBu'))
  pm.gamma <- c(-100, -10, -5, -1, -0.1, 0.1, 1, 5, 10, 100)
  specs <- lapply(pm.gamma, function(GAMMA)
    FreqSpec(my.x, GAMMA, 1))
  
  pdf('figures/freqSpecs.pdf',
      width = 10,
      height = 6)
  plot(
    NA,
    ylim = c(1e-8, max(unlist(specs))),
    xlim = c(0, 1),
    log = 'y',
    bty = 'n',
    xlab = '',
    ylab = ''
  )
  for (i in 1:length(specs)) {
    lines(my.x,
          specs[[i]],
          lwd = 2,
          col = my.cols[i])
  }
  dev.off()
  
  
  
  my.alphas <- seq(1e-8, qnorm(1 - 1e-4), length.out = 1e4)
  my.prevs <- 10 ^ -(4:2)
  my.ts <- qnorm(1 - my.prevs)
  my.piis <- lapply(my.ts,
                    function(MYT)
                      pnorm(MYT - my.alphas, lower.tail = FALSE) - pnorm(MYT, lower.tail =
                                                                           FALSE))
  
  
  
  three.cols <- brewer.pal(3, 'YlGnBu')
  plot(
    NA,
    xlim = c(0, max(my.alphas)),
    ylim = c(0, max(unlist(my.piis))),
    lwd = 2,
    bty = 'n'
  )
  for (i in 1:length(my.piis)) {
    lines(my.alphas,
          my.piis[[i]],
          lwd = 2,
          col = three.cols[i])
  }
  
  
  ## skewed distributions
  t <-
    get(load('talks/2019.03.18_labmeeting/threeRisks.tapprox.Robj'))
  
  
  
  
  
  
  
  
  pdf('figures/shiftPiiChanges.pdf',
      height = 10,
      width = 6)
  plot(
    NA,
    lwd = 2,
    type = 'l',
    xlim = c(0, 0.9),
    ylim = c(0, 0.9),
    bty = 'n',
    xlab = '',
    ylab = ''
  )
  lines(my.piis[[1]],
        my.piis[[3]],
        lwd = 2,
        col = three.cols[1])
  lines(my.piis[[2]],
        my.piis[[3]],
        lwd = 2,
        col = three.cols[2])
  abline(a = 0,
         b = 1,
         lty = 2)
  dev.off()
  
  
}
