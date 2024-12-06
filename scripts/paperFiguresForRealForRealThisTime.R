setwd('~/Dropbox/msdbPaper')
##########################
### Figure 1 schematic ###
##########################
{
  rm(list=ls())
  h2 <- 0.5
  source('scripts/figureFuncs.R')
  my.cols <- c('forestgreen', 'firebrick4','dodgerblue4')
  my.cols <- c(wes_palette('Rushmore1')[c(3,4)],wes_palette('FantasticFox1')[5])[c(1,3,2)]
  png(
    'figures/paperFiguresForRealForRealThisTime/SchematicFigure1.png',
    height = 9 ,
    width = 22,
    units = 'cm',
    res = 500
  )
  nf <- layout(matrix(c(1, 2), nrow = 1, byrow = TRUE))
  
  ## 1A
  ## mid zoom liability distribution
  op1 <- par(mar = c(2, 2, 1, 1.3))
  {
    makePhenLi(
      xlim = c(-3, 4),
      prev = 0.01,
      shade = 'disease',
      w.thr = FALSE,
      cols = my.cols
    )
    makeGenLi(
      h2 = h2,
      xlim = c(-3, 4),
      plot.x.axis = FALSE,
      plot.y.axis = FALSE,
      add = TRUE,
      y.modifier = dnorm(0) / dnorm(0, sd = sqrt(h2)),
      my.color = my.cols[1]
    )
  }
  
  ## lines drawn on the plot
  {
    abline(v = 0 , lty = 2 , lwd = 2, col = my.cols[1])
    lines(
      x = rep(qnorm(0.99), 2),
      y = c(0, dnorm(0) * 0.81) ,
      lty = 1 ,
      lwd = 2 ,
      col = my.cols[2]
    )
    arrows(
      x0 = 0,
      x1 = -h2 * 2,
      y0 = dnorm(0) * 0.35 ,
      y1 = dnorm(0) * 0.35 ,
      length = 0.05
    )
    arrows(
      x0 = 0 ,
      x1 = -1.53,
      y0 = dnorm(0) * 0.3 ,
      y1 = dnorm(0) * 0.3 ,
      length = 0.05
    )
  }
  
  ## text drawn on the plot
  {
    text(
      x = -0.45 ,
      y = dnorm(0) * 0.4 ,
      labels = expression(sqrt(V[A])) ,
      cex = 0.86
    )
    text(
      x = -0.6 ,
      y = dnorm(0) * 0.25 ,
      labels = expression(sqrt(V[A] + V[E])) ,
      cex = 0.86
    )
    text(
      x = -1.4 ,
      y = dnorm(0) * 1.1 ,
      labels = expression(paste('mean liability, ', bar(G), sep = '')),
      cex = 1
    )
    text(
      x = 3.3 ,
      y = dnorm(0) * 0.78 ,
      labels = 'threshold, T',
      cex = 1 ,
      col = my.cols[2]
    )
    text(
      x = 3.3 ,
      y = dnorm(0) * 0.27 ,
      labels = 'individuals' ,
      cex = 0.7 ,
      col = my.cols[2]
    )
    text(
      x = 3.4 ,
      y = dnorm(0) * 0.22 ,
      labels = 'that develop' ,
      cex = 0.7 ,
      col = my.cols[2]
    )
    text(
      x = 3.37 ,
      y = dnorm(0) * 0.18 ,
      labels = 'the disease' ,
      cex = 0.7 ,
      col = my.cols[2]
    )
    arrows(
      x0 = 2.9 ,
      y0 = dnorm(0) * 0.15 ,
      x1 = 2.77 ,
      y = dnorm(0)  * 0.04 ,
      col = my.cols[2] ,
      length = 0.1
    )
  }
  
  ## axis labels and  legend
  {
    mtext(
      side = 1 ,
      text = 'Liability',
      line = 1,
      cex = 1.2
    )
    mtext(
      side = 2 ,
      text = 'Density',
      line = 1,
      cex = 1.2
    )
    mtext(
      side = 3 ,
      text = 'A' ,
      line = -0.4,
      at = -3.5,
      cex = 1.5
    )
    legend(
      'topright' ,
      lty = 1 ,
      lwd = 1.5 ,
      col = c(my.cols[3], my.cols[1]) ,
      legend = c('total liability' , 'genetic liability') ,
      bty = 'n'
    )
  }
  
  par(op1)
  
  {## 1B
    
    this.col = my.cols[2]
    y.max.factor = 2
    op4 <- par(mar = c(2, 2, 1, 1))
    makeEmptyFrame(
      xlim = c(-3, 6.2) ,
      y.max.factor = 2,
      plot.x.axis = FALSE
    )
    axis(side = 1 , at = seq(-3,6,by=1) , labels = rep('',6-(-3)+1))
    this.ymax <- dnorm(0, sd = sqrt(h2)) * y.max.factor
    abline(
      h = (this.ymax / 1.2) / (1 - 0.01 * 0.5) ,
      lty = 2 ,
      col = my.cols[2]
    )
    lines(
      x = c(-4,-0.43) ,
      y = rep((0.5 * this.ymax / 1.2) / (1 - 0.01 * 0.5) , 2 ) ,
      lty = 2 ,
      col = my.cols[2]
    )
    lines(
      x = c(0.44,7) ,
      y = rep((0.5 * this.ymax / 1.2) / (1 - 0.01 * 0.5) , 2 ) ,
      lty = 2 ,
      col = my.cols[2]
    )
    text(
      x = 3.58 ,
      y = this.ymax*1.11/1.2 ,
      labels = 'threshold, T',
      cex = 1 ,
      col = my.cols[2]
    )
    text(
      x = -1 ,
      y = this.ymax*1.06/1.2 ,
      labels = expression(paste('1/', (1-bar(R)*C),sep = '') ),
      cex = 1 ,
      col = my.cols[2]
    )
    text(
      x = 5.17 ,
      y = this.ymax*0.56/1.2 ,
      labels = expression(paste('(1-C)/', (1-bar(R)*C),sep = '') ),
      cex = 1 ,
      col = my.cols[2]
    )
    {
      
      this.ymax <-
        makeGenLi(
          xlim = c(-3, 4),
          h2 = h2,
          plot.y.axis = FALSE ,
          plot.x.axis = FALSE ,
          y.max.factor = 2,
          add = TRUE ,
          my.color=my.cols[1]
        )
      lines(
        x = rep(qnorm(1 - 0.01),2),
        y = c(0,this.ymax*1.14/1.2),
        lty = 1 ,
        lwd = 2 ,
        col = my.cols[2],
      )
      addRiskCurve(
        h2 = h2,
        prev = 0.01,
        xlim = c(-3, 6) ,
        y.max = (this.ymax / 1.2) / (1 - 0.01 * 0.5),
        my.col = this.col,
        fit.surface = TRUE
      )
      mtext(side = 1 ,
            text = 'Genetic Liability',
            line = 2.2)
      mtext(
        side = 2 ,
        text = 'Expected fitness',
        line = 1,
        las = 0
      )
      axis(
        side = 2,
        at = seq(0, 2 * this.ymax / 1.2, length.out = 3),
        labels = seq(0, 2, by = 1),
        las = 1
      )
      mtext(
        side = 3 ,
        text = 'B' ,
        line = -0.4,
        at = -3,
        cex = 1.5
      )
    }
    par(op4)
  }
  dev.off()
}

##########################
### Figure 2 schematic ###
##########################
{
  rm(list=ls())
  png(
    'figures/paperFiguresForRealForRealThisTime/SchematicFigure2.png',
    height = 18 ,
    width = 22,
    units = 'cm',
    res = 500
  )
  {
    my.cols <- c(
      wes_palette('Rushmore1')[4] ,
      wes_palette('FantasticFox1')[c(4,5)]
    )
    
    op1 <- par(mar = c(1.3, 1.7, 1, 0) + 0.1)
    source('scripts/figureFuncs.R')
    ## y.max <- 0.6
    t.pos <- qnorm(0.99)
    alphaS <- 0.09
    my.xlim <- c(1.5, 2.6)
    v.offset <- 0.007
    my.out <- makePhenWEffect(
      xlim = my.xlim,
      alphaS = alphaS,
      t.pos = t.pos,
      v.offset = v.offset ,
      xaxt = 'n',
      return.stuff = TRUE ,
      cols = my.cols
    )
    y.max <- my.out$y.max
    effect.x <- my.out[[3]]
    effect.y <- my.out[[4]]
  }
  par(op1)
  dev.off()
}

################################
### Figure 3 equlibrium bias ###
################################
{
  rm(list=ls())
  library(wesanderson)
  png(
    'figures/paperFiguresForRealForRealThisTime/Figure3.png',
    height = 7 ,
    width = 17.75,
    units = 'cm',
    res = 600
  )
  layout(matrix(c(1, 2, 3),nrow=1))
  
  {
    ############
    #### 3A ####
    ############
    {
      op3 <- par(mar = c(3.1, 2.6, 0, 0.1),
                 mgp = c(3, 0.6, 0))
      my.gamma <-
        10 ^ seq(log(0.01, 10), log(100, 10), length.out = 1000)
      bias <-
        ifelse(my.gamma < 700, (exp(my.gamma) - 1) / (exp(my.gamma) + 1), 1)
      plot(
        x = my.gamma ,
        y = bias ,
        log = 'x' ,
        type = 'n' ,
        bty = 'n' ,
        lwd = 2 ,
        xlab = '' ,
        ylab = '',
        xaxt = 'n',
        yaxt = 'n'
      )
      mtext(
        text = 'Fixation bias',
        side = 2,
        line = 1.4,
        cex = 0.7
      )
      mtext(
        text = 'Population scaled selection coefficient',
        side = 1,
        line = 1.5,
        cex = 0.7,
        at = 0.833666
      )
      axis(
        side = 1 ,
        at = c(0.01, 0.1, 1, 10, 100),
        labels = c('0.01', '0.1', '1', '10', '100'),
        cex.axis = 0.7
      )
      axis(
        side = 2 ,
        at = 0:5 / 5,
        labels = c('0', '0.2', '0.4', '0.6', '0.8', '1'),
        cex.axis = 0.7,
        las = 2
      )
      
      
      # regime markings and labels
      {
        polygon(
          x = c(0.01, 0.2, 0.2, 0.01) ,
          y = c(1, 1, 0, 0) ,
          col = adjustcolor('grey89') ,
          border = NA
        )
        polygon(
          x = c(5, 100, 100, 5) ,
          y = c(1, 1, 0, 0) ,
          col = adjustcolor('grey89') ,
          border = NA
        )
        text(
          labels = 'effectively',
          x = 0.04,
          y = 0.2 ,
          cex = 1
        )
        text(
          labels = 'neutral',
          x = 0.04,
          y = 0.15,
          cex = 1
        )
        text(
          labels = 'weakly',
          x = 0.7,
          y = 0.9 ,
          cex = 1
        )
        text(
          labels = 'selected',
          x = 0.7,
          y = 0.85,
          cex = 1
        )
        text(
          labels = 'strongly',
          x = 20,
          y = 0.9 ,
          cex = 1
        )
        text(
          labels = 'selected',
          x = 20,
          y = 0.85,
          cex = 1
        )
        
      }
      
      lines(x = my.gamma,
            y = bias,
            lwd = 2)
      
      mtext(text = 'A',side = 2, line = 0.4 , las = 2 , at = 0.9)
      
      par(op3)
    }
      
    ############
    #### 3B ####
    ############
    {
      freqSpec <- function(x, y) {
        2 * exp(-2 * y * x) / ((1 + exp(-2 * y)) * x * (1 - x))
      }
      my.x <- seq(1/40000,1-1/40000,length.out=1000000)
      my.y <- c(0.1, 1, 10,100)
      tmp.freq.specs <- lapply(my.y, function(Y)
        freqSpec(my.x, Y))
      prob.dens <- TRUE
      prob.freq.specs <- lapply(tmp.freq.specs , function (X ) X / sum(X))
      my.cdfs <- lapply(prob.freq.specs , cumsum )
      plot.these <- lapply ( my.cdfs , function(X) X < 0.999999)
      if(prob.dens){
        ymax <- freqSpec(0.01,1)/sum(tmp.freq.specs[[2]])
        my.freq.specs <- prob.freq.specs
      } else {
        ymax <- freqSpec(0.01,1)
        my.freq.specs <- tmp.freq.specs
      }
      
      
      my.cols <- wes_palette("AsteroidCity3")
      op4 <- par(mar = c(3.1, 2.1, 0, 0.4),
                 mgp = c(3, 0.5, 0))
      plot(
        NA,
        type = 'n',
        xlim = c(0,1),
        ylim = c(0,ymax) ,
        xaxt = 'n' ,
        yaxt = 'n' ,
        bty = 'n'
      )
      for (i in 1:length(my.freq.specs)) {
        lines(my.x[plot.these[[i]]] ,
              my.freq.specs[[i]][plot.these[[i]]],
              col = my.cols[i],
              lwd = 1.7)
      }
      axis(1, cex.axis = 0.75)
      axis(
        2,
        at = c(0,ymax) ,
        labels = FALSE ,
        cex.axis = 0.75,
        las = 2
      )
      mtext(
        text = 'Frequency of liability increasing allele',
        side = 1,
        line = 1.5,
        cex = 0.7
      )
      mtext(
        text = 'Density',
        side = 2,
        line = 0.5,
        cex = 0.7
      )
      legend(
        'top' ,
        lwd = 2.7 ,
        lty = 1,
        col = my.cols ,
        legend = c('1/10', '1', '10', '100'),
        title = expression(paste('Scaled coef. (', gamma, ')' , sep = '')) ,
        bty = 'n'
      )
      mtext(text = 'B',side = 2, line = 0.4 , las = 2, at = 0.9*ymax)
      par(op4)
    }
      
    ############
    #### 3C ####
    ############
    {
      op5 <- par(mar = c(3.1, 2.5, 0, 0.4),
                 mgp = c(3, 0.8, 0))
      my.gamma <- exp(seq(log(0.01), log(100), length.out = 10000))
      bias <-
        ifelse(my.gamma < 700, (exp(my.gamma) - 1) / (exp(my.gamma) +
                                                        1), 1)
      het <- 2 * bias / my.gamma
      col1 <- wes_palette('FantasticFox1')[1]
      col2 <- wes_palette('FantasticFox1')[3]
      
      
      plot(
        my.gamma ,
        het ,
        log = 'xy' ,
        type = 'n' ,
        ylim = c(1 / 20, 1.4) ,
        xlim = c(0.01, 100) ,
        xaxt = 'n' ,
        yaxt = 'n' ,
        bty = 'n' ,
        xlab = '',
        ylab = ''
      )
      polygon(
        x = c(0.01, 0.2, 0.2, 0.01) ,
        y = c(1.4, 1.4, 1/20, 1/20) ,
        col = adjustcolor('grey89') ,
        border = NA
      )
      polygon(
        x = c(5, 100, 100, 5) ,
        y = c(1.4, 1.4, 1/20, 1/20) ,
        col = adjustcolor('grey89') ,
        border = NA
      )
      lines(
        x = c(0.01,100) ,
        y = c(1,1) ,
        lty = 3 ,
        lwd = 2.3 ,
        col = my.cols[1]
      )
      lines(
        x = my.gamma,
        y = 2 / my.gamma,
        lty = 3 ,
        lwd = 2.3,
        col = my.cols[4]
      )
      lines(my.gamma ,
            het ,
            lwd = 3)
      text(
        x = 0.045,
        y = 1.14 ,
        cex = 1.2 ,
        labels = expression("E[h(x)]"%~~% theta ) ,
        col = my.cols[1]
      )
      text(
        x = 16.5,
        y = 0.195 ,
        cex = 1.2 ,
        srt = 360 - 74 ,
        labels = expression("E[h(x)]"%~~% theta/gamma == 2*u/s) ,
        col = my.cols[4]
      )
      axis(
        side = 1 ,
        at = c(0.01, 0.1, 1, 10, 100) ,
        labels = c('0.01', '0.1', '1', '10', '100') ,
        cex.axis = 0.7
      )
      axis(
        side = 2 ,
        at = c(1 / 20, 1.4) ,
        labels = FALSE ,
        las = 2
      )
      mtext(
        text = 'log(heterozygosity)',
        side = 2,
        line = 0.5,
        cex = 0.7
      )
      mtext(
        text = 'Population scaled selection coefficient',
        side = 1,
        line = 1.5,
        cex = 0.7 ,
        at = 0.833666
      )
      mtext(text = 'C',side = 2, line = 0.4 , las = 2, at = 0.9*1.2)
      par(op5)
    }
  }
  dev.off()
}



##################
#### Figure 4 ####
##################
{
  rm(list=ls())
  library(wesanderson)
  sim.results <- get(load('figures/smallEffectInsensitivityResultsTable.Rdata'))
  png(
    'figures/paperFiguresForRealForRealThisTime/Figure4.png',
    height = 7 ,
    width = 17.75,
    units = 'cm',
    res = 600
  )
  layout(matrix(c(1, 2, 3),nrow=1),widths = c(0.8,1,1))
  
  ############
  #### 4A ####
  ############
  {
    sim.results$costFactor <- as.factor(sim.results$cost)
    my.cols <- wes_palette("Darjeeling1")[c(2, 4)]
    my.bt <- seq(1e-6, 1 - 1e-6 , length.out = 1e5)
    my.gamma <- 0.5 * log((1 + my.bt) / (1 - my.bt))
    these.costs <- unique(sim.results$cost)
    op3 <- par(mar = c(2.9, 2.6, 0, 0.1),
               mgp = c(3, 0.6, 0))
    plot(
      NA,
      xlim = c(0, 1) ,
      ylim = c(0, 6) ,
      bty = 'n' ,
      ylab = '' ,
      xlab = ''
    )
    lines(x = my.bt ,
          y = my.gamma)
    mtext(
      side = 1 ,
      text = expression(paste('Relative threshold position, ' , b[T] , sep = '')) ,
      line = 1.5,
      cex = 0.7,
    )
    mtext(
      side = 2 ,
      text = 'Population scaled selection coefficient' ,
      line = 1.5,
      cex = 0.7,
    )
    points(
      x = sim.results$b ,
      y = 2 * sim.results$Ne * sim.results$cost * sim.results$sim.deltaR ,
      xlim = c(0, 1) ,
      pch = c(4, 20)[sim.results$costFactor] ,
      col = my.cols[sim.results$costFactor] ,
      cex = 2
    )
    legend(
      x = 0.5,
      y = 6,
      pch = c(4, 20) ,
      col = my.cols ,
      bty = 'n' ,
      title = 'Fitness cost' ,
      legend = these.costs ,
      pt.cex = 2
    )
    text(x = 0.05, y = 5.9,labels='A',cex = 2)
    par(op3)
  }
  
  ############
  #### 4B ####
  ############
  {
    #######################
    ### plotting params ###
    #######################
    {
      op4 <- par(mar = c(2.9, 0, 0, 0),
                 mgp = c(3, 0.6, 0))
      x.min <- -4
      x.max <- 8
      y.min <- 0
      y.max <- 1.4 * dnorm(0)
      lia.max <- 4
      my.x <- seq(x.min, lia.max , length.out = 1e5)
      y1 <- dnorm(my.x)
      thr <- qnorm(0.99)
      y1ft <- dnorm(thr)
      mean.shift <-
        uniroot(function(X)
          y1ft * these.costs[1] / these.costs[2] - dnorm(thr + X) ,
          interval = c(0, 1))$root
      y2 <- dnorm(my.x , mean = -mean.shift)
      aspect.ratio <- (y.max - y.min) / (x.max - x.min)
      box.hw <- 0.4
    }
    
    #################
    ### main plot ###
    #################
    {
      plot(
        NA,
        xlim = c(x.min, x.max) ,
        ylim = c(y.min, y.max) ,
        bty = 'n' ,
        xaxt = 'n' ,
        yaxt = 'n'
      )
      text(x = x.min + (x.max - x.min)*0.05, y = y.min + (y.max - y.min)*0.99,labels='B',cex = 2)
      axis(side = 1,
           at = seq(x.min, lia.max),
           labels = FALSE)
      mtext(
        side = 1 ,
        text = 'Liability' ,
        line = 1.2,
        cex = 0.7,
        at = 0
      )
      lines(x = my.x [y1 > y2] ,
            y = y1[y1 > y2] ,
            col = my.cols[1])
      lines(x = my.x [y1 <= y2] ,
            y = y1[y1 <= y2] ,
            col = adjustcolor(my.cols[1] , alpha.f = 0.3))
      polygon(
        x = c(my.x [y1 > y2], x.max , rev(my.x [y1 > y2])) ,
        y = c(y1[y1 > y2], 0 , rev(y2[y1 > y2])) ,
        col = adjustcolor(my.cols[1], alpha.f = 0.08),
        border = NA
      )
      lines(
        x = c(0, 0),
        y = c(0, dnorm(0,-mean.shift)),
        col = adjustcolor(my.cols[1], alpha.f = 0.3),
        lty = 2
      )
      lines(
        x = c(0, 0),
        y = c(dnorm(0,-mean.shift), dnorm(0) * 1.2),
        col = my.cols[1],
        lty = 2
      )
      lines(x = my.x ,
            y = y2 ,
            col = my.cols[2])
      polygon(
        x = c(x.min, my.x, x.max) ,
        y = c(0, y2, 0) ,
        col = adjustcolor(my.cols[2], alpha.f = 0.08),
        border = NA
      )
      lines(
        x = c(-mean.shift,-mean.shift),
        y = c(0, dnorm(0) * 1.2),
        col = my.cols[2],
        lty = 2
      )
      lines(x = rep(thr, 2),
            y = c(0, dnorm(0) * 0.8))
    }
    
    #################
    ### small box ###
    #################
    {
      lines(x = thr + c(-box.hw, box.hw) ,
            y = rep(2 * box.hw, 2) * aspect.ratio)
      lines(x = thr + rep(-box.hw, 2) ,
            y = c(0, 2 * box.hw) * aspect.ratio)
      lines(x = thr + rep(box.hw, 2) ,
            y = c(0, 2 * box.hw) * aspect.ratio)
      lines(x = thr + c(-box.hw, box.hw) ,
            y = rep(0, 2))
    }
    
    #################################
    ### parameters for bigger box ###
    #################################
    {
      vo <- 0.11
      ho <- 3.4
      mult <- 6.6
      inset.ymin <- vo
      inset.ymax <- 2 * box.hw * aspect.ratio * mult + vo
      inset.xmin <- ho + thr - box.hw * mult
      inset.xmax <- ho + thr + box.hw * mult
    }
    
    #####################
    ### zoom in lines ###
    #####################
    {
      lines(
        x = c(thr - box.hw, inset.xmin) ,
        y = c(2 * box.hw * aspect.ratio , inset.ymax) ,
        lty = 2 ,
        lwd = 0.8 ,
        col = adjustcolor('black' , alpha.f = 0.4)
      )
      lines(
        x = c(thr + box.hw, inset.xmax) ,
        y = c(0 , inset.ymin) ,
        lty = 2 ,
        lwd = 0.8 ,
        col = adjustcolor('black' , alpha.f = 0.4)
      )
    }
    
    ##############################
    ### plot inside bigger box ###
    ##############################
    {
      n.pt <- 100000
      plot.x <- seq(inset.xmin, inset.xmax, length.out = n.pt)
      silent.x1 <- seq(thr - box.hw, thr + box.hw, length.out = n.pt)
      silent.x2 <-
        seq(thr - box.hw + mean.shift,
            thr + box.hw + mean.shift,
            length.out = n.pt)
      new.y1 <- inset.ymin + dnorm(silent.x1) * mult
      new.y2 <- inset.ymin + dnorm(silent.x2) * mult
      x.in.max <- max(plot.x[new.y1 < inset.ymax])
      
      
      ## shading inside box
      polygon(
        x = c(plot.x[new.y1 < inset.ymax] ,
              rep(inset.xmax, 2) ,
              rev(plot.x[new.y2 < inset.ymax]),
              inset.xmin) ,
        y = c(new.y1[new.y1 < inset.ymax] ,
              c(tail(new.y1, 1), tail(new.y2, 1)),
              rev(new.y2[new.y2 < inset.ymax]),
              inset.ymax) ,
        col = adjustcolor(my.cols[1], alpha.f = 0.08) ,
        border = NA
      )
      
      polygon(
        x = c(plot.x ,
              rep(inset.xmax, 2) ,
              rev(plot.x)) ,
        y = c(new.y2 ,
              c(tail(new.y2, 1), inset.ymin),
              rep(inset.ymin, length(plot.x))) ,
        col = adjustcolor(my.cols[2], alpha.f = 0.08) ,
        border = NA
      )
      
      ## density lines inside box
      lines(x = plot.x[new.y1 < inset.ymax] ,
            y = new.y1[new.y1 < inset.ymax] ,
            col = my.cols[1])
      lines(x = plot.x[new.y2 < inset.ymax] ,
            y = new.y2[new.y2 < inset.ymax] ,
            col = my.cols[2])
      
      ## threshold inside box
      lines(
        x = rep(ho + thr, 2),
        y = c(inset.ymin, inset.ymax),
        lwd = 1
      )
    }
    
    #######################
    ### draw bigger box ###
    #######################
    {
      ## top
      lines(x = thr + c(-box.hw, box.hw) * mult + ho ,
            y = rep(inset.ymax, 2))
      ## left
      lines(x = thr + rep(-box.hw, 2) * mult + ho ,
            y = c(vo, inset.ymax))
      ## right
      lines(x = thr + rep(box.hw, 2) * mult + ho ,
            y = c(vo, inset.ymax))
      ## bottom
      lines(x = thr + c(-box.hw, box.hw) * mult + ho ,
            y = rep(vo, 2))
    }
    
    #####################################
    ### annotations inside bigger box ###
    #####################################
    {
      ## text inside bigger box
      {
        text(
          x = ho + thr - 1.4 ,
          y = inset.ymin + (inset.ymax - inset.ymin) * 0.135 ,
          labels = 'f(T|C=0.75)' ,
          col = my.cols[2] ,
          cex = 0.76
        )
        text(
          x = ho + thr - 1.25 ,
          y = inset.ymin + (inset.ymax - inset.ymin) * 0.67 ,
          labels = 'f(T|C=0.25)' ,
          col = my.cols[1] ,
          cex = 0.76
        )
      }
      
      ## lines marking threshold density heights
      {
        # high cost
        arrows(
          x0 = ho + thr - 0.2 ,
          x1 = ho + thr - 0.2 ,
          y0 = inset.ymin + 0.002 ,
          y1 = inset.ymin + dnorm(thr,-mean.shift) * mult ,
          col = my.cols[2] ,
          angle = 90 ,
          length = 0.022
        )
        arrows(
          x0 = ho + thr - 0.2 ,
          x1 = ho + thr - 0.2 ,
          y1 = inset.ymin + 0.002 ,
          y0 = inset.ymin + dnorm(thr,-mean.shift) * mult ,
          col = my.cols[2] ,
          angle = 90 ,
          length = 0.022
        )
        
        # low cost
        arrows(
          x0 = ho + thr + 0.2 ,
          x1 = ho + thr + 0.2 ,
          y0 = inset.ymin + 0.002 ,
          y1 = inset.ymin + dnorm(thr) * mult ,
          col = my.cols[1] ,
          angle = 90 ,
          length = 0.022
        )
        arrows(
          x0 = ho + thr + 0.2 ,
          x1 = ho + thr + 0.2 ,
          y1 = inset.ymin + 0.002 ,
          y0 = inset.ymin + dnorm(thr,-mean.shift) * mult ,
          col = my.cols[1] ,
          angle = 90 ,
          length = 0.022
        )
      }
    }
    
    ################################
    ### text labels in main plot ###
    ################################
    {
      text.label.ho <- 0.6
      text(
        x = ho + thr - 1.7 + text.label.ho ,
        y = inset.ymax + 0.1 ,
        labels = 'f(T|C=0.75)' ,
        col = my.cols[2] ,
        cex = 0.9
      )
      text(
        x = ho + thr - 1.7 + text.label.ho,
        y = inset.ymax + 0.07 ,
        labels = 'f(T|C=0.25)' ,
        col = my.cols[1] ,
        cex = 0.9
      )
      lines(
        x = ho + thr + c(-3.15,-0.3) + text.label.ho ,
        y = rep(inset.ymax + 0.086, 2)
      )
      
      text(
        x = ho + thr + text.label.ho,
        y = inset.ymax + 0.086 ,
        labels = '=' ,
        cex = 1
      )
      text(
        x = ho + thr + 0.9 + text.label.ho,
        y = inset.ymax + 0.1 ,
        labels = '0.25' ,
        col = my.cols[1] ,
        cex = 0.9
      )
      text(
        x = ho + thr + 0.9 + text.label.ho,
        y = inset.ymax + 0.07 ,
        labels = '0.75' ,
        col = my.cols[2] ,
        cex = 0.9
      )
      lines(
        x = ho + thr + 0.3 + c(0, 1.15) + text.label.ho,
        y = rep(inset.ymax + 0.086, 2)
      )
    }
    par(op4)
  }
  
  ############
  #### 4C ####
  ############
  solve.model.for.4C <- T
  op4 <- par(mar = c(2.9, 3.1, 0.1, 0.5),
             mgp = c(3, 0.6, 0))
  if (solve.model.for.4C) {
    source('scripts/solveTwoEffect.R')
    getSolnList <-
      function(Y, Z)
        sapply(Y, function(W)
          sapply(W, function(X)
            X[Z]))
    ba <- function(Y)
      ifelse(Y>50,1,(exp(2 * Y) - 1) / (exp(2 * Y) + 1))
    ya <- function(B) 0.5 * log((1 + B) / (1 - B))
    yt <- 1
    bt <- ba(yt)
    Ne = 20000
    C = 0.1
    yt <- ya(bt)
    ft <- yt / (2 * Ne * C)
    L <- c(1.5e6, 5e6, 1.5e7, 5e7, 1.5e8)
    u <- 1e-8
    as <- 1
    soln <- list()
    prev <- list()
    mean.nl <- list()
    deltal <- list()
    Vt <- list()
    al <- list()
    bs <- list()
    raw.ft <- list()
    h2l <- list()
    n.pts <- 500
    pl <- 0.0005
    max.nl <- 3
    
    rf <- F
    ## fixed h2
    for (j in seq_along(L)) {
      thetaL <- 4 * Ne * L[j] * u
      se.vg <- thetaL * as ^ 2 * bt / yt
      min.ya <- ya(1 - 1/(L[j]))
      min.al <- min.ya / (2 * Ne * ft * C)
      al[[j]] <-
        exp(seq(log(min.al), log(1100), length.out =
              1000))
      soln[[j]] <- list()
      last.bs <- NULL
      for (i in 1:length(al[[j]])) {
        if (i > 1)
          last.bs <- soln[[j]][[i - 1]]['bs']
        soln[[j]][[i]] <- solveTwoEffect3D(
          bt = bt,
          bs.guess = ifelse(i == 1, bt, last.bs),
          Ne = 20000,
          as = as,
          al = al[[j]][i],
          L = L[j],
          gs = 1 - pl,
          h2 = NULL ,
          Ve = se.vg ,
          u = u,
          C = C,
          Bval = 1,
          init.deltal = NULL,
          init.tstar = NULL,
          norm.deltal = NULL,
          norm.init.tstar = NULL
        )
      }
    }
    al.mat <- do.call(cbind, al)
    prev.mat <- getSolnList(soln, 'prev')
    deltal.mat <- getSolnList(soln, 'deltal')
    mean.nl.mat <- getSolnList(soln, 'mean.nl')
    Vt.mat <- getSolnList(soln, 'raw.Vt')
    raw.ft.mat <- getSolnList(soln, 'raw.ft')
    h2l.mat <- getSolnList(soln, 'h2l')
    bs.mat <- getSolnList(soln, 'bs')
    code.mat <- getSolnList(soln, 'code')
    se.approx.deltal.mat <- al.mat * raw.ft.mat
    
    
    ## two small effects solve
    small.al <- seq(1, 1000, length.out = 3000)
    # small.ft <- numeric()
    # for (i in 1:length(small.al)) {
    #   # FT <- seq(my.range[1],my.range[2],length.out=100)
    #   # bmean <- (1 - pl) * as * ba(2 * Ne * as * FT * C) + pl * small.al[i] * ba(2 * Ne * small.al[i] * FT * C) 
    #   # amean <- (1 - pl) * as + pl * small.al[i]
    #   # bt - bmean / amean
    #   small.ft[i] <- nleqslv(
    #     x = ft ,
    #     fn = function(FT){
    #         bmean <- (1 - pl) * as * ba(2 * Ne * as * FT * C) + pl * small.al[i] * ba(2 * Ne * small.al[i] * FT * C) 
    #         amean <- (1 - pl) * as + pl * small.al[i]
    #         return(bt - bmean / amean)
    #       }
    #   )$x
    # }
    # small.deltal <- small.al * small.ft
    
    small.deltal <- small.al*ft
    
  }
  {
    my.cols <- wes_palette('Darjeeling2')
    small.xlims <- c(0,40)
    big.xlims <- c(0,1000)
    small.ylims <- c(0, 0.01)
    big.ylims <- c(0, 0.26)
    plot(
      NA,
      xlim = big.xlims,
      ylim = big.ylims,
      bty = 'n',
      lty = 1 ,
      col = my.cols ,
      bty = 'n',
      xlab = '',
      ylab = '',
      xaxt = 'n'
    )
    text(x = big.xlims[1] + (big.xlims[2] - big.xlims[1])*0.03, y = big.ylims[1] + (big.ylims[2] - big.ylims[1])*0.99,labels='C',cex = 2)
    axis(side = 1, at = c(0,400,800), labels = c('0','400','800') )
    mtext(side = 1,
          text = expression(paste('Liability effect of large effect sites, ', a[L], sep = '')),
          line = 1.6,
          cex = 0.7
          )
    mtext(side = 2,
          text = expression(paste('Risk effect of large effect sites, ', delta^{R}, (a[L]), sep = '')),
          line = 1.6,
          cex = 0.7)
    for (i in 1:ncol(al.mat)) {
      plot.these <- code.mat[, i] == 1 & al.mat[,i] < 800 & mean.nl.mat[,i] < 3
      lines(al.mat[plot.these, i],
            deltal.mat[plot.these, i],
            col = my.cols[i])
    }
    
    ## dashed line marking extension of small effect linear solution
    # plot.these <- code.mat[, i] == 1 & al.mat[,i] < 1200
    # lines(al.mat[plot.these, i],
    #       se.approx.deltal.mat[plot.these, i],
    #       col = 'black' ,
    #       lty = 2)
    
    ## dotted line marking two effect solution assuming both are small
    lines(
      small.al ,
      small.deltal ,
      lty = 3,
      col = 'black'
    )
    
    #############
    ### inset ###
    #############
    ## draw little box
    aspect.ratio <- big.ylims[2]/big.xlims[2]
    box.width <- 40
    box.height <- aspect.ratio*box.width
    lines(x=c(0,box.width),y=rep(0,2),lty = 1 )
    lines(x=rep(0,2),y=c(0,box.height),lty = 1 )
    lines(x=c(0,box.width),y=rep(box.height,2),lty = 1 )
    lines(x=rep(box.width,2),y=c(0,box.height),lty = 1 )
    
    
    
    ## plot lines inside inset
    ho <- 570
    vo <- 0
    mult <- 11
    big.box.width <- box.width*mult
    big.box.height <- box.height*mult
    big.box.xmax <- ho + big.box.width
    big.box.ymax <- vo + big.box.height
    for (i in 1:ncol(al.mat)) {
      plot.these1 <- code.mat[, i] == 1 & al.mat[,i] < 40
      tmp.x <- ho + mult*al.mat[plot.these1, i]
      tmp.y <- vo + mult*deltal.mat[plot.these1, i]
      plot.these2 <- tmp.x < big.box.xmax & tmp.y < big.box.ymax
      plot.x <- tmp.x[plot.these2]
      plot.y <- tmp.y[plot.these2]
      lines(x = plot.x,
            y = plot.y,
            col = my.cols[i],
            lty = 1)
    }

    
    ## draw big box
    lines(x=ho+c(0,big.box.width),y=rep(vo,2),lty = 1 )
    lines(x=rep(ho,2),y=vo+c(0,big.box.height),lty = 1 )
    lines(x=ho+c(0,big.box.width),y=vo+rep(big.box.height,2),lty = 1 )
    lines(x=ho+rep(big.box.width,2),y=vo+c(0,big.box.height),lty = 1 )

    
    
    ## draw lines connecting boxes
    lines(
      x = c(box.width, ho) ,
      y = c(0, vo),
      lty = 2,
      col = adjustcolor('black' , alpha.f = 0.4)
    )
    lines(
      x = c(box.width, ho) ,
      y = c(box.height, vo+big.box.height),
      lty = 2,
      col = adjustcolor('black' , alpha.f = 0.4)
    )
    
    ## dashed line marking application of linearity assumption from 
    ## two effect solution to large effect sites
    # plot.these <- code.mat[, i] == 1 & al.mat[,i] < 40
    # lines(ho+mult*al.mat[plot.these, i],
    #       vo+mult*se.approx.deltal.mat[plot.these, i],
    #       col = 'black' ,
    #       lty = 2)
    
    ## dotted line marking two effect solution assuming both are small
    tmp.plot.y <- vo+mult*small.deltal
    tmp.plot.x <- ho+mult*small.al
    plot.y <- tmp.plot.y [ tmp.plot.y < big.box.ymax & tmp.plot.x < big.box.xmax ]
    plot.x <- tmp.plot.x [ tmp.plot.y < big.box.ymax & tmp.plot.x < big.box.xmax ]
    lines(
      plot.x ,
      plot.y ,
      lty = 3,
      col = 'black'
    )
    
    ## legend
    legend(
      x = 715 ,
      y = 0.192,
      legend = c(
        expression(1.5 %*% 10 ^ 6) ,
        expression(5 %*% 10 ^ 6) ,
        expression(1.5 %*% 10 ^ 7) ,
        expression(5 %*% 10 ^ 7) ,
        expression(1.5 %*% 10 ^ 8)
      ), 
      col = my.cols,
      lty = 1,
      title = 'Target size, L' ,
      cex = 0.7 ,
      bty = 'n'
    )
  }
  dev.off()
  
  ######################
  ## exploratory plot ##
  ######################
  exploratory.plots <- F
  if(exploratory.plots) {
    my.cols <- wes_palette("FantasticFox1")
    std.al.mat <- al.mat / sqrt(Vt.mat)
    op4 <- par(mar = c(3.1, 3.1, 0, 0),
               mgp = c(3, 0.6, 0))
    layout(t(matrix(1:4,nrow=2)))
    plot(
      NA ,
      xlim = c(0, 4) ,
      ylim = c(0, 1) ,
      type = 'l',
      lty = 1 ,
      col = my.cols ,
      bty = 'n' ,
      xlab = '',
      ylab = ''
    )
    mtext(side=1,text=expression(a[L]/sqrt(V[T])), line = 2)
    mtext(side=2,text=expression(delta[L]), line = 2)
    for(i in 1:ncol(std.al.mat)) {
      plot.these <- code.mat[,i]==1 & mean.nl.mat[,i]<max.nl
      lines(
        std.al.mat[plot.these,i],
        deltal.mat[plot.these==1,i], 
        col = my.cols[i]
      )
    }
    for(i in 1:ncol(std.al.mat)) {
      plot.these <- code.mat[,i]==1 & mean.nl.mat[,i]<max.nl
      lines(
        std.al.mat[plot.these,i],
        se.approx.deltal.mat[plot.these,i], 
        col = my.cols[i] ,
        lty = 2
      )
    }
    
    plot(
      NA ,
      xlim = c(0, 1500) ,
      ylim = c(0, 1) ,
      type = 'l',
      lty = 1 ,
      col = my.cols ,
      bty = 'n',
      xlab = '',
      ylab = ''
    )
    mtext(side=1,text=expression(a[L]), line = 2)
    mtext(side=2,text=expression(delta[L]), line = 2)
    for(i in 1:ncol(al.mat)) {
      plot.these <- code.mat[, i] == 1 & mean.nl.mat[, i] < max.nl
      lines(
        al.mat[code.mat[,i]==1,i],
        deltal.mat[code.mat[,i]==1,i], 
        col = my.cols[i]
      )
    }
    
    for(i in 1:ncol(al.mat)) {
      plot.these <- code.mat[,i]==1 & mean.nl.mat[,i]<max.nl
      lines(
        al.mat[code.mat[,i]==1,i],
        se.approx.deltal.mat[code.mat[,i]==1,i], 
        col = my.cols[i] ,
        lty = 2
      )
    }
    
    
    plot(
      NA ,
      xlim = c(0, 1000) ,
      ylim = c(0, 0.2) ,
      type = 'l',
      lty = 1 ,
      col = my.cols ,
      bty = 'n',
      xlab = '',
      ylab = ''
    )
    mtext(side=1,text=expression(a[L]), line = 2)
    mtext(side=2,text='Prevalence', line = 2)
    for(i in 1:ncol(al.mat)) {
      plot.these <- code.mat[,i]==1 & mean.nl.mat[,i]<max.nl
      lines(
        al.mat[code.mat[,i]==1,i],
        prev.mat[code.mat[,i]==1,i], 
        col = my.cols[i]
      )
    }
    
    plot(
      NA ,
      xlim = c(0, 1500) ,
      ylim = c(0, 1) ,
      type = 'l',
      lty = 1 ,
      col = my.cols ,
      bty = 'n',
      xlab = '',
      ylab = ''
    )
    mtext(side=1,text=expression(a[L]), line = 2)
    mtext(side=2,text=expression(h[L]^2), line = 2)
    for(i in 1:ncol(al.mat)) {
      plot.these <- code.mat[,i]==1 & mean.nl.mat[,i]<max.nl
      lines(
        al.mat[code.mat[,i]==1,i],
        h2l.mat[code.mat[,i]==1,i], 
        col = my.cols[i]
      )
    }
    legend(
      'bottomright' ,
      col = my.cols ,
      lty = 1 ,
      legend = L ,
      bty = 'n'
    )
  }
}


#### Figure 5 ####
{
  rm(list = ls())
  {
    source('scripts/figureFuncs.R')
    getSolnList <-
      function(Y, Z)
        sapply(Y, function(W)
          sapply(W, function(X)
            X[Z]))
    ba <- function(Y)
      ifelse(Y > 50, 1, (exp(2 * Y) - 1) / (exp(2 * Y) + 1))
    ya <- function(B)
      0.5 * log((1 + B) / (1 - B))
    library(wesanderson)
    my.cols <- wes_palette('Darjeeling2')
    png(
      'figures/paperFiguresForRealForRealThisTime/Figure5.png',
      height = 7 ,
      width = 17.75,
      units = 'cm',
      res = 600
    )
    layout(matrix(c(1, 2, 3), nrow = 1))
    
    #### 5A ####
    {
      op3 <- par(mar = c(3.1, 2.6, 0.4, 0.1),
                 mgp = c(3, 0.6, 0))
      my.bt <- seq(0.01, 0.999, by = 0.001)
      Ne <- 2e4
      u <- 1e-8
      h2 <- 0.5
      C <- 0.1
      L <- c(1.5e6, 5e6, 1.5e7, 5e7, 1.5e8)
      cparam.vals <- sqrt(L * u / (2 * Ne * h2)) / C
      bt.term <- sqrt(my.bt * log((1 + my.bt) / (1 - my.bt)))
      std.dens <- sapply(X = cparam.vals , function(X)
        X * bt.term)
      prevs <- 1 - pnorm(dnorminv(std.dens))
      matplot(
        x = my.bt ,
        y = prevs ,
        log = 'y' ,
        type = 'l' ,
        lty = 1 ,
        bty = 'n' ,
        ylim = c(1e-4, 1e-1) ,
        col = my.cols ,
        yaxt = 'n'
      )
      axis(side = 2 ,
           at = 10 ^ seq(-4,-1, 1))
      mtext(
        side = 1 ,
        text = expression(paste(
          'Relative threshold position, ' , b[T] , sep = ''
        )) ,
        line = 1.5,
        cex = 0.7,
      )
      mtext(
        side = 2 ,
        text = 'Prevalence' ,
        line = 1.5,
        cex = 0.7,
      )
      legend(
        'bottomright' ,
        lty = 1,
        col = my.cols ,
        legend = c(
          expression(1.5 %*% 10 ^ 6) ,
          expression(5 %*% 10 ^ 6) ,
          expression(1.5 %*% 10 ^ 7) ,
          expression(5 %*% 10 ^ 7) ,
          expression(1.5 %*% 10 ^ 8)
        ) ,
        bty = 'n' ,
        title = 'Target size, L'
      )
      legend(
        'topleft' ,
        legend = c(
          expression(N == 20000),
          expression(u == 10 ^ {
            -8
          }) ,
          expression(C == 0.1) ,
          expression(h ^ 2 == 0.5)
        ),
        bty = 'n'
      )
      par(op3)
    }
    
    #### 5B ####
    solve.model.for.B <- T
    if (solve.model.for.B) {
      {
        source('scripts/solveTwoEffect.R')
        bt <- 0.5
        Ne = 20000
        C = 0.1
        #ft <- log((1 + bt) / (1 - bt)) / (2 * Ne * C)
        L <- c(1.5e6, 5e6, 1.5e7, 5e7, 1.5e8)
        u <- 1e-8
        as <- 1
        al <- 100
        yt <- 0.5 * log((1 + bt) / (1 - bt))
        soln <- list()
        prev <- list()
        mean.nl <- list()
        bs <- list()
        raw.ft <- list()
        raw.Vt <- list()
        std.ft <- list()
        n.pts <- 200
        my.gs <- matrix(NA, nrow = n.pts + 1, ncol = length(L))
        max.pl <- 0.008
        max.nl <- 1
      }
      
      rf <- FALSE
      ## fixed Ve
      fixed.Ve <- numeric()
      for (j in seq_along(L)) {
        my.gs[, j] <-
          sort(unique(c(
            1 - max.pl / 2, seq (1 - 1 / L[j] , 1 - max.pl , length.out = n.pts)
          )),decreasing=TRUE)
        thetaL <- 4 * Ne * L[j] * u
        fixed.Ve[j] <- thetaL * as ^ 2 * bt / yt
        soln[[j]] <- list()
        for (i in 1:nrow(my.gs)) {
          if(i == 1 ) {
            init.bs <- NULL
            init.deltal <- NULL
            init.tstar <- NULL
          }
          soln[[j]][[i]] <- solveTwoEffect3D(
            bt = 0.5,
            Ne = 20000,
            as = as,
            al = al,
            L = L[j],
            gs = my.gs[i, j],
            #h2 = 1/2 ,
            Ve = fixed.Ve[j] ,
            u = u,
            C = C,
            Bval = 1,
            init.bs = init.bs,
            init.deltal = init.deltal,
            init.tstar = init.tstar,
            norm.deltal = NULL,
            norm.init.tstar = NULL
          )
          init.bs <- unname(soln[[1]][[1]]['bs'])
          init.deltal <- unname(soln[[1]][[1]]['deltal'])
          init.tstar <- unname(soln[[1]][[1]]['tstar'])
        }
        prev[[j]] <- sapply(soln[[j]], function(X)
          X['prev'])
        mean.nl[[j]] <- sapply(soln[[j]], function(X)
          X['mean.nl'])
      }
      fixed.h2.prev <- list()
      fixed.h2.mean.nl <- list()
      fixed.h2 = 1 / 2
      fixed.h2.soln <- list()
      for (j in seq_along(L)) {
        ## my.gs[, j] <- sort(unique(c(0.025,seq (1 - 1 / L[j] , 1-pl , length.out = n.pts))))
        thetaL <- 4 * Ne * L[j] * u
        se.vg <- thetaL * as ^ 2 * bt / yt
        fixed.h2.soln[[j]] <- list()
        i = 1
        for (i in 1:nrow(my.gs)) {
          if(i == 1 ) {
            init.bs <- NULL
            init.deltal <- NULL
            init.tstar <- NULL
          }
          fixed.h2.soln[[j]][[i]] <- solveTwoEffect3D(
            bt = 0.5,
            Ne = 20000,
            as = as,
            al = al,
            L = L[j],
            gs = my.gs[i, j],
            h2 = fixed.h2 ,
            Ve = NULL ,
            u = u,
            C = C,
            Bval = 1,
            init.bs = init.bs ,
            init.deltal = init.deltal ,
            init.tstar = init.tstar ,
            norm.deltal = NULL,
            norm.init.tstar = NULL
          )
          init.bs <- unname(soln[[1]][[1]]['bs'])
          init.deltal <- unname(soln[[1]][[1]]['deltal'])
          init.tstar <- unname(soln[[1]][[1]]['tstar'])
        }
        fixed.h2.prev[[j]] <- sapply(fixed.h2.soln[[j]], function(X)
          X['prev'])
        fixed.h2.mean.nl[[j]] <-
          sapply(fixed.h2.soln[[j]], function(X)
            X['mean.nl'])
      }
    }
    
    ## fixed h2 plot
    # {
    #   ## my.cols <- wes_palette('Darjeeling2')
    #   op3 <- par(mar = c(3.1, 2.6, 0.4, 0.1),
    #              mgp = c(3, 0.6, 0))
    #   y.max <-
    #     max(fixed.h2.prev[[length(L)]][fixed.h2.mean.nl[[length(L)]] < max.nl])
    #   plot(
    #     NA ,
    #     type = 'l',
    #     bty = 'n' ,
    #     xlim = c(0, max.pl) ,
    #     ylim = c(1e-3, 0.1) ,
    #     log = 'y' ,
    #     yaxt = 'n' ,
    #     xaxt = 'n'
    #   )
    #   axis(2,
    #        at = c(0.001, 0.01, 0.1),
    #        labels = c('0.001', '0.01', '0.1'))
    #   axis(1,
    #        at = c(0, 0.01, 0.02),
    #        labels = c('0', '0.01', '0.02'))
    #   mtext(
    #     side = 2 ,
    #     text = 'Prevalence' ,
    #     line = 1.5,
    #     cex = 0.7,
    #   )
    #   mtext(
    #     side = 1 ,
    #     text = 'Fraction of sites with large effects' ,
    #     line = 1.5,
    #     cex = 0.7,
    #   )
    #   for (j in seq_along(prev)) {
    #     lines(1 - my.gs[fixed.h2.mean.nl[[j]] < max.nl , j] ,
    #           fixed.h2.prev[[j]][fixed.h2.mean.nl[[j]] < max.nl] ,
    #           col = my.cols[j] ,
    #           lty = 1)
    #   }
    #   # for (j in seq_along(prev)) {
    #   #   lines(1 - my.gs[mean.nl[[j]] < max.nl , j] ,
    #   #         prev[[j]][mean.nl[[j]] < max.nl] ,
    #   #         col = my.cols[j] ,
    #   #         lty = 1)
    #   # }
    # }
    
    ## fixed Ve plot
    {
      ## my.cols <- wes_palette('Darjeeling2')
      op3 <- par(mar = c(3.1, 2.6, 0.4, 0.1),
                 mgp = c(3, 0.6, 0))
      y.max <-
        max(prev[[length(L)]][mean.nl[[length(L)]] < max.nl])
      plot(
        NA ,
        type = 'l',
        bty = 'n' ,
        xlim = c(0, max.pl) ,
        ylim = c(1e-3, 0.1) ,
        log = 'y' ,
        yaxt = 'n' ,
        xaxt = 'n'
      )
      axis(2,
           at = c(0.001, 0.01, 0.1),
           labels = c('0.001', '0.01', '0.1'))
      axis(1,
           at = c(0, 0.004, 0.008),
           labels = c('0', '0.004', '0.008'))
      mtext(
        side = 2 ,
        text = 'Prevalence' ,
        line = 1.5,
        cex = 0.7,
      )
      mtext(
        side = 1 ,
        text = 'Fraction of sites with large effects' ,
        line = 1.5,
        cex = 0.7,
      )
      for (j in seq_along(prev)) {
        lines(1 - my.gs[mean.nl[[j]] < max.nl , j] ,
              prev[[j]][mean.nl[[j]] < max.nl] ,
              col = my.cols[j] ,
              lty = 1)
      }
      # for (j in seq_along(prev)) {
      #   lines(1 - my.gs[mean.nl[[j]] < max.nl , j] ,
      #         prev[[j]][mean.nl[[j]] < max.nl] ,
      #         col = my.cols[j] ,
      #         lty = 1)
      # }
      
      ## add points for 5C
      j <- 2
      next.cols <- wes_palette('Rushmore1')[c(1, 3, 4)]
      my.pls <- c(1/2,1)*max(1-my.gs [mean.nl[[j]] < max.nl,j])
      #my.pls <- 1 - c(0.004, 0.008)
      these <- sapply(1-my.pls,function(X) which.min(abs(my.gs[,j] - X) ) )
      points(
        x = c(0,my.pls) ,
        y = c(prevs[my.bt==0.5,j],prev[[j]][these]) , 
        pch = 21 , 
        col ='black' ,
        bg = next.cols
      )
    }
    
    par(op3)
    #### 5C ####
    {
      j <- 2
      my.pls <- c(1/2,1)*max(1-my.gs [mean.nl[[j]] < max.nl,j])
      #my.pls <- 1 - c(0.004, 0.008)
      these <- sapply(1-my.pls,function(X) which.min(abs(my.gs[,j] - X) ) )
      #these <- my.gs[, j] %in% my.pls
      
      
      
      ## normal
      norm.var <- 4 * L[j] * Ne * u * bt / yt + fixed.Ve[j]
      ft.norm <- yt / (2 * Ne * C)
      std.ft.norm <- ft.norm * sqrt (norm.var)
      t.star.norm <- dnorminv(std.ft.norm) * sqrt (norm.var)
      my.xlims <- c(-40, 80)
      my.x <- seq(my.xlims[1], my.xlims[2], length.out = 100000)
      norm.dens <- dnorm(my.x, -t.star.norm, sqrt(norm.var))
      
      
      # pois conv
      these.solns <- fixed.h2.soln[[j]][these]
      
      
      op3 <- par(mar = c(3.1, 2.6, 0.4, 0.3),
                 mgp = c(3, 0.6, 0))
      my.cols <- wes_palette('Rushmore1')[c(1, 3, 4)]
      y.max <- ft.norm * 3
      plot(
        NA,
        xlim = my.xlims ,
        ylim = c(0, y.max),
        bty = 'n',
        xaxt = 'n',
        yaxt = 'n'
      )
      axis(side = 1,
           at = seq(-40, 80, by = 20),
           labels = FALSE)
      axis(side = 2,
           at = seq(0, 0.00035, by = 0.00005),
           labels = FALSE)
      pois.dens <- list()
      
      ## first pois conv
      t.star.pois <-
        these.solns[[2]]['tstar'] * sqrt(these.solns[[2]]['raw.Vt'])
      this.mean.nl <- these.solns[[2]]['mean.nl']
      al <- these.solns[[2]]['al']
      norm.sd <-
        sqrt(these.solns[[2]]['raw.Vas'] + these.solns[[2]]['raw.Ve'])
      pois.dens[[2]] <-
        sapply(my.x, function(X)
          dPoisConv(X + t.star.pois, this.mean.nl, norm.sd, al))
      lines(my.x[pois.dens[[2]] > norm.dens],
            pois.dens[[2]][pois.dens[[2]] > norm.dens],
            col = my.cols[3])
      lines(my.x[pois.dens[[2]] <= norm.dens],
            pois.dens[[2]][pois.dens[[2]] <= norm.dens],
            col = adjustcolor(my.cols[3], alpha.f = 0.3))
      
      ## second pois conv
      t.star.pois <-
        these.solns[[1]]['tstar'] * sqrt(these.solns[[1]]['raw.Vt'])
      this.mean.nl <- these.solns[[1]]['mean.nl']
      al <- these.solns[[1]]['al']
      norm.sd <-
        sqrt(these.solns[[1]]['raw.Vas'] + these.solns[[1]]['raw.Ve'])
      pois.dens[[1]] <-
        sapply(my.x, function(X)
          dPoisConv(X + t.star.pois, this.mean.nl, norm.sd, al))
      lines(my.x[pois.dens[[1]] > pois.dens[[2]]],
            pois.dens[[1]][pois.dens[[1]] > pois.dens[[2]]],
            col = my.cols[2])
      lines(my.x[pois.dens[[1]] <= pois.dens[[2]]],
            pois.dens[[1]][pois.dens[[1]] <= pois.dens[[2]]],
            col = adjustcolor(my.cols[2], alpha.f = 0.3))
      
      lines(x = my.x ,
            y = norm.dens,
            col = my.cols[1])
      
      norm.polygon.x <-
        c(my.xlims[1], my.x, my.xlims[2], my.xlims[1])
      norm.polygon.y <- c(y.max, norm.dens, 0, 0)
      polygon(
        norm.polygon.x ,
        norm.polygon.y ,
        col = adjustcolor(my.cols[1] , alpha.f = 0.08),
        border = NA
      )
      
      ## sliver between normal and lowest pois conv
      pois.polygon1.x <-
        c(my.x[pois.dens[[2]] > norm.dens], rev(my.x[pois.dens[[2]] > norm.dens]))
      pois.polygon1.y <-
        c(pois.dens[[2]][pois.dens[[2]] > norm.dens], rev(norm.dens[pois.dens[[2]] > norm.dens]))
      polygon(
        pois.polygon1.x ,
        pois.polygon1.y ,
        col = adjustcolor(my.cols[3] , alpha.f = 0.08),
        border = NA
      )
      
      
      ## sliver between the two pois conv
      pois.polygon2.x <-
        c(my.x[pois.dens[[1]] > pois.dens[[2]]], rev(my.x[pois.dens[[1]] > pois.dens[[2]]]))
      pois.polygon2.y <-
        c(pois.dens[[1]][pois.dens[[1]] > pois.dens[[2]]], rev(pois.dens[[2]][pois.dens[[1]] > pois.dens[[2]]]))
      polygon(
        pois.polygon2.x ,
        pois.polygon2.y ,
        col = adjustcolor(my.cols[2] , alpha.f = 0.08),
        border = NA
      )
      abline(v = 0 , lwd = 1.5)
    }
    dev.off()
    
  }
  
  #### Figure 5-1 ####
  png(
    'figures/paperFiguresForRealForRealThisTime/SuppFigure5-1.png',
    height = 7 ,
    width = 17.75,
    units = 'cm',
    res = 600
  )
  my.cols <- wes_palette('Darjeeling2')
  op3 <- par(mar = c(3.1, 2.6, 0.4, 0.1),
             mgp = c(3, 0.6, 0))
  layout(t(matrix(1:2,nrow=2)))
  my.pl <- 1 - my.gs
  my.deltal <- getSolnList(soln,'deltal')
  my.h2l <- getSolnList(soln,'h2l')
  plot(
    NA ,
    type = 'l',
    lty = 1 ,
    col = my.cols,
    xlim = c(0,max.pl) ,
    ylim = c(0,max(my.deltal)) ,
    xlab = '' ,
    ylab = '' ,
    bty = 'n'
  )
  mtext(side=1,text='Fraction of sites with large effects',line=1.5,cex=0.7)
  mtext(side=2,text=expression(paste('Risk effect of large effect sites, ', delta^{R}, (a[L]), sep = '')),line=1.5,cex=0.7)
  for( j in 1:ncol(my.deltal)){
    plot.these <- mean.nl[[j]] < max.nl
    lines(
      x = my.pl[plot.these,j] ,
      y = my.deltal[plot.these,j] ,
      col = my.cols[j]
    )
  }
  plot(
    NA ,
    type = 'l',
    lty = 1 ,
    col = my.cols,
    xlim = c(0,max.pl) ,
    ylim = c(0,1) ,
    xlab = '' ,
    ylab = '' ,
    bty = 'n'
  )
  mtext(side=1,text='Fraction of sites with large effects',line=1.5,cex=0.7)
  mtext(side=2,text='Heritability due to large effect sites',line=1.5,cex=0.7)
  for( j in 1:ncol(my.deltal)){
    plot.these <- mean.nl[[j]] < max.nl
    lines(
      x = my.pl[plot.these,j] ,
      y = my.h2l[plot.these,j] ,
      col = my.cols[j]
    )
  }
  legend(
    'topleft',
    legend = L, 
    col = my.cols,
    lty = 1,
    title = 'Target size, L' ,
    bty = 'n' ,
    cex = 0.7
  )
  dev.off()
  
  #### Figure 5-2 ####
  {
    my.cols <- wes_palette('Darjeeling2')
    ## 1 - bs
    my.bs <- getSolnList(soln,'bs')
    ## 2 - raw.ft
    my.ft <- getSolnList(soln,'raw.ft')
    ## 3 - raw.Vt
    my.Vt <- getSolnList(soln,'raw.Vt')
    ## 4 - std.ft
    my.std.ft <- getSolnList(soln,'std.ft')
    
    png(
      'figures/paperFiguresForRealForRealThisTime/SuppFigure5-2.png',
      height = 17.75 ,
      width = 17.75,
      units = 'cm',
      res = 600
    )
    
    layout(t(matrix(1:4,nrow=2)))
    op3 <- par(mar = c(3.1, 2.6, 0.4, 0.1),
               mgp = c(3, 0.6, 0))
    ## 5-2 A
    plot(
      NA ,
      type = 'l',
      lty = 1 ,
      col = my.cols,
      xlim = c(0,max.pl) ,
      ylim = c(0,0.5) ,
      xlab = '' ,
      ylab = '' ,
      bty = 'n'
    )
    mtext(side=1,text='Fraction of sites with large effects',line=1.5,cex=1)
    mtext(side=2,text='Fraction of small effect sites fixed for + allele',line=1.5,cex=1)
    lines(x = my.pl[, 1] ,
          y = my.bs[, 1])
    
    ## 5-2 B
    plot(
      NA ,
      type = 'l',
      lty = 1 ,
      col = my.cols,
      xlim = c(0,max.pl) ,
      ylim = c(0,max(unlist(my.ft))) ,
      xlab = '' ,
      ylab = '' ,
      bty = 'n'
    )
    mtext(side=1,text='Fraction of sites with large effects',line=1.5,cex=1)
    mtext(side=2,text='Raw threshold density',line=1.5,cex=1)
    for( j in 1:ncol(my.ft)){
    ##  plot.these <- mean.nl[[j]] < max.nl
      lines(
        x = my.pl[,j] ,
        y = my.ft[,j]
      )
    }
    
    ## 5-2 C
    plot(
      NA ,
      type = 'l',
      lty = 1 ,
      col = my.cols,
      xlim = c(0,max.pl) ,
      ylim = c(0,max(sqrt(unlist(my.Vt)))) ,
      xlab = '' ,
      ylab = '' ,
      bty = 'n'
    )
    mtext(side=1,text='Fraction of sites with large effects',line=1.5,cex=1)
    mtext(side=2,text='Std. dev. of liability',line=1.5,cex=1)
    for( j in 1:ncol(my.Vt)){
      ##  plot.these <- mean.nl[[j]] < max.nl
      lines(
        x = my.pl[,j] ,
        y = sqrt(my.Vt[,j]) ,
        col = my.cols[j]
      )
    }
    
    ## 5-2 D
    plot(
      NA ,
      type = 'l',
      lty = 1 ,
      col = my.cols,
      xlim = c(0,max.pl) ,
      ylim = c(0,max(unlist(my.std.ft))) ,
      xlab = '' ,
      ylab = '' ,
      bty = 'n'
    )
    mtext(side=1,text='Fraction of sites with large effects',line=1.5,cex=1)
    mtext(side=2,text='Standardized threshold density',line=1.5,cex=1)
    for( j in 1:ncol(my.Vt)){
      ##  plot.these <- mean.nl[[j]] < max.nl
      lines(
        x = my.pl[,j] ,
        y = my.std.ft[,j] ,
        col = my.cols[j]
      )
    }
    dev.off()
  }
    
  #### Figure 5-3 ####
  my.amean <- getSolnList(soln,'amean')
  mean.std.ft <- my.ft / my.amean
  mean.std.std.dev <- sqrt(my.Vt)/my.amean

  matplot(my.pl,my.ft,col = my.cols,lty =1,type ='l')  
  matplot(my.pl,mean.std.ft,col = my.cols,lty =1,type ='l')  
  matplot(my.pl,mean.std.std.dev,col = my.cols,lty =1,type ='l')  
  
}






##################
#### Figure 6 ####
##################
my.dir <- 'figures/smilePlotData'
my.filenames <- list.files(my.dir)
my.data <-
  lapply(my.filenames, function(x)
    read.delim(paste(my.dir, x, sep = '/'), sep = '\t', header = T))
my.cols <- wes_palette("AsteroidCity3")[c(2,4)]

{
  png(
    'figures/paperFiguresForRealForRealThisTime/Figure6.png',
    height = 17.75 ,
    width = 17.75,
    units = 'cm',
    res = 600
  )
  layout(matrix(c(1,1,1, 2, 3, 6,4,5,6), nrow = 3),widths = c(1,12,12),heights = c(
  12,12,1))
  
  par(mar=c(0,0,0,0))
  plot.new()
  text(x=0.5,y=0.5,labels='Effect Size',srt=90,cex = 2.5)
  op6 <- par(mar = c(2.9, 3, 2, 0),
             mgp = c(3, 0.6, 0))
  plot(
    x = my.data[[1]][, 1] ,
    y = my.data[[1]][, 2],
    ylim = c(0.005,0.05),
    xlim = c(0,1),
    pch = 21 ,
    bty = 'n',
    xlab = '',
    ylab = '' ,
    log = 'y' ,
    main = 'CVD' ,
    bg = my.cols[1] ,
    col = my.cols[2]
  )
  
  
  
  plot(
    x = my.data[[3]][, 1] ,
    y = my.data[[3]][, 2],
    ylim = c(0.05,0.5),
    xlim = c(0,1),
    pch = 21 ,
    bty = 'n',
    xlab = '',
    ylab = '' ,
    log = 'y',
    main = 'IBD',
    bg = my.cols[1] ,
    col = my.cols[2]
  )
  
  
  plot(
    x = my.data[[4]][, 1] ,
    y = my.data[[4]][, 2],
    ylim = c(0.045,0.25),
    xlim = c(0,1),
    pch = 21 ,
    bty = 'n',
    xlab = '',
    ylab = '' ,
    log = 'y',
    main = 'SCZ',
    bg = my.cols[1] ,
    col = my.cols[2]
  )
  
  
  plot(
    x = my.data[[5]][, 1] ,
    y = my.data[[5]][, 2],
    ylim = c(0.035,0.3),
    xlim = c(0,1),
    pch = 21 ,
    bty = 'n',
    xlab = '',
    ylab = '' ,
    log = 'y',
    yaxt = 'n',
    main = 'T2D',
    bg = my.cols[1] ,
    col = my.cols[2]
  )
  axis(side = 2, at = c(0.04,0.1,0.25))
  par(op6)
  par(mar=c(0,0,0,0))
  plot.new()
  text(x=0.5,y=0.5,labels='Frequency of risk increasing allele',cex = 2.5)
  dev.off()
}

