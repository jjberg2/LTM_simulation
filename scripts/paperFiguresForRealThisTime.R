#setwd('/Users/jeremyberg/Documents/academics/LTM_simulation/')

##########################
### Figure 1 schematic ###
##########################
{
  h2 <- 0.6
  source('scripts/figureFuncs.R')
  png(
    'figures/paperFiguresForRealThisTime/SchematicFigure1.png',
    height = 18 ,
    width = 22,
    units = 'cm',
    res = 500
  )
  nf <- layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))
  
  ## 1A
  ## mid zoom liability distribution
  op1 <- par(mar = c(3, 3, 2, 4))
  makePhenLi(xlim = c(-3, 4), prev = 0.01)
  mtext(side = 1 ,
        text = 'Liability',
        line = 1.5)
  mtext(side = 2 ,
        text = 'Density on Liability',
        line = 1.5)
  mtext(
    side = 3 ,
    text = 'A' ,
    line = 0,
    at = -3,
    cex = 1.5
  )
  par(op1)
  
  
  ## 1B
  this.col = 'darkorchid4'
  op2 <- par(mar = c(3, 4, 2, 4))
  this.ymax <- makeGenLi(xlim = c(-3, 4), h2 = 0.6)
  abline(v = qnorm(1-0.01),lty = 2 )
  addRiskCurve(
    h2 = h2,
    prev = 0.01,
    y.max = this.ymax,
    my.col = this.col
  )
  axis(
    side = 4,
    at = seq(0, this.ymax, by = 0.2 * this.ymax),
    labels = seq(0, 1, by = 0.2),
    las = 1
  )
  mtext(side = 1 ,
        text = 'Genetic Liability',
        line = 1.5)
  mtext(side = 2 ,
        text = 'Density on Genetic Liability',
        line = 2.2)
  mtext(
    side = 4 ,
    text = 'Genetic Risk',
    line = 2.5,
    col = this.col,
    las = 0
  )
  mtext(
    side = 3 ,
    text = 'B' ,
    line = 0,
    at = -3,
    cex = 1.5
  )
  par(op2)
  
  
  prev <- 0.01
  t.pos <- qnorm(1 - prev, mean = 0, sd = 1)
  my.risks <- exp(seq (log(1e-9) , log(1 - 1e-6) , length.out = 1e5))
  my.g <- qnorm(my.risks , mean = t.pos , sd = sqrt (1 - h2))
  pg <- dnorm(my.g, 0, h2)
  dgdr <-
    1 / dnorm (qnorm (my.risks, t.pos , sd = sqrt (1 - h2)) , 0, h2)
  pr <- pg / dgdr
  plot.risks <- c(0, my.risks, 0)
  plot.pr <- c(0, pr, 0)
  
  this.col = 'darkorchid4'
  op3 <- par(mar = c(4, 4, 1, 4))
  plot(
    NA,
    type = 'l',
    xlim = c (0 , 0.2),
    ylim = c(0, 0.5),
    ylab = '',
    xlab = '',
    bty = 'n'
  )
  mtext(side = 1 ,
        text = 'Genetic Risk',
        line = 2.2)
  mtext(side = 2 ,
        text = 'Density on Genetic Risk',
        line = 2.2)
  polygon(
    x = plot.risks ,
    y = plot.pr,
    col = adjustcolor(this.col, alpha.f = 0.05),
    border = NA
  )
  polygon(
    x = plot.risks ,
    y = plot.pr,
    col = adjustcolor(this.col, alpha.f = 0.4),
    border = NA,
    density = 120 ,
    angle = 315
  )
  mtext(
    side = 3 ,
    text = 'C' ,
    line = 0,
    at = 0,
    cex = 1.5
  )
  par(op3)
  
  this.col = 'firebrick4'
  op4 <- par(mar = c(4, 3, 1, 4))
  this.ymax <- makeGenLi(xlim = c(-3, 4), h2 = 0.6)
  abline(v = qnorm(1-0.01),lty = 2 )
  addRiskCurve(
    h2 = h2,
    prev = 0.01,
    y.max = this.ymax,
    my.col = this.col,
    fit.surface = TRUE
  )
  mtext(side = 1 ,
        text = 'Genetic Liability',
        line = 2.2)
  mtext(side = 2 ,
        text = 'Density on Genetic Liability',
        line = 1.5)
  mtext(
    side = 4 ,
    text = 'Expected fitness',
    line = 2.5,
    col = this.col,
    las = 0
  )
  axis(
    side = 4,
    at = seq(0, this.ymax, by = 0.2 * this.ymax),
    labels = seq(0, 1, by = 0.2),
    las = 1
  )
  mtext(
    side = 3 ,
    text = 'D' ,
    line = 0,
    at = -3,
    cex = 1.5
  )
  par(op4)
  dev.off()
}

##########################
### Figure 2 schematic ###
##########################
{
  png(
    'figures/paperFiguresForRealThisTime/SchematicFigure2.png',
    height = 18 ,
    width = 22,
    units = 'cm',
    res = 500
  )
  op1 <- par(mar = c(3.3, 3.3, 2, 0) + 0.1)
  source('scripts/figureFuncs.R')
  y.max <- 0.06
  t.pos <- qnorm(0.99)
  alpha <- 0.15
  my.out <- makePhenWEffect(
    xlim = c(2, 2.6),
    y.max = y.max,
    prev = 0.1,
    alpha = alpha,
    t.pos = t.pos,
    xaxt = 'n',
    return.stuff = TRUE
  )
  effect.x <- my.out[[3]]
  effect.y <- my.out[[4]]
  axis(side = 1 , labels = FALSE)
  ##axis(side = 2 , labels = FALSE)
  ## polygon(
  ##   x = c(t.pos - alpha , t.pos - alpha , t.pos , t.pos) ,
  ##   y = c(0 , dnorm(t.pos), dnorm(t.pos), 0),
  ##   col = rgb(1, 0, 1, 0.8),
  ##   border = NA
  ## )
  x.pos <- 2.58
  y.pos <- 0.03
  div <- 4
  polygon(
    x = c(x.pos - alpha / div , x.pos - alpha / div , x.pos , x.pos) ,
    y = c(
      y.pos ,
      y.pos + dnorm(t.pos) / div,
      y.pos + dnorm(t.pos) / div,
      y.pos
    ),
    col = rgb(1, 0, 1, 0.8),
    border = NA
  )
  x.offset <- x.pos - alpha / div - effect.x[1]
  y.offset <- 0.04
  min.x <- min(effect.x)
  min.y <- min(effect.y)
  small.effect.x <- x.offset + (effect.x - min.x)/4 + min.x
  small.effect.y <- y.offset + (effect.y - min.y)/4 + min.y
  polygon(
    x = small.effect.x ,
    y = small.effect.y,
    col = rgb(1, 0, 1, 0.05),
    border = NA
  )
  polygon(
    x = small.effect.x ,
    y = small.effect.y,
    col = rgb(1, 0, 1, 0.4),
    border = NA,
    angle = 0 ,
    density = 120
  )
  text(
    x = x.pos - 0.09,
    y = y.pos + 0.0026 ,
    labels = 'a f(T) = ',
    cex = 2.5
  )
  text(
    x = x.pos - 0.09,
    y = y.offset + 0.0026 ,
    labels = expression(paste(delta(a) ,' = ', sep = '')),
    cex = 2.5
  )
  mtext("T",
        at = t.pos,
        cex = 3,
        family = 'Arial')
  mtext("T-a",
        at = t.pos - alpha,
        cex = 3,
        family = 'Arial')
  mtext('Liability' ,
        side = 1,
        line = 2,
        cex = 3)
  lines(
    x = c(2, t.pos),
    y = rep(dnorm(t.pos), 2),
    lty = 3,
    lwd = 3
  )
  mtext(
    side = 2 ,
    at = dnorm(t.pos),
    text = 'f(T)',
    line = -1,
    las = 1,
    cex = 3
  )
  ##text(x = t.pos - alpha/2, y = y.max / 2.7, '{', srt = 270, cex = 8, family = 'Helvetica Neue UltraLight')
  par(op1)
  dev.off()
}

################################
### Figure 3 equlibrium bias ###
################################
{
png('figures/paperFiguresForRealThisTime/AsymmetryFigure3.png', height = 10 , width = 11, units = 'cm', res = 500)
op3 <- par(mar=c(3.3,3.3,2,0.4)+0.1)
op4 <- options(scipen=400)
my.gamma <- 10^seq(-3,3,length.out=1000)
bias <- ifelse(my.gamma<700,(exp(my.gamma)-1)/(exp(my.gamma)+1),1)
plot(
    x = my.gamma ,
    y = bias ,
    log = 'x' ,
    type = 'n' ,
    bty = 'n' ,
    lwd = 2 ,
    xlab = '' ,
    ylab = '',
    xaxt = 'n' 
)
mtext(text='Fixation Asymmetry', side = 2, line = 2)
mtext(text='Population Scaled Selection Coefficient', side = 1, line = 2)
axis(side = 1 , at = c(0.001,0.01,0.1,1,10,100,1000),labels = c('0.001','0.01','0.1','1','10','100','1000'))
polygon(
    x=c(1e1,1e1,1e3,1e3),
    y=c(0,1,1,0),
    border=NA,
    col=adjustcolor('grey',alpha.f=0.6)
)
polygon(
    x=c(1e-1,1e-1,1e-3,1e-3),
    y=c(0,1,1,0),
    border=NA,
    col=adjustcolor('grey',alpha.f=0.6)
)
##mtext(text='Effectively',side = 3 , line = 1,at = 0.01)
##mtext(text='neutral',side = 3 , line = 0,at = 0.01)
text(labels='Effectively', x = 0.01, y=0.2)
text(labels='neutral', x = 0.01, y=0.1)
##mtext(text='Weakly',side = 3 , line = 1,at = 1)
##mtext(text='selected',side = 3 , line = 0,at = 1)
text(labels='Weakly', x = 0.4, y=0.9)
text(labels='selected', x = 0.4, y=0.8)
##mtext(text='Strongly',side = 3 , line = 1,at = 100)
##mtext(text='selected',side = 3 , line = 0,at = 100)
text(labels='Strongly', x = 100, y=0.9)
text(labels='selected', x = 100, y=0.8)
lines(
    x=my.gamma,
    y=bias,
    lwd=2
)


par(op3)
options(op4)
dev.off()
}



png('figures/paperFiguresForRealThisTime/MeanEffectFigure4.png', height = 10 , width = 11, units = 'cm', res = 500)
op4 <- par(mar=c(3.3,3.3,2,0.4)+0.1)
my.gamma <- 10^seq(-2,1,length.out=1000)
bias <- (exp(my.gamma)-1)/(exp(my.gamma)+1)
plot(
    x=bias,
    y=my.gamma,
    type='n',
    bty='n',
    lwd=2,
    xlab='',
    ylab=''
)
mtext(expression(paste('Relative threshold position, ', b[T] ,sep = '')),side = 1 , line = 2)
mtext(expression(paste('Population scaled selection coefficient, ', gamma, '(a)', sep = '')),side = 2 , line = 2)
lines(
    x=bias,
    y=my.gamma,
    lwd=2
)
par(op4)
dev.off()





#########################################
## Figure 4 liability variance scaling ##
#########################################

## bt = 0.6479407
## bt*ay = 1


rm(list = ls())
library("MetBrewer")
library('numDeriv')
source('scripts/simpleVarFuncs.R')
dnorminv<-function(y) sqrt(-2*log(sqrt(2*pi)*y))
prevs = c(0.03, 0.01, 0.003, 0.001)
my.cols = met.brewer('Hiroshige', length(prevs))
thr = sapply(prevs, function(PREV)
  qnorm(1 - PREV))
phits = dnorm(thr)
max.a = 5
min.a = rep(0.001, 4)
my.a = mapply(
  function(MAXA, MINA)
    seq(MINA, MAXA, length.out = 10000),
  MINA = min.a,
  MAXA = max.a,
  SIMPLIFY = FALSE
)
Ne = 10000
L = 1e7
u = 1e-8
U = 2*L*u
aybt <- 1
thetaU = 4*Ne*U

this.cost <- 1/10
h2 = 1/10

tmp <- numeric()
for(i in 1:length(prevs)){
    tmp[i] <- 1/uniroot(function(X) prevs[i] - pnorm(dnorminv(sqrt(2*L*u*aybt / (2*Ne*(1/X)))/(1/X)),lower.tail = FALSE), lower = 1,upper=20)$root
}

## my.ay <- lapply(my.a,function(A) A*sqrt(thetaU / h2))
my.ay <- mapply(function(A,PHIT) 2*Ne*this.cost*PHIT*A,
                A = my.a,
                PHIT = phits,
                SIMPLIFY=FALSE)
norm.del <- mapply(function(THR, A)
  normDel(A, THR),
  THR = thr,
  A = my.a)
li.vars <- mapply(function(THR, AS,AY)
  sigmaa2(AS, AY, THR, Ne * this.cost),
  THR = thr,
  AS = my.a,
  AY = my.ay)
risk.vars <- mapply(function(THR, A)
  sigmaRisk(A, THR, Ne * this.cost),
  THR = thr,
  A = my.a)
small.ay.seq <- seq(0.01,2000,length.out=10000)
small.approx.vars = sigmaaSmall2(small.ay.seq)
plot.a <- do.call (cbind, my.a)
plot.ay <- do.call (cbind, my.ay)
my.cex <- 1.5

png(
    'figures/paperFiguresForRealThisTime/effectSizeVarianceRelationshipFigure4.png',
    units = 'in',
    res = 400,
  width = 12,
  height = 6
)
par(mfrow = c(1, 2),mar = c(5,4.5,4,1)+0.1)
matplot(
  x = plot.ay ,
  y = li.vars,
  xlim = c(0, 6),
  ylim = c(0, 12),
  type = 'l',
  lty = 1 ,
  lwd = 3 ,
  col = my.cols,
  xlab = '' ,
  ylab = 'Liability variance',
  cex.lab = my.cex,
  cex.axis = my.cex
)
legend(
  'topleft',
  legend = prevs,
  lty = 1,
  lwd = 3,
  col = my.cols,
  bty = 'n',
  title = "Prevalence",
  cex = 1
)
legend(
  x = -1/5, 
  y = 9,
  legend = c('Small effect approx','y = 2x'),
  lty = c(2,3),
  lwd = c(3,1),
  col = 'black',
  bty = 'n',
  cex = 1
)
mtext(
  text = 'Liability effect size',
  side = 1,
  line = 3,
  cex = my.cex
)
## mtext(
##   text = '(Population scaled fitness units)',
##   side = 1, 
##   line = 3,
##   cex = my.cex
## )
curve(
  sigmaaSmall2,
  from = 0 , 
  to = 20,
  type = 'l',
  lty = 2 , 
  lwd = 3 ,
  col = 'black',
  add = T
)
abline(
  a = 0,
  b = 2,
  lty = 3 ,
  lwd = 1)
mtext('A', side = 3 , at = 0, cex = 3)
matplot(
  x = plot.ay ,
  y = li.vars,
  xlim = c(0, 1000),
  ylim = c(0, 2000),
  type = 'l',
  lty = 1 ,
  lwd = 3 ,
  col = my.cols,
  xlab = '' ,
  ylab = 'Liability variance',
  cex.lab = my.cex,
  cex.axis = my.cex
)
mtext(
  text = 'Liability effect size',
  side = 1,
  line = 3,
  cex = my.cex
)
## mtext(
##   text = '(Population scaled fitness units)',
##   side = 1, 
##   line = 3,
##   cex = my.cex
## )
curve(
  sigmaaSmall2,
  from = 0 , 
  to = 1000,
  type = 'l',
  lty = 2 , 
  lwd = 3 ,
  col = 'black',
  add = T
)
abline(
  a = 0,
  b = 2,
  lty = 3 ,
  lwd = 1)
mtext('B', side = 3 , at = 0, cex = 3)
dev.off()






###########################
### Figure 5 var dist ###
###########################
{
  library('wesanderson')
  n <- 10000
  my.x <- 1:(n-1)/n
  my.gamma <- c(0.3,1,3,10)
  ##my.var.specs <- sapply(my.gamma, function(G) 2*G*(exp(2*G) + 1) /  exp(2*G)* exp(-2*G*my.x))
  my.var.specs <- sapply(my.gamma, function(G) (1-exp(-2*G*my.x)) / (1-exp(-2*G)) )
  my.der.probs <- sapply(my.gamma, function(G) ( exp(2*G) - exp(2*G*my.x) ) / (exp(2*G) - 1 ) )
  cex.lab = 1.3
  cex.axis = 1.2
  png(
    'figures/paperFiguresForRealThisTime/SegAsymmetryFigure5.png',
    height = 10 ,
    width = 22,
    units = 'cm',
    res = 500
  )
  par(mfrow = c(1, 2))
  op3 <- par(mar = c(3.3, 3.3, 2, 0.4) + 0.1)
  my.cols <- wes_palette('GrandBudapest1', length(my.gamma))
  plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, 1),
    bty = 'n',
    ylab = '',
    xlab = ''
  )
  mtext(
    text = 'A',
    side = 3,
    line = 0.5,
    at = 0,
    cex = cex.lab
  )
  mtext(
    text = 'Proportion of Variance',
    side = 2,
    line = 2.1,
    cex = cex.lab
  )
  mtext(
    text = 'Frequency of Risk Increasing Allele',
    side = 1,
    line = 2.1,
    cex = cex.lab
  )
  abline(a = 0 , b = 1 , lty = 2 )
  matplot(
    my.x,
    my.var.specs,
    type = 'l',
    lty = 1 ,
    lwd = 2 ,
    add = TRUE,
    col = my.cols
  )
  legend(
    'bottomright',
    legend = c(0,my.gamma) ,
    col = c('black',my.cols),
    lty = c(2,rep(1,times=length(my.cols)) ),
    lwd = c(1,rep(2,times=length(my.cols)) ) ,
    bty = 'n',
    title = 'Scaled selection coefficient' ,
    cex = 0.8
  )
  par(op3)
  op4 <- par(mar = c(3.3, 4.3, 2, 0.4) + 0.1)
  plot(NA,
       xlim = c(0, 1),
       ylim = c(0, 1),
       bty = 'n',
       xlab = '',
       ylab = '')
  abline(a = 1 , b = -1 , lty = 2 )
  matplot(
    my.x,
    my.der.probs,
    type = 'l',
    lty = 1 ,
    lwd = 2 ,
    add = TRUE,
    col = my.cols
  )
  mtext(
    text = 'B',
    side = 3,
    line = 0.5,
    at = 0,
    cex = cex.lab
  )
  mtext(
    text = 'Probability risk increasing',
    side = 2,
    line = 3.1,
    cex = cex.lab
  )
  mtext(
    text = 'allele is derived',
    side = 2,
    line = 2.1,
    cex = cex.lab
  )
  mtext(
    text = 'Frequency of risk increasing allele',
    side = 1,
    line = 2.1,
    cex = cex.lab
  )
  par(op4)
  dev.off()
}










########################################
### Figure 6 mean scaled coefficient ###
########################################
{
    library(wesanderson)
  balpha = function(a)
    ifelse(a < 100, (exp(2 * a) - 1) / (exp(2 * a) + 1), 1)
  bt.diff2 = function(my.mean, my.cv, my.bt,my.q=1e-8) {
    my.shape = 1 / my.cv ^ 2
    my.rate = 1 / (my.cv ^ 2 * my.mean)
    my.lower = qgamma(my.q, my.shape, my.rate)
    my.upper = qgamma(1 - my.q, my.shape, my.rate)
    ## balpha(my.mean)
    myFunc <- function(a)
      dgamma(a, my.shape, my.rate) * a * balpha(a)
    my.bt - (1 / my.mean) * integrate(myFunc, lower = my.lower, upper = my.upper)$value
  }
  my.bts = seq(1e-3, 0.99999, length.out = 1000)
  my.cvs = sqrt(c(0.1, 1, 3, 10))
  se.gammas = 0.5 * log((1 + my.bts) / (1 - my.bts))
  Ne = 1e4
  U = 0.2
  h2 = 1 / 2
  cost = 1 / 2
  tmp = list()
  my.gammas = list()
  for (j in 1:length(my.cvs)) {
    tmp[[j]] = list()
    my.gammas[[j]] = rep(NA, length(my.bts))
    for (i in 1:length(my.bts)) {
      my.lower = 0.01 * se.gammas[i]
      my.upper =  2000 * se.gammas[i]
      tmp[[j]][[i]] = uniroot(
        function(X)
          bt.diff2(X, my.cv = my.cvs[j], my.bts[i]),
        lower = my.lower,
        upper = my.upper
      )
      if (tmp[[j]][[i]]$root < 0.002)
        next
      my.gammas[[j]][i] = tmp[[j]][[i]]$root
    }
  }
  {
    cex.lab = 1.3
    cex.axis = 1.2
    png(
      'figures/paperFiguresForRealThisTime/MeanCoefficientFigure6.png',
      height = 10 ,
      width = 22,
      units = 'cm',
      res = 500
    )
    layout(matrix(c(1, 2), nrow = 1))
    op = par(mar = c(3.6, 4.6, 1, 0.5) + 0.1)
    ## 1
    plot(
      x = my.bts,
      y = se.gammas,
      type = 'l' ,
      xlab = '',
      ylab = '',
      cex.axis = cex.axis,
      cex.lab = cex.lab,
      ylim = c(0, 6),
      lwd = 2,
      lty = 1,
      bty = 'n'
    )
      mtext(
    text = 'A',
    side = 3,
    line = 0,
    at = 0,
    cex = cex.lab
  )
    mtext(
      side = 2,
      text = 'Mean population scaled',
      cex = cex.lab,
      line = 3.5
    )
    mtext(
      side = 2,
      text = expression(paste(
        'selection coefficient, ', bar(gamma(alpha)), sep = ''
      )),
      cex = cex.lab,
      line = 1.9
    )
    mtext(
      side = 1,
      text = expression(paste(
        'Relative threshold position, ', b[T], sep = ''
      )),
      cex = cex.lab,
      line = 2.4
    )
    my.cols <- wes_palette('GrandBudapest1', length(my.gammas))
    # for (j in 1:length(my.gammas)) {
    #   plot.these = !is.na(my.gammas[[j]]) & my.gammas[[j]] < 10
    #   lines(
    #     x = my.bts[plot.these],
    #     y = my.gammas[[j]][plot.these],
    #     lty = 2,
    #     lwd = 2,
    #     col = my.cols[j]
    #   )
    # }
    ##abline(v = bt.inflect, lty = 2 , lwd = 2)
    par(op)
    op = par(mar = c(3.6, 4.6, 1, 0.5) + 0.1)
    plot(
      NA,
      xlab = '',
      ylab = '',
      cex.axis = cex.axis,
      cex.lab = cex.lab,
      xlim = c(0, 1),
      ylim = c(0, 1.5),
      lwd = 2,
      lty = 1,
      bty = 'n'
    )
    mtext(
      side = 2,
      text = expression(paste(
        'Fold change in ', bar(gamma(alpha)) , ' relative', sep = ''
      )),
      cex = cex.lab,
      line = 3.3
    )
      mtext(
    text = 'B',
    side = 3,
    line = 0,
    at = 0,
    cex = cex.lab
  )
    mtext(
      side = 2,
      text = 'to single effect model',
      cex = cex.lab,
      line = 2.1
    )
    mtext(
      side = 1,
      text = expression(paste(
        'Relative threshold position, ', b[T], sep = ''
      )),
      cex = cex.lab,
      line = 2.4
    )
    for (i in seq_along(my.cvs)) {
      abline(
        a = 1 / (1 + my.cvs[i] ^ 2) ,
        b = 0,
        col = my.cols[i],
        lwd = 0.3,
        lty = 3
      )
    }
    abline(a = 1 , b = 0 , lty = 2)
    for (j in 1:length(my.gammas)) {
      plot.these = !is.na(my.gammas[[j]]) & my.gammas[[j]] < 10
      lines(
        x = my.bts[plot.these],
        y = my.gammas[[j]][plot.these] / se.gammas[plot.these],
        lty = 2,
        lwd = 2,
        col = my.cols[j]
      )
    }
    legend(
      'topleft',
      legend = c(0, my.cvs ^ 2),
      col = c('black', my.cols),
      lty = c(1, rep(2, length(my.cvs))),
      lwd = c(2, rep(2, length(my.cvs))),
      bty = 'n',
      title = expression(CV[alpha] ^ 2),
      ##'Effect dist\'n coef. of var.',
      cex = 0.75
    )
    par(op)
    dev.off()
  }
}

###############################
### Figure 5 genetic deaths ###
###############################
{
  rm(list = ls())
  library(wesanderson)
  Ne <- 20000
  dels <- seq(1 / (2 * Ne) , 1 - 1 / (2 * Ne), by = 0.1 / (2 * Ne))
  prevs <- c (0.001, 0.003, 0.01 , 0.03)
  deaths <- sapply(prevs, function(P)
    dels / (dels + P))
  png(
    'figures/paperFiguresForRealThisTime/GeneticDeathsFigure5.png',
    height = 10 * 1.6 ,
    width = 11 * 1.6,
    units = 'cm',
    res = 500
  )
  op <- par(mar = c(3.4, 3, 2, 1) + 0.1)
  my.cols <- wes_palette('GrandBudapest1', ncol(deaths))
  plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, 1) ,
    bty = 'n',
    xlab = "" ,
    ylab = ""
  )
  abline(a = 1,
         b = 0 ,
         lty = 2 ,
         lwd = 1)
  matplot(
    x = dels,
    y = deaths,
    type = 'l',
    lty = 1 ,
    lwd = 3 ,
    col = my.cols,
    add = T
  )
  legend(
    'bottomright' ,
    legend = prevs ,
    col = my.cols,
    lty = 1 ,
    lwd = 3 ,
    title = expression(paste("Prevalence, ", bar(R) , sep = '')),
    bty = 'n' ,
    cex = 1.4
  )
  mtext(
    side = 1 ,
    text = expression(paste('Risk effect, ', delta(alpha) , sep = '')),
    line = 2.5 ,
    cex = 1.3
  )
  mtext(
    side = 2 ,
    text = "Number of genetic deaths per mutation",
    line = 2,
    cex = 1.3
  )
  par(op)
  dev.off()
}


#####################################
### Figure 7 prevalence isoclines ###
#####################################
{
  rm(list = ls())
  dnorminv <- function(y)
    sqrt(-2 * log(sqrt(2 * pi) * y))
  h2_C2 <- c(0.2 * sqrt(0.2), 0.5 * sqrt(0.5), 0.8 * sqrt(0.8)) ^ 2
  N <- 20000
  bt <- seq(0.001, 0.999, by = 0.001)
  ybar <- log((1 + bt) / (1 - bt))
  log10prevs <- c(4, 3.5, 3, 2.5, 2, 1.5)
  prevs <- 10 ^ -log10prevs
  phits <- dnorm(-qnorm(prevs))
  
  Us <-
    lapply(h2_C2, function(H2C2)
      sapply(phits, function(PHIT)
        PHIT ^ 2 * H2C2 * N / (bt * ybar)))
  
  library(RColorBrewer)
  my.cols <- brewer.pal(n = length(prevs), name = "RdBu")
  this.panel <- c('A' ,'B', 'C')
  titles <- c(
    expression(paste(h^2 , ' = 0.2, ', C, ' = 0.2' , sep = '')),
    expression(paste(h^2 , ' = 0.5, ', C, ' = 0.5' , sep = '')),
    expression(paste(h^2 , ' = 0.8, ', C, ' = 0.8' , sep = ''))
  )
  png(
    'figures/paperFiguresForRealThisTime/PrevalenceIsoclinesFigure7.png',
    height = 5 * 1.6 ,
    width = 11 * 1.6,
    units = 'cm',
    res = 500
  )
  layout(matrix(c(1,2,3),nrow=1))
  for (i in 1:length(Us)) {
    op <- par(mar=c(3.3,3.3,3.3,0.1))
    matplot(
      bt ,
      Us[[i]] ,
      type = 'l' ,
      lty = 1 ,
      lwd = 2 ,
      xlim = c (0, 1),
      ylim = c (0 ,1),
      col = my.cols,
      bty = 'n',
      xlab = '',
      ylab = ''
    )
    mtext(side = 1, line = 2.1, text = expression(paste('Relative threshold position, ' , b[T], sep ='')), cex = 0.8)
    mtext(side = 2, line = 2.1, text = expression(paste('Number of mutations per gamete, ' , U, sep ='')), cex = 0.8)
    mtext(side = 3, line = 0.4, at = 0 ,text = this.panel[i])
    mtext(side = 3, line = 0.4 , text = titles[i])
    if (i == 1) {
      legend(
        'topright',
        legend = c(
          expression(10^-4), 
          expression(10^-3.5),
          expression(10^-3),
          expression(10^-2.5),
          expression(10^-2),
          expression(10^-2.5)
        ),
        col = my.cols ,
        lwd = 2 ,
        lty = 1 ,
        bty = 'n',
        cex = 0.8 ,
        title = "Prevalence"
      )
    }
  }
  dev.off()
}




######################################
### Two Effect Prevalence Figure 8 ###
######################################


source('scripts/twoEffectSameAsEnv.FixedDeltalR.R')
source('scripts/plotTwoDist.R')


##################################################
### Large Effect Size Sensitivity Supp Figures ###
##################################################

source('scripts/largeEffectPopSizeSensitivity.R')

#################################
### Small Effect Bound Figure ###
#################################
{
  rm(list = ls())
  prevs <- 10 ^ seq(-5, -1, length.out = 1000)
  tstar <- qnorm(1 - prevs)
  bound <- 2 * tstar / (3 * tstar ^ 2 - 2)
  
  png(
    'figures/paperFiguresForRealThisTime/SmallEffectBoundSFigure.png',
    height = 5 * 1.6 ,
    width = 11 * 1.6,
    units = 'cm',
    res = 500
  )
  op <- par(mar = c(3.4, 3.4, 1, 0.4))
  plot(
    prevs,
    bound,
    type = 'l',
    lty = 1 ,
    log = 'x',
    ylim = c(0, max(bound)) ,
    xlab = '',
    ylab = '',
    bty = 'n'
  )
  mtext(side = 1 , line = 2.4 , text = 'Prevalence')
  mtext(side = 2 , line = 2.4 , text = 'Bound')
  par(op)
  dev.off()
}





#################################
### Small Effect Bound Figure ###
#################################
{
  rm(list = ls())
  prevs <- 10 ^ seq(log(5e-5), log(0.05,10), length.out = 1000)
  tstar <- qnorm(1 - prevs)
  bound <- (1 / dnorm(tstar)) * (tstar^3 - 3*tstar) / (tstar^2 - 2)
  easierBound <- exp(tstar^2/2)*sqrt(pi/2)*tstar
  
  
  png(
    'figures/paperFiguresForRealThisTime/AllAdditiveBoundSFigure.png',
    height = 5 * 1.6 ,
    width = 11 * 1.6,
    units = 'cm',
    res = 500
  )
  op <- par(mar = c(3.4, 3.4, 1, 0.4))
  plot(
    prevs,
    bound,
    type = 'l',
    lty = 1 ,
    log = 'xy',
    xlim=c(5e-5,5e-2),
    ylim = c(10, 10000) ,
    xlab = '',
    ylab = '',
    bty = 'n',
    xaxt = 'n',
    yaxt = 'n'
  )
  lines(
    prevs,
    easierBound,
    lty=2
  )
  mtext(side = 1 , line = 2.4 , text = 'Prevalence')
  mtext(side = 2 , line = 2.4 , text = 'Lower Bound')
  axis(side=2, at = 10^seq(1,4,by=1))
  axis(side=1, at = 10^seq(log(5e-5,10),log(5e-2,10),by=1))
  par(op)
  dev.off()
  
  png(
    'figures/paperFiguresForRealThisTime/AllAdditiveBoundSFigure2.png',
    height = 5 * 1.6 ,
    width = 11 * 1.6,
    units = 'cm',
    res = 500
  )
  op <- par(mar = c(3.4, 3.4, 1, 0.4))
  plot(
    prevs,
    bound,
    type = 'l',
    lty = 1 ,
    log = 'x',
    xlim=c(5e-3,1e-1),
    ylim = c(-10, 100) ,
    xlab = '',
    ylab = '',
    bty = 'n'
  )
  lines(
    prevs,
    easierBound,
    lty=2
  )
  mtext(side = 1 , line = 2.4 , text = 'Prevalence')
  mtext(side = 2 , line = 2.4 , text = 'Lower Bound')
  #axis(side=2, at = 10^seq(1,4,by=1))
  #axis(side=1, at = 10^seq(log(5e-5,10),log(5e-2,10),by=1))
  abline(h=0,lty = 3)
  par(op)
  dev.off()
  
  
}


#############################
### Berry-Esseen figure 1 ###
#############################

{
  rm(list = ls())
  b <- function(y)
    ifelse(y > 100, 1, (exp(2 * y) - 1) / (exp(2 * y) + 1))
  sig <- function(y)
    sqrt(2 * b(y) * y)
  rhosig <- function(y)
    (2*b(y) + y*(sig(y)^2-2)) / (sig(y)^3)
  
  Ne <- 20000
  Ne2 <- 200000
  L <- c(1e4, 1e5, 1e6, 1e7, 1e8)
  u <- 1e-8
  theta <- 2 * Ne * L * u
  theta2 <- 2 * Ne2 * L * u
  my.gamma <- 10 ^ seq(-2, 3, length.out = 1000)
  my.rhosig <- rhosig(my.gamma)
  my.bound <-
    lapply(theta, function(TH)
      0.56 / sqrt(TH) * my.rhosig)
  my.bound2 <-
    lapply(theta2, function(TH)
      0.56 / sqrt(TH) * my.rhosig)
  my.cols <- wes_palette('FantasticFox1', length(my.bound))
  
  
  png(
    'figures/paperFiguresForRealThisTime/BerryEsseenMaxAbsError.png',
    height = 7 * 1.6 ,
    width = 11 * 1.6,
    units = 'cm',
    res = 500
  )
  op <- par(mar = c(3.4, 3.4, 1, 0.4))
  plot(
    NA,
    type = 'l',
    bty = 'n' ,
    log = 'xy' ,
    ylim = c(1e-4, 1e1) ,
    xlim = c(1e-2,1e3),
    ylab = '' ,
    xlab = ''
  )
  mtext(text = expression(paste( 'Maximum absolute error, ', K / sqrt(theta) %.% rho / sigma ^ 3)),
        side = 2,
        line = 2)
  mtext(text = expression(paste(
    'Population scaled selection coefficient, ', gamma
  )),
  side = 1,
  line = 2.5)
  for (i in 1:length(my.bound)) {
    lines(my.gamma,
          my.bound[[i]],
          col = my.cols[i],
          lty = 1,
          lwd = 2)
  }
  for (i in 1:length(my.bound)) {
    lines(my.gamma,
          my.bound2[[i]],
          col = my.cols[i],
          lty = 2,
          lwd = 2)
  }
  for (i in seq(-4, 1, 1)) {
    abline(h = 10 ^ i , lty = 2 , lwd = 0.2)
  }
  legend(
    'topleft',
    legend = L,
    col = my.cols ,
    lty = 1,
    lwd = 2 ,
    bty = 'n',
    title = 'L'
  )
  legend(
    x=1.5e-1,
    y=15,
    legend = c('20,000','200,000'),
    col = 'black' ,
    lty = c(1,2),
    lwd = 2 ,
    bty = 'n',
    title = 'N'
  )
  par(op)
  dev.off()
  
  
  
  #############################
  ### Berry-Esseen figure 2 ###
  #############################
  
  dnorminv <- function(y)
    sqrt(-2 * log(sqrt(2 * pi) * y))
  cost <- 1 / 2
  ## h2 <- 1 / 2
  U <- L * u
  phit <-
    lapply(U, function(U)
      cost^(-1) * sqrt (2 * U * b(my.gamma) * my.gamma / (2 * Ne )))
  phit2 <-
    lapply(U, function(U)
      cost^(-1) * sqrt (2 * U * b(my.gamma) * my.gamma / (2 * Ne2)))
  
  tstar <-
    lapply(phit, function(PHIT)
      ifelse(PHIT > dnorm(qnorm(0.7)), NA, dnorminv(PHIT)))
  tstar2 <-
    lapply(phit2, function(PHIT)
      ifelse(PHIT > dnorm(qnorm(0.7)), NA, dnorminv(PHIT)))
  prevs <- lapply(tstar, function(TSTAR)
    1 - pnorm(TSTAR))
  prevs2 <- lapply(tstar2, function(TSTAR)
    1 - pnorm(TSTAR))
  
  max.rel.error <-
    mapply(function(X, Y)
      X / Y , X = my.bound, Y = prevs)
  max.rel.error2 <-
    mapply(function(X, Y)
      X / Y , X = my.bound2, Y = prevs2)
  png(
    'figures/paperFiguresForRealThisTime/BerryEsseenMaxRelError.png',
    height = 7 * 1.6 ,
    width = 11 * 1.6,
    units = 'cm',
    res = 500
  )
  op <- par(mar = c(3.4, 3.4, 1, 0.4))
  plot(
    NA,
    type = 'l',
    bty = 'n' ,
    log = 'xy' ,
    ylim = c(1e-1, 2e6) ,
    xlim = c(1e-2,1e3),
    ylab = '' ,
    xlab = ''
  )
  # mtext(text = expression(paste(C, h ^ 3 / sqrt(theta) %.% rho / sigma ^ 3)),
  #       side = 2,
  #       line = 2)
  mtext(text = 'Maximum relative error',
        side = 2,
        line = 2)
  mtext(text = expression(paste(
    'Population scaled selection coefficient, ', gamma
  )),
  side = 1,
  line = 2.5)
  for (i in 1:ncol(max.rel.error)) {
    lines(my.gamma,
          max.rel.error[,i],
          col = my.cols[i],
          lty = 1,
          lwd = 2)
  }
  for (i in 1:ncol(max.rel.error)) {
    lines(my.gamma,
          max.rel.error2[,i],
          col = my.cols[i],
          lty = 2,
          lwd = 2)
  }
  for (i in seq(-1, 6, 1)) {
    abline(h = 10 ^ i , lty = 2 , lwd = 0.2)
  }
  legend(
    x=1e1,
    y=1.5e6,
    legend = c('20,000','200,000'),
    col = 'black' ,
    lty = c(1,2),
    lwd = 2 ,
    bty = 'n',
    title = 'N'
  )
  legend(
    'topright',
    legend = L,
    col = my.cols ,
    lty = 1,
    lwd = 2 ,
    bty = 'n',
    title = 'L'
  )
  par(op)
  dev.off()

  
  png(
    'figures/paperFiguresForRealThisTime/BerryEsseenPrevVsError.png',
    height = 7 * 1.6 ,
    width = 11 * 1.6,
    units = 'cm',
    res = 500
  )
  op <- par(mar = c(3.4, 3.4, 1, 0.4))
  plot(
    NA,
    type = 'l',
    bty = 'n' ,
    log = 'xy' ,
    ylim = c(1e-1, 1e6) ,
    xlim = c(1e-7,3e-1),
    ylab = '' ,
    xlab = ''
  )
  # mtext(text = expression(paste(C, h ^ 3 / sqrt(theta) %.% rho / sigma ^ 3)),
  #       side = 2,
  #       line = 2)
  mtext(text = 'Maximum relative error',
        side = 2,
        line = 2)
  mtext(text = 'Prevalence',
  side = 1,
  line = 2.5)
  for (i in 1:ncol(max.rel.error)) {
    lines(prevs[[i]],
          max.rel.error[,i],
          col = my.cols[i],
          lty = 1,
          lwd = 2)
  }
  for (i in 1:ncol(max.rel.error)) {
    lines(prevs2[[i]],
          max.rel.error2[,i],
          col = my.cols[i],
          lty = 2,
          lwd = 2)
  }
  for (i in seq(-1, 6, 1)) {
    abline(h = 10 ^ i , lty = 2 , lwd = 0.2)
  }
  legend(
    'topright',
    legend = L,
    col = my.cols ,
    lty = 1,
    lwd = 2 ,
    bty = 'n',
    title = 'L'
  )
  legend(
    x=1e-3,
    y=1.5e6,
    legend = c('20,000','200,000'),
    col = 'black' ,
    lty = c(1,2),
    lwd = 2 ,
    bty = 'n',
    title = 'N'
  )
  par(op)
  dev.off()
  
}

