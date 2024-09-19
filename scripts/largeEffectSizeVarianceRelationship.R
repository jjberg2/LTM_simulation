

#################################
## Four panel variance figures ##
#################################

rm(list = ls())
library("MetBrewer")
library('numDeriv')
source('scripts/simpleVarFuncs.R')
prevs = c(0.03, 0.01, 0.003, 0.001)
my.cols = met.brewer('Hiroshige', length(prevs))
thr = sapply(prevs, function(PREV)
  qnorm(1 - PREV))
phits = dnorm(thr)
max.a = 5
min.a = rep(0.001, 4)
my.a = mapply(
  function(MAXA, MINA)
    seq(MINA, MAXA, length.out = 1000),
  MINA = min.a,
  MAXA = max.a,
  SIMPLIFY = FALSE
)
Ne = 10000
this.cost <- 1
norm.del <- mapply(function(THR, A)
  normDel(A, THR),
  THR = thr,
  A = my.a)
li.vars <- mapply(function(THR, A)
  sigmaa(A, THR, Ne * this.cost),
  THR = thr,
  A = my.a)
risk.vars <- mapply(function(THR, A)
  sigmaRisk(A, THR, Ne * this.cost),
  THR = thr,
  A = my.a)
small.approx.vars = mapply(function(THR, A)
  sigmaaSmall(A, THR, Ne * this.cost),
  THR = thr,
  A = my.a)
plot.a <- do.call (cbind, my.a)
my.cex <- 1.5
pdf(
  'figures/suppFigures/largeEffectSizeVarianceRelationship_2by2.pdf',
  width = 12,
  height = 12
)
par(mfrow = c(2, 2))
matplot(
  x = norm.del ,
  y = risk.vars,
  xlim = c(0, max(unlist(norm.del))),
  ylim = c(0, max(risk.vars)),
  type = 'l',
  lty = 1 ,
  lwd = 3 ,
  col = my.cols,
  xlab = 'Risk scale effect size' ,
  ylab = 'Risk scale variance',
  cex.lab = my.cex,
  cex.axis = my.cex
)
mtext('A', side = 3 , at = 0, cex = 3)
matplot(
  x = plot.a ,
  y = norm.del,
  xlim = c(0, max(unlist(my.a))),
  ylim = c(0, 1),
  type = 'l',
  lty = 1 ,
  lwd = 3 ,
  col = my.cols,
  xlab = 'Standardized liability scale effect size' ,
  ylab = 'Risk scale effect size',
  cex.lab = my.cex,
  cex.axis = my.cex
)
mtext('B', side = 3 , at = 0, cex = 3)
matplot(
  x = norm.del ,
  y = li.vars,
  xlim = c(0, max(unlist(norm.del))),
  ylim = c(0, max(li.vars)),
  type = 'l',
  lty = 1 ,
  lwd = 3 ,
  col = my.cols,
  xlab = 'Risk scale effect size' ,
  ylab = 'Liability scale variance',
  cex.lab = my.cex,
  cex.axis = my.cex
)

legend(
  'topright',
  legend = prevs,
  lty = 1,
  lwd = 2,
  col = my.cols,
  bty = 'n',
  title = "Prevalence",
  cex = my.cex
)
mtext('C', side = 3 , at = 0, cex = 3)
matplot(
  x = plot.a ,
  y = li.vars,
  xlim = c(0, max(unlist(my.a))),
  ylim = c(0, max(li.vars)),
  type = 'l',
  lty = 1 ,
  lwd = 3 ,
  col = my.cols,
  xlab = 'Standardized liability scale effect size' ,
  ylab = 'Liability scale variance',
  cex.lab = my.cex,
  cex.axis = my.cex
)
mtext('D', side = 3 , at = 0, cex = 3)
dev.off()




#######################################################
## Four panel variance figure with a_gamma on x axis ##
#######################################################

## bt = 0.6479407
## bt*ay = 1


rm(list = ls())
library("MetBrewer")
library('numDeriv')
source('scripts/simpleVarFuncs.R')
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
h2 = 1/2
aybt <- 1
thetaU = 4*Ne*U
this.cost <- 1/2
## my.ay <- lapply(my.a,function(A) A*sqrt(thetaU / h2))
my.ay <- mapply(function(A,PHIT) 4*Ne*this.cost*PHIT*A,
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
pdf(
  'figures/suppFigures/largeEffectSizeVarianceRelationship_ay.pdf',
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
  legend = 'Small effect approx',
  lty = 2,
  lwd = 3,
  col = 'black',
  bty = 'n',
  cex = 1
)
mtext(
  text = 'Liability effect size',
  side = 1,
  line = 2,
  cex = my.cex
)
mtext(
  text = '(Population scaled fitness units)',
  side = 1, 
  line = 3,
  cex = my.cex
)
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
  line = 2,
  cex = my.cex
)
mtext(
  text = '(Population scaled fitness units)',
  side = 1, 
  line = 3,
  cex = my.cex
)
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














par(mfrow = c(1, 2))
matplot(
  x = plot.a ,
  y = plot.a^2,
  xlim = c(0, max(plot.a)),
  ylim = c(0, max(plot.a^2)),
  type = 'l',
  lty = 1 ,
  lwd = 3 ,
  col = my.cols,
  xlab = 'Risk scale effect size' ,
  ylab = 'Risk scale variance',
  cex.lab = my.cex,
  cex.axis = my.cex
)
mtext('A', side = 3 , at = 0, cex = 3)
matplot(
  x = plot.a ,
  y = 1/norm.del,
  xlim = c(0, max(unlist(my.a))),
  type = 'l',
  lty = 1 ,
  lwd = 3 ,
  col = my.cols,
  xlab = 'Standardized liability scale effect size' ,
  ylab = 'Risk scale effect size',
  cex.lab = my.cex,
  cex.axis = my.cex,
  log = 'y'
)
mtext('B', side = 3 , at = 0, cex = 3)

###############################
## Ratios of number of sites ##
###############################

this.cost <- 1
Ne <- 3000
baseline.del <- 2.94 / (4*Ne*this.cost)


#####################
## Variance ratios ##
#####################


delseq.ratio <- norm.del^2 / baseline.del^2
baseline.a <- sapply(thr,function(THR) 2.94 / (4*Ne*dnorm(THR)*this.cost)  )
asq.ratio <- t(t(plot.a^2) / baseline.a^2)

## simple.var.ratio <- norm.del^2 / plot.a^2

par(mfrow=c(1,2))
## 1
matplot(
  plot.a,
  asq.ratio/delseq.ratio,
  type = 'l',
  lty = 1 ,
  lwd = 3 ,
  col = my.cols,
  xlab = 'Standardized effect size',
  ylab = 'Relative contribution to liability scale variance'
)
## 2
matplot(
  norm.del,
  asq.ratio/delseq.ratio,
  type = 'l',
  lty = 1 ,
  lwd = 3 ,
  col = my.cols,
  xlab = 'Risk scale effect size',
  ylab = 'Relative contribution to liability scale variance'
)


phit2 <- dnorm(thr)^2

scaled.norm.del.sq <- t(t(norm.del)^2 / phit2)
simple.var.ratio <- scaled.norm.del.sq / plot.a^2

matplot(
  plot.a,
  simple.var.ratio,
  type = 'l',
  lty = 1 ,
  lwd = 3 ,
  col = my.cols,
  xlab = 'Standardized effect size',
  ylab = 'Variance contribution ratio',
  log = 'y'
)




#################################
## how many sites do you need? ##
#################################

## need to deal with the fact that I added an rm(list=ls()) above

rm(list=ls())
library("MetBrewer")
library('numDeriv')
source('scripts/simpleVarFuncs.R')
prevs = c(0.03, 0.01, 0.003, 0.001)
thr = sapply(prevs, function(PREV)
  qnorm(1 - PREV))
my.cols = met.brewer('Hiroshige', length(prevs))
Ne = 3000
my.ratios <- lapply(thr,
                    function(THR) {
                      tmp.bt <- 0.9
                      min.y <- log((1 + tmp.bt) / (1 - tmp.bt))
                      min.del <- min.y / (4 * Ne)
                      my.prev <- 1 - pnorm(THR)
                      min.a <- uniroot(
                        f = function(A)
                          normDel(A, THR) - min.del ,
                        interval = c(0, 5)
                      )$root
                      max.a <- qnorm(0.99) - qnorm(my.prev)
                      these.as <- seq(min.a, max.a, length.out = 1000)
                      these.dels <- normDel(these.as, THR)
                      delsq.ratio <- min.del ^ 2 / these.dels ^ 2
                      asq.ratio <- these.as ^ 2 / min.a ^ 2
                      bgamma.ratio <-
                        bgamma(4 * Ne * these.dels) / bgamma(4 * Ne * min.del)
                      h2l.ratio <- delsq.ratio * asq.ratio
                      phit <- dnorm(THR)
                      new.h2l.ratio <- (these.as^2*phit^2) / these.dels ^ 2 
                      list(these.dels,
                           h2l.ratio,
                           these.as,
                           delsq.ratio,
                           asq.ratio,
                           bgamma.ratio,
                           new.h2l.ratio)
                    })
my.dels <- lapply(my.ratios, function(X)
  X[[1]])
my.h2ratios <- lapply(my.ratios, function(X)
  X[[2]])
my.as <- lapply(my.ratios, function(X)
  X[[3]])
delsq.ratio <- lapply(my.ratios, function(X)
  X[[4]])
new.h2.ratios <- lapply(my.ratios, function(X)
  X[[7]])


pdf(
  'figures/suppFigures/largeEffectSizeRiskVsLiabilityVariance1.pdf',
  width = 12,
  height = 6
)
par(mfrow = c(1, 2))
plot(
  NA,
  xlim = c(0, 1),
  ylim = c(1, max(1/unlist(my.h2ratios))),
  type = 'l',
  xlab = 'Risk effect of large effect alleles',
  ylab = '',
  bty = 'n',
  log = 'y'
)
mtext('Factor by which large effect variance contributions are increased',
      side = 2 ,
      line = 3)
mtext('on the risk scale relative to that of small effect contributions',
      side = 2 ,
      line = 2)
for (j in 1:length(my.dels)) {
  lines(
    my.dels[[j]],
    1/my.h2ratios[[j]],
    lty = 1,
    lwd = 3,
    col = my.cols[j]
  )
}
plot(
  NA,
  xlim = c(0, max(unlist(my.as))),
  ylim = c(1, max(1/unlist(my.h2ratios))),
  type = 'l',
  xlab = 'Standardized liability effect of large effect alleles',
  ylab = '',
  bty = 'n',
  log = 'y'
)
mtext('Factor by which large effect variance contributions are increased',
      side = 2 ,
      line = 3)
mtext('on the risk scale relative to that of small effect contributions',
      side = 2 ,
      line = 2)
for (j in 1:length(my.dels)) {
  lines(my.as[[j]],
        1/my.h2ratios[[j]],
        lty = 1,
        lwd = 3,
        col = my.cols[j])
}
legend(
  'topleft',
  bty = 'n',
  lty = 1 ,
  lwd = 2,
  col = my.cols,
  legend = prevs,
  title = 'Prevalence'
)
dev.off()



pdf(
  'figures/suppFigures/largeEffectSizeRiskVsLiabilityVariance2.pdf',
  width = 12,
  height = 6
)
par(mfrow = c(1, 2))
plot(
  NA,
  xlim = c(0, 1),
  ylim = c(1e-4, 1),
  type = 'l',
  xlab = 'Risk effect of large effect alleles',
  ylab = '',
  bty = 'n',
  log = 'y'
)
mtext('Relative contribution to liability scale variance',
      side = 2 ,
      line = 3)
mtext('if risk scale contributions are equal',
      side = 2 ,
      line = 2)
for (j in 1:length(my.dels)) {
  lines(
    my.dels[[j]],
    my.h2ratios[[j]],
    lty = 1,
    lwd = 3,
    col = my.cols[j]
  )
}
legend(
  'topright',
  bty = 'n',
  lty = 1 ,
  lwd = 2,
  col = my.cols,
  legend = prevs,
  title = 'Prevalence'
)
plot(
  NA,
  xlim = c(0, max(unlist(my.as))),
  ylim = c(1e-4, 1),
  type = 'l',
  xlab = 'Standardized liability effect of large effect alleles',
  ylab = '',
  bty = 'n',
  log = 'y'
)
mtext('Relative contribution to liability scale variance',
      side = 2 ,
      line = 3)
mtext('if risk scale contributions are equal',
      side = 2 ,
      line = 2)
for (j in 1:length(my.dels)) {
  lines(my.as[[j]],
        my.h2ratios[[j]],
        lty = 1,
        lwd = 3,
        col = my.cols[j])
}
dev.off()


pdf(
  'figures/suppFigures/largeEffectSizeRiskVsLiabilityVariance3.pdf',
  width = 6,
  height = 6
)
par(mfrow = c(1, 1))
plot(
  NA,
  xlim = c(0, 1),
  ylim = c(1e-4, 1),
  type = 'l',
  xlab = 'Risk effect of large effect alleles',
  ylab = '',
  bty = 'n',
  log = 'y'
)
mtext('Relative contribution to liability scale variance',
      side = 2 ,
      line = 2.5)
for (j in 1:length(my.dels)) {
  lines(
    my.dels[[j]],
    new.h2.ratios[[j]],
    lty = 1,
    lwd = 3,
    col = my.cols[j]
  )
}
legend(
  'topright',
  bty = 'n',
  lty = 1 ,
  lwd = 2,
  col = my.cols,
  legend = prevs,
  title = 'Prevalence'
)
dev.off()





## how does Poisson shape change things?
sigmaa2 = function(a, del, nc) {
  a ^ 2 * bgamma(4 * nc * del) / (2 * nc * del)
}
source('scripts/solveTwoEffect.R')

soln.in <-
  read.table("twoEffectPrevInsensitivityParamsTable.txt", header = T)
source('scripts/makeTwoEffectPrevInit.R')
keep <- c(2, which.min(abs(soln.in$deltal - 1 / 2)))
my.solns <- soln.in[keep, ]


Ve.mult = c(3, 3, 0.9, 0.65)
tuned.costs = c(1, 0.131, 1, 0.131)
tmp.solns <- list()
for (j in 1:2) {
  tmp.solns[[j]] <- solveTwoEffect(
    bs = my.solns[1, 'bs'],
    bt = my.solns[1, 'bt'],
    Ne = my.solns[1, 'Ne'],
    Ls = my.solns[1, 'Ls'],
    Ll = my.solns[1, 'Ll'],
    as = my.solns[1, 'as'],
    al = my.solns[1, 'al'],
    L.init = L.init,
    r2n = NULL,
    Lmeana = Lmeana,
    Ve = my.solns[1, 'Ve'] * Ve.mult[j],
    u = u,
    C = tuned.costs[j],
    LL.soln = TRUE,
    var.ratio = 4,
    equalize.observed.vars = TRUE
  )
}
for (j in 3:4) {
  tmp.solns[[j]] <- solveTwoEffect(
    bs = my.solns[2, 'bs'],
    bt = my.solns[2, 'bt'],
    Ne = my.solns[2, 'Ne'],
    Ls = my.solns[2, 'Ls'],
    Ll = my.solns[2, 'Ll'],
    as = my.solns[2, 'as'],
    al = my.solns[2, 'al'],
    L.init = L.init,
    r2n = NULL,
    Lmeana = Lmeana,
    Ve = my.solns[2, 'Ve'] * Ve.mult[j],
    u = u,
    C = tuned.costs[j],
    LL.soln = TRUE,
    var.ratio = 4,
    equalize.observed.vars = TRUE
  )
}
new.solns <- as.data.frame(do.call(rbind,
                                   lapply(tmp.solns,
                                          makeOutput)))
new.solns$prev

tol <- 1e-4

deltals <- list()
my.dists <- list()
for (i in 1:4) {
  ## solve for min an max genetic liabilities
  ## that we need to plot distribution
  min.gl <- uniroot(function(X)
    tol - (1 - pPoisConv(X, new.solns[i, 'mean.nl'], new.solns[i, 'norm.sd'], alphal =
                           new.solns[i, 'al'])),
    interval = c(-10 * new.solns[i, 'tstar'], 10 * new.solns[i, 'tstar']))$root
  max.gl <- uniroot(
    function(X)
      tol - pPoisConv(X, new.solns[i, 'mean.nl'], new.solns[i, 'norm.sd'], alphal =
                        new.solns[i, 'al']),
    interval = c(-10 * new.solns[i, 'tstar'], 10 * new.solns[i, 'tstar'])
  )$root
  seq.li <- seq(min.gl, max.gl, length.out = 1000)
  li.dense <-
    sapply(seq.li, function(G)
      dPoisConv(G, new.solns[i, 'mean.nl'], new.solns[i, 'norm.sd'], alphal = new.solns[i, 'al']))
  
  
  Vt <- new.solns[i, 'Vg'] + new.solns[i, 'Ve']
  std.max.gl <- max.gl / sqrt (Vt)
  alpha.seq <-
    seq(1e-8, new.solns[i, 'tstar'] - min.gl, length.out = 1000)
  pens <- sapply(alpha.seq,
                 function(A)
                   pPoisConv(
                     t = new.solns[i, 'tstar'] - A,
                     lambda = new.solns[i, 'mean.nl'],
                     norm.sd = new.solns[i, 'norm.sd'],
                     alphal = new.solns[i, 'al']
                   ))
  deltals[[i]] <- pens - new.solns[i, 'prev']
  tmp.vars <-
    sigmaa2(alpha.seq / sqrt(Vt), deltals[[i]], Ne * tuned.costs[i])
  my.dists[[i]] <-
    list(
      seq.li = seq.li,
      li.dense = li.dense,
      alpha.seq = alpha.seq,
      std.alpha.seq = alpha.seq / sqrt(Vt),
      deltals = deltals[[i]],
      vars = tmp.vars
    )
}

##prevs = new.solns$prev
prevs = c(0.002, 0.02, 0.002, 0.02)
thr = sapply(prevs, function(PREV)
  qnorm(1 - PREV))
max.a = thr + max(thr)
min.a = rep(0.001, 4)
N = 3000
my.a = mapply(
  function(MAXA, MINA)
    seq(MINA, MAXA, length.out = 1000),
  MINA = min.a,
  MAXA = max.a,
  SIMPLIFY = FALSE
)
alt.del = mapply(
  function(THR, A)
    normDel(A, THR),
  THR = thr,
  A = my.a,
  SIMPLIFY = FALSE
)
alt.vars = mapply(
  function(DEL, A, COST)
    sigmaa2(A, DEL, Ne * COST),
  DEL = alt.del,
  A = my.a,
  COST = tuned.costs
)





pdf(
  'figures/suppFigures/largeEffectSizeVarianceRelationshipPoisson.pdf',
  width = 12,
  height = 12
)
par(mfrow = c(1, 2))
op = par(mar = c(5, 5, 4, 4) + 0.1)


## plot risk effects
plot(
  NA,
  type = 'l',
  xlim = c(0, 6),
  ylim = c(0, 1),
  col = my.cols[1],
  lwd = 2
)
lines(my.dists[[1]]$std.alpha.seq ,
      deltals[[1]],
      col = my.cols[1],
      lwd = 2)
lines(my.dists[[2]]$std.alpha.seq,
      deltals[[2]],
      col = my.cols[2],
      lwd = 2)
lines(my.dists[[3]]$std.alpha.seq,
      deltals[[3]],
      col = my.cols[3],
      lwd = 2)
lines(my.dists[[4]]$std.alpha.seq,
      deltals[[4]],
      col = my.cols[4],
      lwd = 2)


for (i in 1:ncol(alt.vars)) {
  lines(
    x = my.a[[i]],
    y = alt.del[[i]],
    col = my.cols[i],
    lwd = 2,
    lty = 2
  )
}


if(FALSE){
  
  
  
}


if (FALSE) {
  par(op)
  ## plot variance contributions
  plot(
    NA,
    type = 'l',
    xlim = c(0, 6),
    ylim = c(0, max(c(
      my.dists[[3]]$vars, alt.vars
    ))),
    col = my.cols[1],
    lwd = 2
  )
  lines(my.dists[[1]]$std.alpha.seq ,
        my.dists[[1]]$vars,
        col = my.cols[1],
        lwd = 2)
  lines(my.dists[[2]]$std.alpha.seq,
        my.dists[[2]]$vars,
        col = my.cols[2],
        lwd = 2)
  lines(my.dists[[3]]$std.alpha.seq,
        my.dists[[3]]$vars,
        col = my.cols[3],
        lwd = 2)
  lines(my.dists[[4]]$std.alpha.seq,
        my.dists[[4]]$vars,
        col = my.cols[4],
        lwd = 2)
  for (i in 1:ncol(alt.vars)) {
    lines(
      x = my.a[[i]],
      y = alt.vars[, i],
      col = my.cols[i],
      lwd = 2,
      lty = 2
    )
  }
  
  dev.off()
  
  
  
  
  ## plot liability distributions
  plot(my.dists[[1]]$seq.li,
       my.dists[[1]]$li.dense,
       type = 'l')
  lines(my.dists[[2]]$seq.li,
        my.dists[[2]]$li.dense,
        col = 'red')
  lines(my.dists[[3]]$seq.li,
        my.dists[[3]]$li.dense,
        lty = 2,
        col = 'black')
  lines(my.dists[[3]]$seq.li,
        my.dists[[3]]$li.dense,
        lty = 2,
        col = 'red')
  
  
}
