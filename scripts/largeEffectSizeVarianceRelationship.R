library("MetBrewer")
library('numDeriv')
prevs = c(0.03,0.01,0.003,0.001)
thr = sapply(prevs, function(PREV) qnorm(1-PREV))
max.a = thr + max(thr)
min.a = rep(0.001,4)
my.a = mapply(function(MAXA,MINA) seq(MINA,MAXA,length.out=1000), MINA = min.a, MAXA = max.a,SIMPLIFY = FALSE)
N = 3000
cost = 1


bgamma = function(gamma) {
    ifelse(gamma>20,1,(exp(gamma)-1)/(exp(gamma)+1))
}

normDel= function(a,thr){
    pnorm(thr) - pnorm(thr - a)
}
sigmaa = function(a,thr,nc){
    del = normDel(a,thr)
    a^2 * bgamma(4*nc*del) / (2*nc*del)
}
sigmaa2 = function(a,del,nc){
  a^2 * bgamma(4*nc*del) / (2*nc*del)
}

my.del = mapply( function(THR,A) normDel(A,THR), THR = thr, A = my.a)
my.vars = mapply(function(THR,A) sigmaa(A,THR,N*cost), THR = thr, A = my.a)
# my.derivs = list()
# for( j in 1:length(thr)){
#     my.derivs[[j]] = matrix(nrow=length(my.a[[j]]),ncol=2)
#     for( i in 1:length(my.a[[j]])){
#         my.derivs[[j]][i,] = genD(
#             func = function(X) sigmaa(X,thr[j],N*cost),
#             x = my.a[[j]][i]
#         )$D
#     }
# }


cex.lab = 1.4
pdf('figures/suppFigures/largeEffectSizeVarianceRelationship1.pdf',width=12,height=12)
par(mfrow=c(2,2))
op = par(mar=c(5,5,4,4)+0.1)

##1 
plot(
    NA,
    xlim = c(0,max(max.a)),
    ylim = c(0,1),
    type = 'l',
    xlab = '',
    ylab = ''
)
for(i in 1:length(my.a)){
    lines(
        my.a[[i]],
        my.del[,i],
        lty = i
    )
}
mtext(
    text=expression(paste('Standardized effect size (', alpha[std] ,')',sep='')),
    side=1,
    line=3,
    cex=cex.lab
)
mtext(
    text=expression(paste('Risk scale effect size (', delta[R](alpha) ,')',sep='')),
    side=2,
    line=3,
    cex=cex.lab
)
text(
    x = 0.8,
    y = 0.9,
    labels = 'Prevalence'
)
legend(
    x = -0.1,
    y = 0.9,
    legend = prevs,
    lty = 1:4,
    bty = 'n'
)


my.costs <- c(1, 2 / 3, 1 / 3)
for (j in 1:length(my.costs)) {
  my.vars = mapply(function(THR, A)
    sigmaa(A, THR, N * my.costs[j]),
    THR = thr,
    A = my.a)
  ##2
  if(j%in%c(1,3)) {
     par(op)
  } else {
    op = par(mar=c(5,5,4,4)+0.1)
  }
  plot(
    NA,
    xlim = c(0, max(max.a)),
    ylim = c(0, 0.035),
    type = 'l',
    xlab = '',
    ylab = ''
  )
  for (i in 1:length(my.a)) {
    lines(my.a[[i]],
          my.vars[, i],
          lty = i)
  }
  mtext(
    text = expression(paste(
      'Standardized effect size (', alpha[std] , ')', sep = ''
    )),
    side = 1,
    line = 3,
    cex = cex.lab
  )
  mtext(
    text = 'Liability scale contribution to variance',
    side = 2,
    line = 3.75,
    cex = cex.lab
  )
  mtext(
    text = expression(alpha ^ 2 * b(alpha) / (2 * N * delta[R](alpha) * C)),
    side = 2,
    line = 2,
    cex = cex.lab
  )
  if(j!=3){
    blah = substitute(paste("Cost = ", a / 3),list(a=j))
  } else {
    blah = "Cost = 1"
  }
  mtext(
    text = blah,
    side = 3,
    line = 1,
    cex = cex.lab
  )
}
dev.off()


## second normal figure
four.costs <- seq(0.1, 1, length.out = 4)
init.prev = 0.001
thr = rep(qnorm(1 - init.prev),4)
##thr <- dnorminv(dnorm(init.thr) / four.costs)
##prevs <- 1 - pnorm(thr)
max.a = thr + max(thr)
min.a = rep(0.001, 4)
my.a = mapply(
  function(MAXA, MINA)
    seq(MINA, MAXA, length.out = 1000),
  MINA = min.a,
  MAXA = max.a,
  SIMPLIFY = FALSE
)
my.vars.four.costs001 = mapply(
  function(THR, A, COST)
    sigmaa(A, THR, N * COST),
  THR = thr,
  A = my.a,
  COST = four.costs,
  SIMPLIFY = FALSE
)


cex.lab = 1.4
pdf('figures/suppFigures/largeEffectSizeVarianceRelationship2.pdf',width=12,height=6)
par(mfrow=c(1,2))
op = par(mar=c(5,5,4,4)+0.1)
plot(
  NA,
  xlim = c(0,max(max.a)),
  ylim = c(0,0.12),
  type = 'l',
  xlab ='',
  ylab =''
)
j=1
for(i in 1:length(my.a)){
  lines(
    my.a[[i]],
    my.vars.four.costs001[[i]],
    lty = i
  )
}
mtext(
  text=expression(paste('Standardized effect size (', alpha[std] ,')',sep='')),
  side=1,
  line=3,
  cex=cex.lab
)
mtext(
  text = 'Liability scale contribution to variance',
  side = 2,
  line = 3.75,
  cex = cex.lab
)
mtext(
  text = expression(alpha ^ 2 * b(alpha) / (2 * N * delta[R](alpha) * C)),
  side = 2,
  line = 2,
  cex = cex.lab
)
mtext(
  text = "Prevalence = 0.001",
  side = 3,
  line = 1,
  cex = cex.lab
)


par(op)
init.prev = 0.01
thr = rep(qnorm(1 - init.prev),4)
my.vars.four.costs01 = mapply(
  function(THR, A, COST)
    sigmaa(A, THR, N * COST),
  THR = thr,
  A = my.a,
  COST = four.costs,
  SIMPLIFY = FALSE
)
plot(
  NA,
  xlim = c(0,max(max.a)),
  ylim = c(0,0.12),
  type = 'l',
  xlab ='',
  ylab =''
)
j=1
for(i in 1:length(my.a)){
  lines(
    my.a[[i]],
    my.vars.four.costs01[[i]],
    lty = i
  )
}
mtext(
  text=expression(paste('Standardized effect size (', alpha[std] ,')',sep='')),
  side=1,
  line=3,
  cex=cex.lab
)
mtext(
  text = 'Liability scale contribution to variance',
  side = 2,
  line = 3.75,
  cex = cex.lab
)
mtext(
  text = expression(alpha ^ 2 * b(alpha) / (2 * N * delta[R](alpha) * C)),
  side = 2,
  line = 2,
  cex = cex.lab
)
text(
  x = 0.8,
  y = 0.12,
  labels = 'Fitness cost'
)
legend(
  x = -0.1,
  y = 0.12,
  legend = four.costs,
  lty = 1:4,
  bty = 'n'
)
mtext(
  text = "Prevalence = 0.01",
  side = 3,
  line = 1,
  cex = cex.lab
)
dev.off()







sigmaa2 = function(a,del,nc){
  a^2 * bgamma(4*nc*del) / (2*nc*del)
}
source('scripts/solveTwoEffect.R')

soln.in <- read.table("twoEffectPrevInsensitivityParamsTable.txt",header = T)
source('scripts/makeTwoEffectPrevInit.R')
keep <- c(2,which.min(abs(soln.in$deltal - 1/2)))
my.solns <- soln.in[keep,]


Ve.mult = c(3,3,0.9,0.65)
tuned.costs = c(1,0.131,1,0.131)
tmp.solns <- list()
for(j in 1:2) {
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
for(j in 3:4) {
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
  min.gl <- uniroot(
    function(X)
      tol - (
        1 - pPoisConv(X, new.solns[i,'mean.nl'], new.solns[i, 'norm.sd'], alphal =
                        new.solns[i, 'al'])
      ),
    interval = c(-10 * new.solns[i,'tstar'], 10 * new.solns[i,'tstar'])
  )$root
  max.gl <- uniroot(
    function(X)
      tol - pPoisConv(X, new.solns[i,'mean.nl'], new.solns[i,'norm.sd'], alphal =
                        new.solns[i,'al']),
    interval = c(-10 * new.solns[i,'tstar'], 10 * new.solns[i,'tstar'])
  )$root
  seq.li <- seq(min.gl, max.gl, length.out = 1000)
  li.dense <-
    sapply(seq.li, function(G)
      dPoisConv(G, new.solns[i,'mean.nl'], new.solns[i,'norm.sd'], alphal = new.solns[i,'al']))
  
  
  Vt <- new.solns[i, 'Vg'] + new.solns[i, 'Ve']
  std.max.gl <- max.gl / sqrt (Vt)
  alpha.seq <- seq(1e-8, new.solns[i, 'tstar'] - min.gl, length.out = 1000)
  pens <- sapply(alpha.seq,
                 function(A)
                   pPoisConv(
                     t = new.solns[i, 'tstar'] - A,
                     lambda = new.solns[i, 'mean.nl'],
                     norm.sd = new.solns[i, 'norm.sd'],
                     alphal = new.solns[i, 'al']
                   ))
  deltals[[i]] <- pens - new.solns[i,'prev']
  tmp.vars <- sigmaa2(alpha.seq/sqrt(Vt),deltals[[i]],N*tuned.costs[i])
  my.dists[[i]] <- list(seq.li = seq.li, li.dense = li.dense, alpha.seq = alpha.seq, std.alpha.seq = alpha.seq/sqrt(Vt), deltals = deltals[[i]], vars = tmp.vars)
}    

##prevs = new.solns$prev
prevs = c(0.002,0.02,0.002,0.02)
my.cols = met.brewer('Hiroshige',4)[c(1,3,2,4)]
thr = sapply(prevs, function(PREV) qnorm(1-PREV))
max.a = thr + max(thr)
min.a = rep(0.001,4)
N = 3000
my.a = mapply(function(MAXA,MINA) seq(MINA,MAXA,length.out=1000), MINA = min.a, MAXA = max.a,SIMPLIFY = FALSE)
alt.del = mapply( function(THR,A) normDel(A,THR), THR = thr, A = my.a,SIMPLIFY=FALSE)
alt.vars = mapply(function(DEL,A,COST) sigmaa2(A,DEL,N*COST), DEL = alt.del, A = my.a, COST = tuned.costs)



## plot variance contributions
plot(
  NA,
  type='l',
  xlim=c(0,6),
  ylim=c(0,max(c(my.dists[[3]]$vars,alt.vars))),
  col = my.cols[1],
  lwd=2
)
lines(
  my.dists[[1]]$std.alpha.seq , 
  my.dists[[1]]$vars,
  col = my.cols[1],
  lwd=2
)
lines(
  my.dists[[2]]$std.alpha.seq,
  my.dists[[2]]$vars,
  col = my.cols[2],
  lwd=2
)
lines(
  my.dists[[3]]$std.alpha.seq,
  my.dists[[3]]$vars,
  col = my.cols[3],
  lwd=2
)
lines(
  my.dists[[4]]$std.alpha.seq,
  my.dists[[4]]$vars,
  col = my.cols[4],
  lwd=2
)
for( i in 1:ncol(alt.vars)){
  lines(
    x = my.a[[i]],
    y = alt.vars[,i],
    col = my.cols[i],
    lwd = 2,
    lty = 2
  )
}



## plot risk effects
plot(
  NA,
  type='l',
  xlim=c(0,6),
  ylim=c(0,1),
  col = my.cols[1],
  lwd=2
)
lines(
  my.dists[[1]]$std.alpha.seq , 
  deltals[[1]],
  col = my.cols[1],
  lwd=2
)
lines(
  my.dists[[2]]$std.alpha.seq,
  deltals[[2]],
  col = my.cols[2],
  lwd=2
)
lines(
  my.dists[[3]]$std.alpha.seq,
  deltals[[3]],
  col = my.cols[3],
  lwd=2
)
lines(
  my.dists[[4]]$std.alpha.seq,
  deltals[[4]],
  col = my.cols[4],
  lwd=2
)


for( i in 1:ncol(alt.vars)){
  lines(
    x = my.a[[i]],
    y = alt.del[[i]],
    col = my.cols[i],
    lwd = 2,
    lty = 2
  )
}





## plot liability distributions
plot(
  my.dists[[1]]$seq.li,
  my.dists[[1]]$li.dense,
  type='l'
)
lines(
  my.dists[[2]]$seq.li,
  my.dists[[2]]$li.dense,
  col='red'
)
lines(
  my.dists[[3]]$seq.li,
  my.dists[[3]]$li.dense,
  lty=2,
  col='black'
)
lines(
  my.dists[[3]]$seq.li,
  my.dists[[3]]$li.dense,
  lty=2,
  col='red'
)


