library("MetBrewer")
library('numDeriv')
prevs = c(0.03,0.01,0.003,0.001)
thr = sapply(prevs, function(PREV) qnorm(1-PREV))
phits = dnorm(thr)
max.a = 5
min.a = rep(0.001,4)
my.a = mapply(function(MAXA,MINA) seq(MINA,MAXA,length.out=1000), MINA = min.a, MAXA = max.a,SIMPLIFY = FALSE)
N = 3000
cost = 1
four.cols = met.brewer('Hiroshige',4)[c(1,3,2,4)]


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
sigmaaSmall = function(a,thr,nc){
    approxDel = a*dnorm(thr)
    a^2 * bgamma(4*nc*approxDel) / (2*nc*approxDel)
}
sigmaa2 = function(a,del,nc){
  a^2 * bgamma(4*nc*del) / (2*nc*del)
}

my.del = mapply( function(THR,A) normDel(A,THR), THR = thr, A = my.a)
## my.vars = mapply(function(THR,A) sigmaa(A,THR,N*cost), THR = thr, A = my.a)

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
pdf('figures/suppFigures/largeEffectSizeVarianceRelationship.pdf',width=12,height=12)
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
        lty = 1,
        lwd = 2,
        col = four.cols[i]
    )
    lines(
        my.a[[i]],
        my.a[[i]]*phits[i],
        lty = 2 ,
        lwd = 1/2,
        col = four.cols[i]
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
    lty = 1,
    lwd = 2,
    col = four.cols,
    bty = 'n'
)
legend(
    x = -0.1,
    y = 0.66,
    legend = c('Normal', 'Small effect approx'),
    lty = c(1,2),
    lwd = c(2,1/2),
    col = 1,
    bty = 'n'
)



my.costs <- c(1, 2 / 3, 1 / 3)
for (j in 1:length(my.costs)) {
  my.vars = mapply(function(THR, A)
    sigmaa(A, THR, N * my.costs[j]),
    THR = thr,
    A = my.a)
  small.approx.vars = mapply(function(THR,A)
      sigmaaSmall(A, THR, N * my.costs[j]),
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
          small.approx.vars[,i],
          lty = 2,
          lwd = 1/2,
          col = four.cols[i])
      lines(my.a[[i]],
          my.vars[, i],
          lty = 1,
          lwd = 2,
          col = four.cols[i])
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


#################################
## how many sites do you need? ##
#################################

##my.costs <- c(1, 2 / 3, 1 / 3)
i = 1
Ne = 3000
my.ratios <- lapply(thr,
       function(THR) {
         tmp.bt <- 0.9
         min.y <- log((1+tmp.bt)/(1-tmp.bt))
         min.del <- min.y / (4 * Ne)
         my.prev <- 1-pnorm(THR)
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
         h2l.ratio <- delsq.ratio * asq.ratio * bgamma.ratio
         list(these.dels,h2l.ratio,these.as,delsq.ratio,asq.ratio,bgamma.ratio)
       }
)
my.dels <- lapply(my.ratios,function(X)X[[1]])
my.h2ratios <- lapply(my.ratios,function(X)X[[2]])
my.as <- lapply(my.ratios,function(X)X[[3]])
par(mfrow=c(1,2))
plot(
  NA,
  xlim = c(0, 1),
  ylim = c(0, 1.2),
  type = 'l',
  xlab = 'PAR of large effect allele',
  ylab = 'Relative contribution to liability scale h2',
  bty = 'n'
)
for(j in 1:length(my.dels)){
lines(my.dels[[j]],
      my.h2ratios[[j]],
      lty = 1,
      lwd = 2,
      col = four.cols[j])
}
legend(
  'topright',
  bty = 'n',
  lty = 1 ,
  lwd = 2,
  col = four.cols,
  legend = prevs,
  title = 'Prevalence'
)
plot(
  NA,
  xlim = c(0, max(unlist(my.as))),
  ylim = c(0, 1.2),
  type = 'l',
  xlab = 'Standardized liability effect of large effect alleles',
  ylab = 'Relative contribution to liability scale h2',
  bty = 'n'
)
for(j in 1:length(my.dels)){
  lines(my.as[[j]],
        my.h2ratios[[j]],
        lty = 1,
        lwd = 2,
        col = four.cols[j])
}









## how does Poisson shape change things?
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
thr = sapply(prevs, function(PREV) qnorm(1-PREV))
max.a = thr + max(thr)
min.a = rep(0.001,4)
N = 3000
my.a = mapply(function(MAXA,MINA) seq(MINA,MAXA,length.out=1000), MINA = min.a, MAXA = max.a,SIMPLIFY = FALSE)
alt.del = mapply( function(THR,A) normDel(A,THR), THR = thr, A = my.a,SIMPLIFY=FALSE)
alt.vars = mapply(function(DEL,A,COST) sigmaa2(A,DEL,N*COST), DEL = alt.del, A = my.a, COST = tuned.costs)





pdf('figures/suppFigures/largeEffectSizeVarianceRelationshipPoisson.pdf',width=12,height=12)
par(mfrow=c(1,2))
op = par(mar=c(5,5,4,4)+0.1)


## plot risk effects
plot(
  NA,
  type='l',
  xlim=c(0,6),
  ylim=c(0,1),
  col = four.cols[1],
  lwd=2
)
lines(
  my.dists[[1]]$std.alpha.seq , 
  deltals[[1]],
  col = four.cols[1],
  lwd=2
)
lines(
  my.dists[[2]]$std.alpha.seq,
  deltals[[2]],
  col = four.cols[2],
  lwd=2
)
lines(
  my.dists[[3]]$std.alpha.seq,
  deltals[[3]],
  col = four.cols[3],
  lwd=2
)
lines(
  my.dists[[4]]$std.alpha.seq,
  deltals[[4]],
  col = four.cols[4],
  lwd=2
)


for( i in 1:ncol(alt.vars)){
  lines(
    x = my.a[[i]],
    y = alt.del[[i]],
    col = four.cols[i],
    lwd = 2,
    lty = 2
  )
}


par(op)
## plot variance contributions
plot(
  NA,
  type='l',
  xlim=c(0,6),
  ylim=c(0,max(c(my.dists[[3]]$vars,alt.vars))),
  col = four.cols[1],
  lwd=2
)
lines(
  my.dists[[1]]$std.alpha.seq , 
  my.dists[[1]]$vars,
  col = four.cols[1],
  lwd=2
)
lines(
  my.dists[[2]]$std.alpha.seq,
  my.dists[[2]]$vars,
  col = four.cols[2],
  lwd=2
)
lines(
  my.dists[[3]]$std.alpha.seq,
  my.dists[[3]]$vars,
  col = four.cols[3],
  lwd=2
)
lines(
  my.dists[[4]]$std.alpha.seq,
  my.dists[[4]]$vars,
  col = four.cols[4],
  lwd=2
)
for( i in 1:ncol(alt.vars)){
  lines(
    x = my.a[[i]],
    y = alt.vars[,i],
    col = four.cols[i],
    lwd = 2,
    lty = 2
  )
}

dev.off()




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


