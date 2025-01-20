rm(list = ls())
setwd('/Users/jeremyberg/Documents/academics/LTM_simulation/')
library('wesanderson')

dnorminv = function(y)
  sqrt(-2 * log(sqrt(2 * pi) * y))
gammabt = function(bt)
  log((1 + bt) / (1 - bt))
balpha = function(a)
  ifelse(a < 100, (exp(a) - 1) / (exp(a) + 1), 1)
my.q = 1e-7
d0fa = function(a)
  a * (exp(2*a) - 1) / (exp(2*a) + 1)
d1fa = function(a)
  (4 * a * exp(2*a) + exp(4 * a) - 1) / ((1 + exp(2*a)) ^ 2)
d2fa = function(a)
  - 8 * exp(2*a) * (a * (exp(2*a) - 1) - exp(2*a) - 1 ) / ((1 + exp(2*a)) ^
                                                            3)


bt.diff = function(my.mean, my.sd, my.bt) {
  my.shape = my.mean ^ 2 / my.sd ^ 2
  my.rate = my.mean / my.sd ^ 2
  ## qgamma(my.q,my.shape,my.rate)
  ## balpha(my.mean)
  myFunc <- function(a)
    dgamma(a, my.shape, my.rate) * a * balpha(a)
  my.bt - (1 / my.mean) * integrate(
    myFunc,
    lower = qgamma(my.q, my.shape, my.rate),
    upper = qgamma(1 - my.q, my.shape, my.rate)
  )$value
}

bt.diff2 = function(my.mean, my.cv, my.bt) {
  my.shape = 1 / my.cv ^ 2
  my.rate = 1 / (my.cv ^ 2 * my.mean)
  my.lower = qgamma(my.q, my.shape, my.rate)
  my.upper = qgamma(1 - my.q, my.shape, my.rate)
  ## balpha(my.mean)
  myFunc <- function(a)
    dgamma(a, my.shape, my.rate) * a * balpha(a)
  my.bt - (1 / my.mean) * integrate(myFunc, lower = my.lower, upper = my.upper)$value
}




my.bts = seq(1e-3, 0.999, length.out = 1000)
my.cvs = c(1e-1, 0.3, 1, 2, 3)

se.gammas = log((1 + my.bts) / (1 - my.bts))


Ne = 1e4
U = 0.2
h2 = 1 / 2
cost = 1 / 2
tmp = list()
my.gammas = list()
new.vars = list()
fold.var.change = list()
new.astd = list()
new.std.dens = list()
new.prev = list()
fold.prev.change = list()
fold.gamma.change = list()
se.vars = 8 * Ne * U * my.bts / se.gammas
se.astd = sqrt (h2 / se.vars)
se.std.dens = U * my.bts * se.astd / (h2 * cost)
se.prev = 1 - pnorm(dnorminv(se.std.dens))
for (j in 1:length(my.cvs)) {
  tmp[[j]] = list()
  my.gammas[[j]] = rep(NA, length(my.bts))
  new.vars[[j]] = rep(NA, length(my.bts))
  new.astd[[j]] = rep(NA, length(my.bts))
  new.std.dens[[j]] = rep(NA, length(my.bts))
  new.prev[[j]] = rep(NA, length(my.bts))
  fold.gamma.change[[j]] = rep(NA, length(my.bts))
  fold.var.change[[j]] = rep(NA, length(my.bts))
  fold.prev.change[[j]] = rep(NA, length(my.bts))
  for (i in 1:length(my.bts)) {
    my.lower = 0.1 * se.gammas[i]
    my.upper =  1000 * se.gammas[i]
    ## bt.diff2(my.lower,my.cv = my.cvs[j],my.bts[i])
    ## bt.diff2(my.upper,my.cv = my.cvs[j],my.bts[i])
    tmp[[j]][[i]] = uniroot(
      function(X)
        bt.diff2(X, my.cv = my.cvs[j], my.bts[i]),
      lower = my.lower,
      upper = my.upper
    )
    #tmp[[j]][[i]] = uniroot(function(X) bt.diff(X,my.sd = my.sds[j],my.bts[i]),lower=my.lowers[j],upper=12)
    if (tmp[[j]][[i]]$root < 0.002)
      next
    my.gammas[[j]][i] = tmp[[j]][[i]]$root
    new.vars[[j]][i] = 8 * Ne * U * my.bts[i] / my.gammas[[j]][i]
    fold.var.change[[j]][i] = new.vars[[j]][i] / se.vars[i]
    new.astd[[j]][i] = sqrt (h2 / new.vars[[j]][i])
    new.std.dens[[j]][i] = U * my.bts[i] * new.astd[[j]][i] / (h2 * cost)
    new.prev[[j]][i] = 1 - pnorm(dnorminv(new.std.dens[[j]][i]))
    fold.gamma.change[[j]][i] = new.vars[[j]][i] / se.vars[i]
    fold.prev.change[[j]][i] = new.prev[[j]][i] / se.prev[i]
  }
}


{
  this.j <- 5
  this.i <- 990
  this.mean <- my.gammas[[this.j]][[this.i]]
  this.cv <- my.cvs[this.j]
  
  this.shape = 1 / this.cv ^ 2
  this.rate = 1 / (this.cv ^ 2 * this.mean)
  
  a <- seq(0, 6, length.out = 10000)
  par(mfrow = c(2, 3))
  plot(NA,
       xlim = c(0, 6),
       ylim = c(0, 2))
  lines(a,
        a * balpha(a),
        col = 'red',
        lty = 1,
        lwd = 3)
  lines(
    a,
    dgamma(a, this.shape, this.rate),
    lty = 1,
    lwd = 3,
    col = 'black'
  )
  abline(h = 0, lty = 3)
  
  plot(NA,
       xlim = c(0, 6),
       ylim = c(0, 2))
  lines(a,
        d1fa(a),
        col = 'red',
        lty = 1,
        lwd = 3)
  lines(
    a,
    dgamma(a, this.shape, this.rate),
    lty = 1,
    lwd = 3,
    col = 'black'
  )
  abline(h = 0, lty = 3)
  
  plot(NA,
       xlim = c(0, 6),
       ylim = c(0, 2))
  lines(a,
        d2fa(a),
        col = 'red',
        lty = 1,
        lwd = 3)
  lines(
    a,
    dgamma(a, this.shape, this.rate),
    lty = 1,
    lwd = 3,
    col = 'black'
  )
  abline(h = 0, lty = 3)
  
  plot(NA,
       xlim = c(0, 6),
       ylim = c(0, 1))
  lines(a,
        a * balpha(a),
        col = 'red',
        lty = 1,
        lwd = 3)
  lines(
    a,
    pgamma(a, this.shape, this.rate),
    lty = 1,
    lwd = 3,
    col = 'black'
  )
  abline(h = 0, lty = 3)
  abline(h = 1, lty = 3)
  
  plot(NA,
       xlim = c(0, 6),
       ylim = c(0, 1))
  lines(a,
        d1fa(a),
        col = 'red',
        lty = 1,
        lwd = 3)
  lines(
    a,
    pgamma(a, this.shape, this.rate),
    lty = 1,
    lwd = 3,
    col = 'black'
  )
  abline(h = 0, lty = 3)
  abline(h = 1, lty = 3)
  
  plot(NA,
       xlim = c(0, 6),
       ylim = c(0, 1))
  lines(a,
        d2fa(a),
        col = 'red',
        lty = 1,
        lwd = 3)
  lines(
    a,
    pgamma(a, this.shape, this.rate),
    lty = 1,
    lwd = 3,
    col = 'black'
  )
  abline(h = 0, lty = 3)
  abline(h = 1, lty = 3)
  
  
}


gamma.inflect = uniroot(d2fa,interval=c(1,2))$root
bt.inflect = uniroot(function(x) gamma.inflect - 0.5*log((1+x)/(1-x)),lower=0,upper=1)$root

cex.axis = 1.9
cex.lab = 1.7

pdf('figures/suppFigures/meanGammaVsbT_distn.pdf',width=20,height=15)
par(mfrow=c(2,3))
op = par(mar=c(5,6,4,1)+0.1)
## 1
plot(
  x = my.bts,
  y = se.gammas,
  type = 'l' ,
  xlab = '',
  ylab = '',
  cex.axis = cex.axis,
  cex.lab = cex.lab,
  ylim = c(0,8),
  lwd = 2, 
  lty = 1
)
mtext(side=2,text='Mean Population Scaled Selection Coefficient',cex=cex.lab,line=3.3)
mtext(side=1,text='Mutational Asymmetry (b_T)',cex=cex.lab,line=3)
my.cols <- wes_palette('Zissou1',length(my.gammas))
for ( j in 1:length(my.gammas) ){
  plot.these = !is.na(my.gammas[[j]]) & my.gammas[[j]] < 10
  lines(
    x = my.bts[plot.these],
    y = my.gammas[[j]][plot.these],
    lty = 2,
    lwd = 3,
    col = my.cols[j]
  )
}
legend(
  x = 0.05,
  y=7,
  legend = c(0,my.cvs),
  col = c('black',my.cols),
  lty = c(1,rep(2,length(my.cvs))),
  lwd = c(2,rep(3,length(my.cvs))),
  bty = 'n',
  title = 'Effect dist\'n coef. of var.',
  cex = cex.lab + 0.4
)
abline(v = bt.inflect, lty = 2 , lwd = 2)
par(op)


## 2
op = par(mar=c(5,6,4,1)+0.1)
y.max <- max(unlist(new.vars),na.rm=T)
plot(
  x = my.bts,
  y = se.vars ,
  type = 'l' ,
  xlab = '',
  ylab = '',
  cex.axis = cex.axis,
  cex.lab = cex.lab,
  ylim = c(0,y.max),
  lwd = 2
)
mtext(side=2,text='Genetic Variance',cex=cex.lab,line=3.3)
mtext(side=1,text='Mutational Asymmetry (b_T)',cex=cex.lab,line=3)
my.cols <- wes_palette('Zissou1',length(my.gammas))
for ( j in 1:length(my.gammas) ){
  plot.these = !is.na(my.gammas[[j]]) & my.gammas[[j]] < 10
  lines(
    x = my.bts[plot.these],
    y = new.vars[[j]][plot.these],
    lty = 2,
    lwd = 3,
    col = my.cols[j]
  )
}
abline(v = bt.inflect, lty = 2 , lwd = 2)
par(op)



## 3
my.ymax = 0.006
plot(
  my.bts,
  se.prev,
  xlim = c(0,1),
  ylim = c(0,my.ymax),
  type = 'l' ,
  xlab = '',
  ylab = '',
  cex.axis = cex.axis,
  cex.lab = cex.lab
)
mtext(side=2,text='Prevalence',cex=cex.lab,line=3.3)
mtext(side=1,text='Mutational Asymmetry (b_T)',cex=cex.lab,line=3)
my.cols <- wes_palette('Zissou1',length(my.gammas))
for ( j in 1:length(my.gammas) ){
  plot.these = !is.na(my.gammas[[j]]) & my.gammas[[j]] < 10
  lines(
    x = my.bts[plot.these],
    y = new.prev[[j]][plot.these],
    lty = 2,
    lwd = 3,
    col = my.cols[j]
  )
}
abline(v = bt.inflect, lty = 2 , lwd = 2)



op = par(mar=c(5,6,4,1)+0.1)
## 4
plot(
  NA,
  xlim = c(0,1),
  ylim = c(0,10),
  type = 'l' ,
  xlab = '',
  ylab = '',
  cex.axis = cex.axis,
  cex.lab = cex.lab
)
mtext(side=2,text='Fold change in mean scaled selection',cex=cex.lab,line=4.2)
mtext(side=2,text='coefficient relative to single effect model',cex=cex.lab,line=2.4)
mtext(side=1,text='Mutational Asymmetry (b_T)',cex=cex.lab,line=3)
my.cols <- wes_palette('Zissou1',length(my.gammas))
abline(h = 1, lty = 3 , lwd = 1)
for ( j in 1:length(my.gammas) ){
  plot.these = !is.na(my.gammas[[j]]) & my.gammas[[j]] < 10
  lines(
    x = my.bts[plot.these],
    y = fold.gamma.change[[j]][plot.these],
    lty = 2,
    lwd = 3,
    col = my.cols[j]
  )
}
abline(v = bt.inflect, lty = 2 , lwd = 2)
par(op)


op = par(mar=c(5,6,4,1)+0.1)
## 5
plot(
  NA,
  xlim = c(0,1),
  ylim = c(0,10),
  type = 'l' ,
  xlab = '',
  ylab = '',
  cex.axis = cex.axis,
  cex.lab = cex.lab
)
mtext(side=2,text='Fold change in genetic variance',cex=cex.lab,line=4.3)
mtext(side=2,text='relative to single effect model',cex=cex.lab,line=2.5)
mtext(side=1,text='Mutational Asymmetry (b_T)',cex=cex.lab,line=3)
my.cols <- wes_palette('Zissou1',length(my.gammas))
abline(h = 1, lty = 3 , lwd = 1)
for ( j in 1:length(my.gammas) ){
  plot.these = !is.na(my.gammas[[j]]) & my.gammas[[j]] < 10
  lines(
    x = my.bts[plot.these],
    y = fold.var.change[[j]][plot.these],
    lty = 2,
    lwd = 3,
    col = my.cols[j]
  )
}
abline(v = bt.inflect, lty = 2 , lwd = 2)
par(op)



op = par(mar=c(5,6,4,1)+0.1)
## 6
plot(
  NA,
  xlim = c(0,1),
  ylim = c(0,1.5),
  type = 'l' ,
  xlab = '',
  ylab = '',
  cex.axis = cex.axis,
  cex.lab = cex.lab
)
mtext(side=2,text='Fold change in prevalence',cex=cex.lab,line=4.3)
mtext(side=2,text='relative to single effect model',cex=cex.lab,line=2.5)
mtext(side=1,text='Mutational Asymmetry (b_T)',cex=cex.lab,line=3)
abline(h = 1, lty = 3 , lwd = 1)
my.cols <- wes_palette('Zissou1',length(my.gammas))
for ( j in 1:length(my.gammas) ){
  plot.these = !is.na(my.gammas[[j]]) & my.gammas[[j]] < 10
  lines(
    x = my.bts[plot.these],
    y = fold.prev.change[[j]][plot.these],
    lty = 2,
    lwd = 3,
    col = my.cols[j]
  )
}
abline(v = bt.inflect, lty = 2 , lwd = 2)
par(op)
dev.off()


if(FALSE){
  
  gamma.inflect = uniroot(d2fa,interval=c(1,2))$root
  a <- seq(0,6,length.out=1000)
  par(mfrow=c(1,3))  
  plot(a,
       d0fa(a),
       type='l',
       bty = 'n',
       xlim = c(0,6)
       )
  abline(v = gamma.inflect, lty = 2)
  plot(a,
       d1fa(a),
       type='l',
       bty = 'n',
       xlim = c(0,6)
  )
  abline(v = gamma.inflect, lty = 2)
  plot(a,
       d2fa(a),
       type='l',
       bty = 'n',
       xlim = c(0,6)
  )
  abline(v = gamma.inflect, lty = 2)
  
  
small.one = TRUE
if(small.one){
  this.one = min(which(my.gammas[[4]] > 0.002))
  this.mean = my.gammas[[4]][this.one]
  this.shape = this.mean^2/tail.sd^2
  this.rate = this.mean / tail.sd^2
} else {
  this.one = max(which(my.gammas[[4]] > 0.002))
  this.mean = my.gammas[[4]][this.one]
  this.shape = this.mean^2/tail.sd^2
  this.rate = this.mean / tail.sd^2
}




plot.gammas = seq(0.01,10,length.out = 1000)
plot(
  plot.gammas ,
  dgamma(plot.gammas,shape = this.shape, rate = this.rate),
  type='l'
)
integrate(
  f=function(a) dgamma(a,shape = this.shape, rate = this.rate)*a*(exp(a)-1)/(exp(a) + 1),
  lower=1e-6,
  upper=15
)$value/this.mean
(exp(this.mean)-1)/(exp(this.mean) + 1)




  
  
  plot(
    plot.gammas, 
    plot.gammas*balpha(plot.gammas),
    type='l'
  )
  
  plot(
    plot.gammas, 
    balpha(plot.gammas),
    se.gammas, 
    balpha(se.gammas),
    type='l'
  )
}
