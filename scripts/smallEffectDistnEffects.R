library('wesanderson')

gammabt = function(bt) log((1+bt)/(1-bt))
balpha = function(a) (exp(a)-1)/(exp(a)+1)
my.q = 1e-6
d0fa = function(a) a*(exp(a)-1)/(exp(a)+1)
d1fa = function(a) (2*a*exp(a) + exp(2*a) - 1) / ((1 + exp(a))^2)
d2fa = function(a) - 2*exp(a)*(a*(exp(a)-1) - 2*(1 + exp(a) ) ) / ((1 + exp(a))^3)

bt.diff = function(my.mean,my.sd,my.bt){
  my.shape = my.mean^2/my.sd^2
  my.rate = my.mean / my.sd^2
  ## qgamma(my.q,my.shape,my.rate)
  ## balpha(my.mean)
  my.bt-(1/my.mean)*integrate(function(a) dgamma(a,my.shape,my.rate)*a*balpha(a),lower = qgamma(my.q,my.shape,my.rate),upper=qgamma(1-my.q,my.shape,my.rate))$value
}






my.bts = seq(1e-3,0.99,length.out=1000)
my.sds = c(1e-2,1e-1,0.5,0.8)
my.lowers = c(1e-6,4e-6,1e-4,1e-4)

se.gammas = log((1+my.bts)/(1-my.bts))
se.vars = 2*se.gammas*my.bts
this.gamma = tail(se.gammas,1)
tail.sd = tail(my.sds,1)
this.shape = this.gamma^2/tail.sd^2
this.rate = this.gamma / tail.sd^2



tmp = list()
my.gammas = list()
for ( j in 1:length(my.sds)){
  tmp[[j]] = list()
  my.gammas[[j]] = numeric()
  for ( i in 1:length(my.bts)){
    tmp[[j]][[i]] = uniroot(function(X) bt.diff(X,my.sd = my.sds[j],my.bts[i]),lower=my.lowers[j],upper=12)
    my.gammas[[j]][i] = tmp[[j]][[i]]$root
  }
}

qgamma(my.q,this.shape,this.rate)
qgamma(1-my.q,this.shape,this.rate)


my.vars = list()
normed.vars = list()
for ( j in 1:length(my.sds)){
  my.vars[[j]] = numeric()
  normed.vars[[j]] = numeric()
  for ( i in 1:length(my.bts)){
    this.shape = my.gammas[[j]][i]^2 / my.sds[j]^2
    this.rate = my.gammas[[j]][i] / my.sds[j]^2
    my.vars[[j]][i] = integrate(
      f = function(a){
        dgamma(a,this.shape,this.rate)*2*a*balpha(a)
      },
      lower = qgamma(my.q,this.shape,this.rate),
      upper = qgamma(1-my.q,this.shape,this.rate),
      rel.tol = .Machine$double.eps^0.5
    )$value
  }
  normed.vars[[j]] = my.vars[[j]]/se.vars
}

gamma.inflect = uniroot(d2fa,interval=c(2,4))$root
bt.inflect = uniroot(function(x) gamma.inflect - log((1+x)/(1-x)),lower=0,upper=1)$root

cex.axis = 1.5
cex.lab = 1.4

pdf('figures/suppFigures/meanGammaVsbT_distn.pdf',width=20,height=15)
par(mfrow=c(2,2))

## 1
plot(
  x = my.bts,
  y = se.gammas,
  type = 'l' ,
  xlab = 'Mutational Asymmetry (b_T)',
  ylab = 'Mean Population Scaled Selection Coefficient',
  cex.axis = cex.axis,
  cex.lab = cex.lab
)
my.cols <- wes_palette('Zissou1',length(my.gammas))
for ( j in 1:length(my.gammas) ){
  plot.these = my.gammas[[j]] > 0.002
  lines(
    x = my.bts[plot.these],
    y = my.gammas[[j]][plot.these],
    lty = 2,
    lwd = 3,
    col = my.cols[j]
  )
}
text(
  x = 0.1,
  y = 5.25,
  labels = "Standard deviation of"
)
text(
  x = 0.1,
  y = 5.1,
  labels = "effect size distribution"
)
legend(
  x = 0.05,
  y=5,
  legend = my.sds,
  col = my.cols,
  lty = 2,
  lwd = 3,
  bty = 'n'
)
abline(v = bt.inflect, lty = 2 , lwd = 2)


## 2
plot(
  x = my.bts,
  y = 2*gammabt(my.bts)*my.bts,
  type = 'l' ,
  xlab = 'Mutational Asymmetry (b_T)',
  ylab = 'Genetic Variance',
  cex.axis = cex.axis,
  cex.lab = cex.lab
)
my.cols <- wes_palette('Zissou1',length(my.gammas))
for ( j in 1:length(my.gammas) ){
  plot.these = my.gammas[[j]] > 0.002
  lines(
    x = my.bts[plot.these],
    y = my.vars[[j]][plot.these],
    lty = 2,
    lwd = 3,
    col = my.cols[j]
  )
}
abline(v = bt.inflect, lty = 2 , lwd = 2)



## 3
plot(
  NA,
  xlim = c(0,1),
  ylim = c(0,1.05),
  type = 'l' ,
  xlab = 'Mutational Asymmetry (b_T)',
  ylab = 'Fold change in genetic variance relative to single effect model',
  cex.axis = cex.axis,
  cex.lab = cex.lab
)
my.cols <- wes_palette('Zissou1',length(my.gammas))
for ( j in 1:length(my.gammas) ){
  plot.these = my.gammas[[j]] > 0.002
  lines(
    x = my.bts[plot.these],
    y = normed.vars[[j]][plot.these],
    lty = 2,
    lwd = 3,
    col = my.cols[j]
  )
}
abline(v = bt.inflect, lty = 2 , lwd = 2)

## 4
plot(
  NA,
  xlim = c(0,1),
  ylim = c(0,2),
  type = 'l' ,
  xlab = 'Mutational Asymmetry (b_T)',
  ylab = 'Fold change in threshold density relative to single effect model',
  cex.axis = cex.axis,
  cex.lab = cex.lab
)
my.cols <- wes_palette('Zissou1',length(my.gammas))
for ( j in 1:length(my.gammas) ){
  plot.these = my.gammas[[j]] > 0.002
  lines(
    x = my.bts[plot.these],
    y = normed.vars[[j]][plot.these]^(-1/2),
    lty = 2,
    lwd = 3,
    col = my.cols[j]
  )
}
abline(v = bt.inflect, lty = 2 , lwd = 2)
dev.off()


if(FALSE){
  
  
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
