library('wesanderson')

gammabt = function(bt) log((1+bt)/(1-bt))
balpha = function(a) (exp(a)-1)/(exp(a)+1)
my.q = 1e-6
## my.mean = 1.098

bt.diff = function(my.mean,my.sd,my.bt){
  my.shape = my.mean^2/my.sd^2
  my.rate = my.mean / my.sd^2
  ## qgamma(my.q,my.shape,my.rate)
  ## balpha(my.mean)
  my.bt-(1/my.mean)*integrate(function(a) dgamma(a,my.shape,my.rate)*a*balpha(a),lower = qgamma(my.q,my.shape,my.rate),upper=qgamma(1-my.q,my.shape,my.rate))$value
}

bt.diff2 = function(my.mean,dist,my.bt){
    ## qgamma(my.q,my.shape,my.rate)
    ## balpha(my.mean)
    low = my.mean - dist/2
    high = my.mean + dist/2
    var = low*balpha(low)/2 + high*balpha(high)/2
    my.bt-(1/my.mean)*var
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


tmp = list()
my.gammas = list()
for ( j in 1:length(my.sds)){
  tmp[[j]] = list()
  my.gammas[[j]] = numeric()
  for ( i in 1:length(my.bts)){
      tmp[[j]][[i]] = uniroot(function(X) bt.diff2(X,dist = 1e-2,my.bts[i]),lower=se.gammas[i]/100,upper=10*se.gammas[i])
      my.gammas[[j]][i] = tmp[[j]][[i]]$root
  }
}

my.gammas[[1]]/se.gammas









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


pdf('figures/suppFigures/meanGammaVsbT_distn.pdf',width=20,height=15)
par(mfrow=c(2,2))

## 1
plot(
  x = my.bts,
  y = se.gammas,
  type = 'l' ,
  xlab = 'Mutational Asymmetry',
  ylab = 'Mean Population Scaled Selection Coefficient'
)
my.cols <- wes_palette('Zissou1',length(my.gammas))
for ( j in 1:length(my.gammas) ){
  plot.these = my.gammas[[j]] > 0.002
  lines(
    x = my.bts[plot.these],
    y = my.gammas[[j]][plot.these],
    lty = 2,
    lwd = 2,
    col = my.cols[j]
  )
}
legend(
  'topleft',
  legend = my.sds,
  col = my.cols,
  lty = 2,
  lwd = 2
)

## 2
plot(
  x = my.bts,
  y = 2*gammabt(my.bts)*my.bts,
  type = 'l' ,
  xlab = 'Mutational Asymmetry',
  ylab = 'Genetic Variance'
)
my.cols <- wes_palette('Zissou1',length(my.gammas))
for ( j in 1:length(my.gammas) ){
  plot.these = my.gammas[[j]] > 0.002
  lines(
    x = my.bts[plot.these],
    y = my.vars[[j]][plot.these],
    lty = 2,
    lwd = 2,
    col = my.cols[j]
  )
}


## 3
plot(
  NA,
  xlim = c(0,1),
  ylim = c(0,1.05),
  type = 'l' ,
  xlab = 'Mutational Asymmetry',
  ylab = 'Fold change in genetic variance relative to single effect model'
)
my.cols <- wes_palette('Zissou1',length(my.gammas))
for ( j in 1:length(my.gammas) ){
  plot.these = my.gammas[[j]] > 0.002
  lines(
    x = my.bts[plot.these],
    y = normed.vars[[j]][plot.these],
    lty = 2,
    lwd = 2,
    col = my.cols[j]
  )
}

plot(
  NA,
  xlim = c(0,1),
  ylim = c(-0.1,0.1),
  type = 'l' ,
  xlab = 'Mutational Asymmetry',
  ylab = 'Fold change in genetic variance relative to single effect model'
)
for ( j in 1:length(my.gammas) ){
  plot.these = my.gammas[[j]] > 0.002
  lines(
    x = my.bts[plot.these],
    y = 2*my.gammas[[j]][plot.these]*my.bts[plot.these] - my.vars[[j]][plot.these],
    lty = 2,
    lwd = 2,
    col = my.cols[j]
  )
}
dev.off()





my.gammas = seq(0.01,100,length.out = 1000)
plot(
  my.gammas ,
  dgamma(my.gammas,shape = this.shape, rate = this.rate),
  type='l'
)






if(FALSE){
  
  
  
  
  plot(
    my.bts, 
    se.gammas*balpha(se.gammas),
    type='l'
  )
  
  plot(
    se.gammas, 
    balpha(se.gammas),
    type='l'
  )
}
