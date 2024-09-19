source('scripts/solveSingleEffect.R')

mu <- 1e-6
Ne <- 3000
theta <- 4*Ne*mu
cost <- c(0.25,0.75)
h2 <- 0.5
my.L <- 1e5
tmp.bt <- seq(0.01,0.99,length.out = 9)
bt <- sort(c(tmp.bt,mean(tail(tmp.bt,2)),0.999))
my.thrs <- my.L * (1 - bt)


my.table <- expand.grid(cost,my.thrs,my.L,theta,h2,Ne)
colnames(my.table) <- c('cost','thr','target.size','theta','h2','Ne')
my.rho <- my.table[,'thr']/(2*my.table[,'target.size'])
my.b <- 1-2*my.rho
my.gamma <- log((1+my.b)/(1-my.b))
env.sd <- sqrt(theta*my.L*my.b/my.gamma*(1-h2)/h2)
tmp.table <- cbind(my.table,'rho'=my.rho,'b'=my.b,'gamma'=my.gamma,'env.sd'=env.sd)


## U <- 2*tmp.table[,'target.size']*mu
## BigGamma <- 2*Ne*tmp.table[,'cost']

norm.prev <- eqNormPrev(my.L,mu,my.b,Ne,h2,cost)
new.table <- cbind(tmp.table,'norm.prev'=norm.prev)
new.table





write.table(new.table,file='smallEffectInsensitivityParamTable.txt')



