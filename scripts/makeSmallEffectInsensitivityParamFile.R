source('scripts/solveSingleEffect.R')

mu <- 1e-6
Ne <- 500
theta <- 4*Ne*mu
cost <- c(seq(0.01,0.97,by=0.06),1)
h2 <- 0.5
my.L <- 1e5
my.thrs <- c(exp(seq(log(4000),log(2*my.L*0.4),length.out=6)))


my.table <- expand.grid(cost,my.thrs,my.L,theta,h2,Ne)
colnames(my.table) <- c('cost','thr','target.size','theta','h2','Ne')
my.rho <- my.table[,'thr']/(2*my.table[,'target.size'])
my.b <- 1-2*my.rho
my.gamma <- log((1+my.b)/(1-my.b))
env.sd <- sqrt(theta*my.L*my.b/my.gamma*(1-h2)/h2)
tmp.table <- cbind(my.table,'rho'=my.rho,'b'=my.b,'gamma'=my.gamma,'env.sd'=env.sd)


U <- 2*tmp.table[,'target.size']*mu
BigGamma <- 2*Ne*tmp.table[,'cost']

norm.prev <- eqNormPrev(U*my.b,BigGamma,my.gamma,h2)
new.table <- cbind(tmp.table,'norm.prev'=norm.prev)
new.table





write.table(new.table,file='smallEffectInsensitivityParamTable.txt')



