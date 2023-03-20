source('scripts/solveSingleEffect.R')
eqNormPrev <- function(deltaU,bigGamma,littleGamma,h2){
    1 - pnorm(dnorminv(bigGamma^-1 * sqrt(deltaU/2 * littleGamma/h2)))
}


mu <- 1e-6
Ne <- 5000
theta <- 4*Ne*mu
cost <- seq(0.01,0.99,length.out=6)
h2 <- 0.5
my.L <- 1e5
my.thrs <- c(exp(seq(log(2000),log(36186.802),length.out=4)),2*my.L*0.35)


my.table <- expand.grid(cost,my.thrs,my.L,theta,h2,Ne)

colnames(my.table) <- c('cost','thr','target.size','theta','h2','Ne')
my.rho <- my.table[,'thr']/(2*my.table[,'target.size'])
my.b <- 1-2*my.rho
my.gamma <- log((1+my.b)/(1-my.b))
env.sd <- sqrt(theta*my.L*my.b*(1-h2)/h2)
tmp.table <- cbind(my.table,'rho'=my.rho,'b'=my.b,'gamma'=my.gamma,'env.sd'=env.sd)
norm.prev <- eqNormPrev(2*tmp.table[,'target.size']*mu*my.b,2*Ne*tmp.table[,'cost'],tmp.table[,'gamma'],h2)
new.table <- cbind(tmp.table,'norm.prev'=norm.prev)
new.table





write.table(new.table,file='costInsensitivityParamTable.txt')



