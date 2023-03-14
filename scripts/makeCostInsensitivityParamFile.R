source('scripts/solveSingleEffect.R')
eqNormPrev <- function(deltaU,bigGamma,littleGamma,h2){
    1 - pnorm(dnorminv(bigGamma^-1 * sqrt(deltaU/2 * littleGamma/h2)))
}


mu <- 1e-6
Ne <- 5000
theta <- 4*Ne*mu
cost <- c(0.01,0.1,0.3,0.5,0.7,0.9,0.99)
h2 <- 0.5
my.L <- 1e5
my.thrs <- exp(seq(log(2000),log(my.L*0.95),length.out=6))
my.rho <- my.thrs/(2*my.L)
my.b <- 1-2*my.rho
my.gamma <- log((1+my.b)/(1-my.b))
my.table <- expand.grid(cost,my.thrs1,my.L,my.rho,my.gamma,theta,h2)
colnames(my.table) <- c('cost','thr','target.size','rho','gamma','theta','h2')

norm.prev <- eqNormPrev(2*my.table[,'target.size']*mu*my.b,2*Ne*my.table[,'cost'],my.table[,'gamma'],h2)
new.table <- cbind(my.table,'norm.prev'=norm.prev)


write.table(new.table,file='costInsensitivityParamFile.txt')



