rm(list=ls())
source('scripts/solveTwoEffect.R')
i=1
j=800
L <- 1e5
my.gs <- 1-seq(10/L,0.1,length.out=100)
## param.table <- read.table("bgsStuff/BGS_twoEffectPrevInsensitivityParamsTable.txt",header=TRUE)##[c(1,20),]
as <- 1


##these.ones <- c(12,29,34)
##param.table <- param.table[these.ones,]


## load(file='tmp.solns.Robj')
solns <- get(load(file='tmp.fixed.deltal.solns.Robj'))
solns[[i]][[j]]
Ne <- 5000
cost <- 1/2
u <- 1e-7

al <- solns[[i]][[j]]['al']
Ll <- L * (1-my.gs[j]) / my.gs[j]
Ls <- L
bs <- param.table$bs[i]

pois.deltal <- solns[[i]][[j]]['deltal']
norm.deltal <- solns[[i]][[j]]['norm.deltal']
ys <- 0.5*log((1 + bs) / (1 - bs))
raw.Vas <- 4 * Ne * u * Ls * as ^ 2 * bs / ys
raw.ft <- 1 / (2*Ne*cost) * ys

mean.nl <- 2 * Ll * u / (pois.deltal*cost)
pois.raw.Val <- al^2 * mean.nl
norm.raw.Val <- 2 * al^2 * Ll * u / (norm.deltal*cost)
Ve <- param.table$Ve[i]/3
norm.raw.Vt <- as.vector(raw.Vas + norm.raw.Val + Ve)
pois.raw.Vt <- as.vector(raw.Vas + pois.raw.Val + Ve)
norm.sd <- sqrt(1 - pois.raw.Val / pois.raw.Vt)
std.al <- al / sqrt( pois.raw.Vt )
pois.tstar <- solns[[i]][[j]]['tstar']
pois.std.ft <- raw.ft * sqrt( pois.raw.Vt )


dPoisConv(pois.tstar,
          mean.nl,
          norm.sd,
          alphal = std.al,
          risk.allele = FALSE)


my.xs <- seq(-4,5,length.out=2000)
pois.dens <- sapply (my.xs,
        function(X)
          dPoisConv(X,
                    mean.nl,
                    norm.sd,
                    alphal = std.al,
                    risk.allele = FALSE))


norm.dens <- dnorm(my.xs)
norm.tstar <- solns[[i]][[j]]['norm.tstar']

offset <- norm.tstar*sqrt(norm.raw.Vt) - pois.tstar*sqrt( pois.raw.Vt )

plot(my.xs*sqrt( pois.raw.Vt )+offset,
     pois.dens/sqrt( pois.raw.Vt ),
     type='l',
     xlim = c(-70,70))
lines(my.xs*sqrt(norm.raw.Vt),
      norm.dens/sqrt(norm.raw.Vt),
      col='blue'
)
abline(v = norm.tstar*sqrt(norm.raw.Vt) , lty = 2, col='blue')
##abline(v = pois.tstar*sqrt( pois.raw.Vt ) , lty = 2 )
abline(h = raw.ft , col = 'red')









