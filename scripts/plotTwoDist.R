rm(list=ls())
source('scripts/solveTwoEffect.R')
i=1
L <- 1e5
##my.gs <- 1-seq(10/L,0.1,length.out=100)
## param.table <- read.table("bgsStuff/BGS_twoEffectPrevInsensitivityParamsTable.txt",header=TRUE)##[c(1,20),]
as <- 1


##these.ones <- c(12,29,34)
##param.table <- param.table[these.ones,]


## load(file='tmp.solns.Robj')
solns <- get(load(file='tmp.fixed.deltal.solns.Robj'))
##j <- length(solns[[i]])

j <- which.max(sapply(solns[[1]],function(X) X['prev']))
## j <- length(solns[[1]])
j <- 8000
my.soln <- unlist(solns[[i]][[j]])
Ne <- 5000
cost <- 1/2
u <- 1e-7
al <- my.soln['al']
Ll <- L * (1-my.soln['gs']) / my.soln['gs']
Ls <- L
bs <- my.soln['bs']
ps <- (1 - bs) / 2
my.t <- 2*Ls*ps
pois.deltal <- my.soln['deltal']
norm.deltal <- my.soln['deltal']
ys <- 0.5*log((1 + bs) / (1 - bs))
raw.Vas <- 4 * Ne * u * Ls * as ^ 2 * bs / ys
raw.ft <- 1 / (2*Ne*cost) * ys
mean.nl <- 2 * Ll * u / (pois.deltal*cost)
pois.raw.Val <- al^2 * mean.nl
norm.raw.Val <- 2 * al^2 * Ll * u / (norm.deltal*cost)
Ve <- my.soln['Ve']
norm.raw.Vt <- as.vector(raw.Vas + norm.raw.Val + Ve)
pois.raw.Vt <- as.vector(raw.Vas + pois.raw.Val + Ve)
norm.sd <- sqrt(1 - pois.raw.Val / pois.raw.Vt)
std.al <- al / sqrt( pois.raw.Vt )
pois.tstar <- my.soln['tstar']
pois.std.ft <- raw.ft * sqrt( pois.raw.Vt )
my.soln['prev']
my.soln['norm.prev']


##tiff('figures/paperFiguresForRealThisTime/PoissonFigure8.tiff', width = 6.5 , height = 6.5*2/3 , units = 'in' , res = 600 )
pdf('figures/paperFiguresForRealThisTime/PoissonFigure8.pdf', width = 18 , height = 18*2/3 )
par(mfrow=c(2,3))
op <- par(mar = c(5,5,3,1) + 0.2)
large.xs <- 0:3
pois.pmf <- dpois(large.xs,lambda = mean.nl )
plot(
    x = large.xs*al ,
    y = pois.pmf ,
    type = 'n',
    xlab = '',
    ylab = '',
    xaxt = 'n',
    yaxt = 'n',
    ylim = c(0,1)
)
mtext('A', side = 3 , line = 1 , at = 0 , cex = 1.7 )
axis(1, cex.axis = 1.8 )
axis(2, cex.axis = 1.8 )
points(
    x = large.xs*al ,
    y = pois.pmf ,
    pch = 21 ,
    bg = 'darkgreen',
    col = 'black',
    cex = 1.6
)
mtext('Probability',side = 2 , line = 3, cex = 1.6 )
mtext('Large Effect Liability',side = 1 , line = 3 , cex = 1.6 )
abline(h = 0 , lty = 2 )
par(op)



my.xs <- seq(-4,7,length.out=2000)
pois.dens <- sapply (my.xs,
        function(X)
          dPoisConv(X,
                    mean.nl,
                    norm.sd,
                    alphal = std.al,
                    risk.allele = FALSE))
norm.dens <- dnorm(my.xs)
norm.tstar <- my.soln['norm.tstar']
my.norm.mean <- uniroot(function(X) dnorm(my.t, X , sqrt( pois.raw.Vt ) ) - raw.ft , lower = 34900, upper = 35000 )$root

pois.raw.t <- pois.tstar*sqrt( pois.raw.Vt )
norm.raw.t <- norm.tstar*sqrt(norm.raw.Vt)
pois.raw.xs <- my.xs *sqrt(pois.raw.Vt)
norm.plotting.offset <- (my.t - my.norm.mean)
norm.plotting.mean <- pois.raw.t - norm.plotting.offset


norm.raw.dens <- dnorm(pois.raw.xs, norm.plotting.mean,sqrt(pois.raw.Vt))

norm.raw.xs <- my.xs *sqrt(norm.raw.Vt)
pois.raw.dens <- pois.dens / sqrt(pois.raw.Vt)

## norm.raw.dens <- norm.dens / sqrt(norm.raw.Vt)
offset <- pois.raw.t - norm.raw.t







op <- par(mar = c(5,7,3,1) + 0.2)
plot.x.at <- seq(-150,50,by=50) + pois.raw.t 
plot.x.labs <- seq(-150,50,by=50) + my.t
plot(NA,
     xlim = c(pois.raw.t - 150,pois.raw.t + 50),
     ylim = c(0,max(c(norm.raw.dens,pois.raw.dens))),
     xaxt = 'n',
     yaxt = 'n',
     ylab = '',
     xlab = ''
     )
mtext('B', side = 3 , line = 1 , at = pois.raw.t - 150 , cex = 1.7 )
axis(1, at = plot.x.at , labels = plot.x.labs , cex.axis = 1.6 )
axis(2, cex.axis = 1.6 )
lines(pois.raw.xs,
      pois.raw.dens,
      col = adjustcolor('black' , alpha.f = 1 ), 
      lwd = 2
      )
lines(pois.raw.xs,
      norm.raw.dens,
      col = adjustcolor( 'black' , alpha.f = 0.4 ) ,
      lwd = 1
      )
abline( v = pois.raw.t , col = 'black' , lty = 3 )
mtext('Probability Density',side = 2 , line = 4, cex = 1.6 )
mtext('Total Liability',side = 1 , line = 3 , cex = 1.6 )
par(op)




op <- par(mar = c(5,7,3,1) + 0.2)
plot.x.at <- seq(-10,20,by=5) + pois.raw.t 
plot.x.labs <- seq(-10,20,by=5) + my.t
plot(NA,
     xlim = c(pois.raw.t - 10,pois.raw.t + 20),
     ylim = c(0,0.001),
     xaxt = 'n',
     yaxt = 'n',
     xlab = '',
     ylab = ''
     )
mtext('C', side = 3 , line = 1 , at = pois.raw.t - 10 , cex = 1.7 )
axis(1, at = plot.x.at , labels = plot.x.labs , cex.axis = 1.6)
axis(2, cex.axis = 1.6)
lines(pois.raw.xs,
      pois.raw.dens,
      col = adjustcolor('black' , alpha.f = 1 ),
      lwd = 2 ,
      
      )
lines(pois.raw.xs,
      norm.raw.dens,
      col = adjustcolor( 'black' , alpha.f = 0.4 ) ,
      lwd = 1
      )
mtext('Probability Density',side = 2 , line = 4, cex = 1.6 )
mtext('Total Liability',side = 1 , line = 3 , cex = 1.6 )
abline( v = pois.raw.t , col = 'black' , lty = 3 )
abline( h = raw.ft , lty = 2 )
par(op) 
     


output <- get(load(file='fixed.deltal.output.Robj'))
op <- par(mar = c(8,5,1,1) + 0.2)
i=1
plot(
    output[[i]][,'Ll']*u,
    output[[i]][,'prev'],
    type = 'l',
    col = 'black',
    lwd = 2,
    xaxt = 'n',
    yaxt = 'n',
    ylim = c(0,0.003),
    ylab = '',
    xlab = ''
)
mtext('D', side = 3 , line = 1 , at = 0 , cex = 1.7 )
axis(1 , cex.axis = 1.6 )
axis(2 , cex.axis = 1.6 , at = 0:3/1000 )
lines(
    output[[i]][,'Ll']*u,
    output[[i]][,'naive.norm.prev'],
    col=adjustcolor('black', alpha.f = 0.4 ),
    lwd = 1 
)
mtext('Prevalence',side = 2 , line = 3, cex = 1.6 )
mtext('Total haploid mutation rate for large effects',side = 1 , line = 4 , cex = 1.6 )
abline( v = output[[i]][j,'Ll']*u , lty = 2 , lwd = 1.5 )
## points(
##     x = output[[i]][j,'Ll']*u,
##     y = output[[i]][j,'prev'],
##     pch = 21 ,
##     bg = 'darkblue',
##     col = 'black' ,
##     cex = 1.6
##     )
par(op)




op <- par(mar = c(8,7,1,1) + 0.2)
plot(
    output[[i]][,'Ll']*u,
    output[[i]][,'pgal'],
    type = 'l',
    ylim = c(0,1),
    xlab = '',
    ylab = '',
    xaxt = 'n',
    yaxt = 'n' ,
    lwd = 2 ,
    col = adjustcolor( 'black' , alpha.f = 1 )
)
axis(side = 1 , cex.axis = 1.6 )
axis(side = 2 , at = 0:5/5, cex.axis = 1.6 )
lines(
    output[[i]][,'Ll']*u,
    output[[i]][,'pgal'],
    col=adjustcolor('black', alpha.f = 0.4 )
)
abline( v = output[[i]][j,'Ll']*u , lty = 2 , lwd = 1.5 )
## points(
##     x = output[[i]][j,'Ll']*u,
##     y = output[[i]][j,'pgal'],
##     pch = 21 ,
##     bg = 'darkblue',
##     col = 'black' ,
##     cex = 1.6
## )
mtext('E', side = 3 , line = 1 , at = 0 , cex = 1.7 )
mtext('Proportion of additive genetic variance' , side = 2 , line = 5 , cex = 1.6 )
mtext('for liability that is due to large effect loci' , side = 2 , line = 3 , cex = 1.6 )
mtext('Total haploid mutation rate for large effects',side = 1 , line = 4 , cex = 1.6 )
par(op)


op <- par(mar = c(8,7,1,1) + 0.2)
plot(
    output[[i]][,'Ll']*u,
    output[[i]][,'pgol'],
    type = 'l' ,
    xlab = '',
    ylab = '',
    xaxt = 'n',
    yaxt = 'n',
    lwd = 2
)
abline( v = output[[i]][j,'Ll']*u , lty = 2 , lwd = 1.5 )
## points(
##     x = output[[i]][j,'Ll']*u,
##     y = output[[i]][j,'pgol'],
##     pch = 21 ,
##     bg = 'darkblue',
##     col = 'black' ,
##     cex = 1.6
## )
mtext('F', side = 3 , line = 1 , at = 0 , cex = 1.7 )
mtext('Total haploid mutation rate for large effects',side = 1 , line = 4 , cex = 1.6 )
mtext('Proportion of additive genetic variance' , side = 2 , line = 5 , cex = 1.6 )
mtext('for risk that is due to large effect loci' , side = 2 , line = 3 , cex = 1.6 )
axis(1 , cex.axis = 1.6 )
axis(2 , cex.axis = 1.6 )
par(op)

dev.off()
