rm(list = ls())
library('MetBrewer')
source('scripts/solveTwoEffect.R')
## main single effect params
L <- 1e6
Ne <- 5e3
u <- 1e-8
cost <- 1 / 2 
as <- 1
my.als <- exp(seq(log(10), log(300), length.out = 1000))
theta <- 4 * Ne * u
my.deltal <- 0.05
## param.table <- read.table("bgsStuff/BGS_twoEffectPrevInsensitivityParamsTable.txt",header=TRUE)##[c(1,20),]
out <- numeric()

## work computer
##these.ones <- c(16,29,34)

## laptop
## these.ones <- c(12,29,34)
my.bs <- c(0.65, 0.75, 0.76)
my.Ve <- c( 531.58, 321.82, 313.79 )/3

## param.table <- param.table[these.ones,]

recover.flag <- F
fixed.deltal.solns <- list()
output <- list()
## i=1
for(i in 1:3){
    my.gs <- 1 - seq(2/L,0.15,length.out=10000)
    fixed.deltal.solns[[i]] <- list()
    last.al <- NULL
    last.tstar <- NULL
    last.norm.al <- NULL
    last.norm.tstar <- NULL
    for(j in 1:length(my.gs)){
        fixed.deltal.solns[[i]][[j]] <- solveTwoEffect2DFixedDeltal(
            bt = my.bs[i],
            bs = my.bs[i],
            Ne = 5000,
            as = 1,
            deltal = my.deltal,
            L = L/my.gs[j],
            gs = my.gs[j],
            ## h2 = 0.5,
            Ve = my.Ve[i],
            u = u,
            C = 1/2,
            Bval=1,
            init.al=last.al,
            init.tstar=last.tstar,
            norm.al=last.al,
            norm.init.tstar=last.norm.tstar
        )
        if(fixed.deltal.solns[[i]][[j]]["code"]==3) break
        last.al <- fixed.deltal.solns[[i]][[j]]["al"]
        last.tstar <- fixed.deltal.solns[[i]][[j]]["tstar"]
        last.norm.al <- fixed.deltal.solns[[i]][[j]]["norm.al"]
        last.norm.tstar <- fixed.deltal.solns[[i]][[j]]["norm.tstar"]
    }
    output[[i]] <- cbind(do.call(rbind,fixed.deltal.solns[[i]]),L=L/my.gs[1:j],Ll=L/my.gs[1:j]-L,gs=my.gs[1:j])
    ##if(i==1) break
}
save(fixed.deltal.solns,file='tmp.fixed.deltal.solns.Robj')
save(output,file='fixed.deltal.output.Robj')


png('figures/largeEffectGrid.png', width = 15, height = 10 , units = 'in', res = 300)
par(mfrow=c(2,3))
i=1
plot(
    output[[i]][,'h2s'],
    output[[i]][,'prev'],
    type = 'l',
    xlim = c(0,1),
    ylim = c(0, max(c(output[[i]][,'prev'],output[[i]][,'naive.norm.prev']))) ,
    xlab = 'Heritability due to small effects' ,
    ylab = 'Prevalence'
)
lines(
    output[[i]][,'h2s'],
    output[[i]][,'naive.norm.prev'],
    col='red'
)
legend(
    'topright',
    col = c('black','red') ,
    lty = 1,
    legend = c('Poisson convolution model','Normal model') 
)
plot(
    output[[i]][,'Ll'],
    output[[i]][,'prev'],
    type = 'l',
    ylim = c(0, max(c(output[[i]][,'prev'],output[[i]][,'naive.norm.prev']))),
    xlab = 'Number of large effect sites' ,
    ylab = 'Prevalence'
)
lines(
    output[[i]][,'Ll'],
    output[[i]][,'naive.norm.prev'],
    col='red'
)
plot(
    output[[i]][,'Ll'],
    output[[i]][,'pgol'],
    type = 'l' ,
    xlab = 'Number of large effect sites' ,
    ylab = 'Proportion of genetic variance in risk that is due to large effects'
)
plot(
    output[[i]][,'pgol'],
    output[[i]][,'h2s'],
    type = 'l',
    xlim = c(0,1),
    xlab = 'Proportion of genetic variance in risk that is due to large effects',
    ylab = 'Heritability of liability due to small effects'
)
plot(
    output[[i]][,'Ll'],
    output[[i]][,'h2s'],
    type = 'l' ,
    xlab = 'Number of large effect sites' ,
    ylab = 'Heritability of liability due to small effects' ,
    )
plot(
    output[[i]][,'pgol'],
    output[[i]][,'prev'],
    type = 'l',
    xlim = c(0,1),
    ylim = c(0, max(c(output[[i]][,'prev'],output[[i]][,'naive.norm.prev']))) ,
    xlab = 'Proportion of genetic variance in risk that is due to large effects' ,
    ylab = 'Prevalence'
)
lines(
    output[[i]][,'pgol'],
    output[[i]][,'naive.norm.prev'],
    col='red'
)
dev.off()

if(FALSE){

    plot(
        output[[i]][,'h2s'],
        1-output[[i]][,'pgol'],
        type = 'l',
        xlim = c(0,1),
        ylim = c(0, 1) ,
        xlab = 'Liability scale heritability due to small effects' ,
        ylab = 'Risk scale heritability due to small effects'
    )

    plot(
        1-my.gs ,
        output[[i]][,'prev'],
        type = 'l',
        xlim = c(0,0.1),
        ylim = c(0, max(c(output[[i]][,'prev'],output[[i]][,'naive.norm.prev']))),
        xlab = 'Fraction of sites with large effects' ,
        ylab = "Prevalence"
    )
    lines(
        1-my.gs ,
        output[[i]][,'naive.norm.prev'],
        col='red'
    )
    lines(
        1-my.gs ,
        output[[i]][,'norm.prev'],
        col='blue'
    )


    plot(
        1-my.gs ,
        output[[i]][,'al'],
        type = 'l',
        xlim = c(0,0.15)
    )

}
