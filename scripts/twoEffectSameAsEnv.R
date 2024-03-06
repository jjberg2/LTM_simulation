rm(list = ls())
library('MetBrewer')
source('scripts/solveTwoEffect.R')
recover.flag <- F
## main single effect params
L <- 1e5
Ne <- 5e3
u <- 1e-7
cost <- 1 / 2 
my.bt <- seq(0.3, 0.9, length.out = 2)  ##c(0.6,0.6,0.8,0.8)
h2 <- c(1 / 6, 5 / 6)
## ratio of large effect observed var to small effect
var.ratio <- c(1 / 3, 1, 3)
as <- 1
my.als <- exp(seq(log(10), log(300), length.out = 1000))
theta <- 4 * Ne * u


Bvalue <- 0.82

param.table <- read.table("bgsStuff/BGS_twoEffectPrevInsensitivityParamsTable.txt",header=TRUE)##[c(1,20),]
out <- numeric()

these.ones <- c(16,29,34)

param.table <- param.table[these.ones,]
 

solns <- list()
output <- list()
for(i in 1:nrow(param.table)){
    init.Ll <- param.table$Ll[i]
    my.Lls <- seq(param.table$Ll[i]/30,param.table$Ll[i]*30 , length.out = 100)
    Ll.diffs <- my.Lls - init.Ll
    solns[[i]] <- list()
    for(j in 1:length(my.Lls)){
        solns[[i]][[j]] <- solveTwoEffect2D(
            bt = param.table$bs[i],
            bs = param.table$bs[i],
            Ne = param.table$Ne[i],
            as = param.table$as[i],
            al = param.table$al[i],
            L = L+Ll.diffs[j],
            gs = param.table$Ls[i] / ( param.table$Ls[i] + param.table$Ll[i] + Ll.diffs[j] ),
            ## h2 = 0.5,
            Ve = param.table$Ve[i],
            u = u,
            C = param.table$C[i],
            Bval=1
        )
    }
    output[[i]] <- do.call(rbind,solns[[i]])
}


plot(
    output[[1]][,'h2s.noBGS'],
    output[[1]][,'prev.noBGS'],
    type = 'l',
    xlim = c(0,0.5),
    ylim = c(0, max(c(output[[1]][,'prev.noBGS'],output[[1]][,'naive.norm.prev.noBGS'])))
)
lines(
    output[[1]][,'h2s.noBGS'],
    output[[1]][,'naive.norm.prev.noBGS'],
    col='red'
)








for(i in 1:nrow(param.table)){

    soln <- solveTwoEffect(
        bt = param.table$bt[i],
        bs = param.table$bs[i],
        Ne = param.table$Ne[i],
        as = param.table$as[i],
        al = param.table$al[i],
        L = L,
        gs = param.table$Ls[i] / ( param.table$Ls[i] + param.table$Ll[i] ),
        h2 = 0.5,
##        Ve = param.table$Ve[i],
        u = u,
        C = param.table$C[i],
        LL.soln = TRUE,
        var.ratio = 2,
        equalize.observed.vars = TRUE
    )
    my.output <- makeOutput(soln)

    
    
    
    blah <- solveTwoEffect2D(
        bt = my.output['bt'],
        bs = my.output['bs'],
        Ne = my.output['Ne'],
        as = my.output['as'],
        al = my.output['al'],
        L = my.output['Ls'] + my.output['Ll'],
        gs = my.output['Ls'] / ( my.output['Ls'] + my.output['Ll']),
        h2 = 0.5,
        Ve = my.output['Ve'],
        u = u,
        C = my.output['C'],
        Bval=Bvalue
        )

    out[i] <- solveTwoEffect2D(
        bt = param.table$bt[i],
        bs = param.table$bs[i],
        Ne = param.table$Ne[i],
        as = param.table$as[i],
        al = param.table$al[i],
        L = param.table$Ls[i] + param.table$Ll[i],
        gs = param.table$Ls[i] / ( param.table$Ls[i] + param.table$Ll[i] ),
        h2 = 0.5,
        Ve = param.table$Ve[i],
        u = u,
        C = param.table$C[i],
        Bval=Bvalue
        )
    
}
