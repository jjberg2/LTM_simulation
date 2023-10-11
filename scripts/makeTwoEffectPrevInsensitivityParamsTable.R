source('scripts/solveTwoEffect.R')
recover.flag <- FALSE


## main single effect params
L <- 1e7
Ne <- 1e3
u <- 1e-8
C <- 1
my.bt <- 0.8
h2 <- 1/2

## ratio of large effect observed var to small effect
var.ratio <- 2

as <- 1 
my.als <- exp(seq(log(8),log(200),length.out=1000))
theta <- 4*Ne*u
last.Ll <- numeric()
last.Ll[1] <- 90000
last.bs <- numeric()
last.bs[1] <- 0.7

solns <- list()
tmp.output <- list()
for( i in seq_along(my.als)){
    ## if(i == 1){
    ##   last.tstar = NULL
    ## } else{
    ##   last.tstar = tmp.output[[i-1]]['tstar']
    ## }
    solns[[i]] <- solveTwoEffect(
        bs=last.bs[i],
        bt=my.bt,
        Ne=Ne,
        as=as,
        al=my.als[i],
#        Ls=last.Ls[i],
        L=L,
        Ll=last.Ll[i],
        last.tstar=NULL,
        h2=h2,
#        Lmeana=Lmeana,
        u=u,
        C=C,
        LL.soln=TRUE,
        var.ratio=var.ratio,
        equalize.observed.vars=TRUE
    )
    tmp.output[[i]] <- makeOutput(solns[[i]]) 
    last.Ll[i+1] <- tmp.output[[i]]['Ll']
    last.bs[i+1] <- tmp.output[[i]]['bs']
    if(i %% 10 == 0) print(i)
}
output <- as.data.frame(do.call(rbind,tmp.output))
my.max <- max(which(output$deltal < 100/(4*Ne*C)))
keep1 <- round(seq(1,my.max,length.out=20))
targets <- seq(100/(4*Ne*C),tail(output$deltal,1),length.out=100)
keep <- list()
for(i in seq_along(targets)){
    keep[[i]] <- which.min(abs(targets[i] - output$deltal))
}
these.ones <- sort(unique(c(unlist(keep),keep1)))
norm.ft <- 1/(4*Ne*C)*log((1+my.bt/(1-my.bt)))
norm.Vg <- 8*Ne*L*u*my.bt / log((1+my.bt/(1-my.bt)))
norm.Vt <- norm.Vg/h2
norm.astd <- 1 / sqrt ( norm.Vt )
norm.dens <- 2*L*u*my.bt*norm.astd / (h2*C) 
norm.prev <- 1-pnorm(dnorminv(norm.dens))

written.output <- output[these.ones,]

write.table(written.output, "twoEffectPrevInsensitivityParamsTable.txt", col.names=T, row.names=F, quote=F)
#a=read.table("params_test_twoeffect_Ls1e6.txt")
