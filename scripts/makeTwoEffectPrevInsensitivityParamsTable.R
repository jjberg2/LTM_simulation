source('scripts/solveTwoEffect.R')

recover.flag <- FALSE

## params
L.init <- 1e7
Lmeana <- 1e7
Ne <- 5e3
u <- 1e-8
as <- 1
C <- 0.1
r2n <- 1/2
theta <- 4*Ne*u



my.bt <- 0.4
my.als <- exp(seq(log(12),log(300),length.out=1000))

last.Ll <- numeric()
last.Ll[1] <- 90000
last.bs <- numeric()
last.bs[1] <- 0.32
last.Ls <- numeric()
last.Ls[1] <- L.init - last.Ll

my.cols <- c('as','Ls','bs','rhos','al','Ll','bl','rhol','bt','rhot','Ne','u','C','Ve','h2s','h2l','prev')
output <- data.frame(matrix(ncol = length(my.cols), nrow = 0))
colnames(output) <- my.cols

solns <- list()
for( i in seq_along(my.als)){

    solns[[i]] <- solveTwoEffect(
        bs=last.bs[i],
        bt=my.bt,
        Ne=Ne,
        Ls=last.Ls[i],
        Ll=last.Ll[i],
        as=as,
        al=my.als[i],
        r2n=r2n,
        Lmeana=Lmeana,
        Ve=NULL,
        u=u,
        C=C,
        LL.soln=TRUE,
        var.ratio=1,
        equalize.observed.vars=TRUE
    )

    output[i,'as'] <- solns[[i]]$as  ## small effect size
    output[i,'Ls'] <- solns[[i]]$Ls  ## numer of small effect loci
    output[i,'bs'] <- solns[[i]]$bs  ## b for small effect loci
    output[i,'rhos'] <- 1/2-output[i,'bs']/2  ## rho for small effect loci
    output[i,'al'] <- solns[[i]]$al  ## large effect size
    output[i,'Ll'] <- round(solns[[i]]$Ll,0)  ## numer of large effect loci
    output[i,'bl'] <- 1              ## b for large effect loci
    output[i,'rhol'] <- 0            ## rho for large effect loci
    output[i,'bt'] <- solns[[i]]$bt  ## total b
    output[i,'rhot'] <-  1/2 - solns[[i]]$bt/2  ## total rho
    
    output[i,'Ne'] <- solns[[i]]$Ne  ## population size
    output[i,'u'] <- solns[[i]]$u    ## mutation rate
    output[i,'C'] <- solns[[i]]$C    ## cost of disease
    output[i,'Ve'] <- solns[[i]]$Ve  ## environmental variance
    output[i,'h2s'] <- solns[[i]]$h2s  ## small effect h2 on liability scale
    output[i,'h2l'] <- solns[[i]]$h2l  ## large effect h2 on liability scale
    output[i,'h2os'] <- solns[[i]]$h2os  ## small effect h2 on liability scale
    output[i,'h2ol'] <- solns[[i]]$h2ol  ## large effect h2 on liability scale
    output[i,'prev'] <- solns[[i]]$prev ## prevalence
    output[i,'deltas'] <- solns[[i]]$deltas ## risk small effect
    output[i,'deltal'] <- solns[[i]]$deltal ## risk small effect
    output[i,'maxG'] <- 2*output[i,'Ls']* output[i,'as'] + 2*output[i,'Ll']* output[i,'al']
    output[i,'bigU'] <- output[i,'maxG']*u*my.bt
    output[i,'thr'] <- 2*output[i,'rhos']*output[i,'Ls']*output[i,'as'] + 1

    
    last.Ls[i+1] <- output[i,'Ls']
    last.Ll[i+1] <- output[i,'Ll']
    last.bs[i+1] <- output[i,'bs']
    if(i %% 10 == 0) print(i)
}


targets <- seq(head(output$deltal,1),tail(output$deltal,1),length.out=100)
keep <- list()
for(i in seq_along(targets)){
    keep[[i]] <- which.min(abs(targets[i] - output$deltal))
}
these.ones <- unlist(keep)

my.gamma <- log((1+my.bt)/(1-my.bt))
env.sd <- sqrt(theta*Lmeana*my.bt/my.gamma*(1-r2n)/r2n)
U <- 2*Lmeana*u
BigGamma <- 4*Ne*C
norm.prev <- eqNormPrev(U*my.bt,BigGamma,my.gamma,r2n)

written.output <- output[these.ones,]



write.table(written.output, "twoEffectPrevInsensitivityParamsTable.txt", col.names=T, row.names=F, quote=F)
#a=read.table("params_test_twoeffect_Ls1e6.txt")
