source('scripts/solveTwoEffect.R')
C <- c(seq(0.05,0.97,by=0.06),1)

my.table <- get(load("twoEffectCostInsensitivityPreliminaryParamsTable.Rdata"))

targets <- c(0.1,0.3,0.5,0.7,0.9)
c.table <- list()
for(j in seq_along(my.table)){
    keep <- numeric()
    for(i in seq_along(targets)){
        keep[i] <- which.min(abs(targets[i] - my.table[[j]]$deltal))
    }
    these.ones <- unlist(keep)
    a.table <- my.table[[j]][these.ones,]
    b.table <- a.table[rep(1:length(targets),each=length(C)),]
    b.table$C <- rep(C,times=length(targets))
    c.table[[j]] <- b.table[,!colnames(b.table) %in% c('h2s','h2l','prev','h2os','h2ol','deltas','deltal')]
}



out.table <- do.call(rbind, c.table)
save(c.table,file="twoEffectCostInsensitivityParamsTable.Rdata")
write.table(out.table, "twoEffectCostInsensitivityParamsTable.txt", col.names=T, row.names=F, quote=F)

















if(FALSE){

recover.flag <- TRUE

## params
L.init <- 1e7
Lmeana <- 1e7
Ne <- 5e2
u <- 1e-8
as <- 1

r2n <- 1/2
theta <- 4*Ne*u



my.bt <- 0.4
al <- 15

last.Ll <- numeric()
Ls <-  9150253
Ll <-  59501
L <- Ls + Ll
meana <- (Ls + al*Ll)/L
Lmeana <- L*meana
last.bs <- numeric()
bs <- 0.38
rhos <- 1/2 - bs/2
rhot <- (rhos * 1*Ls + al*Ll)/Lmeana
last.Ls <- numeric()
last.Ls[1] <- L.init - last.Ll









my.cols <- c('as','Ls','bs','rhos','al','Ll','bl','rhol','bt','rhot','Ne','u','C','Ve','h2s','h2l','prev')
output <- data.frame(matrix(ncol = length(my.cols), nrow = 0))
colnames(output) <- my.cols

solns <- list()
for( i in seq_along(C)){

    solns[[i]] <- solveTwoEffect(
        bs=bs,
        bt=my.bt,
        Ne=Ne,
        Ls=Ls,
        Ll=Ll,
        as=as,
        al=al,
        L.init=L.init,
        r2n=r2n,
        Lmeana=Lmeana,
        Ve=NULL,
        u=u,
        C=C[i],
        LL.soln=FALSE,
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




#a=read.table("params_test_twoeffect_Ls1e6.txt")
}
