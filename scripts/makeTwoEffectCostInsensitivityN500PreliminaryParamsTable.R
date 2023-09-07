source('scripts/solveTwoEffect.R')

recover.flag <- FALSE

## params
L.init <- 1e7
Lmeana <- 1e7
Ne <- 5e2
u <- 1e-8
as <- 1
C <- 0.4
r2n <- 1/2
theta <- 4*Ne*u



my.bt <- c(0.15,0.5,0.9)

my.als <- exp(seq(log(12),log(100),length.out=1000))
written.output <- list()

last.Ll <- lapply(seq_along(my.bt),function(X) numeric())
last.bs <- lapply(seq_along(my.bt),function(X) numeric())
last.Ls <- lapply(seq_along(my.bt),function(X) numeric())
last.Ll[[1]][1] <- 26230
last.Ll[[2]][1] <- 55419
last.Ll[[3]][1] <- 136000

last.bs[[1]][1] <- 0.1223
last.bs[[2]][1] <- 0.46
last.bs[[3]][1] <- 0.88

last.Ls[[1]][1] <- L.init - last.Ll[[1]][1]
last.Ls[[2]][1] <- L.init - last.Ll[[2]][1]
last.Ls[[3]][1] <- L.init - last.Ll[[3]][1]

for( j in 1:length(my.bt)){

    

    
    my.cols <- c('as','Ls','bs','rhos','al','Ll','bl','rhol','bt','rhot','Ne','u','C','Ve','h2s','h2l','prev')
    output <- data.frame(matrix(ncol = length(my.cols), nrow = 0))
    colnames(output) <- my.cols

    solns <- list()
    for( i in seq_along(my.als)){

        solns[[i]] <- solveTwoEffect(
            bs=last.bs[[j]][i],
            bt=my.bt[j],
            Ne=Ne,
            Ls=last.Ls[[j]][i],
            Ll=last.Ll[[j]][i],
            as=as,
            al=my.als[i],
            L.init=L.init,
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
        output[i,'bigU'] <- output[i,'maxG']*u*my.bt[j]

        
        output[i,'thr'] <- 2*output[i,'rhos']*output[i,'Ls']*output[i,'as'] + 1


        last.Ls[[j]][i+1] <- output[i,'Ls']
        last.Ll[[j]][i+1] <- output[i,'Ll']
        last.bs[[j]][i+1] <- output[i,'bs']
        if(i %% 10 == 0) print(i)
    }


    targets <- seq(head(output$deltal,1),tail(output$deltal,1),length.out=100)
    keep <- list()
    for(i in seq_along(targets)){
        keep[[i]] <- which.min(abs(targets[i] - output$deltal))
    }
    these.ones <- unlist(keep)

    my.gamma <- log((1+my.bt[j])/(1-my.bt[j]))
    env.sd <- sqrt(theta*Lmeana*my.bt[j]/my.gamma*(1-r2n)/r2n)
    U <- 2*Lmeana*u
    BigGamma <- 4*Ne*C
    norm.prev <- eqNormPrev(U*my.bt[j],BigGamma,my.gamma,r2n)

    written.output[[j]] <- output[these.ones,]

}

save(written.output, file="twoEffectCostInsensitivityPreliminaryParamsTable.Rdata")
#a=read.table("params_test_twoeffect_Ls1e6.txt")
