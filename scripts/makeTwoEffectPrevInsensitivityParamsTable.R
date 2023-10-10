source('scripts/solveTwoEffect.R')


recover.flag <- FALSE

## params
L.init <- 1e7
Lmeana <- 1e7
Ne <- 3e3
u <- 1e-8
as <- 1 
C <- 1
Ve <- 385
##r2n <- 1/2
theta <- 4*Ne*u



my.bt <- 0.8
my.als <- exp(seq(log(12),log(130),length.out=1000))

last.Ll <- numeric()
last.Ll[1] <- 90000
last.bs <- numeric()
last.bs[1] <- 0.42
last.Ls <- numeric()
last.Ls[1] <- L.init - last.Ll

# my.cols <- c('as','Ls','bs','rhos','al','Ll','bl','rhol','bt','rhot','Ne','u','C','Ve','h2s','h2l','prev')
# output <- data.frame(matrix(ncol = length(my.cols), nrow = 0))
# colnames(output) <- my.cols

solns <- list()
tmp.output <- list()
for( i in seq_along(my.als)){
    
    if(i == 1){
      last.tstar = NULL
    } else{
      last.tstar = tmp.output[[i-1]]['tstar']
    }
  
    solns[[i]] <- solveTwoEffect(
        bs=last.bs[i],
        bt=my.bt,
        Ne=Ne,
        Ls=last.Ls[i],
        Ll=last.Ll[i],
        as=as,
        al=my.als[i],
        L.init=L.init,
        last.tstar=NULL,
        r2n=NULL,
        Lmeana=Lmeana,
        Ve=Ve,
        u=u,
        C=C,
        LL.soln=TRUE,
        var.ratio=4,
        equalize.observed.vars=TRUE
    )

    tmp.output[[i]] <- makeOutput(solns[[i]]) 

    last.Ls[i+1] <- tmp.output[[i]]['Ls']
    last.Ll[i+1] <- tmp.output[[i]]['Ll']
    last.bs[i+1] <- tmp.output[[i]]['bs']
    if(i %% 10 == 0) print(i)
}
output <- as.data.frame(do.call(rbind,tmp.output))


keep1 = which(output$deltal < 100/(4*Ne*C))
targets <- seq(100/(4*Ne*C),tail(output$deltal,1),length.out=100)
keep <- list()
for(i in seq_along(targets)){
    keep[[i]] <- which.min(abs(targets[i] - output$deltal))
}
these.ones <- sort(unique(c(unlist(keep),keep1)))


norm.ft <- 1/(4*N*C)*log((1+my.bt/(1-my.bt)))
norm.Vg <- 4*N*L.init*u*my.bt / log((1+my.bt/(1-my.bt)))
norm.Vt <- norm.Vg + Ve
h2 <- norm.Vg / norm.Vt
norm.astd <- 1 / sqrt ( norm.Vt )
norm.dens <- 2*L.init*u*my.bt*norm.astd / (h2*C) 
1-pnorm(dnorminv(norm.dens))


written.output <- output[these.ones,]



write.table(written.output, "twoEffectPrevInsensitivityParamsTable.txt", col.names=T, row.names=F, quote=F)
#a=read.table("params_test_twoeffect_Ls1e6.txt")
