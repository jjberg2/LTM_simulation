rm(list=ls())
setwd('~/Dropbox/msdbPaper')
source('scripts/newSolveTwoEffect.R')
recover.flag <- FALSE


## params
theory.L <- c(1.5e6, 5e6, 1.5e7, 5e7, 1.5e8)
theory.u <- 1e-8
theory.Ne <- 20000
theory.al <- 100

theory.cost <- 0.1
bt <- 0.5
pt <-  bt/2
h2 <- 1/2
as <- 1 

sim.Ne <- 1000
Ne.rt.ratio <- sqrt (sim.Ne /theory.Ne )
sim.cost <- theory.cost / Ne.rt.ratio
sim.al <- theory.al * Ne.rt.ratio
sim.L <- theory.L * 0.01
sim.u <- theory.u / 0.01

al <- sim.al
cost <- sim.cost
Ne <- sim.Ne
L <- sim.L
u <- sim.u

setwd('~/Documents/academics/LTM_simulation')
my.gs <- get(load(file = 'paramFiles/twoEffect.gs.Robj'))
setwd('~/Dropbox/msdbPaper')

gs.params <- unlist(my.gs)
L.params.list <- mapply(
  function (GS,MYL) 
    rep(MYL,times=length(GS)) , 
  GS = my.gs , 
  MYL = L )
L.params <- unlist(
  L.params.list
)
Ls.params.list <- mapply(
  function (GS,MYL) 
    MYL*GS , 
  GS = my.gs , 
  MYL = L )
Ls.params <- unlist(
  Ls.params.list
)
Ll.params.list <- mapply(
  function (GS,MYL) 
    MYL*(1-GS) , 
  GS = my.gs , 
  MYL = L )
Ll.params <- unlist(
  Ll.params.list
)

a.mean.list <- mapply(
  function(LS,LL,L){
    (LS * as + LL * al) / L
  } ,
  LS = Ls.params.list ,
  LL = Ll.params.list ,
  L = L.params.list
)
a.mean <- unlist(a.mean.list)
bs <- (a.mean * bt - (1-gs.params)*al) / ( gs.params * as )
rhos.params <- bs / 2

thr.params <- unlist(
  mapply(
    function (GS,MYL,AMEAN) 
      2*MYL*pt*AMEAN , 
    GS = my.gs , 
    MYL = L ,
    AMEAN = a.mean.list
  )
)



soln <- list()
for (j in seq_along(L)) {
  {
    init.ft <- (4*Ne*cost*as)^(-1) * log ((1+bt)/(1-bt))
    init.Vas <- 2*L[j]*my.gs[[j]]*u*as*bt / (cost*init.ft)
    init.Vtot <- init.Vas/h2
    init.tstar <- dnorminv(init.ft,s=sqrt(init.Vtot))
    init.Ft <- pnorm(init.tstar,sd=sqrt(init.Vtot),lower.tail=FALSE)
    init.Fta <- pnorm(init.tstar - al,sd=sqrt(init.Vtot),lower.tail=FALSE)
    init.dl <- init.Fta - init.Ft
  }
  soln[[j]] <- list()
  for (i in 1:length(my.gs[[j]])) {
    soln[[j]][[i]] <- solveTwoEffect1D(
      bt = bt,
      Ne = Ne,
      as = as,
      al = al,
      L = L[j],
      gs = my.gs[[j]][i],
      h2 = h2 ,
      Ve = NULL ,
      u = u,
      C = cost,
      Bval = 1,
      init.dl = ifelse(i==1,init.dl,prev.dl),
      init.tstar = ifelse(i==1,init.tstar,prev.tstar)
    )
    if(is.na(soln[[j]][[i]]['dl'])) break
    prev.dl <- soln[[j]][[i]]['dl']
    prev.tstar <- soln[[j]][[i]]['tstar']
  }
}


Ve.params <- unlist(sapply(soln,function(X) sapply(X,function(Y)Y['Ve'])))

my.params <- data.frame(
  'rhot' = (1-bt)/2 , 
  'rhos' = rhos.params ,
  'thr' = thr.params , 
  'C'= cost ,
  'al' = al ,
  'Ne' = Ne ,
  'Ve' = Ve.params , 
  'as' = as ,
  'Ls' = floor(Ls.params) ,
  'Ll' = ceiling(Ll.params) ,
  'u' = u
  )

2*my.params$Ls * my.params$rhos


setwd('~/Documents/academics/LTM_simulation')
write.table(my.params, "paramFiles/twoEffectPrevParamTableN1000.txt", col.names=T, row.names=F, quote=F)
