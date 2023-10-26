rm(list = ls())
library('MetBrewer')
source('scripts/solveTwoEffect.R')
recover.flag <- F
## main single effect params
L <- 1e7
Ne <- 1e3
u <- 1.3e-8
cost <- 1 / 2 
my.bt <- seq(0.3, 0.9, length.out = 2)  ##c(0.6,0.6,0.8,0.8)
h2 <- c(1 / 6, 5 / 6)
## ratio of large effect observed var to small effect
var.ratio <- c(1 / 3, 1, 3)
as <- 1
my.als <- exp(seq(log(10), log(300), length.out = 1000))
theta <- 4 * Ne * u
tmp.output <- list()
nulls <- list()
output <- list()
solns <- list()
output.for.sims <- data.frame()
for (l in seq_along(my.bt)) {
  output[[l]] <- list()
  solns[[l]] <- list()
  for (k in seq_along(h2)) {
    output[[l]][[k]] <- list()
    solns[[l]][[k]] <- list()
    for (j in seq_along(var.ratio)) {
      last.gs <- numeric()
      last.gs[1] <- 0.99
      last.bs <- numeric()
      last.bs[1] <- my.bt[l] - 0.05
      solns[[l]][[k]][[j]] <- list()
      for (i in seq_along(my.als)) {
        if(l==2 && k==2 && j==1 && i==1) recover.flag <- F
        solns[[l]][[k]][[j]][[i]] <- solveTwoEffect(
          bs = last.bs[i],
          bt = my.bt[l],
          Ne = Ne,
          as = as,
          al = my.als[i],
          L = L,
          gs = last.gs[i],
          last.tstar = NULL,
          h2 = h2[k],
          u = u,
          C = cost,
          LL.soln = TRUE,
          var.ratio = var.ratio[j],
          equalize.observed.vars = TRUE
        )
        tmp.output[[i]] <- makeOutput(solns[[l]][[k]][[j]][[i]])
        last.gs[i + 1]  <- tmp.output[[i]]['gs']
        last.bs[i + 1] <- tmp.output[[i]]['bs']
        if (i %% 10 == 0)
          print(paste(l, k, j, i, sep = ':'))
      }
      output[[l]][[k]][[j]] <- as.data.frame(do.call(rbind, tmp.output))
      targets <- seq(0.01,0.99,length.out=99)
      keep <- list()
      for(i in seq_along(targets)){
        keep[[i]] <- which.min(abs(targets[i] - output[[l]][[k]][[j]]$deltal))
      }
      these.ones <- sort(unlist(keep))
      output.for.sims <- rbind(output.for.sims,output[[l]][[k]][[j]][these.ones,])
    }
  }
}


save(solns, file = 'solutions/twoEffectPrevalenceSolutions_solutions.Robj')
save(output, file = 'solutions/twoEffectPrevalenceSolutions_output.Robj')
save(output.for.sims, file = 'solutions/twoEffectPrevalenceSolutions_outputForSims.Robj')

