source('scripts/solveSingleEffect.R')


args = commandArgs(trailingOnly=TRUE)
liaMutRate = as.numeric(args[1])
THR = as.numeric(args[2])
fitCost = as.numeric(args[3])
popSize = as.numeric(args[4])
liabilitySize = as.numeric(args[5])
target_h2 = as.numeric(args[6])
theta = 4 * popSize * liaMutRate

if (THR <1000){
   temp = SolvePoissonGivenThr(theta,THR, fitCost,popSize, liabilitySize, target_h2,verbose.output=T)
   poissonGamma = temp$my.gamma
   if (poissonGamma <50){
      b = (exp(poissonGamma)-1)/(exp(poissonGamma)+1)
   }else {b=1}
   envSD = sqrt(2.0 * liabilitySize * liaMutRate/temp$my.sel * b * (1-target_h2)/target_h2)
}else{
   rho= THR/(2*liabilitySize)
   gamma = log((1-rho)/rho)
   envSD = sqrt(2.0 * liabilitySize * liaMutRate/(gamma/(4*popSize)) * (exp(gamma)-1)/(exp(gamma)+1) * (1-target_h2)/target_h2)
}
print(envSD)
