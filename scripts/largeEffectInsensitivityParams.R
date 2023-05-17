my.thr <- c(1,2,3,4,5,10,25) - 0.001
my.Ne <- 5000
L <- 2500
h2 <- 0.99999
mu <- 1e-6
my.theta <- 4*my.Ne*mu
sim.costs <- seq(0.01,1,by=0.03)
target.prev <- 0.01
anchor.cost <- 0.1
