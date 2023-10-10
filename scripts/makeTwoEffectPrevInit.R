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
last.bs[1] <- 0.32
last.Ls <- numeric()
last.Ls[1] <- L.init - last.Ll

my.cols <- c('as','Ls','bs','rhos','al','Ll','bl','rhol','bt','rhot','Ne','u','C','Ve','h2s','h2l','prev')
output <- data.frame(matrix(ncol = length(my.cols), nrow = 0))