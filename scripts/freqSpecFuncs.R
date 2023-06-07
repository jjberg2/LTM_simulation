#dnorminv<-function(y) sqrt(-2*log(sqrt(2*pi)*y))

GetPi <- function ( alpha , R , x = 0 , diploid = TRUE ) {
    ## recover()
    if ( diploid ) {
        Raa <- pnorm ( qnorm ( 1-R ) - 2*alpha*(1-x) , lower.tail = FALSE )
        RAa <- pnorm ( qnorm ( 1-R ) - 2*alpha*(1/2-x) , lower.tail = FALSE )
        RAA <- pnorm ( qnorm ( 1-R ) + 2*alpha*x , lower.tail = FALSE )

        Ra <- x * Raa + ( 1 - x ) * RAa
        RA <- x * RAa + ( 1 - x ) * RAA
    } else {
        Ra <- pnorm ( qnorm ( 1-R ) - alpha*(1-x) , lower.tail = FALSE )
        RA <- pnorm ( qnorm ( 1-R ) + alpha*x , lower.tail = FALSE )
    }
    Ra - RA
}
GetAlpha <- function ( pi , P ) {
    qnorm ( P , lower.tail = FALSE ) - qnorm ( P + pi , lower.tail = FALSE )
}
FreqSpec <- function ( x , gamma , theta = 1000 ) {
    ## if gamma is negative, the mutation is deleterious, if positive, it is beneficial
    one <- theta / ( x * ( 1 - x ) )
    two <- ( 1 - exp ( -2*gamma*(1-x))) /( 1 - exp ( -2*gamma ) )
    one * two
}
meanDerFreq <- function(fourNS,theta,Ne){
    theta <- 10
    num <- integrate(
        function(X) X*FreqSpec(X,fourNS,theta) ,
        lower=1/(4*Ne),
        upper=1-1/(4*Ne)
    )$value
    denom <- integrate(
        function(X) FreqSpec(X,fourNS,theta) ,
        lower=1/(4*Ne),
        upper=1-1/(4*Ne)
    )$value
    num/denom
}
SiteSpec <- function ( x , gamma , theta , normalize=FALSE ) {

    ## if gamma is negative, the mutation is deleterious, if positive, it is beneficial
    one <- ( x * ( 1 - x ) )^(theta-1)
    two <- exp ( -2*gamma*x )
    if ( !normalize ) {
        one * two
    } else if ( normalize ) {
        denom <- integrate(
            function(x) SiteSpec(x,gamma,theta),
            lower=0,
            upper=1
        )$value
        one*two/denom
    }
}
SiteSpecInt <- function(gamma4NS,theta,x1,x2){
    myInt <- integrate(
        function(x) SiteSpec(x,gamma4NS,theta),
        lower=x1,
        upper=x2
    )$value
    return(myInt)
}
VarSpecInt <- function(gamma4NS,theta,x1,x2){
    myInt <- integrate(
        function(x) x*(1-x)*SiteSpec(x,gamma4NS,theta),
        lower=x1,
        upper=x2
    )$value
    return(myInt)
}
SiteSpecFullInt <- function(gamma4NS,theta){
    gamma4NS <- abs(gamma4NS)
    if ( gamma4NS == 0 ) {
        gamma(theta)^2/gamma(2*theta)
    } else {
        exp(-gamma4NS/2)*sqrt(pi)*gamma4NS^(0.5-theta)*besselI(gamma4NS/2,theta-0.5)*gamma(theta)
    }
}
IntSegRange <- function(gamma,theta,epsilon,my.range=NULL){
    if(is.null(my.range)){
        my.range <- c(epsilon,1-epsilon)
    }

    num <- mapply(
        FUN=function(THETA,GAMMA){
            integrate(
                f = function ( x )
                    SiteSpec(x,GAMMA,THETA),
                lower=my.range[1],
                upper=my.range[2]
            )$value
        },
        THETA=theta,
        GAMMA=gamma
    )
    num

}
InfSiteSpec <- function(x,gamma,theta){

    one <- theta/(x*(1-x))
    two <- exp(2*gamma*x)/(1+exp(2*gamma))
    one*two

}
SegSites <- function ( gamma , Ne , theta = 1000 ) {
    ## if gamma is negative, the mutation is deleterious, if positive, it is beneficial
    integrate (
        f = function ( x ) FreqSpec ( x , gamma , theta  ) ,
        lower = 1 / ( 2 * Ne ) ,
        upper = 1 - 1 / ( 2 * Ne )
    )$value
}
VarSpec <- function ( x, gamma , theta = 1000 , approx=FALSE ) {

    if(gamma < -500 & !approx){
        theta*( 1 - exp ( -gamma*(1-x)))/( 1 - exp ( -gamma ) )
    } else {
        theta*exp(-gamma*x)
    }
}
ChiSqNonCentral <- function(N,Q,R,x,pii){
    R.A <- R-pii*x
    R.a <- R+pii*(1-x)
    x.case <- x*R.a/R
    x.control <- x*(1-R.a)/(1-R)
    x.hat <- x.case*Q+x.control*(1-Q)
    chisq <- N*Q*(1-Q)*(x.case-x.control)^2/(x.hat*(1-x.hat))
    chisq
}
FiniteSiteVarSpec <- function ( x, gamma , theta = 1000 , approx=FALSE ) {

##    if(gamma < -500 & !approx){
    num <- exp(-gamma*x)/(x*(1-x))^(theta)
    denom <- exp(-gamma/2)*sqrt(pi)*gamma^(0.5-theta)*besselI(gamma/2,theta-0.5)*gamma(theta)
    num/denom
##    } else {
##        theta*exp(-gamma*x)
##    }
}
MutBias <- function ( gamma ) {
    exp ( -gamma ) / ( exp ( gamma ) + exp ( - gamma ) )
}
EqFreqSpec <- function ( x , gamma , theta = 1000 ) {
    MutBias ( -gamma ) * FreqSpec ( x , gamma , theta )
}
EqTotVar <- function ( gamma , theta = 1000 , Ne = 1e4 , lower = NULL , upper = NULL ) {
    if ( is.null ( lower ) )
        lower <- 1 / ( 2 * Ne )
    if ( is.null ( upper ) )
        upper <- 1 - 1 / ( 2 * Ne )
    integrate (
        f = function ( x ) MutBias ( -gamma ) * gamma^2*VarSpec ( x , gamma , theta ) ,
        lower = lower,
        upper = upper
    )
}
EqVarSpec <- function ( x , gamma , theta = 1000 ) {
    MutBias ( -gamma ) * VarSpec ( x , gamma , theta )
}
EqFreqSpecBoth <- function ( x , gamma , theta = 1000 ) {
    FreqSpec ( x , gamma , theta ) +
        EqFreqSpec ( x , -gamma , theta )
}
EqVarSpecBoth <- function ( x , gamma , theta = 1000 ) {
    EqVarSpec ( x , gamma , theta ) +
        EqVarSpec ( x , -gamma , theta )
}
MapToRiskScale <- function ( alpha , P , x = 0 ) {

    recover()
    t <- qnorm ( 1 - P )

    t1 <- t - 2*alpha * ( 1 - x )
    t2 <- t - 2*alpha * ( 1/2 - x )
    t3 <- t + 2*alpha * x

    A <- pnorm ( t1 , lower.tail = FALSE ) - pnorm ( t3 , lower.tail = FALSE )
    D <- A/2 - ( pnorm ( t1 , lower.tail = FALSE ) - pnorm ( t2 , lower.tail = FALSE ) )

    A/2 + D * ( 1 - 2*x )
}
FixProb <- function ( Ne , s , N = Ne ) {
    one <- 1-exp(-4*Ne*s/(2*N))
    two <- (1-exp(-4*Ne*s))
    one / two
}

## library(Rmpfr)
library(expint)
ProbSeg <- function(theta,gamma,Ne,inf.sites=TRUE) {

    if(!inf.sites){
        stopifnot(gamma<1480)
        num <- mapply(
            FUN=function(THETA,GAMMA){
                integrate(
                    f = function ( x )
                        SiteSpec(x,GAMMA,THETA),
                    lower=1/(4*Ne),
                    upper=1-1/(4*Ne)
                )$value
            },
            THETA=theta,
            GAMMA=gamma
        )
        denom <- exp(-gamma/2)*sqrt(pi)*gamma^(0.5-theta)*besselI(gamma/2,theta-0.5)*gamma(theta)
        return(num/denom)
    } else {
        num <- expint_Ei(gamma*(1/(2*Ne)-1))-expint_Ei(-gamma/(2*Ne))
        return(theta*num)
    }

}

SegFreq <- function(myTheta,myGamma,Ne){

    recover()
    mu <- myTheta/(4*Ne)
    my.s <- myGamma/(4*Ne)
    det.freq <- mu/my.s
    seg.probs <- ProbSeg(myTheta,myGamma,Ne)
    denom <- IntSegRange(myGamma,myTheta,1/(4*Ne))
    num <- mapply(
        FUN = function(THETA,GAMMA){
            integrate(
                f = function(X)
                    X*SiteSpec(X,GAMMA,THETA),
                lower=1/(4*Ne),
                upper=1-1/(4*Ne)
            )$value
        },
        THETA=myTheta,
        GAMMA=new.gammas
    )
    num/denom
}
