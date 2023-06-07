library(numDeriv,quietly=TRUE,warn.conflicts=FALSE)
## library(BB,quietly=TRUE)
library(nleqslv)
source('scripts/freqSpecFuncs.R')
recoversOn <- F
WrightProbFixed <- function(myTheta=4e-4,myGamma=1e1,Ne=1e4,ploidy=2){

    prob.seg <- ProbSeg(myTheta,myGamma,Ne)
    (1-prob.seg)*exp(-myGamma)/(exp(-myGamma)+exp(myGamma))
}
#################################################################################
## CGF and first two derivatives for contribution from genetic component sites ##
#################################################################################
GetLFixed <- function(myTheta,myGamma,Ne,L){
    mapply(
        function(MYL,GAMMA)
            MYL*WrightProbFixed(myTheta=myTheta,myGamma=GAMMA,Ne=Ne,ploidy=ploidy),
        MYL=L,
        GAMMA=myGamma
    )
}
GetWeights <- function(gamma,theta,Ne,ploidy=2){
    k <- 0:ploidy
    tmp.weights <- lapply(gamma,function(GAMMA) sapply(k,function(x) WrightDiploPDF(x,theta,GAMMA,Ne)))
    lapply(tmp.weights,function(TMP)TMP/sum(TMP))
}
binomCGF <- function(t,x,k=2,alpha=1){
    k*log(1+x*(exp(alpha*t)-1))
}
binomCGFprime <- function(t,x,k=2,alpha=1){
    ## t = dummy var
    ## x = frequency
    ## k = ploidy
    ## alpha = effect size
    num <- k*exp(alpha*t)*x
    denom <- 1+x*(exp(alpha*t)-1)
    num/denom
}
binomCGFdprime <- function(t,x,k=2,alpha=1){
    num <- exp(alpha*t)*x
    denom <- 1+x*(exp(alpha*t)-1)
    return(k*alpha^2*(num/denom-(num/denom)^2))
}
binomCGFtprime <- function(t,x,k=2,alpha=1){
    num <- exp(alpha*t)*x
    denom <- 1+x*(exp(alpha*t)-1)
    return(k*alpha^3*(num/denom-(num/denom)^2+(num/denom)^3))
}
WrightGenoCGF <- function(s=0,myTheta=4e-4,myGamma=1e1,alpha=1,Ne=1e4,L=1e4,ploidy=2,epsilon1=1/(4*Ne),epsilon2=1-1/(4*Ne),return.sum=TRUE){
    Lseg <- mapply(
        function(GAMMA,TARGET)
            integrate(
                f=function(X) InfSiteSpec(X,-GAMMA,myTheta*TARGET),
                lower=epsilon1,
                upper=epsilon2
            )$value,
        GAMMA=myGamma,
        TARGET=L
    )
    binomInt <- mapply(
        function(GAMMA,ALPHA){
            integrate(
                f=function(p)
                    binomCGF(s,p,ploidy,ALPHA)*SiteSpec(p,gamma=GAMMA,theta=myTheta),
                lower=epsilon1,
                upper=epsilon2,
                subdivisions=4*Ne
            )$value
        },
        GAMMA=myGamma,
        ALPHA=alpha
    )
    denom <- sapply(myGamma,function(GAMMA)SiteSpecInt(GAMMA,myTheta,epsilon1,epsilon2))
    cgfs <- binomInt/denom
    Lfixed <- GetLFixed(myTheta,myGamma,Ne,L)
    if(return.sum)
        sum(Lseg*cgfs,2*alpha*s*Lfixed)
    else
        list(seg=Lseg*cgfs,fixed=2*alpha*s*Lfixed)
}
WrightGenoCGFprime <- function(s=0,myTheta=4e-4,myGamma=1e1,alpha=1,Ne=1e4,L=1e4,ploidy=2,return.sum=TRUE){
    ## recover()
    Lseg <- mapply(
        function(GAMMA,TARGET)
            integrate(
                f=function(X) InfSiteSpec(X,-GAMMA,myTheta*TARGET),
                lower=1/(4*Ne),
                upper=1-1/(4*Ne),
                subdivisions=4*Ne
            )$value,
        GAMMA=myGamma,
        TARGET=L
    )
    binomInt <- mapply(
        function(GAMMA,ALPHA){
            integrate(
                f=function(p)
                    binomCGFprime(s,p,ploidy,ALPHA)*SiteSpec(p,gamma=GAMMA,theta=myTheta),
                lower=1/(4*Ne),
                upper=1-1/(4*Ne),
                subdivisions=4*Ne
            )$value
        },
        GAMMA=myGamma,
        ALPHA=alpha
    )

    denom <- sapply(myGamma,function(GAMMA)SiteSpecInt(GAMMA,myTheta,1/(4*Ne),1-1/(4*Ne)))
    cgfps <- binomInt/denom
    Lfixed <- GetLFixed(myTheta,myGamma,Ne,L)
    if(return.sum)
        sum(Lseg*cgfps,2*Lfixed*alpha)
    else
        list(seg=Lseg*cgfps,fixed=2*Lfixed*alpha)
}
WrightGenoCGFdoubleprime <- function(s=0,myTheta=4e-4,myGamma=1e1,alpha=1,Ne=1e4,L=1e4,ploidy=2,return.sum=TRUE){

    Lseg <- mapply(
        function(GAMMA,TARGET)
            integrate(
                f=function(X)
                    InfSiteSpec(X,-GAMMA,myTheta*TARGET),
                lower=1/(4*Ne),
                upper=1-1/(4*Ne)
            )$value,
        GAMMA=myGamma,
        TARGET=L
    )

    binomInt <- mapply(
        function(GAMMA,ALPHA){
            integrate(
                f=function(p)
                    binomCGFdprime(s,p,ploidy,ALPHA)*SiteSpec(p,gamma=GAMMA,theta=myTheta),
                lower=1/(4*Ne),
                upper=1-1/(4*Ne),
                subdivisions=4*Ne
            )$value
        },
        GAMMA=myGamma,
        ALPHA=alpha
    )

    denom <- sapply(myGamma,function(GAMMA)SiteSpecInt(GAMMA,myTheta,1/(4*Ne),1-1/(4*Ne)))
    cgfdps <- binomInt/denom
    if(return.sum)
        sum(Lseg*cgfdps)
    else
        c(Lseg*cgfdps)
}
WrightGenoCGFtripleprime <- function(s=0,myTheta=4e-4,myGamma=1e1,alpha=1,Ne=1e4,L=1e4,ploidy=2,return.sum=TRUE){

    Lseg <- mapply(
        function(GAMMA,TARGET)
            integrate(
                f=function(X)
                    InfSiteSpec(X,-GAMMA,theta*TARGET),
                lower=1/(4*Ne),
                upper=1-1/(4*Ne)
            )$value,
        GAMMA=myGamma,
        TARGET=L
    )

    binomInt <- mapply(
        function(GAMMA,ALPHA){
            integrate(
                f=function(p)
                    binomCGFtprime(s,p,ploidy,ALPHA)*SiteSpec(p,gamma=GAMMA,theta=myTheta),
                lower=1/(4*Ne),
                upper=1-1/(4*Ne),
                subdivisions=4*Ne
            )$value
        },
        GAMMA=myGamma,
        ALPHA=alpha
    )

    denom <- sapply(myGamma,function(GAMMA)SiteSpecInt(GAMMA,myTheta,1/(4*Ne),1-1/(4*Ne)))
    cgftps <- binomInt/denom
    if(return.sum)
        sum(Lseg*cgftps)
    else
        c(Lseg*cgftps)
}

## expression for bound
## not currently using cause bound is uselessly large
bns <- function(s,myTheta,myGamma,Ne,L,ploidy=2){
    b <- WrightDiploAbsThr(s,myTheta,myGamma,Ne,ploidy)
    Lseg <- ProbSeg(myTheta,myGamma,Ne)*L
    num <- 2*0.56*Lseg*b
    denom <- WrightGenoCGFdoubleprime(s,myTheta,myGamma,Ne,L,ploidy)^(3/2)
    num/denom
}

#################################################################################
## CGF and first two derivatives for contribution from genetic component sites ##
#################################################################################
EnvVarWarning <- function(h2=NULL,env.var=NULL){
    if(is.null(h2)+is.null(env.var)==2){
        stop('must specify either h2 or env.var')
    } else if(is.null(h2)+is.null(env.var)==0){
        stop('cannot specify both h2 and env.var. Pick one')
    }
}
NormalCGF <- function(s,mu,sigma2,prime=0){
    if(prime==0){
        s^2*sigma2/2 + s*mu
    } else if(prime==1){
        mu+s*sigma2
    } else if (prime==2){
        sigma2
    } else if(prime>2){
        0
    } else {
        stop('prime must be positive')
    }
}
LogisticCGF <- function(s,mu,sigma2,l.scale=NULL,prime=0){
    ## note: if you go on wikipedia, or other resources, you will see 's' parameterizing the
    ## scale of the logistic distribution. Here, 's' is the dummy variable in the
    ## cumulant generating function, and the scale is parameterized by 'l.scale'
    if(is.null(l.scale))
        l.scale <- sqrt(sigma2*3/pi)
    else
        stop('code not yet configured to handle parameterization of logistic scale directly')
    if(prime==0){
        s*mu+beta(1-l.scale*s,1+l.scale*s)
    } else if(prime==1){
        mu + l.scale*(digamma(1+s*l.scale)-digamma(1-s*l.scale))
    } else if (prime==2){
        l.scale^2*(trigamma(1+s*l.scale)+trigamma(1-s*l.scale))
    } else {
        stop('variable prime must by = 0, 1, or 2')
    }
}
WrightPhenoCGF <- function(s=0,myTheta=4e-4,myGamma=1e1,alpha=1,Ne=1e4,L=1e4,ploidy=2,h2=NULL,env.var=NULL,env.mean=0,return.sum=TRUE,env.shape='normal'){
    EnvVarWarning(h2,env.var)
    gen.var <- WrightGenoCGFdoubleprime(s=0,myTheta,myGamma,alpha,Ne,L,ploidy,return.sum=TRUE)
    if(is.null(env.var)){
        env.var <- gen.var*(1-h2)/h2
    }
    geno.cgf <- WrightGenoCGF(s,myTheta,myGamma,alpha,Ne,L,ploidy,return.sum=FALSE)
    if(env.shape=='normal'){
        env.cgf <- NormalCGF(s,env.mean,env.var,prime=0)
    } else if(env.shape=='logistic'){
        env.cgf <- LogisticCGF(s,env.mean,env.var,prime=0)
    }
    if(return.sum)
        sum(c(unlist(geno.cgf),env.cgf))
    else
        c(geno.cgf,env=env.cgf)
}
WrightPhenoCGFprime <- function(s=0,myTheta=4e-4,myGamma=1e1,alpha=1,Ne=1e4,L=1e4,ploidy=2,h2=NULL,env.var=NULL,env.mean=0,return.sum=TRUE,env.shape='normal'){
    EnvVarWarning(h2,env.var)
    geno.cgfps <- WrightGenoCGFprime(s,myTheta,myGamma,alpha,Ne,L,ploidy,return.sum=FALSE)
    gen.var <- WrightGenoCGFdoubleprime(0,myTheta,myGamma,alpha,Ne,L,ploidy,return.sum=TRUE)
    if(is.null(env.var)){
        env.var <- gen.var*(1-h2)/h2
    }
    if(env.shape=='normal'){
        env.cgfps <- NormalCGF(s,env.mean,env.var,prime=1)
    } else if(env.shape=='logistic'){
        env.cgfps <- LogisticCGF(s,env.mean,env.var,prime=1)
    }
    if(return.sum)
        sum(c(unlist(geno.cgfps),env.cgfps))
    else
        c(geno.cgfps,env=env.cgfps)
}
WrightPhenoCGFdoubleprime <- function(s=0,myTheta=4e-4,myGamma=1e1,alpha=1,Ne=1e4,L=1e4,ploidy=2,h2=NULL,env.var=NULL,return.sum=TRUE,analytical=TRUE,env.shape='normal'){
    EnvVarWarning(h2,env.var)
    geno.cgfdps <- WrightGenoCGFdoubleprime(s,myTheta,myGamma,alpha,Ne,L,ploidy,return.sum=F)
    gen.var <- WrightGenoCGFdoubleprime(0,myTheta,myGamma,alpha,Ne,L,ploidy,return.sum=T)
    if(is.null(env.var)){
        env.var <- gen.var*(1-h2)/h2
    }
    if(env.shape=='normal'){
        env.cgfdps <- NormalCGF(s,env.mean,env.var,prime=2)
    } else if(env.shape=='logistic'){
        env.cgfdps <- LogisticCGF(s,env.mean,env.var,prime=2)
    }
    if(return.sum)
        sum(c(unlist(geno.cgfdps),env.cgfdps))
    else
        c(geno.cgfdps,env=env.cgfdps)
}
WrightPhenoCGFtripleprime <- function(s=0,myTheta=4e-4,myGamma=1e1,alpha=1,Ne=1e4,L=1e4,ploidy=2,h2=NULL,env.var=NULL,return.sum=TRUE,analytical=TRUE,env.shape='normal'){
    ##EnvVarWarning(h2,env.var)
    geno.cgftps <- WrightGenoCGFtripleprime(s,myTheta,myGamma,alpha,Ne,L,ploidy,return.sum=F)
    gen.var <- WrightGenoCGFdoubleprime(0,myTheta,myGamma,alpha,Ne,L,ploidy,return.sum=T)
    if(env.shape=='normal'){
        env.cgftps <- 0
    } else if(env.shape=='logistic'){
        env.cgftps <- LogisticCGF(s,env.mean,env.var,prime=3)
    }
    if(return.sum)
        sum(c(unlist(geno.cgftps),env.cgftps))
    else
        c(geno.cgftps,env=env.cgftps)
}
WrightPhenoCGFprimeInverse <- function(thr,eta,myTheta=4e-4,myGamma=1e1,alpha=1,Ne=1e4,L=1e4,ploidy=2,h2=NULL,env.var=NULL,env.mean=0,return.sum=TRUE,env.shape='normal'){

    ## WIP a bit confused
    recover()
    nleqslv(
        x=thr,
        fn=function(x,...)WrightPhenoCGFprime(x,...)-eta,
        myTheta=myTheta,
        myGamma=myGamma,
        alpha=alpha,
        Ne=Ne,
        L=L,
        ploidy=ploidy,
        h2=h2,
        env.var=env.var,
        env.mean=env.mean,
        return.sum=TRUE,
        env.shape=env.shape
    )



}
getTailProbs <- function(my.s,myTheta=4e-4,myGamma=1e1,alpha=1,Ne=1e4,L=1e4,ploidy=2,pheno=TRUE,h2=0.5,env.var=NULL,give.bound=FALSE,env.mean=0,env.shape='normal'){
    ## recover()

    if(FALSE){

        myGamma <- c(1,1)
        alpha <- c(1,1)
        L <- c(250000,250000)


        myGamma <- 1
        alpha <- 1
        L <- 500000

    }



    if(!pheno){
        my.mean <- WrightGenoCGFprime(
            s=0,
            myTheta=myTheta,
            myGamma=myGamma,
            alpha=alpha,
            Ne=Ne,
            L=L,
            ploidy=ploidy
        )
        t0 <- sapply(
            my.s,
            function(s)
                WrightGenoCGF(
                    s=s,
                    myTheta=myTheta,
                    myGamma=myGamma,
                    alpha=alpha,
                    Ne=Ne,
                    L=L,
                    ploidy=ploidy
                )
        )
        t <- sapply(
            my.s,
            function(s)
                WrightGenoCGFprime(
                    s=s,
                    myTheta=myTheta,
                    myGamma=myGamma,
                    alpha=alpha,
                    Ne=Ne,
                    L=L,
                    ploidy=ploidy
                )
        )
        tt <- sapply(
            my.s,
            function(s)
                WrightGenoCGFdoubleprime(
                    s=s,
                    myTheta=myTheta,
                    myGamma=myGamma,
                    alpha=alpha,
                    Ne=Ne,
                    L=L,
                    ploidy=ploidy
                )
        )
        c(my.mean,t0,t,tt)
    } else {
        my.mean <- WrightPhenoCGFprime(
            s=0,
            myTheta=myTheta,
            myGamma=myGamma,
            alpha=alpha,
            Ne=Ne,
            L=L,
            ploidy=ploidy,
            h2=h2,
            env.var=env.var,
            env.mean=env.mean,
            env.shape=env.shape
        )
        t0 <- sapply(
            my.s,
            function(s)
                WrightPhenoCGF(
                    s=s,
                    myTheta=myTheta,
                    myGamma=myGamma,
                    alpha=alpha,
                    Ne=Ne,
                    L=L,
                    ploidy=ploidy,
                    h2=h2,
                    env.var=env.var,
                    env.mean=env.mean,
                    env.shape=env.shape
                )
        )
        t <- sapply(
            my.s,
            function(s)
                WrightPhenoCGFprime(
                    s=s,
                    myTheta=myTheta,
                    myGamma=myGamma,
                    alpha=alpha,
                    Ne=Ne,
                    L=L,
                    ploidy=ploidy,
                    h2=h2,
                    env.var=env.var,
                    env.mean=env.mean,
                    env.shape=env.shape
                )
        )
        tt <- sapply(
            my.s,
            function(s)
                WrightPhenoCGFdoubleprime(
                    s,
                    myTheta=myTheta,
                    myGamma=myGamma,
                    alpha=alpha,
                    Ne=Ne,
                    L=L,
                    ploidy=ploidy,
                    h2=h2,
                    env.var=env.var,
                    env.shape=env.shape
                )
        )
        c(my.mean,t0,t,tt)
    }
    my.log.Q <- pnorm(sqrt(my.s^2*tt),lower.tail=FALSE,log.p=TRUE)
    my.log.term <- t0-my.s*t+my.s^2*tt/2
    my.prod <- exp(my.log.Q+my.log.term)
    tp <- ifelse(
        t>c(my.mean) ,
        my.prod,
        1-(ifelse(my.s>0,1,0)-sign(my.s)*my.prod)
    )

    if(my.s==0){
        return(
            list(
                my.t=t,
                tail.prob=0.5
            )
        )
    }

    if(give.bound){
        ## this bound is basically useless
        ## have not been keeping this bit of code up to date
        exp.bound.term <- exp(t0-my.s*t)
        bns.term <- sapply(my.s,function(S) bns(S,myTheta,myGamma,Ne,L,ploidy))
        bound.term <- exp.bound.term*bns.term
        my.range <- pmax(c(tp-bound.term,tp+bound.term),0)
        return(
            list(
                my.t=t,
                tail.prob=tp,
                range=my.range
            )
        )
    } else {
        return(
            list(
                my.t=t,
                tail.prob=tp
            )
        )
    }
}
getDensity <- function(my.s,myTheta=4e-4,myGamma=1e1,alpha=1,Ne=1e4,L=1e4,ploidy=2,pheno=TRUE,h2=0.5,give.bound=FALSE,env.mean=0,env.shape='normal'){

    recover()
    if(FALSE){

        myGamma <- c(1,1)
        alpha <- c(1,1)
        L <- c(250000,250000)


        myGamma <- 1
        alpha <- 1
        L <- 500000

    }

    if(!pheno){
        my.mean <- WrightGenoCGFprime(
            s=0,
            myTheta=myTheta,
            myGamma=myGamma,
            alpha=alpha,
            Ne=Ne,
            L=L,
            ploidy=ploidy
        )
        t0 <- sapply(
            my.s,
            function(s)
                WrightGenoCGF(
                    s=s,
                    myTheta=myTheta,
                    myGamma=myGamma,
                    alpha=alpha,
                    Ne=Ne,
                    L=L,
                    ploidy=ploidy
                )
        )
        t <- sapply(
            my.s,
            function(s)
                WrightGenoCGFprime(
                    s=s,
                    myTheta=myTheta,
                    myGamma=myGamma,
                    alpha=alpha,
                    Ne=Ne,
                    L=L,
                    ploidy=ploidy
                )
        )
        tt <- sapply(
            my.s,
            function(s)
                WrightGenoCGFdoubleprime(
                    s=s,
                    myTheta=myTheta,
                    myGamma=myGamma,
                    alpha=alpha,
                    Ne=Ne,
                    L=L,
                    ploidy=ploidy
                )
        )
        t3 <- sapply(
            my.s,
            function(s)
                WrightGenoCGFtripleprime(
                    s=s,
                    myTheta=myTheta,
                    myGamma=myGamma,
                    alpha=alpha,
                    Ne=Ne,
                    L=L,
                    ploidy=ploidy
                )
        )
        c(my.mean,t0,t,tt,t3)
    } else {
        my.mean <- WrightPhenoCGFprime(
            s=0,
            myTheta=myTheta,
            myGamma=myGamma,
            alpha=alpha,
            Ne=Ne,
            L=L,
            ploidy=ploidy,
            h2=h2,
            env.mean=env.mean,
            env.shape=env.shape
        )
        t0 <- sapply(
            my.s,
            function(s)
                WrightPhenoCGF(
                    s=s,
                    myTheta=myTheta,
                    myGamma=myGamma,
                    alpha=alpha,
                    Ne=Ne,
                    L=L,
                    ploidy=ploidy,
                    h2=h2,
                    env.mean=env.mean,
                    env.shape=env.shape
                )
        )
        t <- sapply(
            my.s,
            function(s)
                WrightPhenoCGFprime(
                    s=s,
                    myTheta=myTheta,
                    myGamma=myGamma,
                    alpha=alpha,
                    Ne=Ne,
                    L=L,
                    ploidy=ploidy,
                    h2=h2,
                    env.mean=env.mean,
                    env.shape=env.shape
                )
        )
        tt <- sapply(
            my.s,
            function(s)
                WrightPhenoCGFdoubleprime(
                    s,
                    myTheta=myTheta,
                    myGamma=myGamma,
                    alpha=alpha,
                    Ne=Ne,
                    L=L,
                    ploidy=ploidy,
                    h2=h2,
                    env.shape=env.shape
                )
        )
        t3 <- sapply(
            my.s,
            function(s)
                WrightPhenoCGFtripleprime(
                    s,
                    myTheta=myTheta,
                    myGamma=myGamma,
                    alpha=alpha,
                    Ne=Ne,
                    L=L,
                    ploidy=ploidy,
                    h2=h2,
                    env.shape=env.shape
                )
        )
        c(my.mean,t0,t,tt,t3)
    }
    my.Q <- pnorm(sqrt(my.s^2*tt),lower.tail=FALSE)
    my.exp.term <- exp(t0-my.s*t+my.s^2*tt/2)
    one <- my.s^2*my.Q*sqrt(tt)*t3
    two <- dnorm(eta*sqrt(tt))*(2*tt+my.s*t3)
    denom <- 2*sqrt(tt)

    DdDs <- my.exp.term*(one+two)/denom
    tmp <- WrightPhenoCGFprimeInverse(
        thr=t,eta=my.s,myTheta=4e-4,myGamma=1e1,alpha=1,
        Ne=1e4,L=1e4,ploidy=2,h2=h2,env.var=NULL,env.mean=0,
        return.sum=TRUE,env.shape='normal'
    )
    DsDt <-
    return(my.dens)
    ## tp <- ifelse(
    ##     t>c(my.mean) ,
    ##     my.prod,
    ##     1-(ifelse(my.s>0,1,0)-sign(my.s)*my.prod)
    ## )


    ## return(
    ##     if(give.bound){
    ##         ## this bound is basically useless
    ##         ## have not been keeping this bit of code up to date
    ##         exp.bound.term <- exp(t0-my.s*t)
    ##         bns.term <- sapply(my.s,function(S) bns(S,myTheta,myGamma,Ne,L,ploidy))
    ##         bound.term <- exp.bound.term*bns.term
    ##         my.range <- pmax(c(tp-bound.term,tp+bound.term),0)
    ##         list(
    ##             my.t=t,
    ##             tail.prob=tp,
    ##             range=my.range
    ##         )
    ##     } else {
    ##         list(
    ##             my.t=t,
    ##             tail.prob=tp
    ##         )
    ##     }
    ## )
}
getRiskAlleleS <- function(my.s,myTheta,myGamma,Ne,L,ploidy,h2,alpha=1){
    base.t <- WrightPhenoCGFprime(my.s,myTheta,myGamma,Ne,L,ploidy,h2=h2)
    s.range <- seq(my.s-100,my.s,length.out=20)
    t.range <- sapply(s.range,function(S)WrightPhenoCGFprime(S,myTheta,myGamma,Ne,L,ploidy,h2=h2))
    narrow.s.range <- c(
        max(s.range[(base.t-alpha-t.range)>0]),
        min(s.range[(base.t-alpha-t.range)<0])
    )
    mut.s <- uniroot(
        f=function(s) (base.t - alpha) - WrightPhenoCGFprime(s,myTheta,myGamma,Ne,L,ploidy,h2=h2) ,
        lower=narrow.s.range[1],
        upper=narrow.s.range[2]
    )$root
    return(mut.s)
}
getPiDiffFromS <- function(start.s,this.pi,myTheta,myGamma,Ne,L,ploidy,h2,r=FALSE){
    if(r)
        recover()
    base.t <- WrightPhenoCGFprime(start.s,myTheta,myGamma,Ne,L,ploidy,h2=h2)
    mut.s <- getRiskAlleleS(start.s,myTheta,myGamma,Ne,L,ploidy,h2)
    my.probs <- getTailProbs(
        c(mut.s,start.s),
        myTheta=myTheta,
        myGamma=myGamma,
        Ne=Ne,
        L=L,
        ploidy=ploidy,
        h2=h2
    )$tail.prob
    guess.pi <- my.probs[1]-my.probs[2]
    list(
        pi.diff=this.pi-guess.pi,
        guess.pi=guess.pi,
        probs=my.probs,
        base.t=base.t
    )
}
findPrevApprox <- function(myTheta,myGamma,fitCost,Ne,L,ploidy,h2) {

    my.pi <- myGamma/(4*Ne*fitCost)
    pi.diffs <- numeric()
    ss.tried <- numeric()
    pi.guess <- numeric()
    cur.s <- 0
    i <- 1
    j <- 1
    s.inc <- 1/4



    while ( TRUE ) {
        if(FALSE){
            ## optim(
            ##     par=0.001,
            ##     function(s)
            ##         abs(getPiDiffFromS(s,my.pi,myTheta,myGamma,Ne,L,ploidy,h2)$pi.diff),
            ##     method='CG',
            ##     control=list('abstol'=1e-10)
            ## )
            grad(
                func=function(s)
                    getPiDiffFromS(s,my.pi,myTheta,myGamma,Ne,L,ploidy,h2)$pi.diff,
                x=0.01
            )
        }
        out <- getPiDiffFromS(cur.s,my.pi,myTheta,myGamma,Ne,L,ploidy,h2,r=FALSE)
        pi.diffs[i] <- out$pi.diff
        pi.guess[i] <- out$guess.pi
        ss.tried[i] <- cur.s
        if ( length(pi.diffs)>1 ){
            if (pi.diffs[i] > 0 & pi.diffs[i-1] < 0)
                break
            if(pi.diffs[i]==my.pi){
                i <- 1
                pi.diffs <- numeric()
                ss.tried <- numeric()
                cur.s <- 0
                i <- 1
                s.inc <- s.inc/4
                j <- j+1
                next
            }
            if(j>6)
                stop()
        }
        cur.s <- cur.s + s.inc
        i <- i+1
    }
    if(diff(tail(pi.diffs,2))<1e-8){
        min.s <- cur.s - 2*s.inc
        max.s <- cur.s - s.inc
    } else {
        min.s <- cur.s - s.inc
        max.s <- cur.s
    }

    out <- uniroot(
        f = function(S){
            out <- getPiDiffFromS(S,my.pi,myTheta,myGamma,Ne,L,ploidy,h2)
            out$pi.diff
        },
        lower=min.s,
        upper=max.s
    )
    the.s <- out$root
    my.prev <- getTailProbs(the.s,myTheta,myGamma,Ne,L,ploidy,pheno=TRUE,h2)
    return(
        list(
            the.s,
            my.prev
        )
    )
}

getTiltedCDF <- function(interval,length.out=100,lower=interval[1],upper=interval[2],myTheta=4e-4,myGamma=1e1,Ne=1e4,L=1e4,ploidy=2,h2=NULL){
    right.s <- 0
    left.s <- 0
    while(TRUE){
        this.prob <- getTailProbs(my.s=right.s,myTheta,myGamma,Ne,L,ploidy,TRUE,h2)$tail.prob
        if(this.prob<lower)
            break
        right.s <- right.s +  1/2
    }
    while(TRUE){
        this.prob <- getTailProbs(my.s=left.s,myTheta,myGamma,Ne,L,ploidy,TRUE,h2)$tail.prob
        if(this.prob>upper)
            break
        left.s <- left.s +  1/2
    }
    my.s <- seq(from=left.s,to=right.s,length.out=length.out)
    tmp <- sapply(my.s,function(S) getTailProbs(S,myTheta,myGamma,Ne,L,ploidy,TRUE,h2))
    return(rbind(tmp,my.s=my.s))
}

#############################################
## helper functions to solve for threshold ##
#############################################

riskObjective <- function(eta,beta,this.pi,myTheta,myGamma,alpha,Ne,L,ploidy,h2,env.var=NULL){
    ## recover()
    ds <- sapply(
        c(0,beta),
        function(BETA)
            getTailProbs(
                my.s=eta-BETA,
                myTheta=myTheta,
                myGamma=myGamma,
                alpha=alpha,
                Ne=Ne,
                L=L,
                ploidy=2,
                h2=h2,
                env.var=env.var
            )$tail.prob
    )
    ds[2]-ds[1]-this.pi
}

prevObjective <- function(gamma,theta,alpha,Ne,L,h2,fitCost,prev,ploidy=2,env.mean=0,renew.weights=TRUE,env.shape='normal'){
    ##recover()
    solnEta <- uniroot(
        f=function(ETA){
            getTailProbs(
                my.s=ETA,
                myTheta=theta,
                myGamma=gamma,
                alpha=alpha,
                Ne=Ne,
                L=L,
                ploidy=ploidy,
                h2=h2,
                env.mean=env.mean,
                renew.weights=TRUE,
                env.shape=env.shape
            )$tail.prob-prev
        },
        interval=c(0,40),
        tol=1e-13
    )
    myEta <<- solnEta$root
    this.pi <- gamma/(4*Ne*fitCost)
    solnBeta <- uniroot(
        f=function(BETA){
            getTailProbs(
                my.s=myEta-BETA,
                myTheta=theta,
                myGamma=gamma,
                alpha=alpha,
                Ne=Ne,
                L=L,
                ploidy=ploidy,
                h2=h2,
                env.mean=env.mean,
                renew.weights=TRUE,
                env.shape=env.shape
            )$tail.prob-this.pi-prev
        },
        interval=c(0,40),
        tol=1e-13
    )
    myBeta <<- solnBeta$root

    effectSizeObjective(myEta,myBeta,theta,gamma,alpha,Ne,L,ploidy,h2,renew.weights)

}
TwoObjectiveGammaFixed <- function(X,my.pi,myTheta,myGamma,alpha,Ne,L,ploidy,h2){
    ## X[1] = eta
    ## X[2] = beta
    if(any(X<0)){
        ifelse(X<0,X*rep(1e7,2),X)
    } else {
        thr <- WrightPhenoCGFprime(X[1],myTheta,myGamma,alpha,Ne,L,ploidy,h2)
        c(
            riskObjective(X[1],X[2],my.pi,myTheta,myGamma,alpha,Ne,L,ploidy,h2),
            #WrightPhenoCGFprime(X[1],myTheta,myGamma,alpha,Ne,L,ploidy,h2)-X[3],
            WrightPhenoCGFprime(X[1]-X[2],myTheta,myGamma,alpha,Ne,L,ploidy,h2) - (thr - alpha)
        )
    }
}
ThreeObjectiveThrFixed <- function(X,thr,myTheta,alpha,Ne,L,ploidy,h2,fitCost){
    ## X[1] = eta
    ## X[2] = beta
    ## X[3] = gamma

    eta <- X[1]
    beta <- X[2]
    gamma <- X[3]
    this.pi <- gamma/(4*Ne*fitCost)

    out <- c(
        riskObjective(eta,beta,this.pi,myTheta,gamma,alpha,Ne,L,ploidy,h2),
        WrightPhenoCGFprime(eta,myTheta,gamma,alpha,Ne,L,ploidy,h2)-thr,
        WrightPhenoCGFprime(eta-beta,myTheta,gamma,alpha,Ne,L,ploidy,h2)-(thr-alpha)
    )

    return(ifelse(X<0,X*rep(1e50,3),out))

}
MultEffObGammaFixed <- function(X,my.pi,myTheta,myGamma,alpha,rho,Ne,L,ploidy,h2,rt=1:length(X),eta=NULL,betas=NULL,big.alphas=NULL,logit.ps=NULL){
    ## X[1] = eta
    ## X[2:N] = betas
    ## X[(M+1):N] = alphas
    ## X[(N+1):O] = ps
    if(recoversOn)
        recover()

    n.var <- length(X)
    n.effects <- ifelse(n.var==2,1,(n.var+1)/3)

    eta.ind <- 1
    beta.inds <- seq(2,n.effects+1)
    alpha.inds <- seq(max(beta.inds)+1,max(beta.inds)+n.effects-1)
    p.inds <- seq(max(alpha.inds)+1,max(alpha.inds)+(n.effects-1))
    ## thr.ind <- n.var

    if(is.null(eta))
        eta <- X[eta.ind]
    if(is.null(betas))
        betas <- X[beta.inds]
    if(is.null(big.alphas))
        big.alphas <- X[alpha.inds]
    all.alphas <- c(alpha,big.alphas)
    if(is.null(logit.ps)){
        tmp.ps <- plogis(X[p.inds])
        ps <- c(tmp.ps,1-sum(tmp.ps)) #c(X[p.inds],1-sum(X[p.inds]))
    }

    thr <- WrightPhenoCGFprime(eta,myTheta,myGamma,all.alphas,Ne,ps*L,ploidy,h2,return.sum=TRUE)

    out1 <- mapply(
        function(BETA,PI)
            riskObjective(eta,BETA,PI,myTheta,myGamma,all.alphas,Ne,ps*L,ploidy,h2),
        BETA=betas,
        PI=my.pi
    )
    tmp.out2 <- sapply(
        betas,
        function(BETA)
            WrightPhenoCGFprime(eta-BETA,myTheta,myGamma,all.alphas,Ne,ps*L,ploidy,h2,return.sum=TRUE)
    )
    out2 <- tmp.out2-(thr-all.alphas)

    tmp.out3 <- WrightGenoCGFdoubleprime(0,myTheta,myGamma,all.alphas,Ne,ps*L,ploidy,return.sum=FALSE)
    ## out3 <- head(tmp.out3/sum(tmp.out3)-rho,-1)
    out3 <- head((tmp.out3/sum(tmp.out3))/rho-1,-1)

    ## ifelse(X<0,X*rep(1e7,n.var),X)
    my.out <- c(out1/my.pi,out2/all.alphas,out3)
    my.out[rt]

}
MultEffObThrFixed <- function(X,thr,delta,myTheta,alpha,rho,Ne,L,ploidy,h2,rt=1:length(X),eta=NULL,betas=NULL,big.alphas=NULL,logit.ps=NULL,small.gamma=NULL){
    ## X[1] = eta
    ## X[2:N] = betas
    ## X[(M+1):N] = alphas
    ## X[(N+1):O] = ps
    ## X[(O+1):P] = gammas
    if(recoversOn)
        recover()

    n.var <- length(X)
    n.effects <- ifelse(n.var==2,1,(n.var)/3)

    eta.ind <- 1
    beta.inds <- seq(2,n.effects+1)
    alpha.inds <- seq(max(beta.inds)+1,max(beta.inds)+n.effects-1)
    p.inds <- seq(max(alpha.inds)+1,max(alpha.inds)+(n.effects-1))
    gamma.inds <- seq(max(alpha.inds)+(n.effects),max(p.inds)+(n.effects-1))

    if(is.null(eta))
        eta <- X[eta.ind]
    if(is.null(betas))
        betas <- X[beta.inds]
    if(is.null(big.alphas))
        big.alphas <- X[alpha.inds]
    all.alphas <- c(alpha,big.alphas)
    if(is.null(logit.ps)){
        tmp.ps <- plogis(X[p.inds])
        ps <- c(tmp.ps,1-sum(tmp.ps)) #c(X[p.inds],1-sum(X[p.inds]))
    }
    if(is.null(small.gamma)){
        small.gamma <- X[gamma.inds]
        gammas <- c(1,1+delta)*small.gamma
    }
    my.pi <- gammas/(4*Ne*fitCost)

    ## thr <- WrightPhenoCGFprime(eta,myTheta,myGamma,all.alphas,Ne,ps*L,ploidy,h2,return.sum=TRUE)

    out1 <- mapply(
        function(BETA,PI)
            riskObjective(eta,BETA,PI,myTheta,gammas,all.alphas,Ne,ps*L,ploidy,h2),
        BETA=betas,
        PI=my.pi
    )
    tmp.out2 <- sapply(
        c(0,betas),
        function(BETA)
            WrightPhenoCGFprime(eta-BETA,myTheta,gammas,all.alphas,Ne,ps*L,ploidy,h2,return.sum=TRUE)
    )
    out2 <- tmp.out2-(thr-c(0,all.alphas))

    tmp.out3 <- WrightGenoCGFdoubleprime(0,myTheta,gammas,all.alphas,Ne,ps*L,ploidy,return.sum=FALSE)
    ## out3 <- head(tmp.out3/sum(tmp.out3)-rho,-1)
    out3 <- head((tmp.out3/sum(tmp.out3))/rho-1,-1)

    ## ifelse(X<0,X*rep(1e7,n.var),X)
    my.out <- c(out1/my.pi,out2/c(alpha,all.alphas),out3)
    my.out[rt]

}
