library('numDeriv')
prevs = c(0.03,0.01,0.003,0.001)
thr = sapply(prevs, function(PREV) qnorm(1-PREV))
max.a = thr + max(thr)
min.a = rep(0.001,4)
N = 10000
cost = 0.5

my.a = mapply(function(MAXA,MINA) seq(MINA,MAXA,length.out=1000), MINA = min.a, MAXA = max.a,SIMPLIFY = FALSE)

bgamma = function(gamma) {
    ifelse(gamma>20,1,(exp(gamma)-1)/(exp(gamma)+1))
}

normDel= function(a,thr){
    pnorm(thr) - pnorm(thr - a)
}
sigmaa = function(a,thr,nc){
    del = normDel(a,thr)
    a^2 * bgamma(4*nc*del) / (2*nc*del)
}

my.del = mapply( function(THR,A) normDel(A,THR), THR = thr, A = my.a)
my.vars = mapply(function(THR,A) sigmaa(A,THR,N*cost), THR = thr, A = my.a)
my.derivs = list()
for( j in 1:length(thr)){
    my.derivs[[j]] = matrix(nrow=length(my.a[[j]]),ncol=2)
    for( i in 1:length(my.a[[j]])){
        my.derivs[[j]][i,] = genD(
            func = function(X) sigmaa(X,thr[j],N*cost),
            x = my.a[[j]][i]
        )$D
    }
}



cex.lab = 1.4

pdf('figures/suppFigures/largeEffectSizeVarianceRelationship.pdf',width=12,height=6)
par(mfrow=c(1,2))
op = par(mar=c(5,5,4,4)+0.1)

##1 
plot(
    NA,
    xlim = c(0,max(max.a)),
    ylim = c(0,1),
    type = 'l',
    xlab = '',
    ylab = ''
)
for(i in 1:length(my.a)){
    lines(
        my.a[[i]],
        my.del[,i],
        lty = i
    )
}
mtext(
    text=expression(paste('Standardized effect size (', alpha[std] ,')',sep='')),
    side=1,
    line=3,
    cex=cex.lab
)
mtext(
    text=expression(paste('Risk scale effect size (', delta[R](alpha) ,')',sep='')),
    side=2,
    line=3,
    cex=cex.lab
)
text(
    x = 0.8,
    y = 0.9,
    labels = 'Prevalence'
)
legend(
    x = -0.1,
    y = 0.9,
    legend = prevs,
    lty = 1:4,
    bty = 'n'
)

##2
par(op)
plot(
    NA,
    xlim = c(0,max(max.a)),
    ylim = c(0,max(unlist(my.vars))*1.05),
    type = 'l',
    xlab ='',
    ylab =''
)
for(i in 1:length(my.a)){
    lines(
        my.a[[i]],
        my.vars[,i],
        lty = i
    )
}
mtext(
    text=expression(paste('Standardized effect size (', alpha[std] ,')',sep='')),
    side=1,
    line=3,
    cex=cex.lab
)
mtext(
    text='Liability scale contribution to variance',
    side=2,
    line=3.75,
    cex=cex.lab
)
mtext(
    text=expression(alpha^2*b(alpha)/(2*N*delta[R](alpha)*C)),
    side=2,
    line=2,
    cex=cex.lab
)
dev.off()





source('scripts/solveTwoEffect.R')
soln.in <- read.table("twoEffectPrevInsensitivityParamsTable.txt",header = T)

keep <- c(2,which.min(abs(soln.in$deltal - 1/2)))
my.solns <- soln.in[keep,]

tol <- 1e-5

my.dists <- list()
for( i in 1:2){
    min.gl <- uniroot(
        function(X)
            tol-(1-pPoisConv(X,my.solns[i,]$mean.nl,my.solns[i,]$norm.sd,alphal=my.solns[i,]$al)),
        interval=c(-10*my.solns[i,]$tstar,10*my.solns[i,]$tstar)
    )$root
    max.gl <- uniroot(
        function(X)
            tol-pPoisConv(X,my.solns[i,]$mean.nl,my.solns[i,]$norm.sd,alphal=my.solns[i,]$al),
        interval=c(-10*my.solns[i,]$tstar,10*my.solns[i,]$tstar)
    )$root
    seq.li <- seq(min.gl,max.gl,length.out=1000)
    li.dense <- sapply(seq.li,function(G) dPoisConv(G,my.solns[i,]$mean.nl,my.solns[i,]$norm.sd,alphal=my.solns[i,]$al))
    my.dists[[i]] <- list(seq.li=seq.li,li.dense=li.dense)
}    


plot(
    my.dists[[1]]$seq.li,
    my.dists[[1]]$li.dense,
    type='l'
)
lines(
    my.dists[[2]]$seq.li,
    my.dists[[2]]$li.dense,
    type='l',
    col='red'
)


min.gl <- uniroot(
    function(X)
        tol-(1-pPoisConv(X,my.solns[1,]$mean.nl,my.solns[1,]$norm.sd,alphal=my.solns[1,]$al)),
    interval=c(-10*my.solns[1,]$tstar,10*my.solns[1,]$tstar)
)$root
