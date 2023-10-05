library('numDeriv')
thr = qnorm(0.99)
N = 10000
cost = 0.5
my.a = seq(2.5e-2,thr*2,length.out=1000)

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

my.del = normDel(my.a,thr)
my.vars = sigmaa(my.a,thr,N*cost)
my.derivs = sapply(
  my.a,
    function(Y){
      genD(
        func = function(X) sigmaa(X,thr,N*cost),
        x = Y
      )$D
    }
)

pdf('figures/suppFigures/largeEffectSizeVarianceRelationship.pdf',width=12,height=12)
par(mfrow=c(2,2))
op = par(mar=c(5,5,4,4)+0.1)

#1 
plot(
  my.a,
  my.del,
  type = 'l',
  xlab = '',
  ylab = ''
)
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

#2
par(op)
plot(
  my.a,
  my.vars,
  type = 'l',
  xlab ='',
  ylab =''
)
mtext(
  text=expression(paste('Standardized effect size (', alpha[std] ,')',sep='')),
  side=1,
  line=3,
  cex=cex.lab
)
mtext(
  text='Variance on liability scale',
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

#3
plot(
  my.a,
  my.derivs[1,],
  type = 'l',
  xlab = '',
  ylab = ''
)
abline(h=0,lty=2,lwd=2)

#4
plot(
  my.a,
  my.derivs[2,],
  type = 'l',
  xlab = '',
  ylab = ''
)
abline(h=0,lty=2,lwd=2)

dev.off()
