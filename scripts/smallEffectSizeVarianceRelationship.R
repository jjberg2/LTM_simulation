
c = 0.5
Ne = 5000
bt = 0.5
my.gamma = log((1+bt)/(1-bt))
ft = my.gamma/(4*Ne*c)

d0fa = function(a) a*(exp(a)-1)/(exp(a)+1)
d1fa = function(a) (2*a*exp(a) + exp(2*a) - 1) / ((1 + exp(a))^2)
d2fa = function(a) - 2*exp(a)*(a*(exp(a)-1) - 2*(1 + exp(a) ) ) / ((1 + exp(a))^3)

gamma.inflect = uniroot(d2fa,interval=c(2,4))$root
bt.inflect = uniroot(function(x) gamma.inflect - log((1+x)/(1-x)),lower=0,upper=1)$root


cex.axis = 1.5
cex.lab = 1.4

pdf('figures/suppFigures/smallEffectSizeVarianceRelationship.pdf',width=12,height=5)
par(mfrow=c(1,3))
op = par(mar=c(5,5,4,4)+0.1)
my.a = seq(0,6,length.out=1000)
plot(
  my.a ,
  d0fa(my.a),
  type = 'l',
  xlab = '',
  ylab =  '',
  cex.axis = cex.axis,
  cex.lab = cex.lab
)
mtext(
  text=expression(paste('Scaled selection coefficient (', gamma ,')',sep='')),
  side=1,
  line=3,
  cex=cex.lab
)
mtext(
  text=expression(paste('Contribution to variance (', sigma^2 ,'(', gamma ,'))',sep='')),
  side=2,
  line=2.1,
  cex=cex.lab
)

par(op)
plot(
  my.a ,
  d1fa(my.a),
  type = 'l',
  xlab = '',
  ylab = '',
  cex.axis = cex.axis,
  cex.lab = cex.lab
)
mtext(
  text=expression(paste('Scaled selection coefficient (', gamma ,')',sep='')),
  side=1,
  line=3,
  cex=cex.lab
)
mtext(
  text='First derivative',
  side=2,
  line=2.9,
  cex=cex.lab
)
abline(h=1,lty=2)

plot(
  my.a ,
  d2fa(my.a),
  type = 'l',
  xlab = '',
  ylab = '',
  cex.axis = cex.axis,
  cex.lab = cex.lab
)
mtext(
  text='Second derivative',
  side=2,
  line=2.9,
  cex=cex.lab
)
abline(h=0,lty=2)
mtext(
  text=expression(paste('Scaled selection coefficient (', gamma ,')',sep='')),
  side=1,
  line=3,
  cex=cex.lab
)
dev.off()
