thr = seq(qnorm(0.98),6,length=1000)

ft = dnorm(thr)
prev = 1-pnorm(thr)

xs = range(ft)
ys = range(prev)

my.axis.cex = 1.8
pdf(paste('figures/suppFigures/prevalenceDensityRelationship.pdf', sep = ''),width=20,height=12)
par(mar=c(5,4,4,2)+0.1)
plot(
  ft,
  prev,
  ylim = c(0,0.02),
  type = "l",
  xlab = '',
  ylab = '',
  cex.axis = my.axis.cex,
  bty = 'n',
  lwd = 3
)
lines(
  x = xs,
  y = ys, 
  lty = 2,
  lwd = 1.5
)
mtext(
  side = 1 ,
  text = 'Threshold density of standard Normal',
  line = 2.7,
  cex = my.axis.cex
)
mtext(
  side = 2 ,
  text = 'Upper tail probability of standard Normal',
  line = 2.7,
  cex = my.axis.cex
)
dev.off()
