

### Figure 1 schematic

h2 <- 0.6
source('scripts/figureFuncs.R')
pdf('figures/paperFiguresForRealThisTime/SchematicFigure1.pdf', height = 18/2.54 , width = 22/2.54)
nf <- layout( matrix(c(1,1,2,3,4,5), nrow=3, byrow=TRUE), heights = c(1,4,4))

## zoomed out liability distribution
op1 <- par(mar=c(0,0,1,0))
plot.new()
par(op1)

## mid zoom liability distribution
op2 <- par(mar=c(3,3,1,4))
makePhenLi(xlim=c(-3,4),prev=0.01)
mtext(side = 1 , text = 'Liability', line = 1.5)
mtext(side = 2 , text = 'Density on Liability', line = 1.5)
mtext(side = 3 , text = 'B' , line = 0, at = -3, cex = 1.5)
par(op2)


this.col = 'darkorchid4'
op3 <- par(mar=c(3,3,1,4))
this.ymax <- makeGenLi(xlim=c(-3,4),h2=0.6)
addRiskCurve(h2=h2,prev=0.01,y.max = this.ymax,my.col = this.col)
axis(side=4,at=seq(0,this.ymax,by=0.2*this.ymax),labels = seq(0,1,by = 0.2),las=1)
mtext(side = 1 , text = 'Genetic Liability', line = 1.5)
mtext(side = 2 , text = 'Density on Genetic Liability', line = 1.5)
mtext(side = 4 , text = 'Genetic Risk', line = 2.5, col = this.col,las = 0)
mtext(side = 3 , text = 'C' , line = 0, at = -3, cex = 1.5)
par(op3)


this.col = 'firebrick4'
op4 <- par(mar=c(4,3,1,4))
this.ymax <- makeGenLi(xlim=c(-3,4),h2=0.6)
addRiskCurve(h2=h2,prev=0.01,y.max = this.ymax,my.col = this.col,fit.surface = TRUE)
mtext(side = 1 , text = 'Genetic Liability', line = 1.5)
mtext(side = 2 , text = 'Density on Genetic Liability', line = 1.5)
mtext(side = 4 , text = 'Expected fitness', line = 2.5, col = this.col,las = 0)
axis(side=4,at=seq(0,this.ymax,by=0.2*this.ymax),labels = seq(0,1,by = 0.2),las=1)
mtext(side = 3 , text = 'D' , line = 0, at = -3, cex = 1.5)



this.col = 'darkorchid4'
op4 <- par(mar=c(4,3,1,4))
prev <- 0.01
t.pos <- qnorm(1-prev,mean=0,sd=1)

my.risks <- exp(seq ( log(1e-9) , log(1-1e-6) , length.out = 1e5 ))

my.g <- qnorm(my.risks , mean = t.pos , sd = sqrt ( 1 - h2 ))



pg <- dnorm(my.g,0,h2)
dgdr <- 1 / dnorm ( qnorm (my.risks, t.pos , sd = sqrt ( 1 -h2 )) , 0, h2 )

pr <- pg / dgdr

plot.risks <- c(0,my.risks,0)
plot.pr <- c(0,pr,0)

plot(
    NA,
    type = 'l',
    xlim = c ( 0 , 0.2),
    ylim = c(0,0.5),
    bty = 'n',
    yaxt = 'n',
    xlab = ''
)
polygon(
    x = plot.risks ,
    y = plot.pr,
    col = adjustcolor(this.col,alpha.f=0.05),
    border = NA
)
polygon(
    x = plot.risks ,
    y = plot.pr,
    col = adjustcolor(this.col,alpha.f=0.4),
    border = NA,
    density = 60 ,
    angle = 315
)
axis(side = 2 , at = seq(0,1,0.2)*0.5,labels = FALSE)
mtext(side = 1 , text = 'Genetic Risk', line = 2.5,las = 0)
mtext(side = 2 , text = 'Density on Genetic Risk', line = 1.5)
mtext(side = 3 , text = 'E' , line = 0, at = 0, cex = 1.5)
dev.off()

























## equilibrium bias figure

pdf('figures/paperFiguresForRealThisTime/AsymmetryFigure2.pdf',width = 12,height = 6)
my.gamma <- 10^seq(-2,2,length.out=1000)
bias <- (exp(my.gamma)-1)/(exp(my.gamma)+1)


par(mfrow = c(1,2))
plot(
    x = my.gamma ,
    y = bias ,
    log = 'x' ,
    type ='n' ,
    bty = 'n' ,
    lwd = 2 ,
    xlab = 'Population Scaled Selection Coefficient' ,
    ylab = 'Fixation Asymmetry'
)
polygon(
    x=c(1e1,1e1,1e2,1e2),
    y=c(0,1,1,0),
    border=NA,
    col='grey'
)
polygon(
    x=c(1e-1,1e-1,1e-2,1e-2),
    y=c(0,1,1,0),
    border=NA,
    col='grey'
)
lines(
    x=my.gamma,
    y=bias,
    lwd=2
)


my.gamma <- 10^seq(-2,1,length.out=1000)
bias <- (exp(my.gamma)-1)/(exp(my.gamma)+1)



plot(
    x=bias,
    y=my.gamma,
    type='n',
    bty='n',
    lwd=2,
    xlab='Relative Threshold Position',
    ylab='Population Scaled Selection Coefficient'
)
lines(
    x=bias,
    y=my.gamma,
    lwd=2
)
dev.off()





