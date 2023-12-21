#setwd('/Users/jeremyberg/Documents/academics/LTM_simulation/')

##########################
### Figure 1 schematic ###
##########################

h2 <- 0.6
source('scripts/figureFuncs.R')
png('figures/paperFiguresForRealThisTime/SchematicFigure1.png', height = 18 , width = 22, units = 'cm', res = 500)
nf <- layout( matrix(c(1,2,3,4), nrow=2, byrow=TRUE))


## mid zoom liability distribution
op2 <- par(mar=c(3,3,2,4))
makePhenLi(xlim=c(-3,4),prev=0.01)
mtext(side = 1 , text = 'Liability', line = 1.5)
mtext(side = 2 , text = 'Density on Liability', line = 1.5)
mtext(side = 3 , text = 'A' , line = 0, at = -3, cex = 1.5)
par(op2)


this.col = 'darkorchid4'
op3 <- par(mar=c(3,4,2,4))
this.ymax <- makeGenLi(xlim=c(-3,4),h2=0.6)
addRiskCurve(h2=h2,prev=0.01,y.max = this.ymax,my.col = this.col)
axis(side=4,at=seq(0,this.ymax,by=0.2*this.ymax),labels = seq(0,1,by = 0.2),las=1)
mtext(side = 1 , text = 'Genetic Liability', line = 1.5)
mtext(side = 2 , text = 'Density on Genetic Liability', line = 2.2)
mtext(side = 4 , text = 'Genetic Risk', line = 2.5, col = this.col,las = 0)
mtext(side = 3 , text = 'B' , line = 0, at = -3, cex = 1.5)
par(op3)


this.col = 'firebrick4'
op4 <- par(mar=c(4,3,1,4))
this.ymax <- makeGenLi(xlim=c(-3,4),h2=0.6)
addRiskCurve(h2=h2,prev=0.01,y.max = this.ymax,my.col = this.col,fit.surface = TRUE)
mtext(side = 1 , text = 'Genetic Liability', line = 2.2)
mtext(side = 2 , text = 'Density on Genetic Liability', line = 1.5)
mtext(side = 4 , text = 'Expected fitness', line = 2.5, col = this.col,las = 0)
axis(side=4,at=seq(0,this.ymax,by=0.2*this.ymax),labels = seq(0,1,by = 0.2),las=1)
mtext(side = 3 , text = 'C' , line = 0, at = -3, cex = 1.5)




prev <- 0.01
t.pos <- qnorm(1-prev,mean=0,sd=1)
my.risks <- exp(seq ( log(1e-9) , log(1-1e-6) , length.out = 1e5 ))
my.g <- qnorm(my.risks , mean = t.pos , sd = sqrt ( 1 - h2 ))
pg <- dnorm(my.g,0,h2)
dgdr <- 1 / dnorm ( qnorm (my.risks, t.pos , sd = sqrt ( 1 -h2 )) , 0, h2 )
pr <- pg / dgdr
plot.risks <- c(0,my.risks,0)
plot.pr <- c(0,pr,0)

this.col = 'darkorchid4'
op4 <- par(mar=c(4,4,1,4))
plot(
    NA,
    type = 'l',
    xlim = c ( 0 , 0.2),
    ylim = c(0,0.5),
    ylab = '',
    xlab = '',
    bty = 'n'
)
mtext(side = 1 , text = 'Genetic Risk', line = 2.2)
mtext(side = 2 , text = 'Density on Genetic Risk', line = 2.2)
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
    density = 120 ,
    angle = 315
)
mtext(side = 3 , text = 'D' , line = 0, at = 0, cex = 1.5)
dev.off()



##########################
### Figure 2 schematic ###
##########################


png('figures/paperFiguresForRealThisTime/SchematicFigure2.png', height = 18 , width = 22, units = 'cm', res = 500)
op1 <- par(mar=c(3.3,3.3,2,0)+0.1)
source('scripts/figureFuncs.R')
y.max <- 0.06
t.pos <- qnorm(0.99)
alpha <- 0.1
makePhenWEffect(xlim=c(2,2.6),y.max=y.max,prev=0.1,alpha=alpha,t.pos=t.pos,xaxt = 'n')
axis(side = 1 , labels = FALSE)
##axis(side = 2 , labels = FALSE)
polygon( 
  x = c(t.pos-alpha , t.pos-alpha , t.pos , t.pos ) ,
  y = c(0 , dnorm(t.pos), dnorm(t.pos),0),
  col = rgb(1,0,1,0.8),
  border = NA
  )
x.pos <- 2.58
y.pos <- 0.03
div <- 4
polygon( 
  x = c(x.pos-alpha/div , x.pos-alpha/div , x.pos , x.pos ) ,
  y = c(y.pos , y.pos+dnorm(t.pos)/div, y.pos+dnorm(t.pos)/div,y.pos),
  col = rgb(1,0,1,0.8),
  border = NA
)
text(x = x.pos - 0.08, y = y.pos + 0.0026 , labels = 'a f(T) = ',cex = 2.5)
mtext("T", at = t.pos, cex = 3, family = 'Arial')
mtext("T-a", at = t.pos - alpha, cex = 3, family = 'Arial')
mtext('Liability' , side = 1, line = 2, cex = 3)
lines(x=c(2,t.pos),y = rep(dnorm(t.pos),2),lty = 3,lwd=3)
mtext(side = 2 , at = dnorm(t.pos), text = 'f(T)', line = -1, las = 1, cex = 3)
##text(x = t.pos - alpha/2, y = y.max / 2.7, '{', srt = 270, cex = 8, family = 'Helvetica Neue UltraLight')
par(op1)
dev.off()











##########################
### Figure 3 schematic ###
##########################

## equilibrium bias figure
{
png('figures/paperFiguresForRealThisTime/AsymmetryFigure3.png', height = 10 , width = 11, units = 'cm', res = 500)
op3 <- par(mar=c(3.3,3.3,2,0.4)+0.1)
op4 <- options(scipen=400)
my.gamma <- 10^seq(-3,3,length.out=1000)
bias <- ifelse(my.gamma<700,(exp(my.gamma)-1)/(exp(my.gamma)+1),1)
plot(
    x = my.gamma ,
    y = bias ,
    log = 'x' ,
    type ='n' ,
    bty = 'n' ,
    lwd = 2 ,
    xlab = '' ,
    ylab = '',
    xaxt = 'n' 
)
mtext(text='Fixation Asymmetry', side = 2, line = 2)
mtext(text='Population Scaled Selection Coefficient', side = 1, line = 2)
axis(side = 1 , at = c(0.001,0.01,0.1,1,10,100,1000),labels = c('0.001','0.01','0.1','1','10','100','1000'))
polygon(
    x=c(1e1,1e1,1e3,1e3),
    y=c(0,1,1,0),
    border=NA,
    col=adjustcolor('grey',alpha.f=0.6)
)
polygon(
    x=c(1e-1,1e-1,1e-3,1e-3),
    y=c(0,1,1,0),
    border=NA,
    col=adjustcolor('grey',alpha.f=0.6)
)
##mtext(text='Effectively',side = 3 , line = 1,at = 0.01)
##mtext(text='neutral',side = 3 , line = 0,at = 0.01)
text(labels='Effectively', x = 0.01, y=0.2)
text(labels='neutral', x = 0.01, y=0.1)
##mtext(text='Weakly',side = 3 , line = 1,at = 1)
##mtext(text='selected',side = 3 , line = 0,at = 1)
text(labels='Weakly', x = 0.4, y=0.9)
text(labels='selected', x = 0.4, y=0.8)
##mtext(text='Strongly',side = 3 , line = 1,at = 100)
##mtext(text='selected',side = 3 , line = 0,at = 100)
text(labels='Strongly', x = 100, y=0.9)
text(labels='selected', x = 100, y=0.8)
lines(
    x=my.gamma,
    y=bias,
    lwd=2
)


par(op3)
options(op4)
dev.off()
}



png('figures/paperFiguresForRealThisTime/MeanEffectFigure4.png', height = 10 , width = 11, units = 'cm', res = 500)
op4 <- par(mar=c(3.3,3.3,2,0.4)+0.1)
my.gamma <- 10^seq(-2,1,length.out=1000)
bias <- (exp(my.gamma)-1)/(exp(my.gamma)+1)
plot(
    x=bias,
    y=my.gamma,
    type='n',
    bty='n',
    lwd=2,
    xlab='',
    ylab=''
)
mtext(expression(paste('Relative threshold position, ', b[T] ,sep = '')),side = 1 , line = 2)
mtext(expression(paste('Population scaled selection coefficient, ', gamma, '(a)', sep = '')),side = 2 , line = 2)
lines(
    x=bias,
    y=my.gamma,
    lwd=2
)
par(op4)
dev.off()





png('figures/paperFiguresForRealThisTime/SelectiveDeathsPer2NFigure5.png', height = 10 , width = 11, units = 'cm', res = 500)
op4 <- par(mar=c(3.3,4.3,2,0.4)+0.1)
my.gamma <- 10^seq(-2,log(6,10),length.out=1000)
bias <- (exp(my.gamma)-1)/(exp(my.gamma)+1)
plot(
    x=my.gamma,
    y=bias*my.gamma,
    type='n',
    bty='n',
    lwd=2,
    xlab='',
    ylab=''
)
mtext('Selective deaths due to additive effects',side = 2 , line = 3 )
mtext(expression(paste('per 2N carriers of the derived allele, b(a)', gamma, '(a)',sep = '')), side = 2 , line = 2)
mtext(expression(paste('Population scaled selection coefficient, ', gamma, '(a)', sep = '')),side = 1 , line = 2)
lines(
    x=my.gamma,
    y=bias*my.gamma,
    lwd=2
)
abline(a = 0 , b = 1 , lty = 3 )
par(op4)
dev.off()

