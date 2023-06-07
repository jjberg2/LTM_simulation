rm(list=ls())
load(file='costInsensitivityResultsTable.Rdata')










pdf('figures/SingleEffectCostInsensitivityFigure.pdf',width=12,height=5)
par(mfrow=c(1,2))
max.prev <- max(results$prev_emp)
my.b <- unique(results$b)
cols <- seq_along(my.b)
plot(
    NA,
    xlim = c(0,1),
    ylim = c(0,0.05),
    bty='n',
    xlab='',
    ylab=''
)
abline(h=1/2,lty=2)
mtext('Fitness Cost of Disease',side=1,line=2.5,cex=2)
mtext('Disease Prevalence',side=2,line=2.5,cex=2)
for ( i in seq_along(my.b)){
    b = my.b[i]
    costs <- results[results[,'b']==b,'cost']
    prev_theory <- results[results[,'b']==b,'norm.prev']
    prev_emp <- results[results[,'b']==b,'prev_emp']
    points(
        x=costs,
        y=prev_emp,
        pch=20,
        col=i,
        cex=1.6
    )
    lines(
        x=costs,
        y=prev_emp,
        pch=20,
        col=i,
        cex=1.6
    )
     legend(
        x=0.8,
        y=0.044,
        legend=round(my.b,2),
        pch=20,
        col=cols,
        bty='n',
        cex=1.6
    )
    text(
        x=0.87,
        y=0.048,
        labels='Mutational',
        cex=1.6
    )
    text(
        x=0.87,
        y=0.044,
        labels='Asymmetry',
        cex=1.6
    )
}






plot(
    NA,
    xlim = c(0,1),
    ylim = c(0,1),
    bty='n',
    xlab='',
    ylab=''
)
abline(h=1/2,lty=2)
mtext('Fitness Cost of Disease',side=1,line=2.5,cex=2)
mtext('Heritability (liability scale)',side=2,line=2.5,cex=2)

for ( i in seq_along(my.b)){
    b = my.b[i]
    costs <- results[results[,'b']==b,'cost']
    h2s <- results[results[,'b']==b,'sim_h2']
    points(
        x=costs,
        y=h2s,
        pch=20,
        col=i,
        cex=1.6
    )
    lines(
        x=costs,
        y=h2s,
        pch=20,
        col=i
    )
   
}
dev.off()
