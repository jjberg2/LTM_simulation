input <- get(load('smallEffectInsensitivityResultsTable.Rdata'))


table <- input[input$sim.mean<input$thr,]
split.table <- split(table,table$thr)


par(mfrow=c(2,4))
thrs <- sort(unique(table$thr))
b <- sort(unique(table$b),decreasing=TRUE)
for(i in 1:length(split.table)){


    this.cost <- split.table[[i]]$cost
    this.mean <- split.table[[i]]$sim.mean
    this.fixed <- split.table[[i]]$sim.fixed


    plot(
        this.cost,
        this.mean,
        pch=20,
        ylim=c(floor(min(this.mean,this.fixed)),thrs[i]),
        main=paste('b = ' , b[i], sep ='')
    )
    points(
        this.cost,
        this.fixed,
        pch=20,
        col='darkgreen'
    )
    abline(
        h=thrs[i],
        lty=2
    )

}
