## functions
merge_into_paramtable <- function(params.table,file.roots,file.exts,my.path){
    ## recover()
    param.rows <- nrow(params.table)
    n.sims <- length(file.roots)
    if(param.rows!=n.sims)
        stop("number of parameter combinations in parameter table must match files in directory")
    n.exts <- length(file.exts)
    sim.results <- matrix(NA,nrow=n.sims,ncol=n.exts)
    colnames(sim.results) <- sapply(file.exts,function(X) paste("sim",X,sep="."))
    for (i in 1:nrow(params.table) ) {
        myNe <- format(round(params.table[i,'Ne'],3),3)
        myL <- format(round(params.table[i,'target.size'],3),3)
        myrho <- format(round(params.table[i,'rho'],5),nsmall=5)
        mycost  <- format(round(params.table[i,'cost'],2),2)
        myenvSD  <- format(round(params.table[i,'env.sd'],3),3)
        temp_prefix = paste(my.path,"/PopSize", myNe, "_LiaSize", myL, "_rho", myrho, "_cost", mycost, "_envsd", myenvSD, "_all", sep="")
        files.to.read <- sapply(file.exts,function(X) paste(temp_prefix,X,sep="."))
        sim.results[i,] <- sapply(files.to.read,function(X) colMeans(read.table(X)))
    }
    results <- cbind(params.table,sim.results)
    return(results)
}

options(scipen=400)
if(interactive()){
    my.args <- c(
        "largeEffectInsensitivityParamTable.txt",
        sapply(dir(my.path),function(X) paste(my.path,X,sep="/")),
        "largeEffectInsensitivity/all",
        "largeEffectInsensitivityResultsTable.Rdata"
    )
} else {
    ## read in command line arguments
    my.args <- commandArgs(trailingOnly=TRUE)
}


input.table.filename <- head(my.args,1)
output.table.filename <- tail(my.args,1)
my.path <- head(tail(my.args,2),1)
sim.filenames <- tail(head(my.args,-2),-1)

## load the parameter table
param.table <- read.table(input.table.filename)

## get unique sims and all the file extensions
files <- unique(sapply(sim.filenames,function(X) strsplit(X,'_all.')[[1]][1]))
exts <-  unique(sapply(sim.filenames,function(X) strsplit(X,'_all.')[[1]][2]))


results.table <- merge_into_paramtable(param.table,files,exts,my.path)
save(results.table,file=output.table.filename)


