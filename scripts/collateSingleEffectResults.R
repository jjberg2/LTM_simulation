## functions
merge_into_paramtable <- function(params.table,file.roots,sim.exts,site.exts,my.path){
    ## recover()
    param.rows <- nrow(params.table)
    n.sims <- length(file.roots)
    if(param.rows!=n.sims)
        stop("number of parameter combinations in parameter table must match files in directory")
    ## things where there's one number per generation
    n.sim.exts <- length(sim.exts)
    sim.results <- matrix(NA,nrow=n.sims,ncol=n.sim.exts)
    colnames(sim.results) <- sapply(sim.exts,function(X) paste("sim",X,sep="."))
    ## things where there's one number per site per generation
    n.site.exts <- length(site.exts)
    site.results <- matrix(NA,nrow=n.sims,ncol=n.site.exts)
    colnames(site.results) <- sapply(site.exts,function(X) paste("sim",X,sep="."))
    for (i in 1:nrow(params.table) ) {
        myNe <- format(round(params.table[i,'Ne'],3),3)
        myL <- floor(params.table[i,'target.size'])
        myrho <- format(round(params.table[i,'rho'],5),nsmall=5)
        mycost  <- format(round(params.table[i,'cost'],2),nsmall=2)
        myenvSD  <- format(round(params.table[i,'env.sd'],5),nsmall=5)
        temp_prefix = paste(my.path,"/PopSize", myNe, "_LiaSize", myL, "_rho", myrho, "_cost", mycost, "_envsd", myenvSD, "_all", sep="")
        sim.files <- sapply(sim.exts,function(X) paste(temp_prefix,X,sep="."))
        sim.results[i,] <- sapply(sim.files,function(X) colMeans(read.table(X)))
        site.files <- sapply(site.exts,function(X) paste(temp_prefix,X,sep="."))
        for(j in 1:length(site.files)){
            site.results[i,j] <- mean(unlist(sapply(readLines(site.files[j]),function(X) as.numeric(strsplit(X,',')[[1]]))))
        }
    }
    results <- cbind(params.table,sim.results,site.results)
    return(results)
}

options(scipen=400)
if(interactive()){
    my.path <- "largeEffectInsensitivity/all"
    my.args <- c(
        "largeEffectInsensitivityParamTable.txt",
        sapply(dir(my.path),function(X) paste(my.path,X,sep="/")),
        my.path,
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

site.exts <- exts[grepl('Freq',exts) | grepl('siteVar',exts) ]
sim.exts  <- exts[!(exts %in% site.exts)]

results.table <- merge_into_paramtable(param.table,files,sim.exts,site.exts,my.path)
save(results.table,file=output.table.filename)


