## functions
merge_into_paramtable <- function(params.table,file.roots,sim.exts,site.exts,my.path){
    ## recover()
    param.rows <- nrow(params.table)
    n.sims <- length(file.roots)
    if(param.rows!=n.sims)
        stop(paste("number of parameter combinations in parameter table must match files in directory. There are ", param.rows, " rows, but there are ", n.sims, " files.", sep=""))
    ## things where there's one number per generation
    n.sim.exts <- length(sim.exts)
    sim.results <- matrix(NA,nrow=n.sims,ncol=n.sim.exts)
    colnames(sim.results) <- sapply(sim.exts,function(X) paste("sim",X,sep="."))
    ## things where there's one number per site per generation
    n.site.exts <- length(site.exts)
    site.results <- matrix(NA,nrow=n.sims,ncol=n.site.exts)
    colnames(site.results) <- sapply(site.exts,function(X) paste("sim",X,sep="."))
    for (i in 1:nrow(params.table) ) {
        myNe <- format(round(params.table[i,'Ne'],3),nsmall=0)
        myaL <- format(round(params.table[i,'al'],3),nsmall=3)
        myLl <- format(round(params.table[i,'Ll'],3),nsmall=3)
        mythr <- format(round(params.table[i,'thr'],1),nsmall=1)
        mycost  <- format(round(params.table[i,'C'],2),nsmall=2)
        myenvSD  <- format(round(sqrt(params.table[i,'Ve']),3),nsmall=3)
        temp_prefix = paste(my.path,"/PopSize", myNe, "_aL", myaL, "_thr", mythr, "_envSD", myenvSD,  "_cost", mycost, "_all", sep="")
        sim.files <- sapply(sim.exts,function(X) paste(temp_prefix,X,sep="."))
        sim.results[i,] <- sapply(
            sim.files,
            function(X) {
                my.file <- file(X,'r')
                data <- read.table(X)
                close(my.file)
                mean(as.numeric(data[[1]]),na.rm=TRUE)
            }
        )
	## sapply(sim.files,function(X) close(file(X)))
        ## was causing memory issues on cluster so temporarilily deleted
        site.files <- sapply(site.exts,function(X) paste(temp_prefix,X,sep="."))
        ## print(i)
        for(j in 1:length(site.files)){
            my.file <- file(site.files[j],'r')
            lines <- list()
            my.mean <- 0
            k <- 1
            n.tot <- 0
            while(TRUE){
                line <- readLines(my.file,n=1)
                if(length(line)==0){
                    close(my.file)
                    break
                }
                if(line=="NA"){
                    next
                }
                obs <- as.numeric(strsplit(line,',')[[1]])
                this.mean <- mean(obs)
                n.obs <- length(obs)
                new.tot <- n.tot + n.obs
                my.mean <- my.mean*n.tot/(new.tot) + this.mean*n.obs/(new.tot)
                n.tot <- new.tot
                k <- k+1
            }
            if(k==1){
                site.results[i,j] <- NA
            } else {
                site.results[i,j] <- my.mean
            }
        }
        if(i %% 10 ==0 ) print(i)
    }
    results <- cbind(params.table,sim.results,site.results)
    return(results)
}

options(scipen=400)
if(interactive()){
    my.path <- "Users/jeremyberg/Documents/academics/LTM_simulation/twoEffectPrevN1000/all"
    my.args <- c(
        "paramFiles/twoEffectPrevParamTableN1000.txt",
        "twoEffectPrevN1000/all/filenames.txt",
        my.path,
        "resultsFiles/twoEffectPrevResultsTableN1000.Rdata"##,
##        "resultsFiles/twoEffectPrevDerProbsN1000.Rdata"
    )
} else {
    ## read in command line arguments
    my.args <- commandArgs(trailingOnly=TRUE)
}


input.table.filename <- my.args[1]
filenames.file <- my.args[2]
my.path <- my.args[3]
output.table.filename <- my.args[4]


sim.filenames <- readLines(filenames.file)

## load the parameter table
param.table <- read.table(input.table.filename,header=T)

## get unique sims and all the file extensions
files <- unique(sapply(sim.filenames,function(X) strsplit(X,'_all.')[[1]][1]))
exts <-  unique(sapply(sim.filenames,function(X) strsplit(X,'_all.')[[1]][2]))

site.exts <- exts[grepl('Freq',exts) | grepl('siteVar',exts) ]
sim.exts  <- exts[!(exts %in% site.exts)]

results <- merge_into_paramtable(param.table,files,sim.exts,site.exts,my.path)
##sim.results <- results[[2]]
save(results,file=output.table.filename)


