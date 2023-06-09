## functions
merge_into_paramtable <- function(params.table,file.roots,sim.exts,site.exts,my.path){
    recover()
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
    probDerTable <- list()
    for (i in 1:nrow(params.table) ) {
        myNe <- format(round(params.table[i,'Ne'],3),3)
        myL <- floor(params.table[i,'target.size'])
        myrho <- format(round(params.table[i,'rho'],5),nsmall=5)
        mythr <- format(round(params.table[i,'thr'],2),nsmall=2)
        mycost  <- format(round(params.table[i,'cost'],2),nsmall=2)
        myenvSD  <- format(round(params.table[i,'env.sd'],5),nsmall=5)
        temp_prefix = paste(my.path,"/PopSize", myNe, "_LiaSize", myL, "_thr", mythr, "_cost", mycost, "_envsd", myenvSD, "_all", sep="")
        sim.files <- sapply(sim.exts,function(X) paste(temp_prefix,X,sep="."))
        sim.results[i,] <- colMeans(sapply(sim.files,function(X) as.numeric(read.table(X)[[1]])),na.rm=TRUE)
        ## was causing memory issues on cluster so temporarilily deleted
        site.files <- sapply(site.exts,function(X) paste(temp_prefix,X,sep="."))
        for(j in 1:length(site.files)){
            my.file <- file(site.files[j],'r')
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
        ## this is inefficient as it currently loops through the risk and der freq files twice
        ## but don't think the cost is worth dealing with
        der.file <- file(site.files[grep('der',site.files)],'r')
        risk.file <- file(site.files[grep('risk',site.files)],'r')
        Ne <- as.numeric(myNe)
        possibleFreqs <- 1:(2*Ne-1)/(2*Ne)
        riskFreqCounts <- rep(0,2*Ne-1)
        riskIsDerCounts <- rep(0,2*Ne-1)
        while(TRUE){
            derLine <- readLines(der.file,n=1)
            riskLine <- readLines(risk.file,n=1)
            if(length(derLine)==0){
                close(der.file)
                close(risk.file)
                break
            }
            if(derLine=="NA"){
                next
            }
            derFreqs <- factor(strsplit(derLine,',')[[1]],levels=possibleFreqs)
            riskFreqs <- factor(strsplit(riskLine,',')[[1]],levels=possibleFreqs)
            riskIsDerFreqs <- riskFreqs[riskFreqs==derFreqs]
            riskFreqCounts <- riskFreqCounts + table(riskFreqs)
            riskIsDerCounts <- riskIsDerCounts + table(riskIsDerFreqs)
        }
        probDerTable[[i]] <- list(
            params.table[i,],
            rbind(riskFreqCounts,riskIsDerCounts,riskIsDerCounts/riskFreqCounts)
        )
        if(i %% 10 ==0 ) print(i)
    }
    results <- list(
        cbind(params.table,sim.results,site.results),
        probDerTable
    )
    return(results)
}

options(scipen=400)
if(interactive()){
    my.path <- "smallEffectInsensitivity/all"
    my.args <- c(
        "smallEffectInsensitivityParamTable.txt",
        sapply(dir(my.path),function(X) paste(my.path,X,sep="/")),
        my.path,
        "smallEffectInsensitivityResultsTable.Rdata",
        "smallEffectInsensitivityDerProbs.Rdata"
    )
} else {
    ## read in command line arguments
    my.args <- commandArgs(trailingOnly=TRUE)
}


input.table.filename <- head(my.args,1)
output.table.filename1 <- tail(my.args,2)[1]
output.table.filename2 <- tail(my.args,2)[2]
my.path <- head(tail(my.args,3),1)
sim.filenames <- tail(head(my.args,-3),-1)

## load the parameter table
param.table <- read.table(input.table.filename)

## get unique sims and all the file extensions
files <- unique(sapply(sim.filenames,function(X) strsplit(X,'_all.')[[1]][1]))
exts <-  unique(sapply(sim.filenames,function(X) strsplit(X,'_all.')[[1]][2]))

site.exts <- exts[grepl('Freq',exts) | grepl('siteVar',exts) ]
sim.exts  <- exts[!(exts %in% site.exts)]

out <- merge_into_paramtable(param.table,files,sim.exts,site.exts,my.path)
results.table <- out[[1]]
der.probs <- out[[2]]
save(results.table,file=output.table.filename1)
save(der.probs,file=output.table.filename2)


