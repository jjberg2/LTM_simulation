## set results directory  
##setwd("~/Documents/PhD/Project/Rotation3_Berg/JeremySimulation/prev_results/set4_nucleotide/")
##setwd("~/Documents/academics/liability-model/LTM_simulation/CostInsensitivity")
## load the parameter table for simulation or input the table 
#load("../param.table.Rdata")  
new.table <- read.table("CostInsensitivity/costInsensitivityParamTable.txt",header=T)


## the input is the parameter table used for simulation 
## rho is the all the rho range simulated 
merge_into_paramtable <- function(params_table, rho){
    #recover()
    for (i in 1:nrow(params_table))
    {
        #index = which(params_table$rho== as.numeric(r))
        #e = round(params_table$env.sd[index], 3) ## find the environmental variance in the table 
        myrho  <- round(params_table[i,'rho'],3)
        mycost  <- params_table[i,'cost']
        
        temp_prefix = paste("CostInsensitivity/PopSize5000_LiaSize100000_rho",myrho,"_cost", mycost, "_all", sep="")
        temp_h2 = read.table(paste(temp_prefix,".h2", sep="")) 
        temp_prev = read.table(paste(temp_prefix,".prev", sep="")) 
        params_table$sim_h2[i] = mean(temp_h2$V1)
        params_table$sim_h2_sd[i] = sd(temp_h2$V1)/ sqrt(dim(temp_h2)[1])
        params_table$prev_emp[i] = mean(temp_prev$V1)
        params_table$prev_sd[i] = sd(temp_prev$V1)/ sqrt(dim(temp_prev)[1])
    }
    results = as.data.frame(params_table)
    return(results)
}


rho <- unique(new.table$rho)


results = merge_into_paramtable(new.table, rho)

save(results,file='costInsensitivity/costInsensitivityResultsTable.Rdata')



## ggplot() + geom_line(data= results, aes(x=thr, y =prev, col="Theory")) + 
##   geom_line(data=results, aes(x=thr,y=prev_emp, col="Simulations")) + geom_point(data=results, aes(x=thr,y=prev_emp)) +
##   scale_x_continuous(trans='log10') + xlab('Threshold Position') + ylab('Prevalence')
## ggsave('../figures/poisPrev.pdf',width=8,height=6,units='in')


## ggplot() + geom_line(data= results, aes(x=rho, y =h2, col="h2"))  + geom_point(data=results, aes(x=rho,y=h2)) +
##   scale_x_continuous(trans='log10')
