import sys
print (sys.version)
print (sys.path)

import numpy as np
import pandas as pd


## global parameter ( doesn't change)
smallCyc = 200
largeCyc = 1600
mu=1e-6
sampleInt = 25
toyRun=0
if(toyRun==1):
    print("Warning: the toyRun flag is on!")


rep = list(np.arange(0,4))




## make parameter tables
rule makeLargeEffectParamTable:
  input:
    params="scripts/largeEffectInsensitivityParams.R",
    script="scripts/makeLargeEffectInsensitivityParamFile.R"
  params:
    time="36:00:00",
    partition="broadwl",
    mem="2Gb"
  output:
    paramTable= "largeEffectInsensitivityParamTable.txt"
  shell:
    """Rscript {input.script}"""


# rule makeSmallEffectParamTable:
#   input:
#     script="scripts/makeSmallEffectInsensitivityParamFile.R"
#   params:
#     time="36:00:00",
#     partition="broadwl",
#     mem="2Gb"
#   output:
#     paramTable= "smallEffectInsensitivityParamTable.txt"
#   shell:
#     """Rscript {input.script}"""

rule makeSmallEffectParamTable:
  input:
    script="scripts/{name}ParamFile.R"
  params:
    time="36:00:00",
    partition="broadwl",
    mem="2Gb"
  output:
    paramTable= "{name}ParamTable.txt"
  shell:
    """Rscript {input.script}"""


paramTables = [
    "largeEffectInsensitivityParamTable.txt",
    "smallEffectInsensitivityParamTable.txt",
    "smallEffectVarianceInsensParamTable.txt"
]

    
rule makeAllParamTables:
  input:
     paramTables






##################################################
###### large effect cost insensitivity sims ######
##################################################     

## read large parameter tables
input_table_filename_large = "largeEffectInsensitivityParamTable.txt"
output_table_filename_large = "largeEffectInsensitivityResultsTable.Rdata"
params_table_large = pd.read_csv(input_table_filename_large, delim_whitespace=True)


## simulation variable
liaSizesLarge = np.array((params_table_large["target.size"]).astype(int))
tmpCostLarge = np.round((params_table_large["cost"]),2)
costLarge = np.array(['{:.2f}'.format(r) for r in tmpCostLarge], dtype=np.str)
NLarge = np.array(params_table_large["Ne"].astype(int))
## thr = np.array(params_table_large["thr"].astype(int))
tmpEnvSDsLarge = np.array(params_table_large["env.sd"])
envsdLarge = np.array(['{:.5f}'.format(r) for r in tmpEnvSDsLarge], dtype=np.str)
tmpRhosLarge = np.array(params_table_large["rho"])
rhosLarge = np.array(['{:.5f}'.format(r) for r in tmpRhosLarge], dtype=np.str)




rule allLargeEffect:
  input:
    input_table_filename_large,
    expand("largeEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.fixed",zip, N=NLarge, liaSizes=liaSizesLarge, rhos=rhosLarge, cost=costLarge, envsd=envsdLarge),
    expand("largeEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.mean",zip, N=NLarge, liaSizes=liaSizesLarge, rhos=rhosLarge, cost=costLarge, envsd=envsdLarge),
    expand("largeEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.h2",zip, N=NLarge, liaSizes=liaSizesLarge, rhos=rhosLarge, cost=costLarge, envsd=envsdLarge),    
    expand("largeEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.prev",zip, N=NLarge, liaSizes=liaSizesLarge, rhos=rhosLarge, cost=costLarge, envsd=envsdLarge),
    expand("largeEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.genVar",zip, N=NLarge, liaSizes=liaSizesLarge, rhos=rhosLarge, cost=costLarge, envsd=envsdLarge),
    expand("largeEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.nSeg",zip, N=NLarge, liaSizes=liaSizesLarge, rhos=rhosLarge, cost=costLarge, envsd=envsdLarge),
    expand("largeEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.deltaR",zip, N=NLarge, liaSizes=liaSizesLarge, rhos=rhosLarge, cost=costLarge, envsd=envsdLarge),
    expand("largeEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.riskFreq",zip, N=NLarge, liaSizes=liaSizesLarge, rhos=rhosLarge, cost=costLarge, envsd=envsdLarge),
    expand("largeEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.derFreq",zip, N=NLarge, liaSizes=liaSizesLarge, rhos=rhosLarge, cost=costLarge, envsd=envsdLarge),
    expand("largeEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.siteVar",zip, N=NLarge, liaSizes=liaSizesLarge, rhos=rhosLarge, cost=costLarge, envsd=envsdLarge)
  params:
     time="36:00:00",
     partition="broadwl",
     mem="4Gb",
     path="largeEffectInsensitivity/all"
  output:
    protected(output_table_filename_large)
  shell:
    """Rscript scripts/collateSingleEffectResults.R {input} {params.path} {output}"""
   


    
##################################################
###### small effect cost insensitivity sims ######
##################################################

## read small parameter tables
input_table_filename_small = "smallEffectInsensitivityParamTable.txt"
output_table_filename_small = "smallEffectInsensitivityResultsTable.Rdata"
params_table_small = pd.read_csv(input_table_filename_small, delim_whitespace=True)


## simulation variable
liaSizesSmall = np.array((params_table_small["target.size"]).astype(int))
tmpCostSmall = np.round((params_table_small["cost"]),2)
costSmall = np.array(['{:.2f}'.format(r) for r in tmpCostSmall], dtype=np.str)
NSmall = np.array(params_table_small["Ne"].astype(int))
## thr = np.array(params_table_small["thr"].astype(int))
tmpEnvSDsSmall = np.array(params_table_small["env.sd"])
envsdSmall = np.array(['{:.5f}'.format(r) for r in tmpEnvSDsSmall], dtype=np.str)
tmpRhosSmall = np.array(params_table_small["rho"])
rhosSmall = np.array(['{:.5f}'.format(r) for r in tmpRhosSmall], dtype=np.str)
    
rule allSmallEffectCost:
  input:
    input_table_filename_small,
    expand("smallEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.fixed",zip, N=NSmall, liaSizes=liaSizesSmall, rhos=rhosSmall, cost=costSmall, envsd=envsdSmall),
    expand("smallEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.mean",zip, N=NSmall, liaSizes=liaSizesSmall, rhos=rhosSmall, cost=costSmall, envsd=envsdSmall),
    expand("smallEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.h2",zip, N=NSmall, liaSizes=liaSizesSmall, rhos=rhosSmall, cost=costSmall, envsd=envsdSmall),    
    expand("smallEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.prev",zip, N=NSmall, liaSizes=liaSizesSmall, rhos=rhosSmall, cost=costSmall, envsd=envsdSmall),
    expand("smallEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.genVar",zip, N=NSmall, liaSizes=liaSizesSmall, rhos=rhosSmall, cost=costSmall, envsd=envsdSmall),
    expand("smallEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.nSeg",zip, N=NSmall, liaSizes=liaSizesSmall, rhos=rhosSmall, cost=costSmall, envsd=envsdSmall),
    expand("smallEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.deltaR",zip, N=NSmall, liaSizes=liaSizesSmall, rhos=rhosSmall, cost=costSmall, envsd=envsdSmall),
    expand("smallEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.riskFreq",zip, N=NSmall, liaSizes=liaSizesSmall, rhos=rhosSmall, cost=costSmall, envsd=envsdSmall),
    expand("smallEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.derFreq",zip, N=NSmall, liaSizes=liaSizesSmall, rhos=rhosSmall, cost=costSmall, envsd=envsdSmall),
    expand("smallEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.siteVar",zip, N=NSmall, liaSizes=liaSizesSmall, rhos=rhosSmall, cost=costSmall, envsd=envsdSmall)
  params:
     time="36:00:00",
     partition="broadwl",
     mem="4Gb",
     path="smallEffectInsensitivity/all"
  output:
     "smallEffectInsensitivityResultsTable.Rdata"
  shell:
    """Rscript scripts/collateSingleEffectResults.R {input} {params.path} {output}"""
   


    
##################################################
#### small effect variance insensitivity sims ####
##################################################

## read small parameter tables
input_table_filename_smallVe = "smallEffectVarianceInsensParamTable.txt"
output_table_filename_smallVe = "smallEffectVarianceInsensResultsTable.Rdata"
params_table_smallVe = pd.read_csv(input_table_filename_smallVe, delim_whitespace=True)


## simulation variable
liaSizesSmallVe = np.array((params_table_smallVe["target.size"]).astype(int))
tmpCostSmallVe = np.round((params_table_smallVe["cost"]),2)
costSmallVe = np.array(['{:.2f}'.format(r) for r in tmpCostSmallVe], dtype=np.str)
NSmallVe = np.array(params_table_smallVe["Ne"].astype(int))
## thr = np.array(params_table_smallVe["thr"].astype(int))
tmpEnvSDsSmallVe = np.array(params_table_smallVe["env.sd"])
envsdSmallVe = np.array(['{:.5f}'.format(r) for r in tmpEnvSDsSmallVe], dtype=np.str)
tmpRhosSmallVe = np.array(params_table_smallVe["rho"])
rhosSmallVe = np.array(['{:.5f}'.format(r) for r in tmpRhosSmallVe], dtype=np.str)
    
rule allSmallEffectVariance:
  input:
    input_table_filename_smallVe,
    expand("smallEffectVarianceInsens/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.fixed",zip, N=NSmallVe, liaSizes=liaSizesSmallVe, rhos=rhosSmallVe, cost=costSmallVe, envsd=envsdSmallVe),
    expand("smallEffectVarianceInsens/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.mean",zip, N=NSmallVe, liaSizes=liaSizesSmallVe, rhos=rhosSmallVe, cost=costSmallVe, envsd=envsdSmallVe),
    expand("smallEffectVarianceInsens/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.h2",zip, N=NSmallVe, liaSizes=liaSizesSmallVe, rhos=rhosSmallVe, cost=costSmallVe, envsd=envsdSmallVe),    
    expand("smallEffectVarianceInsens/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.prev",zip, N=NSmallVe, liaSizes=liaSizesSmallVe, rhos=rhosSmallVe, cost=costSmallVe, envsd=envsdSmallVe),
    expand("smallEffectVarianceInsens/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.genVar",zip, N=NSmallVe, liaSizes=liaSizesSmallVe, rhos=rhosSmallVe, cost=costSmallVe, envsd=envsdSmallVe),
    expand("smallEffectVarianceInsens/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.nSeg",zip, N=NSmallVe, liaSizes=liaSizesSmallVe, rhos=rhosSmallVe, cost=costSmallVe, envsd=envsdSmallVe),
    expand("smallEffectVarianceInsens/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.deltaR",zip, N=NSmallVe, liaSizes=liaSizesSmallVe, rhos=rhosSmallVe, cost=costSmallVe, envsd=envsdSmallVe),
    expand("smallEffectVarianceInsens/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.riskFreq",zip, N=NSmallVe, liaSizes=liaSizesSmallVe, rhos=rhosSmallVe, cost=costSmallVe, envsd=envsdSmallVe),
    expand("smallEffectVarianceInsens/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.derFreq",zip, N=NSmallVe, liaSizes=liaSizesSmallVe, rhos=rhosSmallVe, cost=costSmallVe, envsd=envsdSmallVe),
    expand("smallEffectVarianceInsens/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.siteVar",zip, N=NSmallVe, liaSizes=liaSizesSmallVe, rhos=rhosSmallVe, cost=costSmallVe, envsd=envsdSmallVe)
  params:
     time="36:00:00",
     partition="broadwl",
     mem="4Gb",
     path="smallEffectVarianceInsens/all"
  output:
     "smallEffectVarianceInsensResultsTable.Rdata"
  shell:
    """Rscript scripts/collateSingleEffectResults.R {input} {params.path} {output}"""
   






rule allSmallEffect:
  input:
    "smallEffectInsensitivityResultsTable.Rdata",
    "smallEffectVarianceInsensResultsTable.Rdata"
    

    

    
rule slim_simulate_small:
  input:
    slim_script="LTM_prev_nucleotide.slim",
    paramTable="smallEffect{suffix}ParamTable.txt"
  params:
    mu=mu,
    cyc=smallCyc,
    sampleInt = sampleInt,
    toyRun = toyRun,
    time="36:00:00",
    partition="broadwl",
    mem="4Gb"
  log:
    "logs/smallEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.log"
  output:
    fixed="smallEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.fixed",
    mean="smallEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.mean",
    h2="smallEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.h2",
    prev="smallEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.prev",
    genVar="smallEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.genVar",
    nSeg="smallEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.nSeg",
    deltaR="smallEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.deltaR",
    riskFreq="smallEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.riskFreq",
    derFreq="smallEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.derFreq",
    siteVar="smallEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.siteVar",
    tmp="smallEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.tmp"
  shell:
    """thr=`awk 'BEGIN {{print {wildcards.liaSizes}*2*{wildcards.rhos}}}'`;
    set +u; slim  -d mu={params.mu} -d  rho_input={wildcards.rhos} -d p={wildcards.N} -d liaSize={wildcards.liaSizes} -d f={wildcards.cost}  -d e={wildcards.envsd} -d cyc={params.cyc} -d sampleInt={params.sampleInt} -d rep={wildcards.rep} -d "fixedOut='{output.fixed}'" -d "meanOut='{output.mean}'" -d "h2Out='{output.h2}'" -d "prevOut='{output.prev}'" -d "genVarOut='{output.genVar}'" -d "nSegOut='{output.nSeg}'" -d "deltaROut='{output.deltaR}'" -d "riskFreqOut='{output.riskFreq}'" -d "derFreqOut='{output.derFreq}'" -d "siteVarOut='{output.siteVar}'" -d toyRun={params.toyRun} {input.slim_script} > {output.tmp}; set -u"""




    
rule slim_simulate_large:
  input:
    slim_script="LTM_prev_nucleotide.slim",
    paramTable="largeEffect{suffix}ParamTable.txt"
  params:
    mu=mu,
    cyc=largeCyc,
    sampleInt = sampleInt,
    toyRun = toyRun,
    time="36:00:00",
    partition="broadwl",
    mem="4Gb"
  log:
    "logs/largeEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.log"
  output:
    fixed="largeEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.fixed",
    mean="largeEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.mean",
    h2="largeEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.h2",
    prev="largeEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.prev",
    genVar="largeEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.genVar",
    nSeg="largeEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.nSeg",
    deltaR="largeEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.deltaR",
    riskFreq="largeEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.riskFreq",
    derFreq="largeEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.derFreq",
    siteVar="largeEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.siteVar",
    tmp="largeEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.tmp"
  shell:
    """thr=`awk 'BEGIN {{print {wildcards.liaSizes}*2*{wildcards.rhos}}}'`;
    set +u; slim  -d mu={params.mu} -d  rho_input={wildcards.rhos} -d p={wildcards.N} -d liaSize={wildcards.liaSizes} -d f={wildcards.cost}  -d e={wildcards.envsd} -d cyc={params.cyc} -d sampleInt={params.sampleInt} -d rep={wildcards.rep} -d "fixedOut='{output.fixed}'" -d "meanOut='{output.mean}'" -d "h2Out='{output.h2}'" -d "prevOut='{output.prev}'" -d "genVarOut='{output.genVar}'" -d "nSegOut='{output.nSeg}'" -d "deltaROut='{output.deltaR}'" -d "riskFreqOut='{output.riskFreq}'" -d "derFreqOut='{output.derFreq}'" -d "siteVarOut='{output.siteVar}'" -d toyRun={params.toyRun} {input.slim_script} > {output.tmp}; set -u"""








    
    
rule result_combined_small:
   input:
     fixed=expand("smallEffect{{suffix}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.fixed", rep=rep),
     mean=expand("smallEffect{{suffix}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.mean", rep=rep),
     h2=expand("smallEffect{{suffix}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.h2", rep=rep),
     prev=expand("smallEffect{{suffix}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.prev", rep=rep),
     genVar=expand("smallEffect{{suffix}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.genVar", rep=rep),
     nSeg=expand("smallEffect{{suffix}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.nSeg", rep=rep),
     deltaR=expand("smallEffect{{suffix}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.deltaR", rep=rep),
     riskFreq=expand("smallEffect{{suffix}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.riskFreq", rep=rep),
     derFreq=expand("smallEffect{{suffix}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.derFreq", rep=rep),
     siteVar=expand("smallEffect{{suffix}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.siteVar", rep=rep),
   params:
     time="36:00:00",
     partition="broadwl",
     mem="2Gb"
   log:
     "smallEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.log"
   output:
     fixed="smallEffect{suffix}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.fixed",
     mean="smallEffect{suffix}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.mean",
     h2="smallEffect{suffix}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.h2",
     prev="smallEffect{suffix}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.prev",
     genVar="smallEffect{suffix}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.genVar",
     nSeg="smallEffect{suffix}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.nSeg",
     deltaR="smallEffect{suffix}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.deltaR",
     riskFreq="smallEffect{suffix}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.riskFreq",
     derFreq="smallEffect{suffix}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.derFreq",
     siteVar="smallEffect{suffix}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.siteVar"          
   shell:
     """cat {input.fixed} >> {output.fixed}; cat {input.mean} >> {output.mean}; cat {input.h2} >> {output.h2}; cat {input.prev} >> {output.prev}; cat {input.genVar} >> {output.genVar}; cat {input.nSeg} >> {output.nSeg}; cat {input.deltaR} >> {output.deltaR}; cat {input.riskFreq} >> {output.riskFreq}; cat {input.derFreq} >> {output.derFreq}; cat {input.siteVar} >> {output.siteVar}"""






     

    
    
rule result_combined_large:
   input:
     fixed=expand("largeEffect{suffix}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.fixed", rep=rep),
     mean=expand("largeEffect{suffix}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.mean", rep=rep),
     h2=expand("largeEffect{suffix}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.h2", rep=rep),
     prev=expand("largeEffect{suffix}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.prev", rep=rep),
     genVar=expand("largeEffect{suffix}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.genVar", rep=rep),
     nSeg=expand("largeEffect{suffix}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.nSeg", rep=rep),
     deltaR=expand("largeEffect{suffix}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.deltaR", rep=rep),
     riskFreq=expand("largeEffect{suffix}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.riskFreq", rep=rep),
     derFreq=expand("largeEffect{suffix}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.derFreq", rep=rep),
     siteVar=expand("largeEffect{suffix}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.siteVar", rep=rep),
   params:
     time="36:00:00",
     partition="broadwl",
     mem="2Gb"
   log:
     "largeEffect{suffix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.log"
   output:
     fixed="largeEffect{suffix}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.fixed",
     mean="largeEffect{suffix}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.mean",
     h2="largeEffect{suffix}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.h2",
     prev="largeEffect{suffix}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.prev",
     genVar="largeEffect{suffix}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.genVar",
     nSeg="largeEffect{suffix}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.nSeg",
     deltaR="largeEffect{suffix}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.deltaR",
     riskFreq="largeEffect{suffix}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.riskFreq",
     derFreq="largeEffect{suffix}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.derFreq",
     siteVar="largeEffect{suffix}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.siteVar"          
   shell:
     """cat {input.fixed} >> {output.fixed}; cat {input.mean} >> {output.mean}; cat {input.h2} >> {output.h2}; cat {input.prev} >> {output.prev}; cat {input.genVar} >> {output.genVar}; cat {input.nSeg} >> {output.nSeg}; cat {input.deltaR} >> {output.deltaR}; cat {input.riskFreq} >> {output.riskFreq}; cat {input.derFreq} >> {output.derFreq}; cat {input.siteVar} >> {output.siteVar}"""






