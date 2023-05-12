import numpy as np
import pandas as pd


## global parameter ( doesn't change)
mu=1e-6
cyc = 200
sampleInt = 25
toyRun=0
if(toyRun==1):
    print("Warning: the toyRun flag is on!")


rep = list(np.arange(0,1))


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


rule makeSmallEffectParamTable:
  input:
    script="scripts/makeSmallEffectInsensitivityParamFile.R"
  params:
    time="36:00:00",
    partition="broadwl",
    mem="2Gb"
  output:
    paramTable= "smallEffectInsensitivityParamTable.txt"
  shell:
    """Rscript {input.script}"""






## read large parameter tables
input_table_filename_large = "largeEffectInsensitivityParamTable.txt"
output_table_filename_large = "largeEffectInsensitivityResultsTable.Rdata"
params_table_large = pd.read_csv(input_table_filename_large, delim_whitespace=True)


## simulation variable
liaSizesLarge = np.array((params_table_large["target.size"]).astype(int))
costLarge = np.round((params_table_large["cost"]),3).astype(str)
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
     mem="2Gb",
     path="largeEffectInsensitivity/all"
  output:
    protected(output_table_filename_large)
  shell:
    """Rscript scripts/collateSingleEffectResults.R {input} {params.path} {output}"""
   





## read small parameter tables
input_table_filename_small = "smallEffectInsensitivityParamTable.txt"
output_table_filename_small = "smallEffectInsensitivityResultsTable.Rdata"
params_table_small = pd.read_csv(input_table_filename_small, delim_whitespace=True)


## simulation variable
liaSizesSmall = np.array((params_table_small["target.size"]).astype(int))
costSmall = np.round((params_table_small["cost"]),3).astype(str)
NSmall = np.array(params_table_small["Ne"].astype(int))
## thr = np.array(params_table_small["thr"].astype(int))
envsdSmall = np.round(np.array(params_table_small["env.sd"]),3).astype(str)
tmpRhosSmall = np.array(params_table_small["rho"])
rhosSmall = np.array(['{:.5f}'.format(r) for r in tmpRhosSmall], dtype=np.str)
    
rule allSmallEffect:
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
     mem="2Gb",
     path="smallEffectInsensitivity/all"
  output:
    protected(output_table_filename_small)
  shell:
    """Rscript scripts/collateSingleEffectResults.R {input} {params.path} {output}"""
   


    

    

    
rule slim_simulate_withsegregating:
  input:
    slim_script="LTM_prev_nucleotide.slim",
    paramTable="{path}ParamTable.txt"
  params:
    mu=mu,
    cyc=cyc,
    sampleInt = sampleInt,
    toyRun = toyRun,
    time="36:00:00",
    partition="broadwl",
    mem="2Gb"
  log:
    "logs/{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.log"
  output:
    fixed="{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.fixed",
    mean="{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.mean",
    h2="{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.h2",
    prev="{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.prev",
    genVar="{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.genVar",
    nSeg="{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.nSeg",
    deltaR="{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.deltaR",
    riskFreq="{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.riskFreq",
    derFreq="{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.derFreq",
    siteVar="{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.siteVar",
    tmp="{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.tmp"
  shell:
    """thr=`awk 'BEGIN {{print {wildcards.liaSizes}*2*{wildcards.rhos}}}'`;
    set +u; slim  -d mu={params.mu} -d  rho_input={wildcards.rhos} -d p={wildcards.N} -d liaSize={wildcards.liaSizes} -d f={wildcards.cost}  -d e={wildcards.envsd} -d cyc={params.cyc} -d sampleInt={params.sampleInt} -d rep={wildcards.rep} -d "fixedOut='{output.fixed}'" -d "meanOut='{output.mean}'" -d "h2Out='{output.h2}'" -d "prevOut='{output.prev}'" -d "genVarOut='{output.genVar}'" -d "nSegOut='{output.nSeg}'" -d "deltaROut='{output.deltaR}'" -d "riskFreqOut='{output.riskFreq}'" -d "derFreqOut='{output.derFreq}'" -d "siteVarOut='{output.siteVar}'" -d toyRun={params.toyRun} {input.slim_script} > {output.tmp}; set -u"""


    
rule result_combined:
   input:
     fixed=expand("{{path}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.fixed", rep=rep),
     mean=expand("{{path}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.mean", rep=rep),
     h2=expand("{{path}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.h2", rep=rep),
     prev=expand("{{path}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.prev", rep=rep),
     genVar=expand("{{path}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.genVar", rep=rep),
     nSeg=expand("{{path}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.nSeg", rep=rep),
     deltaR=expand("{{path}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.deltaR", rep=rep),
     riskFreq=expand("{{path}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.riskFreq", rep=rep),
     derFreq=expand("{{path}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.derFreq", rep=rep),
     siteVar=expand("{{path}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.siteVar", rep=rep),
   params:
     time="36:00:00",
     partition="broadwl",
     mem="2Gb"
   log:
     "{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.log"
   output:
     fixed="{path}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.fixed",
     mean="{path}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.mean",
     h2="{path}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.h2",
     prev="{path}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.prev",
     genVar="{path}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.genVar",
     nSeg="{path}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.nSeg",
     deltaR="{path}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.deltaR",
     riskFreq="{path}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.riskFreq",
     derFreq="{path}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.derFreq",
     siteVar="{path}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.siteVar"          
   shell:
     """cat {input.fixed} >> {output.fixed}; cat {input.mean} >> {output.mean}; cat {input.h2} >> {output.h2}; cat {input.prev} >> {output.prev}; cat {input.genVar} >> {output.genVar}; cat {input.nSeg} >> {output.nSeg}; cat {input.deltaR} >> {output.deltaR}; cat {input.riskFreq} >> {output.riskFreq}; cat {input.derFreq} >> {output.derFreq}; cat {input.siteVar} >> {output.siteVar}"""
