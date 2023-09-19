#! /usr/bin/python3
import sys
print(sys.version)
import numpy as np 
import pandas as pd
import subprocess as sp

def trim_trailingzero(x):
    return(x.astype(str).strip("0").strip("."))


#subprocess.call ("Rscript --vanilla scripts/makeTwoEffectPrevInsensitivityParamsTable.R", shell=True)
## sp.call ("Rscript --vanilla scripts/makeTwoEffectCostInsensitivityParamsTable.R", shell=True)


input_table_costInsens = "twoEffectCostInsensitivityN1000ParamsTable.txt"
params_table = pd.read_csv(input_table_costInsens, delim_whitespace=True)
#params_table = params_table.iloc[[0,49,99,149,199],:].reset_index(drop=True)

## global parameter ( doesn't change) 
cyc = 800
sampleInt = 50
reps = 40
toyRun = 0
if(toyRun==1):
    print("Warning: the toyRun flag is on!")



    

    
## simulation variable 
rep = list(np.arange(0,reps))
rhot = np.round(np.array(params_table["rhot"]),3)
tmpThr = np.round(np.array(params_table["thr"]),1)
thr = np.array(['{:.1f}'.format(r) for r in tmpThr], dtype=str)
tmpCost = np.round(np.array(params_table["C"]),2)
cost = np.array(['{:.2f}'.format(r) for r in tmpCost], dtype=str)
tmpAlphaLarge = np.round(np.array(params_table["al"]),3)
alphaLarge = np.array(['{:.3f}'.format(r) for r in tmpAlphaLarge], dtype=str)
#alphaLarge = [trim_trailingzero(x) for x in alphaLarge]
# print(cost)
tmpN = np.array(params_table["Ne"].astype(int))
N = np.array(['{:.0f}'.format(r) for r in tmpN], dtype=str)
tmpEnvSD = np.round(np.sqrt(params_table["Ve"]),3)
envSD = np.array(['{:.3f}'.format(r) for r in tmpEnvSD], dtype=str)

print(params_table)


rule make_results_table:
  input:
    input_table_costInsens,
    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.fixedSmall", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost),
    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.fixedLarge", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost),
    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.meanSmall", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost),
    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.meanLarge", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost),
    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.mean", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost),
    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.prev", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost),
    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.h2l", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost),
    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.h2s", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost),
#    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.h2ol", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost),
#    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.h2os", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost),
#    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.h2o", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost),        
    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.genVar", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost),
    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.nSegSmall", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost),
    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.nSegLarge", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost),
    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.deltaRSmall", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost),
    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.deltaRLarge", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost),
    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.riskFreqSmall", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost),
    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.riskFreqLarge", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost),
    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.derFreqSmall", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost),
    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.derFreqLarge", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost),
    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.siteVarSmall", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost),
    expand("{{prefix}}/all/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.siteVarLarge", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost)
  params:
     time="36:00:00",
     partition="broadwl",
     mem="4Gb",
     path="{prefix}/all"
  output:
     "{prefix}ResultsTable.Rdata"
  shell:
    """Rscript scripts/collateTwoEffectResults.R {input} {params.path} {output}"""



rule N500_results_table:
  input:
    "twoEffectCostInsensitivityN500ResultsTable.Rdata"



rule N1000_results_table:
  input:
    "twoEffectCostInsensitivityN1000ResultsTable.Rdata"


















    

    
rule twoEffectBackup:
  input:
    "{paramTableName}ResultsTable.Rdata"
  params:
    time="36:00:00",
    partition="broadwl",
    mem="2Gb"
  output:
    protected("{paramTableName}ResultsTableBackup.Rdata")
  shell:
    """cp {input} {output}"""    
    



    
def find_index(wildcards, col):                                                                                                                                          
    index = np.where( (np.round(params_table["al"],3) == float(wildcards.alphaLarge)) & (params_table["C"] == float(wildcards.cost)))[0][0]             
    return(params_table[col][index])




rule slim_simulate_withsegregating:
  input:
    slim_script="scripts/LTM_prev_nucleotide_multieffect_diffstart.slim"
  output:
    fixedSmall="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.fixedSmall",
    fixedLarge="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.fixedLarge",
    meanSmall="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.meanSmall",
    meanLarge="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.meanLarge",
    mean="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.mean",
    h2="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.h2",
    prev="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.prev",
    h2s="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.h2s",
    h2l="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.h2l",
    h2os="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.h2os",
    h2ol="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.h2ol",
    h2o="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.h2o",
    genVar="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.genVar",
    nSegSmall="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.nSegSmall",
    nSegLarge="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.nSegLarge",
    deltaRSmall="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.deltaRSmall",
    deltaRLarge="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.deltaRLarge",
    riskFreqSmall="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.riskFreqSmall",
    riskFreqLarge="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.riskFreqLarge",
    derFreqSmall="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.derFreqSmall",
    derFreqLarge="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.derFreqLarge",
    siteVarSmall="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.siteVarSmall",
    siteVarLarge="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.siteVarLarge",
    tmp="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.tmp"
  params:
    mu = lambda wildcards: find_index(wildcards, col="u"),
    fitCost= lambda wildcards: find_index(wildcards, col="C"),
    alphaSmall = lambda wildcards: find_index(wildcards, col="as"),
    liaSmall = lambda wildcards: find_index(wildcards, col="Ls"),
    liaLarge = lambda wildcards: find_index(wildcards, col="Ll"),
    rhos = lambda wildcards:find_index(wildcards, col="rhos"),
    thr = lambda wildcards:find_index(wildcards, col="thr"), 
    cyc=cyc,
    sampleInt = sampleInt,
    toyRun=toyRun,
    time="36:00:00",
    partition="broadwl",
    mem="4Gb"
  log:
    "{prefix}/logs/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.log"
  shell:
    """set +u; slim -d mu={params.mu} -d thr={params.thr} -d rhos={params.rhos} -d p={wildcards.N}  -d f={params.fitCost}  -d e={wildcards.envSD} -d cyc={params.cyc} -d sampleInt={params.sampleInt} -d rep={wildcards.rep} -d aS={params.alphaSmall} -d aL={wildcards.alphaLarge} -d liaSmall={params.liaSmall} -d liaLarge={params.liaLarge} -d "fixedSmallOut='{output.fixedSmall}'" -d "fixedLargeOut='{output.fixedLarge}'" -d "meanSmallOut='{output.meanSmall}'" -d "meanLargeOut='{output.meanLarge}'" -d "meanOut='{output.mean}'" -d "h2Out='{output.h2}'" -d "h2sOut='{output.h2s}'" -d "h2lOut='{output.h2l}'" -d "h2osOut='{output.h2os}'" -d "h2olOut='{output.h2ol}'" -d "h2oOut='{output.h2o}'" -d "prevOut='{output.prev}'" -d "genVarOut='{output.genVar}'" -d "nSegSmallOut='{output.nSegSmall}'"  -d "nSegLargeOut='{output.nSegLarge}'" -d "deltaRSmallOut='{output.deltaRSmall}'" -d "deltaRLargeOut='{output.deltaRLarge}'" -d "riskFreqSmallOut='{output.riskFreqSmall}'" -d "riskFreqLargeOut='{output.riskFreqLarge}'" -d "derFreqSmallOut='{output.derFreqSmall}'" -d "derFreqLargeOut='{output.derFreqLarge}'" -d "siteVarSmallOut='{output.siteVarSmall}'" -d "siteVarLargeOut='{output.siteVarLarge}'" -d toyRun={params.toyRun} {input.slim_script} > {output.tmp}; set -u; """



rule result_combined_two_small:
  input:
    expand("{{prefix}}/PopSize{{N}}_aL{{liaSizes}}_thr{{thr}}_envSD{{envsd}}_cost{{cost}}_rep{rep}.{{ext}}", rep=rep)
  params:
     time="36:00:00",
     partition="broadwl",
     mem="2Gb"
  log:
    "{prefix}/PopSize{N}_aL{liaSizes}_thr{thr}_envSD{envsd}_cost{cost}_all_{ext}.log"
  output:
    "{prefix}/all/PopSize{N}_aL{liaSizes}_thr{thr}_envSD{envsd}_cost{cost}_all.{ext}"
  shell:
     """cat {input} >> {output}"""




    

    
# rule result_combined: 
#    input:
#      fixedSmall=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.fixedSmall", rep=rep),
#      fixedLarge=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.fixedLarge", rep=rep),
#      meanSmall=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.meanSmall", rep=rep),
#      meanLarge=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.meanLarge", rep=rep),
#      mean=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.mean", rep=rep),
#      h2=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.h2", rep=rep),
#      prev=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.prev", rep=rep),
#      h2l=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.h2l", rep=rep),
#      h2s=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.h2s", rep=rep),
#      genVar=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.genVar", rep=rep),
#      nSegSmall=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.nSegSmall", rep=rep),
#      nSegLarge=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.nSegLarge", rep=rep)
#    output:
#      fixedSmall="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.fixedSmall",
#      fixedLarge="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.fixedLarge",
#      meanSmall="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.meanSmall",
#      meanLarge="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.meanLarge",
#      mean="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.mean",
#      h2="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.h2",
#      prev="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.prev",
#      h2l="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.h2l",
#      h2s="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.h2s",
#      genVar="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.genVar",
#      nSegSmall="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.nSegSmall",
#      nSegLarge="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.nSegLarge"
#    shell: 
#      """cat {input.fixedSmall} >> {output.fixedSmall}; cat {input.fixedLarge} >> {output.fixedLarge}; cat {input.meanSmall} >> {output.meanSmall}; cat {input.meanLarge} >> {output.meanLarge}; cat {input.mean} >> {output.mean}; cat {input.h2} >> {output.h2}; cat {input.prev} >> {output.prev}; cat {input.h2l} >> {output.h2l}; cat {input.h2s} >> {output.h2s}; cat {input.genVar} >> {output.genVar}; cat {input.nSegSmall} >> {output.nSegSmall}; cat {input.nSegLarge} >> {output.nSegLarge}""" 
 
