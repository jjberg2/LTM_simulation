import numpy as np 
import pandas as pd

def trim_trailingzero(x):
    return(x.astype(str).strip("0").strip("."))

params_table = pd.read_csv("twoEffectPrevInsensitivityParamsTable.txt", delim_whitespace=True)
#params_table = params_table.iloc[[0,49,99,149,199],:].reset_index(drop=True)

## global parameter ( doesn't change) 
cyc = 200
sampleInt = 50
reps = 1
toyRun = 0
if(toyRun==1):
    print("Warning: the toyRun flag is on!")



    

    
## simulation variable 
rep = list(np.arange(0,reps))
rhot = np.round(np.array(params_table["rhot"]),3)
thr = np.round(np.array(params_table["thr"]),3)
cost = np.round(np.array(params_table["C"]),3)
alphaLarge = np.round(np.array(params_table["al"]),3)
#alphaLarge = [trim_trailingzero(x) for x in alphaLarge]
# print(cost)

N = np.array(params_table["Ne"].astype(int))
envSD = np.round(np.sqrt(params_table["Ve"]),3)

rule all:
  input: 
    expand(expand("twoEffectInsensitivity/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{{rep}}.prev",zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost), rep=rep),
    expand("twoEffectInsensitivity/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.prev", zip, N=N, alphaLarge = alphaLarge, thr=thr, envSD=envSD, cost=cost)



    
def find_index(wildcards, col):                                                                                                                                          
    index = np.where( np.round(params_table["al"],3) == float(wildcards.alphaLarge))[0][0]             
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
    h2l="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.h2l",
    h2s="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.h2s",
    genVar="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.genVar",
    nSegSmall="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.nSegSmall",
    nSegLarge="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.nSegLarge",
    deltaRSmall="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.deltaRSmall",
    deltaRLarge="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.deltaRLarge",
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
    '"{prefix}/logs/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_rep{rep}.log",'
  shell:
    """set +u; slim -d mu={params.mu} -d thr={params.thr} -d rhos={params.rhos} -d p={wildcards.N}  -d f={params.fitCost}  -d e={wildcards.envSD} -d cyc={params.cyc} -d sampleInt={params.sampleInt} -d rep={wildcards.rep} -d aS={params.alphaSmall} -d aL={wildcards.alphaLarge} -d liaSmall={params.liaSmall} -d liaLarge={params.liaLarge} -d "fixedSmallOut='{output.fixedSmall}'" -d "fixedLargeOut='{output.fixedLarge}'" -d "meanSmallOut='{output.meanSmall}'" -d "meanLargeOut='{output.meanLarge}'" -d "meanOut='{output.mean}'" -d "h2Out='{output.h2}'" -d "prevOut='{output.prev}'" -d "h2lOut='{output.h2l}'" -d "h2sOut='{output.h2s}'" -d "genVarOut='{output.genVar}'" -d "nSegSmallOut='{output.nSegSmall}'"  -d "nSegLargeOut='{output.nSegLarge}'" -d "deltaRSmallOut='{output.deltaRSmall}'" -d "deltaRLargeOut='{output.deltaRLarge}'" -d toyRun={params.toyRun} {input.slim_script} > {output.tmp}; set -u; """


rule result_combined: 
   input:
     fixedSmall=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.fixedSmall", rep=rep),
     fixedLarge=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.fixedLarge", rep=rep),
     meanSmall=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.meanSmall", rep=rep),
     meanLarge=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.meanLarge", rep=rep),
     mean=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.mean", rep=rep),
     h2=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.h2", rep=rep),
     prev=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.prev", rep=rep),
     h2l=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.h2l", rep=rep),
     h2s=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.h2s", rep=rep),
     genVar=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.genVar", rep=rep),
     nSegSmall=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.nSegSmall", rep=rep),
     nSegLarge=expand("{{prefix}}/PopSize{{N}}_aL{{alphaLarge}}_thr{{thr}}_envSD{{envSD}}_cost{{cost}}_rep{rep}.nSegLarge", rep=rep)
   output:
     fixedSmall="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.fixedSmall",
     fixedLarge="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.fixedLarge",
     meanSmall="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.meanSmall",
     meanLarge="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.meanLarge",
     mean="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.mean",
     h2="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.h2",
     prev="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.prev",
     h2l="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.h2l",
     h2s="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.h2s",
     genVar="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.genVar",
     nSegSmall="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.nSegSmall",
     nSegLarge="{prefix}/PopSize{N}_aL{alphaLarge}_thr{thr}_envSD{envSD}_cost{cost}_all.nSegLarge"
   shell: 
     """cat {input.fixedSmall} >> {output.fixedSmall}; cat {input.fixedLarge} >> {output.fixedLarge}; cat {input.meanSmall} >> {output.meanSmall}; cat {input.meanLarge} >> {output.meanLarge}; cat {input.mean} >> {output.mean}; cat {input.h2} >> {output.h2}; cat {input.prev} >> {output.prev}; cat {input.h2l} >> {output.h2l}; cat {input.h2s} >> {output.h2s}; cat {input.genVar} >> {output.genVar}; cat {input.nSegSmall} >> {output.nSegSmall}; cat {input.nSegLarge} >> {output.nSegLarge}""" 
 
