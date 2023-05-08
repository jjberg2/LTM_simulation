import numpy as np 
import pandas as pd

input_table_filename = "smallEffectInsensitivityParamTable.txt"
results_table_filename = "smallEffectInsensitivityResultsTable.Rdata"

params_table = pd.read_csv(input_table_filename, delim_whitespace=True)


## global parameter ( doesn't change) 
mu=1e-6
cyc = 200
sampleInt = 50
toyRun=0

if(toyRun==1):
    print("Warning: the toyRun flag is on!")

## simulation variable 
rep = list(np.arange(0,2))
liaSizes = np.array((params_table["target.size"]).astype(int))
cost = np.round((params_table["cost"]),3)
N = np.array(params_table["Ne"].astype(int))
## thr = np.array(params_table["thr"].astype(int))
envsd = np.round(np.array(params_table["env.sd"]),3)
tmpRhos = np.array(params_table["rho"])
rhos = np.array(['{:.5f}'.format(r) for r in tmpRhos], dtype=np.str)

rule all:
  input:
    input_table_filename,
    expand("smallEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.fixed",zip, N=N, liaSizes=liaSizes, rhos=rhos, cost=cost, envsd=envsd),
    expand("smallEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.mean",zip, N=N, liaSizes=liaSizes, rhos=rhos, cost=cost, envsd=envsd),
    expand("smallEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.h2",zip, N=N, liaSizes=liaSizes, rhos=rhos, cost=cost, envsd=envsd),    
    expand("smallEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.prev",zip, N=N, liaSizes=liaSizes, rhos=rhos, cost=cost, envsd=envsd),
    expand("smallEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.genVar",zip, N=N, liaSizes=liaSizes, rhos=rhos, cost=cost, envsd=envsd),
    expand("smallEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.nSeg",zip, N=N, liaSizes=liaSizes, rhos=rhos, cost=cost, envsd=envsd),
    expand("smallEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.deltaR",zip, N=N, liaSizes=liaSizes, rhos=rhos, cost=cost, envsd=envsd)
  params:
     time="36:00:00",
     partition="broadwl",
     mem="2Gb",
     path="smallEffectInsensitivity/all"
  output:
    protected(results_table_filename)
  shell:
    """Rscript scripts/collateSingleEffectResults.R {input} {params.path} {output}"""


rule slim_simulate_withsegregating:
  input:
    slim_script="LTM_prev_nucleotide.slim"
  params:
    mu=mu,
    cyc=cyc,
    sampleInt = sampleInt,
    toyRun = toyRun,
    time="36:00:00",
    partition="broadwl",
    mem="4Gb"
  log:
    "logs/{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.log"
  output:
    fixed=temp("{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.fixed"),
    mean=temp("{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.mean"),
    h2=temp("{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.h2"),
    prev=temp("{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.prev"),
    genVar=temp("{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.genVar"),
    nSeg=temp("{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.nSeg"),
    deltaR=temp("{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.deltaR"),
    tmp=temp("{path}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.tmp")
  shell:
    """thr=`awk 'BEGIN {{print {wildcards.liaSizes}*2*{wildcards.rhos}}}'`;
    set +u; slim  -d mu={params.mu} -d  rho_input={wildcards.rhos} -d p={wildcards.N} -d liaSize={wildcards.liaSizes} -d f={wildcards.cost}  -d e={wildcards.envsd} -d cyc={params.cyc} -d sampleInt={params.sampleInt} -d rep={wildcards.rep} -d "fixedOut='{output.fixed}'" -d "meanOut='{output.mean}'" -d "h2Out='{output.h2}'" -d "prevOut='{output.prev}'" -d "genVarOut='{output.genVar}'" -d "nSegOut='{output.nSeg}'" -d "deltaR='{output.deltaR}'" -d toyRun={params.toyRun} {input.slim_script} > {output.tmp}; set -u"""



rule result_combined:
   input:
     fixed=expand("{{path}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.fixed", rep=rep),
     mean=expand("{{path}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.mean", rep=rep),
     h2=expand("{{path}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.h2", rep=rep),
     prev=expand("{{path}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.prev", rep=rep),
     genVar=expand("{{path}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.genVar", rep=rep),
     nSeg=expand("{{path}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.nSeg", rep=rep),
     deltaR=expand("{{path}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.deltaR", rep=rep)
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
     deltaR="{path}/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.deltaR"          
   shell:
     """cat {input.fixed} >> {output.fixed}; cat {input.mean} >> {output.mean}; cat {input.h2} >> {output.h2}; cat {input.prev} >> {output.prev}; cat {input.genVar} >> {output.genVar}; cat {input.nSeg} >> {output.nSeg}; cat {input.deltaR} >> {output.deltaR}"""

