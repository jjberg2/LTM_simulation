import numpy as np
import pandas as pd

params_table = pd.read_csv("largeEffectInsensitivityParamTable.txt", delim_whitespace=True)


## global parameter ( doesn't change)
mu=1e-6
cyc = 200
sampleInt = 50

# print(params_table)

## simulation variable
rep = list(np.arange(0,2))
#rhos = np.array(params_table["rho"])
liaSizes = np.array((params_table["target.size"]).astype(int))
cost = np.round((params_table["cost"]),3).astype(str)
rhos = np.round(np.array(params_table["rho"]),3).astype(str)
N = np.array(params_table["Ne"].astype(int))
thr = np.array(params_table["thr"].astype(int))
envsd = np.round(np.array(params_table["env.sd"]),3).astype(str)

print(envsd)

rule all:
  input:
    expand("largeEffectInsensitivity/all/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.prev",zip, N=N, liaSizes=liaSizes, rhos=rhos, cost=cost, envsd=envsd)
    #expand("PopSize{N}_LiaSize{liaSizes}_rho{rhos}_all.h2", zip, N=N, liaSizes=liaSizes, rhos=rhos)

rule slim_simulate_withsegregating:
  input:
    slim_script="LTM_prev_nucleotide_largeEffectOutput.slim"
  params:
    mu=mu,
    cyc=cyc,
    sampleInt = sampleInt,
    time="36:00:00",
    partition="broadwl",
    mem="4Gb"
  log:
    "logs/{prefix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.log"
  output:
    "{prefix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.prev",
    "{prefix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.h2"
  shell:
    """thr=`awk 'BEGIN {{print 1e5*2*{wildcards.rhos}}}'`;
    set +u; slim  -d mu={params.mu} -d -d "rho_input='{wildcards.rhos}'" -d p={wildcards.N} -d liaSize={wildcards.liaSizes} -d f={wildcards.cost}  -d "e='{wildcards.envsd}'" -d cyc={params.cyc} -d sampleInt={params.sampleInt} -d rep={wildcards.rep} -d "pre='{wildcards.prefix}'" {input.slim_script} > {wildcards.prefix}/PopSize{wildcards.N}_LiaSize{wildcards.liaSizes}_rho{wildcards.rhos}_cost{wildcards.cost}_envsd{wildcards.envsd}_rep{wildcards.rep}.temp; set -u;
    set +u; rm {wildcards.prefix}/PopSize{wildcards.N}_LiaSize{wildcards.liaSizes}_rho{wildcards.rhos}_cost{wildcards.cost}_envsd{wildcards.envsd}_rep{wildcards.rep}.temp; set -u"""

rule result_combined:
   input:
     h2=expand("{{prefix}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.h2", rep=rep),
     prev=expand("{{prefix}}/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.prev", rep=rep)
   params:
     time="36:00:00",
     partition="broadwl",
     mem="4Gb"
   log:
     "{prefix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.log"
   output:
     h2="{prefix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.h2",
     prev="{prefix}/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.prev"
   shell:
     """cat {input.h2} >> {output.h2}; cat {input.prev} >> {output.prev}"""
