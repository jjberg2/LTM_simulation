import numpy as np 
import pandas as pd

params_table = pd.read_csv("costInsensitivityParamTable.txt", delim_whitespace=True)


## global parameter ( doesn't change) 
mu=1e-6
cyc = 100
sampleInt = 50

print(params_table)

## simulation variable 
rep = list(np.arange(0,1))
rhos = np.array(params_table["rho"])
liaSizes = np.array((params_table["target.size"]).astype(int))
cost = np.round((params_table["cost"]),3)
rhos = np.round(np.array(params_table["rho"]),3)
N = np.array(params_table["Ne"].astype(int))

print(cost)

rule all:
  input: 
    expand(expand("CostInsensitivity/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_rep{{rep}}.prev",zip, N=N, liaSizes=liaSizes, rhos=rhos, cost=cost), rep=rep),
    #expand("PopSize{N}_LiaSize{liaSizes}_rho{rhos}_all.h2", zip, N=N, liaSizes=liaSizes, rhos=rhos)

rule slim_simulate_withsegregating:
  input:
    slim_script="LTM_prev_nucleotide.slim"
  params:
    mu=mu,
    cyc=cyc, 
    sampleInt = sampleInt, 
    time="50:00:00",
##    partition="jjberg",
    mem="5Gb"
  output:
    "CostInsensitivity/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_rep{rep}.prev"
  shell:
    """thr=`awk 'BEGIN {{print 1e5*2*{wildcards.rhos}}}'`;
    env=`Rscript --vanilla scripts/getEnvSD.R 1e-6 ${{thr}} 0.9 5000 1e5 0.5 | awk '{{print $2}}'`;
    set +u; slim  -d mu={params.mu} -d rho_input={wildcards.rhos} -d p={wildcards.N} -d liaSize={wildcards.liaSizes} -d f={wildcards.cost}  -d e=${{env}} -d cyc={params.cyc} -d sampleInt={params.sampleInt} -d rep={wildcards.rep} {input.slim_script} > PopSize{wildcards.N}_LiaSize{wildcards.liaSizes}_rho{wildcards.rhos}_rep{wildcards.rep}.temp; set -u;""" 
    #set +u; rm PopSize{wildcards.N}_LiaSize{wildcards.liaSizes}_rho{wildcards.rhos}_rep{wildcards.rep}.temp; set -u"""

rule result_combined: 
   input: 
     h2=expand("CostInsensitivity/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_rep{rep}.h2", rep=rep),
     prev=expand("CostInsensitivity/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_rep{rep}.prev", rep=rep)
   output:
     h2="CostInsensitivity/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_all.h2",
     prev="CostInsensitivity/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_all.prev"
   shell: 
     """cat {input.h2} >> {output.h2}; cat {input.prev} >> {output.prev}""" 
 
