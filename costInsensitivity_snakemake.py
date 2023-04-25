import numpy as np 
import pandas as pd

params_table = pd.read_csv("CostInsensitivity/costInsensitivityParamTable.txt", delim_whitespace=True)


## global parameter ( doesn't change) 
mu=1e-6
cyc = 200
sampleInt = 50

#print(params_table)

## simulation variable 
rep = list(np.arange(0,3))
rhos = np.array(params_table["rho"])
liaSizes = np.array((params_table["target.size"]).astype(int))
cost = np.round((params_table["cost"]),1)
rhos = np.round(np.array(params_table["rho"]),3)
N = np.array(params_table["Ne"].astype(int))
h2 = np.round(params_table["h2"],1)
envsd = np.round(np.array(params_table["env.sd"]),3)
## heritability

print(cost)

#rule all:
#  input: 
#    expand("CostInsensitivity/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_all.prev",zip, N=N, liaSizes=liaSizes, rhos=rhos, cost=cost),

    
rule all:
  input: 
    expand("CostInsensitivity/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.prev",zip, N=N, liaSizes=liaSizes, rhos=rhos, cost=cost, envsd=envsd),
    expand("CostInsensitivity/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.h2",zip, N=N, liaSizes=liaSizes, rhos=rhos, cost=cost, envsd=envsd),
    "CostInsensitivity/costInsensitivityParamTable.txt"
  output:
    "CostInsensitivity/costInsensitivityResultsTable.Rdata"
  shell:
    """Rscript scripts/costInsensitivityResults.R"""


    
rule slim_simulate_withsegregating:
  input:
    slim_script="LTM_prev_nucleotide.slim"
  params:
    mu=mu,
    cyc=cyc, 
    sampleInt = sampleInt, 
    time="50:00:00",
    partition="broadwl",
    mem="5Gb"
  output:
    "CostInsensitivity/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.prev",
    "CostInsensitivity/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_rep{rep}.h2"
  shell:
    """thr=`awk 'BEGIN {{print 1e5*2*{wildcards.rhos}}}'`;
    set +u; slim -d mu={params.mu} -d rho_input={wildcards.rhos} -d p={wildcards.N} -d liaSize={wildcards.liaSizes} -d f={wildcards.cost}  -d e={wildcards.envsd} -d cyc={params.cyc} -d sampleInt={params.sampleInt} -d rep={wildcards.rep} {input.slim_script} > CostInsensitivity/PopSize{wildcards.N}_LiaSize{wildcards.liaSizes}_rho{wildcards.rhos}_rep{wildcards.rep}_cost{wildcards.cost}_envsd{wildcards.envsd}.temp; set -u;""" 
    #set +u; rm PopSize{wildcards.N}_LiaSize{wildcards.liaSizes}_rho{wildcards.rhos}_rep{wildcards.rep}.temp; set -u"""

rule result_combined: 
   input: 
     h2=expand("CostInsensitivity/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.h2", rep=rep),
     prev=expand("CostInsensitivity/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_cost{{cost}}_envsd{{envsd}}_rep{rep}.prev", rep=rep)
   output:
     h2="CostInsensitivity/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.h2",
     prev="CostInsensitivity/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_cost{cost}_envsd{envsd}_all.prev"
   shell: 
     """cat {input.h2} >> {output.h2}; cat {input.prev} >> {output.prev}""" 
 
