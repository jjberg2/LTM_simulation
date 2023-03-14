import numpy as np 
import pandas as pd

params_table = pd.read_csv("costInsensitivityParamTable.txt", delim_whitespace=True)



## global parameter ( doesn't change) 
mu=1e-6
cyc = 100
sampleInt = 100

## simulation variable 
rep = list(np.arange(0,2))
liaSizes = np.array((params_table["target.size"]).astype(int))
rhos = np.round(np.array(params_table["rho"]),3)
N = np.array(params_table["Ne"].astype(int))
envSD = np.round(np.array(params_table["env.sd"]),3)
cost = np.array(params_table["cost"])

#envSD = np.repeat(0,len(N))

rule all:
  input: 
    expand(expand("CostInsensitive/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_envSD{envSD}_cost{cost}__rep{{rep}}.prev",zip, N=N, liaSizes=liaSizes, rhos=rhos, envSD=envSD, cost=cost), rep=rep),
    expand("CostInsensitive/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_envSD{envSD}_cost{cost}_all.h2", zip, N=N, liaSizes=liaSizes, rhos=rhos, envSD=envSD, cost=cost)
    

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
    "CostInsensitive/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_envSD{envSD}_cost{cost}__rep{rep}.prev",
    "CostInsensitive/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_envSD{envSD}_cost{cost}__rep{rep}.h2"
  shell:
    """set +u; slim  -d mu={params.mu} -d rho_input={wildcards.rhos} -d p={wildcards.N} -d liaSize={wildcards.liaSizes} -d f={wildcards.cost}  -d e={wildcards.envSD} -d cyc={params.cyc} -d sampleInt={params.sampleInt} -d rep={wildcards.rep} {input.slim_script} > CostInsensitive/PopSize{wildcards.N}_LiaSize{wildcards.liaSizes}_rho{wildcards.rhos}_envSD{wildcards.envSD}_cost{wildcards.cost}__rep{wildcards.rep}.temp; set -u; 
    set +u; rm CostInsensitive/PopSize{wildcards.N}_LiaSize{wildcards.liaSizes}_rho{wildcards.rhos}_envSD{wildcards.envSD}_cost{wildcards.cost}__rep{wildcards.rep}.temp; set -u; """


rule result_combined: 
   input: 
     h2=expand("CostInsensitive/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_envSD{{envSD}}_cost{{cost}}__rep{rep}.h2", rep=rep),
     prev=expand("CostInsensitive/PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_envSD{{envSD}}_cost{{cost}}__rep{rep}.prev", rep=rep)
   output:
     h2="CostInsensitive/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_envSD{envSD}_cost{cost}_all.h2",
     prev="CostInsensitive/PopSize{N}_LiaSize{liaSizes}_rho{rhos}_envSD{envSD}_cost{cost}_all.prev"
   shell: 
     """cat {input.h2} >> {output.h2}; cat {input.prev} >> {output.prev}""" 
 
