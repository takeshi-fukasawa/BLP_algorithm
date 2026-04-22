import numpy as np
import csv
output_path="C:/Users/fukas/Dropbox/BLP/pyblp/"

#exec(open("test_nevo.py").read())
#exec(open("test_BLP.py").read())

import pandas as pd


n_trial=3
for trial_id in range(0,n_trial):
    #exec(open("test_nevo.py").read())
    exec(open("run_Nevo_est.py").read())
    exec(open("run_BLP_est.py").read())


#########
## Summarise Nevo results over trials
n_trial=3
n_param_Nevo=13
nevo_results_list = [
    np.loadtxt(f'{output_path}Nevo_est_results_trial_{i}.csv', delimiter=',')
    for i in range(n_trial)
]
Nevo_est_results_summary = sum(nevo_results_list) / n_trial
Nevo_est_results_summary=np.round(Nevo_est_results_summary[:,[0,2,5,6,3]],2)
Nevo_est_results_summary[8:11,[1,2]]=Nevo_est_results_summary[8:11,[1,2]]*(2*n_param_Nevo+1)# Because PyBLP itself does report only the contraction evaluations / objective evaluations excluding the ones used to compute numerical derivatives, we adjust those values here by multiplying by (2*n_param_BLP+1).


with open(output_path+'Nevo_est_results_summary.csv', 'w',newline="") as f:
    writer = csv.writer(f)
    writer.writerows(Nevo_est_results_summary)

results = pd.read_csv(output_path + "Nevo_est_results_summary.csv", header=None)
results = results.rename(columns={4: "comp_time", 1: "obj_eval"})

nonnumerical_time_per_obs = 0.098# cf. code profiling (profile_results_Nevo_numer_diff.csv)
results['comp_time_numerical'] = results['comp_time'] - (nonnumerical_time_per_obs * results['obj_eval'])


results.to_csv(output_path + "Nevo_est_results_summary_subtract_nonnumerical.csv", index=False)

#########
## Summarise BLP results over trials
n_param_BLP=6
blp_results_list = [
    np.loadtxt(f'{output_path}BLP_est_results_trial_{i}.csv', delimiter=',')
    for i in range(n_trial)
]
BLP_est_results_summary = sum(blp_results_list) / n_trial
BLP_est_results_summary=np.round(BLP_est_results_summary[:,[0,2,5,6,3]],2)
BLP_est_results_summary[8:10,[1,2]]=BLP_est_results_summary[8:10,[1,2]]*(2*n_param_BLP+1)# Because PyBLP itself does report only the contraction evaluations / objective evaluations excluding the ones used to compute numerical derivatives, we adjust those values here by multiplying by (2*n_param_BLP+1).

with open(output_path+'BLP_est_results_summary.csv', 'w',newline="") as f:
    writer = csv.writer(f)
    writer.writerows(BLP_est_results_summary)
