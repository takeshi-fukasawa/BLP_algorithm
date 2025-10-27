import numpy as np
import csv
output_path="C:/Users/fukas/Dropbox/BLP/pyblp/"

#exec(open("test_nevo.py").read())
#exec(open("test_BLP.py").read())

from turtle import pd


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
