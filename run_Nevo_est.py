import pyblp_test
import numpy as np
import pandas as pd
import pickle
import csv
pyblp_test.options.verbose = True  ##output status updates?
pyblp_test.options.digits=2
pyblp_test.__version__

output_path="C:/Users/fukas/Dropbox/BLP/pyblp/"

product_data=pd.read_csv(pyblp_test.data.NEVO_PRODUCTS_LOCATION)
agent_data = pd.read_csv(pyblp_test.data.NEVO_AGENTS_LOCATION)


## Random coefficient
X1_formulation = pyblp_test.Formulation('0 + prices',absorb='C(product_ids)')
#X1_formulation = pyblp_test.Formulation('0 + prices',absorb='C(market_ids)')
X2_formulation = pyblp_test.Formulation('1 + prices + sugar + mushy')
product_formulations = (X1_formulation, X2_formulation)

agent_formulation = pyblp_test.Formulation('0 + income + income_squared + age + child')

problem = pyblp_test.Problem(product_formulations, product_data,agent_formulation,agent_data)

initial_sigma = np.diag([0.3302, 2.4526, 0.0163, 0.2441])
initial_pi = np.array([
    [ 5.4819, 0, 0.2037, 0 ],
    [15.8935, -1.2000, 0, 2.6342],
    [-0.2506, 0, 0.0511, 0 ],
    [ 1.2650, 0, -0.8091, 0 ]
])

opt_setting = pyblp_test.Optimization('l-bfgs-b', {'gtol': 1e-4,'ftol':1e-4})


#####################

## Anderson-1
iteration=pyblp_test.Iteration('Anderson_acceleration',{'atol':1e-14,'scheme':2,'mem_size':5},new_delta_mapping=True)

results_Anderson_1 = problem.solve(initial_sigma,
  initial_pi,iteration=iteration,optimization=opt_setting,method='1s',finite_differences=False,check_optimality='gradient')

with open(output_path+'Nevo_est_results_Anderson_1_trial_'+str(trial_id)+'.pickle', mode='wb') as fo:
  pickle.dump(results_Anderson_1, fo)

## Anderson-1-mumer-diff
iteration=pyblp_test.Iteration('Anderson_acceleration',{'atol':1e-14,'scheme':2,'mem_size':5},new_delta_mapping=True)

results_Anderson_1_numer_diff = problem.solve(initial_sigma,
  initial_pi,iteration=iteration,optimization=opt_setting,method='1s',finite_differences=True,check_optimality='gradient')

with open(output_path+'Nevo_est_results_Anderson_1_numer_diff_trial_'+str(trial_id)+'.pickle', mode='wb') as fo:
  pickle.dump(results_Anderson_1_numer_diff, fo)

## Spectral-1
iteration=pyblp_test.Iteration('df-sane',{'ftol':0,'fatol':1e-14,'line_search':"no"},new_delta_mapping=True)

results_spectral_1 = problem.solve(initial_sigma,
  initial_pi,iteration=iteration,optimization=opt_setting,method='1s',finite_differences=False,check_optimality='gradient')

with open(output_path+'Nevo_est_results_spectral_1_trial_'+str(trial_id)+'.pickle', mode='wb') as fo:
  pickle.dump(results_spectral_1, fo)

## SQUAREM-1
iteration=pyblp_test.Iteration('squarem',{'atol':1e-14},new_delta_mapping=True)

results_squarem_1 = problem.solve(initial_sigma,
    initial_pi,iteration=iteration,optimization=opt_setting,method='1s',finite_differences=False,check_optimality='gradient')

with open(output_path+'Nevo_est_results_squarem_1_trial_'+str(trial_id)+'.pickle', mode='wb') as fo:
  pickle.dump(results_squarem_1, fo)

## Simple-1
iteration=pyblp_test.Iteration('simple',{'atol':1e-14},new_delta_mapping=True)

results_simple_1 = problem.solve(initial_sigma,
                                 initial_pi,iteration=iteration,optimization=opt_setting,method='1s',finite_differences=False,check_optimality='gradient')

with open(output_path+'Nevo_est_results_simple_1_trial_'+str(trial_id)+'.pickle', mode='wb') as fo:
  pickle.dump(results_simple_1, fo)

## Simple-1-numer-diff
iteration=pyblp_test.Iteration('simple',{'atol':1e-14},new_delta_mapping=True)

results_simple_1_numer_diff = problem.solve(initial_sigma,
                                 initial_pi,iteration=iteration,optimization=opt_setting,method='1s',finite_differences=True,check_optimality='gradient')

with open(output_path+'Nevo_est_results_simple_1_numer_diff_trial_'+str(trial_id)+'.pickle', mode='wb') as fo:
  pickle.dump(results_simple_1_numer_diff, fo)

########################
## Anderson-0
iteration=pyblp_test.Iteration('Anderson_acceleration',{'atol':1e-14,'scheme':2,'mem_size':5},new_delta_mapping=False)

results_Anderson_0 = problem.solve(initial_sigma,
  initial_pi,iteration=iteration,optimization=opt_setting,method='1s',finite_differences=False,check_optimality='gradient')

with open(output_path+'Nevo_est_results_Anderson_0_trial_'+str(trial_id)+'.pickle', mode='wb') as fo:
  pickle.dump(results_Anderson_0, fo)

## Spectral-0
iteration=pyblp_test.Iteration('df-sane',{'ftol':0,'fatol':1e-14,'line_search':"no"},new_delta_mapping=False)

results_spectral_0 = problem.solve(initial_sigma,
    initial_pi,iteration=iteration,optimization=opt_setting,method='1s',finite_differences=False,check_optimality='gradient')

with open(output_path+'Nevo_est_results_spectral_0_trial_'+str(trial_id)+'.pickle', mode='wb') as fo:
  pickle.dump(results_spectral_0, fo)

## SQUAREM-0
iteration=pyblp_test.Iteration('squarem',{'atol':1e-14},new_delta_mapping=False)

results_squarem_0 = problem.solve(initial_sigma,
    initial_pi,iteration=iteration,optimization=opt_setting,method='1s',finite_differences=False,check_optimality='gradient')

with open(output_path+'Nevo_est_results_squarem_0_trial_'+str(trial_id)+'.pickle', mode='wb') as fo:
  pickle.dump(results_squarem_0, fo)

## SQUAREM-0-numer-diff
iteration=pyblp_test.Iteration('squarem',{'atol':1e-14},new_delta_mapping=False)

results_squarem_0_numer_diff = problem.solve(initial_sigma,
    initial_pi,iteration=iteration,optimization=opt_setting,method='1s',finite_differences=True,check_optimality='gradient')

with open(output_path+'Nevo_est_results_squarem_0_numer_diff_trial_'+str(trial_id)+'.pickle', mode='wb') as fo:
  pickle.dump(results_squarem_0_numer_diff, fo)

## Simple-0
iteration=pyblp_test.Iteration('simple',{'atol':1e-14},new_delta_mapping=False)

results_simple_0 = problem.solve(initial_sigma,
                                 initial_pi,iteration=iteration,optimization=opt_setting,method='1s',finite_differences=False,check_optimality='gradient')

with open(output_path+'Nevo_est_results_simple_0_trial_'+str(trial_id)+'.pickle', mode='wb') as fo:
  pickle.dump(results_simple_0, fo)

## LM
#iteration=pyblp_test.Iteration('lm',{'ftol':1e-14},new_delta_mapping=False)
#results_LM_0 = problem.solve(initial_sigma,
#    initial_pi,iteration=iteration,optimization=opt_setting,method='1s')

#with open(output_path+'Nevo_est_results_LM_0.pickle', mode='wb') as fo:
#  pickle.dump(results_LM_0, fo)

####################
with open(output_path+'Nevo_est_results_Anderson_1_trial_'+str(trial_id)+'.pickle', mode='br') as fi:
  results_Anderson_1 = pickle.load(fi)

with open(output_path+'Nevo_est_results_Anderson_1_numer_diff_trial_'+str(trial_id)+'.pickle', mode='br') as fi:
  results_Anderson_1_numer_diff = pickle.load(fi)

with open(output_path+'Nevo_est_results_spectral_1_trial_'+str(trial_id)+'.pickle', mode='br') as fi:
  results_spectral_1 = pickle.load(fi)

with open(output_path+'Nevo_est_results_squarem_1_trial_'+str(trial_id)+'.pickle', mode='br') as fi:
  results_squarem_1 = pickle.load(fi)


with open(output_path+'Nevo_est_results_Anderson_0_trial_'+str(trial_id)+'.pickle', mode='br') as fi:
  results_Anderson_0 = pickle.load(fi)

with open(output_path+'Nevo_est_results_spectral_0_trial_'+str(trial_id)+'.pickle', mode='br') as fi:
  results_spectral_0 = pickle.load(fi)

with open(output_path+'Nevo_est_results_squarem_0_trial_'+str(trial_id)+'.pickle', mode='br') as fi:
  results_squarem_0 = pickle.load(fi)

with open(output_path+'Nevo_est_results_squarem_0_numer_diff_trial_'+str(trial_id)+'.pickle', mode='br') as fi:
  results_squarem_0_numer_diff = pickle.load(fi)

with open(output_path+'Nevo_est_results_simple_1_trial_'+str(trial_id)+'.pickle', mode='br') as fi:
  results_simple_1 = pickle.load(fi)

with open(output_path+'Nevo_est_results_simple_1_numer_diff_trial_'+str(trial_id)+'.pickle', mode='br') as fi:
  results_simple_1_numer_diff = pickle.load(fi)

with open(output_path+'Nevo_est_results_simple_0_trial_'+str(trial_id)+'.pickle', mode='br') as fi:
  results_simple_0 = pickle.load(fi)

import results_functions

results_table=np.array([
  results_functions.results_func(results_Anderson_1),
  results_functions.results_func(results_spectral_1),
  results_functions.results_func(results_squarem_1),
  results_functions.results_func(results_simple_1),
  results_functions.results_func(results_Anderson_0),
  results_functions.results_func(results_spectral_0),
  results_functions.results_func(results_squarem_0),
  results_functions.results_func(results_simple_0),
  results_functions.results_func(results_simple_1_numer_diff),
  results_functions.results_func(results_Anderson_1_numer_diff),
  results_functions.results_func(results_squarem_0_numer_diff)])

np.set_printoptions(suppress=True)
results_table=np.round(results_table,3)

with open(output_path+'Nevo_est_results_trial_'+str(trial_id)+'.csv', 'w',newline="") as f:
    writer = csv.writer(f)
    writer.writerows(results_table)

## Convergence properties
import scipy.io as sio
choice_prob=results_Anderson_1.compute_probabilities()

sio.savemat(output_path+'data_Nevo.mat', {'market_ids':np.array(product_data.market_ids),
                             'shares':np.array(product_data.shares),'weights':np.array(agent_data.weights),'choice_prob':np.array(choice_prob)})
 