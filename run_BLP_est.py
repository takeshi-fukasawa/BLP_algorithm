import pyblp_test
import numpy as np
import pandas as pd
import pickle
import csv
pyblp_test.options.verbose = True  ##output status updates?
pyblp_test.options.digits=2
pyblp_test.__version__

tight_tol=False
numerical_diff=False

if tight_tol==True :
  output_path="C:/Users/fukas/Dropbox/BLP/pyblp/tight_tol/"
  optimization = pyblp_test.Optimization('l-bfgs-b', {'ftol': 0, 'gtol': 1e-8})
else:
  output_path="C:/Users/fukas/Dropbox/BLP/pyblp/"
  optimization = pyblp_test.Optimization('l-bfgs-b', {'ftol': 0, 'gtol': 1e-4})

if numerical_diff==True :
  output_path="C:/Users/fukas/Dropbox/BLP/pyblp/"
  optimization = pyblp_test.Optimization('l-bfgs-b', {'ftol': 0, 'gtol': 1e-4},'compute_gradient',False)
  


product_data=pd.read_csv(pyblp_test.data.BLP_PRODUCTS_LOCATION)
agent_data=pd.read_csv(pyblp_test.data.BLP_AGENTS_LOCATION)


product_formulations = (
pyblp_test.Formulation('1 + hpwt + air + mpd + space'),
pyblp_test.Formulation('1 + prices + hpwt + air + mpd + space'),
pyblp_test.Formulation('1 + log(hpwt) + air + log(mpg) + log(space) + trend')
)

agent_formulation = pyblp_test.Formulation('0 + I(1 / income)')

problem = pyblp_test.Problem(product_formulations, product_data, agent_formulation, agent_data, costs_type='log')

initial_sigma = np.diag([3.612, 0, 4.628, 1.818, 1.050, 2.056])
initial_pi = np.c_[[0, -43.501, 0, 0, 0, 0]]

iteration=pyblp_test.Iteration('anderson')###### slow??? ###
iteration=pyblp_test.Iteration('squarem',new_delta_mapping=False)###### slower than squarem?
iteration=pyblp_test.Iteration('simple',new_delta_mapping=False)######
#iteration=pyblp_test.Iteration('lm',new_delta_mapping=False)######
iteration=pyblp_test.Iteration('df-sane',new_delta_mapping=False)


#########################
## Anderson-1
iteration=pyblp_test.Iteration('Anderson_acceleration',{'atol':1e-14,'scheme':2,'mem_size':5},new_delta_mapping=True)

results_Anderson_1 = problem.solve(initial_sigma,
    initial_pi,iteration=iteration,optimization=optimization,costs_bounds=(0.001, None),
    W_type='clustered', se_type='clustered',initial_update=True)

with open(output_path+'BLP_est_results_Anderson_1.pickle', mode='wb') as fo:
  pickle.dump(results_Anderson_1, fo)

## Spectral-1
iteration=pyblp_test.Iteration('df-sane',{'ftol':0,'fatol':1e-14,'line_search':"no"},new_delta_mapping=True)

results_spectral_1 = problem.solve(initial_sigma,
    initial_pi,iteration=iteration,optimization=optimization,costs_bounds=(0.001, None),
    W_type='clustered', se_type='clustered',initial_update=True)

with open(output_path+'BLP_est_results_spectral_1.pickle', mode='wb') as fo:
  pickle.dump(results_spectral_1, fo)

## SQUAREM-1
iteration=pyblp_test.Iteration('squarem',{'atol':1e-14},new_delta_mapping=True)

results_squarem_1 = problem.solve(initial_sigma,
    initial_pi,iteration=iteration,optimization=optimization,costs_bounds=(0.001, None),
    W_type='clustered', se_type='clustered',initial_update=True)

with open(output_path+'BLP_est_results_squarem_1.pickle', mode='wb') as fo:
  pickle.dump(results_squarem_1, fo)

## Simple-1
iteration=pyblp_test.Iteration('simple',{'atol':1e-14},new_delta_mapping=True)


results_simple_1 = problem.solve(initial_sigma,
    initial_pi,iteration=iteration,optimization=optimization,costs_bounds=(0.001, None),
    W_type='clustered', se_type='clustered',initial_update=True)

with open(output_path+'BLP_est_results_simple_1.pickle', mode='wb') as fo:
  pickle.dump(results_simple_1, fo)

########################
## Anderson-0
iteration=pyblp_test.Iteration('Anderson_acceleration',{'atol':1e-14,'scheme':2,'mem_size':5},new_delta_mapping=False)

results_Anderson_0 = problem.solve(initial_sigma,
    initial_pi,iteration=iteration,optimization=optimization,costs_bounds=(0.001, None),
    W_type='clustered', se_type='clustered',initial_update=True)

with open(output_path+'BLP_est_results_Anderson_0.pickle', mode='wb') as fo:
  pickle.dump(results_Anderson_0, fo)

## Spectral-0
iteration=pyblp_test.Iteration('df-sane',{'ftol':0,'fatol':1e-14,'line_search':"no"},new_delta_mapping=False)


results_spectral_0 = problem.solve(initial_sigma,
    initial_pi,iteration=iteration,optimization=optimization,costs_bounds=(0.001, None),
    W_type='clustered', se_type='clustered',initial_update=True)

with open(output_path+'BLP_est_results_spectral_0.pickle', mode='wb') as fo:
  pickle.dump(results_spectral_0, fo)

## SQUAREM-0
iteration=pyblp_test.Iteration('squarem',{'atol':1e-14},new_delta_mapping=False)


results_squarem_0 = problem.solve(initial_sigma,
    initial_pi,iteration=iteration,optimization=optimization,costs_bounds=(0.001, None),
    W_type='clustered', se_type='clustered',initial_update=True)

with open(output_path+'BLP_est_results_squarem_0.pickle', mode='wb') as fo:
  pickle.dump(results_squarem_0, fo)

## Simple-0
iteration=pyblp_test.Iteration('simple',{'atol':1e-14},new_delta_mapping=False)

results_simple_0 = problem.solve(initial_sigma,
    initial_pi,iteration=iteration,optimization=optimization,costs_bounds=(0.001, None),
    W_type='clustered', se_type='clustered',initial_update=True)

with open(output_path+'BLP_est_results_simple_0.pickle', mode='wb') as fo:
  pickle.dump(results_simple_0, fo)

####################
with open(output_path+'BLP_est_results_Anderson_1.pickle', mode='br') as fi:
  results_Anderson_1 = pickle.load(fi)

with open(output_path+'BLP_est_results_spectral_1.pickle', mode='br') as fi:
  results_spectral_1 = pickle.load(fi)

with open(output_path+'BLP_est_results_squarem_1.pickle', mode='br') as fi:
  results_squarem_1 = pickle.load(fi)

with open(output_path+'BLP_est_results_simple_1.pickle', mode='br') as fi:
  results_simple_1 = pickle.load(fi)

with open(output_path+'BLP_est_results_Anderson_0.pickle', mode='br') as fi:
  results_Anderson_0 = pickle.load(fi)

with open(output_path+'BLP_est_results_spectral_0.pickle', mode='br') as fi:
  results_spectral_0 = pickle.load(fi)

with open(output_path+'BLP_est_results_squarem_0.pickle', mode='br') as fi:
  results_squarem_0 = pickle.load(fi)

with open(output_path+'BLP_est_results_simple_0.pickle', mode='br') as fi:
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
         results_functions.results_func(results_simple_0)])

np.set_printoptions(suppress=True)
results_table=np.round(results_table,3)

with open(output_path+'BLP_est_results.csv', 'w',newline="") as f:
    writer = csv.writer(f)
    writer.writerows(results_table)

## Convergence properties
import scipy.io as sio
choice_prob=results_Anderson_1.compute_probabilities()
            
sio.savemat(output_path+'data_BLP.mat', {'market_ids':np.array(product_data.market_ids),
                             'shares':np.array(product_data.shares),'weights':np.array(agent_data.weights),'choice_prob':np.array(choice_prob)})


