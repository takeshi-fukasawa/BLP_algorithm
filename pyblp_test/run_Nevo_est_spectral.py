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

tighter_bfgs = pyblp_test.Optimization('bfgs', {'gtol': 1e-5})
#nevo_results = problem.solve(
#    initial_sigma,initial_pi,optimization=tighter_bfgs,method='1s'
#)

#####################
#########################
## Spectral-1-S1
iteration=pyblp_test.Iteration('df-sane',{'ftol':0,'fatol':1e-14,'line_search':"no",'scheme':1},new_delta_mapping=True)

results_spectral_1_S1 = problem.solve(initial_sigma,
  initial_pi,iteration=iteration,optimization=tighter_bfgs,method='1s')

with open(output_path+'Nevo_est_results_spectral_1_S1.pickle', mode='wb') as fo:
  pickle.dump(results_spectral_1_S1, fo)

  ## Spectral-1-S2
iteration=pyblp_test.Iteration('df-sane',{'ftol':0,'fatol':1e-14,'line_search':"no",'scheme':2},new_delta_mapping=True)

results_spectral_1_S2 = problem.solve(initial_sigma,
  initial_pi,iteration=iteration,optimization=tighter_bfgs,method='1s')

with open(output_path+'Nevo_est_results_spectral_1_S2.pickle', mode='wb') as fo:
  pickle.dump(results_spectral_1_S2, fo)

####
## Spectral-1-S1-global
iteration=pyblp_test.Iteration('df-sane',{'ftol':0,'fatol':1e-14,'line_search':"cruz",'scheme':1},new_delta_mapping=True)

results_spectral_1_S1_global = problem.solve(initial_sigma,
  initial_pi,iteration=iteration,optimization=tighter_bfgs,method='1s')

with open(output_path+'Nevo_est_results_spectral_1_S1_global.pickle', mode='wb') as fo:
  pickle.dump(results_spectral_1_S1_global, fo)

  ## Spectral-1-S2
iteration=pyblp_test.Iteration('df-sane',{'ftol':0,'fatol':1e-14,'line_search':"cruz",'scheme':2},new_delta_mapping=True)

results_spectral_1_S2_global = problem.solve(initial_sigma,
  initial_pi,iteration=iteration,optimization=tighter_bfgs,method='1s')

with open(output_path+'Nevo_est_results_spectral_1_S2_global.pickle', mode='wb') as fo:
  pickle.dump(results_spectral_1_S2_global, fo)

## Spectral-1-S3
iteration=pyblp_test.Iteration('df-sane',{'ftol':0,'fatol':1e-14,'line_search':"cruz",'scheme':3},new_delta_mapping=True)

results_spectral_1_S3_global = problem.solve(initial_sigma,
  initial_pi,iteration=iteration,optimization=tighter_bfgs,method='1s')

with open(output_path+'Nevo_est_results_spectral_1_S3_global.pickle', mode='wb') as fo:
  pickle.dump(results_spectral_1_S3_global, fo)


## SQUAREM-1
iteration=pyblp_test.Iteration('squarem',{'atol':1e-14},new_delta_mapping=True,scheme=1)

results_squarem_1_S1 = problem.solve(initial_sigma,
    initial_pi,iteration=iteration,optimization=tighter_bfgs,method='1s')

with open(output_path+'Nevo_est_results_squarem_1_S1.pickle', mode='wb') as fo:
  pickle.dump(results_squarem_1_S1, fo)

iteration=pyblp_test.Iteration('squarem',{'atol':1e-14},new_delta_mapping=True,scheme=2)

results_squarem_1_S2 = problem.solve(initial_sigma,
    initial_pi,iteration=iteration,optimization=tighter_bfgs,method='1s')

with open(output_path+'Nevo_est_results_squarem_1_S2.pickle', mode='wb') as fo:
  pickle.dump(results_squarem_1_S2, fo)

####################
with open(output_path+'Nevo_est_results_spectral_1_S1.pickle', mode='br') as fi:
  results_spectral_1_S1 = pickle.load(fi)

with open(output_path+'Nevo_est_results_spectral_1_S2.pickle', mode='br') as fi:
  results_spectral_1_S2 = pickle.load(fi)

with open(output_path+'Nevo_est_results_spectral_1.pickle', mode='br') as fi:
  results_spectral_1 = pickle.load(fi)

with open(output_path+'Nevo_est_results_spectral_1_S1_global.pickle', mode='br') as fi:
  results_spectral_1_S1_global = pickle.load(fi)

with open(output_path+'Nevo_est_results_spectral_1_S2_global.pickle', mode='br') as fi:
  results_spectral_1_S2_global = pickle.load(fi)

with open(output_path+'Nevo_est_results_spectral_1_S3_global.pickle', mode='br') as fi:
  results_spectral_1_S3_global = pickle.load(fi)

with open(output_path+'Nevo_est_results_squarem_1_S1.pickle', mode='br') as fi:
  results_squarem_1_S1 = pickle.load(fi)

with open(output_path+'Nevo_est_results_squarem_1_S2.pickle', mode='br') as fi:
  results_squarem_1_S2 = pickle.load(fi)

with open(output_path+'Nevo_est_results_squarem_1.pickle', mode='br') as fi:
  results_squarem_1 = pickle.load(fi)

import results_functions

results_table=np.array([results_functions.results_func(results_spectral_1_S1),
         results_functions.results_func(results_spectral_1_S2),
         results_functions.results_func(results_spectral_1),
         results_functions.results_func(results_spectral_1_S1_global),
         results_functions.results_func(results_spectral_1_S2_global),
         results_functions.results_func(results_spectral_1_S3_global),
         results_functions.results_func(results_squarem_1_S1),
         results_functions.results_func(results_squarem_1_S2),
         results_functions.results_func(results_squarem_1)
         ])

np.set_printoptions(suppress=True)
results_table=np.round(results_table,3)


with open(output_path+'Nevo_est_results_spectral.csv', 'w',newline="") as f:
    writer = csv.writer(f)
    writer.writerows(results_table)