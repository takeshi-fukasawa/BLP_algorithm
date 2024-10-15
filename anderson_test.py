#iteration=pyblp_test.Iteration('anderson',{'w0':0.01,'line_search':None},new_delta_mapping=True)

#results_spectral_anderson = problem.solve(initial_sigma,
#  initial_pi,iteration=iteration,optimization=tighter_bfgs,method='1s')

############## Nevo #############
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

iteration=pyblp_test.Iteration('Anderson_acceleration',{'atol':1e-14,'scheme':2,'mem_size':5},new_delta_mapping=True)

results_Anderson_1 = problem.solve(initial_sigma,initial_pi,iteration=iteration,optimization=tighter_bfgs,method='1s')

comp_time_delta_Anderson=sum(sum(results_Anderson_1.cumulative_comp_time_solve_delta))

iteration=pyblp_test.Iteration('df-sane',{'ftol':0,'fatol':1e-14,'line_search':"no"},new_delta_mapping=True)

results_spectral_1 = problem.solve(initial_sigma,initial_pi,iteration=iteration,optimization=tighter_bfgs,method='1s')

comp_time_delta_spectral=sum(sum(results_spectral_1.cumulative_comp_time_solve_delta))
