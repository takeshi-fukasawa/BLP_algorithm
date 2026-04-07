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

## Differentiation IV

#formulation=pyblp_test.Formulation('1 + sugar + meshy')
#local_instruments = pyblp_test.build_differentiation_instruments(
#  formulation,
#  product_data
#)

#quadratic_instruments = pyblp_test.#build_differentiation_instruments(
#  formulation,
#  product_data,
#  version='quadratic'
#)

#for i in range(8):
#  del product_data[f'demand_instruments{i}']

#for i, column in enumerate(quadratic_instruments.T):
#  product_data[f'demand_instruments{100+i}'] = column

#for i, column in enumerate(local_instruments.T):
#  product_data[f'demand_instruments{200+i}'] = column


#product_data_dict=product_data

#product_data_dict = {k: product_data[k] for k in ['market_ids', 'city_ids','quarter','product_ids','firm_ids','brand_ids','shares', 'prices', 'sugar', 'mushy']}
#product_data_dict['demand_instruments'] = local_instruments

problem = pyblp_test.Problem(product_formulations, product_data,agent_formulation,agent_data)

initial_sigma = np.diag([0.3302, 2.4526, 0.0163, 0.2441])
initial_pi = np.array([
    [ 5.4819, 0, 0.2037, 0 ],
    [15.8935, -1.2000, 0, 2.6342],
    [-0.2506, 0, 0.0511, 0 ],
    [ 1.2650, 0, -0.8091, 0 ]
])

opt_setting = pyblp_test.Optimization('l-bfgs-b', {'gtol': 1e-4,'ftol':1e-4})


## Anderson-1
iteration=pyblp_test.Iteration('Anderson_acceleration',{'atol':1e-14,'scheme':2,'mem_size':5},new_delta_mapping=True)

#import Anderson_acceleration_functions
#iteration=pyblp_test.Iteration(Anderson_acceleration_functions.Anderson_acceleration_iterator,{'max_evaluations':5000,'norm':Anderson_acceleration_functions.infinity_norm,
#    'atol':1e-14,'rtol':0,'scheme':2,'mem_size':5},new_delta_mapping=True)


results_Anderson_1 = problem.solve(initial_sigma,
                                   initial_pi,iteration=iteration,optimization=opt_setting,method='1s',finite_differences=False,check_optimality='gradient')

### SQUAREM-0
#iteration=pyblp_test.Iteration('squarem',{'atol':1e-14},new_delta_mapping=False)

#results_squarem_0 = problem.solve(initial_sigma,
#    initial_pi,iteration=iteration,optimization=opt_setting,finite_differences=False,method='1s',check_optimality='gradient')

