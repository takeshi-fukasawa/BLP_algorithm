import pyblp_test
import numpy as np
import pandas as pd
import pickle
import csv
pyblp_test.options.verbose = True  ##output status updates?
pyblp_test.options.digits=2
pyblp_test.__version__

tight_tol=False
if tight_tol==True :
  output_path="C:/Users/fukas/Dropbox/BLP/pyblp/tight_tol/"
  optimization = pyblp_test.Optimization('l-bfgs-b', {'ftol': 0, 'gtol': 1e-8})
else:
  output_path="C:/Users/fukas/Dropbox/BLP/pyblp/"
  optimization = pyblp_test.Optimization('l-bfgs-b', {'ftol': 0, 'gtol': 1e-4})

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
## Spectral-1
#iteration=pyblp_test.Iteration('df-sane',{'ftol':0,'fatol':1e-14,'line_search':"no"},new_delta_mapping=True)

iteration=pyblp_test.Iteration('Anderson_acceleration',{'atol':1e-14,'scheme':2,'mem_size':5},new_delta_mapping=True)

results_Anderson_1 = problem.solve(initial_sigma,
    initial_pi,iteration=iteration,optimization=optimization,costs_bounds=(0.001, None),
    W_type='clustered', se_type='clustered',initial_update=True)

####
iteration=pyblp_test.Iteration('df-sane',{'ftol':0,'fatol':1e-14,'line_search':"no"},new_delta_mapping=True)

results_spectral_1 = problem.solve(initial_sigma,
    initial_pi,iteration=iteration,optimization=optimization,costs_bounds=(0.001, None),
    W_type='clustered', se_type='clustered',initial_update=True)
