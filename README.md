This repository contains replication codes of Fukasawa (2024), investigating computationally efficient inner-loop algorithms for estimating static/dynamic BLP models. 
It is composed of two folders. The first one is "static_dynamic_BLP_MATLAB", and it contains replication codes of the paper using MATLAB. The second one is pyblp_test, and it slightly modifies the original pyBLP package (Version 0.13.0; Conlon and Gortmaker (2020)) so that we can introduce the new mapping $\delta$-(1) and flexible specifications regarding the spectral algorithm.

## static_dynamic_BLP_MATLAB
We can obtain the static BLP Monte Carlo simulation results by running the file "run_static_RCNL_main.m".
Regarding dynamic BLP, we can obtain the dynamic BLP Monte Carlo simulation results by running the file "run_dynamic_BLP_main.m".
Appendix A.2. of the paper shows numerical results for the case with large consumer heterogeneity. We can obtain the results by running "run_RCNL_main_hetero_test.m".

## pyblp_test
Codes in the folder are built on the PyBLP package (version 0.13.0). Differences with the original codes are:
* The package name is pyblp-test, rather than pyblp.  
* Introduction of new mapping $\delta$-(1):   
By specifying new_delta_mapping_mapping=TRUE, we can use the new mapping $\delta$-(1). If new_delta_mapping_mapping=FALSE, we can use the standard BLP contraction mapping. We should specify it in pyblp_test.Iteration.
I mainly modified markets/market.py so as to introduce the new mapping.   

* Spectral algorithm  
The original PyBLP package relies on the dfsane function in SciPy package (version 1.12.0), if we want to use the spectral algorithm. Nevertheless, the function only allows the step size S2, and does not allow the algorithm without globalization strategies. Hence, I introduced an alternative for the original dfsane function, which is in configurations/spectral.py. In the new function, we can specify different step sizes (S1,S2,S3) by changing the value of "scheme". The default value is 3. Also, if we want to try the spectral algorithm without globalization steps, we should specify 'line_search':"no". For instance, if we want to try the algorithm using the new mapping combined with the spectral algorithm with step size S1, we should specify:  
```
iteration=pyblp_test.Iteration('df-sane',{'ftol':0,'fatol':1e-14,'line_search':"no",'scheme':1},new_delta_mapping=True)
```

We can obtain the results of Fukasawa (2024), by running the following codes:
* "run_Nevo_est.py": Estimation using Nevo (2001)'s dataset
* "run_BLP_est.py": Estimation using Berry et al. (1995,1999)'s dataset
* "run_Nevo_est_spectral.py": Estimation using Nevo (2001)'s dataset, under different specifications of the spectral algorithm

To run these codes, we should set the current directory to the folder where BLP repository is located.
We can run "run_Nevo_est.py" by 
```
python run_Nevo_est.py
```

## References
* Berry, A., Levinsohn, J., & Pakes, A. (1995). Automobile Prices in Market Equilibrium. Econometrica, 63(4), 841-890.
* Berry, S., Levinsohn, J., & Pakes, A. (1999). Voluntary export restraints on automobiles: Evaluating a trade policy. American Economic Review, 89(3), 400-431.  
* Conlon, C., & Gortmaker, J. (2020). Best practices for differentiated products demand estimation with pyblp. The RAND Journal of Economics, 51(4), 1108-1161.  
* Fukasawa, T. (2024). Fast and simple inner-loop algorithms of static/dynamic BLP estimations. arXiv preprint arXiv:2404.04494.  
* Nevo, A. (2001). Measuring market power in the ready-to-eat cereal industry. Econometrica, 69(2), 307-342.
