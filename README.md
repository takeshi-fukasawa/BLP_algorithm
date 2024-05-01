This repository contains replication codes of Fukasawa (2024), investigating computationally efficient inner-loop algorithms for estimating static/dynamic BLP models. 
It is composed of two folders. The first one is "static_dynamic_BLP_MATLAB", and it contains replication codes of the paper using MATLAB codes. The second one is pyblp_test, and it slightly modifies the original pyBLP package (Version 0.13.0; Conlon and Gortmaker (2020)) so that we can introduce the new mapping $\delta$-(1) and flexible specifications regarding the spectral algorithm.

## static_dynamic_BLP_MATLAB
We can obtain the static BLP Monte Carlo simulation results by running the file "XXX".
Regarding dynamic BLP, we can obtain the dynamic BLP Monte Carlo simulation results by running the file "".
Appendix A.2. of the paper shows numerical results under the case with large consumer heterogeneity. We can obtain the results by runnning XXX.

## pyblp_test
Codes in the folder build on the PyBLP package (version 0.13.0). Differences with the original codes are:
• The package name is pyblp-test, rather than pyblp.
• Introduction of new mapping $\delta$-(1):
By specifying new_delta_mapping_mapping=TRUE, we can use the new mapping $\delta$-(1). If new_delta_mapping_mapping=FALSE, we can use the standard BLP contraction mapping. We should specify it in pyblp_test.Iteration.
I modified XXX so as to introduce the new mapping.
• Spectral algorithm
The original PyBLP package relies on the dfsane function in SciPy package, if we want to use the spectral algorithm. Nevertheless, the function only allows the step size S2, and does not allow the algorithm without globalization strategies. Hence, I introduced an alternative for the original dfsane function, which is in configurations/spectral.py. In  the new function, we can specify different step sizes (S1,S2,S3). For instance, if we want to specify the step size S1, we should specify 
XXX in pyblp_test.
The default value of "scheme" is 3. Also, if we want to try the spectral algorithm without globalization steps, we should specify 'line_search':"no".

We can obtain the results of Fukasawa (2024), by running the following codes:
*XXX

## References
• Fukasawa
• Conlon
