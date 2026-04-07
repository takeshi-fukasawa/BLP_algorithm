# Replication code of "Fast and simple inner-loop algorithms of static / dynamic BLP estimations" (Fukasawa, 2024)
This repository contains replication codes of the paper titled "Fast and simple inner-loop algorithms of static / dynamic BLP estimations" by Takeshi Fukasawa, investigating computationally efficient inner-loop algorithms for estimating static/dynamic BLP models. 
It is composed of two folders. The first is "static_dynamic_BLP_algorithm", and it contains replication codes of the paper using MATLAB. The second is pyblp_test, and it slightly modifies the original pyBLP package (Version 0.13.0; Conlon and Gortmaker, 2020) so that we can introduce the new mapping $\delta$-(1) and flexible specifications regarding the spectral algorithm.

Note that the code relies on my other repository takeshi-fukasawa/spectral. If you run the replication code, please also download the spectral repository, and modify the path in static_dynamic_BLP_algorithm/preliminary_path.m.

## static_dynamic_BLP_algorithm
We can obtain the static BLP Monte Carlo simulation results by running the file "run_static_RCNL_main.m".
Regarding dynamic BLP, we can obtain the dynamic BLP Monte Carlo simulation results by running the file "run_dynamic_BLP_main.m".
Appendix A.2 of the paper shows numerical results for the case with large consumer heterogeneity. We can obtain the results by running "run_RCNL_main_hetero_test.m".

## pyblp_test
Codes in the folder are built on the PyBLP package (version 0.13.0). Differences with the original codes are:
* The package name is pyblp-test, rather than pyblp.  
* Introduction of a new mapping $\delta$-(1):   
By specifying new_delta_mapping=TRUE, we can use the new mapping $\delta$-(1). If new_delta_mapping=FALSE, we can use the standard BLP contraction mapping. We should specify it in pyblp_test.Iteration.
I mainly modified markets/market.py so as to introduce the new mapping.   

* Spectral algorithm  
The original PyBLP package relies on the dfsane function in SciPy package (version 1.12.0), if we want to use the spectral algorithm. Nevertheless, the function only allows the step size S2, and does not allow the algorithm without globalization strategies. Hence, I introduced an alternative for the original dfsane function, which is in pyblp_test/configurations/spectral.py. In the new function, we can specify different step sizes (S1,S2,S3) by changing the value of "scheme". The default value is 3 (S3). Also, if we want to try the spectral algorithm without globalization steps, we should specify `'line_search':"no"`. For instance, if we want to try the algorithm using the new mapping $\delta$-(1) combined with the spectral algorithm with step size S1, we should specify:  
```
iteration=pyblp_test.Iteration('df-sane',{'ftol':0,'fatol':1e-14,'line_search':"no",'scheme':1},new_delta_mapping=True)
```

* Anderson acceleration  
Because Anderson acceleration method for fixed point iterations is not directly available in PyBLP, I added a code on Anderson acceleration. For instance, if we want to try the algorithm using the new mapping $\delta$-(1) combined with the Anderson acceleration with memory size $m=5$ abd scheme 2 (Type-II Anderson acceleration), we should specify:  
```
iteration=pyblp_test.Iteration('Anderson_acceleration',{'atol':1e-14,'scheme':2,'mem_size':5},new_delta_mapping=True)
```

We can obtain the results of Fukasawa (2024), by running the following code:
* "run_Nevo_est.py": Estimation using Nevo (2001)'s dataset
* "run_BLP_est.py": Estimation using Berry et al. (1995,1999)'s dataset
* "run_Nevo_est_spectral.py": Estimation using Nevo (2001)'s dataset, under different specifications of the spectral algorithm

To run these codes, we should set the current directory to the folder where BLP repository is located.
We can run "run_Nevo_est.py" by 
```
python run_Nevo_est.py
```

## Use of Anderson acceleration in PyBLP
While the current study significantly modified the original PyBLP for numerical experiments, Anderson acceleration can be easily implemented in PyBLP without modifying the core source code by utilizing a custom Iteration configuration. The pseudo-code is provided below:

```
import Anderson_acceleration_functions
iteration=pyblp_test.Iteration(Anderson_acceleration_functions.Anderson_acceleration_iterator,{'max_evaluations':5000,'norm':Anderson_acceleration_functions.infinity_norm,'atol':1e-14,'rtol':0,'scheme':2,'mem_size':5},new_delta_mapping=True)

results = problem.solve(initial_sigma,initial_pi,iteration=iteration)

```
To run the code, save the following Python file (Anderson_acceleration_functions.py) in your working directory:

```
import numpy as np
def Anderson_acceleration_iterator(
        initial, contraction, iteration_callback, max_evaluations,
        atol, rtol, norm, scheme, mem_size):
    """Apply the Anderson acceleration method for fixed point iteration."""
    m = int(mem_size)
    N = initial.size  # dimension of the variable


    x = initial
    failed = False
    k = 0

    # Initialize matrices to store past values

    resid_past_mat = np.empty((N, max_evaluations))
    fun_past_mat = np.empty((N, max_evaluations))
    x_past_mat = np.empty((N, max_evaluations))

    while True:
        # first step
        x0, (x, weights) = x, contraction(x)[:2]
        if not all_finite(x, weights):
            x = x0
            failed = True
            break

        g0 = x - x0
        #print(max(abs(g0)))

        # check for convergence
        if k >= max_evaluations or termination_check(x, g0, weights, atol, rtol, norm):
            break

        resid_k_vec=x-x0
        x_k_vec=x0
        fun_k_vec=x

        resid_past_mat[:,k] = resid_k_vec
        fun_past_mat[:,k] = fun_k_vec
        x_past_mat[:,k] = x_k_vec

        if k >= 1:
            m_k = min(m, k)

            Z = resid_past_mat[:, k]
            
            DF = np.diff(resid_past_mat[:, k - m_k:k+1])#resid_past_mat(:, k - m_k+1:k+1) in MATLAB

            if scheme == 1:
                # Type I Anderson (Corresponding to Good Broyden update)

                DX = np.diff(x_past_mat[:, k - m_k :k+1])
                gamma = np.linalg.solve(DX.T @ DF, DX.T @ Z)# Singular matrix error??
                #gamma=np.linalg.lstsq(DX.T @ DF,DX.T @ Z,rcond=None)[0]

            else:  # scheme == 2
                # Type II Anderson (Corresponding to Bad Broyden update)
                gamma = np.linalg.solve(DF.T @ DF, DF.T @ Z) # Singular matrix error??
                #gamma=np.linalg.lstsq(DF.T @ DF,DF.T @ Z,rcond=None)[0]
            

            alpha_vec = np.empty(m_k + 1) # m_k = len(gamma)
            
            # alpha_0 = gamma_0
            alpha_vec[0] = gamma[0]
            
            # alpha_i = gamma_i - gamma_{i-1} for 1 <= i <= m_k - 1
            alpha_vec[1:-1] = gamma[1:] - gamma[:-1]
            
            # alpha_{m_k} = 1 - gamma_{m_k - 1}
            alpha_vec[-1] = 1 - gamma[-1]
            
            x = fun_past_mat[:, k - m_k:k + 1] @ alpha_vec # Update x

        # record the completion of a major iteration
        iteration_callback() ######

        k += 1

    # determine whether there was convergence
    converged = not failed and k < max_evaluations
    return x, k

def infinity_norm(x) -> float:
    """Compute the infinity norm of a vector."""
    return np.abs(x).max()

def all_finite(*arrays):
    """Validate that multiple arrays are either None or all finite."""
    return all(a is None or np.isfinite(a).all() for a in arrays)


def termination_check(
        x, residual, weights, atol, rtol,
        norm):
    """Check whether the residual indicates that iteration should be terminated."""
    tol = atol
    if rtol > 0:
        tol += rtol * norm(weight(x, weights))
    return norm(weight(residual, weights)) < tol

def weight(x, weights):
    """Optionally weight an array."""
    if weights is None:
        return x
    return weights * x

```

## References
* Berry, S., Levinsohn, J., & Pakes, A. (1995). Automobile Prices in Market Equilibrium. Econometrica, 63(4), 841-890.
* Berry, S., Levinsohn, J., & Pakes, A. (1999). Voluntary export restraints on automobiles: Evaluating a trade policy. American Economic Review, 89(3), 400-431.  
* Conlon, C., & Gortmaker, J. (2020). Best practices for differentiated products demand estimation with pyblp. The RAND Journal of Economics, 51(4), 1108-1161.  
* Fukasawa, T. (2024). Fast and simple inner-loop algorithms of static/dynamic BLP estimations. arXiv preprint arXiv:2404.04494.  
* Nevo, A. (2001). Measuring market power in the ready-to-eat cereal industry. Econometrica, 69(2), 307-342.
