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
