The script `check_Iaria_Wang_code.m` in the same folder is provided to validate and profile the replication code ("GenerateMarketShare_parametric" function) used in **Iaria and Wang (2025)**, because the efficiency of the implementation can significantly affect the relative performance of fixed-point iteration-based methods and Jacobian-based methods, as discussed in Supplemental Appendix S.9 of the current paper.

In the original MATLAB implementation, the choice probabilities were computed as follows:

```matlab
% Original implementation
ms_all = exp_utility_all * diag(1./(sum(exp_utility_all, 1) + 1));
```
However, the expression diag(...) entails allocating an $I \times I$ matrix, which becomes a major computational bottleneck and is highly memory-intensive for large $I$. The following implementation is significantly more efficient as it eliminates this overhead by utilizing element-wise division:


```matlab
ms_all = exp_utility_all./(sum(exp_utility_all,1)+1);
```

Code profiling (refer to profile_results/file2.html) indicates that this optimization alone provides a tenfold speedup in the market share calculation, fundamentally altering the assessment of the algorithm's overall efficiency.

## References
* Iaria, A., & Wang, A. (2025). Real analytic discrete choice models of demand: Theory and implications. Econometric Theory, 41(5), 1080-1128.
