def results_func(results):
   return [sum(sum(results.cumulative_contraction_evaluations))/
           (results.cumulative_objective_evaluations*results.problem.T),### mean contraction evaluation per market 
           results.cumulative_optimization_time/
                 (results.cumulative_objective_evaluations*results.problem.T),## Mean comp time per market
            results.cumulative_objective_evaluations,
           results.cumulative_optimization_time,
            sum(sum(results.cumulative_comp_time_solve_delta)),
           sum(sum(results.cumulative_contraction_evaluations)),
           results.objective[0][0]
           ]

