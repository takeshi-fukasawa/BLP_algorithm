function IV_grid=construct_IV_grid_func(IV_obs_pt,n_grid_IV)
      min_ratio=0.7;
      diff_ratio=1-min_ratio;
      IV_grid=IV_obs_pt*(min_ratio+2*diff_ratio*reshape([0:n_grid_IV-1],1,1,1,1,n_dim_V)/(n_grid_IV-1));%1*ns*1*T*n_grid_V

end
