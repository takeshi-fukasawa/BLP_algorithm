function IV=IVS_compute_IV_func(IV_obs_pt,n_grid_V)

    IV_grid=construct_IV_grid_func(IV_obs_pt,n_grid_V);%1*ns*G*T*n_grid_V
   IV=cat(5,IV_obs_pt,IV_grid);%1*ns*G*T*n_dim_V
end
