function IV=IVS_compute_IV_func(IV_obs_pt,n_grid_V)

      n_grid_IV=n_grid_V;

      min_ratio=0.9;%%%%%%
      diff_ratio=1-min_ratio;
      IV_grid=IV_obs_pt*(min_ratio+2*diff_ratio*...
          reshape([0:n_grid_IV-1],1,1,1,1,n_grid_IV)/(n_grid_IV-1));%1*ns*1*T*n_grid_V

   IV=cat(5,IV_obs_pt,IV_grid);%1*ns*G*T*(n_dim_V=n_grid_IV+1)
end
