function IV=IVS_compute_IV_func(IV_obs_pt,n_grid_V)

    [ns,T]=size(IV_obs_pt,[2,4]);
    n_grid_IV=n_grid_V;

    diff_ratio=0.5;%%%%
    n_dim_Chebyshev=n_grid_IV;

    Chebyshev_extrema=cos([0:n_dim_Chebyshev-1]*pi/(n_dim_Chebyshev-1));

    %%%%%%%%%%%% large range important for convergence 
    IV_min=-4+min(IV_obs_pt,[],4);%1*ns*1*1
    IV_max=4+max(IV_obs_pt,[],4);%1*ns*1*1
    
    %%%% Fixed IV grid %%%%
    %%IV_min=IV_min*0-16;%% Performance rarely change??
    %%IV_max=IV_max*0+6;%% Performance rarely change??
    %%%%%%%%%%%%

    IV_grid=IV_min+...
        (IV_max-IV_min).*reshape((Chebyshev_extrema+1)/2,1,1,1,1,n_grid_IV);%1*ns*1*1*n_grid_IV
    IV_grid=repmat(IV_grid,[1,1,1,T,1]);%1*ns*1*T*n_grid_IV

   IV=cat(5,IV_obs_pt,IV_grid);%1*ns*G*T*(n_dim_V=n_grid_IV+1)
end
