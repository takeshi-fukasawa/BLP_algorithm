function IV=IVS_compute_IV_func(IV_obs_pt,n_grid_V)

    [ns,T]=size(IV_obs_pt,[2,4]);
      n_grid_IV=n_grid_V;

    diff_ratio=0.2;%%%%
    n_dim_Chebyshev=n_grid_IV;


    if 1==0
     Chebyshev_extrema=cos([0:n_dim_Chebyshev-1]*pi/(n_dim_Chebyshev-1));

     temp=repmat(mean(IV_obs_pt,[2,4]),[1,ns,1,T,1]);
      IV_grid=repmat(mean(IV_obs_pt,[2,4]),[1,ns,1,T,1])+...
        diff_ratio*reshape(Chebyshev_extrema,1,1,1,1,n_grid_IV);%1*ns*1*T*n_grid_V; common grid for all t=1,...,T, i=1,...,ns
    end

      min_ratio=0.9;%%%%%%
      diff_ratio=1-min_ratio;
      IV_grid=repmat(mean(IV_obs_pt,[2,4]),[1,ns,1,T,1])*(min_ratio+2*diff_ratio*...
          reshape([0:n_grid_IV-1],1,1,1,1,n_grid_IV)/(n_grid_IV-1));%1*ns*1*T*n_grid_V; common grid for all t=1,...,T, i=1,...,ns

   IV=cat(5,IV_obs_pt,IV_grid);%1*ns*G*T*(n_dim_V=n_grid_IV+1)
end
