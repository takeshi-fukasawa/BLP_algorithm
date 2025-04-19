function [numer_1,denom_1,IV,EV,coef_V,coef_0_AR1,coef_1_AR1,sigma_AR1]=compute_IV_EV_func(...
            V,u_ijt_tilde,beta_C,rho,weight_V,x_V);

    n_dim_V=size(V,4);
    T=size(u_ijt_tilde,4);
    n_grid_IV=n_dim_V-T;

    coef_V=[];coef_0_AR1=[];coef_1_AR1=[];sigma_AR1=[];

    numer_1=exp(u_ijt_tilde./(1-rho));%J*ns*G*T
    denom_1=sum(numer_1,1);%1*ns*G*T
    
    IV_obs_pt=(1-rho).*log(denom_1);%1*ns*G*T

    if n_grid_IV==0
        IV=IV_obs_pt;
    else
        IV=IVS_compute_IV_func(IV_obs_pt,n_dim_V-T);%1*ns*1*n_dim_V
    end

    if beta_C>0
        [EV,coef_V,coef_0_AR1,coef_1_AR1,sigma_AR1]=...
     compute_EV_func(V,IV,weight_V,x_V,T);%1*ns*1*n_dim_V
    else
        EV=0;
    end
    
end
    