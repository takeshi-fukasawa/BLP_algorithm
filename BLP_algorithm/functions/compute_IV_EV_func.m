function [numer_1,denom_1,IV,EV]=compute_IV_EV_func(...
            V,u_ijt_tilde,beta_C,rho,weight_V,x_V);

    n_dim_V=size(V,5);


    if n_dim_V==1
        numer_1=exp(u_ijt_tilde./(1-rho));%J*ns*G*T*n_dim_V
        denom_1=sum(numer_1,1);%1*ns*G*T*n_dim_V
        IV=(1-rho).*log(denom_1);%1*ns*G*T*n_dim_V

    if n_dim_V>=2
        IV_obs_pt=IV(:,:,:,:,1);%1*ns*G*T
        IV=IVS_compute_IV_func(IV_obs_pt,n_dim_V-1);%1*ns*1*T*n_dim_V
    end

    if beta_C>0
        EV=compute_EV_func(V,IV,weight_V,x_V);%1*ns*1*T*n_dim_V
    else
        EV=0;
    end
    
end
    