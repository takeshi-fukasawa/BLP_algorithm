function resid_V=Bellman_update_func(...
    V_initial,delta,mu_ijt,beta_C,L,rho,weight)

global IV_temp0
%%% rho>0 allowed; Bellman update

    n_dim_V=size(V_initial,5);

    u_ijt_tilde=delta+mu_ijt+(beta_C^L).*V_initial;%J*I*G*T*n_dim_V

    numer_1=exp(u_ijt_tilde./(1-rho));%J*ns*G*T*n_dim_V
    denom_1=sum(numer_1,1);%1*ns*G*T*n_dim_V
    IV=(1-rho).*log(denom_1);%1*ns*G*T*n_dim_V
    IV_obs_pt=IV(:,:,:,:,1);%1*ns*G*T

    if n_dim_V==1
       IV=IV_obs_pt;
    else
        IV=IVS_compute_IV_func(IV_obs_pt,n_dim_V-1);%1*ns*1*T*n_dim_V
    end

    weight_V=[];
    x_V=[];
    EV=compute_EV_func(V_initial,IV,weight_V,x_V);%1*ns*1*T*n_dim_V
    u_i0t_tilde=beta_C*EV;%J*I*1*T*n_dim_V


    IV=log(sum(exp(u_ijt_tilde),1));%1*I*1*T*n_dim_V
    V_updated=log(exp(u_i0t_tilde)+sum(exp(IV),3));%1*ns*1*T*n_dim_V

    IV_temp0=IV;

resid_V={V_initial-V_updated};

end
