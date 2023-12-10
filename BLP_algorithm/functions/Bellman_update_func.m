function resid_V=Bellman_update_func(...
    V_initial,delta,mu_ijt,beta_C,L,rho,weight)

    n_dim_V=size(V_initial,5);

    u_ijt_tilde=delta+mu_ijt+(beta_C^L).*V_initial;%J*I*G*T

    numer_1=exp(u_ijt_tilde./(1-rho));%J*ns*G*T
    denom_1=sum(numer_1,1);%1*ns*G*T
    IV_obs_pt=(1-rho).*log(denom_1);%1*ns*G*T

    if n_dim_V==1
       IV=IV_obs_pt;
    else
        IV=IVS_compute_IV_func(IV_obs_pt,n_dim_V-1);%1*ns*1*T*n_dim_V
    end

    weight_V=[];
    EV=compute_EV_func(V_initial,IV,weight_V);%1*ns*1*T*n_dim_V
    u_i0t_tilde=beta_C*EV;%J*I*1*T*n_dim_V


    IV=log(sum(exp(u_ijt_tilde),[1,3]));
    V_updated=exp(u_i0t_tilde)+sum(exp(IV),3);%1*ns*1*T*n_dim_V


resid_V={V_initial-V_updated};

end
