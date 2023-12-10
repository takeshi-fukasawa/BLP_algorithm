
function output=...
    V_update_func(...
    V_initial,weight,mu_ijt,rho,...
    S_jt_data,S_0t_data,numer_1_without_delta,...
    beta_C,L,tune_param,Newton_spec)
 
    %%% G==1 case only
    [J,ns,~,T]=size(mu_ijt);
    n_dim_V=size(V_initial,5);
   
   
    V_initial_obs_pt=V_initial(:,:,:,:,1);%1*ns*1*T; V at obs data pts
      
    numer_1=exp(mu_ijt);%J*ns*1*T
    denom_1=sum(reshape(weight,1,ns,1,1).*reshape(numer_1,J,ns,1,T).*...
        reshape(exp(-V_initial_obs_pt),1,ns,1,T),2);%J*1*1*T
    exp_IV_new=sum(S_jt_data.*numer_1./denom_1,1);%1*ns*1*T
     
if n_dim_V==1
    IV_new=log(exp_IV_new);
else
%%% Inclusive value sufficiency
IV_new_obs_pt=log(exp_IV_new);%J*1*1*T
IV_new=IVS_compute_IV_func(IV_new_obs_pt,n_dim_V-1);%1*ns*1*T*n_dim_V
exp_IV_new=exp(IV_new);%1*ns*1*T*n_dim_V
end

    weight_V=[];
    EV=compute_EV_func(V_initial,IV_new,weight_V);%1*ns*1*T*n_dim_V
    v_i0t=beta_C*EV;%1*ns*1*T*n_dim_V

    v_i0t_obs_pt=v_i0t(:,:,:,:,1);%1*ns*1*T*n_dim_V; v_i0t at obs data pts

    numer_2=sum(reshape(weight,1,ns).*reshape(exp(v_i0t_obs_pt-V_initial_obs_pt),1,ns,1,T),2);%1*1*1*T
    S_0_ratio=numer_2./reshape(S_0t_data,1,1,1,T);%1*1*1*T


        V_updated=log(exp(v_i0t)+exp_IV_new.*(S_0_ratio.^tune_param));%1*ns*1*T*n_dim_V

    if Newton_spec==1 & n_grid_V==1
        %%% Newton iteration
        prob_0=reshape(exp(v_i0t),1,ns,1,T)./reshape(exp(V_initial),1,ns,1,T);%1*ns*1*T
        diff=log(exp(v_i0t)+temp_1.*(temp_2))-V_initial;
        V_updated=V_initial+diff./(1-beta_C*prob_0);%1*ns*1*T
    end
    
    resid=V_initial-V_updated;

    output={resid};
end
