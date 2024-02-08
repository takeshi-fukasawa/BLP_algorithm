
function [output,other_vars]=...
    V_update_func(...
    V_initial,weight,mu_ijt,rho,...
    S_jt_data,S_0t_data,...
    weight_V,x_V,beta_C,L,tune_param,Newton_spec)
 
    %%% G==1 & rho==0 case only
    [J,ns,~,T]=size(mu_ijt);
    n_dim_V=size(V_initial,5);
    n_state=size(weight,5);
   
    if n_state==1 & n_dim_V>1 %% IVS spec
       V_initial_obs_pt=V_initial(:,:,:,:,1);%1*ns*1*T; V at obs data pts
    else
       V_initial_obs_pt=V_initial;%1*ns*1*T*n_state
    end

    v_ijt_tilde=mu_ijt;%J*ns*1*T*n_state %%%%%%%% add beta_C*EV %%%%%%%
    exp_v_ijt_tilde=exp(v_ijt_tilde);%J*ns*1*T
    denom_1=sum(reshape(weight,1,ns,1,1,n_state).*reshape(exp_v_ijt_tilde,J,ns,1,T).*...
        reshape(exp(-V_initial_obs_pt),1,ns,1,T,n_state),[2,5]);%J*1*1*T
    exp_delta_jt=S_jt_data./denom_1;%J*1*J*T
    exp_IV_new=sum(exp_delta_jt.*exp_v_ijt_tilde,1);%1*ns*1*T*n_state

if n_state>1 & n_dim_V==1 % Setting other than IVS 
    IV_new=log(exp_IV_new);%1*ns*1*T*n_state
else
%%% Inclusive value sufficiency
    IV_new_obs_pt=log(exp_IV_new);%1*1*G*T
    IV_new=IVS_compute_IV_func(IV_new_obs_pt,n_dim_V-1);%1*ns*1*T*n_dim_V
    exp_IV_new=exp(IV_new);%1*ns*1*T*n_dim_V
end

    EV=compute_EV_func(V_initial,IV_new,weight_V,x_V);%1*ns*1*T*n_dim_V
    v_i0t_tilde=beta_C*EV;%1*ns*1*T*n_dim_V

    v_i0t_tilde_obs_pt=v_i0t_tilde(:,:,:,:,1);%1*ns*1*T*n_dim_V; v_i0t_tilde at obs data pts

    s_i0t_ccp=reshape(exp(v_i0t_tilde_obs_pt-V_initial_obs_pt),1,ns,1,T,n_state);
    s_0t_predict=sum(reshape(weight,1,ns,1,1,n_state).*s_i0t_ccp,[2,5]);%1*1*1*T %%%Pr0 not incorporated???

    S_0_ratio=s_0t_predict./reshape(S_0t_data,1,1,1,T);%1*1*1*T


        V_updated=log(exp(v_i0t_tilde)+exp_IV_new.*(S_0_ratio.^tune_param));%1*ns*1*T*n_dim_V

    if Newton_spec==1 & n_grid_V==1
        %%% Newton iteration
        prob_0=reshape(exp(v_i0t_tilde),1,ns,1,T)./reshape(exp(V_initial),1,ns,1,T);%1*ns*1*T
        diff=log(exp(v_i0t_tilde)+temp_1.*(temp_2))-V_initial;
        V_updated=V_initial+diff./(1-beta_C*prob_0);%1*ns*1*T
    end
    
    resid=V_initial-V_updated;

    output={resid};
    other_vars.s_i0t_ccp=s_i0t_ccp;

end
