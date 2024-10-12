
function [output,other_vars]=...
    V_EV_update_func(...
    V_EV_initial,weight,mu_ijt,...
    S_jt_data,S_0t_data,...
    weight_V,x_V,beta_C,tune_param,Newton_spec)
 
    %%% G==1 & rho==0 case only
    [J,ns,~,T]=size(mu_ijt);
    n_dim_V_EV=size(V_EV_initial,4);
    n_grid_IV=n_dim_V_EV-T*2;

    V_initial_obs_pt=V_EV_initial(:,:,:,1:T);%1*ns*1*T; V at obs data pts
    V_initial=V_EV_initial(:,:,:,1:end-T);%1*ns*1*(T+n_grid_IV);
    EV_initial_obs_pt=V_EV_initial(:,:,:,end-T+1:end);%1*ns*1*T; EV at obs

    v_i0t_tilde_obs_pt=beta_C*EV_initial_obs_pt;

    
    %% Compute delta
    [delta_jt,Pr0]=compute_delta_from_V_func(...
       mu_ijt,weight,S_jt_data,V_initial_obs_pt,[],...
       []);

    if tune_param>0
        s_i0t_ccp_obs_pt=reshape(exp(v_i0t_tilde_obs_pt-V_initial_obs_pt),1,ns,1,T);

        s_0t_predict=...
            sum(reshape(weight,1,ns,1,1).*(1-Pr0+Pr0.*s_i0t_ccp_obs_pt),[2,5]);%1*1*1*T 
        S_0_ratio=s_0t_predict./reshape(S_0t_data,1,1,1,T);%1*1*1*T
    else
        s_i0t_ccp_obs_pt=[];
    end

    %% Update V
    rho=0;
    v_ijt_tilde=delta_jt+mu_ijt;

    [numer_1,denom_1,IV_new,EV]=compute_IV_EV_func(...
            V_initial,v_ijt_tilde,beta_C,rho,weight_V,x_V);

    v_i0t_tilde=beta_C*EV;

    if n_grid_IV==0
        V_updated=log(exp(v_i0t_tilde)+exp(IV_new));%1*ns*1*n_dim_V
    else
        if tune_param==0
            V_obs_pt=log(exp(v_i0t_tilde(:,:,:,1:T))+...
                exp(IV_new(:,:,:,1:T)));%1*ns*1*T
        else
            V_obs_pt=log(exp(v_i0t_tilde(:,:,:,1:T))+...
                exp(IV_new(:,:,:,1:T).*(S_0_ratio.^tune_param)));%1*ns*1*T
        end

        V_grid=log(exp(v_i0t_tilde(:,:,:,T+1:end))+...
            exp(IV_new(:,:,:,T+1:end)));%1*ns*1*n_grid_IV
        V_updated=cat(4,V_obs_pt,V_grid);%1*ns*1*n_dim_V
    end

        EV_updated_obs_pt=EV(:,:,:,end-T+1:end);%1*ns*1*T

        V_EV_updated=cat(4,V_updated,EV_updated_obs_pt);

    output={V_EV_updated};
    other_vars.s_i0t_ccp=s_i0t_ccp_obs_pt;
    other_vars.v_i0t_tilde=v_i0t_tilde;
    other_vars.IV=IV_new;
    other_vars.delta_jt=delta_jt;
    other_vars.v_ijt_tilde=v_ijt_tilde;
    
end
