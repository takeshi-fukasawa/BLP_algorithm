
function [output,other_vars]=...
    V_update_func(...
    V_initial,weight,mu_ijt,...
    S_jt_data,S_0t_data,...
    weight_V,x_V,beta_C,tune_param,Newton_spec)
 
    %%% G==1 & rho==0 case only
    [J,ns,~,T]=size(mu_ijt);
    n_dim_V=size(V_initial,5);
    n_state=size(weight,5);
   
    if n_state==1 & n_dim_V>1 %% IVS spec
       V_initial_obs_pt=V_initial(:,:,:,:,1);%1*ns*1*T; V at obs data pts
    else
       V_initial_obs_pt=V_initial;%1*ns*1*T*n_state
    end

    %% Compute delta
    exp_mu_ijt=exp(mu_ijt);
    s_jt_predict_up_to_scale=sum(reshape(weight,1,ns,1,1,n_state).*...
        reshape(exp_mu_ijt,J,ns,1,T).*...
        reshape(exp(-V_initial_obs_pt),1,ns,1,T,n_state),[2,5]);%J*1*1*T
    exp_delta_jt=S_jt_data./s_jt_predict_up_to_scale;%J*1*J*T
    

    %% Update V
    rho=0;

    [numer_1,denom_1,IV_new,EV]=compute_IV_EV_func(...
            V_initial,log(exp_delta_jt)+mu_ijt,beta_C,rho,weight_V,x_V);

    v_i0t_tilde=beta_C*EV;

    s_i0t_ccp=reshape(exp(v_i0t_tilde(:,:,:,:,1)-V_initial(:,:,:,:,1)),1,ns,1,T,n_state);
    s_0t_predict=sum(reshape(weight,1,ns,1,1,n_state).*s_i0t_ccp,[2,5]);%1*1*1*T %%%Pr0 not incorporated???

    S_0_ratio=s_0t_predict./reshape(S_0t_data,1,1,1,T);%1*1*1*T

    V_updated=log(exp(v_i0t_tilde)+exp(IV_new).*(S_0_ratio.^tune_param));%1*ns*1*T*n_dim_V

%%%%%%%%%%%%%%%%%%%%%%%
    if Newton_spec==1 & n_grid_V==1
        %%% Newton iteration
        prob_0=reshape(exp(v_i0t_tilde),1,ns,1,T)./reshape(exp(V_initial),1,ns,1,T);%1*ns*1*T
        diff=log(exp(v_i0t_tilde)+temp_1.*(temp_2))-V_initial;
        V_updated=V_initial+diff./(1-beta_C*prob_0);%1*ns*1*T
    end
    %%%%%%%%%%%%%%%%%%%%%%%%

    resid=V_initial-V_updated;

    output={resid};
    other_vars.s_i0t_ccp=s_i0t_ccp;
    other_vars.v_i0t_tilde=v_i0t_tilde;
    other_vars.IV=IV_new;
    other_vars.exp_delta_jt=exp_delta_jt;
    
end
